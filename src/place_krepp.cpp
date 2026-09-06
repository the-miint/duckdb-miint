#include "place_krepp.hpp"

#include "alignment_functions_internal.hpp"
#include "id_column_utils.hpp"
#include "miint_log.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

#include <algorithm>
#include <functional>
#include <stdexcept>

namespace duckdb {

namespace {

// Reads an optional named parameter, leaving `target` alone when absent.
template <typename T>
void ReadOptional(const named_parameter_map_t &params, const char *key, T &target,
                  const std::function<T(const Value &)> &convert) {
	auto it = params.find(key);
	if (it != params.end() && !it->second.IsNull()) {
		target = convert(it->second);
	}
}

} // namespace

unique_ptr<FunctionData> PlaceKreppTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                       vector<LogicalType> &return_types, vector<std::string> &names) {
	auto data = make_uniq<Data>();

	auto query_param = input.named_parameters.find("query_table");
	if (query_param == input.named_parameters.end() || query_param->second.IsNull()) {
		throw BinderException("place_krepp requires a query_table parameter");
	}
	data->query_table = query_param->second.ToString();

	auto index_param = input.named_parameters.find("index_path");
	if (index_param == input.named_parameters.end() || index_param->second.IsNull()) {
		throw BinderException("place_krepp requires an index_path parameter");
	}
	data->index_path = index_param->second.ToString();

	ReadOptional<std::string>(input.named_parameters, "newick_path", data->newick_path,
	                          [](const Value &v) { return v.ToString(); });
	ReadOptional<uint32_t>(input.named_parameters, "hdist_th", data->config.hdist_th,
	                       [](const Value &v) { return v.GetValue<uint32_t>(); });
	ReadOptional<uint32_t>(input.named_parameters, "tau", data->config.tau,
	                       [](const Value &v) { return v.GetValue<uint32_t>(); });
	ReadOptional<double>(input.named_parameters, "chisq", data->config.chisq,
	                     [](const Value &v) { return v.GetValue<double>(); });
	ReadOptional<bool>(input.named_parameters, "multi", data->config.multi,
	                   [](const Value &v) { return BooleanValue::Get(v); });
	ReadOptional<bool>(input.named_parameters, "filter", data->config.filter,
	                   [](const Value &v) { return BooleanValue::Get(v); });

	// Both of these re-do a validation krepp performs in its CLI layer, which
	// this build excludes (ext/krepp/src/krepp.hpp holds
	// validate_configuration_place, ext/krepp/src/krepp.cpp:844 the --chisq
	// CLI::PositiveNumber check). Without them the library path is entered with
	// arguments krepp itself would have refused.
	//
	// tau > hdist_th is memory-unsafe, not merely useless: Minfo sizes
	// hdisthist_v to hdist_th + 1 (query.hpp:144-146) and get_leq_tau sums
	// hdisthist_v[0..tau] with unchecked operator[] (query.hpp:200-206).
	if (data->config.hdist_th < data->config.tau) {
		throw BinderException("place_krepp: tau (%d) must not exceed hdist_th (%d)", data->config.tau,
		                      data->config.hdist_th);
	}
	// hdist_th has no upper bound of its own in krepp - the CLI checks only
	// NonNegativeNumber - and the library path is memory-unsafe above k.
	// HDistHistLLH's constructor (ext/krepp/src/hdhistllh.hpp:58-67) sizes
	// binom_coef_k to k + 1 and then reads binom_coef_k[i] for i up to
	// hdist_th, and sizes binom_coef_hnk to hdist_th + 1, which wraps to 0 at
	// UINT32_MAX and makes the very next write out of bounds. Measured: 64
	// returns a different row count off the end of the array, and 4294967295
	// segfaults the DuckDB process outright.
	//
	// 31 is krepp's own hard maximum for k (krepp.hpp:72-74 refuses to build an
	// index above it), so nothing past that can ever be valid and it is
	// rejected here rather than after loading what may be tens of gigabytes.
	// The exact bound needs the index's k and is applied in InitGlobal.
	if (data->config.hdist_th > 31) {
		throw BinderException("place_krepp: hdist_th (%d) exceeds the largest k-mer length krepp supports (31); "
		                      "it must not exceed the index's k",
		                      data->config.hdist_th);
	}
	// Written as !(x > 0) so NaN is rejected too. A non-positive cutoff makes
	// krepp's `chisq < chisq_value` unsatisfiable, so every candidate edge is
	// discarded and the function returns nothing - indistinguishable from "these
	// reads do not place".
	if (!(data->config.chisq > 0.0)) {
		throw BinderException("place_krepp: chisq must be a positive number (got %s)",
		                      Value::DOUBLE(data->config.chisq).ToString());
	}

	// Same contract as the aligners: read_id plus sequence1, BIGINT ids allowed.
	data->schema = ValidateSequenceTableSchema(context, data->query_table, /*allow_bigint=*/true);

	// fragment mirrors the query relation's read_id type, as align_minimap2 and
	// align_sortmerna do. Pinning it to VARCHAR would make a BIGINT read_id come
	// back as text and force every join back to the query table through a cast.
	data->types[0] = data->schema.id_type;

	return_types = data->types;
	names = data->names;
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> PlaceKreppTableFunction::InitGlobal(ClientContext &context,
                                                                         TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	// Loaded once and shared: the query path only reads it, and indexes are
	// routinely far too large to hold per thread.
	//
	// KreppPlacer is built into the C++ unit-test binary, which does not link
	// DuckDB, so it raises miint::InvalidInputException and std::runtime_error
	// instead of DuckDB's types. Translating here keeps the error classes a SQL
	// user sees meaningful.
	try {
		gstate->index = std::make_shared<miint::SharedKreppIndex>(bind_data.index_path, bind_data.newick_path);
	} catch (const miint::InvalidInputException &e) {
		throw InvalidInputException(std::string(e.what()));
	} catch (const std::runtime_error &e) {
		throw IOException(std::string(e.what()));
	}
	// The exact bound, now that k is known. Above it krepp reads past the end of
	// binom_coef_k; see the note on the bind-time check.
	const uint32_t k = gstate->index->kmer_length();
	if (bind_data.config.hdist_th > k) {
		throw InvalidInputException("place_krepp: hdist_th (%d) must not exceed the index's k-mer length (%d)",
		                            bind_data.config.hdist_th, k);
	}
	// A backbone tip the index does not know contributes no placements, and
	// krepp says nothing about it. Sharing none of them is rejected in the
	// SharedKreppIndex constructor; a partial overlap is legitimate - a backbone
	// covering more than this index - but it silently shrinks the space a read
	// can be placed into, so it is worth one line. InitGlobal runs once per
	// scan, so this needs no warn-once guard.
	const size_t tips_total = gstate->index->backbone_tips_total();
	const size_t tips_matched = gstate->index->backbone_tips_matched();
	if (tips_total > 0 && tips_matched < tips_total) {
		miint::EmitWarning(context, "place_krepp: only " + std::to_string(tips_matched) + " of " +
		                                std::to_string(tips_total) +
		                                " backbone tree tips are present in the index; placements can only land "
		                                "on the matched part of the tree");
	}
	gstate->stream = make_uniq<QuerySequenceStream>(context, bind_data.query_table, bind_data.schema);
	gstate->max_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	return std::move(gstate);
}

unique_ptr<LocalTableFunctionState> PlaceKreppTableFunction::InitLocal(ExecutionContext &context,
                                                                       TableFunctionInitInput &input,
                                                                       GlobalTableFunctionState *gstate_p) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto &gstate = gstate_p->Cast<GlobalState>();
	auto lstate = make_uniq<LocalState>();
	lstate->placer = make_uniq<miint::KreppPlacer>(gstate.index, bind_data.config);
	return std::move(lstate);
}

void PlaceKreppTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();
	auto &lstate = data_p.local_state->Cast<LocalState>();

	// Refill from the query stream until there is something to emit, or the
	// stream runs dry. A batch can place nothing at all, so this loops rather
	// than returning an empty chunk that DuckDB would read as end-of-input.
	while (lstate.emitted >= lstate.pending.size()) {
		auto batch = gstate.stream->FetchSubBatch();
		if (batch.size() == 0) {
			output.SetCardinality(0);
			return;
		}
		std::vector<miint::KreppQuery> queries;
		queries.reserve(batch.size());
		for (size_t i = 0; i < batch.size(); ++i) {
			queries.push_back({batch.read_ids[i], batch.sequences1[i]});
		}
		lstate.pending.clear();
		lstate.emitted = 0;
		size_t skipped = 0;
		try {
			skipped = lstate.placer->place(queries, lstate.pending);
		} catch (const miint::InvalidInputException &e) {
			throw InvalidInputException(std::string(e.what()));
		} catch (const std::runtime_error &e) {
			throw IOException(std::string(e.what()));
		}
		// Silently dropping these would leave a user whose reads are all
		// shorter than k staring at an empty result with nothing to explain it.
		//
		// No count is quoted. `skipped` is this batch's tally, and other workers
		// are still placing when the first one warns, so any number printed here
		// would be a floor presented as a total - worse than no number at all.
		if (skipped > 0) {
			bool already_warned = false;
			if (gstate.warned_short.compare_exchange_strong(already_warned, true)) {
				miint::EmitWarning(context, "place_krepp: skipping query sequence(s) shorter than the index "
				                            "k-mer length (" +
				                                std::to_string(gstate.index->kmer_length()) +
				                                "); these cannot be placed");
			}
		}
	}

	const idx_t count = std::min<idx_t>(STANDARD_VECTOR_SIZE, lstate.pending.size() - lstate.emitted);
	// Column-wise, like the aligner output paths in align_common.hpp: SetValue
	// boxes every cell into a Value first. EmitIdCell is the shared id codec, so
	// a BIGINT or UUID fragment is written the same way align_minimap2 writes it.
	auto &fragment_out = output.data[0];
	auto edge_num_out = FlatVector::GetData<int64_t>(output.data[1]);
	auto likelihood_out = FlatVector::GetData<double>(output.data[2]);
	auto lwr_out = FlatVector::GetData<double>(output.data[3]);
	auto distal_out = FlatVector::GetData<double>(output.data[4]);
	auto pendant_out = FlatVector::GetData<double>(output.data[5]);
	auto distance_out = FlatVector::GetData<double>(output.data[6]);
	for (idx_t row = 0; row < count; ++row) {
		const auto &placement = lstate.pending[lstate.emitted + row];
		EmitIdCell(fragment_out, row, placement.fragment, bind_data.schema.id_type);
		edge_num_out[row] = placement.edge_num;
		likelihood_out[row] = placement.likelihood;
		lwr_out[row] = placement.like_weight_ratio;
		distal_out[row] = placement.distal_length;
		pendant_out[row] = placement.pendant_length;
		distance_out[row] = placement.distance;
	}
	lstate.emitted += count;
	output.SetCardinality(count);
}

void PlaceKreppTableFunction::Register(ExtensionLoader &loader) {
	TableFunction place("place_krepp", {}, Execute, Bind, InitGlobal, InitLocal);
	place.named_parameters["query_table"] = LogicalType::VARCHAR;
	place.named_parameters["index_path"] = LogicalType::VARCHAR;
	place.named_parameters["newick_path"] = LogicalType::VARCHAR;
	place.named_parameters["hdist_th"] = LogicalType::UINTEGER;
	place.named_parameters["tau"] = LogicalType::UINTEGER;
	place.named_parameters["chisq"] = LogicalType::DOUBLE;
	place.named_parameters["multi"] = LogicalType::BOOLEAN;
	place.named_parameters["filter"] = LogicalType::BOOLEAN;
	loader.RegisterFunction(place);
}

} // namespace duckdb
