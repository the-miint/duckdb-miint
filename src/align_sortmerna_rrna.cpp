#include "align_sortmerna_rrna.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

namespace duckdb {

namespace {

// Parse a LIST(VARCHAR) named parameter into std::vector<std::string>.
// Throws BinderException on missing / empty / wrong-type inputs.
std::vector<std::string> ParseRefPaths(const named_parameter_map_t &params) {
	auto it = params.find("ref_paths");
	if (it == params.end() || it->second.IsNull()) {
		throw BinderException("align_sortmerna_rrna requires ref_paths (LIST of FASTA paths)");
	}
	auto &children = ListValue::GetChildren(it->second);
	if (children.empty()) {
		throw BinderException("align_sortmerna_rrna: ref_paths must be a non-empty list");
	}
	std::vector<std::string> result;
	result.reserve(children.size());
	for (auto &child : children) {
		if (child.IsNull()) {
			throw BinderException("align_sortmerna_rrna: ref_paths contains NULL entry");
		}
		result.push_back(child.ToString());
	}
	return result;
}

void ParseSortMeRNAConfig(const named_parameter_map_t &params, miint::SortMeRNAConfig &cfg) {
	auto set_i32 = [&](const char *name, int32_t &dst) {
		auto it = params.find(name);
		if (it != params.end() && !it->second.IsNull()) {
			dst = it->second.GetValue<int32_t>();
		}
	};
	auto set_u32 = [&](const char *name, uint32_t &dst) {
		auto it = params.find(name);
		if (it != params.end() && !it->second.IsNull()) {
			dst = it->second.GetValue<uint32_t>();
		}
	};
	auto set_bool = [&](const char *name, bool &dst) {
		auto it = params.find(name);
		if (it != params.end() && !it->second.IsNull()) {
			dst = it->second.GetValue<bool>();
		}
	};
	auto set_double = [&](const char *name, double &dst) {
		auto it = params.find(name);
		if (it != params.end() && !it->second.IsNull()) {
			dst = it->second.GetValue<double>();
		}
	};

	set_i32("num_threads", cfg.num_threads);
	set_i32("match", cfg.match);
	set_i32("mismatch", cfg.mismatch);
	set_i32("gap_open", cfg.gap_open);
	set_i32("gap_ext", cfg.gap_ext);
	set_i32("score_N", cfg.score_N);
	set_double("evalue", cfg.evalue);
	set_u32("seed_win_len", cfg.seed_win_len);
	set_u32("num_alignments", cfg.num_alignments);
	set_bool("best", cfg.best);
	set_bool("paired", cfg.paired);
	set_bool("forward_only", cfg.forward_only);
	set_bool("reverse_only", cfg.reverse_only);
	set_bool("full_search", cfg.full_search);
}

} // namespace

unique_ptr<FunctionData> AlignSortMeRNARRNATableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                               vector<LogicalType> &return_types,
                                                               vector<std::string> &names) {
	auto data = make_uniq<Data>();

	if (input.inputs.empty() || input.inputs[0].IsNull()) {
		throw BinderException("align_sortmerna_rrna: query_table is required");
	}
	data->query_table = input.inputs[0].ToString();
	data->ref_paths = ParseRefPaths(input.named_parameters);
	ParseSortMeRNAConfig(input.named_parameters, data->config);

	data->query_schema = ValidateSequenceTableSchema(context, data->query_table);
	if (data->query_schema.has_sequence2 != data->config.paired) {
		throw BinderException("align_sortmerna_rrna: query table paired-ness (sequence2 %s) does not match "
		                      "paired=%s; set the paired parameter to match or reshape the query",
		                      data->query_schema.has_sequence2 ? "present" : "absent",
		                      data->config.paired ? "true" : "false");
	}

	for (const auto &n : data->names)
		names.emplace_back(n);
	for (const auto &t : data->types)
		return_types.emplace_back(t);
	return data;
}

unique_ptr<GlobalTableFunctionState> AlignSortMeRNARRNATableFunction::InitGlobal(ClientContext &context,
                                                                                 TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	miint::SortMeRNAConfig cfg = data.config;
	if (cfg.num_threads <= 0) {
		cfg.num_threads = NumericCast<int32_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	}
	gstate->aligner = std::make_unique<miint::SortMeRNAAligner>(cfg, data.ref_paths);
	gstate->query_stream = std::make_unique<QuerySequenceStream>(context, data.query_table, data.query_schema);
	// Execute holds gstate.lock across align() on the assumption that
	// MaxThreads() is 1. If that invariant ever changes, the shared
	// result_buffer's deposit path needs a rethink before this body is run
	// concurrently.
	D_ASSERT(gstate->MaxThreads() == 1);
	return gstate;
}

unique_ptr<LocalTableFunctionState>
AlignSortMeRNARRNATableFunction::InitLocal(ExecutionContext &, TableFunctionInitInput &, GlobalTableFunctionState *) {
	return make_uniq<LocalState>();
}

void AlignSortMeRNARRNATableFunction::Execute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();
	auto &bind_data = data_p.bind_data->Cast<Data>();

	// MaxThreads() == 1 means only one DuckDB thread ever reaches this body;
	// the lock is defensive against future changes. Holding it across the
	// align() call is deliberate — the streaming invariant (one sub-batch in
	// flight at a time, shared buffer consumed in order) requires it, and
	// sortmerna's g_run_mutex would serialize concurrent align() calls anyway.
	std::lock_guard<std::mutex> lock(gstate.lock);

	while (true) {
		// 1. Emit buffered rows first.
		idx_t available = gstate.result_buffer.size() - gstate.buffer_offset;
		if (available > 0) {
			idx_t count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			OutputSortMeRNARRNABatch(output, gstate.result_buffer, gstate.buffer_offset, count);
			gstate.buffer_offset += count;
			return;
		}

		// 2. Buffer exhausted — pull next sub-batch from the query stream.
		gstate.result_buffer.clear();
		gstate.buffer_offset = 0;

		auto query_batch = gstate.query_stream->FetchSubBatch();
		if (query_batch.empty()) {
			output.SetCardinality(0);
			return;
		}

		// 3. Convert SequenceRecordBatch → SortMeRNAQueryBatch.
		miint::SortMeRNAQueryBatch queries;
		queries.read_ids = std::move(query_batch.read_ids);
		queries.sequences = std::move(query_batch.sequences1);
		if (bind_data.config.paired) {
			queries.sequences2 = std::move(query_batch.sequences2);
		}

		// 4. Align. Rethrow library errors as IOException so DuckDB's
		//    table-function error path handles them cleanly rather than
		//    letting std::runtime_error cross the ABI boundary.
		try {
			gstate.aligner->align(queries, gstate.result_buffer);
		} catch (const std::exception &e) {
			throw IOException("align_sortmerna_rrna: %s", e.what());
		}
	}
}

TableFunction AlignSortMeRNARRNATableFunction::GetFunction() {
	auto tf = TableFunction("align_sortmerna_rrna", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["ref_paths"] = LogicalType::LIST(LogicalType::VARCHAR);
	tf.named_parameters["num_threads"] = LogicalType::INTEGER;
	tf.named_parameters["match"] = LogicalType::INTEGER;
	tf.named_parameters["mismatch"] = LogicalType::INTEGER;
	tf.named_parameters["gap_open"] = LogicalType::INTEGER;
	tf.named_parameters["gap_ext"] = LogicalType::INTEGER;
	tf.named_parameters["score_N"] = LogicalType::INTEGER;
	tf.named_parameters["evalue"] = LogicalType::DOUBLE;
	tf.named_parameters["seed_win_len"] = LogicalType::UINTEGER;
	tf.named_parameters["num_alignments"] = LogicalType::UINTEGER;
	tf.named_parameters["best"] = LogicalType::BOOLEAN;
	tf.named_parameters["paired"] = LogicalType::BOOLEAN;
	tf.named_parameters["forward_only"] = LogicalType::BOOLEAN;
	tf.named_parameters["reverse_only"] = LogicalType::BOOLEAN;
	tf.named_parameters["full_search"] = LogicalType::BOOLEAN;

	// Alignment output order depends on sortmerna's internal batch traversal,
	// not on the input row order. Declare NO_ORDER so DuckDB doesn't try to
	// preserve insertion order for CTAS pipelines.
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void AlignSortMeRNARRNATableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
