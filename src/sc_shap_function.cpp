#include "sc_shap_function.hpp"

#include "catalog_utils.hpp"
#include "miint_log.hpp"
#include "sc_common.hpp"
#include "sc_coo_builder.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace duckdb {

namespace {

//! Default budget of attributions one batch may hold: samples x classes x
//! features.
//!
//! sc computes every attribution in a batch as a dense double before anything is
//! filtered, so this bounds real memory rather than output size: 10M
//! attributions is ~80 MB. That buffer is the only full copy -- sc writes each
//! sample into it in place, hands it to Arrow without copying, and rows are read
//! straight out of it -- so 80 MB plus per-thread scratch is the peak however
//! many samples are explained. At 5,000 features that is 2,000 samples a batch;
//! at 200,000 features, 50.
constexpr int64_t kDefaultMaxAttributions = 10000000;

struct ScShapData : public TableFunctionData {
	string data_relation;
	string model_relation;
	string model_name;
	bool classification = false;
	bool predicted_class_only = true;
	//! 0 means every feature.
	int64_t top_k = 0;
	int64_t max_attributions = kDefaultMaxAttributions;
	bool max_attributions_set = false;
	//! Samples per batch; 0 means derive it from max_attributions.
	int64_t batch_size = 0;
	int32_t n_threads = 0;
	//! The data's id types. Its feature_id type equals the model's, which bind
	//! enforces, so it also types the feature_id output.
	miint::ScIdTypes id_types;
};

struct ScShapGlobalState : public GlobalTableFunctionState {
	// Held for the whole scan: every batch is explained by the same context and
	// model against the same scanned input. Declared in dependency order, so the
	// batcher is destroyed before the table it indexes.
	miint::ScContext ctx;
	miint::ScModel model;
	//! Every sample; sample ids and coverage are read from here.
	std::unique_ptr<miint::ScCooTable> table;
	std::unique_ptr<miint::ScCooBatcher> batcher;
	std::vector<std::string> feature_ids;
	//! Empty for a regressor.
	std::vector<std::string> classes;
	//! One per output. Comes from the trees alone, so it is the same every batch.
	std::vector<double> base_values;
	//! Predicted output per sample; unused for a regressor.
	std::vector<uint32_t> predicted;
	size_t n_outputs = 1;
	size_t n_features = 0;
	size_t batch_size = 1;

	// The batch being emitted. Only one batch's attributions exist at a time.
	//! sc's attribution export for the current batch, read in place.
	std::unique_ptr<miint::OwnedArrowArray> shap_array;
	//! Row-major batch_count x n_outputs x n_features, output-major within a
	//! sample. Points into shap_array.
	const double *values = nullptr;
	size_t batch_first = 0;
	size_t batch_count = 0;

	// Emission cursor over all samples. Rows are produced on demand rather than
	// materialised, so the output costs nothing beyond the current batch.
	size_t cur_sample = 0;
	size_t cur_output = 0;
	std::vector<uint32_t> chosen;
	size_t next_chosen = 0;
	bool have_row = false;
	std::vector<uint32_t> scratch_pos;
	std::vector<uint32_t> scratch_neg;

	bool loaded = false;
};

unique_ptr<FunctionData> ScShapBind(ClientContext &context, TableFunctionBindInput &input,
                                    vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<ScShapData>();
	data->data_relation = input.inputs[0].GetValue<string>();
	data->model_relation = input.inputs[1].GetValue<string>();
	if (data->data_relation.empty() || data->model_relation.empty()) {
		throw InvalidInputException("shap: data and model relation names must not be empty");
	}
	for (auto &kv : input.named_parameters) {
		const auto &k = kv.first;
		const auto &v = kv.second;
		if (v.IsNull()) {
			throw InvalidInputException("shap: named parameter '%s' must not be NULL", k);
		}
		if (StringUtil::CIEquals(k, "name")) {
			data->model_name = v.GetValue<string>();
		} else if (StringUtil::CIEquals(k, "predicted_class_only")) {
			data->predicted_class_only = v.GetValue<bool>();
		} else if (StringUtil::CIEquals(k, "top_k")) {
			data->top_k = v.GetValue<int64_t>();
			if (data->top_k <= 0) {
				throw InvalidInputException("shap: top_k must be > 0 (got %lld); omit it to return every feature",
				                            (int64_t)data->top_k);
			}
		} else if (StringUtil::CIEquals(k, "max_attributions")) {
			data->max_attributions = v.GetValue<int64_t>();
			data->max_attributions_set = true;
			if (data->max_attributions <= 0) {
				throw InvalidInputException("shap: max_attributions must be > 0 (got %lld)",
				                            (int64_t)data->max_attributions);
			}
		} else if (StringUtil::CIEquals(k, "batch_size")) {
			data->batch_size = v.GetValue<int64_t>();
			if (data->batch_size <= 0) {
				throw InvalidInputException("shap: batch_size must be > 0 (got %lld)", (int64_t)data->batch_size);
			}
		} else if (StringUtil::CIEquals(k, "n_threads")) {
			data->n_threads = v.GetValue<int32_t>();
		}
	}
	if (data->batch_size > 0 && data->max_attributions_set) {
		throw InvalidInputException("shap: pass batch_size or max_attributions, not both -- batch_size fixes the "
		                            "samples per batch, max_attributions derives it from a memory budget");
	}
	{
		auto conn = MakeReadOnlyHelperConnection(context);
		data->id_types = miint::ProbeIdTypes(conn, data->data_relation, "shap");
		data->classification =
		    miint::ReadModelTask(conn, data->model_relation, data->model_name, "shap") == "classification";
		miint::RequireModelFeatureIdType(
		    miint::ReadModelFeatureIdType(conn, data->model_relation, data->model_name, "shap"),
		    data->id_types.feature_id, "shap", data->model_relation, data->model_name, data->data_relation);
	}

	names = {"sample_id", "class", "feature_id", "shap_value", "base_value", "sample_coverage"};
	return_types = {data->id_types.sample_id, LogicalType::VARCHAR, data->id_types.feature_id,
	                LogicalType::DOUBLE,      LogicalType::DOUBLE,  LogicalType::DOUBLE};
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ScShapInitGlobal(ClientContext &, TableFunctionInitInput &) {
	return make_uniq<ScShapGlobalState>();
}

void ScanForShap(Connection &conn, const ScShapData &bind, miint::ScCooBuilder &builder) {
	const auto q = KeywordHelper::WriteOptionallyQuoted(bind.data_relation);
	// The casts guarantee the physical layout the buffer reads below assume;
	// see the note in sc_fit_function.cpp.
	auto result = conn.Query("SELECT sample_id::VARCHAR, feature_id::VARCHAR, value::DOUBLE FROM " + q);
	if (result->HasError()) {
		throw InvalidInputException("shap: Data relation '%s' does not match the required COO triplet schema.\n"
		                            "  Expected columns : sample_id, feature_id, value\n"
		                            "  Engine error     : %s\n\n"
		                            "Remedy:\n"
		                            "  Wrap it in an aliased view before explaining:\n"
		                            "    CREATE VIEW my_counts AS\n"
		                            "      SELECT your_sample_col  AS sample_id,\n"
		                            "             your_feature_col AS feature_id,\n"
		                            "             your_count_col   AS value\n"
		                            "      FROM %s;",
		                            bind.data_relation, result->GetError(), bind.data_relation);
	}
	while (auto chunk = result->Fetch()) {
		const idx_t n = chunk->size();
		UnifiedVectorFormat sf, ff, vf;
		chunk->data[0].ToUnifiedFormat(n, sf);
		chunk->data[1].ToUnifiedFormat(n, ff);
		chunk->data[2].ToUnifiedFormat(n, vf);
		const auto *samples = UnifiedVectorFormat::GetData<string_t>(sf);
		const auto *features = UnifiedVectorFormat::GetData<string_t>(ff);
		const auto *values = UnifiedVectorFormat::GetData<double>(vf);
		for (idx_t row = 0; row < n; row++) {
			const auto si = sf.sel->get_index(row);
			const auto fi = ff.sel->get_index(row);
			const auto vi = vf.sel->get_index(row);
			if (!sf.validity.RowIsValid(si) || !ff.validity.RowIsValid(fi) || !vf.validity.RowIsValid(vi)) {
				throw InvalidInputException(
				    "shap: NULL in data relation '%s' (sample_id/feature_id/value must all be non-NULL)",
				    bind.data_relation);
			}
			// favor 8 byte alias addr over 16 byte struct value copy onto the stack
			// with string_t& s reads from the chunk's buffer, not a copy, so the builder copies it onto the heap for
			// later use double is 8 bytes so copy is cheap and no pointer indirection is needed
			const string_t &s = samples[si];
			const string_t &f = features[fi];
			builder.Append(std::string_view(s.GetData(), s.GetSize()), std::string_view(f.GetData(), f.GetSize()),
			               values[vi]);
		}
	}
}

//! Take ownership of one array sc wrote inside sc_shap_result_t.
//!
//! The C Data Interface permits moving these structs: copy the bytes, then mark
//! the source released so nothing frees it twice. Done immediately after the
//! call, before the status is checked, so an error path still cleans up.
void TakeArray(ArrowArray &src, ArrowSchema &src_schema, miint::OwnedArrowArray &dst) {
	*dst.array() = src;
	*dst.schema() = src_schema;
	src.release = nullptr;
	src_schema.release = nullptr;
}

//! Choose which features of one (sample, output) attribution row to emit, and
//! in what order.
//!
//! Without top_k: every feature. With it: the strongest positive and strongest
//! negative attributions, ceil(k/2) and floor(k/2), the other side filling in
//! when one runs short. Exact zeros push neither way and are never picked.
//!
//! Either way the result is in waterfall order -- shap_value descending -- with
//! ties on column index, so both the choice and the order are deterministic.
void SelectFeatures(const double *row, size_t n_features, int64_t top_k, std::vector<uint32_t> &chosen,
                    std::vector<uint32_t> &pos, std::vector<uint32_t> &neg) {
	const auto descending = [row](uint32_t a, uint32_t b) {
		return row[a] != row[b] ? row[a] > row[b] : a < b;
	};
	chosen.clear();
	if (top_k <= 0) {
		chosen.resize(n_features);
		for (size_t f = 0; f < n_features; f++) {
			chosen[f] = static_cast<uint32_t>(f);
		}
		std::sort(chosen.begin(), chosen.end(), descending);
		return;
	}
	pos.clear();
	neg.clear();
	for (size_t f = 0; f < n_features; f++) {
		if (row[f] > 0.0) {
			pos.push_back(static_cast<uint32_t>(f));
		} else if (row[f] < 0.0) {
			neg.push_back(static_cast<uint32_t>(f));
		}
	}
	const auto k = static_cast<size_t>(top_k);
	size_t take_pos = std::min((k + 1) / 2, pos.size());
	size_t take_neg = std::min(k / 2, neg.size());
	size_t spare = k - take_pos - take_neg;
	const size_t more_pos = std::min(spare, pos.size() - take_pos);
	take_pos += more_pos;
	spare -= more_pos;
	take_neg += std::min(spare, neg.size() - take_neg);

	const auto pos_end = pos.begin() + static_cast<std::ptrdiff_t>(take_pos);
	const auto neg_end = neg.begin() + static_cast<std::ptrdiff_t>(take_neg);
	std::partial_sort(pos.begin(), pos_end, pos.end(), descending);
	std::partial_sort(neg.begin(), neg_end, neg.end(),
	                  [row](uint32_t a, uint32_t b) { return row[a] != row[b] ? row[a] < row[b] : a < b; });
	chosen.insert(chosen.end(), pos.begin(), pos_end);
	chosen.insert(chosen.end(), neg.begin(), neg_end);
	std::sort(chosen.begin(), chosen.end(), descending);
}

//! Everything every batch shares: context, model, the scanned input and each
//! sample's predicted class. Runs once, before the first row.
void LoadInput(ClientContext &context, const ScShapData &bind, ScShapGlobalState &gstate) {
	auto conn = MakeReadOnlyHelperConnection(context);

	sc_config_t config {};
	config.n_threads = bind.n_threads;
	if (auto st = sc_context_new(&config, &gstate.ctx.ptr); st != SC_OK) {
		miint::ThrowSc("shap", nullptr, st);
	}
	miint::LoadModelFromRelation(conn, bind.model_relation, bind.model_name, "shap", gstate.ctx.ptr, gstate.model);

	miint::OwnedArrowArray vocab;
	if (auto st = sc_model_feature_ids(gstate.model.ptr, vocab.array(), vocab.schema()); st != SC_OK) {
		miint::ThrowSc("sc_model_feature_ids", gstate.ctx.ptr, st);
	}
	gstate.feature_ids = vocab.ReadUtf8("sc_model_feature_ids");

	miint::ScCooBuilder builder;
	builder.SetFeatureVocabulary(gstate.feature_ids);
	ScanForShap(conn, bind, builder);

	const auto dropped = builder.DroppedCells();
	gstate.table = builder.Finalize();
	if (!gstate.table) {
		throw InvalidInputException("shap: data relation '%s' produced no samples", bind.data_relation);
	}
	if (dropped > 0 && gstate.table->NumNonZeros() == 0) {
		throw InvalidInputException("shap: none of the %llu cells in '%s' use a feature this model was trained on; "
		                            "the data and the model do not share a feature vocabulary",
		                            (uint64_t)dropped, bind.data_relation);
	}
	if (dropped > 0) {
		// A sample left with no known features still gets a full, additive
		// explanation -- of an all-zero row. Nothing in the numbers says so.
		size_t empty_samples = 0;
		for (auto c : gstate.table->SampleCoverage()) {
			if (c == 0.0) {
				empty_samples++;
			}
		}
		miint::EmitWarning(context,
		                   "shap: dropped %llu cell(s) from '%s' whose feature the model was not trained on%s. "
		                   "See the sample_coverage column.",
		                   (uint64_t)dropped, bind.data_relation.c_str(),
		                   empty_samples > 0 ? (" -- " + std::to_string(empty_samples) +
		                                        " sample(s) retained no features at all; their attributions explain "
		                                        "an all-zero row")
		                                           .c_str()
		                                     : "");
	}
	gstate.n_features = gstate.feature_ids.size();
	const size_t n_samples = gstate.table->SampleIds().size();

	// A classifier's class count sizes the batches, and the predicted class
	// filters the output. Both come from one proba call over every sample: its
	// output is only samples x classes, and one call spares each batch a second
	// pass over the trees. argmax keeps the first maximum, matching np.argmax --
	// which is what sklearn's and sc's hard prediction is.
	gstate.predicted.assign(n_samples, 0);
	gstate.n_outputs = 1;
	if (bind.classification) {
		miint::OwnedArrowArray proba, classes;
		if (auto st = sc_predict_proba(gstate.ctx.ptr, gstate.model.ptr, gstate.table->get(), proba.array(),
		                               proba.schema(), classes.array(), classes.schema());
		    st != SC_OK) {
			miint::ThrowSc("sc_predict_proba", gstate.ctx.ptr, st);
		}
		gstate.classes = classes.ReadUtf8("sc_predict_proba classes");
		gstate.n_outputs = gstate.classes.size();
		int64_t width = 0;
		const auto p = proba.ReadFixedSizeListFloat64("sc_predict_proba", width);
		for (size_t s = 0; s < n_samples; s++) {
			uint32_t best = 0;
			for (size_t c = 1; c < gstate.n_outputs; c++) {
				if (p[s * gstate.n_outputs + c] > p[s * gstate.n_outputs + best]) {
					best = static_cast<uint32_t>(c);
				}
			}
			gstate.predicted[s] = best;
		}
	}

	// Batches bound memory: sc computes every class of every sample in a batch as
	// a dense double, whatever top_k or predicted_class_only later keeps. Each
	// sample is explained on its own, so the grouping changes no value. Long
	// double so the product cannot overflow on its way to the comparison.
	const long double per_sample =
	    static_cast<long double>(gstate.n_outputs) * static_cast<long double>(gstate.n_features);
	if (bind.batch_size > 0) {
		gstate.batch_size = static_cast<size_t>(bind.batch_size);
	} else if (per_sample > static_cast<long double>(bind.max_attributions)) {
		// A sample cannot be split, so it still runs -- alone.
		gstate.batch_size = 1;
		miint::EmitWarning(context,
		                   "shap: one sample needs %llu attributions (%llu %s x %llu features), above "
		                   "max_attributions %lld; explaining one sample at a time, which needs about %.1f MB at "
		                   "peak.",
		                   (uint64_t)per_sample, (uint64_t)gstate.n_outputs,
		                   gstate.n_outputs == 1 ? "output" : "classes", (uint64_t)gstate.n_features,
		                   (int64_t)bind.max_attributions, static_cast<double>(per_sample * 8 / 1e6));
	} else {
		gstate.batch_size = static_cast<size_t>(static_cast<long double>(bind.max_attributions) / per_sample);
	}
	gstate.batcher = std::make_unique<miint::ScCooBatcher>(*gstate.table);
}

//! Make samples [first, first + batch_size) the current batch.
//!
//! The previous batch's export is released BEFORE sc allocates the next, so
//! only one batch's attributions ever exist at once.
void LoadBatch(const ScShapData &bind, ScShapGlobalState &gstate, size_t first) {
	gstate.values = nullptr;
	gstate.shap_array = std::make_unique<miint::OwnedArrowArray>();
	const size_t n_samples = gstate.table->SampleIds().size();
	const size_t count = std::min(gstate.batch_size, n_samples - first);
	auto batch = gstate.batcher->Batch(first, count);

	sc_shap_result_t res {};
	const auto st = sc_shap(gstate.ctx.ptr, gstate.model.ptr, batch->get(), &res);
	miint::OwnedArrowArray base_values, shap_classes;
	TakeArray(res.shap_values, res.shap_values_schema, *gstate.shap_array);
	TakeArray(res.base_values, res.base_values_schema, base_values);
	TakeArray(res.classes, res.classes_schema, shap_classes);
	if (st != SC_OK) {
		miint::ThrowSc("sc_shap", gstate.ctx.ptr, st);
	}

	// A model trained on a feature-selected subset explains fewer columns than
	// its vocabulary. The fit functions never select features, so this only
	// fires on a model table assembled some other way.
	if (static_cast<size_t>(res.n_features) != gstate.n_features ||
	    static_cast<size_t>(res.n_outputs) != gstate.n_outputs) {
		throw InvalidInputException(
		    "shap: the model in '%s' explains %lld features x %lld outputs, but its vocabulary has %llu features "
		    "and %llu outputs; feature-selected models are not supported by shap",
		    bind.model_relation, (int64_t)res.n_features, (int64_t)res.n_outputs, (uint64_t)gstate.n_features,
		    (uint64_t)gstate.n_outputs);
	}
	gstate.base_values = base_values.ReadFloat64("sc_shap base_values");
	int64_t width = 0;
	// Read in place: this buffer is the one copy of the batch's attributions, and
	// gstate keeps it alive until the batch's last row is emitted.
	gstate.values = gstate.shap_array->FixedSizeListFloat64Data("sc_shap", width);
	const auto rows = static_cast<size_t>(gstate.shap_array->array()->length);
	if (static_cast<size_t>(width) != gstate.n_outputs * gstate.n_features || rows != count ||
	    gstate.base_values.size() != gstate.n_outputs) {
		throw InternalException("sc_shap: attribution array has width %lld and %llu rows for %llu samples x "
		                        "%llu outputs x %llu features",
		                        (int64_t)width, (uint64_t)rows, (uint64_t)count, (uint64_t)gstate.n_outputs,
		                        (uint64_t)gstate.n_features);
	}
	// The output axis of the attributions must be the same class order proba
	// reported, or every attribution would be labelled with the wrong class.
	if (bind.classification && shap_classes.ReadUtf8("sc_shap classes") != gstate.classes) {
		throw InternalException("sc_shap: class order differs between sc_shap and sc_predict_proba");
	}
	gstate.batch_first = first;
	gstate.batch_count = count;
}

void ScShapExecute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<ScShapGlobalState>();
	const auto &bind = input.bind_data->Cast<ScShapData>();
	if (!gstate.loaded) {
		gstate.loaded = true;
		LoadInput(context, bind, gstate);
	}

	const auto &sample_ids = gstate.table->SampleIds();
	const auto &coverage = gstate.table->SampleCoverage();
	const size_t n_samples = sample_ids.size();
	const bool filter_output = bind.classification && bind.predicted_class_only;
	auto advance = [&gstate]() {
		gstate.have_row = false;
		if (++gstate.cur_output == gstate.n_outputs) {
			gstate.cur_output = 0;
			gstate.cur_sample++;
		}
	};

	idx_t n = 0;
	while (n < STANDARD_VECTOR_SIZE) {
		if (!gstate.have_row) {
			while (gstate.cur_sample < n_samples && filter_output &&
			       gstate.cur_output != gstate.predicted[gstate.cur_sample]) {
				advance();
			}
			if (gstate.cur_sample >= n_samples) {
				break;
			}
			if (gstate.cur_sample >= gstate.batch_first + gstate.batch_count) {
				LoadBatch(bind, gstate, gstate.cur_sample);
			}
			// Output-major within a sample: [o0f0, o0f1, ..., o1f0, ...].
			const double *row =
			    gstate.values +
			    ((gstate.cur_sample - gstate.batch_first) * gstate.n_outputs + gstate.cur_output) * gstate.n_features;
			SelectFeatures(row, gstate.n_features, bind.top_k, gstate.chosen, gstate.scratch_pos, gstate.scratch_neg);
			gstate.next_chosen = 0;
			gstate.have_row = true;
		}
		if (gstate.next_chosen == gstate.chosen.size()) {
			advance();
			continue;
		}
		const auto f = gstate.chosen[gstate.next_chosen++];
		const auto base =
		    ((gstate.cur_sample - gstate.batch_first) * gstate.n_outputs + gstate.cur_output) * gstate.n_features;
		output.SetValue(0, n, Value(sample_ids[gstate.cur_sample]).DefaultCastAs(bind.id_types.sample_id));
		output.SetValue(1, n,
		                bind.classification ? Value(gstate.classes[gstate.cur_output]) : Value(LogicalType::VARCHAR));
		output.SetValue(2, n, Value(gstate.feature_ids[f]).DefaultCastAs(bind.id_types.feature_id));
		output.SetValue(3, n, Value::DOUBLE(gstate.values[base + f]));
		output.SetValue(4, n, Value::DOUBLE(gstate.base_values[gstate.cur_output]));
		output.SetValue(5, n, Value::DOUBLE(coverage[gstate.cur_sample]));
		n++;
	}
	output.SetCardinality(n);
}

} // namespace

void ScShapFunction::Register(ExtensionLoader &loader) {
	TableFunction fn("shap", {LogicalType::VARCHAR, LogicalType::VARCHAR}, ScShapExecute, ScShapBind, ScShapInitGlobal);
	fn.named_parameters["name"] = LogicalType::VARCHAR;
	fn.named_parameters["predicted_class_only"] = LogicalType::BOOLEAN;
	fn.named_parameters["top_k"] = LogicalType::BIGINT;
	fn.named_parameters["max_attributions"] = LogicalType::BIGINT;
	fn.named_parameters["batch_size"] = LogicalType::BIGINT;
	fn.named_parameters["n_threads"] = LogicalType::INTEGER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
