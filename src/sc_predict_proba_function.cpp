#include "sc_predict_proba_function.hpp"

#include "catalog_utils.hpp"
#include "miint_log.hpp"
#include "id_column_utils.hpp"
#include "sc_common.hpp"
#include "sc_rf_common.hpp"
#include "coo_builder.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <string>
#include <string_view>
#include <vector>

namespace duckdb {

namespace {

struct ScProbaData : public TableFunctionData {
	string data_relation;
	string model_relation;
	string model_name;
	int32_t n_threads = 0;
	//! Mirrors the data relation; see ScTrainingInput.
	LogicalType sample_id_type = LogicalType::VARCHAR;
	//! From the model: what its class labels were.
	LogicalType target_type = LogicalType::VARCHAR;
};

struct ScProbaGlobalState : public GlobalTableFunctionState {
	std::vector<std::string> sample_ids;
	std::vector<std::string> classes;
	std::vector<double> coverage;
	//! Row-major, n_samples * n_classes.
	std::vector<double> proba;
	//! One output row per (sample, class) pair.
	idx_t emitted = 0;
	bool loaded = false;
};

unique_ptr<FunctionData> ScProbaBind(ClientContext &context, TableFunctionBindInput &input,
                                     vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<ScProbaData>();
	data->data_relation = input.inputs[0].GetValue<string>();
	data->model_relation = input.inputs[1].GetValue<string>();
	if (data->data_relation.empty() || data->model_relation.empty()) {
		throw InvalidInputException("sc_predict_proba: data and model relation names must not be empty");
	}
	for (auto &kv : input.named_parameters) {
		if (kv.second.IsNull()) {
			continue;
		}
		if (StringUtil::CIEquals(kv.first, "n_threads")) {
			data->n_threads = kv.second.GetValue<int32_t>();
		} else if (StringUtil::CIEquals(kv.first, "name")) {
			data->model_name = kv.second.GetValue<string>();
		}
	}

	// Reject a regressor at bind, before any scanning. sc would refuse too, but
	// only after the whole COO had been marshalled.
	{
		auto conn = MakeReadOnlyHelperConnection(context);
		const auto task = miint::ReadModelTask(conn, data->model_relation, data->model_name, "sc_predict_proba");
		if (task != "classification") {
			throw InvalidInputException(
			    "sc_predict_proba: The model in '%s' is a regressor, which has no classes to give "
			    "probabilities over.\n\n"
			    "Remedy:\n"
			    "  Use sc_predict for a regressor -- it returns the predicted value directly:\n"
			    "    SELECT * FROM sc_predict('%s', '%s');",
			    data->model_relation, data->data_relation, data->model_relation);
		}
		const auto id_types = sc_rf::DetectCooIdTypes(conn, data->data_relation, "sc_predict_proba");
		data->sample_id_type = id_types.sample_id_type;
		data->target_type =
		    miint::ReadModelTypes(conn, context, data->model_relation, data->model_name, "sc_predict_proba")
		        .target_type;
	}

	names = {"sample_id", "class", "probability", "sample_coverage"};
	return_types = {data->sample_id_type, data->target_type, LogicalType::DOUBLE, LogicalType::DOUBLE};
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ScProbaInitGlobal(ClientContext &, TableFunctionInitInput &) {
	return make_uniq<ScProbaGlobalState>();
}

void ScanForProba(Connection &conn, const ScProbaData &bind, miint::CooBuilder &builder) {
	const auto q = KeywordHelper::WriteOptionallyQuoted(bind.data_relation);
	// The casts guarantee the physical layout the buffer reads below assume;
	// see the note in sc_fit_function.cpp.
	auto result = conn.Query("SELECT sample_id::VARCHAR, feature_id::VARCHAR, value::DOUBLE FROM " + q);
	if (result->HasError()) {
		sc_rf::ThrowNotCooTriplet(bind.data_relation, result->GetError(), "sc_predict_proba");
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
				throw InvalidInputException("sc_predict_proba: NULL in data relation '%s' "
				                            "(sample_id/feature_id/value must all be non-NULL)",
				                            bind.data_relation);
			}
			const auto &s = samples[si];
			const auto &f = features[fi];
			builder.Append(std::string_view(s.GetData(), s.GetSize()), std::string_view(f.GetData(), f.GetSize()),
			               values[vi]);
		}
	}
}

void ScProbaExecute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<ScProbaGlobalState>();
	const auto &bind = input.bind_data->Cast<ScProbaData>();

	if (!gstate.loaded) {
		gstate.loaded = true;
		auto conn = MakeReadOnlyHelperConnection(context);

		sc_config_t config {};
		config.n_threads = bind.n_threads;
		miint::ScContext ctx;
		if (auto st = sc_context_new(&config, &ctx.ptr); st != SC_OK) {
			miint::ThrowSc("sc_predict_proba", nullptr, st);
		}
		miint::ScModel model;
		miint::LoadModelFromRelation(conn, bind.model_relation, bind.model_name, "sc_predict_proba", ctx.ptr, model);

		miint::OwnedArrowArray vocab;
		if (auto st = sc_model_feature_ids(model.ptr, vocab.array(), vocab.schema()); st != SC_OK) {
			miint::ThrowSc("sc_model_feature_ids", ctx.ptr, st);
		}
		miint::CooBuilder builder;
		builder.SetFeatureVocabulary(vocab.ReadUtf8("sc_model_feature_ids"));
		ScanForProba(conn, bind, builder);

		const auto dropped = builder.DroppedCells();
		// Finalize() resets the builder, so take the diagnostic sample first.
		const auto dropped_examples = builder.DroppedExamples();
		auto table = builder.Finalize();
		if (!table) {
			throw InvalidInputException("sc_predict_proba: data relation '%s' produced no samples", bind.data_relation);
		}
		if (dropped > 0 && table->NumNonZeros() == 0) {
			throw InvalidInputException(
			    "sc_predict_proba: none of the %llu cells in '%s' use a feature this model was trained on; "
			    "the data and the model do not share a feature vocabulary%s",
			    (unsigned long long)dropped, bind.data_relation,
			    miint::VocabularyMismatchHint(dropped_examples, table->FeatureIds(), bind.data_relation));
		}
		if (dropped > 0) {
			miint::EmitWarning(context,
			                   "sc_predict_proba: dropped %llu cell(s) from '%s' whose feature the model was not "
			                   "trained on. See the sample_coverage column.",
			                   (unsigned long long)dropped, bind.data_relation.c_str());
		}
		gstate.sample_ids = table->SampleIds();
		gstate.coverage = table->SampleCoverage();

		const auto sc_table = miint::AsScTable(*table);
		miint::OwnedArrowArray proba, classes;
		if (auto st = sc_predict_proba(ctx.ptr, model.ptr, &sc_table, proba.array(), proba.schema(), classes.array(),
		                               classes.schema());
		    st != SC_OK) {
			miint::ThrowSc("sc_predict_proba", ctx.ptr, st);
		}
		gstate.classes = classes.ReadUtf8("sc_predict_proba classes");

		int64_t width = 0;
		gstate.proba = proba.ReadFixedSizeListFloat64("sc_predict_proba", width);
		// The list width IS the class count -- one probability per class, in the
		// order the classes array reports.
		if (static_cast<size_t>(width) != gstate.classes.size()) {
			throw InternalException("sc_predict_proba: %lld probabilities per sample but %llu classes",
			                        (long long)width, (unsigned long long)gstate.classes.size());
		}
		if (gstate.proba.size() != gstate.sample_ids.size() * gstate.classes.size()) {
			throw InternalException("sc_predict_proba: %llu probabilities for %llu samples x %llu classes",
			                        (unsigned long long)gstate.proba.size(),
			                        (unsigned long long)gstate.sample_ids.size(),
			                        (unsigned long long)gstate.classes.size());
		}
	}

	// One row per (sample, class): index i is sample i / n_classes, class
	// i % n_classes -- the same row-major order sc flattened them in.
	const idx_t n_classes = gstate.classes.size();
	const idx_t total = gstate.sample_ids.size() * n_classes;
	const idx_t remaining = total - gstate.emitted;
	const idx_t n = remaining < STANDARD_VECTOR_SIZE ? remaining : STANDARD_VECTOR_SIZE;
	output.SetCardinality(n);
	for (idx_t i = 0; i < n; i++) {
		const auto at = gstate.emitted + i;
		const auto sample = at / n_classes;
		const auto klass = at % n_classes;
		EmitIdCell(output.data[0], i, gstate.sample_ids[sample], bind.sample_id_type);
		output.SetValue(1, i, Value(gstate.classes[klass]).DefaultCastAs(bind.target_type));
		output.SetValue(2, i, Value::DOUBLE(gstate.proba[at]));
		output.SetValue(3, i, Value::DOUBLE(gstate.coverage[sample]));
	}
	gstate.emitted += n;
}

} // namespace

void ScPredictProbaFunction::Register(ExtensionLoader &loader) {
	TableFunction fn("sc_predict_proba", {LogicalType::VARCHAR, LogicalType::VARCHAR}, ScProbaExecute, ScProbaBind,
	                 ScProbaInitGlobal);
	fn.named_parameters["n_threads"] = LogicalType::INTEGER;
	fn.named_parameters["name"] = LogicalType::VARCHAR;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
