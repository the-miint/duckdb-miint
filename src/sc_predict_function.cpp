#include "sc_predict_function.hpp"

#include "catalog_utils.hpp"
#include "sc_common.hpp"
#include "sc_coo_builder.hpp"
#include "miint_log.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <string>
#include <string_view>
#include <vector>

namespace duckdb {

namespace {

struct ScPredictData : public TableFunctionData {
	string data_relation;
	string model_relation;
	//! Empty selects the whole relation, which must then be one row.
	string model_name;
	bool classification = false;
	int32_t n_threads = 0;
	miint::ScIdTypes id_types;
};

struct ScPredictGlobalState : public GlobalTableFunctionState {
	// std::vector, not duckdb::vector -- these are handed over from
	// ScCooBuilder / OwnedArrowArray, which are outside duckdb's namespace.
	std::vector<std::string> sample_ids;
	std::vector<std::string> labels; // classification
	std::vector<double> values;      // regression
	std::vector<double> coverage;
	idx_t emitted = 0;
	bool loaded = false;
};

unique_ptr<FunctionData> ScPredictBind(ClientContext &context, TableFunctionBindInput &input,
                                       vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<ScPredictData>();
	data->data_relation = input.inputs[0].GetValue<string>();
	data->model_relation = input.inputs[1].GetValue<string>();
	if (data->data_relation.empty() || data->model_relation.empty()) {
		throw InvalidInputException("predict: data and model relation names must not be empty");
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
	// The model knows its own task; ask it rather than making the caller assert
	// one. This is why the prediction column can be typed at all.
	{
		auto conn = MakeReadOnlyHelperConnection(context);
		data->id_types = miint::ProbeIdTypes(conn, data->data_relation, "predict");
		data->classification =
		    miint::ReadModelTask(conn, data->model_relation, data->model_name, "predict") == "classification";
		miint::RequireModelFeatureIdType(
		    miint::ReadModelFeatureIdType(conn, data->model_relation, data->model_name, "predict"),
		    data->id_types.feature_id, "predict", data->model_relation, data->model_name, data->data_relation);
	}

	names = {"sample_id", "prediction", "sample_coverage"};
	// A classifier predicts one of its training labels; a regressor a number.
	// sample_coverage is matched/observed for that sample -- see ScCooTable.
	return_types = {data->id_types.sample_id, data->classification ? LogicalType::VARCHAR : LogicalType::DOUBLE,
	                LogicalType::DOUBLE};
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ScPredictInitGlobal(ClientContext &, TableFunctionInitInput &) {
	return make_uniq<ScPredictGlobalState>();
}

//! Scan the prediction data into `builder`, which already carries the model's
//! vocabulary. Reads through UnifiedVectorFormat so a cell costs no allocation.
void ScanForPrediction(Connection &conn, const ScPredictData &bind, miint::ScCooBuilder &builder) {
	const auto q = KeywordHelper::WriteOptionallyQuoted(bind.data_relation);
	// The casts are load-bearing, not cosmetic. Reading a vector's buffer
	// directly assumes its physical type, and DuckDB infers `42.0` as
	// DECIMAL(3,1) (physical INT16), an INTEGER count column as INT32, and
	// woltka's feature ids as BIGINT or UUID. Vector::GetValue used to convert
	// on the way out; a raw buffer read cannot. Casting in SQL makes DuckDB do
	// the conversion and guarantees the layout this loop reads.
	auto result = conn.Query("SELECT sample_id::VARCHAR, feature_id::VARCHAR, value::DOUBLE FROM " + q);
	if (result->HasError()) {
		throw InvalidInputException("predict: Data relation '%s' does not match the required COO triplet schema.\n"
		                            "  Expected columns : sample_id, feature_id, value\n"
		                            "  Engine error     : %s\n\n"
		                            "Remedy:\n"
		                            "  Prediction data must be shaped exactly like the data the model was fit on.\n"
		                            "  If your relation uses different column names, wrap it in an aliased view:\n"
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
				    "predict: NULL in data relation '%s' (sample_id/feature_id/value must all be non-NULL)",
				    bind.data_relation);
			}
			const auto &s = samples[si];
			const auto &f = features[fi];
			builder.Append(std::string_view(s.GetData(), s.GetSize()), std::string_view(f.GetData(), f.GetSize()),
			               values[vi]);
		}
	}
}

void ScPredictExecute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<ScPredictGlobalState>();
	const auto &bind = input.bind_data->Cast<ScPredictData>();

	if (!gstate.loaded) {
		gstate.loaded = true;
		auto conn = MakeReadOnlyHelperConnection(context);

		sc_config_t config {};
		config.n_threads = bind.n_threads;
		miint::ScContext ctx;
		if (auto st = sc_context_new(&config, &ctx.ptr); st != SC_OK) {
			miint::ThrowSc("predict", nullptr, st);
		}
		miint::ScModel model;
		miint::LoadModelFromRelation(conn, bind.model_relation, bind.model_name, "predict", ctx.ptr, model);

		// The vocabulary comes out of the model, in the model's column order.
		// This is the whole reason sc_model_feature_ids exists.
		miint::OwnedArrowArray vocab_array;
		if (auto st = sc_model_feature_ids(model.ptr, vocab_array.array(), vocab_array.schema()); st != SC_OK) {
			miint::ThrowSc("sc_model_feature_ids", ctx.ptr, st);
		}
		auto vocab = vocab_array.ReadUtf8("sc_model_feature_ids");

		miint::ScCooBuilder builder;
		builder.SetFeatureVocabulary(std::move(vocab));
		ScanForPrediction(conn, bind, builder);

		const auto builder_dropped = builder.DroppedCells();
		auto table = builder.Finalize();
		if (!table) {
			throw InvalidInputException("predict: data relation '%s' produced no samples", bind.data_relation);
		}
		gstate.sample_ids = table->SampleIds();
		gstate.coverage = table->SampleCoverage();

		// Every cell dropped means the model shares no vocabulary at all with
		// this data -- a different reference database or pipeline, not a
		// borderline case. Predicting from an all-zero matrix would return
		// confident nonsense.
		if (builder_dropped > 0 && table->NumNonZeros() == 0) {
			throw InvalidInputException(
			    "predict: none of the %llu cells in '%s' use a feature this model was trained on; "
			    "the data and the model do not share a feature vocabulary",
			    (uint64_t)builder_dropped, bind.data_relation);
		}
		size_t empty_samples = 0;
		for (auto c : gstate.coverage) {
			if (c == 0.0) {
				empty_samples++;
			}
		}
		if (builder_dropped > 0) {
			miint::EmitWarning(context,
			                   "predict: dropped %llu cell(s) from '%s' whose feature the model was not "
			                   "trained on%s. See the sample_coverage column.",
			                   (uint64_t)builder_dropped, bind.data_relation.c_str(),
			                   empty_samples > 0
			                       ? (" -- " + std::to_string(empty_samples) +
			                          " sample(s) retained no features at all and are predicted from an all-zero row")
			                             .c_str()
			                       : "");
		}

		miint::OwnedArrowArray pred;
		if (auto st = sc_predict(ctx.ptr, model.ptr, table->get(), pred.array(), pred.schema()); st != SC_OK) {
			miint::ThrowSc("sc_predict", ctx.ptr, st);
		}
		// sc's output type follows the MODEL's task, not the function the caller
		// reached for: "u" for a classifier's labels, "g" for a regressor's
		// numbers. Checking it here turns an Arrow format mismatch -- which
		// would surface as "expected 'g', got 'u'" -- into advice.
		// The bound return type came from the relation's `task` column; this is
		// what sc actually produced. A disagreement means the column does not
		// describe the blob beside it -- a hand-assembled model table.
		const char *fmt = pred.Format();
		const bool model_is_classifier = fmt && std::string(fmt) == "u";
		if (model_is_classifier != bind.classification) {
			throw InvalidInputException(
			    "predict: model relation '%s' says task '%s', but the stored model is a %s; the task column "
			    "and the model blob do not match",
			    bind.model_relation, bind.classification ? "classification" : "regression",
			    model_is_classifier ? "classifier" : "regressor");
		}
		if (bind.classification) {
			gstate.labels = pred.ReadUtf8("sc_predict");
		} else {
			gstate.values = pred.ReadFloat64("sc_predict");
		}

		const size_t got = bind.classification ? gstate.labels.size() : gstate.values.size();
		if (got != gstate.sample_ids.size()) {
			throw InternalException("sc_predict: %llu predictions for %llu samples", (uint64_t)got,
			                        (uint64_t)gstate.sample_ids.size());
		}
	}

	const idx_t total = gstate.sample_ids.size();
	const idx_t remaining = total - gstate.emitted;
	const idx_t n = remaining < STANDARD_VECTOR_SIZE ? remaining : STANDARD_VECTOR_SIZE;
	output.SetCardinality(n);
	for (idx_t i = 0; i < n; i++) {
		const auto at = gstate.emitted + i;
		output.SetValue(0, i, Value(gstate.sample_ids[at]).DefaultCastAs(bind.id_types.sample_id));
		if (bind.classification) {
			output.SetValue(1, i, Value(gstate.labels[at]));
		} else {
			output.SetValue(1, i, Value::DOUBLE(gstate.values[at]));
		}
		output.SetValue(2, i, Value::DOUBLE(gstate.coverage[at]));
	}
	gstate.emitted += n;
}

} // namespace

void ScPredictFunction::Register(ExtensionLoader &loader) {
	TableFunction fn("predict", {LogicalType::VARCHAR, LogicalType::VARCHAR}, ScPredictExecute, ScPredictBind,
	                 ScPredictInitGlobal);
	fn.named_parameters["n_threads"] = LogicalType::INTEGER;
	fn.named_parameters["name"] = LogicalType::VARCHAR;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
