#include "sc_feature_importances_function.hpp"

#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "sc_common.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/main/connection.hpp"

#include <string>
#include <vector>

namespace duckdb {

namespace {

struct ScImportancesData : public TableFunctionData {
	string model_relation;
	//! Empty selects the whole relation, which must then be one row.
	string model_name;
	//! The type these ids had when the model was fit; see ScModelTypes.
	LogicalType feature_id_type = LogicalType::VARCHAR;
};

struct ScImportancesGlobalState : public GlobalTableFunctionState {
	vector<string> feature_ids;
	vector<double> importances;
	idx_t emitted = 0;
	bool loaded = false;
};

unique_ptr<FunctionData> ScImportancesBind(ClientContext &context, TableFunctionBindInput &input,
                                           vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<ScImportancesData>();
	data->model_relation = input.inputs[0].GetValue<string>();
	if (data->model_relation.empty()) {
		throw InvalidInputException("sc_feature_importances: model relation name must not be empty");
	}
	for (auto &kv : input.named_parameters) {
		if (!kv.second.IsNull() && StringUtil::CIEquals(kv.first, "name")) {
			data->model_name = kv.second.GetValue<string>();
		}
	}
	{
		auto conn = MakeReadOnlyHelperConnection(context);
		data->feature_id_type =
		    miint::ReadModelTypes(conn, context, data->model_relation, data->model_name, "sc_feature_importances")
		        .feature_id_type;
	}
	names = {"feature_id", "importance"};
	return_types = {data->feature_id_type, LogicalType::DOUBLE};
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ScImportancesInitGlobal(ClientContext &, TableFunctionInitInput &) {
	return make_uniq<ScImportancesGlobalState>();
}

void ScImportancesExecute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<ScImportancesGlobalState>();
	const auto &bind = input.bind_data->Cast<ScImportancesData>();

	if (!gstate.loaded) {
		gstate.loaded = true;
		auto conn = MakeReadOnlyHelperConnection(context);

		miint::ScContext ctx;
		sc_config_t config {};
		if (auto st = sc_context_new(&config, &ctx.ptr); st != SC_OK) {
			miint::ThrowSc("sc_feature_importances", nullptr, st);
		}
		miint::ScModel model;
		miint::LoadModelFromRelation(conn, bind.model_relation, bind.model_name, "sc_feature_importances", ctx.ptr,
		                             model);

		// Two exports, both owned by us from here. The RAII wrapper releases
		// them even if the second call throws.
		miint::OwnedArrowArray ids;
		if (auto st = sc_model_feature_ids(model.ptr, ids.array(), ids.schema()); st != SC_OK) {
			miint::ThrowSc("sc_model_feature_ids", ctx.ptr, st);
		}
		miint::OwnedArrowArray weights;
		if (auto st = sc_feature_importances(model.ptr, weights.array(), weights.schema()); st != SC_OK) {
			miint::ThrowSc("sc_feature_importances", ctx.ptr, st);
		}

		gstate.feature_ids = ids.ReadUtf8("sc_model_feature_ids");
		gstate.importances = weights.ReadFloat64("sc_feature_importances");
		// Both are indexed by model column, so a length mismatch means the two
		// have drifted apart inside sc -- not something a caller can cause, but
		// silently zipping mismatched arrays would mislabel every row.
		if (gstate.feature_ids.size() != gstate.importances.size()) {
			throw InternalException("sc_feature_importances: %llu feature ids but %llu importances",
			                        (unsigned long long)gstate.feature_ids.size(),
			                        (unsigned long long)gstate.importances.size());
		}
	}

	const idx_t remaining = gstate.feature_ids.size() - gstate.emitted;
	const idx_t n = remaining < STANDARD_VECTOR_SIZE ? remaining : STANDARD_VECTOR_SIZE;
	output.SetCardinality(n);
	for (idx_t i = 0; i < n; i++) {
		EmitIdCell(output.data[0], i, gstate.feature_ids[gstate.emitted + i], bind.feature_id_type);
		output.SetValue(1, i, Value::DOUBLE(gstate.importances[gstate.emitted + i]));
	}
	gstate.emitted += n;
}

} // namespace

void ScFeatureImportancesFunction::Register(ExtensionLoader &loader) {
	TableFunction fn("sc_feature_importances", {LogicalType::VARCHAR}, ScImportancesExecute, ScImportancesBind,
	                 ScImportancesInitGlobal);
	fn.named_parameters["name"] = LogicalType::VARCHAR;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
