#include "sc_model_features_function.hpp"

#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "sc_common.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/main/connection.hpp"

#include <string>
#include <vector>

namespace duckdb {

namespace {

struct ScModelFeaturesData : public TableFunctionData {
	string model_relation;
	//! Empty selects the whole relation, which must then be one row.
	string model_name;
	//! The type these ids had when the model was fit; see ScModelTypes.
	LogicalType feature_id_type = LogicalType::VARCHAR;
};

struct ScModelFeaturesGlobalState : public GlobalTableFunctionState {
	std::vector<std::string> feature_ids;
	idx_t emitted = 0;
	bool loaded = false;
};

unique_ptr<FunctionData> ScModelFeaturesBind(ClientContext &context, TableFunctionBindInput &input,
                                             vector<LogicalType> &return_types, vector<string> &names) {
	// calls ScModelFeaturesData constructor and allocates on the heap, returning a unique_ptr to it.  The unique_ptr
	// will automatically free the memory when it goes out of scope, so we don't have to worry about memory leaks. do
	// make_uniq for RAII.  It makes a smart pointer. this is C++ lifecycle management for heap data
	auto data = make_uniq<ScModelFeaturesData>();
	data->model_relation = input.inputs[0].GetValue<string>();
	if (data->model_relation.empty()) {
		// if we didn't have a smart ptr then this path would cause a memory leak
		throw InvalidInputException("sc_model_features: model relation name must not be empty");
	}
	for (auto &kv : input.named_parameters) {
		if (!kv.second.IsNull() && StringUtil::CIEquals(kv.first, "name")) {
			data->model_name = kv.second.GetValue<string>();
		}
	}
	// set the output column names and types.  The first column is the feature id, which is a string.  The second column
	// is the column index, which is an integer.
	{
		auto conn = MakeReadOnlyHelperConnection(context);
		data->feature_id_type =
		    miint::ReadModelTypes(conn, context, data->model_relation, data->model_name, "sc_model_features")
		        .feature_id_type;
	}
	names = {"feature_id", "column_index"};
	return_types = {data->feature_id_type, LogicalType::BIGINT};
	// return metadata about the table function to the engine.  The engine will use this metadata to create the output
	// table.  The engine will call ScModelFeaturesExecute to fill in the output table.
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ScModelFeaturesInitGlobal(ClientContext &, TableFunctionInitInput &) {
	return make_uniq<ScModelFeaturesGlobalState>();
}

void ScModelFeaturesExecute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<ScModelFeaturesGlobalState>();
	const auto &bind = input.bind_data->Cast<ScModelFeaturesData>();

	if (!gstate.loaded) {
		gstate.loaded = true;
		auto conn = MakeReadOnlyHelperConnection(context);

		miint::ScContext ctx;
		sc_config_t config {};
		if (auto st = sc_context_new(&config, &ctx.ptr); st != SC_OK) {
			miint::ThrowSc("sc_model_features", nullptr, st);
		}
		miint::ScModel model;
		miint::LoadModelFromRelation(conn, bind.model_relation, bind.model_name, "sc_model_features", ctx.ptr, model);

		miint::OwnedArrowArray ids;
		if (auto st = sc_model_feature_ids(model.ptr, ids.array(), ids.schema()); st != SC_OK) {
			miint::ThrowSc("sc_model_feature_ids", ctx.ptr, st);
		}
		gstate.feature_ids = ids.ReadUtf8("sc_model_feature_ids");
	}

	const idx_t remaining = gstate.feature_ids.size() - gstate.emitted;
	const idx_t n = remaining < STANDARD_VECTOR_SIZE ? remaining : STANDARD_VECTOR_SIZE;
	output.SetCardinality(n);
	for (idx_t i = 0; i < n; i++) {
		const auto at = gstate.emitted + i;
		EmitIdCell(output.data[0], i, gstate.feature_ids[at], bind.feature_id_type);
		// Position in the model's matrix, which is what fit and predict must
		// agree on.
		output.SetValue(1, i, Value::BIGINT(static_cast<int64_t>(at)));
	}
	gstate.emitted += n;
}

} // namespace

void ScModelFeaturesFunction::Register(ExtensionLoader &loader) {
	TableFunction fn("sc_model_features", {LogicalType::VARCHAR}, ScModelFeaturesExecute, ScModelFeaturesBind,
	                 ScModelFeaturesInitGlobal);
	fn.named_parameters["name"] = LogicalType::VARCHAR;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
