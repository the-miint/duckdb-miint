#include "sc_fit_function.hpp"

#include "catalog_utils.hpp"
#include "coo_builder.hpp"
#include "sc_common.hpp"
#include "sc_rf_common.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "sc.h"

namespace duckdb {

namespace {

// Everything fit and cross-validate share -- input scanning, validation,
// forest-parameter parsing -- lives in sc_rf_common.
using namespace sc_rf;

//! The two relations plus the forest parameters, resolved at bind time.
//!
//! Inherits ScTrainingInput for the relation names the shared helpers read.
struct ScFitData : public TableFunctionData, public ScTrainingInput {
	//! Required. A model is selected by name, so an unnamed one is reachable
	//! only by position or by whatever hyperparameters happen to differ.
	string name;
	//! The algorithm, echoed into the output so a row says what produced it.
	string model = "random_forest";
	bool classification = false;
	int32_t n_threads = 0;
	// Defaults are sklearn's RandomForest{Classifier,Regressor} defaults, so an
	// unparameterised call matches what a user would get from scikit-learn.
	// `criterion` and `max_features` differ by task, so they are resolved in
	// Bind once `classification` is known.
	sc_rf_params_t params {};
};

// ---------------------------------------------------------------------------
// Bind / Execute
// ---------------------------------------------------------------------------

unique_ptr<FunctionData> ScFitBind(ClientContext &context, TableFunctionBindInput &input,
                                   vector<LogicalType> &return_types, vector<string> &names, bool classification) {
	auto data = make_uniq<ScFitData>();
	data->classification = classification;
	data->data_relation = input.inputs[0].GetValue<string>();
	data->metadata_relation = input.inputs[1].GetValue<string>();
	if (data->data_relation.empty() || data->metadata_relation.empty()) {
		throw InvalidInputException("sc_fit: data and metadata relation names must not be empty");
	}

	ApplyDefaults(data->params, classification);
	auto &params = data->params;
	for (auto &kv : input.named_parameters) {
		const auto &k = kv.first;
		const auto &v = kv.second;
		if (v.IsNull()) {
			throw InvalidInputException("sc_fit: named parameter '%s' must not be NULL", k);
		}
		if (StringUtil::CIEquals(k, "name")) {
			data->name = v.GetValue<string>();
			if (data->name.empty()) {
				throw InvalidInputException("sc_fit: name must not be empty; omit it to leave the model unnamed");
			}
		} else if (StringUtil::CIEquals(k, "target_column")) {
			data->target_column = v.GetValue<string>();
		} else if (StringUtil::CIEquals(k, "n_threads")) {
			data->n_threads = v.GetValue<int32_t>();
			params.n_threads = data->n_threads;
		} else if (StringUtil::CIEquals(k, "model")) {
			const auto m = v.GetValue<string>();
			if (!StringUtil::CIEquals(m, "random_forest")) {
				throw InvalidInputException("sc_fit: model must be 'random_forest' (got '%s')", m);
			}
			data->model = "random_forest"; // normalised, so the column is stable
		} else if (StringUtil::CIEquals(k, "n_estimators")) {
			params.n_estimators = v.GetValue<int64_t>();
			if (params.n_estimators <= 0) {
				throw InvalidInputException("sc_fit: n_estimators must be > 0 (got %lld)",
				                            (long long)params.n_estimators);
			}
		} else if (StringUtil::CIEquals(k, "random_state")) {
			params.random_state = static_cast<uint64_t>(v.GetValue<int64_t>());
		} else if (StringUtil::CIEquals(k, "max_depth")) {
			// <= 0 is sklearn's None (unbounded), which sc reads the same way.
			params.max_depth = v.GetValue<int64_t>();
		} else if (StringUtil::CIEquals(k, "criterion")) {
			ParseCriterion(v, classification, data->caller, params.criterion);
		} else if (StringUtil::CIEquals(k, "max_features")) {
			ParseMaxFeatures(v, data->caller, params.max_features);
		} else if (StringUtil::CIEquals(k, "min_samples_split")) {
			ParseMinSamples(v, "min_samples_split", 2, data->caller, params.min_samples_split);
		} else if (StringUtil::CIEquals(k, "min_samples_leaf")) {
			ParseMinSamples(v, "min_samples_leaf", 1, data->caller, params.min_samples_leaf);
		} else if (StringUtil::CIEquals(k, "min_weight_fraction_leaf")) {
			params.min_weight_fraction_leaf = v.GetValue<double>();
		} else if (StringUtil::CIEquals(k, "min_impurity_decrease")) {
			params.min_impurity_decrease = v.GetValue<double>();
		} else if (StringUtil::CIEquals(k, "bootstrap")) {
			params.bootstrap = v.GetValue<bool>();
		} else if (StringUtil::CIEquals(k, "max_samples")) {
			ParseMaxSamples(v, data->caller, params.max_samples);
		}
	}
	{
		auto conn = MakeReadOnlyHelperConnection(context);
		// Only probe for the target when the caller did not say. An explicit
		// target_column short-circuits that part entirely.
		if (data->target_column.empty()) {
			data->target_column = ResolveTargetColumn(conn, data->metadata_relation, data->caller);
		}
		// Captured here so the model can hand ids and labels back as the types
		// they arrived as. The sample type is not stored: every function that
		// returns sample ids takes its own data relation and mirrors that.
		// Validates both id columns; only the feature type is kept, since a fit
		// returns no sample ids and the model is what later calls read.
		data->feature_id_type = DetectCooIdTypes(conn, data->data_relation, data->caller).feature_id_type;
		data->target_type = DetectColumnType(conn, data->metadata_relation, data->target_column, data->caller);
	}

	// Required rather than optional. A name is not an input to the fit -- the
	// same seed with and without one produces a byte-identical model -- but
	// making it mandatory removes a state from the system: `name` is never NULL,
	// so `sc_predict(..., name := ...)` always works and nobody has to reason
	// about the unnamed case. Optional names get skipped, and then a registry
	// needs ALTER TABLE and rowid archaeology to become selectable.
	if (data->name.empty()) {
		throw InvalidInputException(
		    "sc_fit: name is required -- `name := 'my_model'`. It is how a model is selected later, e.g. "
		    "sc_predict(data, models, name := 'my_model'). Any label will do if you are just exploring.");
	}

	// sklearn rejects max_samples without bootstrap rather than silently ignoring
	// it; sc has no opinion, so catch it here where the message can be specific.
	if (!params.bootstrap && params.max_samples.kind != SC_MAX_SAMPLES_ALL) {
		throw InvalidInputException("sc_fit: max_samples is only meaningful with bootstrap := true");
	}

	// `name` is required, so this column is never NULL and a registry is
	// selectable from the moment it is created. Nothing is ever invented to fill
	// it -- a generated 'm0' would collide on the very next call, since a table
	// function cannot see what already exists in the table its row is headed for.
	//
	// `random_state` is echoed because it is the one thing needed to reproduce a
	// fit that the caller would otherwise have to remember and re-type. Richer
	// provenance stays the caller's: `SELECT 'x' AS target, * FROM sc_fit_...`.
	// `model` names the algorithm and mirrors the `model :=` parameter; the
	// serialized forest lives in `model_blob`. They were one column called
	// `model`, which meant `model := 'random_forest'` went in and bytes came
	// out under the same name -- and left a row unable to say which algorithm
	// produced it once there is more than one.
	// feature_id_type and target_type travel with the model because the functions
	// that return features or class labels -- sc_feature_importances,
	// sc_model_features, sc_shap, sc_predict -- have only the model to go on.
	names = {"name",    "model", "model_blob",   "n_samples",       "n_features",
	         "n_trees", "task",  "random_state", "feature_id_type", "target_type"};
	return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::BLOB,    LogicalType::BIGINT,
	                LogicalType::BIGINT,  LogicalType::BIGINT,  LogicalType::VARCHAR, LogicalType::BIGINT,
	                LogicalType::VARCHAR, LogicalType::VARCHAR};
	return std::move(data);
}

struct ScFitGlobalState : public GlobalTableFunctionState {
	bool done = false;
};

unique_ptr<GlobalTableFunctionState> ScFitInitGlobal(ClientContext &, TableFunctionInitInput &) {
	return make_uniq<ScFitGlobalState>();
}

/*
    Entry point called by DuckDB execution engine to pull output chunks.
    receives query context, function input state wrappers, and destination chunk
    the output we write to. In this case a single row is produced, that is the fitted model.
*/
void ScFitExecute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	// global state to tell us if we are finished
	auto &gstate = input.global_state->Cast<ScFitGlobalState>();
	if (gstate.done) {
		output.SetCardinality(0); // turn off the engine
		return;
	}
	// write to first byte in gstate (global duck db state) so that next time we come through here we know we are done
	// output.SetCardinality(0); is a contract for table functions to tell the engine that we are done producing output,
	// so we don't produce any more rows
	// executes ScFitExecute exactly once
	gstate.done = true;

	const auto &bind = input.bind_data->Cast<ScFitData>();
	auto conn = MakeReadOnlyHelperConnection(context);

	// instantiate the builder class on the function stack, intake triplets, canoncialise them, and produce a sparse
	// matrix.
	miint::CooBuilder builder;
	ScanCounts(conn, bind, builder);
	if (builder.NumNonZeros() == 0) {
		throw InvalidInputException("sc_fit: data relation '%s' produced no cells", bind.data_relation);
	}
	// remove duplicate triplets using sorted flat array of bitwise packed values
	// prefetchable, cache-friendly and less memory than a hash table. The builder's internal state is now a non
	// duplicated sparse matrix in COO format, ready to be packaged to arrow and passed to sc
	RequireNoDuplicateCells(builder, bind);

	// Targets are positional: element i must be the label for sample index i.
	// The order comes from the builder's sorted dictionary, looked up through a
	// hash -- never from a second ORDER BY, whose collation need not match
	// std::string's byte ordering.
	TargetArray targets;
	std::vector<std::string> data_samples;
	auto table = std::unique_ptr<miint::CooTable> {};

	if (bind.classification) {
		auto labels = ScanTargets<string>(conn, bind, "VARCHAR");
		table = builder.Finalize();
		// sorted sample ids after Finalize()
		data_samples = table->SampleIds();
		RequireSameSamples(data_samples, labels, bind);
		targets.labels.reserve(data_samples.size());
		for (const auto &s : data_samples) {
			targets.labels.push_back(labels.at(s));
		}
	} else {
		auto values = ScanTargets<double>(conn, bind, "DOUBLE");
		table = builder.Finalize();
		// sorted sample ids after Finalize()
		data_samples = table->SampleIds();
		RequireSameSamples(data_samples, values, bind);
		targets.numbers.reserve(data_samples.size());
		for (const auto &s : data_samples) {
			targets.numbers.push_back(values.at(s));
		}
	}
	BuildTargets(targets, bind.classification);

	// Fully resolved at bind time, so a bad parameter fails the query before any
	// scanning happens.
	const sc_rf_params_t &params = bind.params;

	// multi-threaded execution context for sc, with the number of threads specified by the user. If n_threads is 0, sc
	// will use all available threads.
	sc_config_t config {};
	config.n_threads = bind.n_threads;
	miint::ScContext ctx;
	if (auto st = sc_context_new(&config, &ctx.ptr); st != SC_OK) {
		miint::ThrowSc("sc_fit", nullptr, st);
	}

	// hand over the arrow data to sc, which will take ownership of the buffers and free them when done. The table is
	// now owned by sc and must not be freed by the caller. Fit the model and serialize it into a blob
	miint::ScModel model;
	// sc borrows these arrays, so the view may be a local: `table` stays the owner.
	const auto sc_table = miint::AsScTable(*table);
	const auto fit = bind.classification ? sc_fit_classifier : sc_fit_regressor;
	if (auto st = fit(ctx.ptr, &sc_table, &targets.array, &targets.schema, &params, &model.ptr); st != SC_OK) {
		miint::ThrowSc(bind.classification ? "sc_fit_classifier" : "sc_fit_regressor", ctx.ptr, st);
	}

	uint8_t *blob = nullptr;
	size_t blob_len = 0;
	// hand over **blob, one more level of indirection, to copy over the value of the memory addr mapped to the variable
	// blob and capture that inside the scope of the C function so it can dereference it and write to it.
	if (auto st = sc_model_serialize(model.ptr, &blob, &blob_len); st != SC_OK) {
		miint::ThrowSc("sc_model_serialize", ctx.ptr, st);
	}
	// sc owns this buffer until sc_buffer_free; copy it into DuckDB's heap first. Then free it from sc's heap. This is
	// a one-time copy, so the model blob is now owned by DuckDB.
	duckdb::Value blob_value = Value::BLOB(blob, blob_len);
	sc_buffer_free(blob, blob_len);

	// DuckDB's table function output is a single row with the fitted model, so we set the cardinality to 1 and fill in
	// the columns with the model's metadata and serialized blob.
	output.SetCardinality(1); // one row returned, which is the fitted model
	// The output columns are:
	output.SetValue(0, 0, Value(bind.name));
	output.SetValue(1, 0, Value(bind.model));
	output.SetValue(2, 0, blob_value);
	output.SetValue(3, 0, Value::BIGINT(table->NumSamples()));
	output.SetValue(4, 0, Value::BIGINT(table->NumFeatures()));
	output.SetValue(5, 0, Value::BIGINT(params.n_estimators));
	output.SetValue(6, 0, Value(bind.classification ? "classification" : "regression"));
	output.SetValue(7, 0, Value::BIGINT(static_cast<int64_t>(params.random_state)));
	output.SetValue(8, 0, Value(bind.feature_id_type.ToString()));
	output.SetValue(9, 0, Value(bind.target_type.ToString()));
}
/**
 * Creates a new table function for fitting a model based on sc rf
 */
TableFunction MakeFitFunction(const char *name, table_function_bind_t bind) {
	// declare a variable named fn and initialize an instance of the TableFunction class with the constructor call. The
	// TableFunction constructor takes the following parameters:
	TableFunction fn(name, {LogicalType::VARCHAR, LogicalType::VARCHAR}, ScFitExecute, bind, ScFitInitGlobal);
	fn.named_parameters["name"] = LogicalType::VARCHAR;
	fn.named_parameters["target_column"] = LogicalType::VARCHAR;
	fn.named_parameters["n_estimators"] = LogicalType::BIGINT;
	fn.named_parameters["random_state"] = LogicalType::BIGINT;
	fn.named_parameters["n_threads"] = LogicalType::INTEGER;
	fn.named_parameters["model"] = LogicalType::VARCHAR;
	fn.named_parameters["max_depth"] = LogicalType::BIGINT;
	fn.named_parameters["criterion"] = LogicalType::VARCHAR;
	fn.named_parameters["min_weight_fraction_leaf"] = LogicalType::DOUBLE;
	fn.named_parameters["min_impurity_decrease"] = LogicalType::DOUBLE;
	fn.named_parameters["bootstrap"] = LogicalType::BOOLEAN;
	// ANY, not a fixed type: these are sklearn's tagged unions, where an integer
	// means a count and a float means a fraction. SQL literal types carry that
	// distinction already, so ParseMaxFeatures / ParseMinSamples / ParseMaxSamples
	// dispatch on what the user wrote.
	fn.named_parameters["max_features"] = LogicalType::ANY;
	fn.named_parameters["min_samples_split"] = LogicalType::ANY;
	fn.named_parameters["min_samples_leaf"] = LogicalType::ANY;
	fn.named_parameters["max_samples"] = LogicalType::ANY;
	return fn;
}

unique_ptr<FunctionData> ScFitClassifierBind(ClientContext &context, TableFunctionBindInput &input,
                                             vector<LogicalType> &return_types, vector<string> &names) {
	return ScFitBind(context, input, return_types, names, /*classification=*/true);
}
unique_ptr<FunctionData> ScFitRegressorBind(ClientContext &context, TableFunctionBindInput &input,
                                            vector<LogicalType> &return_types, vector<string> &names) {
	return ScFitBind(context, input, return_types, names, /*classification=*/false);
}

} // namespace

void ScFitFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(MakeFitFunction("sc_fit_classifier", ScFitClassifierBind));
	loader.RegisterFunction(MakeFitFunction("sc_fit_regressor", ScFitRegressorBind));
}

} // namespace duckdb
