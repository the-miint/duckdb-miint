#include "sc_cross_validate_function.hpp"

#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "sc_common.hpp"
#include "coo_builder.hpp"
#include "sc_rf_common.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/main/connection.hpp"

#include <memory>
#include <string>
#include <vector>

#include "sc.h"

namespace duckdb {

namespace {

// Everything shared with sc_fit -- input scanning, validation, forest-parameter
// parsing -- lives in sc_rf_common.
using namespace sc_rf;

//! q2's `cv` default, and sklearn's for cross_val_score.
constexpr int64_t kDefaultCvFolds = 5;

struct ScCvData : public TableFunctionData, public ScTrainingInput {
	bool classification = false;
	int32_t n_threads = 0;
	sc_rf_params_t params {};
};

struct ScCvGlobalState : public GlobalTableFunctionState {
	std::vector<std::string> sample_ids;
	//! OOF prediction per sample, in sample_ids order. One of the two is filled.
	std::vector<std::string> labels;
	std::vector<double> values;
	//! The targets the folds were scored against, same order. Already read to
	//! hand to sc, so carrying them costs nothing and spares every caller a join
	//! back to the metadata relation to score or plot the result.
	std::vector<std::string> actual_labels;
	std::vector<double> actual_values;
	//! OOF class probabilities, row-major n_samples x n_classes, in `classes`
	//! order. Empty for a regressor. Flat rather than one Value per sample so the
	//! rows are built only as they are emitted.
	std::vector<double> proba;
	size_t n_classes = 0;
	//! Constant for the whole run, so every row repeats them.
	Value classes;
	Value fold_scores;
	double mean_score = 0;
	double std_score = 0;

	size_t cursor = 0;
	bool done_setup = false;
};

unique_ptr<FunctionData> ScCvBind(ClientContext &context, TableFunctionBindInput &input,
                                  vector<LogicalType> &return_types, vector<string> &names, bool classification) {
	auto data = make_uniq<ScCvData>();
	data->caller = classification ? "sc_cross_validate_classifier" : "sc_cross_validate_regressor";
	data->classification = classification;
	data->data_relation = input.inputs[0].GetValue<string>();
	data->metadata_relation = input.inputs[1].GetValue<string>();
	if (data->data_relation.empty() || data->metadata_relation.empty()) {
		throw InvalidInputException("%s: data and metadata relation names must not be empty", data->caller);
	}

	ApplyDefaults(data->params, classification);
	auto &params = data->params;
	// `cv` is what q2, sklearn's cross_val_score and sc's own C struct call this,
	// so it keeps working; `n_folds` says what the number actually is. Accepting
	// both but never silently resolving a disagreement.
	bool n_folds_given = false;
	bool cv_given = false;
	params.cv = kDefaultCvFolds;
	for (auto &kv : input.named_parameters) {
		const auto &k = kv.first;
		const auto &v = kv.second;
		if (v.IsNull()) {
			throw InvalidInputException("%s: named parameter '%s' must not be NULL", data->caller, k);
		}
		if (StringUtil::CIEquals(k, "target_column")) {
			data->target_column = v.GetValue<string>();
		} else if (StringUtil::CIEquals(k, "n_folds") || StringUtil::CIEquals(k, "cv")) {
			(StringUtil::CIEquals(k, "cv") ? cv_given : n_folds_given) = true;
			params.cv = v.GetValue<int64_t>();
			// One fold has nothing to hold out: the model would be scored on the
			// data it was fit on, which is the thing cross-validation exists to
			// avoid.
			if (params.cv < 2) {
				throw InvalidInputException("%s: n_folds must be >= 2 (got %lld)", data->caller, (long long)params.cv);
			}
		} else if (StringUtil::CIEquals(k, "parameter_tuning")) {
			params.parameter_tuning = v.GetValue<bool>();
		} else if (StringUtil::CIEquals(k, "n_threads")) {
			data->n_threads = v.GetValue<int32_t>();
			params.n_threads = data->n_threads;
		} else if (StringUtil::CIEquals(k, "n_estimators")) {
			params.n_estimators = v.GetValue<int64_t>();
			if (params.n_estimators <= 0) {
				throw InvalidInputException("%s: n_estimators must be > 0 (got %lld)", data->caller,
				                            (long long)params.n_estimators);
			}
		} else if (StringUtil::CIEquals(k, "random_state")) {
			params.random_state = static_cast<uint64_t>(v.GetValue<int64_t>());
		} else if (StringUtil::CIEquals(k, "max_depth")) {
			params.max_depth = v.GetValue<int64_t>();
		} else if (StringUtil::CIEquals(k, "criterion")) {
			ParseCriterion(v, classification, "sc_cross_validate", params.criterion);
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
	if (n_folds_given && cv_given) {
		throw InvalidInputException("%s: pass n_folds or cv, not both -- cv is the q2/sklearn spelling of n_folds",
		                            data->caller);
	}
	{
		auto conn = MakeReadOnlyHelperConnection(context);
		if (data->target_column.empty()) {
			data->target_column = ResolveTargetColumn(conn, data->metadata_relation, data->caller);
		}
		// No model is involved, so the types come straight from this call's own
		// relations -- ids from the data, labels from the metadata column.
		// Validates both id columns; cross-validation returns sample ids only.
		data->sample_id_type = DetectCooIdTypes(conn, data->data_relation, data->caller).sample_id_type;
		data->target_type = DetectColumnType(conn, data->metadata_relation, data->target_column, data->caller);
	}
	// sklearn rejects max_samples without bootstrap rather than ignoring it.
	if (!params.bootstrap && params.max_samples.kind != SC_MAX_SAMPLES_ALL) {
		throw InvalidInputException("%s: max_samples is only meaningful with bootstrap := true", data->caller);
	}

	// No model comes out, so there is nothing to name and nothing to store: the
	// per-fold forests exist only to be scored. The prediction column follows the
	// task, exactly as in sc_predict.
	names = {"sample_id", "prediction", "actual"};
	// A classifier hands back the labels it was given; a regressor predicts a
	// continuous value, so both its columns are DOUBLE whatever the column was.
	const auto value_type = classification ? data->target_type : LogicalType::DOUBLE;
	return_types = {data->sample_id_type, value_type, value_type};
	if (classification) {
		// Positional, not a map: probabilities[i] is the probability of classes[i],
		// and `classes` is the same list on every row, so a whole column reads as
		// one matrix. The regressor has no classes to be probable, and the task is
		// fixed by which function was called, so its schema simply lacks these
		// rather than carrying two always-NULL columns.
		names.push_back("probabilities");
		return_types.push_back(LogicalType::LIST(LogicalType::DOUBLE));
		names.push_back("classes");
		return_types.push_back(LogicalType::LIST(data->target_type));
	}
	for (const auto &n : {"metric", "mean_score", "std_score", "fold_scores"}) {
		names.push_back(n);
	}
	return_types.push_back(LogicalType::VARCHAR);
	return_types.push_back(LogicalType::DOUBLE);
	return_types.push_back(LogicalType::DOUBLE);
	return_types.push_back(LogicalType::LIST(LogicalType::DOUBLE));
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ScCvInitGlobal(ClientContext &, TableFunctionInitInput &) {
	return make_uniq<ScCvGlobalState>();
}

//! Take ownership of one array sc wrote into sc_cv_result_t.
//!
//! The C Data Interface permits moving these structs: copy the bytes, then mark
//! the source released so nothing frees it twice. Done for every field before
//! the status is checked, so an error path still cleans up.
void TakeArray(ArrowArray &src, ArrowSchema &src_schema, miint::OwnedArrowArray &dst) {
	*dst.array() = src;
	*dst.schema() = src_schema;
	src.release = nullptr;
	src_schema.release = nullptr;
}

//! Scan both relations, run the folds, and keep what the rows need.
void RunCrossValidation(ClientContext &context, const ScCvData &bind, ScCvGlobalState &gstate) {
	auto conn = MakeReadOnlyHelperConnection(context);

	miint::CooBuilder builder;
	ScanCounts(conn, bind, builder);
	if (builder.NumNonZeros() == 0) {
		throw InvalidInputException("%s: data relation '%s' produced no cells", bind.caller, bind.data_relation);
	}
	RequireNoDuplicateCells(builder, bind);

	// Targets are positional: element i is the target of sample index i, in the
	// builder's sorted dictionary order.
	TargetArray targets;
	std::unique_ptr<miint::CooTable> table;
	if (bind.classification) {
		auto labels = ScanTargets<string>(conn, bind, "VARCHAR");
		table = builder.Finalize();
		gstate.sample_ids = table->SampleIds();
		RequireSameSamples(gstate.sample_ids, labels, bind);
		targets.labels.reserve(gstate.sample_ids.size());
		for (const auto &s : gstate.sample_ids) {
			targets.labels.push_back(labels.at(s));
		}
	} else {
		auto values = ScanTargets<double>(conn, bind, "DOUBLE");
		table = builder.Finalize();
		gstate.sample_ids = table->SampleIds();
		RequireSameSamples(gstate.sample_ids, values, bind);
		targets.numbers.reserve(gstate.sample_ids.size());
		for (const auto &s : gstate.sample_ids) {
			targets.numbers.push_back(values.at(s));
		}
	}
	BuildTargets(targets, bind.classification);

	const auto n_samples = gstate.sample_ids.size();
	// Caught here rather than in sc, where the message would name a matrix rather
	// than the relation the user passed.
	if (static_cast<size_t>(bind.params.cv) > n_samples) {
		throw InvalidInputException("%s: n_folds is %lld but '%s' has only %llu samples; a fold cannot be empty",
		                            bind.caller, (long long)bind.params.cv, bind.data_relation,
		                            (unsigned long long)n_samples);
	}

	sc_config_t config {};
	config.n_threads = bind.n_threads;
	miint::ScContext ctx;
	if (auto st = sc_context_new(&config, &ctx.ptr); st != SC_OK) {
		miint::ThrowSc(bind.caller, nullptr, st);
	}

	const auto sc_table = miint::AsScTable(*table);
	sc_cv_result_t res {};
	const auto st = sc_cross_validate(ctx.ptr, &sc_table, &targets.array, &targets.schema, &bind.params, &res);
	// Every field, including the probabilities this function does not return:
	// unclaimed arrays would leak.
	miint::OwnedArrowArray predictions, proba, classes, fold_scores;
	TakeArray(res.predictions, res.predictions_schema, predictions);
	TakeArray(res.proba, res.proba_schema, proba);
	TakeArray(res.classes, res.classes_schema, classes);
	TakeArray(res.fold_scores, res.fold_scores_schema, fold_scores);
	if (st != SC_OK) {
		miint::ThrowSc(bind.caller, ctx.ptr, st);
	}

	// Safe only now that sc has read them: the Arrow arrays it borrowed point
	// into these buffers.
	gstate.actual_labels = std::move(targets.labels);
	gstate.actual_values = std::move(targets.numbers);

	if (bind.classification) {
		gstate.labels = predictions.ReadUtf8("sc_cross_validate predictions");
		// The OOF probability behind each prediction, which is what ROC curves,
		// calibration and "which samples was it least sure about" need. sc aligns
		// the class array to the probability columns; that order is carried into
		// every row rather than re-derived here.
		const auto class_labels = classes.ReadUtf8("sc_cross_validate classes");
		int64_t width = 0;
		gstate.proba = proba.ReadFixedSizeListFloat64("sc_cross_validate proba", width);
		gstate.n_classes = class_labels.size();
		if (static_cast<size_t>(width) != gstate.n_classes ||
		    gstate.proba.size() != gstate.sample_ids.size() * gstate.n_classes) {
			throw InternalException("sc_cross_validate: %llu probabilities of width %lld for %llu samples x %llu "
			                        "classes",
			                        (unsigned long long)gstate.proba.size(), (long long)width,
			                        (unsigned long long)gstate.sample_ids.size(), (unsigned long long)gstate.n_classes);
		}
		duckdb::vector<Value> class_values;
		class_values.reserve(class_labels.size());
		for (const auto &c : class_labels) {
			class_values.push_back(Value(c).DefaultCastAs(bind.target_type));
		}
		gstate.classes = Value::LIST(bind.target_type, std::move(class_values));
	} else {
		gstate.values = predictions.ReadFloat64("sc_cross_validate predictions");
	}
	const auto scores = fold_scores.ReadFloat64("sc_cross_validate fold_scores");
	const auto n_predictions = bind.classification ? gstate.labels.size() : gstate.values.size();
	if (n_predictions != n_samples || scores.size() != static_cast<size_t>(bind.params.cv)) {
		throw InternalException("sc_cross_validate: %llu predictions and %llu fold scores for %llu samples over "
		                        "%lld folds",
		                        (unsigned long long)n_predictions, (unsigned long long)scores.size(),
		                        (unsigned long long)n_samples, (long long)bind.params.cv);
	}
	duckdb::vector<Value> fold_values;
	fold_values.reserve(scores.size());
	for (const auto s : scores) {
		fold_values.push_back(Value::DOUBLE(s));
	}
	gstate.fold_scores = Value::LIST(LogicalType::DOUBLE, std::move(fold_values));
	gstate.mean_score = res.mean_score;
	gstate.std_score = res.std_score;
}

void ScCvExecute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<ScCvGlobalState>();
	const auto &bind = input.bind_data->Cast<ScCvData>();
	if (!gstate.done_setup) {
		gstate.done_setup = true;
		RunCrossValidation(context, bind, gstate);
	}

	const auto metric = bind.classification ? "accuracy" : "mse";
	const size_t total = gstate.sample_ids.size();
	idx_t n = 0;
	while (n < STANDARD_VECTOR_SIZE && gstate.cursor < total) {
		const auto i = gstate.cursor++;
		EmitIdCell(output.data[0], n, gstate.sample_ids[i], bind.sample_id_type);
		output.SetValue(1, n,
		                bind.classification ? Value(gstate.labels[i]).DefaultCastAs(bind.target_type)
		                                    : Value::DOUBLE(gstate.values[i]));
		output.SetValue(2, n,
		                bind.classification ? Value(gstate.actual_labels[i]).DefaultCastAs(bind.target_type)
		                                    : Value::DOUBLE(gstate.actual_values[i]));
		// The probability columns exist only for a classifier, so everything after
		// them shifts by two.
		idx_t col = 3;
		if (bind.classification) {
			duckdb::vector<Value> row;
			row.reserve(gstate.n_classes);
			for (size_t c = 0; c < gstate.n_classes; c++) {
				row.push_back(Value::DOUBLE(gstate.proba[i * gstate.n_classes + c]));
			}
			output.SetValue(col++, n, Value::LIST(LogicalType::DOUBLE, std::move(row)));
			output.SetValue(col++, n, gstate.classes);
		}
		output.SetValue(col++, n, Value(metric));
		output.SetValue(col++, n, Value::DOUBLE(gstate.mean_score));
		output.SetValue(col++, n, Value::DOUBLE(gstate.std_score));
		output.SetValue(col, n, gstate.fold_scores);
		n++;
	}
	output.SetCardinality(n);
}

TableFunction MakeCvFunction(const char *name, table_function_bind_t bind) {
	TableFunction fn(name, {LogicalType::VARCHAR, LogicalType::VARCHAR}, ScCvExecute, bind, ScCvInitGlobal);
	fn.named_parameters["target_column"] = LogicalType::VARCHAR;
	fn.named_parameters["n_folds"] = LogicalType::BIGINT;
	fn.named_parameters["cv"] = LogicalType::BIGINT;
	fn.named_parameters["parameter_tuning"] = LogicalType::BOOLEAN;
	fn.named_parameters["n_estimators"] = LogicalType::BIGINT;
	fn.named_parameters["random_state"] = LogicalType::BIGINT;
	fn.named_parameters["n_threads"] = LogicalType::INTEGER;
	fn.named_parameters["max_depth"] = LogicalType::BIGINT;
	fn.named_parameters["criterion"] = LogicalType::VARCHAR;
	fn.named_parameters["min_weight_fraction_leaf"] = LogicalType::DOUBLE;
	fn.named_parameters["min_impurity_decrease"] = LogicalType::DOUBLE;
	fn.named_parameters["bootstrap"] = LogicalType::BOOLEAN;
	// ANY for sklearn's tagged unions: an integer means a count, a float a
	// fraction. See MakeFitFunction.
	fn.named_parameters["max_features"] = LogicalType::ANY;
	fn.named_parameters["min_samples_split"] = LogicalType::ANY;
	fn.named_parameters["min_samples_leaf"] = LogicalType::ANY;
	fn.named_parameters["max_samples"] = LogicalType::ANY;
	return fn;
}

unique_ptr<FunctionData> ScCvClassifierBind(ClientContext &context, TableFunctionBindInput &input,
                                            vector<LogicalType> &return_types, vector<string> &names) {
	return ScCvBind(context, input, return_types, names, /*classification=*/true);
}
unique_ptr<FunctionData> ScCvRegressorBind(ClientContext &context, TableFunctionBindInput &input,
                                           vector<LogicalType> &return_types, vector<string> &names) {
	return ScCvBind(context, input, return_types, names, /*classification=*/false);
}

} // namespace

void ScCrossValidateFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(MakeCvFunction("sc_cross_validate_classifier", ScCvClassifierBind));
	loader.RegisterFunction(MakeCvFunction("sc_cross_validate_regressor", ScCvRegressorBind));
}

} // namespace duckdb
