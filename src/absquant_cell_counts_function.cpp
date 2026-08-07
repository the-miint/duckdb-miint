#include "absquant.hpp"

#include "absquant_readers.hpp"
#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "miint_log.hpp"
#include "unifrac_function_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace duckdb {

namespace {

using absquant_internal::ReadKeyedColumns;
using absquant_internal::ReadLongFormValues;
using miint::absquant::CellCountsMetric;
using miint::absquant::CellCountsOptions;
using miint::absquant::CellCountsResult;
using miint::absquant::ComputeCellCounts;
using miint::absquant::CountObservation;
using miint::absquant::CoverageObservation;
using miint::absquant::DenominatorColumnName;
using miint::absquant::FeatureLength;
using miint::absquant::FeatureTableValue;
using miint::absquant::FormatIdList;
using miint::absquant::MetricName;
using miint::absquant::ParseCellCountsMetric;
using miint::absquant::SampleCellParams;
using miint::absquant::SampleRegression;
using unifrac_internal::ReadFeatureTable;
using unifrac_internal::ResolveSampleIdOutputType;

constexpr const char *kCallerName = "absquant_cell_counts";

// "'cells_per_g_of_gdna', 'cells_per_g_of_sample', ..." for the rejection
// message. Rendered from the enum rather than spelled out, so a fifth metric
// cannot be accepted by the parser and still be missing from the list of what
// this function accepts.
std::string AcceptedMetricList() {
	std::string out;
	for (const auto metric : miint::absquant::kAllCellCountsMetrics) {
		if (!out.empty()) {
			out += ", ";
		}
		out += std::string("'") + MetricName(metric) + "'";
	}
	return out;
}

struct AbsQuantCellCountsBindData : public TableFunctionData {
	std::string counts_table;
	std::string models_table;
	std::string coverage_table;
	std::string lengths_table;
	std::string params_table;
	CellCountsOptions options;
	// Output types for the two id columns, mirrored from the counts relation
	// (BIGINT/UUID preserved, else VARCHAR).
	LogicalType sample_id_type = LogicalType::VARCHAR;
	LogicalType feature_id_type = LogicalType::VARCHAR;
};

// The whole result is computed in InitGlobal, because the core needs every row
// of every relation before it can decide which samples survive the gates. What
// is left for Execute is a copy loop, so it stays single-threaded -- the same
// shape as absquant_fit_models and community_distances.
struct AbsQuantCellCountsGlobalState : public GlobalTableFunctionState {
	std::vector<FeatureTableValue> values;
	LogicalType sample_id_type = LogicalType::VARCHAR;
	LogicalType feature_id_type = LogicalType::VARCHAR;
	idx_t cursor = 0;
	idx_t MaxThreads() const override {
		return 1;
	}
};

unique_ptr<FunctionData> AbsQuantCellCountsBind(ClientContext &context, TableFunctionBindInput &input,
                                                vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<AbsQuantCellCountsBindData>();
	for (idx_t i = 0; i < 7; ++i) {
		if (input.inputs[i].IsNull()) {
			throw BinderException("%s: argument %llu must not be NULL", kCallerName, static_cast<uint64_t>(i + 1));
		}
	}
	data->counts_table = input.inputs[0].GetValue<string>();
	data->models_table = input.inputs[1].GetValue<string>();
	data->coverage_table = input.inputs[2].GetValue<string>();
	data->lengths_table = input.inputs[3].GetValue<string>();
	data->params_table = input.inputs[4].GetValue<string>();
	data->options.min_coverage = input.inputs[6].GetValue<double>();
	for (const auto &kv : input.named_parameters) {
		if (StringUtil::Lower(kv.first) == "min_rsquared") {
			if (kv.second.IsNull()) {
				throw BinderException("%s: min_rsquared must not be NULL", kCallerName);
			}
			data->options.min_rsquared = kv.second.GetValue<double>();
		}
	}

	// The metric IS checked here, unlike the two numeric thresholds. Their range
	// checks live in the pure core, where they are unit-tested against pysyndna's
	// bounds and InitGlobal re-throws them with this function's name. An
	// unparseable metric never gets as far as the core -- it decides which
	// columns to read -- so bind is the only place that can reject it.
	//
	// ParseCellCountsMetric folds case itself, so the raw argument goes in and
	// the raw argument comes back out in the message.
	if (!ParseCellCountsMetric(input.inputs[5].GetValue<string>(), data->options.metric)) {
		throw BinderException("%s: unsupported metric '%s'; accepted: %s", kCallerName,
		                      input.inputs[5].GetValue<string>(), AcceptedMetricList());
	}

	// All five relations are resolved here so a typo in ANY of them is a
	// bind-time error, matching absquant_fit_models' handling of its three. Only
	// the counts table's columns are needed at bind, for the id-type mirroring
	// below; the other four are looked up purely to fail early, and their column
	// contracts are enforced by the LIMIT 0 probes in InitGlobal, whose errors
	// carry DuckDB's own candidate-binding hints.
	//
	// "table", not "relation", in the entity strings: GetTableOrViewColumns
	// renders them as "<entity> or view '<name>' does not exist".
	auto cols = GetTableOrViewColumns(context, data->counts_table, "feature counts table");
	GetTableOrViewColumns(context, data->models_table, "models table");
	GetTableOrViewColumns(context, data->coverage_table, "coverage table");
	GetTableOrViewColumns(context, data->lengths_table, "feature lengths table");
	GetTableOrViewColumns(context, data->params_table, "sample parameters table");

	// Both id columns are mirrored, and both from the COUNTS relation, which is
	// the one relation guaranteed to name every id that can appear in the output.
	// ResolveSampleIdOutputType says "sample" but the rule it encodes is
	// id-generic -- BIGINT and UUID pass through, everything else becomes VARCHAR
	// -- and it is deliberately permissive, so it never rejects a column the
	// ::VARCHAR-casting reader would have happily read.
	for (idx_t i = 0; i < cols.names.size(); ++i) {
		const auto lname = StringUtil::Lower(cols.names[i]);
		if (lname == "sample_id") {
			data->sample_id_type = ResolveSampleIdOutputType(cols.types[i]);
		} else if (lname == "feature_id") {
			data->feature_id_type = ResolveSampleIdOutputType(cols.types[i]);
		}
	}

	names = {"sample_id", "feature_id", "value"};
	return_types = {data->sample_id_type, data->feature_id_type, LogicalType::DOUBLE};
	return std::move(data);
}

// Every sample the caller handed us that produced no cells is reported, and so
// is every feature dropped for coverage. A missing row is the one outcome a user
// cannot debug from the output alone.
//
// pysyndna logs the first two of these and is silent about the other three: a
// wholly uncovered sample never reaches its per-sample loop, it has no notion of
// a model that is present but unusable, and its dense output writes the
// all-zero sample out rather than omitting it. The values agree either way;
// those three lines are miint saying out loud what pysyndna leaves the user to
// infer from an absent row.
void EmitDiagnostics(ClientContext &context, const AbsQuantCellCountsBindData &data, const CellCountsResult &result) {
	if (!result.low_coverage_feature_ids.empty()) {
		miint::EmitWarning(
		    context, std::string(kCallerName) + ": dropped " + std::to_string(result.low_coverage_feature_ids.size()) +
		                 " feature(s) whose coverage fell below " + std::to_string(data.options.min_coverage) +
		                 " in at least one sample: " + FormatIdList(result.low_coverage_feature_ids));
	}
	if (!result.filtered_sample_ids.empty()) {
		miint::EmitWarning(context, std::string(kCallerName) + ": dropped " +
		                                std::to_string(result.filtered_sample_ids.size()) +
		                                " sample(s) with a NULL, negative or non-finite required parameter: " +
		                                FormatIdList(result.filtered_sample_ids));
	}
	if (!result.uncovered_sample_ids.empty()) {
		miint::EmitWarning(
		    context, std::string(kCallerName) + ": dropped " + std::to_string(result.uncovered_sample_ids.size()) +
		                 " sample(s) with no feature at or above " + std::to_string(data.options.min_coverage) +
		                 " coverage: " + FormatIdList(result.uncovered_sample_ids));
	}
	if (!result.samples_without_models.empty()) {
		miint::EmitWarning(context, std::string(kCallerName) + ": no usable model in '" + data.models_table + "' for " +
		                                std::to_string(result.samples_without_models.size()) +
		                                " sample(s): " + FormatIdList(result.samples_without_models));
	}
	if (!result.low_rsquared_sample_ids.empty()) {
		miint::EmitWarning(
		    context, std::string(kCallerName) + ": dropped " + std::to_string(result.low_rsquared_sample_ids.size()) +
		                 " sample(s) whose model r^2 fell below " + std::to_string(data.options.min_rsquared) + ": " +
		                 FormatIdList(result.low_rsquared_sample_ids));
	}
	if (!result.zero_valued_sample_ids.empty()) {
		// These passed every gate and still contribute no row, because the
		// output is sparse and every one of their cells is zero. Naming the
		// denominator is the actionable part: the cause is almost always a zero
		// extraction concentration or elution volume, which is exactly the case
		// the NULL/negative screen lets through.
		miint::EmitWarning(
		    context, std::string(kCallerName) + ": " + std::to_string(result.zero_valued_sample_ids.size()) +
		                 " sample(s) produced only zero-valued cells for metric '" + MetricName(data.options.metric) +
		                 "' and appear in no output row: " + FormatIdList(result.zero_valued_sample_ids));
	}
}

unique_ptr<GlobalTableFunctionState> AbsQuantCellCountsInitGlobal(ClientContext &context,
                                                                  TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<AbsQuantCellCountsBindData>();
	auto gstate = make_uniq<AbsQuantCellCountsGlobalState>();
	gstate->sample_id_type = data.sample_id_type;
	gstate->feature_id_type = data.feature_id_type;

	// ReadFeatureTable's zero-drop is what makes the log10 in the core safe, and
	// what makes the output sparse where pysyndna's is a dense zero (D10). The
	// coverage relation must NOT go through it -- zero is a real coverage -- so it
	// uses the zero-preserving reader instead. See absquant_readers.hpp.
	const auto count_rows = ReadFeatureTable(context, data.counts_table, kCallerName);
	const auto model_rows = ReadKeyedColumns(context, data.models_table, "sample_id", {"slope", "intercept", "rvalue"},
	                                         "models relation", kCallerName);
	const auto coverage_rows =
	    ReadLongFormValues(context, data.coverage_table, "coverage", "coverage relation", kCallerName);
	const auto length_rows = ReadKeyedColumns(context, data.lengths_table, "feature_id", {"ogu_len_in_bp"},
	                                          "feature lengths relation", kCallerName);
	// The requested metric's denominator joins the three base columns, and only
	// it does. Asking for cells_per_g_of_gdna must not require a sample_volume_ul
	// column to exist -- the LIMIT 0 probe inside ReadKeyedColumns would reject
	// the whole relation for a column this query never reads.
	const char *denominator_column = DenominatorColumnName(data.options.metric);
	std::vector<const char *> param_columns = {"sequenced_sample_gdna_mass_ng", "extracted_gdna_concentration_ng_ul",
	                                           "vol_extracted_elution_ul"};
	if (denominator_column != nullptr) {
		param_columns.push_back(denominator_column);
	}
	const auto param_rows = ReadKeyedColumns(context, data.params_table, "sample_id", param_columns,
	                                         "sample parameters relation", kCallerName);

	std::vector<CountObservation> counts;
	counts.reserve(count_rows.size());
	for (const auto &row : count_rows) {
		counts.push_back({row.sample_id, row.feature_id, row.count});
	}
	std::vector<SampleRegression> models;
	models.reserve(model_rows.keys.size());
	const auto &slope = model_rows.values.at("slope");
	const auto &intercept = model_rows.values.at("intercept");
	const auto &rvalue = model_rows.values.at("rvalue");
	for (size_t i = 0; i < model_rows.keys.size(); ++i) {
		models.push_back({model_rows.keys[i], slope[i], intercept[i], rvalue[i]});
	}
	std::vector<CoverageObservation> coverage;
	coverage.reserve(coverage_rows.size());
	for (const auto &row : coverage_rows) {
		coverage.push_back({row.sample_id, row.feature_id, row.value});
	}
	std::vector<FeatureLength> lengths;
	lengths.reserve(length_rows.keys.size());
	const auto &ogu_len_in_bp = length_rows.values.at("ogu_len_in_bp");
	for (size_t i = 0; i < length_rows.keys.size(); ++i) {
		lengths.push_back({length_rows.keys[i], ogu_len_in_bp[i]});
	}
	std::vector<SampleCellParams> params;
	params.reserve(param_rows.keys.size());
	const auto &gdna_mass_ng = param_rows.values.at("sequenced_sample_gdna_mass_ng");
	const auto &gdna_conc = param_rows.values.at("extracted_gdna_concentration_ng_ul");
	const auto &elution_ul = param_rows.values.at("vol_extracted_elution_ul");
	// Null for cells_per_g_of_gdna, which leaves sample_denominator at its zero
	// default -- the core does not read it and its validation skips it.
	const std::vector<double> *denominator =
	    denominator_column == nullptr ? nullptr : &param_rows.values.at(denominator_column);
	for (size_t i = 0; i < param_rows.keys.size(); ++i) {
		params.push_back({param_rows.keys[i], gdna_mass_ng[i], gdna_conc[i], elution_ul[i],
		                  denominator == nullptr ? 0.0 : (*denominator)[i]});
	}

	// The pure core carries the whole rejection layer and throws
	// std::invalid_argument with no function name attached, exactly as
	// BuildDenseDistanceMatrix does for ReadDistanceTable. Prefixing here is what
	// turns those into SQL errors a user can act on.
	CellCountsResult result;
	try {
		result = ComputeCellCounts(counts, models, coverage, lengths, params, data.options);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s: %s", kCallerName, e.what());
	}

	EmitDiagnostics(context, data, result);
	gstate->values = std::move(result.values);
	return std::move(gstate);
}

void AbsQuantCellCountsExecute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &g = data_p.global_state->Cast<AbsQuantCellCountsGlobalState>();
	const idx_t total = g.values.size();
	const idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - g.cursor);

	auto &sample_id = output.data[0];
	auto &feature_id = output.data[1];
	auto value = FlatVector::GetData<double>(output.data[2]);

	for (idx_t r = 0; r < count; ++r) {
		const auto &cell = g.values[g.cursor + r];
		EmitIdCell(sample_id, r, cell.sample_id, g.sample_id_type);
		EmitIdCell(feature_id, r, cell.feature_id, g.feature_id_type);
		value[r] = cell.value;
	}
	g.cursor += count;
	output.SetCardinality(count);
}

} // namespace

void RegisterAbsQuantCellCounts(ExtensionLoader &loader) {
	TableFunction cells("absquant_cell_counts",
	                    {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
	                     LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::DOUBLE},
	                    AbsQuantCellCountsExecute, AbsQuantCellCountsBind, AbsQuantCellCountsInitGlobal);
	cells.named_parameters["min_rsquared"] = LogicalType::DOUBLE;
	cells.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(cells);
}

} // namespace duckdb
