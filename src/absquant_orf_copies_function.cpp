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
using miint::absquant::ComputeOrfCopies;
using miint::absquant::CountObservation;
using miint::absquant::FeatureTableValue;
using miint::absquant::FormatIdList;
using miint::absquant::OrfCoords;
using miint::absquant::OrfCopiesResult;
using miint::absquant::SampleOrfParams;
using unifrac_internal::ReadFeatureTable;
using unifrac_internal::ResolveSampleIdOutputType;

constexpr const char *kCallerName = "absquant_orf_copies";

struct AbsQuantOrfCopiesBindData : public TableFunctionData {
	std::string counts_table;
	std::string coords_table;
	std::string params_table;
	// Output types for the two id columns, mirrored from the counts relation
	// (BIGINT/UUID preserved, else VARCHAR).
	LogicalType sample_id_type = LogicalType::VARCHAR;
	LogicalType feature_id_type = LogicalType::VARCHAR;
};

// The whole result is computed in InitGlobal, because the core needs every row
// of every relation before it can decide which samples survive the filter. What
// is left for Execute is a copy loop, so it stays single-threaded -- the same
// shape as absquant_cell_counts and community_distances.
struct AbsQuantOrfCopiesGlobalState : public GlobalTableFunctionState {
	std::vector<FeatureTableValue> values;
	LogicalType sample_id_type = LogicalType::VARCHAR;
	LogicalType feature_id_type = LogicalType::VARCHAR;
	idx_t cursor = 0;
	idx_t MaxThreads() const override {
		return 1;
	}
};

unique_ptr<FunctionData> AbsQuantOrfCopiesBind(ClientContext &context, TableFunctionBindInput &input,
                                               vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<AbsQuantOrfCopiesBindData>();
	for (idx_t i = 0; i < 3; ++i) {
		if (input.inputs[i].IsNull()) {
			throw BinderException("%s: argument %llu must not be NULL", kCallerName, static_cast<uint64_t>(i + 1));
		}
	}
	data->counts_table = input.inputs[0].GetValue<string>();
	data->coords_table = input.inputs[1].GetValue<string>();
	data->params_table = input.inputs[2].GetValue<string>();

	// No options and no named parameters: pysyndna exposes exactly one metric
	// for this workflow, and there is no threshold to set -- no standard curve,
	// no coverage filter, no r^2 gate. Nothing here needs the bind-time
	// rejection layer absquant_cell_counts has for its metric argument.
	//
	// All three relations are resolved here so a typo in ANY of them is a
	// bind-time error, matching the other two absquant functions. Only the
	// counts table's columns are needed at bind, for the id-type mirroring
	// below; the other two are looked up purely to fail early, and their column
	// contracts are enforced by the LIMIT 0 probes in InitGlobal, whose errors
	// carry DuckDB's own candidate-binding hints.
	//
	// "table", not "relation", in the entity strings: GetTableOrViewColumns
	// renders them as "<entity> or view '<name>' does not exist".
	auto cols = GetTableOrViewColumns(context, data->counts_table, "feature counts table");
	GetTableOrViewColumns(context, data->coords_table, "ORF coordinates table");
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

// Every sample the caller handed us that produced no cells is reported. A
// missing row is the one outcome a user cannot debug from the output alone, and
// this function has exactly two ways to produce one.
//
// pysyndna logs the first of these and is silent about the second: its output is
// dense, so an all-zero sample is written out as zeros rather than omitted. The
// values agree either way; that line is miint saying out loud what the sparse
// form would otherwise leave the user to infer from an absent row.
void EmitDiagnostics(ClientContext &context, const OrfCopiesResult &result) {
	if (!result.filtered_sample_ids.empty()) {
		miint::EmitWarning(context, std::string(kCallerName) + ": dropped " +
		                                std::to_string(result.filtered_sample_ids.size()) +
		                                " sample(s) with a NULL, negative or non-finite required parameter: " +
		                                FormatIdList(result.filtered_sample_ids));
	}
	if (!result.zero_valued_sample_ids.empty()) {
		// These passed the filter and still contribute no row, because the output
		// is sparse and every one of their cells is zero. Naming the likely cause
		// is the actionable part: a zero total_rna_concentration_ng_ul or
		// vol_extracted_elution_ul, which is exactly what the NULL/negative screen
		// lets through.
		miint::EmitWarning(context, std::string(kCallerName) + ": " +
		                                std::to_string(result.zero_valued_sample_ids.size()) +
		                                " sample(s) produced only zero-valued cells -- check "
		                                "total_rna_concentration_ng_ul and vol_extracted_elution_ul -- and appear in "
		                                "no output row: " +
		                                FormatIdList(result.zero_valued_sample_ids));
	}
}

unique_ptr<GlobalTableFunctionState> AbsQuantOrfCopiesInitGlobal(ClientContext &context,
                                                                 TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<AbsQuantOrfCopiesBindData>();
	auto gstate = make_uniq<AbsQuantOrfCopiesGlobalState>();
	gstate->sample_id_type = data.sample_id_type;
	gstate->feature_id_type = data.feature_id_type;

	// ReadFeatureTable's zero-drop is what makes the output sparse where
	// pysyndna's is a dense zero (D10). Unlike absquant_cell_counts nothing here
	// DEPENDS on it -- a zero count is legal in this chain, there being no
	// log10 -- so the drop is purely the storage convention.
	const auto count_rows = ReadFeatureTable(context, data.counts_table, kCallerName);
	const auto coord_rows = ReadKeyedColumns(context, data.coords_table, "feature_id", {"ogu_orf_start", "ogu_orf_end"},
	                                         "ORF coordinates relation", kCallerName);
	const auto param_rows = ReadKeyedColumns(context, data.params_table, "sample_id",
	                                         {"calc_mass_sample_aliquot_input_g", "total_rna_concentration_ng_ul",
	                                          "vol_extracted_elution_ul", "total_biological_reads_r1r2"},
	                                         "sample parameters relation", kCallerName);

	std::vector<CountObservation> counts;
	counts.reserve(count_rows.size());
	for (const auto &row : count_rows) {
		counts.push_back({row.sample_id, row.feature_id, row.count});
	}
	std::vector<OrfCoords> coords;
	coords.reserve(coord_rows.keys.size());
	const auto &start = coord_rows.values.at("ogu_orf_start");
	const auto &end = coord_rows.values.at("ogu_orf_end");
	for (size_t i = 0; i < coord_rows.keys.size(); ++i) {
		coords.push_back({coord_rows.keys[i], start[i], end[i]});
	}
	std::vector<SampleOrfParams> params;
	params.reserve(param_rows.keys.size());
	const auto &aliquot_mass_g = param_rows.values.at("calc_mass_sample_aliquot_input_g");
	const auto &rna_conc = param_rows.values.at("total_rna_concentration_ng_ul");
	const auto &elution_ul = param_rows.values.at("vol_extracted_elution_ul");
	const auto &biological_reads = param_rows.values.at("total_biological_reads_r1r2");
	for (size_t i = 0; i < param_rows.keys.size(); ++i) {
		params.push_back({param_rows.keys[i], aliquot_mass_g[i], rna_conc[i], elution_ul[i], biological_reads[i]});
	}

	// The pure core carries the whole rejection layer and throws
	// std::invalid_argument with no function name attached, exactly as
	// BuildDenseDistanceMatrix does for ReadDistanceTable. Prefixing here is what
	// turns those into SQL errors a user can act on.
	OrfCopiesResult result;
	try {
		result = ComputeOrfCopies(counts, coords, params);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s: %s", kCallerName, e.what());
	}

	EmitDiagnostics(context, result);
	gstate->values = std::move(result.values);
	return std::move(gstate);
}

void AbsQuantOrfCopiesExecute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &g = data_p.global_state->Cast<AbsQuantOrfCopiesGlobalState>();
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

void RegisterAbsQuantOrfCopies(ExtensionLoader &loader) {
	TableFunction orf_copies("absquant_orf_copies", {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR},
	                         AbsQuantOrfCopiesExecute, AbsQuantOrfCopiesBind, AbsQuantOrfCopiesInitGlobal);
	orf_copies.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(orf_copies);
}

} // namespace duckdb
