#include "absquant.hpp"

#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "miint_log.hpp"
#include "unifrac_function_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/main/query_result.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace duckdb {

namespace {

using miint::absquant::CountObservation;
using miint::absquant::FitOptions;
using miint::absquant::FitResult;
using miint::absquant::FitSyndnaModels;
using miint::absquant::FormatIdList;
using miint::absquant::SampleMass;
using miint::absquant::SampleModel;
using miint::absquant::SyndnaConcentration;
using unifrac_internal::ReadFeatureTable;
using unifrac_internal::ResolveSampleIdOutputType;

constexpr const char *kCallerName = "absquant_fit_models";

// A two-column reference relation read into parallel vectors.
struct KeyedValues {
	std::vector<std::string> keys;
	std::vector<double> values;
};

// Read `(key_column, value_column)` from a user-named relation, following the
// separate-connection recipe in docs/internals/reading-tables-views.md.
//
// A NULL value becomes NaN, which is how SQL NULL travels into the pure core --
// for the parameters relation that is a legitimate "skip this sample", and for
// the concentrations relation the core rejects it. A NULL KEY is refused
// outright: an unnamed row cannot be joined to anything, so silently dropping it
// would quietly shrink the pool the mass fractions are computed against.
//
// (ReadFeatureTable, which reads the counts relation, instead drops rows with a
// NULL id -- an established contract shared with unifrac_* and
// community_distances that is not this milestone's to change. The asymmetry is
// defensible on its own terms, a NULL cell in a sparse matrix being absence
// rather than a broken config, but it is worth revisiting across all of them.)
//
// `table_name` is quoted before interpolation and may come from the user.
// `key_column` and `value_column` are NOT quoted and MUST be trusted
// compile-time literals -- wiring a user-supplied column name (a named
// parameter, say) into either is a SQL injection. M3-M5 add per-sample
// parameter columns and will be tempted to reuse this; keep the column names
// baked in and select among them, never pass one through.
KeyedValues ReadKeyedValues(ClientContext &context, const std::string &table_name, const char *key_column,
                            const char *value_column, const char *entity) {
	auto conn = MakeReadOnlyHelperConnection(context);
	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);
	const std::string projection = std::string(key_column) + "::VARCHAR, " + std::string(value_column) + "::DOUBLE";

	// Schema probe first, so a missing column or an impossible cast surfaces
	// before the whole relation is materialized (mirrors ReadFeatureTable).
	auto probe = conn.Query("SELECT " + projection + " FROM " + qname + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException("%s: %s '%s' must expose (%s VARCHAR, %s DOUBLE): %s", kCallerName, entity,
		                            table_name, key_column, value_column, probe->GetError());
	}
	auto result = conn.Query("SELECT " + projection + " FROM " + qname);
	if (result->HasError()) {
		throw InvalidInputException("%s: failed to read %s '%s': %s", kCallerName, entity, table_name,
		                            result->GetError());
	}

	KeyedValues out;
	auto &materialized = result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		const idx_t n = chunk->size();
		if (n == 0) {
			break;
		}
		UnifiedVectorFormat key_u, val_u;
		chunk->data[0].ToUnifiedFormat(n, key_u);
		chunk->data[1].ToUnifiedFormat(n, val_u);
		auto key_data = UnifiedVectorFormat::GetData<string_t>(key_u);
		auto val_data = UnifiedVectorFormat::GetData<double>(val_u);
		for (idx_t i = 0; i < n; ++i) {
			const auto ki = key_u.sel->get_index(i);
			const auto vi = val_u.sel->get_index(i);
			if (!key_u.validity.RowIsValid(ki)) {
				throw InvalidInputException("%s: %s '%s' has a row with a NULL %s", kCallerName, entity, table_name,
				                            key_column);
			}
			out.keys.push_back(key_data[ki].GetString());
			out.values.push_back(val_u.validity.RowIsValid(vi) ? val_data[vi] : std::nan(""));
		}
	}
	return out;
}

struct AbsQuantFitBindData : public TableFunctionData {
	std::string counts_table;
	std::string concentrations_table;
	std::string params_table;
	FitOptions options;
	// Output type for sample_id, mirrored from the counts relation's sample_id
	// column (BIGINT/UUID preserved, else VARCHAR).
	LogicalType sample_id_type = LogicalType::VARCHAR;
};

// The whole result is computed in InitGlobal and emitted single-threaded: one
// row per sample, and a synDNA pool is a handful of spike-ins across at most a
// few thousand samples.
struct AbsQuantFitGlobalState : public GlobalTableFunctionState {
	std::vector<SampleModel> models;
	LogicalType sample_id_type = LogicalType::VARCHAR;
	idx_t cursor = 0;
	idx_t MaxThreads() const override {
		return 1;
	}
};

unique_ptr<FunctionData> AbsQuantFitBind(ClientContext &context, TableFunctionBindInput &input,
                                         vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<AbsQuantFitBindData>();
	for (idx_t i = 0; i < 4; ++i) {
		if (input.inputs[i].IsNull()) {
			throw BinderException("%s: argument %llu must not be NULL", kCallerName, static_cast<uint64_t>(i + 1));
		}
	}
	data->counts_table = input.inputs[0].GetValue<string>();
	data->concentrations_table = input.inputs[1].GetValue<string>();
	data->params_table = input.inputs[2].GetValue<string>();
	data->options.syndna_contributing_fraction = input.inputs[3].GetValue<double>();
	for (const auto &kv : input.named_parameters) {
		if (StringUtil::Lower(kv.first) == "min_syndna_counts") {
			if (kv.second.IsNull()) {
				throw BinderException("%s: min_syndna_counts must not be NULL", kCallerName);
			}
			data->options.min_syndna_counts = kv.second.GetValue<int64_t>();
		}
	}
	// The RANGE checks on both scalars are deliberately NOT repeated here. They
	// live in the pure core, where they are unit-tested against pysyndna's own
	// bounds, and InitGlobal re-throws them with this function's name attached.
	// Nothing is emitted before InitGlobal, so failing there rather than at bind
	// costs the user nothing and keeps one definition of the bounds.

	// All three relations are resolved here so that a typo in ANY of them is a
	// bind-time error, matching phylo_ancestral_parsimony's handling of its three
	// (phylo_ancestral_parsimony.cpp:129-144). Only the counts table's columns
	// are needed at bind, for the id-type mirroring below; the other two are
	// looked up purely to fail early, and their column contracts are enforced by
	// the LIMIT 0 probes in InitGlobal, whose errors carry DuckDB's own
	// candidate-binding hints and so read better than anything repeated here.
	//
	// "table", not "relation", in the entity strings: GetTableOrViewColumns
	// renders them as "<entity> or view '<name>' does not exist".
	auto cols = GetTableOrViewColumns(context, data->counts_table, "synDNA counts table");
	GetTableOrViewColumns(context, data->concentrations_table, "synDNA concentrations table");
	GetTableOrViewColumns(context, data->params_table, "sample parameters table");
	for (idx_t i = 0; i < cols.names.size(); ++i) {
		if (StringUtil::Lower(cols.names[i]) == "sample_id") {
			data->sample_id_type = ResolveSampleIdOutputType(cols.types[i]);
			break;
		}
	}

	names = {"sample_id", "slope", "intercept", "rvalue", "pvalue", "stderr", "intercept_stderr"};
	return_types = {data->sample_id_type, LogicalType::DOUBLE, LogicalType::DOUBLE, LogicalType::DOUBLE,
	                LogicalType::DOUBLE,  LogicalType::DOUBLE, LogicalType::DOUBLE};
	return std::move(data);
}

// Every sample the caller handed us that produced no model is reported. A
// missing row is the one outcome a user cannot debug from the output alone, so
// silence is not an option for any of these four.
void EmitDiagnostics(ClientContext &context, const AbsQuantFitBindData &data, const FitResult &fit) {
	if (!fit.dropped_syndna_ids.empty()) {
		miint::EmitWarning(context,
		                   std::string(kCallerName) + ": dropped " + std::to_string(fit.dropped_syndna_ids.size()) +
		                       " synDNA feature(s) with fewer than " + std::to_string(data.options.min_syndna_counts) +
		                       " total reads across all samples: " + FormatIdList(fit.dropped_syndna_ids));
	}
	if (!fit.filtered_sample_ids.empty()) {
		miint::EmitWarning(context, std::string(kCallerName) + ": dropped " +
		                                std::to_string(fit.filtered_sample_ids.size()) +
		                                " sample(s) whose mass_syndna_input_ng is NULL, negative or not finite: " +
		                                FormatIdList(fit.filtered_sample_ids));
	}
	if (!fit.unfittable_sample_ids.empty()) {
		miint::EmitWarning(context, std::string(kCallerName) + ": no model could be fit for " +
		                                std::to_string(fit.unfittable_sample_ids.size()) +
		                                " sample(s): " + FormatIdList(fit.unfittable_sample_ids));
	}
	if (!fit.samples_without_counts.empty()) {
		// "no nonzero counts", not "no rows": the reader drops zero-valued cells
		// before the core sees them, so a sample present in the counts relation
		// with a zero in every synDNA arrives here indistinguishable from one
		// that is genuinely absent. Saying "no rows" would send a user to look
		// at a table that visibly does contain the sample.
		miint::EmitWarning(context, std::string(kCallerName) + ": " +
		                                std::to_string(fit.samples_without_counts.size()) + " sample(s) in '" +
		                                data.params_table + "' contributed no nonzero synDNA counts from '" +
		                                data.counts_table + "': " + FormatIdList(fit.samples_without_counts));
	}
}

unique_ptr<GlobalTableFunctionState> AbsQuantFitInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<AbsQuantFitBindData>();
	auto gstate = make_uniq<AbsQuantFitGlobalState>();
	gstate->sample_id_type = data.sample_id_type;

	// ReadFeatureTable drops NULL and zero-valued cells, which is exactly right
	// here: pysyndna excludes zero counts before the log10 anyway, and a synDNA
	// that is zero in every sample then simply has no rows -- the sparse case the
	// concentrations-without-counts rule exists to accept. Negative values are
	// NOT dropped by it, so they reach the core and are rejected there.
	const auto rows = ReadFeatureTable(context, data.counts_table, kCallerName);
	const auto pool = ReadKeyedValues(context, data.concentrations_table, "feature_id", "syndna_indiv_ng_ul",
	                                  "synDNA concentrations relation");
	const auto params =
	    ReadKeyedValues(context, data.params_table, "sample_id", "mass_syndna_input_ng", "sample parameters relation");

	std::vector<CountObservation> counts;
	counts.reserve(rows.size());
	for (const auto &row : rows) {
		counts.push_back({row.sample_id, row.feature_id, row.count});
	}
	std::vector<SyndnaConcentration> concentrations;
	concentrations.reserve(pool.keys.size());
	for (size_t i = 0; i < pool.keys.size(); ++i) {
		concentrations.push_back({pool.keys[i], pool.values[i]});
	}
	std::vector<SampleMass> masses;
	masses.reserve(params.keys.size());
	for (size_t i = 0; i < params.keys.size(); ++i) {
		masses.push_back({params.keys[i], params.values[i]});
	}

	// The pure core carries the whole rejection layer and throws
	// std::invalid_argument with no function name attached, exactly as
	// BuildDenseDistanceMatrix does for ReadDistanceTable. Prefixing here is what
	// turns those into SQL errors a user can act on.
	FitResult fit;
	try {
		fit = FitSyndnaModels(counts, concentrations, masses, data.options);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s: %s", kCallerName, e.what());
	}

	EmitDiagnostics(context, data, fit);
	gstate->models = std::move(fit.models);
	return std::move(gstate);
}

void AbsQuantFitExecute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &g = data_p.global_state->Cast<AbsQuantFitGlobalState>();
	const idx_t total = g.models.size();
	const idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - g.cursor);

	auto &sample_id = output.data[0];
	auto slope = FlatVector::GetData<double>(output.data[1]);
	auto intercept = FlatVector::GetData<double>(output.data[2]);
	auto rvalue = FlatVector::GetData<double>(output.data[3]);
	auto pvalue = FlatVector::GetData<double>(output.data[4]);
	auto stderr_ = FlatVector::GetData<double>(output.data[5]);
	auto intercept_stderr = FlatVector::GetData<double>(output.data[6]);

	for (idx_t r = 0; r < count; ++r) {
		const auto &model = g.models[g.cursor + r];
		EmitIdCell(sample_id, r, model.sample_id, g.sample_id_type);
		slope[r] = model.fit.slope;
		intercept[r] = model.fit.intercept;
		rvalue[r] = model.fit.rvalue;
		pvalue[r] = model.fit.pvalue;
		stderr_[r] = model.fit.stderr_;
		intercept_stderr[r] = model.fit.intercept_stderr;
	}
	g.cursor += count;
	output.SetCardinality(count);
}

} // namespace

void RegisterAbsQuant(ExtensionLoader &loader) {
	TableFunction fit("absquant_fit_models",
	                  {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::DOUBLE},
	                  AbsQuantFitExecute, AbsQuantFitBind, AbsQuantFitInitGlobal);
	fit.named_parameters["min_syndna_counts"] = LogicalType::BIGINT;
	fit.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(fit);
}

} // namespace duckdb
