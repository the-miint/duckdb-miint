#include "unifrac_table_functions.hpp"

#include <algorithm>
#include <climits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "NewickTree.hpp"
#include "catalog_utils.hpp"
#include "tree_table_reader.hpp"
#include "unifrac_bptree.hpp"
#include "unifrac_distance.hpp"
#include "unifrac_function_common.hpp"
#include "unifrac_metadata.hpp"
#include "unifrac_omp_scope.hpp"
#include "unifrac_support_biom.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"

// scikit-bio-binaries — PERMANOVA pseudo-F + p-value on a fp32 distance matrix.
#include "distance.h"

namespace duckdb {
namespace {

using unifrac_internal::AcceptedVariantList;
using unifrac_internal::IsValidVariant;
using unifrac_internal::ReadDistanceTable;
using unifrac_internal::ReadFeatureTable;
using unifrac_internal::ReadWideMetadata;
using unifrac_internal::ResolveThreadsParameter;

struct PermanovaRow {
	int32_t iteration;
	std::string variable;
	int32_t n_groups;
	double f_stat;
	double p_value;
	int32_t n_permutations;
};

struct UnifracPermanovaData : public TableFunctionData {
	std::vector<PermanovaRow> rows;
};

struct UnifracPermanovaGlobalState : public GlobalTableFunctionState {
	std::vector<PermanovaRow> rows;
	size_t cursor = 0;
	idx_t MaxThreads() const override {
		return 1;
	}
};

// Run PERMANOVA (skbb_permanova_fp32) on a dense fp32 distance matrix for every
// grouping the metadata factorizes into, appending one PermanovaRow per
// grouping. Shared by unifrac_permanova (once per subsample iteration, over a
// freshly computed UniFrac matrix) and the metric-agnostic `permanova` (a single
// iteration over any distance table), so the grouping + skbb path is identical
// for both — the property the permanova/unifrac_permanova equivalence test pins.
//
// `mat` is n×n row-major fp32; `ids` are the matrix's sample ids (size n, the
// order BuildGroupings aligns groupings to). `seed` is used across all
// variables in the call, guaranteeing the Rule-7 invariant: identical groupings
// under the same distance matrix produce identical pseudo-F. `caller_name`
// prefixes the BuildGroupings error (missing sample / unknown variable).
void RunPermanovaOnMatrix(const float *mat, uint32_t n, const std::vector<std::string> &ids,
                          const std::vector<miint::unifrac::MetadataRow> &metadata_rows,
                          const std::vector<std::string> &requested_variables, uint32_t n_permutations, int seed,
                          int n_threads, int32_t iteration_index, const char *caller_name,
                          std::vector<PermanovaRow> &out_rows) {
	std::vector<miint::unifrac::NamedGrouping> groupings;
	try {
		groupings = miint::unifrac::BuildGroupings(metadata_rows, ids, requested_variables);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s: %s", caller_name, e.what());
	}

	for (const auto &grouping : groupings) {
		// No process-wide lock: skbb's PERMANOVA permutes from a generator array
		// seeded per call, so it is re-entrant and two queries may run at once. The
		// scope pins this thread's OpenMP fan-out and guarantees the non-negative
		// seed that keeps skbb off its process-global generator (see ComputeCallScope).
		float f_stat = 0.0f;
		float p_value = 0.0f;
		{
			miint::unifrac::ComputeCallScope skbb(n_threads, seed);
			skbb_permanova_fp32(n, mat, grouping.labels.data(), n_permutations, skbb.seed(), &f_stat, &p_value);
		}
		PermanovaRow row;
		row.iteration = iteration_index;
		row.variable = grouping.variable;
		row.n_groups = static_cast<int32_t>(grouping.n_groups);
		row.f_stat = static_cast<double>(f_stat);
		row.p_value = static_cast<double>(p_value);
		row.n_permutations = static_cast<int32_t>(n_permutations);
		out_rows.push_back(std::move(row));
	}
}

// Declare the PERMANOVA output schema. Shared by unifrac_permanova and permanova
// so the two functions can never drift apart column-wise — the "identical output
// schema" invariant is enforced structurally rather than by discipline.
void DeclarePermanovaOutputSchema(vector<LogicalType> &return_types, vector<string> &names) {
	names.emplace_back("iteration");
	return_types.emplace_back(LogicalType::INTEGER);
	names.emplace_back("variable");
	return_types.emplace_back(LogicalType::VARCHAR);
	names.emplace_back("n_groups");
	return_types.emplace_back(LogicalType::INTEGER);
	names.emplace_back("f_stat");
	return_types.emplace_back(LogicalType::DOUBLE);
	names.emplace_back("p_value");
	return_types.emplace_back(LogicalType::DOUBLE);
	names.emplace_back("n_permutations");
	return_types.emplace_back(LogicalType::INTEGER);
}

void ComputeOneIteration(const miint::unifrac::UnifracSupportBiomView &biom_view,
                         const miint::unifrac::UnifracBptreeView &bptree_view, const std::string &variant_fp32,
                         bool variance_adjust, double alpha, bool bypass_tips, bool normalize_sample_counts,
                         uint32_t subsample_depth, bool subsample_with_replacement,
                         const std::vector<miint::unifrac::MetadataRow> &metadata_rows,
                         const std::vector<std::string> &requested_variables, int seed_iter, uint32_t n_permutations,
                         int n_threads, int32_t iteration_index, std::vector<PermanovaRow> &out_rows) {
	// UnifracDistanceMatrix::Compute throws std::runtime_error on libssu
	// errors. The tree/feature mismatch case is already caught upstream by
	// ValidateTreeCoversFeatures, so anything reaching here is a libssu
	// internal failure we surface as InvalidInputException.
	miint::unifrac::UnifracDistanceMatrix dist = [&]() {
		try {
			return miint::unifrac::UnifracDistanceMatrix::Compute(
			    biom_view, bptree_view, variant_fp32, variance_adjust, alpha, bypass_tips, normalize_sample_counts,
			    subsample_depth, subsample_with_replacement, seed_iter, n_threads);
		} catch (const std::runtime_error &e) {
			throw InvalidInputException("unifrac_permanova: %s", e.what());
		}
	}();

	RunPermanovaOnMatrix(dist.matrix(), dist.n_samples(), dist.sample_ids(), metadata_rows, requested_variables,
	                     n_permutations, seed_iter, n_threads, iteration_index, "unifrac_permanova", out_rows);
}

unique_ptr<FunctionData> UnifracPermanovaBind(ClientContext &context, TableFunctionBindInput &input,
                                              vector<LogicalType> &return_types, vector<string> &names) {
	const std::string table_name = input.inputs[0].GetValue<string>();
	const std::string tree_name = input.inputs[1].GetValue<string>();
	const std::string metadata_name = input.inputs[2].GetValue<string>();
	if (table_name.empty() || tree_name.empty() || metadata_name.empty()) {
		throw BinderException("unifrac_permanova: all three positional arguments (table, tree, metadata) "
		                      "must be non-empty");
	}

	// All INTEGER-typed named parameters are read as int32_t to match
	// LogicalType::INTEGER (see align_common.hpp:75 for the project
	// convention).
	std::string variant = "weighted_normalized";
	int32_t n_permutations = 999;
	std::vector<std::string> requested_variables;
	bool variance_adjust = false;
	double alpha = 1.0;
	bool bypass_tips = false;
	bool normalize_sample_counts = true;
	int32_t subsample_depth = 0;
	bool subsample_with_replacement = false;
	int32_t n_subsamples = 1;
	int32_t seed = -1;
	int32_t threads = 0; // 0 = follow DuckDB's TaskScheduler::NumberOfThreads()
	for (const auto &kv : input.named_parameters) {
		const auto key = StringUtil::Lower(kv.first);
		if (key == "variant") {
			variant = kv.second.GetValue<string>();
		} else if (key == "n_permutations") {
			n_permutations = kv.second.GetValue<int32_t>();
		} else if (key == "variables") {
			auto &list_children = ListValue::GetChildren(kv.second);
			requested_variables.reserve(list_children.size());
			for (const auto &child : list_children) {
				requested_variables.push_back(child.GetValue<string>());
			}
		} else if (key == "variance_adjust") {
			variance_adjust = kv.second.GetValue<bool>();
		} else if (key == "alpha") {
			alpha = kv.second.GetValue<double>();
		} else if (key == "bypass_tips") {
			bypass_tips = kv.second.GetValue<bool>();
		} else if (key == "normalize_sample_counts") {
			normalize_sample_counts = kv.second.GetValue<bool>();
		} else if (key == "subsample_depth") {
			subsample_depth = kv.second.GetValue<int32_t>();
		} else if (key == "subsample_with_replacement") {
			subsample_with_replacement = kv.second.GetValue<bool>();
		} else if (key == "n_subsamples") {
			n_subsamples = kv.second.GetValue<int32_t>();
		} else if (key == "seed") {
			seed = kv.second.GetValue<int32_t>();
		} else if (key == "threads") {
			threads = kv.second.GetValue<int32_t>();
		}
	}
	const int n_threads = ResolveThreadsParameter(context, threads, "unifrac_permanova");

	if (!IsValidVariant(variant)) {
		throw BinderException("unifrac_permanova: variant '%s' is not valid. Accepted variants: %s", variant,
		                      AcceptedVariantList());
	}
	if (n_permutations < 1) {
		throw BinderException("unifrac_permanova: n_permutations must be >= 1 (got %d)", n_permutations);
	}
	if (n_subsamples < 1) {
		throw BinderException("unifrac_permanova: n_subsamples must be >= 1 (got %d)", n_subsamples);
	}
	if (subsample_depth < 0) {
		throw BinderException("unifrac_permanova: subsample_depth must be >= 0 (got %d)", subsample_depth);
	}
	if (n_subsamples > 1 && subsample_depth == 0) {
		throw BinderException("unifrac_permanova: n_subsamples > 1 requires subsample_depth > 0 (iterations would "
		                      "otherwise be identical)");
	}
	// Guard the seed_iter = seed + i arithmetic against signed int32 overflow.
	if (seed >= 0 &&
	    static_cast<int64_t>(seed) + static_cast<int64_t>(n_subsamples) - 1 > static_cast<int64_t>(INT_MAX)) {
		throw BinderException("unifrac_permanova: seed (%d) + n_subsamples (%d) - 1 exceeds INT_MAX; "
		                      "pick a smaller seed or fewer subsamples",
		                      seed, n_subsamples);
	}

	auto coo_rows = ReadFeatureTable(context, table_name, "unifrac_permanova");
	if (coo_rows.empty()) {
		throw InvalidInputException("unifrac_permanova: feature-table '%s' is empty after dropping NULL/zero rows",
		                            table_name);
	}
	miint::unifrac::UnifracSupportBiomView biom_view = [&]() {
		try {
			return miint::unifrac::UnifracSupportBiomView::FromCoo(std::move(coo_rows));
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("unifrac_permanova: %s", e.what());
		}
	}();

	const auto *biom_struct = biom_view.support_biom();
	std::vector<std::string> feature_ids;
	feature_ids.reserve(biom_struct->n_obs);
	for (int i = 0; i < biom_struct->n_obs; ++i) {
		feature_ids.emplace_back(biom_struct->obs_ids[i]);
	}
	if (biom_struct->n_samples < 2) {
		throw BinderException(
		    "unifrac_permanova: feature-table '%s' has %d sample(s); at least 2 are required for PERMANOVA", table_name,
		    biom_struct->n_samples);
	}

	auto tree_inputs = ReadTreeTable(context, tree_name);
	auto tree = miint::NewickTree::build(tree_inputs);
	try {
		miint::unifrac::ValidateTreeCoversFeatures(tree, feature_ids);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("unifrac_permanova: %s", e.what());
	}
	auto bptree_view = miint::unifrac::UnifracBptreeView::FromNewickTree(tree);

	auto metadata = ReadWideMetadata(context, metadata_name, requested_variables, "unifrac_permanova");
	if (metadata.column_names.empty()) {
		throw BinderException("unifrac_permanova: metadata '%s' has no variable columns (only 'sample_id' present)",
		                      metadata_name);
	}

	const std::string variant_fp32 = variant + "_fp32";

	auto data = make_uniq<UnifracPermanovaData>();
	data->rows.reserve(static_cast<size_t>(n_subsamples) * metadata.column_names.size());
	for (int32_t i = 0; i < n_subsamples; ++i) {
		// seed + i overflow is prevented by the bind-time check above.
		const int seed_iter = (seed >= 0) ? (seed + i) : -1;
		ComputeOneIteration(biom_view, bptree_view, variant_fp32, variance_adjust, alpha, bypass_tips,
		                    normalize_sample_counts, static_cast<uint32_t>(subsample_depth), subsample_with_replacement,
		                    metadata.rows, metadata.column_names, seed_iter, static_cast<uint32_t>(n_permutations),
		                    n_threads, i, data->rows);
	}

	DeclarePermanovaOutputSchema(return_types, names);

	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> UnifracPermanovaInitGlobal(ClientContext &, TableFunctionInitInput &input) {
	auto &data = input.bind_data->CastNoConst<UnifracPermanovaData>();
	auto gstate = make_uniq<UnifracPermanovaGlobalState>();
	gstate->rows = std::move(data.rows);
	return std::move(gstate);
}

void UnifracPermanovaExecute(ClientContext &, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<UnifracPermanovaGlobalState>();
	const idx_t total = gstate.rows.size();
	if (gstate.cursor >= total) {
		output.SetCardinality(0);
		return;
	}
	const idx_t remaining = total - gstate.cursor;
	const idx_t n = std::min<idx_t>(STANDARD_VECTOR_SIZE, remaining);

	auto iter_data = FlatVector::GetData<int32_t>(output.data[0]);
	auto &variable_vec = output.data[1];
	auto variable_data = FlatVector::GetData<string_t>(variable_vec);
	auto n_groups_data = FlatVector::GetData<int32_t>(output.data[2]);
	auto f_stat_data = FlatVector::GetData<double>(output.data[3]);
	auto p_value_data = FlatVector::GetData<double>(output.data[4]);
	auto n_perm_data = FlatVector::GetData<int32_t>(output.data[5]);

	for (idx_t i = 0; i < n; ++i) {
		const auto &r = gstate.rows[gstate.cursor + i];
		iter_data[i] = r.iteration;
		variable_data[i] = StringVector::AddString(variable_vec, r.variable);
		n_groups_data[i] = r.n_groups;
		f_stat_data[i] = r.f_stat;
		p_value_data[i] = r.p_value;
		n_perm_data[i] = r.n_permutations;
	}
	gstate.cursor += n;
	output.SetCardinality(n);
}

// ── permanova(distances, metadata, ...) — metric-agnostic PERMANOVA ───────────
// Decoupled from UniFrac: reads any `(sample_a, sample_b, distance)` relation
// (the unifrac_distances output, a beta_* macro result, or a precomputed
// Bray-Curtis/Jaccard/Euclidean table) into a dense matrix and runs the exact
// same BuildGroupings + skbb_permanova_fp32 path as unifrac_permanova via
// RunPermanovaOnMatrix. Reuses UnifracPermanovaData / UnifracPermanovaGlobalState
// / UnifracPermanovaExecute / UnifracPermanovaInitGlobal and ReadWideMetadata
// unchanged — the output schema is identical, iteration is always 0 (kept for
// parity), and there is no subsampling (a distance table is a fixed matrix).
unique_ptr<FunctionData> PermanovaFromDistancesBind(ClientContext &context, TableFunctionBindInput &input,
                                                    vector<LogicalType> &return_types, vector<string> &names) {
	const std::string table_name = input.inputs[0].GetValue<string>();
	const std::string metadata_name = input.inputs[1].GetValue<string>();
	if (table_name.empty() || metadata_name.empty()) {
		throw BinderException("permanova: both positional arguments (distances, metadata) must be non-empty");
	}

	int32_t n_permutations = 999;
	std::vector<std::string> requested_variables;
	int32_t seed = -1;
	int32_t threads = 0; // 0 = follow DuckDB's TaskScheduler::NumberOfThreads()
	for (const auto &kv : input.named_parameters) {
		const auto key = StringUtil::Lower(kv.first);
		if (key == "n_permutations") {
			n_permutations = kv.second.GetValue<int32_t>();
		} else if (key == "variables") {
			auto &list_children = ListValue::GetChildren(kv.second);
			requested_variables.reserve(list_children.size());
			for (const auto &child : list_children) {
				requested_variables.push_back(child.GetValue<string>());
			}
		} else if (key == "seed") {
			seed = kv.second.GetValue<int32_t>();
		} else if (key == "threads") {
			threads = kv.second.GetValue<int32_t>();
		}
	}
	const int n_threads = ResolveThreadsParameter(context, threads, "permanova");
	if (n_permutations < 1) {
		throw BinderException("permanova: n_permutations must be >= 1 (got %d)", n_permutations);
	}

	auto dist = ReadDistanceTable(context, table_name, "permanova");
	auto metadata = ReadWideMetadata(context, metadata_name, requested_variables, "permanova");
	if (metadata.column_names.empty()) {
		throw BinderException("permanova: metadata '%s' has no variable columns (only 'sample_id' present)",
		                      metadata_name);
	}

	auto data = make_uniq<UnifracPermanovaData>();
	data->rows.reserve(metadata.column_names.size());
	RunPermanovaOnMatrix(dist.matrix.data(), dist.n_samples, dist.sample_ids, metadata.rows, metadata.column_names,
	                     static_cast<uint32_t>(n_permutations), seed, n_threads, /*iteration_index*/ 0, "permanova",
	                     data->rows);

	DeclarePermanovaOutputSchema(return_types, names);

	return std::move(data);
}

} // namespace

void RegisterUnifracPermanova(ExtensionLoader &loader) {
	TableFunction fn("unifrac_permanova", {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR},
	                 UnifracPermanovaExecute, UnifracPermanovaBind, UnifracPermanovaInitGlobal);
	fn.named_parameters["variant"] = LogicalType::VARCHAR;
	fn.named_parameters["n_permutations"] = LogicalType::INTEGER;
	fn.named_parameters["variables"] = LogicalType::LIST(LogicalType::VARCHAR);
	fn.named_parameters["variance_adjust"] = LogicalType::BOOLEAN;
	fn.named_parameters["alpha"] = LogicalType::DOUBLE;
	fn.named_parameters["bypass_tips"] = LogicalType::BOOLEAN;
	fn.named_parameters["normalize_sample_counts"] = LogicalType::BOOLEAN;
	fn.named_parameters["subsample_depth"] = LogicalType::INTEGER;
	fn.named_parameters["subsample_with_replacement"] = LogicalType::BOOLEAN;
	fn.named_parameters["n_subsamples"] = LogicalType::INTEGER;
	fn.named_parameters["seed"] = LogicalType::INTEGER;
	fn.named_parameters["threads"] = LogicalType::INTEGER;
	loader.RegisterFunction(fn);
}

void RegisterPermanovaFromDistances(ExtensionLoader &loader) {
	TableFunction fn("permanova", {LogicalType::VARCHAR, LogicalType::VARCHAR}, UnifracPermanovaExecute,
	                 PermanovaFromDistancesBind, UnifracPermanovaInitGlobal);
	fn.named_parameters["n_permutations"] = LogicalType::INTEGER;
	fn.named_parameters["variables"] = LogicalType::LIST(LogicalType::VARCHAR);
	fn.named_parameters["seed"] = LogicalType::INTEGER;
	fn.named_parameters["threads"] = LogicalType::INTEGER;
	fn.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
