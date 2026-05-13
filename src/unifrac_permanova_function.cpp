#include "unifrac_table_functions.hpp"

#include <algorithm>
#include <climits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "NewickTree.hpp"
#include "tree_table_reader.hpp"
#include "unifrac_bptree.hpp"
#include "unifrac_distance.hpp"
#include "unifrac_function_common.hpp"
#include "unifrac_metadata.hpp"
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
using unifrac_internal::ReadFeatureTable;

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

struct WideMetadata {
	std::vector<std::string> column_names;         // chosen variables, in canonical order
	std::vector<miint::unifrac::MetadataRow> rows; // unpivoted long-form view
};

// Wide-form reader: metadata must have a `sample_id` column
// (case-insensitive); every other column is a variable whose values are
// cast to VARCHAR. `requested_variables` empty → use all non-sample_id
// columns in original column order. Non-empty → exact-match lookup
// (case-insensitive), preserving user-supplied order.
WideMetadata ReadWideMetadata(ClientContext &context, const std::string &table_name,
                              const std::vector<std::string> &requested_variables) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);

	auto probe = conn.Query("SELECT * FROM " + qname + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException("unifrac_permanova: failed to read metadata relation '%s': %s", table_name,
		                            probe->GetError());
	}
	auto &probe_mat = probe->Cast<MaterializedQueryResult>();
	const auto &all_names = probe_mat.names;

	idx_t sample_id_col = DConstants::INVALID_INDEX;
	std::vector<std::string> non_sample_cols;
	std::vector<idx_t> non_sample_indices;
	for (idx_t i = 0; i < all_names.size(); ++i) {
		if (StringUtil::Lower(all_names[i]) == "sample_id") {
			if (sample_id_col != DConstants::INVALID_INDEX) {
				throw BinderException("unifrac_permanova: metadata '%s' has multiple 'sample_id' columns", table_name);
			}
			sample_id_col = i;
		} else {
			non_sample_cols.push_back(all_names[i]);
			non_sample_indices.push_back(i);
		}
	}
	if (sample_id_col == DConstants::INVALID_INDEX) {
		throw BinderException("unifrac_permanova: metadata '%s' must have a 'sample_id' column", table_name);
	}

	std::vector<std::string> chosen_variables;
	std::vector<idx_t> chosen_indices;
	if (requested_variables.empty()) {
		chosen_variables = non_sample_cols;
		chosen_indices = non_sample_indices;
	} else {
		// Lookup is case-insensitive on the table's column names, but the
		// stored canonical name is the actual column name from the relation —
		// not the user-supplied spelling. The user-supplied spelling would
		// drift into both the emitted `variable` column and the
		// MetadataRow.variable field, surprising downstream consumers who
		// expect the column name they see in the metadata table.
		std::unordered_map<std::string, idx_t> lookup;
		for (size_t k = 0; k < non_sample_cols.size(); ++k) {
			lookup[StringUtil::Lower(non_sample_cols[k])] = k;
		}
		for (const auto &v : requested_variables) {
			auto it = lookup.find(StringUtil::Lower(v));
			if (it == lookup.end()) {
				throw BinderException(
				    "unifrac_permanova: variable '%s' not found in metadata '%s' (sample_id column is reserved)", v,
				    table_name);
			}
			chosen_variables.push_back(non_sample_cols[it->second]);
			chosen_indices.push_back(non_sample_indices[it->second]);
		}
	}

	std::string sql = "SELECT " + KeywordHelper::WriteOptionallyQuoted(all_names[sample_id_col]) + "::VARCHAR";
	for (auto col_idx : chosen_indices) {
		sql += ", " + KeywordHelper::WriteOptionallyQuoted(all_names[col_idx]) + "::VARCHAR";
	}
	sql += " FROM " + qname;

	auto result = conn.Query(sql);
	if (result->HasError()) {
		throw InvalidInputException(
		    "unifrac_permanova: failed to read metadata variables from '%s': %s\nGenerated SQL: %s", table_name,
		    result->GetError(), sql);
	}

	WideMetadata out;
	out.column_names = chosen_variables;
	auto &mat = result->Cast<MaterializedQueryResult>();
	while (auto chunk = mat.Fetch()) {
		const idx_t n = chunk->size();
		if (n == 0) {
			break;
		}
		UnifiedVectorFormat sid_u;
		chunk->data[0].ToUnifiedFormat(n, sid_u);
		auto sid_data = UnifiedVectorFormat::GetData<string_t>(sid_u);

		std::vector<UnifiedVectorFormat> var_u(chosen_variables.size());
		std::vector<const string_t *> var_data(chosen_variables.size());
		for (size_t v = 0; v < chosen_variables.size(); ++v) {
			chunk->data[v + 1].ToUnifiedFormat(n, var_u[v]);
			var_data[v] = UnifiedVectorFormat::GetData<string_t>(var_u[v]);
		}
		for (idx_t i = 0; i < n; ++i) {
			const auto si = sid_u.sel->get_index(i);
			if (!sid_u.validity.RowIsValid(si)) {
				continue; // skip rows with NULL sample_id — never join targets
			}
			const std::string sample = sid_data[si].GetString();
			for (size_t v = 0; v < chosen_variables.size(); ++v) {
				const auto vi = var_u[v].sel->get_index(i);
				// NULL values become "" — BuildGroupings treats them as a real
				// (empty-string) value. Users wanting to drop NULL samples
				// should filter the metadata table before passing it in.
				const std::string value =
				    var_u[v].validity.RowIsValid(vi) ? var_data[v][vi].GetString() : std::string();
				out.rows.push_back({sample, chosen_variables[v], value});
			}
		}
	}
	return out;
}

void ComputeOneIteration(const miint::unifrac::UnifracSupportBiomView &biom_view,
                         const miint::unifrac::UnifracBptreeView &bptree_view, const std::string &variant_fp32,
                         bool variance_adjust, double alpha, bool bypass_tips, bool normalize_sample_counts,
                         uint32_t subsample_depth, bool subsample_with_replacement,
                         const std::vector<miint::unifrac::MetadataRow> &metadata_rows,
                         const std::vector<std::string> &requested_variables, int seed_iter, uint32_t n_permutations,
                         int32_t iteration_index, std::vector<PermanovaRow> &out_rows) {
	// UnifracDistanceMatrix::Compute throws std::runtime_error on libssu
	// errors. The tree/feature mismatch case is already caught upstream by
	// ValidateTreeCoversFeatures, so anything reaching here is a libssu
	// internal failure we surface as InvalidInputException.
	miint::unifrac::UnifracDistanceMatrix dist = [&]() {
		try {
			return miint::unifrac::UnifracDistanceMatrix::Compute(
			    biom_view, bptree_view, variant_fp32, variance_adjust, alpha, bypass_tips, normalize_sample_counts,
			    subsample_depth, subsample_with_replacement, seed_iter);
		} catch (const std::runtime_error &e) {
			throw InvalidInputException("unifrac_permanova: %s", e.what());
		}
	}();

	std::vector<miint::unifrac::NamedGrouping> groupings;
	try {
		groupings = miint::unifrac::BuildGroupings(metadata_rows, dist.sample_ids(), requested_variables);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("unifrac_permanova: %s", e.what());
	}

	const uint32_t n_samples = dist.n_samples();
	for (const auto &grouping : groupings) {
		// skbb_permanova_fp32 takes its own per-call seed (no libssu mutex
		// needed). Using seed_iter across all variables in an iteration
		// guarantees the Rule-7 invariant: identical groupings under the
		// same distance matrix produce identical pseudo-F.
		float f_stat = 0.0f;
		float p_value = 0.0f;
		skbb_permanova_fp32(n_samples, dist.matrix(), grouping.labels.data(), n_permutations, seed_iter, &f_stat,
		                    &p_value);
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
		}
	}

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

	auto metadata = ReadWideMetadata(context, metadata_name, requested_variables);
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
		                    metadata.rows, metadata.column_names, seed_iter, static_cast<uint32_t>(n_permutations), i,
		                    data->rows);
	}

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
	loader.RegisterFunction(fn);
}

} // namespace duckdb
