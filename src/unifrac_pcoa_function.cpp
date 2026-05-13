#include "unifrac_table_functions.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "NewickTree.hpp"
#include "tree_table_reader.hpp"
#include "unifrac_bptree.hpp"
#include "unifrac_distance.hpp"
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

// scikit-bio-binaries — randomized PCoA on a libssu fp32 distance matrix.
#include "ordination.h"

namespace duckdb {

namespace {

constexpr std::array<const char *, 5> kAcceptedVariants = {"unweighted", "weighted_normalized", "weighted_unnormalized",
                                                           "unweighted_unnormalized", "generalized"};

bool IsValidVariant(const std::string &v) {
	for (const auto *name : kAcceptedVariants) {
		if (v == name) {
			return true;
		}
	}
	return false;
}

std::string AcceptedVariantList() {
	std::string out;
	for (size_t i = 0; i < kAcceptedVariants.size(); ++i) {
		if (i != 0) {
			out += ", ";
		}
		out += kAcceptedVariants[i];
	}
	return out;
}

struct PcoaRow {
	int32_t iteration;
	std::string sample_id;
	int32_t axis;
	double coordinate;
	double eigenvalue;
	double proportion_explained;
};

struct UnifracPcoaData : public TableFunctionData {
	std::vector<PcoaRow> rows;
};

struct UnifracPcoaGlobalState : public GlobalTableFunctionState {
	std::vector<PcoaRow> rows;
	size_t cursor = 0;
	idx_t MaxThreads() const override {
		return 1;
	}
};

// Read the user-named feature-table relation (must expose
// `(sample_id VARCHAR, feature_id VARCHAR, value DOUBLE)`, matching the
// columns produced by read_biom) into long-form COO rows.
std::vector<miint::unifrac::CooRow> ReadFeatureTable(ClientContext &context, const std::string &table_name) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);

	// Schema probe via LIMIT 0 — surfaces missing columns or unsafe casts as a
	// binder-time error before we materialize the full table.
	auto probe = conn.Query("SELECT sample_id::VARCHAR, feature_id::VARCHAR, value::DOUBLE FROM " + qname + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException(
		    "unifrac_pcoa: feature-table '%s' must expose (sample_id VARCHAR, feature_id VARCHAR, value DOUBLE): %s",
		    table_name, probe->GetError());
	}

	auto result = conn.Query("SELECT sample_id::VARCHAR, feature_id::VARCHAR, value::DOUBLE FROM " + qname);
	if (result->HasError()) {
		throw InvalidInputException("unifrac_pcoa: failed to read feature-table '%s': %s", table_name,
		                            result->GetError());
	}

	std::vector<miint::unifrac::CooRow> rows;
	auto &materialized = result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		const idx_t n = chunk->size();
		if (n == 0) {
			break;
		}
		UnifiedVectorFormat sid_u, fid_u, val_u;
		chunk->data[0].ToUnifiedFormat(n, sid_u);
		chunk->data[1].ToUnifiedFormat(n, fid_u);
		chunk->data[2].ToUnifiedFormat(n, val_u);
		auto sid_data = UnifiedVectorFormat::GetData<string_t>(sid_u);
		auto fid_data = UnifiedVectorFormat::GetData<string_t>(fid_u);
		auto val_data = UnifiedVectorFormat::GetData<double>(val_u);
		for (idx_t i = 0; i < n; ++i) {
			const auto si = sid_u.sel->get_index(i);
			const auto fi = fid_u.sel->get_index(i);
			const auto vi = val_u.sel->get_index(i);
			if (!sid_u.validity.RowIsValid(si) || !fid_u.validity.RowIsValid(fi) || !val_u.validity.RowIsValid(vi)) {
				continue;
			}
			const double v = val_data[vi];
			if (v == 0.0 || std::isnan(v)) {
				continue; // sparse-storage invariant; UnifracSupportBiomView would drop these anyway
			}
			rows.push_back({sid_data[si].GetString(), fid_data[fi].GetString(), v});
		}
	}
	return rows;
}

std::vector<std::string> CollectIds(char **ids, int n) {
	std::vector<std::string> out;
	out.reserve(n);
	for (int i = 0; i < n; ++i) {
		out.emplace_back(ids[i]);
	}
	return out;
}

void ComputeOneIteration(const miint::unifrac::UnifracSupportBiomView &biom_view,
                         const miint::unifrac::UnifracBptreeView &bptree_view, const std::string &variant_fp32,
                         bool variance_adjust, double alpha, bool bypass_tips, bool normalize_sample_counts,
                         uint32_t subsample_depth, bool subsample_with_replacement, int seed_iter, uint32_t n_dims,
                         int32_t iteration_index, std::vector<PcoaRow> &out_rows) {
	miint::unifrac::UnifracDistanceMatrix dist = [&]() {
		try {
			return miint::unifrac::UnifracDistanceMatrix::Compute(
			    biom_view, bptree_view, variant_fp32, variance_adjust, alpha, bypass_tips, normalize_sample_counts,
			    subsample_depth, subsample_with_replacement, seed_iter);
		} catch (const std::runtime_error &e) {
			throw InvalidInputException("unifrac_pcoa: %s", e.what());
		}
	}();

	// Subsampling can drop samples whose total counts fall below
	// subsample_depth, so the distance matrix's n_samples may be smaller
	// than the input feature-table's n_samples; n_dims must still fit.
	const uint32_t actual_n_samples = dist.n_samples();
	if (actual_n_samples < n_dims + 1) {
		throw InvalidInputException(
		    "unifrac_pcoa: after subsampling iteration %d only %u sample(s) survive "
		    "(samples whose total count falls below subsample_depth=%u are dropped); n_dims=%u requires >= %u samples",
		    iteration_index, actual_n_samples, subsample_depth, n_dims, n_dims + 1);
	}

	std::vector<float> eigvals(n_dims);
	std::vector<float> samples(static_cast<size_t>(actual_n_samples) * n_dims);
	std::vector<float> prop(n_dims);

	// skbb_pcoa_fsvd_fp32 uses its own per-call seed; no libssu mutex needed.
	skbb_pcoa_fsvd_fp32(actual_n_samples, dist.matrix(), n_dims, seed_iter, eigvals.data(), samples.data(),
	                    prop.data());

	// samples is laid out (actual_n_samples × n_dims), sample-major. The header
	// comment in ordination.h reads "(n_eighs × n_dims)" but the actual
	// implementation (principal_coordinate_analysis.cpp:574-578 and the
	// preceding transpose_T(n_dims, n_eighs, ...)) writes `samples + row *
	// n_eighs` with row iterating over samples, i.e. sample-major.
	const auto &sample_ids = dist.sample_ids();
	for (uint32_t s = 0; s < actual_n_samples; ++s) {
		for (uint32_t axis = 0; axis < n_dims; ++axis) {
			PcoaRow row;
			row.iteration = iteration_index;
			row.sample_id = sample_ids[s];
			row.axis = static_cast<int32_t>(axis);
			row.coordinate = static_cast<double>(samples[s * n_dims + axis]);
			row.eigenvalue = static_cast<double>(eigvals[axis]);
			row.proportion_explained = static_cast<double>(prop[axis]);
			out_rows.push_back(std::move(row));
		}
	}
}

unique_ptr<FunctionData> UnifracPcoaBind(ClientContext &context, TableFunctionBindInput &input,
                                         vector<LogicalType> &return_types, vector<string> &names) {
	const std::string table_name = input.inputs[0].GetValue<string>();
	const std::string tree_name = input.inputs[1].GetValue<string>();
	if (table_name.empty()) {
		throw BinderException("unifrac_pcoa: feature-table name must not be empty");
	}
	if (tree_name.empty()) {
		throw BinderException("unifrac_pcoa: tree name must not be empty");
	}

	std::string variant = "weighted_normalized";
	int64_t n_dims = 3;
	bool variance_adjust = false;
	double alpha = 1.0;
	bool bypass_tips = false;
	bool normalize_sample_counts = true;
	int64_t subsample_depth = 0;
	bool subsample_with_replacement = false;
	int64_t n_subsamples = 1;
	int64_t seed = -1;
	for (const auto &kv : input.named_parameters) {
		const auto key = StringUtil::Lower(kv.first);
		if (key == "variant") {
			variant = kv.second.GetValue<string>();
		} else if (key == "n_dims") {
			n_dims = kv.second.GetValue<int64_t>();
		} else if (key == "variance_adjust") {
			variance_adjust = kv.second.GetValue<bool>();
		} else if (key == "alpha") {
			alpha = kv.second.GetValue<double>();
		} else if (key == "bypass_tips") {
			bypass_tips = kv.second.GetValue<bool>();
		} else if (key == "normalize_sample_counts") {
			normalize_sample_counts = kv.second.GetValue<bool>();
		} else if (key == "subsample_depth") {
			subsample_depth = kv.second.GetValue<int64_t>();
		} else if (key == "subsample_with_replacement") {
			subsample_with_replacement = kv.second.GetValue<bool>();
		} else if (key == "n_subsamples") {
			n_subsamples = kv.second.GetValue<int64_t>();
		} else if (key == "seed") {
			seed = kv.second.GetValue<int64_t>();
		}
	}

	if (!IsValidVariant(variant)) {
		throw BinderException("unifrac_pcoa: variant '%s' is not valid. Accepted variants: %s", variant,
		                      AcceptedVariantList());
	}
	if (n_dims < 1) {
		throw BinderException("unifrac_pcoa: n_dims must be >= 1 (got %lld)", static_cast<long long>(n_dims));
	}
	if (n_subsamples < 1) {
		throw BinderException("unifrac_pcoa: n_subsamples must be >= 1 (got %lld)",
		                      static_cast<long long>(n_subsamples));
	}
	if (subsample_depth < 0) {
		throw BinderException("unifrac_pcoa: subsample_depth must be >= 0 (got %lld)",
		                      static_cast<long long>(subsample_depth));
	}
	if (n_subsamples > 1 && subsample_depth == 0) {
		throw BinderException(
		    "unifrac_pcoa: n_subsamples > 1 requires subsample_depth > 0 (iterations would otherwise be identical)");
	}

	auto coo_rows = ReadFeatureTable(context, table_name);
	if (coo_rows.empty()) {
		throw InvalidInputException("unifrac_pcoa: feature-table '%s' is empty after dropping NULL/zero rows",
		                            table_name);
	}
	miint::unifrac::UnifracSupportBiomView biom_view = [&]() {
		try {
			return miint::unifrac::UnifracSupportBiomView::FromCoo(std::move(coo_rows));
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("unifrac_pcoa: %s", e.what());
		}
	}();

	const auto *biom_struct = biom_view.support_biom();
	auto ordered_sample_ids = CollectIds(biom_struct->sample_ids, biom_struct->n_samples);
	auto feature_ids = CollectIds(biom_struct->obs_ids, biom_struct->n_obs);
	const auto n_samples = static_cast<uint32_t>(ordered_sample_ids.size());

	if (n_samples < 2) {
		throw BinderException("unifrac_pcoa: feature-table '%s' has %u sample(s); at least 2 are required for PCoA",
		                      table_name, n_samples);
	}
	if (static_cast<uint64_t>(n_dims) > static_cast<uint64_t>(n_samples) - 1) {
		throw BinderException(
		    "unifrac_pcoa: n_dims (%lld) must be <= n_samples - 1 (%u). PCoA loses one dimension to centering.",
		    static_cast<long long>(n_dims), n_samples - 1);
	}

	auto tree_inputs = ReadTreeTable(context, tree_name);
	auto tree = miint::NewickTree::build(tree_inputs);
	try {
		miint::unifrac::ValidateTreeCoversFeatures(tree, feature_ids);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("unifrac_pcoa: %s", e.what());
	}
	auto bptree_view = miint::unifrac::UnifracBptreeView::FromNewickTree(tree);

	// libssu accepts both bare and `_fp32`-suffixed variant strings
	// (unifrac-binaries/src/api.cpp:60-89). We append _fp32 explicitly so the
	// caller's choice is pinned to fp32 even if libssu changes which bare
	// names default to fp32 in a future release.
	const std::string variant_fp32 = variant + "_fp32";

	auto data = make_uniq<UnifracPcoaData>();
	const auto rows_per_iter = static_cast<size_t>(n_samples) * static_cast<size_t>(n_dims);
	data->rows.reserve(static_cast<size_t>(n_subsamples) * rows_per_iter);
	for (int32_t i = 0; i < static_cast<int32_t>(n_subsamples); ++i) {
		const int seed_iter = (seed >= 0) ? static_cast<int>(seed + i) : -1;
		ComputeOneIteration(biom_view, bptree_view, variant_fp32, variance_adjust, alpha, bypass_tips,
		                    normalize_sample_counts, static_cast<uint32_t>(subsample_depth), subsample_with_replacement,
		                    seed_iter, static_cast<uint32_t>(n_dims), i, data->rows);
	}

	names.emplace_back("iteration");
	return_types.emplace_back(LogicalType::INTEGER);
	names.emplace_back("sample_id");
	return_types.emplace_back(LogicalType::VARCHAR);
	names.emplace_back("axis");
	return_types.emplace_back(LogicalType::INTEGER);
	names.emplace_back("coordinate");
	return_types.emplace_back(LogicalType::DOUBLE);
	names.emplace_back("eigenvalue");
	return_types.emplace_back(LogicalType::DOUBLE);
	names.emplace_back("proportion_explained");
	return_types.emplace_back(LogicalType::DOUBLE);

	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> UnifracPcoaInitGlobal(ClientContext &, TableFunctionInitInput &input) {
	auto &data = input.bind_data->CastNoConst<UnifracPcoaData>();
	auto gstate = make_uniq<UnifracPcoaGlobalState>();
	gstate->rows = std::move(data.rows);
	return std::move(gstate);
}

void UnifracPcoaExecute(ClientContext &, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<UnifracPcoaGlobalState>();
	const idx_t total = gstate.rows.size();
	if (gstate.cursor >= total) {
		output.SetCardinality(0);
		return;
	}
	const idx_t remaining = total - gstate.cursor;
	const idx_t n = std::min<idx_t>(STANDARD_VECTOR_SIZE, remaining);

	auto iter_data = FlatVector::GetData<int32_t>(output.data[0]);
	auto &sample_id_vec = output.data[1];
	auto sample_id_data = FlatVector::GetData<string_t>(sample_id_vec);
	auto axis_data = FlatVector::GetData<int32_t>(output.data[2]);
	auto coord_data = FlatVector::GetData<double>(output.data[3]);
	auto eig_data = FlatVector::GetData<double>(output.data[4]);
	auto pe_data = FlatVector::GetData<double>(output.data[5]);

	for (idx_t i = 0; i < n; ++i) {
		const auto &r = gstate.rows[gstate.cursor + i];
		iter_data[i] = r.iteration;
		sample_id_data[i] = StringVector::AddString(sample_id_vec, r.sample_id);
		axis_data[i] = r.axis;
		coord_data[i] = r.coordinate;
		eig_data[i] = r.eigenvalue;
		pe_data[i] = r.proportion_explained;
	}

	gstate.cursor += n;
	output.SetCardinality(n);
}

} // namespace

void RegisterUnifracPcoa(ExtensionLoader &loader) {
	TableFunction fn("unifrac_pcoa", {LogicalType::VARCHAR, LogicalType::VARCHAR}, UnifracPcoaExecute, UnifracPcoaBind,
	                 UnifracPcoaInitGlobal);
	fn.named_parameters["variant"] = LogicalType::VARCHAR;
	fn.named_parameters["n_dims"] = LogicalType::INTEGER;
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
