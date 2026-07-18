#include "unifrac_table_functions.hpp"

#include <algorithm>
#include <climits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "NewickTree.hpp"
#include "id_column_utils.hpp"
#include "tree_table_reader.hpp"
#include "unifrac_bptree.hpp"
#include "unifrac_distance.hpp"
#include "unifrac_function_common.hpp"
#include "unifrac_omp_scope.hpp"
#include "unifrac_support_biom.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"

// scikit-bio-binaries — randomized PCoA on a libssu fp32 distance matrix.
#include "ordination.h"

namespace duckdb {
namespace {

using unifrac_internal::AcceptedVariantList;
using unifrac_internal::IsValidVariant;
using unifrac_internal::ReadDistanceTable;
using unifrac_internal::ReadFeatureTable;
using unifrac_internal::ResolveSampleIdOutputType;
using unifrac_internal::ResolveThreadsParameter;

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
	// Output type for sample_id — mirrors the input sample_id type (BIGINT/UUID)
	// or VARCHAR otherwise. See ResolveSampleIdOutputType.
	LogicalType sample_id_type = LogicalType::VARCHAR;
};

struct UnifracPcoaGlobalState : public GlobalTableFunctionState {
	std::vector<PcoaRow> rows;
	size_t cursor = 0;
	LogicalType sample_id_type = LogicalType::VARCHAR;
	idx_t MaxThreads() const override {
		return 1;
	}
};

std::vector<std::string> CollectIds(char **ids, int n) {
	std::vector<std::string> out;
	out.reserve(n);
	for (int i = 0; i < n; ++i) {
		out.emplace_back(ids[i]);
	}
	return out;
}

// Run randomized PCoA (skbb_pcoa_fsvd_fp32) on a dense fp32 distance matrix and
// append one PcoaRow per (sample, axis). Shared by unifrac_pcoa (once per
// subsample iteration, over a freshly computed UniFrac matrix) and the
// metric-agnostic `pcoa` (a single iteration over any distance table), so the
// ordination + row-emission path is byte-for-byte identical for both — the
// property the pcoa/unifrac_pcoa equivalence test pins.
//
// `mat` is n×n row-major fp32; `ids` are the matrix's sample ids (size n).
// `caller_name` prefixes the too-few-samples error. Throws InvalidInputException
// when n < n_dims + 1 (PCoA loses one dimension to centering).
void RunPcoaOnMatrix(const float *mat, uint32_t n, const std::vector<std::string> &ids, uint32_t n_dims, int seed,
                     int n_threads, int32_t iteration_index, const char *caller_name, std::vector<PcoaRow> &out_rows) {
	if (n < n_dims + 1) {
		throw InvalidInputException("%s: only %u sample(s) available for ordination; n_dims=%u requires at least %u "
		                            "samples (PCoA loses one dimension to centering)",
		                            caller_name, n, n_dims, n_dims + 1);
	}

	std::vector<float> eigvals(n_dims);
	std::vector<float> samples(static_cast<size_t>(n) * n_dims);
	std::vector<float> prop(n_dims);

	// skbb_pcoa_fsvd_fp32 uses its own per-call seed for randomization, but its
	// `#pragma omp parallel for` regions (principal_coordinate_analysis.cpp)
	// still need the process-wide OpenMP serialization so concurrent queries
	// don't race on omp_set_num_threads.
	{
		miint::unifrac::OmpThreadScope omp_scope(n_threads);
		skbb_pcoa_fsvd_fp32(n, mat, n_dims, seed, eigvals.data(), samples.data(), prop.data());
	}

	// samples is laid out (n × n_dims), sample-major. The header comment in
	// ordination.h reads "(n_eighs × n_dims)" but the actual implementation
	// (principal_coordinate_analysis.cpp:574-578 and the preceding
	// transpose_T(n_dims, n_eighs, ...)) writes `samples + row * n_eighs` with
	// row iterating over samples, i.e. sample-major.
	for (uint32_t s = 0; s < n; ++s) {
		for (uint32_t axis = 0; axis < n_dims; ++axis) {
			PcoaRow row;
			row.iteration = iteration_index;
			row.sample_id = ids[s];
			row.axis = static_cast<int32_t>(axis);
			row.coordinate = static_cast<double>(samples[s * n_dims + axis]);
			row.eigenvalue = static_cast<double>(eigvals[axis]);
			row.proportion_explained = static_cast<double>(prop[axis]);
			out_rows.push_back(std::move(row));
		}
	}
}

// Declare the PCoA output schema. Shared by unifrac_pcoa and pcoa so the two
// functions can never drift apart column-wise — the "identical output schema"
// invariant is enforced structurally rather than by discipline. `sample_id_type`
// is the mirrored input id type (see ResolveSampleIdOutputType).
void DeclarePcoaOutputSchema(const LogicalType &sample_id_type, vector<LogicalType> &return_types,
                             vector<string> &names) {
	names.emplace_back("iteration");
	return_types.emplace_back(LogicalType::INTEGER);
	names.emplace_back("sample_id");
	return_types.emplace_back(sample_id_type);
	names.emplace_back("axis");
	return_types.emplace_back(LogicalType::INTEGER);
	names.emplace_back("coordinate");
	return_types.emplace_back(LogicalType::DOUBLE);
	names.emplace_back("eigenvalue");
	return_types.emplace_back(LogicalType::DOUBLE);
	names.emplace_back("proportion_explained");
	return_types.emplace_back(LogicalType::DOUBLE);
}

void ComputeOneIteration(const miint::unifrac::UnifracSupportBiomView &biom_view,
                         const miint::unifrac::UnifracBptreeView &bptree_view, const std::string &variant_fp32,
                         bool variance_adjust, double alpha, bool bypass_tips, bool normalize_sample_counts,
                         uint32_t subsample_depth, bool subsample_with_replacement, int seed_iter, uint32_t n_dims,
                         int n_threads, int32_t iteration_index, std::vector<PcoaRow> &out_rows) {
	// UnifracDistanceMatrix::Compute throws std::runtime_error on libssu
	// errors (unknown_method, table_empty, table_and_tree_do_not_overlap,
	// output_error). The tree/feature mismatch case is already caught
	// upstream by ValidateTreeCoversFeatures, so anything reaching here is
	// a libssu-internal failure we surface as InvalidInputException.
	miint::unifrac::UnifracDistanceMatrix dist = [&]() {
		try {
			return miint::unifrac::UnifracDistanceMatrix::Compute(
			    biom_view, bptree_view, variant_fp32, variance_adjust, alpha, bypass_tips, normalize_sample_counts,
			    subsample_depth, subsample_with_replacement, seed_iter, n_threads);
		} catch (const std::runtime_error &e) {
			throw InvalidInputException("unifrac_pcoa: %s", e.what());
		}
	}();

	// Subsampling can drop samples whose total counts fall below
	// subsample_depth, so the distance matrix's n_samples may be smaller than
	// the input feature-table's n_samples. RunPcoaOnMatrix guards n vs n_dims
	// (the same check the bind-time validation applies to the pre-subsample
	// count) and emits the rows.
	RunPcoaOnMatrix(dist.matrix(), dist.n_samples(), dist.sample_ids(), n_dims, seed_iter, n_threads, iteration_index,
	                "unifrac_pcoa", out_rows);
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

	// All INTEGER-typed named parameters are read as int32_t to match
	// LogicalType::INTEGER (see align_common.hpp:75 for the project
	// convention). Promote to int64_t at use sites that need a wider range.
	std::string variant = "weighted_normalized";
	int32_t n_dims = 3;
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
		} else if (key == "n_dims") {
			n_dims = kv.second.GetValue<int32_t>();
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
	const int n_threads = ResolveThreadsParameter(context, threads, "unifrac_pcoa");

	if (!IsValidVariant(variant)) {
		throw BinderException("unifrac_pcoa: variant '%s' is not valid. Accepted variants: %s", variant,
		                      AcceptedVariantList());
	}
	if (n_dims < 1) {
		throw BinderException("unifrac_pcoa: n_dims must be >= 1 (got %d)", n_dims);
	}
	if (n_subsamples < 1) {
		throw BinderException("unifrac_pcoa: n_subsamples must be >= 1 (got %d)", n_subsamples);
	}
	if (subsample_depth < 0) {
		throw BinderException("unifrac_pcoa: subsample_depth must be >= 0 (got %d)", subsample_depth);
	}
	if (n_subsamples > 1 && subsample_depth == 0) {
		throw BinderException(
		    "unifrac_pcoa: n_subsamples > 1 requires subsample_depth > 0 (iterations would otherwise be identical)");
	}
	// Guard the seed_iter = seed + i arithmetic against signed int32 overflow.
	// Promotes to int64_t for the check; libssu's ssu_set_random_seed and
	// skbb_pcoa_fsvd_fp32 both take 32-bit ints, so values past INT32_MAX
	// can't be honoured anyway.
	if (seed >= 0 &&
	    static_cast<int64_t>(seed) + static_cast<int64_t>(n_subsamples) - 1 > static_cast<int64_t>(INT_MAX)) {
		throw BinderException(
		    "unifrac_pcoa: seed (%d) + n_subsamples (%d) - 1 exceeds INT_MAX; pick a smaller seed or fewer subsamples",
		    seed, n_subsamples);
	}

	LogicalType sample_id_col_type = LogicalType::VARCHAR;
	auto coo_rows = ReadFeatureTable(context, table_name, "unifrac_pcoa", &sample_id_col_type);
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
	if (static_cast<uint32_t>(n_dims) > n_samples - 1) {
		throw BinderException(
		    "unifrac_pcoa: n_dims (%d) must be <= n_samples - 1 (%u). PCoA loses one dimension to centering.", n_dims,
		    n_samples - 1);
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
	data->sample_id_type = ResolveSampleIdOutputType(sample_id_col_type);
	const auto rows_per_iter = static_cast<size_t>(n_samples) * static_cast<size_t>(n_dims);
	data->rows.reserve(static_cast<size_t>(n_subsamples) * rows_per_iter);
	for (int32_t i = 0; i < n_subsamples; ++i) {
		// seed + i overflow is prevented by the bind-time check above.
		const int seed_iter = (seed >= 0) ? (seed + i) : -1;
		ComputeOneIteration(biom_view, bptree_view, variant_fp32, variance_adjust, alpha, bypass_tips,
		                    normalize_sample_counts, static_cast<uint32_t>(subsample_depth), subsample_with_replacement,
		                    seed_iter, static_cast<uint32_t>(n_dims), n_threads, i, data->rows);
	}

	DeclarePcoaOutputSchema(data->sample_id_type, return_types, names);

	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> UnifracPcoaInitGlobal(ClientContext &, TableFunctionInitInput &input) {
	auto &data = input.bind_data->CastNoConst<UnifracPcoaData>();
	auto gstate = make_uniq<UnifracPcoaGlobalState>();
	gstate->rows = std::move(data.rows);
	gstate->sample_id_type = data.sample_id_type;
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
	auto axis_data = FlatVector::GetData<int32_t>(output.data[2]);
	auto coord_data = FlatVector::GetData<double>(output.data[3]);
	auto eig_data = FlatVector::GetData<double>(output.data[4]);
	auto pe_data = FlatVector::GetData<double>(output.data[5]);

	for (idx_t i = 0; i < n; ++i) {
		const auto &r = gstate.rows[gstate.cursor + i];
		iter_data[i] = r.iteration;
		// EmitIdCell mirrors the id type; its ""/"*"→NULL sentinel branch is
		// unreachable here (ReadFeatureTable drops NULL sample_ids).
		EmitIdCell(sample_id_vec, i, r.sample_id, gstate.sample_id_type);
		axis_data[i] = r.axis;
		coord_data[i] = r.coordinate;
		eig_data[i] = r.eigenvalue;
		pe_data[i] = r.proportion_explained;
	}

	gstate.cursor += n;
	output.SetCardinality(n);
}

// ── pcoa(distances, ...) — metric-agnostic PCoA over a condensed distance table ──
// Decoupled from UniFrac: reads any `(sample_a, sample_b, distance)` relation
// (the unifrac_distances output, a beta_* macro result, or a precomputed
// Bray-Curtis/Jaccard/Euclidean table) into a dense matrix and runs the exact
// same skbb_pcoa_fsvd_fp32 + emit path as unifrac_pcoa via RunPcoaOnMatrix. It
// reuses UnifracPcoaData / UnifracPcoaGlobalState / UnifracPcoaExecute /
// UnifracPcoaInitGlobal unchanged — the output schema is identical, iteration is
// always 0 (kept for schema parity), and there is no subsampling (a distance
// table is a fixed matrix).
unique_ptr<FunctionData> PcoaFromDistancesBind(ClientContext &context, TableFunctionBindInput &input,
                                               vector<LogicalType> &return_types, vector<string> &names) {
	const std::string table_name = input.inputs[0].GetValue<string>();
	if (table_name.empty()) {
		throw BinderException("pcoa: distance-table name must not be empty");
	}

	int32_t n_dims = 3;
	int32_t seed = -1;
	int32_t threads = 0; // 0 = follow DuckDB's TaskScheduler::NumberOfThreads()
	for (const auto &kv : input.named_parameters) {
		const auto key = StringUtil::Lower(kv.first);
		if (key == "n_dims") {
			n_dims = kv.second.GetValue<int32_t>();
		} else if (key == "seed") {
			seed = kv.second.GetValue<int32_t>();
		} else if (key == "threads") {
			threads = kv.second.GetValue<int32_t>();
		}
	}
	const int n_threads = ResolveThreadsParameter(context, threads, "pcoa");
	if (n_dims < 1) {
		throw BinderException("pcoa: n_dims must be >= 1 (got %d)", n_dims);
	}

	auto dist = ReadDistanceTable(context, table_name, "pcoa");
	if (static_cast<uint32_t>(n_dims) > dist.n_samples - 1) {
		throw BinderException("pcoa: n_dims (%d) must be <= n_samples - 1 (%u). PCoA loses one dimension to centering.",
		                      n_dims, dist.n_samples - 1);
	}

	auto data = make_uniq<UnifracPcoaData>();
	data->sample_id_type = dist.sample_id_type;
	data->rows.reserve(static_cast<size_t>(dist.n_samples) * static_cast<size_t>(n_dims));
	RunPcoaOnMatrix(dist.matrix.data(), dist.n_samples, dist.sample_ids, static_cast<uint32_t>(n_dims), seed, n_threads,
	                /*iteration_index*/ 0, "pcoa", data->rows);

	DeclarePcoaOutputSchema(data->sample_id_type, return_types, names);

	return std::move(data);
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
	fn.named_parameters["threads"] = LogicalType::INTEGER;
	loader.RegisterFunction(fn);
}

void RegisterPcoaFromDistances(ExtensionLoader &loader) {
	TableFunction fn("pcoa", {LogicalType::VARCHAR}, UnifracPcoaExecute, PcoaFromDistancesBind, UnifracPcoaInitGlobal);
	fn.named_parameters["n_dims"] = LogicalType::INTEGER;
	fn.named_parameters["seed"] = LogicalType::INTEGER;
	fn.named_parameters["threads"] = LogicalType::INTEGER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
