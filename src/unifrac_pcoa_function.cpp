#include "unifrac_table_functions.hpp"

#include <algorithm>
#include <climits>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "NewickTree.hpp"
#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "progressive_pcoa_core.hpp"
#include "tree_table_reader.hpp"
#include "unifrac_bptree.hpp"
#include "unifrac_dense_distance.hpp"
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
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"

// scikit-bio-binaries — randomized PCoA on a libssu fp32 distance matrix.
#include "ordination.h"

namespace duckdb {
namespace {

using unifrac_internal::AcceptedVariantList;
using unifrac_internal::DistanceRelationIds;
using unifrac_internal::EnumerateDistanceIds;
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
//
// `mat` is CONSUMED: the in-place fsvd entry point centers into the caller's
// buffer instead of allocating a second N×N one (skbb's non-in-place overload
// differs only in that it constructs NewCentered<float>(n) — same pcoa_T body,
// so results are unchanged). At 25k samples that copy is 2.5 GB, and it doubled
// the ordination-phase peak. Callers must be done with `mat` afterwards; both
// current ones drop it immediately.
void RunPcoaOnMatrix(float *mat, uint32_t n, const std::vector<std::string> &ids, uint32_t n_dims, int seed,
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
		skbb_pcoa_fsvd_inplace_fp32(n, mat, n_dims, seed, eigvals.data(), samples.data(), prop.data());
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
	RunPcoaOnMatrix(dist.mutable_matrix(), dist.n_samples(), dist.sample_ids(), n_dims, seed_iter, n_threads,
	                iteration_index, "unifrac_pcoa", out_rows);
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

// ── progressive_pcoa_from_distances(distances, ...) — scalable PCoA ──────────────
// Reference-anchored progressive PCoA over a condensed (sample_a, sample_b,
// distance) relation, for sample counts where the dense N×N eigendecomposition is
// infeasible but the condensed distances still fit on disk. A seeded-random anchor
// set defines a common frame (miint::progressive::RunProgressivePcoa); the
// remaining samples stream in batches, each ordinated with the anchors and aligned
// back by partial procrustes. Output is the identical long-format schema as pcoa
// (iteration=0; eigenvalues/proportions are the anchor reference PCoA's — a
// documented caveat) so it reuses UnifracPcoaData / *Execute / *InitGlobal.
//
// The distance blocks are sourced by querying the relation per batch (bounded
// memory — one (k+a)² block plus the id set is held at a time; the dense matrix is
// never materialized). NOTE: each block issues its own scan-and-filter query, so
// anchor rows are re-read every batch; a future optimization caches the
// anchor-touching rows once and reads batch rows via contiguous ranges.

// Partition sorted ids into (anchors, remaining). Anchors: `n_anchors` chosen by a
// seeded partial Fisher-Yates shuffle (random+seed — no percentile heuristic).
// Both halves are returned in sorted order: anchors deterministic given the seed,
// remaining in id order so batches are contiguous-ish ranges (sequential reads).
struct AnchorPartition {
	std::vector<std::string> anchors;
	std::vector<std::string> remaining;
};

AnchorPartition PickAnchors(const std::vector<std::string> &sorted_ids, uint32_t n_anchors, int seed) {
	const uint32_t n = static_cast<uint32_t>(sorted_ids.size());
	std::vector<uint32_t> idx(n);
	std::iota(idx.begin(), idx.end(), 0u);
	std::mt19937_64 rng(seed >= 0 ? static_cast<uint64_t>(seed) : std::random_device {}());
	for (uint32_t i = 0; i < n_anchors && i < n; ++i) {
		std::uniform_int_distribution<uint32_t> pick(i, n - 1);
		std::swap(idx[i], idx[pick(rng)]);
	}
	std::unordered_set<uint32_t> anchor_idx(idx.begin(), idx.begin() + n_anchors);
	AnchorPartition part;
	part.anchors.reserve(n_anchors);
	part.remaining.reserve(n - n_anchors);
	for (uint32_t i = 0; i < n; ++i) {
		if (anchor_idx.count(i)) {
			part.anchors.push_back(sorted_ids[i]); // sorted_ids ascending → anchors sorted
		} else {
			part.remaining.push_back(sorted_ids[i]);
		}
	}
	return part;
}

// Partition sorted ids using an explicit, caller-supplied anchor set (in place of
// the seeded random PickAnchors). Every anchor must be a distinct sample present
// in the distance table. Both halves come out in sorted id order (anchors sorted
// for determinism; remaining in id order so batches stay contiguous-ish ranges).
// This makes the reference frame fully reproducible/controllable — required for
// apples-to-apples parity against an external oracle that fixes its own anchors.
AnchorPartition PartitionWithExplicitAnchors(const std::vector<std::string> &sorted_ids,
                                             const std::vector<std::string> &explicit_anchors) {
	std::unordered_set<std::string> id_set(sorted_ids.begin(), sorted_ids.end());
	std::unordered_set<std::string> anchor_set;
	anchor_set.reserve(explicit_anchors.size());
	for (const auto &a : explicit_anchors) {
		if (!id_set.count(a)) {
			throw InvalidInputException(
			    "progressive_pcoa_from_distances: explicit anchor '%s' is not a sample in the distance-table", a);
		}
		if (!anchor_set.insert(a).second) {
			throw BinderException("progressive_pcoa_from_distances: duplicate explicit anchor '%s'", a);
		}
	}
	AnchorPartition part;
	part.anchors.reserve(anchor_set.size());
	part.remaining.reserve(sorted_ids.size() - anchor_set.size());
	for (const auto &id : sorted_ids) {
		if (anchor_set.count(id)) {
			part.anchors.push_back(id);
		} else {
			part.remaining.push_back(id);
		}
	}
	return part;
}

// Build the dense distance block over exactly `requested` (row/col order =
// requested order) by scanning the relation for pairs with both endpoints in the
// set, then delegating the fill + validation to the same BuildDenseDistanceMatrix
// that ReadDistanceTable/pcoa use — so the block enforces the identical contract
// (negative / non-finite / nonzero-self / conflicting-duplicate / incomplete all
// throw), rather than a weaker hand-rolled reimplementation. A batch needs its
// full (anchors + batch)² block present in the relation; a gap fails loud.
miint::progressive::DistanceBlock QueryDistanceBlock(ClientContext &context, const std::string &qname,
                                                     const std::vector<std::string> &requested) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	// Filter by a literal IN-list. Value::ToSQLString escapes quotes; an id
	// containing an embedded NUL would be truncated in the SQL text and fail to
	// match — but that degrades to the completeness check below (a "missing pair"
	// error), never silent wrong output.
	std::string in_list;
	in_list.reserve(requested.size() * 8);
	for (size_t i = 0; i < requested.size(); ++i) {
		if (i) {
			in_list += ',';
		}
		in_list += Value(requested[i]).ToSQLString();
	}
	const std::string sql = "SELECT sample_a::VARCHAR, sample_b::VARCHAR, distance::DOUBLE FROM " + qname +
	                        " WHERE sample_a::VARCHAR IN (" + in_list + ") AND sample_b::VARCHAR IN (" + in_list + ")";
	auto res = conn.Query(sql);
	if (res->HasError()) {
		throw InvalidInputException("progressive_pcoa_from_distances: block query failed: %s", res->GetError());
	}

	const uint32_t m = static_cast<uint32_t>(requested.size());
	std::unordered_map<std::string, uint32_t> pos;
	pos.reserve(m);
	for (uint32_t i = 0; i < m; ++i) {
		pos.emplace(requested[i], i);
	}

	// Resolve rows to index entries (both endpoints are in the set by the WHERE, so
	// pos.at cannot miss). NULL ids / NULL / NaN distances are "not provided" and
	// dropped — a pair left unfilled surfaces as an incompleteness error below,
	// mirroring ReadDistanceTable. Self-rows are kept so BuildDenseDistanceMatrix
	// can reject a nonzero self-distance.
	std::vector<miint::unifrac::DistanceEntry> entries;
	auto &mat = res->Cast<MaterializedQueryResult>();
	while (auto chunk = mat.Fetch()) {
		const idx_t rn = chunk->size();
		if (rn == 0) {
			break;
		}
		UnifiedVectorFormat a_u, b_u, d_u;
		chunk->data[0].ToUnifiedFormat(rn, a_u);
		chunk->data[1].ToUnifiedFormat(rn, b_u);
		chunk->data[2].ToUnifiedFormat(rn, d_u);
		auto a_data = UnifiedVectorFormat::GetData<string_t>(a_u);
		auto b_data = UnifiedVectorFormat::GetData<string_t>(b_u);
		auto d_data = UnifiedVectorFormat::GetData<double>(d_u);
		for (idx_t i = 0; i < rn; ++i) {
			const auto ai = a_u.sel->get_index(i);
			const auto bi = b_u.sel->get_index(i);
			const auto di = d_u.sel->get_index(i);
			if (!a_u.validity.RowIsValid(ai) || !b_u.validity.RowIsValid(bi) || !d_u.validity.RowIsValid(di)) {
				continue;
			}
			const double dv = d_data[di];
			if (std::isnan(dv)) {
				continue;
			}
			entries.push_back({pos.at(a_data[ai].GetString()), pos.at(b_data[bi].GetString()), dv});
		}
	}

	miint::progressive::DistanceBlock block;
	block.ids = requested;
	try {
		block.matrix = miint::unifrac::BuildDenseDistanceMatrix(entries, m, requested);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("progressive_pcoa_from_distances: %s", e.what());
	}
	return block;
}

unique_ptr<FunctionData> ProgressivePcoaFromDistancesBind(ClientContext &context, TableFunctionBindInput &input,
                                                          vector<LogicalType> &return_types, vector<string> &names) {
	const std::string table_name = input.inputs[0].GetValue<string>();
	if (table_name.empty()) {
		throw BinderException("progressive_pcoa_from_distances: distance-table name must not be empty");
	}

	int32_t n_dims = 3;
	int32_t n_anchors = 100; // PoC default; validated against N below
	int32_t batch_size = 1000;
	int32_t seed = -1;
	int32_t threads = 0;             // 0 = follow DuckDB's TaskScheduler::NumberOfThreads()
	vector<string> explicit_anchors; // if non-empty: override seeded random anchor selection
	for (const auto &kv : input.named_parameters) {
		const auto key = StringUtil::Lower(kv.first);
		if (key == "n_dims") {
			n_dims = kv.second.GetValue<int32_t>();
		} else if (key == "n_anchors") {
			n_anchors = kv.second.GetValue<int32_t>();
		} else if (key == "batch_size") {
			batch_size = kv.second.GetValue<int32_t>();
		} else if (key == "seed") {
			seed = kv.second.GetValue<int32_t>();
		} else if (key == "threads") {
			threads = kv.second.GetValue<int32_t>();
		} else if (key == "anchors") {
			for (const auto &child : ListValue::GetChildren(kv.second)) {
				if (child.IsNull()) {
					throw BinderException("progressive_pcoa_from_distances: anchors list must not contain NULL");
				}
				explicit_anchors.push_back(child.ToString());
			}
		}
	}
	const int n_threads = ResolveThreadsParameter(context, threads, "progressive_pcoa_from_distances");
	if (n_dims < 1) {
		throw BinderException("progressive_pcoa_from_distances: n_dims must be >= 1 (got %d)", n_dims);
	}
	if (batch_size < 1) {
		throw BinderException("progressive_pcoa_from_distances: batch_size must be >= 1 (got %d)", batch_size);
	}

	auto ids = EnumerateDistanceIds(context, table_name, "progressive_pcoa_from_distances");
	const auto n_samples = static_cast<uint32_t>(ids.sorted_ids.size());
	if (n_samples < 2) {
		throw InvalidInputException(
		    "progressive_pcoa_from_distances: distance-table '%s' has %u distinct sample(s); at least 2 are required",
		    table_name, n_samples);
	}
	if (static_cast<uint32_t>(n_dims) > n_samples - 1) {
		throw BinderException("progressive_pcoa_from_distances: n_dims (%d) must be <= n_samples - 1 (%u). PCoA loses "
		                      "one dimension to centering.",
		                      n_dims, n_samples - 1);
	}

	AnchorPartition part;
	if (!explicit_anchors.empty()) {
		// Explicit anchor set: takes precedence over n_anchors/seed. Must be a
		// valid, distinct subset of the table's samples, at least n_dims + 1 of them.
		part = PartitionWithExplicitAnchors(ids.sorted_ids, explicit_anchors);
		if (part.anchors.size() < static_cast<size_t>(n_dims) + 1) {
			throw BinderException(
			    "progressive_pcoa_from_distances: %zu explicit anchor(s) given but n_dims + 1 (%d) are "
			    "required; the reference PCoA and the procrustes fit each need at least that many",
			    part.anchors.size(), n_dims + 1);
		}
	} else {
		if (n_anchors < 1) {
			throw BinderException("progressive_pcoa_from_distances: n_anchors must be >= 1 (got %d)", n_anchors);
		}
		if (static_cast<uint32_t>(n_anchors) > n_samples) {
			throw BinderException("progressive_pcoa_from_distances: n_anchors (%d) exceeds the %u distinct sample(s) "
			                      "in the distance-table",
			                      n_anchors, n_samples);
		}
		if (n_anchors < n_dims + 1) {
			throw BinderException(
			    "progressive_pcoa_from_distances: n_anchors (%d) must be >= n_dims + 1 (%d); the reference "
			    "PCoA and the procrustes fit each need at least that many anchors",
			    n_anchors, n_dims + 1);
		}
		part = PickAnchors(ids.sorted_ids, static_cast<uint32_t>(n_anchors), seed);
	}

	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);
	const miint::progressive::BlockProvider provider = [&context, qname](const std::vector<std::string> &requested) {
		return QueryDistanceBlock(context, qname, requested);
	};

	miint::progressive::ProgressivePcoaResult result;
	try {
		result = miint::progressive::RunProgressivePcoa(part.anchors, part.remaining, static_cast<uint32_t>(n_dims),
		                                                static_cast<uint32_t>(batch_size), seed, n_threads, provider);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("progressive_pcoa_from_distances: %s", e.what());
	}

	auto data = make_uniq<UnifracPcoaData>();
	data->sample_id_type = ids.sample_id_type;
	data->rows.reserve(result.coords.size());
	for (const auto &c : result.coords) {
		const auto axis = static_cast<uint32_t>(c.axis);
		PcoaRow row;
		row.iteration = 0; // kept for schema parity with pcoa / unifrac_pcoa
		row.sample_id = c.sample_id;
		row.axis = c.axis;
		row.coordinate = c.coordinate;
		row.eigenvalue = result.eigvals[axis];
		row.proportion_explained = result.proportion_explained[axis];
		data->rows.push_back(std::move(row));
	}

	DeclarePcoaOutputSchema(data->sample_id_type, return_types, names);
	return std::move(data);
}

// ── progressive_pcoa_from_unifrac(feature_table, tree, ...) — the true-10M path ──
// Same reference-anchored progressive PCoA as progressive_pcoa_from_distances,
// but the distance blocks are UniFrac matrices computed on the fly per batch from
// a sliced sub-biom, so the full N×N UniFrac matrix is never formed. This is
// correct because UniFrac is pairwise-local: the distance between two samples
// depends only on their own abundance vectors and the tree, never on which other
// samples share the biom — so a per-batch UniFrac block equals the corresponding
// slice of the full matrix. Reuses the B1 core + the pcoa emit path; only the
// BlockProvider differs from the _distances variant. subsample_depth is fixed at 0
// (rarefaction + progressive alignment do not compose cleanly — each batch would
// rarefy independently against a different RNG draw).

// Distinct nonzero sample ids (sorted) and distinct features of a feature table,
// plus the resolved sample-id output type. Reads only the id dictionaries (bounded
// by sample/feature counts), never the full COO — the per-batch slices are read
// lazily by the block provider. The nonzero/non-NaN filter mirrors
// ReadFeatureTable's drop rules so an enumerated sample always materializes in its
// block (a sample with no valid feature can't be ordinated).
struct FeatureTableIds {
	std::vector<std::string> sorted_sample_ids;
	std::vector<std::string> feature_ids;
	LogicalType sample_id_type = LogicalType::VARCHAR;
};

std::vector<std::string> CollectStringColumn(QueryResult &result) {
	std::vector<std::string> out;
	auto &mat = result.Cast<MaterializedQueryResult>();
	while (auto chunk = mat.Fetch()) {
		const idx_t rn = chunk->size();
		if (rn == 0) {
			break;
		}
		UnifiedVectorFormat u;
		chunk->data[0].ToUnifiedFormat(rn, u);
		auto data = UnifiedVectorFormat::GetData<string_t>(u);
		for (idx_t i = 0; i < rn; ++i) {
			const auto ii = u.sel->get_index(i);
			if (u.validity.RowIsValid(ii)) {
				out.emplace_back(data[ii].GetString());
			}
		}
	}
	return out;
}

FeatureTableIds EnumerateFeatureTableIds(ClientContext &context, const std::string &table_name,
                                         const std::string &caller) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);
	auto probe = conn.Query("SELECT sample_id::VARCHAR, feature_id::VARCHAR, value::DOUBLE FROM " + qname + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException("%s: feature-table '%s' must expose (sample_id, feature_id, value DOUBLE): %s",
		                            caller, table_name, probe->GetError());
	}

	FeatureTableIds out;
	{
		LogicalType sid_type = LogicalType::VARCHAR;
		auto cols = GetTableOrViewColumns(context, table_name, "feature-table");
		for (idx_t i = 0; i < cols.names.size(); ++i) {
			if (StringUtil::Lower(cols.names[i]) == "sample_id") {
				sid_type = cols.types[i];
				break;
			}
		}
		out.sample_id_type = ResolveSampleIdOutputType(sid_type);
	}

	// A row enters a biom only if it survives ReadFeatureTable's exact drop rules:
	// ALL of sample_id, feature_id, value non-NULL, and value nonzero/non-NaN
	// (dropping any one column). Filter over the cast columns in a subquery — a
	// numeric-string VARCHAR `value` column (the probe accepts it, and
	// ReadFeatureTable handles it via the same ::DOUBLE cast) would not bind
	// `isnan()`/`!= 0` on the raw column. An enumerated sample thus always
	// materializes in its block, so the core's block-coverage check never fires
	// for a row this filter should have dropped.
	const std::string filtered =
	    "(SELECT sample_id::VARCHAR AS sid, feature_id::VARCHAR AS fid, value::DOUBLE AS v FROM " + qname +
	    ") t WHERE t.sid IS NOT NULL AND t.fid IS NOT NULL AND t.v IS NOT NULL AND t.v != 0 AND NOT isnan(t.v)";
	auto sres = conn.Query("SELECT DISTINCT t.sid AS id FROM " + filtered + " ORDER BY id");
	if (sres->HasError()) {
		throw InvalidInputException("%s: failed to enumerate samples of feature-table '%s': %s", caller, table_name,
		                            sres->GetError());
	}
	out.sorted_sample_ids = CollectStringColumn(*sres);

	auto fres = conn.Query("SELECT DISTINCT t.fid AS id FROM " + filtered);
	if (fres->HasError()) {
		throw InvalidInputException("%s: failed to enumerate features of feature-table '%s': %s", caller, table_name,
		                            fres->GetError());
	}
	out.feature_ids = CollectStringColumn(*fres);
	return out;
}

// Read the feature rows of exactly the requested samples (all their features) into
// COO form, applying ReadFeatureTable's NULL/zero/NaN drops.
std::vector<miint::unifrac::CooRow> QueryFeatureRows(ClientContext &context, const std::string &qname,
                                                     const std::vector<std::string> &requested) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	std::string in_list;
	in_list.reserve(requested.size() * 8);
	for (size_t i = 0; i < requested.size(); ++i) {
		if (i) {
			in_list += ',';
		}
		in_list += Value(requested[i]).ToSQLString();
	}
	const std::string sql = "SELECT sample_id::VARCHAR, feature_id::VARCHAR, value::DOUBLE FROM " + qname +
	                        " WHERE sample_id::VARCHAR IN (" + in_list + ")";
	auto res = conn.Query(sql);
	if (res->HasError()) {
		throw InvalidInputException("progressive_pcoa_from_unifrac: feature slice query failed: %s", res->GetError());
	}
	std::vector<miint::unifrac::CooRow> rows;
	auto &mat = res->Cast<MaterializedQueryResult>();
	while (auto chunk = mat.Fetch()) {
		const idx_t rn = chunk->size();
		if (rn == 0) {
			break;
		}
		UnifiedVectorFormat sid_u, fid_u, val_u;
		chunk->data[0].ToUnifiedFormat(rn, sid_u);
		chunk->data[1].ToUnifiedFormat(rn, fid_u);
		chunk->data[2].ToUnifiedFormat(rn, val_u);
		auto sid_data = UnifiedVectorFormat::GetData<string_t>(sid_u);
		auto fid_data = UnifiedVectorFormat::GetData<string_t>(fid_u);
		auto val_data = UnifiedVectorFormat::GetData<double>(val_u);
		for (idx_t i = 0; i < rn; ++i) {
			const auto si = sid_u.sel->get_index(i);
			const auto fi = fid_u.sel->get_index(i);
			const auto vi = val_u.sel->get_index(i);
			if (!sid_u.validity.RowIsValid(si) || !fid_u.validity.RowIsValid(fi) || !val_u.validity.RowIsValid(vi)) {
				continue;
			}
			const double v = val_data[vi];
			if (v == 0.0 || std::isnan(v)) {
				continue;
			}
			rows.push_back({sid_data[si].GetString(), fid_data[fi].GetString(), v});
		}
	}
	return rows;
}

// Compute the UniFrac distance block over exactly the requested samples: slice the
// feature table → sub-biom → one_off UniFrac against the (full) tree. Compute()
// takes the process-wide OmpThreadScope itself, so this adds none.
miint::progressive::DistanceBlock ComputeUnifracBlock(ClientContext &context, const std::string &qname,
                                                      const std::vector<std::string> &requested,
                                                      const miint::unifrac::UnifracBptreeView &bptree_view,
                                                      const std::string &variant_fp32, bool variance_adjust,
                                                      double alpha, bool bypass_tips, bool normalize_sample_counts,
                                                      int seed, int n_threads) {
	auto rows = QueryFeatureRows(context, qname, requested);
	miint::unifrac::UnifracSupportBiomView biom_view = [&]() {
		try {
			return miint::unifrac::UnifracSupportBiomView::FromCoo(std::move(rows));
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("progressive_pcoa_from_unifrac: %s", e.what());
		}
	}();
	miint::unifrac::UnifracDistanceMatrix dist = [&]() {
		try {
			return miint::unifrac::UnifracDistanceMatrix::Compute(
			    biom_view, bptree_view, variant_fp32, variance_adjust, alpha, bypass_tips, normalize_sample_counts,
			    /*subsample_depth=*/0,
			    /*subsample_with_replacement=*/false, seed, n_threads);
		} catch (const std::runtime_error &e) {
			throw InvalidInputException("progressive_pcoa_from_unifrac: %s", e.what());
		}
	}();

	miint::progressive::DistanceBlock block;
	block.ids = dist.sample_ids();
	const uint32_t n = dist.n_samples();
	block.matrix.assign(dist.matrix(), dist.matrix() + static_cast<size_t>(n) * n);
	return block;
}

unique_ptr<FunctionData> ProgressivePcoaFromUnifracBind(ClientContext &context, TableFunctionBindInput &input,
                                                        vector<LogicalType> &return_types, vector<string> &names) {
	const std::string table_name = input.inputs[0].GetValue<string>();
	const std::string tree_name = input.inputs[1].GetValue<string>();
	if (table_name.empty()) {
		throw BinderException("progressive_pcoa_from_unifrac: feature-table name must not be empty");
	}
	if (tree_name.empty()) {
		throw BinderException("progressive_pcoa_from_unifrac: tree name must not be empty");
	}

	std::string variant = "weighted_normalized";
	int32_t n_dims = 3;
	int32_t n_anchors = 100; // PoC default; validated against N below
	int32_t batch_size = 1000;
	int32_t seed = -1;
	int32_t threads = 0; // 0 = follow DuckDB's TaskScheduler::NumberOfThreads()
	bool variance_adjust = false;
	double alpha = 1.0;
	bool bypass_tips = false;
	bool normalize_sample_counts = true;
	for (const auto &kv : input.named_parameters) {
		const auto key = StringUtil::Lower(kv.first);
		if (key == "variant") {
			variant = kv.second.GetValue<string>();
		} else if (key == "n_dims") {
			n_dims = kv.second.GetValue<int32_t>();
		} else if (key == "n_anchors") {
			n_anchors = kv.second.GetValue<int32_t>();
		} else if (key == "batch_size") {
			batch_size = kv.second.GetValue<int32_t>();
		} else if (key == "seed") {
			seed = kv.second.GetValue<int32_t>();
		} else if (key == "threads") {
			threads = kv.second.GetValue<int32_t>();
		} else if (key == "variance_adjust") {
			variance_adjust = kv.second.GetValue<bool>();
		} else if (key == "alpha") {
			alpha = kv.second.GetValue<double>();
		} else if (key == "bypass_tips") {
			bypass_tips = kv.second.GetValue<bool>();
		} else if (key == "normalize_sample_counts") {
			normalize_sample_counts = kv.second.GetValue<bool>();
		}
	}
	const int n_threads = ResolveThreadsParameter(context, threads, "progressive_pcoa_from_unifrac");
	if (!IsValidVariant(variant)) {
		throw BinderException("progressive_pcoa_from_unifrac: variant '%s' is not valid. Accepted variants: %s",
		                      variant, AcceptedVariantList());
	}
	if (n_dims < 1) {
		throw BinderException("progressive_pcoa_from_unifrac: n_dims must be >= 1 (got %d)", n_dims);
	}
	if (batch_size < 1) {
		throw BinderException("progressive_pcoa_from_unifrac: batch_size must be >= 1 (got %d)", batch_size);
	}
	if (n_anchors < 1) {
		throw BinderException("progressive_pcoa_from_unifrac: n_anchors must be >= 1 (got %d)", n_anchors);
	}

	auto ids = EnumerateFeatureTableIds(context, table_name, "progressive_pcoa_from_unifrac");
	const auto n_samples = static_cast<uint32_t>(ids.sorted_sample_ids.size());
	if (n_samples < 2) {
		throw InvalidInputException(
		    "progressive_pcoa_from_unifrac: feature-table '%s' has %u sample(s) with nonzero features; at least 2 are "
		    "required",
		    table_name, n_samples);
	}
	if (static_cast<uint32_t>(n_anchors) > n_samples) {
		throw BinderException(
		    "progressive_pcoa_from_unifrac: n_anchors (%d) exceeds the %u sample(s) in the feature-table", n_anchors,
		    n_samples);
	}
	if (static_cast<uint32_t>(n_dims) > n_samples - 1) {
		throw BinderException(
		    "progressive_pcoa_from_unifrac: n_dims (%d) must be <= n_samples - 1 (%u). PCoA loses one "
		    "dimension to centering.",
		    n_dims, n_samples - 1);
	}
	if (n_anchors < n_dims + 1) {
		throw BinderException("progressive_pcoa_from_unifrac: n_anchors (%d) must be >= n_dims + 1 (%d); the reference "
		                      "PCoA and the procrustes fit each need at least that many anchors",
		                      n_anchors, n_dims + 1);
	}

	// Build the tree once and validate it covers every feature (each batch's
	// features are a subset, so this upfront check makes them all safe).
	auto tree_inputs = ReadTreeTable(context, tree_name);
	auto tree = miint::NewickTree::build(tree_inputs);
	try {
		miint::unifrac::ValidateTreeCoversFeatures(tree, ids.feature_ids);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("progressive_pcoa_from_unifrac: %s", e.what());
	}
	auto bptree_view = miint::unifrac::UnifracBptreeView::FromNewickTree(tree);
	const std::string variant_fp32 = variant + "_fp32";

	auto part = PickAnchors(ids.sorted_sample_ids, static_cast<uint32_t>(n_anchors), seed);

	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);
	const miint::progressive::BlockProvider provider = [&](const std::vector<std::string> &requested) {
		return ComputeUnifracBlock(context, qname, requested, bptree_view, variant_fp32, variance_adjust, alpha,
		                           bypass_tips, normalize_sample_counts, seed, n_threads);
	};

	miint::progressive::ProgressivePcoaResult result;
	try {
		result = miint::progressive::RunProgressivePcoa(part.anchors, part.remaining, static_cast<uint32_t>(n_dims),
		                                                static_cast<uint32_t>(batch_size), seed, n_threads, provider);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("progressive_pcoa_from_unifrac: %s", e.what());
	}

	auto data = make_uniq<UnifracPcoaData>();
	data->sample_id_type = ids.sample_id_type;
	data->rows.reserve(result.coords.size());
	for (const auto &c : result.coords) {
		const auto axis = static_cast<uint32_t>(c.axis);
		PcoaRow row;
		row.iteration = 0; // kept for schema parity with unifrac_pcoa
		row.sample_id = c.sample_id;
		row.axis = c.axis;
		row.coordinate = c.coordinate;
		row.eigenvalue = result.eigvals[axis];
		row.proportion_explained = result.proportion_explained[axis];
		data->rows.push_back(std::move(row));
	}

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

void RegisterProgressivePcoaFromDistances(ExtensionLoader &loader) {
	TableFunction fn("progressive_pcoa_from_distances", {LogicalType::VARCHAR}, UnifracPcoaExecute,
	                 ProgressivePcoaFromDistancesBind, UnifracPcoaInitGlobal);
	fn.named_parameters["n_dims"] = LogicalType::INTEGER;
	fn.named_parameters["n_anchors"] = LogicalType::INTEGER;
	fn.named_parameters["batch_size"] = LogicalType::INTEGER;
	fn.named_parameters["seed"] = LogicalType::INTEGER;
	fn.named_parameters["threads"] = LogicalType::INTEGER;
	fn.named_parameters["anchors"] = LogicalType::LIST(LogicalType::VARCHAR);
	loader.RegisterFunction(fn);
}

void RegisterProgressivePcoaFromUnifrac(ExtensionLoader &loader) {
	TableFunction fn("progressive_pcoa_from_unifrac", {LogicalType::VARCHAR, LogicalType::VARCHAR}, UnifracPcoaExecute,
	                 ProgressivePcoaFromUnifracBind, UnifracPcoaInitGlobal);
	fn.named_parameters["variant"] = LogicalType::VARCHAR;
	fn.named_parameters["n_dims"] = LogicalType::INTEGER;
	fn.named_parameters["n_anchors"] = LogicalType::INTEGER;
	fn.named_parameters["batch_size"] = LogicalType::INTEGER;
	fn.named_parameters["seed"] = LogicalType::INTEGER;
	fn.named_parameters["threads"] = LogicalType::INTEGER;
	fn.named_parameters["variance_adjust"] = LogicalType::BOOLEAN;
	fn.named_parameters["alpha"] = LogicalType::DOUBLE;
	fn.named_parameters["bypass_tips"] = LogicalType::BOOLEAN;
	fn.named_parameters["normalize_sample_counts"] = LogicalType::BOOLEAN;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
