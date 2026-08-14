#include "unifrac_table_functions.hpp"

#include <algorithm>
#include <climits>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <cstdlib>
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

#include "duckdb/common/types/uuid.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/appender.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/storage/buffer_manager.hpp"

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
	// Progressive PCoA only (see BatchDiagnostic): which batch placed this sample
	// and that batch's anchor-overlap disparity. batch < 0 means "not placed by a
	// batch" — the anchor rows, which ARE the reference frame — and emits NULL for
	// both columns. Unused (and unemitted) by pcoa / unifrac_pcoa.
	int32_t batch = -1;
	double batch_anchor_m2 = 0.0;
};

struct UnifracPcoaData : public TableFunctionData {
	std::vector<PcoaRow> rows;
	// Output type for sample_id — mirrors the input sample_id type (BIGINT/UUID)
	// or VARCHAR otherwise. See ResolveSampleIdOutputType.
	LogicalType sample_id_type = LogicalType::VARCHAR;
	// Progressive functions append (batch, batch_anchor_m2); the dense ones don't,
	// which keeps pcoa/unifrac_pcoa's schema exactly as documented.
	bool with_batch_diagnostics = false;
};

struct UnifracPcoaGlobalState : public GlobalTableFunctionState {
	std::vector<PcoaRow> rows;
	size_t cursor = 0;
	LogicalType sample_id_type = LogicalType::VARCHAR;
	bool with_batch_diagnostics = false;
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

	// No process-wide lock: skbb's ordination is re-entrant once its randomization
	// is seeded per call, so two queries may ordinate at the same time. The scope
	// pins this thread's OpenMP fan-out and guarantees the non-negative seed that
	// keeps skbb off its process-global generator (see ComputeCallScope).
	{
		miint::unifrac::ComputeCallScope skbb(n_threads, seed);
		skbb_pcoa_fsvd_inplace_fp32(n, mat, n_dims, skbb.seed(), eigvals.data(), samples.data(), prop.data());
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
//
// `with_batch_diagnostics` APPENDS (batch, batch_anchor_m2) for the progressive
// functions, whose result is an approximation and must therefore carry its own
// quality evidence. Appending (rather than interleaving) keeps the first six
// columns positionally identical to pcoa's, so `SELECT sample_id, axis,
// coordinate` and column-name-based consumers work across all four functions.
void DeclarePcoaOutputSchema(const LogicalType &sample_id_type, vector<LogicalType> &return_types,
                             vector<string> &names, bool with_batch_diagnostics = false) {
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
	if (with_batch_diagnostics) {
		names.emplace_back("batch");
		return_types.emplace_back(LogicalType::INTEGER);
		names.emplace_back("batch_anchor_m2");
		return_types.emplace_back(LogicalType::DOUBLE);
	}
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
	gstate->with_batch_diagnostics = data.with_batch_diagnostics;
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
	int32_t *batch_data = nullptr;
	double *m2_data = nullptr;
	if (gstate.with_batch_diagnostics) {
		batch_data = FlatVector::GetData<int32_t>(output.data[6]);
		m2_data = FlatVector::GetData<double>(output.data[7]);
	}

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
		if (batch_data != nullptr) {
			// Anchor rows (batch < 0) are the reference frame itself, not a fitted
			// batch — NULL rather than a fabricated 0 disparity, which would read as
			// "perfect fit" and quietly flatter the run.
			if (r.batch < 0) {
				FlatVector::SetNull(output.data[6], i, true);
				FlatVector::SetNull(output.data[7], i, true);
			} else {
				batch_data[i] = r.batch;
				m2_data[i] = r.batch_anchor_m2;
			}
		}
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
// in the source relation. Both halves come out in sorted id order (anchors sorted
// for determinism; remaining in id order so batches stay contiguous-ish ranges).
// This makes the reference frame fully reproducible/controllable — required for
// apples-to-apples parity against an external oracle that fixes its own anchors,
// and the way a caller supplies anchors chosen by some rule of its own rather than
// by the seeded random draw.
//
// Shared by both progressive functions, so `caller_name` and `source_desc` carry
// the wording ("distance-table" vs "feature-table") into the errors.
AnchorPartition PartitionWithExplicitAnchors(const std::vector<std::string> &sorted_ids,
                                             const std::vector<std::string> &explicit_anchors, const char *caller_name,
                                             const char *source_desc) {
	std::unordered_set<std::string> id_set(sorted_ids.begin(), sorted_ids.end());
	std::unordered_set<std::string> anchor_set;
	anchor_set.reserve(explicit_anchors.size());
	for (const auto &a : explicit_anchors) {
		if (!id_set.count(a)) {
			throw InvalidInputException("%s: explicit anchor '%s' is not a sample in the %s", caller_name, a,
			                            source_desc);
		}
		if (!anchor_set.insert(a).second) {
			throw BinderException("%s: duplicate explicit anchor '%s'", caller_name, a);
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

miint::progressive::DistanceBlock QueryDistanceBlock(ClientContext &context, const std::string &qname,
                                                     const std::vector<std::string> &requested);

// Fill EVERY block of one wave from a single pass over the relation.
//
// The per-block query below costs one full scan of the relation apiece (an IN-list
// filter has no index to use), so B batches cost B scans — O(B · pairs), i.e.
// O(N³ / batch_size), which is what made this path unusable well before the dense
// path's memory ceiling. Yet a single row is needed by at most one block, with two
// exceptions, so one pass can serve a whole wave:
//
//   both endpoints anchors      → every block's anchor×anchor corner
//   one anchor, one batch id    → that id's block only
//   both ids in the SAME batch  → that block only
//   ids in two DIFFERENT batches → NO block needs it (the common case, skipped)
//
// So the scan is unfiltered — no WHERE, no 2000-element IN-lists to plan — and
// routing is two hash lookups per row. Validation is unchanged: each block is
// filled through DenseDistanceMatrixBuilder, the same builder pcoa uses, so
// negative / non-finite / nonzero-self / conflicting-duplicate / incomplete all
// still throw per block.
//
// Memory is W blocks at once (5·(batch+anchors)² bytes each), which is why the
// caller sizes the wave from its memory budget.
class WaveDistanceBlockSource {
public:
	WaveDistanceBlockSource(ClientContext &context, std::string qname, const std::vector<std::string> &anchors)
	    : context_(context), qname_(std::move(qname)), anchors_(anchors),
	      // The staging tables below are created on a connection that inherits the
	      // caller's TEMP catalog — needed so a TEMP distance relation resolves —
	      // which means they land in the USER's session rather than in a private
	      // catalog that dies with the connection. So the names must be unique per
	      // provider (two progressive runs in one session would otherwise collide on
	      // a fixed name) and each wave must drop them explicitly. Name shape follows
	      // ReadBatchByIds.
	      wave_batch_table_("_miint_wave_batch_" +
	                        StringUtil::Replace(UUID::ToString(UUID::GenerateRandomUUID()), "-", "")),
	      wave_anchor_table_("_miint_wave_anchor_" +
	                         StringUtil::Replace(UUID::ToString(UUID::GenerateRandomUUID()), "-", "")),
	      wave_batch_quoted_(KeywordHelper::WriteOptionallyQuoted(wave_batch_table_)),
	      wave_anchor_quoted_(KeywordHelper::WriteOptionallyQuoted(wave_anchor_table_)) {
	}

	// Fill the wave's blocks in one scan. `requests` are exactly the requests the
	// core will then ask for, each being (batch ids..., all anchors...).
	void Prefetch(const std::vector<std::vector<std::string>> &requests) {
		const uint32_t a = static_cast<uint32_t>(anchors_.size());
		blocks_.clear();
		by_request_.clear();
		batch_len_.assign(requests.size(), 0);
		builders_.clear();
		blocks_.reserve(requests.size());
		for (size_t k = 0; k < requests.size(); ++k) {
			if (requests[k].size() <= a) {
				throw InvalidInputException(
				    "progressive_pcoa_from_distances: internal error — wave request %llu has no non-anchor samples",
				    static_cast<unsigned long long>(k));
			}
			batch_len_[k] = static_cast<uint32_t>(requests[k].size() - a);
			miint::progressive::DistanceBlock blk;
			blk.ids = requests[k];
			blocks_.push_back(std::move(blk));
			by_request_.emplace(requests[k].front(), k);
		}
		// Builders reference blocks_[k].ids for error messages, so build them only
		// once blocks_ has stopped reallocating.
		builders_.reserve(blocks_.size());
		for (size_t k = 0; k < blocks_.size(); ++k) {
			builders_.push_back(make_uniq<miint::unifrac::DenseDistanceMatrixBuilder>(
			    static_cast<uint32_t>(blocks_[k].ids.size()), blocks_[k].ids));
		}

		// Inherits the caller's TEMP catalog so a TEMP distance relation resolves in
		// the routed queries below. The guards drop the two staging tables when this
		// wave finishes — including on the exception paths, which matter here because
		// a failed block fill throws mid-wave and the tables now live in the user's
		// session rather than dying with the connection.
		auto conn = MakeReadOnlyHelperConnection(context_);
		// Armed BEFORE the CREATEs, not after: StageWaveMaps issues two of them, and a
		// failure on the second would otherwise leak the first. DropHelperTempRelation
		// is DROP ... IF EXISTS, so guarding a table that was never created is a no-op.
		HelperTempRelation batch_guard(conn, wave_batch_quoted_);
		HelperTempRelation anchor_guard(conn, wave_anchor_quoted_);
		StageWaveMaps(conn, requests);

		// Routing happens in SQL — DuckDB resolves ids to integer cell positions
		// vectorized and across all cores, and returns only the pairs some block
		// needs. Three queries, split by which case a row falls into, because the
		// obvious single query does not work: anchors belong to EVERY block, so a map
		// holding one row per (anchor, block) makes the second join emit W² rows per
		// anchor×anchor pair before the block filter discards all but W. At W=121
		// that measured 23 GB and 92 s. Keeping anchors out of the batch map means
		// every join key below is 1:1 on at least one side, so nothing explodes.
		//
		// Case 1 — both endpoints in the SAME batch. Rows spanning two batches match
		// no block and are dropped by the join itself.
		RunRoutedQuery(conn,
		               "SELECT b1.block, b1.pos, b2.pos, d.distance::DOUBLE FROM " + qname_ + " d JOIN " +
		                   wave_batch_quoted_ + " b1 ON d.sample_a::VARCHAR = b1.id JOIN " + wave_batch_quoted_ +
		                   " b2 ON d.sample_b::VARCHAR = b2.id"
		                   " WHERE b1.block = b2.block AND d.distance IS NOT NULL AND"
		                   " NOT isnan(d.distance::DOUBLE)",
		               /*a_is_anchor_ord=*/false, /*b_is_anchor_ord=*/false);
		// Case 2/3 — one endpoint an anchor, the other a batch sample: the batch side
		// alone determines the block, so the anchor map needs no block column.
		RunRoutedQuery(conn,
		               "SELECT b.block, an.ord, b.pos, d.distance::DOUBLE FROM " + qname_ + " d JOIN " +
		                   wave_anchor_quoted_ + " an ON d.sample_a::VARCHAR = an.id JOIN " + wave_batch_quoted_ +
		                   " b ON d.sample_b::VARCHAR = b.id"
		                   " WHERE d.distance IS NOT NULL AND NOT isnan(d.distance::DOUBLE)",
		               /*a_is_anchor_ord=*/true, /*b_is_anchor_ord=*/false);
		RunRoutedQuery(conn,
		               "SELECT b.block, b.pos, an.ord, d.distance::DOUBLE FROM " + qname_ + " d JOIN " +
		                   wave_batch_quoted_ + " b ON d.sample_a::VARCHAR = b.id JOIN " + wave_anchor_quoted_ +
		                   " an ON d.sample_b::VARCHAR = an.id"
		                   " WHERE d.distance IS NOT NULL AND NOT isnan(d.distance::DOUBLE)",
		               /*a_is_anchor_ord=*/false, /*b_is_anchor_ord=*/true);
		// Case 4 — anchor×anchor. This corner is IDENTICAL in every block of every
		// wave, so it is read from the relation once for the whole run and replayed
		// into each block from memory (a²/2 triples ≈ 8 MB at 1000 anchors).
		LoadAnchorCornerOnce(conn);
		for (size_t k = 0; k < builders_.size(); ++k) {
			const uint32_t off = batch_len_[k];
			for (const auto &t : anchor_corner_) {
				builders_[k]->Add(off + t.oa, off + t.ob, t.distance);
			}
		}

		try {
			for (size_t k = 0; k < blocks_.size(); ++k) {
				blocks_[k].matrix = builders_[k]->Finish();
			}
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("progressive_pcoa_from_distances: %s", e.what());
		}
		builders_.clear();
	}

	// Serve a request: from the wave cache when it was announced (every batch
	// block), else by a direct query (the one-off anchors-only reference block).
	//
	// Called CONCURRENTLY once the core runs a wave's batches in parallel, which is
	// safe because everything it touches was fixed by Prefetch and is only read
	// here. The direct-query fallback is not concurrent in practice: batch ids are
	// disjoint, so every announced request is found (front id + full id-vector
	// match), and only the anchors-only block — fetched before any fan-out — misses.
	miint::progressive::DistanceBlock Get(const std::vector<std::string> &requested) const {
		if (!requested.empty()) {
			auto it = by_request_.find(requested.front());
			if (it != by_request_.end() && blocks_[it->second].ids == requested) {
				return blocks_[it->second];
			}
		}
		return QueryDistanceBlock(context_, qname_, requested);
	}

private:
	// One anchor×anchor distance, addressed by ordinal in the shared anchor order.
	struct AnchorPair {
		uint32_t oa;
		uint32_t ob;
		double distance;
	};

	// Stage the wave's id→position maps as temp tables (the same pattern
	// sequence_table_reader.cpp uses for _batch_ids, uniquely named and explicitly
	// dropped for the same reason — see the constructor). The batch map holds one
	// row per non-anchor sample of the wave; the anchor map holds one row per
	// anchor, with NO block column — that is what keeps the joins from exploding.
	void StageWaveMaps(Connection &conn, const std::vector<std::vector<std::string>> &requests) {
		auto create = conn.Query("CREATE TEMPORARY TABLE " + wave_batch_quoted_ +
		                         " (id VARCHAR, block INTEGER, pos INTEGER);"
		                         "CREATE TEMPORARY TABLE " +
		                         wave_anchor_quoted_ + " (id VARCHAR, ord INTEGER)");
		if (create->HasError()) {
			throw InvalidInputException("progressive_pcoa_from_distances: failed to stage the wave map: %s",
			                            create->GetError());
		}
		{
			Appender appender(conn, wave_batch_table_);
			for (size_t k = 0; k < requests.size(); ++k) {
				for (uint32_t p = 0; p < batch_len_[k]; ++p) {
					appender.AppendRow(Value(requests[k][p]), Value::INTEGER(static_cast<int32_t>(k)),
					                   Value::INTEGER(static_cast<int32_t>(p)));
				}
			}
			appender.Close();
		}
		{
			Appender appender(conn, wave_anchor_table_);
			for (uint32_t j = 0; j < anchors_.size(); ++j) {
				appender.AppendRow(Value(anchors_[j]), Value::INTEGER(static_cast<int32_t>(j)));
			}
			appender.Close();
		}
	}

	// Run one wave query to completion and hand back its rows.
	//
	// MATERIALIZED, not streamed, and this is the single most important line in the
	// class. A streaming result advances only when the fetching thread executes a
	// task (duckdb executor.cpp: `Executor::ExecuteTask` runs it inline) and may
	// only run ~`streaming_buffer_size` (1 MB by default) ahead of the consumer — so
	// scanning 316 M rows through `SendQuery` ran the joins essentially
	// single-threaded: 28 s wall for 2.7 s of work. Materializing lets DuckDB
	// execute the whole pipeline across all cores first (25 batches: 28.1 s → 6.0 s;
	// 121 batches: 15.2 s → 4.7 s). The per-batch provider always did this — via
	// `conn.Query` in QueryDistanceBlock — which is why that path reached 8.4x
	// effective cores while the wave path sat at 2.1x.
	//
	// BUFFER_MANAGED, not the `conn.Query` default: the result is the wave's largest
	// allocation (~20 bytes × anchors × the wave's samples, ≈ 500 MB at 1000 anchors
	// and a 25-batch wave), and the default IN_MEMORY collector allocates it from
	// Allocator::DefaultAllocator — untracked heap that `memory_limit` cannot see.
	// That is exactly the category that made pcoa() SIGKILL. Buffer-managed, it is
	// counted and spillable instead: the same run under `memory_limit='2GB'`
	// completes at 1.97 GB peak rather than dying.
	static unique_ptr<QueryResult> RunWaveQuery(Connection &conn, const std::string &sql, const char *what) {
		PendingQueryParameters params;
		params.query_parameters.output_type = QueryResultOutputType::FORCE_MATERIALIZED;
		params.query_parameters.memory_type = QueryResultMemoryType::BUFFER_MANAGED;
		auto pending = conn.PendingQuery(sql, params);
		if (pending->HasError()) {
			throw InvalidInputException("progressive_pcoa_from_distances: %s failed: %s", what, pending->GetError());
		}
		auto res = pending->Execute();
		if (res->HasError()) {
			throw InvalidInputException("progressive_pcoa_from_distances: %s failed: %s", what, res->GetError());
		}
		return res;
	}

	// Route one query's rows into the builders. Columns are always
	// (block, position_or_ordinal_a, position_or_ordinal_b, distance); an ordinal is
	// turned into a cell position by adding that block's non-anchor count, since a
	// request is laid out as [batch ids..., anchors...].
	void RunRoutedQuery(Connection &conn, const std::string &sql, bool a_is_anchor_ord, bool b_is_anchor_ord) {
		auto res = RunWaveQuery(conn, sql, "wave scan");
		try {
			while (auto chunk = res->Fetch()) {
				const idx_t rn = chunk->size();
				if (rn == 0) {
					break;
				}
				UnifiedVectorFormat blk_u, pa_u, pb_u, d_u;
				chunk->data[0].ToUnifiedFormat(rn, blk_u);
				chunk->data[1].ToUnifiedFormat(rn, pa_u);
				chunk->data[2].ToUnifiedFormat(rn, pb_u);
				chunk->data[3].ToUnifiedFormat(rn, d_u);
				auto blk_data = UnifiedVectorFormat::GetData<int32_t>(blk_u);
				auto pa_data = UnifiedVectorFormat::GetData<int32_t>(pa_u);
				auto pb_data = UnifiedVectorFormat::GetData<int32_t>(pb_u);
				auto d_data = UnifiedVectorFormat::GetData<double>(d_u);
				for (idx_t i = 0; i < rn; ++i) {
					const auto k = static_cast<size_t>(blk_data[blk_u.sel->get_index(i)]);
					const uint32_t off = batch_len_[k];
					const auto pa = static_cast<uint32_t>(pa_data[pa_u.sel->get_index(i)]);
					const auto pb = static_cast<uint32_t>(pb_data[pb_u.sel->get_index(i)]);
					builders_[k]->Add(a_is_anchor_ord ? off + pa : pa, b_is_anchor_ord ? off + pb : pb,
					                  d_data[d_u.sel->get_index(i)]);
				}
			}
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("progressive_pcoa_from_distances: %s", e.what());
		}
		if (res->HasError()) {
			throw InvalidInputException("progressive_pcoa_from_distances: wave scan failed: %s", res->GetError());
		}
	}

	// Read the anchor×anchor distances once per run. Every block of every wave holds
	// this same corner, so re-reading it per block (as the per-batch provider did)
	// re-scans a×N rows for every batch — the term that made the old path's cost grow
	// with the batch count.
	void LoadAnchorCornerOnce(Connection &conn) {
		if (anchor_corner_loaded_) {
			return;
		}
		auto res =
		    RunWaveQuery(conn,
		                 "SELECT a1.ord, a2.ord, d.distance::DOUBLE FROM " + qname_ + " d JOIN " + wave_anchor_quoted_ +
		                     " a1 ON d.sample_a::VARCHAR = a1.id JOIN " + wave_anchor_quoted_ +
		                     " a2 ON d.sample_b::VARCHAR = a2.id"
		                     " WHERE d.distance IS NOT NULL AND NOT isnan(d.distance::DOUBLE)",
		                 "anchor block scan");
		while (auto chunk = res->Fetch()) {
			const idx_t rn = chunk->size();
			if (rn == 0) {
				break;
			}
			UnifiedVectorFormat oa_u, ob_u, d_u;
			chunk->data[0].ToUnifiedFormat(rn, oa_u);
			chunk->data[1].ToUnifiedFormat(rn, ob_u);
			chunk->data[2].ToUnifiedFormat(rn, d_u);
			auto oa_data = UnifiedVectorFormat::GetData<int32_t>(oa_u);
			auto ob_data = UnifiedVectorFormat::GetData<int32_t>(ob_u);
			auto d_data = UnifiedVectorFormat::GetData<double>(d_u);
			for (idx_t i = 0; i < rn; ++i) {
				anchor_corner_.push_back({static_cast<uint32_t>(oa_data[oa_u.sel->get_index(i)]),
				                          static_cast<uint32_t>(ob_data[ob_u.sel->get_index(i)]),
				                          d_data[d_u.sel->get_index(i)]});
			}
		}
		if (res->HasError()) {
			throw InvalidInputException("progressive_pcoa_from_distances: anchor block scan failed: %s",
			                            res->GetError());
		}
		anchor_corner_loaded_ = true;
	}

	ClientContext &context_;
	std::string qname_;
	std::vector<std::string> anchors_;
	// Per-provider unique names for the two wave staging tables, plus their quoted
	// forms for embedding in SQL. See the constructor for why they are not fixed.
	std::string wave_batch_table_;
	std::string wave_anchor_table_;
	std::string wave_batch_quoted_;
	std::string wave_anchor_quoted_;
	std::vector<miint::progressive::DistanceBlock> blocks_;
	std::vector<std::unique_ptr<miint::unifrac::DenseDistanceMatrixBuilder>> builders_;
	std::vector<uint32_t> batch_len_;                    // per block: non-anchor sample count
	std::unordered_map<std::string, size_t> by_request_; // request.front() → block index
	std::vector<AnchorPair> anchor_corner_;              // anchor×anchor, read once per run
	bool anchor_corner_loaded_ = false;
};

// Build the dense distance block over exactly `requested` (row/col order =
// requested order) by scanning the relation for pairs with both endpoints in the
// set, then delegating the fill + validation to the same BuildDenseDistanceMatrix
// that ReadDistanceTable/pcoa use — so the block enforces the identical contract
// (negative / non-finite / nonzero-self / conflicting-duplicate / incomplete all
// throw), rather than a weaker hand-rolled reimplementation. A batch needs its
// full (anchors + batch)² block present in the relation; a gap fails loud.
miint::progressive::DistanceBlock QueryDistanceBlock(ClientContext &context, const std::string &qname,
                                                     const std::vector<std::string> &requested) {
	auto conn = MakeReadOnlyHelperConnection(context);
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
		part = PartitionWithExplicitAnchors(ids.sorted_ids, explicit_anchors, "progressive_pcoa_from_distances",
		                                    "distance-table");
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
	WaveDistanceBlockSource source(context, qname, part.anchors);
	const miint::progressive::BlockProvider provider = [&source](const std::vector<std::string> &requested) {
		return source.Get(requested);
	};
	const miint::progressive::WavePrefetch prefetch = [&source](const std::vector<std::vector<std::string>> &requests) {
		source.Prefetch(requests);
	};

	// How many batches to serve per relation scan (see ChooseWaveWidth for what a
	// wave costs). The budget is a quarter of what memory_limit currently leaves —
	// a quarter because the blocks are extension heap the buffer manager can neither
	// see nor evict, the same reason pcoa()'s dense matrix is guarded rather than
	// tracked. The scan rows ARE buffer-managed (see RunWaveQuery), so overshooting
	// the estimate spills rather than dies; they are still charged, because spilling
	// a wave's scan result costs far more than running a narrower wave.
	const uint32_t wave_batches = [&]() -> uint32_t {
		auto &buffer_manager = BufferManager::GetBufferManager(context);
		const auto max_memory = buffer_manager.GetMaxMemory();
		const auto used_memory = buffer_manager.GetUsedMemory();
		const uint64_t budget = (max_memory > used_memory ? max_memory - used_memory : 0) / 4;
		const size_t n_batches =
		    (part.remaining.size() + static_cast<size_t>(batch_size) - 1) / static_cast<size_t>(batch_size);
		return miint::progressive::ChooseWaveWidth(part.anchors.size(), static_cast<uint32_t>(batch_size),
		                                           static_cast<uint32_t>(n_threads), budget, n_batches);
	}();

	miint::progressive::ProgressivePcoaResult result;
	try {
		// Waves, always. The per-batch provider re-reads every anchor-touching row for
		// every batch, so its cost grows with the batch count; a wave reads them once.
		// It used to win anyway at low batch counts, purely because its blocks came
		// from a materialized query while the wave scans were streamed one row at a
		// time (see RunWaveQuery) — measured on the 25,145-sample EMP matrix (316 M
		// pairs, 1000 anchors, d=10), 14 cores:
		//
		//   batches   waves (streamed)   waves (materialized)   per-batch queries
		//   25        28.1 s / 58 CPU    6.0 s / 40 CPU         18.1 s / 152 CPU
		//   121       15.2 s / 44 CPU    4.7 s / 36 CPU         67.4 s / 666 CPU
		//
		// With that fixed the wave path wins in both regimes, so the old
		// `n_batches >= 64` crossover and its MIINT_PROGRESSIVE_NO_WAVES escape hatch
		// are gone: one path, no regime to regress.
		//
		// A wave's batches run concurrently, n_threads of them at a time, each
		// ordination pinned to one OpenMP thread — so `threads :=` still bounds total
		// fan-out, it just buys concurrent blocks instead of a wider fsvd. Safe here
		// because a wave's blocks are already in the source's cache (Get is read-only
		// during a wave) and skbb's ordination is re-entrant when seeded per call
		// (see ComputeCallScope). It is worth little on this path (~0.2 s of the 6.0 s; the
		// ordination stage is ~4% of a run at 1000 anchors, ~11% at 3000) — the
		// from_unifrac variant, where a block IS a UniFrac compute, is where it would
		// pay, and that stays serial until libssu's per-compute global `report_status`
		// is thread-local.
		result = miint::progressive::RunProgressivePcoa(part.anchors, part.remaining, static_cast<uint32_t>(n_dims),
		                                                static_cast<uint32_t>(batch_size), seed, n_threads, provider,
		                                                prefetch, wave_batches, static_cast<uint32_t>(n_threads));
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
		row.batch = c.batch;
		if (c.batch >= 0) {
			row.batch_anchor_m2 = result.batches[static_cast<size_t>(c.batch)].anchor_m2;
		}
		data->rows.push_back(std::move(row));
	}

	data->with_batch_diagnostics = true;
	DeclarePcoaOutputSchema(data->sample_id_type, return_types, names, /*with_batch_diagnostics=*/true);
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
	auto conn = MakeReadOnlyHelperConnection(context);
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
	// Called concurrently from block workers, one connection per call. Safe: the
	// inherit only copies the caller's temporary_objects shared_ptr (a read of a
	// pointer nobody is writing) into this thread's own fresh context, and the
	// shared catalog is then only read from.
	auto conn = MakeReadOnlyHelperConnection(context);
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
// pins its own OpenMP width and takes no lock, so blocks may run concurrently.
// The anchor samples' feature rows are identical for every batch, but the slice
// query re-fetched them for each one: with 1000 anchors and batch_size 1000, half
// of every batch's query output was a re-read of rows already seen. Cached once per
// run (bounded by the anchors' own row count, ~20 MB at 1000 anchors), so each
// batch queries only its own samples.
//
// Order does not matter: the rows go to UnifracSupportBiomView::FromCoo, which
// builds its own sorted dictionary, and the core maps block rows back by id.
// Thread safety: the anchor rows are loaded EAGERLY here, in the constructor, on
// the binding thread. They used to be filled lazily on the first RowsFor call
// behind a plain `loaded_` bool, which is a data race the moment blocks run
// concurrently — several workers would write anchor_rows_ while others read it,
// and a vector being resized under a reader is a use-after-free, not a stale
// read. Loading before any fan-out makes every later RowsFor a read-only use of
// shared state, which is what lets all workers share one cache with no lock.
class AnchorFeatureRowCache {
public:
	AnchorFeatureRowCache(ClientContext &context, const std::string &qname, const std::vector<std::string> &anchors)
	    : anchor_set_(anchors.begin(), anchors.end()), anchor_rows_(QueryFeatureRows(context, qname, anchors)) {
	}

	// Safe to call concurrently: reads only, and QueryFeatureRows opens its own
	// Connection per call (see docs/internals/reading-tables-views.md).
	std::vector<miint::unifrac::CooRow> RowsFor(ClientContext &context, const std::string &qname,
	                                            const std::vector<std::string> &requested) const {
		std::vector<std::string> non_anchor;
		non_anchor.reserve(requested.size());
		for (const auto &id : requested) {
			if (!anchor_set_.count(id)) {
				non_anchor.push_back(id);
			}
		}
		if (non_anchor.empty()) {
			return anchor_rows_; // the anchors-only reference block
		}
		auto rows = QueryFeatureRows(context, qname, non_anchor);
		rows.insert(rows.end(), anchor_rows_.begin(), anchor_rows_.end());
		return rows;
	}

private:
	std::unordered_set<std::string> anchor_set_;
	std::vector<miint::unifrac::CooRow> anchor_rows_;
};

miint::progressive::DistanceBlock
ComputeUnifracBlock(ClientContext &context, const std::string &qname, const std::vector<std::string> &requested,
                    const miint::NewickTree &table_tree, const std::string &variant_fp32, bool variance_adjust,
                    double alpha, bool bypass_tips, bool normalize_sample_counts, int seed, int n_threads,
                    const AnchorFeatureRowCache &anchor_cache) {
	auto rows = [&]() {
		return anchor_cache.RowsFor(context, qname, requested);
	}();
	miint::unifrac::UnifracSupportBiomView biom_view = [&]() {
		try {
			return miint::unifrac::UnifracSupportBiomView::FromCoo(std::move(rows));
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("progressive_pcoa_from_unifrac: %s", e.what());
		}
	}();

	// Shear the (already table-scoped) tree to just this block's features. A block
	// is a fraction of the table's samples and so touches a fraction of its
	// features — measured on a 700-sample block of a 197,711-feature table, the
	// block needed 26,741 tips and the compute went 1.83 s -> 0.61 s, bit-identical
	// (244,650 pairs, max abs diff 0.0). The shear and the bptree rebuild are
	// charged below so the trade stays visible if a future table inverts it.
	miint::NewickTree block_tree = [&]() {
		const auto &fids = biom_view.feature_ids();
		std::unordered_set<std::string> keep(fids.begin(), fids.end());
		try {
			return table_tree.shear(keep, /*collapse=*/true, /*ignore_missing=*/false);
		} catch (const std::exception &e) {
			throw InvalidInputException("progressive_pcoa_from_unifrac: failed to shear the tree to a batch's "
			                            "features: %s",
			                            e.what());
		}
	}();
	miint::unifrac::UnifracBptreeView block_bptree = [&]() {
		return miint::unifrac::UnifracBptreeView::FromNewickTree(block_tree);
	}();

	miint::unifrac::UnifracDistanceMatrix dist = [&]() {
		try {
			return miint::unifrac::UnifracDistanceMatrix::Compute(
			    biom_view, block_bptree, variant_fp32, variance_adjust, alpha, bypass_tips, normalize_sample_counts,
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
	vector<string> explicit_anchors; // if non-empty: override seeded random anchor selection
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
		} else if (key == "anchors") {
			for (const auto &child : ListValue::GetChildren(kv.second)) {
				if (child.IsNull()) {
					throw BinderException("progressive_pcoa_from_unifrac: anchors list must not contain NULL");
				}
				explicit_anchors.push_back(child.ToString());
			}
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

	auto ids = [&]() {
		return EnumerateFeatureTableIds(context, table_name, "progressive_pcoa_from_unifrac");
	}();
	const auto n_samples = static_cast<uint32_t>(ids.sorted_sample_ids.size());
	if (n_samples < 2) {
		throw InvalidInputException(
		    "progressive_pcoa_from_unifrac: feature-table '%s' has %u sample(s) with nonzero features; at least 2 are "
		    "required",
		    table_name, n_samples);
	}
	if (static_cast<uint32_t>(n_dims) > n_samples - 1) {
		throw BinderException(
		    "progressive_pcoa_from_unifrac: n_dims (%d) must be <= n_samples - 1 (%u). PCoA loses one "
		    "dimension to centering.",
		    n_dims, n_samples - 1);
	}

	// Build the tree once and validate it covers every feature (each batch's
	// features are a subset, so this upfront check makes them all safe).
	auto tree_inputs = ReadTreeTable(context, tree_name);
	auto tree = [&]() {
		return miint::NewickTree::build(tree_inputs);
	}();
	// Not needed once the tree is built, and shear() below allocates an
	// intermediate copy — release it first (same reasoning as shear_tree.cpp).
	tree_inputs = {};
	{
		try {
			miint::unifrac::ValidateTreeCoversFeatures(tree, ids.feature_ids);
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("progressive_pcoa_from_unifrac: %s", e.what());
		}
	}

	// A reference phylogeny is usually a superset of what one feature table uses,
	// and every batch would otherwise pay for the unused tips: measured, 800,000
	// tips absent from the table cost ~0.49 s on EVERY block. Shear once here so
	// each block's own shear (in ComputeUnifracBlock) starts from a tree already
	// scoped to the table. Validation above guarantees every feature is a tip, so
	// nothing can be missing and the LCA re-rooting drops only branches no sample
	// traverses — distances are unchanged.
	{
		std::unordered_set<std::string> table_features(ids.feature_ids.begin(), ids.feature_ids.end());
		try {
			tree = tree.shear(table_features, /*collapse=*/true, /*ignore_missing=*/false);
		} catch (const std::exception &e) {
			throw InvalidInputException(
			    "progressive_pcoa_from_unifrac: failed to shear the tree to the feature-table's features: %s",
			    e.what());
		}
	}
	const std::string variant_fp32 = variant + "_fp32";

	// Anchor selection: an explicit list takes precedence over n_anchors/seed. Same
	// contract as progressive_pcoa_from_distances — it is what makes the reference
	// frame reproducible against an external oracle, and the way to feed in anchors
	// chosen by some rule other than the seeded random draw (see
	// pick_maxvol_anchors), which matters because anchor QUALITY is the accuracy
	// ceiling of this method.
	AnchorPartition part;
	if (!explicit_anchors.empty()) {
		part = PartitionWithExplicitAnchors(ids.sorted_sample_ids, explicit_anchors, "progressive_pcoa_from_unifrac",
		                                    "feature-table");
		if (part.anchors.size() < static_cast<size_t>(n_dims) + 1) {
			throw BinderException("progressive_pcoa_from_unifrac: %zu explicit anchor(s) given but n_dims + 1 (%d) are "
			                      "required; the reference PCoA and the procrustes fit each need at least that many",
			                      part.anchors.size(), n_dims + 1);
		}
	} else {
		if (n_anchors < 1) {
			throw BinderException("progressive_pcoa_from_unifrac: n_anchors must be >= 1 (got %d)", n_anchors);
		}
		if (static_cast<uint32_t>(n_anchors) > n_samples) {
			throw BinderException(
			    "progressive_pcoa_from_unifrac: n_anchors (%d) exceeds the %u sample(s) in the feature-table",
			    n_anchors, n_samples);
		}
		if (n_anchors < n_dims + 1) {
			throw BinderException(
			    "progressive_pcoa_from_unifrac: n_anchors (%d) must be >= n_dims + 1 (%d); the reference "
			    "PCoA and the procrustes fit each need at least that many anchors",
			    n_anchors, n_dims + 1);
		}
		part = PickAnchors(ids.sorted_sample_ids, static_cast<uint32_t>(n_anchors), seed);
	}

	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);
	AnchorFeatureRowCache anchor_cache(context, qname, part.anchors);

	// How wide each block's libssu compute may fan out. The core pins each worker's
	// ORDINATION to one OpenMP thread itself, but it cannot reach into our provider,
	// so the UniFrac width has to be divided HERE — otherwise W concurrent workers
	// would each pin n_threads and oversubscribe the machine W-fold (196 OpenMP
	// threads at n_threads = 14). Divided by the workers that will actually run, not
	// by n_threads, so a table with fewer batches than threads still gets a wide
	// compute per block instead of being needlessly pinned to one core.
	const auto n_batches =
	    (part.remaining.size() + static_cast<size_t>(batch_size) - 1) / static_cast<size_t>(batch_size);
	const auto workers = static_cast<uint32_t>(n_threads);
	const auto active_workers = std::max<size_t>(1, std::min<size_t>(workers, n_batches));
	const int block_threads = std::max<int>(1, n_threads / static_cast<int>(active_workers));

	const miint::progressive::BlockProvider provider = [&](const std::vector<std::string> &requested) {
		return ComputeUnifracBlock(context, qname, requested, tree, variant_fp32, variance_adjust, alpha, bypass_tips,
		                           normalize_sample_counts, seed, block_threads, anchor_cache);
	};

	// Run the blocks CONCURRENTLY. This is the one path where it pays: a block here
	// is a UniFrac compute (seconds) rather than an fsvd (~9 ms at m=2000), and a
	// single block cannot use more than one core no matter what `threads` says —
	// libssu's parallel degree is ceil(n_samples/2048) and a block is only
	// n_anchors + batch_size samples, so at any sane config that is one stripe.
	// Concurrency across blocks is therefore the only way this path uses the machine.
	//
	// Workers are drawn from one wave, so wave_batches must be set too or
	// batch_workers does nothing (wave_workers = min(batch_workers, wave_count)).
	// Sized equal to the worker count: without a prefetch each block is fetched and
	// released inside its own batch, so live memory is `workers` blocks rather than
	// the whole wave, and a bigger wave would only buy fewer barriers at the cost of
	// holding more batch output. Waves are a barrier, so a straggler can idle its
	// wave's other workers; revisit together with the streaming refactor.
	//
	// Bit-identity is preserved, not merely approximate: each worker pins its
	// ordination to ONE OpenMP thread, which is what the serial path uses at
	// n_threads=1, and skbb's centering reduction sums in thread-count-dependent
	// order (see test_ProgressivePcoa's serial-vs-parallel case).
	miint::progressive::ProgressivePcoaResult result;
	try {
		result = miint::progressive::RunProgressivePcoa(part.anchors, part.remaining, static_cast<uint32_t>(n_dims),
		                                                static_cast<uint32_t>(batch_size), seed, n_threads, provider,
		                                                /*prefetch=*/nullptr, /*wave_batches=*/workers,
		                                                /*batch_workers=*/workers);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("progressive_pcoa_from_unifrac: %s", e.what());
	}

	auto data = make_uniq<UnifracPcoaData>();
	data->sample_id_type = ids.sample_id_type;
	{
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
			row.batch = c.batch;
			if (c.batch >= 0) {
				row.batch_anchor_m2 = result.batches[static_cast<size_t>(c.batch)].anchor_m2;
			}
			data->rows.push_back(std::move(row));
		}
	}

	data->with_batch_diagnostics = true;
	DeclarePcoaOutputSchema(data->sample_id_type, return_types, names, /*with_batch_diagnostics=*/true);
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
	fn.named_parameters["anchors"] = LogicalType::LIST(LogicalType::VARCHAR);
	loader.RegisterFunction(fn);
}

} // namespace duckdb
