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
#include "community_distances.hpp"
#include "id_column_utils.hpp"
#include "miint_log.hpp"
#include "progressive_pcoa_core.hpp"
#include "tree_table_reader.hpp"
#include "unifrac_bptree.hpp"
#include "unifrac_dense_distance.hpp"
#include "unifrac_distance.hpp"
#include "unifrac_function_common.hpp"
#include "unifrac_omp_scope.hpp"
#include "unifrac_support_biom.hpp"

#include "duckdb/common/types/uuid.hpp"
#include "duckdb/common/types/column/column_data_collection.hpp"

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
using unifrac_internal::ProbeDistanceTableIdType;
using unifrac_internal::ProbeFeatureTableIdType;
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

// Inputs for a dense (whole-matrix) PCoA, resolved and validated at bind. Shared
// by pcoa and unifrac_pcoa, which differ only in where the matrix comes from —
// `source` selects that, and the UNIFRAC-only fields below are read accordingly.
//
// It holds INPUTS, not results. The ordination runs per execution, in InitGlobal,
// for the same reasons as the progressive functions (see ProgressivePcoaData) plus
// one specific to these: a result computed in Bind was handed to the scan by
// MOVING it out of the bind data, so a prepared statement — which binds once and
// re-executes that plan — returned rows the first time and NOTHING thereafter.
struct UnifracPcoaData : public TableFunctionData {
	enum class Source { DISTANCES, UNIFRAC };

	Source source = Source::DISTANCES;
	// Output type for sample_id — mirrors the input sample_id type (BIGINT/UUID)
	// or VARCHAR otherwise. See ResolveSampleIdOutputType.
	LogicalType sample_id_type = LogicalType::VARCHAR;
	std::string table_name;
	int32_t n_dims = 3;
	int32_t seed = -1;
	int n_threads = 1;

	// ── UNIFRAC only ──
	std::string tree_name;
	std::string variant_fp32;
	bool variance_adjust = false;
	double alpha = 1.0;
	bool bypass_tips = false;
	bool normalize_sample_counts = true;
	int32_t subsample_depth = 0;
	bool subsample_with_replacement = false;
	int32_t n_subsamples = 1;
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

	auto data = make_uniq<UnifracPcoaData>();
	data->source = UnifracPcoaData::Source::UNIFRAC;
	// Resolved without reading the table, so the schema is known at bind while the
	// ordination is not (see ProbeFeatureTableIdType).
	data->sample_id_type = ProbeFeatureTableIdType(context, table_name, "unifrac_pcoa");
	data->table_name = table_name;
	data->tree_name = tree_name;
	data->n_dims = n_dims;
	data->seed = seed;
	data->n_threads = n_threads;
	// libssu accepts both bare and `_fp32`-suffixed variant strings
	// (unifrac-binaries/src/api.cpp:60-89). We append _fp32 explicitly so the
	// caller's choice is pinned to fp32 even if libssu changes which bare
	// names default to fp32 in a future release.
	data->variant_fp32 = variant + "_fp32";
	data->variance_adjust = variance_adjust;
	data->alpha = alpha;
	data->bypass_tips = bypass_tips;
	data->normalize_sample_counts = normalize_sample_counts;
	data->subsample_depth = subsample_depth;
	data->subsample_with_replacement = subsample_with_replacement;
	data->n_subsamples = n_subsamples;

	DeclarePcoaOutputSchema(data->sample_id_type, return_types, names);

	return std::move(data);
}

// The UniFrac ordination itself: read the table, build and check the tree, then run
// one PCoA per subsample iteration. Runs per execution — see UnifracPcoaData.
void RunUnifracPcoa(ClientContext &context, const UnifracPcoaData &data, std::vector<PcoaRow> &out_rows) {
	auto coo_rows = ReadFeatureTable(context, data.table_name, "unifrac_pcoa");
	if (coo_rows.empty()) {
		throw InvalidInputException("unifrac_pcoa: feature-table '%s' is empty after dropping NULL/zero rows",
		                            data.table_name);
	}
	miint::unifrac::UnifracSupportBiomView biom_view = [&]() {
		try {
			return miint::unifrac::UnifracSupportBiomView::FromCoo(std::move(coo_rows));
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("unifrac_pcoa: %s", e.what());
		}
	}();

	const auto *biom_struct = biom_view.support_biom();
	auto feature_ids = CollectIds(biom_struct->obs_ids, biom_struct->n_obs);
	const auto n_samples = static_cast<uint32_t>(biom_struct->n_samples);

	if (n_samples < 2) {
		throw InvalidInputException(
		    "unifrac_pcoa: feature-table '%s' has %u sample(s); at least 2 are required for PCoA", data.table_name,
		    n_samples);
	}
	if (static_cast<uint32_t>(data.n_dims) > n_samples - 1) {
		throw InvalidInputException(
		    "unifrac_pcoa: n_dims (%d) must be <= n_samples - 1 (%u). PCoA loses one dimension to centering.",
		    data.n_dims, n_samples - 1);
	}

	auto tree_inputs = ReadTreeTable(context, data.tree_name);
	auto tree = miint::NewickTree::build(tree_inputs);
	try {
		miint::unifrac::ValidateTreeCoversFeatures(tree, feature_ids);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("unifrac_pcoa: %s", e.what());
	}
	auto bptree_view = miint::unifrac::UnifracBptreeView::FromNewickTree(tree);

	const auto rows_per_iter = static_cast<size_t>(n_samples) * static_cast<size_t>(data.n_dims);
	out_rows.reserve(static_cast<size_t>(data.n_subsamples) * rows_per_iter);
	for (int32_t i = 0; i < data.n_subsamples; ++i) {
		// seed + i overflow is prevented by the bind-time check.
		const int seed_iter = (data.seed >= 0) ? (data.seed + i) : -1;
		ComputeOneIteration(biom_view, bptree_view, data.variant_fp32, data.variance_adjust, data.alpha,
		                    data.bypass_tips, data.normalize_sample_counts, static_cast<uint32_t>(data.subsample_depth),
		                    data.subsample_with_replacement, seed_iter, static_cast<uint32_t>(data.n_dims),
		                    data.n_threads, i, out_rows);
	}
}

// The metric-agnostic ordination: read the condensed relation into a dense matrix
// and ordinate it. Runs per execution — see UnifracPcoaData.
void RunPcoaFromDistances(ClientContext &context, const UnifracPcoaData &data, std::vector<PcoaRow> &out_rows) {
	auto dist = ReadDistanceTable(context, data.table_name, "pcoa");
	if (static_cast<uint32_t>(data.n_dims) > dist.n_samples - 1) {
		throw InvalidInputException(
		    "pcoa: n_dims (%d) must be <= n_samples - 1 (%u). PCoA loses one dimension to centering.", data.n_dims,
		    dist.n_samples - 1);
	}
	out_rows.reserve(static_cast<size_t>(dist.n_samples) * static_cast<size_t>(data.n_dims));
	RunPcoaOnMatrix(dist.matrix.data(), dist.n_samples, dist.sample_ids, static_cast<uint32_t>(data.n_dims), data.seed,
	                data.n_threads, /*iteration_index*/ 0, "pcoa", out_rows);
}

// The ordination happens HERE, once per execution, rather than in Bind.
//
// Three things follow, all of which were wrong before. A prepared statement binds
// once and re-executes the same plan, so a Bind-computed result had to be either
// moved (second EXECUTE returned nothing — the bug this fixes) or copied (paying
// for the result twice). `EXPLAIN` binds without executing, so it used to run the
// whole ordination: measured in-process on a 4,000-sample distance relation,
// EXPLAIN took 0.644 s against 0.638 s to actually run the query, and now takes
// 0.000 s while the query still takes 0.648 s. And an unseeded run, re-executed,
// should redraw — which it now does, rather than replaying the first execution's
// subsample.
unique_ptr<GlobalTableFunctionState> UnifracPcoaInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<UnifracPcoaData>();
	auto gstate = make_uniq<UnifracPcoaGlobalState>();
	gstate->sample_id_type = data.sample_id_type;
	if (data.source == UnifracPcoaData::Source::UNIFRAC) {
		RunUnifracPcoa(context, data, gstate->rows);
	} else {
		RunPcoaFromDistances(context, data, gstate->rows);
	}
	return std::move(gstate);
}

// Page up to one vector of rows out of `rows`, starting at `cursor` and advancing
// it. Shared by the dense functions (pcoa / unifrac_pcoa, whose `rows` is the whole
// result) and the progressive ones (whose `rows` is refilled one wave at a time),
// so the two can never drift apart in how a row becomes output.
void EmitPcoaChunk(const std::vector<PcoaRow> &rows, size_t &cursor, const LogicalType &sample_id_type,
                   bool with_batch_diagnostics, DataChunk &output) {
	const idx_t remaining = rows.size() - cursor;
	const idx_t n = std::min<idx_t>(STANDARD_VECTOR_SIZE, remaining);

	auto iter_data = FlatVector::GetData<int32_t>(output.data[0]);
	auto &sample_id_vec = output.data[1];
	auto axis_data = FlatVector::GetData<int32_t>(output.data[2]);
	auto coord_data = FlatVector::GetData<double>(output.data[3]);
	auto eig_data = FlatVector::GetData<double>(output.data[4]);
	auto pe_data = FlatVector::GetData<double>(output.data[5]);
	int32_t *batch_data = nullptr;
	double *m2_data = nullptr;
	if (with_batch_diagnostics) {
		batch_data = FlatVector::GetData<int32_t>(output.data[6]);
		m2_data = FlatVector::GetData<double>(output.data[7]);
	}

	for (idx_t i = 0; i < n; ++i) {
		const auto &r = rows[cursor + i];
		iter_data[i] = r.iteration;
		// EmitIdCell mirrors the id type; its ""/"*"→NULL sentinel branch is
		// unreachable here (ReadFeatureTable drops NULL sample_ids).
		EmitIdCell(sample_id_vec, i, r.sample_id, sample_id_type);
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

	cursor += n;
	output.SetCardinality(n);
}

void UnifracPcoaExecute(ClientContext &, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<UnifracPcoaGlobalState>();
	if (gstate.cursor >= gstate.rows.size()) {
		output.SetCardinality(0);
		return;
	}
	EmitPcoaChunk(gstate.rows, gstate.cursor, gstate.sample_id_type, /*with_batch_diagnostics=*/false, output);
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

	auto data = make_uniq<UnifracPcoaData>();
	data->source = UnifracPcoaData::Source::DISTANCES;
	// Resolved without reading the relation, so a mis-shaped or mixed-id-type
	// relation is still rejected at bind while the matrix is not built until the
	// query runs (see ProbeDistanceTableIdType).
	data->sample_id_type = ProbeDistanceTableIdType(context, table_name, "pcoa");
	data->table_name = table_name;
	data->n_dims = n_dims;
	data->seed = seed;
	data->n_threads = n_threads;

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
// documented caveat), plus the two batch-diagnostic columns.
//
// Bind validates and picks the anchors; the run itself is driven from the scan, one
// wave at a time (see ProgressivePcoaGlobalState). The distance blocks are sourced
// by querying the relation a wave at a time — bounded memory, one wave of blocks
// plus the id set, and the dense matrix is never materialized.

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
// Shared by every progressive_pcoa_from_* bind, so `caller_name` and `source_desc` carry
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

// Everything a progressive run needs, resolved and validated at bind time. Both
// progressive functions share it: the parameters are the same apart from how a
// block is obtained, which `source` selects.
//
// It holds INPUTS, not results. The run itself happens at execution time
// (ProgressivePcoaExecute) rather than in Bind, because a Bind that computes the
// whole ordination has to hand it over as one materialized vector — which the
// executor then copies again — and offers nowhere to poll for cancellation. See
// ProgressivePcoaGlobalState.
struct ProgressivePcoaData : public TableFunctionData {
	// Rotate the finished configuration onto its own principal axes before emitting.
	// ON by default, because without it the emitted axes are the ANCHOR block's
	// principal axes and "PC1" is not the leading axis of the output. It costs the
	// ability to stream — see StageRotatedRun — so it is a named parameter rather
	// than unconditional.
	bool global_rotation = true;
	// Which block source serves this run. The UNIFRAC-only and FEATURES-only
	// fields below are read only when the source is that one.
	enum class Source { DISTANCES, UNIFRAC, FEATURES };

	const char *CallerName() const {
		switch (source) {
		case Source::UNIFRAC:
			return "progressive_pcoa_from_unifrac";
		case Source::FEATURES:
			return "progressive_pcoa_from_features";
		default:
			return "progressive_pcoa_from_distances";
		}
	}

	Source source = Source::DISTANCES;
	// Output type for sample_id — mirrors the input id type. See UnifracPcoaData.
	LogicalType sample_id_type = LogicalType::VARCHAR;
	std::string qname; // quoted source relation (distance table or feature table)
	AnchorPartition part;
	uint32_t n_dims = 3;
	uint32_t batch_size = 1000;
	int seed = -1;
	int n_threads = 1;

	// ── UNIFRAC + FEATURES (any source that COMPUTES its blocks) ──
	// How wide one block's compute may fan out, and how many blocks run at once.
	// Resolved together at bind (ResolveBlockConcurrency) so their product cannot
	// oversubscribe, and read by MakeComputedBlockRun for either source.
	int block_threads = 1;
	// Per-block dense-operand allowance, already divided by the concurrent workers —
	// see ResolveBlockConcurrency. 0 would mean "library default", which is per call and
	// would therefore be W times too generous.
	size_t gram_operand_bytes = 0;
	uint32_t workers = 1;

	// ── FEATURES only ──
	// Which community distance a block is made of. Validated at bind as both a
	// known metric AND a pairwise-local one (see IsPairwiseLocalCommunityMetric).
	std::string metric;

	// ── UNIFRAC only ──
	// Built, validated against the feature table, and sheared to its tips once at
	// bind time; every block shears it again to its own features. Held by pointer
	// only because NewickTree is not default-constructible in a way this struct
	// could use.
	unique_ptr<miint::NewickTree> tree;
	std::string variant_fp32;
	bool variance_adjust = false;
	double alpha = 1.0;
	bool bypass_tips = false;
	bool normalize_sample_counts = true;
};

// How a computed-block run splits `threads` between concurrent blocks and the width
// of one block's own compute.
//
// The progressive core pins each worker's ORDINATION to one OpenMP thread itself, but
// it cannot reach into a provider, so the per-block width has to be divided HERE or W
// concurrent workers would each claim n_threads and oversubscribe the machine W-fold.
// Divided by the workers that will actually run rather than by n_threads, so a table
// with fewer batches than threads still gets a wide compute per block instead of
// being needlessly pinned to one core.
//
// Shared by the UniFrac and community-distance binds for the same reason
// MakeComputedBlockRun is shared: the two halves of the oversubscription contract
// must not drift apart.
// The same division applies to MEMORY, not only to threads. community_distances caps
// the dense operand PER CALL and cannot see how many calls are in flight, so a run
// with W concurrent blocks would be entitled to W times that cap. This is the run's
// total allowance for dense operands, split the same way `block_threads` is; a block
// that does not fit its share stays on the sparse merge, which needs only CSR.
constexpr size_t kRunGramOperandBudgetBytes = 512ull << 20;

struct BlockConcurrency {
	uint32_t workers = 1;
	int block_threads = 1;
	size_t gram_operand_bytes = kRunGramOperandBudgetBytes;
};

BlockConcurrency ResolveBlockConcurrency(int n_threads, size_t n_batches) {
	BlockConcurrency c;
	c.workers = static_cast<uint32_t>(n_threads);
	const auto active = std::max<size_t>(1, std::min<size_t>(c.workers, n_batches));
	c.block_threads = std::max<int>(1, n_threads / static_cast<int>(active));
	c.gram_operand_bytes = kRunGramOperandBudgetBytes / active;
	return c;
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
	bool global_rotation = true;     // see ProgressivePcoaData::global_rotation
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
		} else if (key == "global_rotation") {
			global_rotation = kv.second.GetValue<bool>();
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

	auto data = make_uniq<ProgressivePcoaData>();
	data->source = ProgressivePcoaData::Source::DISTANCES;
	data->sample_id_type = ids.sample_id_type;
	data->qname = KeywordHelper::WriteOptionallyQuoted(table_name);
	data->part = std::move(part);
	data->n_dims = static_cast<uint32_t>(n_dims);
	data->batch_size = static_cast<uint32_t>(batch_size);
	data->seed = seed;
	data->n_threads = n_threads;
	data->global_rotation = global_rotation;

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
// slice of the full matrix. Reuses the B1 core and the whole streaming scan; only
// the BlockProvider differs from the _distances variant. subsample_depth is fixed at 0
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

// Warn when the feature table is not stored in sample_id order.
//
// WHY: batches are contiguous ranges of sorted sample_id, so a batch's slice query
// prunes to its own rows when the table's physical order agrees — and reads the
// whole table when it does not, turning the run's slicing cost from one pass into
// n_batches passes. Measured on a 13.1 M-row table with ASV-sequence ids, one
// batch's slice took 8 ms sorted against 140 ms unsorted (see the sort-order note
// in docs/diversity.md); at scale that is the difference between minutes and hours,
// and nothing in the output reveals it. Sort order never changes the coordinates,
// so this is a warning, not an error.
//
// The probe counts DESCENTS — adjacent rows, in physical order, whose sample_id
// goes backwards. A sequence is sorted if and only if it never descends, so
// `descents > 0` is an exact test rather than a heuristic, and it needs one pass
// over the id column alone: measured 0.09 s over 4.7 M rows and 19.6 s over 1.85 B,
// against runs that take minutes to hours on tables of those sizes.
//
// Note what a cheaper threshold would miss. Scaling the bar to the row-group count
// (pruning's granularity) looks reasonable and is wrong: a feature-major table with
// few features descends only once per feature block — a handful of times — while
// every row group still spans the whole id range, so pruning is dead and the run
// pays in full. Exactness costs nothing here, so take it.
void WarnIfFeatureTableUnsorted(ClientContext &context, const char *caller, const std::string &qname,
                                const std::string &table_name, size_t n_batches) {
	auto conn = MakeReadOnlyHelperConnection(context);
	// lag() over the whole relation with no ORDER BY reads it in physical order,
	// which is the order the slice queries will have to prune against. That relies on
	// the scan reaching the window sink in storage order; checked to still hold under
	// `preserve_insertion_order=false` (a sorted 4.7 M-row table reports 0 descents
	// either way), since a single scan feeding one unpartitioned window has nothing to
	// reorder. A future plan that broke that would over-report, i.e. warn about a
	// sorted table — noise, never a wrong result.
	auto res = conn.Query("SELECT count(*), count(*) FILTER (WHERE t.prev IS NOT NULL AND t.sid < t.prev) FROM "
	                      "(SELECT sample_id::VARCHAR AS sid, lag(sample_id::VARCHAR) OVER () AS prev FROM " +
	                      qname + ") t");
	if (res->HasError()) {
		// A diagnostic must not fail a valid run. Say so rather than swallowing it —
		// silence would be indistinguishable from "your table is fine".
		miint::EmitWarning(context,
		                   "%s: could not check whether feature-table '%s' is "
		                   "stored in sample_id order: %s",
		                   caller, table_name, res->GetError());
		return;
	}
	auto &mat = res->Cast<MaterializedQueryResult>();
	auto chunk = mat.Fetch();
	if (!chunk || chunk->size() == 0) {
		return;
	}
	const int64_t rows = chunk->GetValue(0, 0).GetValue<int64_t>();
	const int64_t descents = chunk->GetValue(1, 0).GetValue<int64_t>();
	if (descents <= 0) {
		return;
	}
	miint::EmitWarning(context,
	                   "%s: feature-table '%s' is not stored in sample_id order (%lld of "
	                   "%lld rows step backwards). Batches are contiguous sample_id ranges, so each of the %llu "
	                   "batches reads the whole table instead of only its own rows. Store it sorted — e.g. CREATE "
	                   "TABLE t AS SELECT * FROM ... ORDER BY sample_id::VARCHAR — which changes only what the run "
	                   "reads, never the coordinates it produces.",
	                   caller, table_name, static_cast<long long>(descents), static_cast<long long>(rows),
	                   static_cast<unsigned long long>(n_batches));
}

// Read the feature rows of exactly the requested samples (all their features) into
// COO form, applying ReadFeatureTable's NULL/zero/NaN drops.
std::vector<miint::unifrac::CooRow> QueryFeatureRows(ClientContext &context, const char *caller,
                                                     const std::string &qname,
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
		throw InvalidInputException("%s: feature slice query failed: %s", caller, res->GetError());
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
	AnchorFeatureRowCache(ClientContext &context, const char *caller, const std::string &qname,
	                      const std::vector<std::string> &anchors)
	    : caller_(caller), anchor_set_(anchors.begin(), anchors.end()),
	      anchor_rows_(QueryFeatureRows(context, caller, qname, anchors)) {
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
		auto rows = QueryFeatureRows(context, caller_, qname, non_anchor);
		rows.insert(rows.end(), anchor_rows_.begin(), anchor_rows_.end());
		return rows;
	}

	// The anchors' own cells, for a caller that needs to build something over just
	// them (the anchor-corner cache). By reference: on a dense table this is
	// anchors x features rows and copying it would undo the point of caching.
	const std::vector<miint::unifrac::CooRow> &anchor_rows() const {
		return anchor_rows_;
	}

private:
	const char *caller_;
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

// Compute the community-distance block over exactly the requested samples: slice
// the feature table -> CSR over the block's OWN features ->
// CommunityDistancesCondensedSparse -> the row-major fp32 matrix the core consumes.
//
// CSR rather than a dense sample x feature matrix because the pair loop is the whole
// cost here and a dense one pays for the entire feature space on every pair. Measured
// on the 1.2M-sample reference table, one default-sized block spans 11,018 features
// but averages 89 nonzeros per sample, so the merge does ~62x less arithmetic for a
// bit-identical answer (test_CommunityDistances pins that equality).
//
// Correct because bind admits only pairwise-local metrics. A block carries just
// the features its own samples use, and for such a metric a feature that is zero
// in both samples of a pair contributes nothing to that pair -- so the block's
// distances are exactly the corresponding slice of the full matrix, and the
// missing columns cost nothing. (test_CommunityDistances pins this bitwise, and
// the SQL parity oracle pins it end to end, per metric.)
//
// Block rows are laid out in the caller's `requested` order, NOT sorted. Any order
// is correct — the core maps rows back by id, and a pairwise-local metric's value
// does not depend on where its two rows sit — but matching `requested` makes this
// block bit-identical to the one progressive_pcoa_from_distances builds for the
// same samples. That matters because the ordination is a RANDOMIZED SVD: permuting
// the rows perturbs its result in the last bits, which showed up as a procrustes
// disparity of ~1e-9 between the two paths. Same order, no perturbation.
//
// The FEATURE dictionary is first-seen order over this block's rows, so its column
// order differs from block to block and from community_distances' whole-table
// dictionary. Which features are present is identical either way — that is what
// pairwise locality buys — but the order in which a metric's inner loop accumulates
// them is not, and floating-point addition is not associative. Values that are
// exactly representable (integer counts) are unaffected; others could differ in the
// last ulp of the double, which the fp32 narrowing then absorbs. Measured, not
// assumed — see the dictionary comment below.
//
// A requested sample with no rows would leave an all-zero row rather than vanish
// from the block, which is what the metrics' empty-community conventions expect
// and what keeps the core's block-coverage check from firing spuriously.
//
// Thread safety: reads only, and RowsFor opens its own Connection per call, so
// blocks may run concurrently.
// One block's cells as CSR, in `requested` order. Split out of
// ComputeCommunityBlock so the anchor-corner cache below can build the anchors'
// own CSR with exactly the same code — a second implementation of the dictionary
// and duplicate-coalescing rules would be a second place for them to drift.
struct BlockCsr {
	std::vector<uint32_t> indptr;
	std::vector<uint32_t> indices;
	std::vector<double> values;
	uint32_t n_features = 0;
};

BlockCsr BuildBlockCsr(const std::vector<miint::unifrac::CooRow> &rows, const std::vector<std::string> &requested) {
	const auto n = static_cast<uint32_t>(requested.size());
	std::unordered_map<std::string, uint32_t> s_index;
	s_index.reserve(n);
	for (uint32_t i = 0; i < n; ++i) {
		s_index.emplace(requested[i], i);
	}

	// Feature dictionary in first-seen order, the same convention
	// community_distances_function.cpp uses. Column order is a per-block property, so
	// in principle it changes the ORDER the metric's inner loop accumulates in, and
	// floating-point addition is not associative — but it cannot change what this
	// function returns, because the distances are narrowed to fp32 before the core
	// sees them and the effect is at the last ulp of a double. Checked rather than
	// assumed: a fixture of non-representable abundances (1/k + sqrt(k)/7) placed
	// every sample at max abs difference 0.0 from progressive_pcoa_from_distances
	// under both this order and a sorted one, at two different batch sizes. Sorting
	// the dictionary would therefore buy nothing observable, so it is not done.
	// (Within a ROW the indices still ascend — the CSR contract requires it, and the
	// merge below depends on it.)
	std::unordered_map<std::string, uint32_t> f_index;
	uint32_t f = 0;
	for (const auto &r : rows) {
		if (f_index.emplace(r.feature_id, f).second) {
			++f;
		}
	}

	// Bucket the block's cells by sample, then sort each row by feature index and
	// coalesce duplicate cells (summed, as the dense path summed them) so the CSR
	// rows are strictly ascending.
	std::vector<std::vector<std::pair<uint32_t, double>>> per_sample(n);
	for (const auto &r : rows) {
		const auto it = s_index.find(r.sample_id);
		if (it == s_index.end()) {
			continue; // not part of this block (the anchor cache is shared across blocks)
		}
		per_sample[it->second].emplace_back(f_index.at(r.feature_id), r.count);
	}
	BlockCsr csr;
	csr.n_features = f;
	auto &indptr = csr.indptr;
	auto &indices = csr.indices;
	auto &values = csr.values;
	indptr.reserve(static_cast<size_t>(n) + 1);
	indices.reserve(rows.size());
	values.reserve(rows.size());
	indptr.push_back(0);
	for (uint32_t i = 0; i < n; ++i) {
		auto &row = per_sample[i];
		// STABLE: duplicate (sample, feature) cells are summed just below, and
		// community_distances sums them in the order its scan produced. A non-stable
		// sort could order equal-index cells differently and land on different last
		// bits, which the SQL parity test asserts are equal.
		std::stable_sort(row.begin(), row.end(),
		                 [](const std::pair<uint32_t, double> &a, const std::pair<uint32_t, double> &b) {
			                 return a.first < b.first;
		                 });
		for (size_t c = 0; c < row.size(); ++c) {
			if (c > 0 && row[c].first == row[c - 1].first) {
				values.back() += row[c].second;
			} else {
				indices.push_back(row[c].first);
				values.push_back(row[c].second);
			}
		}
		indptr.push_back(static_cast<uint32_t>(indices.size()));
	}
	return csr;
}

// The anchor x anchor quadrant of every block, computed once.
//
// WHY: the core requests each batch as (batch ids ++ ALL anchors), so that quadrant
// is re-derived in every block of the run, and re-derived to the SAME VALUE — that is
// exactly what IsPairwiseLocalCommunityMetric guarantees, and what
// test_CommunityDistances' "unchanged by the block a pair lands in" case pins — so
// the recomputation buys nothing.
//
// "Same value" here means mathematically, not to the last bit. This block's feature
// dictionary is built over the anchor rows alone, so it is a different permutation of
// the same features than a batch block's (which sees the batch's features first), and
// every metric's inner sum is accumulated in ascending dictionary order. See
// BuildBlockCsr's note on exactly this: the difference is at the last ulp of a double,
// below the fp32 narrowing the core applies before it sees any of it, and it was
// checked rather than assumed. It is a(a-1)/2 of a block's (a+k)(a+k-1)/2
// pairs: 25% of the work at the documented defaults (a = k = 1000), 44% at
// k = 500, and 64% at a = 200, k = 50, which is the small-batch regime the
// accuracy guidance recommends when only the leading axes will be interpreted.
//
// Held as the condensed upper triangle in ANCHOR order and narrowed to fp32 at
// the splice, the same single narrowing the uncached path applies — so a cached
// block and a recomputed one hand the core identical fp32 matrices. 4 MB at
// 1000 anchors.
//
// Read-only after construction, so every concurrent block shares one.
class AnchorCommunityCorner {
public:
	AnchorCommunityCorner(const std::vector<miint::unifrac::CooRow> &anchor_rows,
	                      const std::vector<std::string> &anchors, const std::string &metric, int n_threads)
	    : anchors_(anchors) { // copied, not aliased: a few tens of KB buys one less lifetime rule
		if (anchors_.size() < 2) {
			return; // no pairs to cache
		}
		const auto csr = BuildBlockCsr(anchor_rows, anchors_);
		const auto a = static_cast<uint32_t>(anchors_.size());
		try {
			condensed_ = miint::CommunityDistancesCondensedSparse(
			    csr.indptr, csr.indices, csr.values, a, csr.n_features, metric, static_cast<unsigned>(n_threads));
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("progressive_pcoa_from_features: %s", e.what());
		}
	}

	uint32_t size() const {
		return static_cast<uint32_t>(anchors_.size());
	}

	// True iff `requested` really ends with the anchor list, in anchor order — the
	// only arrangement under which the cached quadrant belongs where it would be
	// spliced. Checked per block rather than assumed: getting it wrong would not
	// raise, it would silently place every sample against the wrong anchors. O(a)
	// string compares against ~a*k distance computations is free.
	bool MatchesTailOf(const std::vector<std::string> &requested) const {
		if (condensed_.empty() || requested.size() < anchors_.size()) {
			return false;
		}
		const size_t head = requested.size() - anchors_.size();
		for (size_t p = 0; p < anchors_.size(); ++p) {
			if (requested[head + p] != anchors_[p]) {
				return false;
			}
		}
		return true;
	}

	// Anchor pair (p, q), p < q, in the condensed upper-triangle layout the pure
	// core uses.
	double At(uint32_t p, uint32_t q) const {
		return condensed_[miint::CondensedRowBase(static_cast<uint32_t>(anchors_.size()), p) + (q - p - 1)];
	}

private:
	const std::vector<std::string> anchors_;
	std::vector<double> condensed_;
};

miint::progressive::DistanceBlock
ComputeCommunityBlock(ClientContext &context, const std::string &qname, const std::vector<std::string> &requested,
                      const std::string &metric, int n_threads, const AnchorFeatureRowCache &anchor_cache,
                      const AnchorCommunityCorner &corner, size_t gram_operand_bytes) {
	const auto rows = anchor_cache.RowsFor(context, qname, requested);
	const auto n = static_cast<uint32_t>(requested.size());
	const auto csr = BuildBlockCsr(rows, requested);
	const uint32_t f = csr.n_features;

	// Skip the anchor quadrant when this request really carries it as its tail;
	// it is spliced back in from the cache below.
	const bool use_corner = corner.MatchesTailOf(requested);
	const uint32_t cached_tail = use_corner ? corner.size() : 0;

	std::vector<double> condensed;
	try {
		condensed =
		    miint::CommunityDistancesCondensedSparse(csr.indptr, csr.indices, csr.values, n, f, metric,
		                                             static_cast<unsigned>(n_threads), cached_tail, gram_operand_bytes);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("progressive_pcoa_from_features: %s", e.what());
	}
	if (use_corner) {
		// Anchor pair (p, q) sits at block indices (head + p, head + q); the condensed
		// slot for (i, j) is RowBase(n, i) + (j - i - 1).
		const uint32_t head = n - cached_tail;
		for (uint32_t p = 0; p + 1 < cached_tail; ++p) {
			const size_t row_base = miint::CondensedRowBase(n, head + p);
			for (uint32_t q = p + 1; q < cached_tail; ++q) {
				condensed[row_base + (q - p - 1)] = corner.At(p, q);
			}
		}
	}

	// Condensed upper triangle -> dense symmetric fp32, zero diagonal. The fp32
	// narrowing is the same one progressive_pcoa_from_distances applies to the
	// DOUBLE it reads back out of a community_distances table, so for exactly
	// representable abundances the two paths hand the core identical blocks (the
	// SQL parity test asserts the coordinates match to the bit); see the note on
	// feature-dictionary order above for the one case that is only identical after
	// the narrowing.
	miint::progressive::DistanceBlock block;
	block.ids = requested;
	block.matrix.assign(static_cast<size_t>(n) * n, 0.0f);
	size_t k = 0;
	for (uint32_t i = 0; i + 1 < n; ++i) {
		for (uint32_t j = i + 1; j < n; ++j, ++k) {
			const auto d = static_cast<float>(condensed[k]);
			block.matrix[static_cast<size_t>(i) * n + j] = d;
			block.matrix[static_cast<size_t>(j) * n + i] = d;
		}
	}
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
	bool global_rotation = true;     // see ProgressivePcoaData::global_rotation
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
		} else if (key == "global_rotation") {
			global_rotation = kv.second.GetValue<bool>();
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

	// See ResolveBlockConcurrency for why the width is divided here. Concretely on
	// this path: without it, W concurrent workers each pin n_threads to libssu — 196
	// OpenMP threads at n_threads = 14.
	const auto n_batches =
	    (part.remaining.size() + static_cast<size_t>(batch_size) - 1) / static_cast<size_t>(batch_size);
	const auto concurrency = ResolveBlockConcurrency(n_threads, n_batches);

	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);
	// Only worth probing — and only worth telling the user about — when there is more
	// than one batch: a single batch reads the table once whatever its order, so
	// sorting it would save nothing and the warning would be noise.
	if (n_batches > 1) {
		WarnIfFeatureTableUnsorted(context, "progressive_pcoa_from_unifrac", qname, table_name, n_batches);
	}

	auto data = make_uniq<ProgressivePcoaData>();
	data->source = ProgressivePcoaData::Source::UNIFRAC;
	data->sample_id_type = ids.sample_id_type;
	data->qname = qname;
	data->part = std::move(part);
	data->n_dims = static_cast<uint32_t>(n_dims);
	data->batch_size = static_cast<uint32_t>(batch_size);
	data->seed = seed;
	data->n_threads = n_threads;
	data->global_rotation = global_rotation;
	data->tree = make_uniq<miint::NewickTree>(std::move(tree));
	data->variant_fp32 = variant_fp32;
	data->variance_adjust = variance_adjust;
	data->alpha = alpha;
	data->bypass_tips = bypass_tips;
	data->normalize_sample_counts = normalize_sample_counts;
	data->block_threads = concurrency.block_threads;
	data->gram_operand_bytes = concurrency.gram_operand_bytes;
	data->workers = concurrency.workers;

	DeclarePcoaOutputSchema(data->sample_id_type, return_types, names, /*with_batch_diagnostics=*/true);
	return std::move(data);
}

// ── progressive_pcoa_from_features(feature_table, metric, ...) — the tree-free path ──
// The same reference-anchored progressive PCoA, with each batch's block computed
// on the fly as a community_distances matrix over (anchors + batch) rather than a
// UniFrac one — so a non-phylogenetic analysis gets the same memory-bounded
// ordination, instead of being capped at whatever N a dense matrix fits in.
//
// Correct for the SAME reason the UniFrac path is: the metric must be pairwise
// local. Here that means two things, and both are load-bearing. d(i,j) must read
// no statistic taken over the other samples, AND it must be unchanged by dropping
// the features that are zero in both i and j — because a block carries only the
// features its own samples use, which differs from block to block. Five of the
// eight community metrics satisfy both; pearson, chisq and gower do not, and are
// refused at bind rather than computed into a silently wrong ordination. The
// classification itself lives with the metric definitions
// (IsPairwiseLocalCommunityMetric), where a future metric's author will meet it.
unique_ptr<FunctionData> ProgressivePcoaFromFeaturesBind(ClientContext &context, TableFunctionBindInput &input,
                                                         vector<LogicalType> &return_types, vector<string> &names) {
	const std::string table_name = input.inputs[0].GetValue<string>();
	const std::string metric = StringUtil::Lower(input.inputs[1].GetValue<string>());
	if (table_name.empty()) {
		throw BinderException("progressive_pcoa_from_features: feature-table name must not be empty");
	}
	if (!miint::IsValidCommunityMetric(metric)) {
		throw BinderException("progressive_pcoa_from_features: unknown metric '%s' (must be one of %s)", metric,
		                      miint::CommunityMetricList());
	}
	if (!miint::IsPairwiseLocalCommunityMetric(metric)) {
		// Refused rather than computed. Each of these reads a statistic over the
		// whole matrix — pearson the feature-space size, chisq the column sums and
		// grand total, gower the per-feature ranges — so a per-block value is a
		// different metric in every block, and nothing in the output would show it.
		throw BinderException(
		    "progressive_pcoa_from_features: metric '%s' cannot be computed progressively. It depends on statistics "
		    "taken over the whole feature table, so each batch would silently measure a different distance and the "
		    "ordination would be wrong with no error raised. Use progressive_pcoa_from_distances over "
		    "community_distances('%s', '%s') instead, which forms the full matrix and is therefore limited to sample "
		    "counts that fit in memory. Progressive metrics: %s",
		    metric, table_name, metric, miint::PairwiseLocalCommunityMetricList());
	}

	int32_t n_dims = 3;
	int32_t n_anchors = 100;
	int32_t batch_size = 1000;
	int32_t seed = -1;
	int32_t threads = 0;             // 0 = follow DuckDB's TaskScheduler::NumberOfThreads()
	bool global_rotation = true;     // see ProgressivePcoaData::global_rotation
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
		} else if (key == "global_rotation") {
			global_rotation = kv.second.GetValue<bool>();
		} else if (key == "anchors") {
			for (const auto &child : ListValue::GetChildren(kv.second)) {
				if (child.IsNull()) {
					throw BinderException("progressive_pcoa_from_features: anchors list must not contain NULL");
				}
				explicit_anchors.push_back(child.ToString());
			}
		}
	}
	const int n_threads = ResolveThreadsParameter(context, threads, "progressive_pcoa_from_features");
	if (n_dims < 1) {
		throw BinderException("progressive_pcoa_from_features: n_dims must be >= 1 (got %d)", n_dims);
	}
	if (batch_size < 1) {
		throw BinderException("progressive_pcoa_from_features: batch_size must be >= 1 (got %d)", batch_size);
	}

	auto ids = EnumerateFeatureTableIds(context, table_name, "progressive_pcoa_from_features");
	const auto n_samples = static_cast<uint32_t>(ids.sorted_sample_ids.size());
	if (n_samples < 2) {
		throw InvalidInputException(
		    "progressive_pcoa_from_features: feature-table '%s' has %u sample(s) with nonzero features; at least 2 are "
		    "required",
		    table_name, n_samples);
	}
	if (static_cast<uint32_t>(n_dims) > n_samples - 1) {
		throw BinderException("progressive_pcoa_from_features: n_dims (%d) must be <= n_samples - 1 (%u). PCoA loses "
		                      "one dimension to centering.",
		                      n_dims, n_samples - 1);
	}

	AnchorPartition part;
	if (!explicit_anchors.empty()) {
		part = PartitionWithExplicitAnchors(ids.sorted_sample_ids, explicit_anchors, "progressive_pcoa_from_features",
		                                    "feature-table");
		if (part.anchors.size() < static_cast<size_t>(n_dims) + 1) {
			throw BinderException(
			    "progressive_pcoa_from_features: %zu explicit anchor(s) given but n_dims + 1 (%d) are "
			    "required; the reference PCoA and the procrustes fit each need at least that many",
			    part.anchors.size(), n_dims + 1);
		}
	} else {
		if (n_anchors < 1) {
			throw BinderException("progressive_pcoa_from_features: n_anchors must be >= 1 (got %d)", n_anchors);
		}
		if (static_cast<uint32_t>(n_anchors) > n_samples) {
			throw BinderException(
			    "progressive_pcoa_from_features: n_anchors (%d) exceeds the %u sample(s) in the feature-table",
			    n_anchors, n_samples);
		}
		if (n_anchors < n_dims + 1) {
			throw BinderException("progressive_pcoa_from_features: n_anchors (%d) must be >= n_dims + 1 (%d); the "
			                      "reference PCoA and the procrustes fit each need at least that many anchors",
			                      n_anchors, n_dims + 1);
		}
		part = PickAnchors(ids.sorted_sample_ids, static_cast<uint32_t>(n_anchors), seed);
	}

	// See ResolveBlockConcurrency.
	const auto n_batches =
	    (part.remaining.size() + static_cast<size_t>(batch_size) - 1) / static_cast<size_t>(batch_size);
	const auto concurrency = ResolveBlockConcurrency(n_threads, n_batches);

	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);
	// Only worth warning about with more than one batch: a single batch reads the
	// table once whatever its order.
	if (n_batches > 1) {
		WarnIfFeatureTableUnsorted(context, "progressive_pcoa_from_features", qname, table_name, n_batches);
	}

	auto data = make_uniq<ProgressivePcoaData>();
	data->source = ProgressivePcoaData::Source::FEATURES;
	data->sample_id_type = ids.sample_id_type;
	data->qname = qname;
	data->part = std::move(part);
	data->n_dims = static_cast<uint32_t>(n_dims);
	data->batch_size = static_cast<uint32_t>(batch_size);
	data->seed = seed;
	data->n_threads = n_threads;
	data->global_rotation = global_rotation;
	data->metric = metric;
	data->block_threads = concurrency.block_threads;
	data->gram_operand_bytes = concurrency.gram_operand_bytes;
	data->workers = concurrency.workers;

	DeclarePcoaOutputSchema(data->sample_id_type, return_types, names, /*with_batch_diagnostics=*/true);
	return std::move(data);
}

// ── Driving a progressive run from the scan ──────────────────────────────────────
// Every progressive_pcoa_from_* shares this: the run is constructed and stepped here,
// one wave at a time, so a wave's rows reach DuckDB while the next wave is still
// being computed and only one wave of rows is ever buffered.
//
// WHY not in Bind, where it used to be: Bind can only hand over a finished result,
// which meant the coordinates existed twice at once (the core's vector plus the
// row vector built from it — ~1.3 GB at 1M samples × d=10, ~13 GB at 10M), nothing
// was emitted until the last batch finished, and there was no point at which a
// multi-hour run could notice it had been cancelled.
struct ProgressivePcoaGlobalState : public GlobalTableFunctionState {
	// Null until the first Execute. The run owns its block source (the provider
	// closures hold it), so this is also what bounds the source's lifetime.
	unique_ptr<miint::progressive::ProgressivePcoaRun> run;
	// The current wave's rows and how far Execute has paged into them. One wave, not
	// the run — that is the whole point. With global_rotation on this holds one
	// STAGED chunk's rows instead of one wave's; either way it is bounded and
	// EmitPcoaChunk drains it the same way.
	std::vector<PcoaRow> rows;
	size_t cursor = 0;
	LogicalType sample_id_type = LogicalType::VARCHAR;

	// ── global_rotation staging ──────────────────────────────────────────────
	// The whole configuration, one row per sample with d coordinate columns, held in
	// a buffer-managed ColumnDataCollection rather than extension heap: that draws
	// against memory_limit and SPILLS to temp_directory, so a run that no longer fits
	// gets slower instead of dying. It has to be constructed from the BufferManager to
	// get that — see StageRotatedRun. The blocks themselves are still fetched and released
	// per batch — what has to be held is the finished coordinates, because applying
	// the rotation means revisiting rows that would otherwise already be gone.
	unique_ptr<ColumnDataCollection> staged;
	ColumnDataScanState staged_scan;
	DataChunk staged_chunk;
	bool staged_scan_started = false;
	unique_ptr<miint::progressive::PrincipalAxisAccumulator> axes;
	std::vector<double> rotation; // d*d row-major
	std::vector<double> centroid; // d
	// Copied out of the run so it can be destroyed before the emit phase, which frees
	// the anchor row cache and distance-corner cache for the whole of it.
	std::vector<double> eigvals;
	std::vector<double> proportions;
	// The run is stepped from one place, in order; a parallel scan would interleave
	// waves and reorder rows.
	idx_t MaxThreads() const override {
		return 1;
	}
};

// Convert one batch's coordinates into output rows. `anchor_m2` is that batch's
// anchor-overlap disparity; it is ignored for the anchor coordinates themselves
// (batch < 0), which report NULL because they define the frame rather than being
// fitted into it.
void AppendPcoaRows(const std::vector<miint::progressive::ProgressiveCoord> &coords, const std::vector<double> &eigvals,
                    const std::vector<double> &proportions, double anchor_m2, std::vector<PcoaRow> &out) {
	out.reserve(out.size() + coords.size());
	for (const auto &c : coords) {
		const auto axis = static_cast<uint32_t>(c.axis);
		PcoaRow row;
		row.iteration = 0; // kept for schema parity with pcoa / unifrac_pcoa
		row.sample_id = c.sample_id;
		row.axis = c.axis;
		row.coordinate = c.coordinate;
		row.eigenvalue = eigvals[axis];
		row.proportion_explained = proportions[axis];
		row.batch = c.batch;
		if (c.batch >= 0) {
			row.batch_anchor_m2 = anchor_m2;
		}
		out.push_back(std::move(row));
	}
}

// Wave sizing and run construction for a provider that COMPUTES each block inside
// its own batch and releases it there — the UniFrac and community-distance paths
// both do. Nothing of a wave is held but the output of the batches already in it,
// so `bytes_per_batch` is all there is to charge.
//
// The wave is sized as WIDE as its buffered output can afford rather than to the
// worker count, because a wave is a BARRIER and one-batch-per-worker is its worst
// case: every barrier waits for the wave's slowest batch, so the run pays max(batch)
// per wave instead of mean(batch), and the ragged final wave leaves the pool idle to
// the end. Widening amortizes both; a wave wide enough to hold the whole run has no
// barrier at all. Width is an I/O choice only — it never changes a coordinate.
//
// A continuous cross-wave work queue was priced and rejected: at a wave of
// budget/bytes_per_batch batches the barrier costs about (max-mean)/width per wave,
// ~2% at defaults, which does not pay for a producer/consumer pool and its
// cancellation and lifetime edges. The 64 MB budget is what keeps this from being
// the old "materialize the whole run" bug on a big table — well under the memory the
// workers' own blocks already hold, and the floor of `workers` batches is exactly
// what this path buffered before. Per-path measurements are at the call sites.
//
// Shared rather than copied because the two paths must stay in step: a wave width
// that drifted between them would silently give one path a different memory
// profile from the other for no stated reason.
unique_ptr<miint::progressive::ProgressivePcoaRun>
MakeComputedBlockRun(const ProgressivePcoaData &data, const miint::progressive::BlockProvider &provider,
                     const miint::progressive::InterruptCheck &interrupt) {
	// A batch's coordinates are charged twice over because they exist twice while
	// a wave is drained: the core's ProgressiveCoord and the PcoaRow copied from it
	// (see AppendPcoaRows).
	const uint64_t bytes_per_batch = static_cast<uint64_t>(data.batch_size) * data.n_dims *
	                                 (sizeof(miint::progressive::ProgressiveCoord) + sizeof(PcoaRow));
	const size_t n_batches = (data.part.remaining.size() + data.batch_size - 1) / data.batch_size;
	const uint32_t wave_batches = miint::progressive::ChooseWaveWidthByOutput(n_batches, data.workers, bytes_per_batch,
	                                                                          /*budget_bytes=*/64ull << 20);
	return make_uniq<miint::progressive::ProgressivePcoaRun>(
	    data.part.anchors, data.part.remaining, data.n_dims, data.batch_size, data.seed, data.n_threads, provider,
	    /*prefetch=*/nullptr, wave_batches, /*batch_workers=*/data.workers, interrupt);
}

// Build the run and its block source. Called once per scan, at execution time,
// because the source issues queries on the context that is running the scan.
unique_ptr<miint::progressive::ProgressivePcoaRun> MakeProgressiveRun(ClientContext &context,
                                                                      const ProgressivePcoaData &data) {
	// Cooperative cancellation. The core polls this before every batch, which is the
	// only thing standing between a user and an uninterruptible multi-hour query:
	// Ctrl-C sets context.interrupted and, until now, nothing on this path ever read
	// it. Polled from worker threads too, hence the atomic read.
	const miint::progressive::InterruptCheck interrupt = [&context]() {
		if (context.interrupted) {
			throw InterruptException();
		}
	};

	if (data.source == ProgressivePcoaData::Source::UNIFRAC) {
		// Shared, not owned by one closure: the cache is read-only after construction
		// (its anchor rows are loaded eagerly, before any fan-out) so every concurrent
		// block reads the same one. The run holds the provider, so the provider owning
		// the cache is what keeps it alive exactly as long as the run.
		auto cache = std::make_shared<AnchorFeatureRowCache>(context, data.CallerName(), data.qname, data.part.anchors);
		const miint::progressive::BlockProvider provider = [&context, &data,
		                                                    cache](const std::vector<std::string> &requested) {
			return ComputeUnifracBlock(context, data.qname, requested, *data.tree, data.variant_fp32,
			                           data.variance_adjust, data.alpha, data.bypass_tips, data.normalize_sample_counts,
			                           data.seed, data.block_threads, *cache);
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
		//
		// The wave is sized as WIDE as its buffered output can afford, not to the
		// worker count, because a wave is a barrier and one-batch-per-worker is the
		// worst way to draw it. Instrumented on 23,814 EMP samples, 14 threads, 23
		// batches: a batch is 94.6% UniFrac compute (the fsvd is 14 ms; the row read,
		// biom build and shear together are 5%) and genuinely uneven — 4.4 s to 9.5 s,
		// max/mean 1.19 — so wave 1's fourteen batches finished between 8.5 s and
		// 12.1 s and wave 2 could not start until 12.1 s. Widened, a worker takes the
		// next batch the moment it frees (batch 14 started at 7.959 s, when batch 13
		// ended at 7.958 s). A/B'd against a stashed build, interleaved, 4 pairs:
		// 20.89 s -> 19.39 s, same checksum to 6 dp over all 238,140 coordinates —
		// wave width has never changed a coordinate, and this confirms it end to end.
		//
		// The gain is smaller than the idle workers suggest because keeping the pool
		// busy makes each block dearer: in the instrumented run the same 2000-sample
		// block averaged 8.45 s while fourteen ran at once and 6.23 s while nine did.
		// The batch phase is memory-bandwidth-bound, so filling an idle worker returns
		// less than the arithmetic promises — the 1→14 thread curve bends for the same
		// reason, 8 → 14 threads buying only 1.04× in this build AND in the one before
		// it. What is left is quantization — 23 batches over 14 workers is two rounds
		// however they are scheduled, and no scheduler fixes that; it fades to under 1%
		// by a few hundred batches. Two things that DO close it were measured and left
		// to the caller: `threads := 28` on this 14-core box ran 17.2 s (a block is one
		// core, so oversubscribing lets the OS split the tail) but would break `threads`
		// as a bound on a machine with more cores than the budget. Re-cutting batch_size
		// to make the batch count a multiple of the pool did NOT help once the barrier
		// was gone (28 batches: 18.6 s vs 18.9 s) — it only helped the barrier build, so
		// do not re-derive it. Judge all of this on WALL time: the CLI timer's CPU column
		// reported 284 vs 141 CPU-seconds for two runs that execute identically.
		//
		// A continuous cross-wave work queue was priced and rejected: at a wave of
		// `budget/bytes_per_batch` batches the barrier costs about (max-mean)/width per
		// wave, ~2% at defaults, which does not pay for a producer/consumer thread pool
		// and its cancellation and lifetime edges. The budget is what keeps this from
		// being the old "materialize the whole run" bug (#11) on a big table: 64 MB of
		// coordinates, well under the ~280 MB of blocks the workers already hold, and
		// the floor of `workers` batches is exactly what this path buffered before.
		//
		// Bit-identity is preserved, not merely approximate: each worker pins its
		// ordination to ONE OpenMP thread, which is what the serial path uses at
		// n_threads=1, and skbb's centering reduction sums in thread-count-dependent
		// order (see test_ProgressivePcoa's serial-vs-parallel case).
		//
		return MakeComputedBlockRun(data, provider, interrupt);
	}

	if (data.source == ProgressivePcoaData::Source::FEATURES) {
		// Same shape as the UniFrac path above and for the same reasons: the anchor
		// feature rows are read once and shared read-only across concurrent blocks,
		// and the run owns the provider, which owns the cache.
		//
		// Blocks run concurrently here too. A block is an O(m^2 * f) pair loop rather
		// than an fsvd, so — as with UniFrac — the machine is used by running several
		// blocks at once; `data.block_threads` is the per-block share bind already
		// divided out of `threads`, so the two levels cannot oversubscribe.
		auto cache = std::make_shared<AnchorFeatureRowCache>(context, data.CallerName(), data.qname, data.part.anchors);
		// Built here, on the calling thread, before any fan-out — like the row cache
		// above it is read-only afterwards, which is what lets every concurrent block
		// share one. It costs one a x a distance computation up front and removes the
		// same computation from every block of the run.
		// data.n_threads, not data.block_threads: nothing else is running yet, so the
		// per-block share would leave the machine idle through the whole serial phase.
		auto corner = std::make_shared<AnchorCommunityCorner>(cache->anchor_rows(), data.part.anchors, data.metric,
		                                                      data.n_threads);
		const miint::progressive::BlockProvider provider = [&context, &data, cache,
		                                                    corner](const std::vector<std::string> &requested) {
			return ComputeCommunityBlock(context, data.qname, requested, data.metric, data.block_threads, *cache,
			                             *corner, data.gram_operand_bytes);
		};
		return MakeComputedBlockRun(data, provider, interrupt);
	}

	auto source = std::make_shared<WaveDistanceBlockSource>(context, data.qname, data.part.anchors);
	const miint::progressive::BlockProvider provider = [source](const std::vector<std::string> &requested) {
		return source->Get(requested);
	};
	const miint::progressive::WavePrefetch prefetch = [source](const std::vector<std::vector<std::string>> &requests) {
		source->Prefetch(requests);
	};

	// How many batches to serve per relation scan (see ChooseWaveWidth for what a
	// wave costs). The budget is a quarter of what memory_limit currently leaves —
	// a quarter because the blocks are extension heap the buffer manager can neither
	// see nor evict, the same reason pcoa()'s dense matrix is guarded rather than
	// tracked. The scan rows ARE buffer-managed (see RunWaveQuery), so overshooting
	// the estimate spills rather than dies; they are still charged, because spilling
	// a wave's scan result costs far more than running a narrower wave.
	auto &buffer_manager = BufferManager::GetBufferManager(context);
	const auto max_memory = buffer_manager.GetMaxMemory();
	const auto used_memory = buffer_manager.GetUsedMemory();
	const uint64_t budget = (max_memory > used_memory ? max_memory - used_memory : 0) / 4;
	const size_t n_batches = (data.part.remaining.size() + data.batch_size - 1) / data.batch_size;
	const uint32_t wave_batches = miint::progressive::ChooseWaveWidth(
	    data.part.anchors.size(), data.batch_size, static_cast<uint32_t>(data.n_threads), budget, n_batches);

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
	// from_unifrac variant, where a block IS a UniFrac compute, is where it pays.
	return make_uniq<miint::progressive::ProgressivePcoaRun>(
	    data.part.anchors, data.part.remaining, data.n_dims, data.batch_size, data.seed, data.n_threads, provider,
	    prefetch, wave_batches, static_cast<uint32_t>(data.n_threads), interrupt);
}

// Refill `rows` with the next piece of the run: the anchor coordinates first (they
// are the reference frame, and the whole batch phase depends on them), then one
// wave of batches per call. Returns false once the run is complete.
bool AdvanceProgressiveRun(ClientContext &context, const ProgressivePcoaData &data,
                           ProgressivePcoaGlobalState &gstate) {
	gstate.rows.clear();
	gstate.cursor = 0;
	try {
		if (!gstate.run) {
			gstate.run = MakeProgressiveRun(context, data);
			AppendPcoaRows(gstate.run->Start(), gstate.run->eigvals(), gstate.run->proportion_explained(),
			               /*anchor_m2=*/0.0, gstate.rows);
			return true;
		}
		std::vector<miint::progressive::BatchOutput> wave = gstate.run->NextWave();
		// Sized for the whole wave before appending any of it: reserving per batch
		// instead reallocates once per batch and copies everything appended so far,
		// which is quadratic in the wave's width (and a wave can be over a hundred
		// batches wide on the _distances path). Capacity survives the clear() above, so
		// only the first wave pays even this.
		size_t wave_rows = 0;
		for (const auto &out : wave) {
			wave_rows += out.coords.size();
		}
		gstate.rows.reserve(wave_rows);
		for (const auto &out : wave) {
			AppendPcoaRows(out.coords, gstate.run->eigvals(), gstate.run->proportion_explained(), out.diag.anchor_m2,
			               gstate.rows);
		}
	} catch (const std::invalid_argument &e) {
		// The core's own guards (see RunProgressivePcoa) — a bad anchor set, a provider
		// that returned the wrong samples. Bind pre-validates the parameter forms, so
		// what reaches here is about the data.
		throw InvalidInputException("%s: %s", data.CallerName(), e.what());
	}
	// A wave that ran any batch always yields rows — a batch has at least one sample
	// and d >= 1 axes — so an empty refill means NextWave found nothing left to run.
	return !gstate.rows.empty();
}

// Types of the staging collection: the sample id as the core produced it (a string;
// EmitIdCell converts to the caller's id type at emit), the batch diagnostics, and
// one DOUBLE per axis. Wide rather than long — one row per sample instead of d —
// because rotating needs a sample's d coordinates together, and a wide row cannot be
// split across a chunk boundary.
vector<LogicalType> StagedRowTypes(uint32_t d) {
	vector<LogicalType> types {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::DOUBLE};
	for (uint32_t a = 0; a < d; ++a) {
		types.push_back(LogicalType::DOUBLE);
	}
	return types;
}

// Fold one wave's long-format rows into staged wide rows, accumulating the second
// moments as they pass. Every sample's d axis rows are consecutive and complete
// within a wave (a wave emits whole batches, and AppendPcoaRows walks coords in
// sample-major order), which is what lets this gather d at a time without a carry
// buffer between calls — asserted rather than assumed.
void StageWave(const ProgressivePcoaData &data, ProgressivePcoaGlobalState &gstate, DataChunk &append_chunk) {
	const auto d = data.n_dims;
	if (gstate.rows.size() % d != 0) {
		throw InternalException("%s: staged %llu rows for d = %u; a sample's axes must arrive together",
		                        data.CallerName(), static_cast<unsigned long long>(gstate.rows.size()), d);
	}
	std::vector<double> coords(d);
	for (size_t r = 0; r < gstate.rows.size(); r += d) {
		if (append_chunk.size() == STANDARD_VECTOR_SIZE) {
			gstate.staged->Append(append_chunk);
			append_chunk.Reset();
		}
		const idx_t row = append_chunk.size();
		const auto &first = gstate.rows[r];
		append_chunk.SetValue(0, row, Value(first.sample_id));
		if (first.batch < 0) {
			append_chunk.SetValue(1, row, Value(LogicalType::INTEGER));
			append_chunk.SetValue(2, row, Value(LogicalType::DOUBLE));
		} else {
			append_chunk.SetValue(1, row, Value::INTEGER(first.batch));
			append_chunk.SetValue(2, row, Value::DOUBLE(first.batch_anchor_m2));
		}
		for (uint32_t a = 0; a < d; ++a) {
			const auto &pr = gstate.rows[r + a];
			if (static_cast<uint32_t>(pr.axis) != a || pr.sample_id != first.sample_id) {
				throw InternalException("%s: staged rows for '%s' are not axis-ordered", data.CallerName(),
				                        first.sample_id.c_str());
			}
			coords[a] = pr.coordinate;
			append_chunk.SetValue(3 + a, row, Value::DOUBLE(pr.coordinate));
		}
		append_chunk.SetCardinality(row + 1);
		gstate.axes->Add(coords.data(), 1);
	}
}

// Run the WHOLE run, staging it, then solve for the rotation. This is where
// global_rotation costs streaming: nothing can be emitted until the last batch has
// been placed, because the transform is a property of the finished configuration.
void StageRotatedRun(ClientContext &context, const ProgressivePcoaData &data, ProgressivePcoaGlobalState &gstate) {
	const auto types = StagedRowTypes(data.n_dims);
	// BufferManager, not Allocator: the Allocator overload selects
	// IN_MEMORY_ALLOCATOR, which is plain heap — uncharged against memory_limit and
	// unable to spill. Since this holds the WHOLE run, that is the difference between
	// a large run getting slower and the process being OOM-killed.
	gstate.staged = make_uniq<ColumnDataCollection>(BufferManager::GetBufferManager(context), types);
	gstate.axes = make_uniq<miint::progressive::PrincipalAxisAccumulator>(data.n_dims);
	DataChunk append_chunk;
	append_chunk.Initialize(Allocator::Get(context), types);

	while (AdvanceProgressiveRun(context, data, gstate)) {
		if (context.interrupted) {
			throw InterruptException();
		}
		StageWave(data, gstate, append_chunk);
	}
	if (append_chunk.size() > 0) {
		gstate.staged->Append(append_chunk);
	}
	gstate.rows.clear();
	gstate.cursor = 0;

	// Held past the run so the emit phase can report them, and captured HERE so the
	// run — and with it the anchor row cache and the anchor-distance corner — can be
	// released before the emit phase allocates anything.
	gstate.eigvals = gstate.run->eigvals();
	gstate.proportions = gstate.run->proportion_explained();
	try {
		gstate.rotation = gstate.axes->Rotation();
		gstate.centroid = gstate.axes->Mean();
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s: %s", data.CallerName(), e.what());
	}
	gstate.axes.reset();
	gstate.run.reset();
}

// One staged chunk, rotated, expanded back into the long-format rows EmitPcoaChunk
// already knows how to drain.
bool RefillFromStaged(const ProgressivePcoaData &data, ProgressivePcoaGlobalState &gstate) {
	const auto d = data.n_dims;
	if (!gstate.staged_scan_started) {
		gstate.staged->InitializeScan(gstate.staged_scan);
		gstate.staged->InitializeScanChunk(gstate.staged_chunk);
		gstate.staged_scan_started = true;
	}
	gstate.rows.clear();
	gstate.cursor = 0;
	if (!gstate.staged->Scan(gstate.staged_scan, gstate.staged_chunk)) {
		return false;
	}
	const idx_t n = gstate.staged_chunk.size();
	gstate.rows.reserve(static_cast<size_t>(n) * d);
	std::vector<double> raw(d), rotated(d);
	for (idx_t r = 0; r < n; ++r) {
		const auto id_val = gstate.staged_chunk.data[0].GetValue(r);
		const auto batch_val = gstate.staged_chunk.data[1].GetValue(r);
		const auto m2_val = gstate.staged_chunk.data[2].GetValue(r);
		for (uint32_t a = 0; a < d; ++a) {
			raw[a] = gstate.staged_chunk.data[3 + a].GetValue(r).GetValue<double>();
		}
		// y' = R (y - centroid): centred as well as rotated, so the emitted
		// configuration is a centred principal-axis one like a real PCoA's.
		for (uint32_t a = 0; a < d; ++a) {
			double acc = 0.0;
			for (uint32_t k = 0; k < d; ++k) {
				acc += gstate.rotation[static_cast<size_t>(a) * d + k] * (raw[k] - gstate.centroid[k]);
			}
			rotated[a] = acc;
		}
		for (uint32_t a = 0; a < d; ++a) {
			PcoaRow row;
			row.iteration = 0;
			row.sample_id = id_val.ToString();
			row.axis = static_cast<int32_t>(a);
			row.coordinate = rotated[a];
			row.eigenvalue = gstate.eigvals[a];
			row.proportion_explained = gstate.proportions[a];
			row.batch = batch_val.IsNull() ? -1 : batch_val.GetValue<int32_t>();
			if (!m2_val.IsNull()) {
				row.batch_anchor_m2 = m2_val.GetValue<double>();
			}
			gstate.rows.push_back(std::move(row));
		}
	}
	return true;
}

unique_ptr<GlobalTableFunctionState> ProgressivePcoaInitGlobal(ClientContext &, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<ProgressivePcoaData>();
	auto gstate = make_uniq<ProgressivePcoaGlobalState>();
	gstate->sample_id_type = data.sample_id_type;
	return std::move(gstate);
}

void ProgressivePcoaExecute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<ProgressivePcoaGlobalState>();
	const auto &data = input.bind_data->Cast<ProgressivePcoaData>();
	if (data.global_rotation) {
		// The whole run happens inside this first call — the rotation is a property of
		// the finished configuration, so there is nothing correct to emit before then.
		if (!gstate.staged) {
			// AdvanceProgressiveRun builds the run and emits the anchor frame on its
			// first call; creating it here instead would skip that and lose the anchors.
			StageRotatedRun(context, data, gstate);
		}
		while (gstate.cursor >= gstate.rows.size()) {
			if (context.interrupted) {
				throw InterruptException();
			}
			if (!RefillFromStaged(data, gstate)) {
				output.SetCardinality(0);
				return;
			}
		}
		EmitPcoaChunk(gstate.rows, gstate.cursor, gstate.sample_id_type, /*with_batch_diagnostics=*/true, output);
		return;
	}
	while (gstate.cursor >= gstate.rows.size()) {
		// Also polled here, not only inside the run: a cancellation that arrives while
		// DuckDB is draining the wave already produced is noticed on the next refill
		// rather than after another wave's worth of work.
		if (context.interrupted) {
			throw InterruptException();
		}
		if (!AdvanceProgressiveRun(context, data, gstate)) {
			output.SetCardinality(0);
			return;
		}
	}
	EmitPcoaChunk(gstate.rows, gstate.cursor, gstate.sample_id_type, /*with_batch_diagnostics=*/true, output);
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
	fn.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(fn);
}

void RegisterProgressivePcoaFromDistances(ExtensionLoader &loader) {
	TableFunction fn("progressive_pcoa_from_distances", {LogicalType::VARCHAR}, ProgressivePcoaExecute,
	                 ProgressivePcoaFromDistancesBind, ProgressivePcoaInitGlobal);
	fn.named_parameters["n_dims"] = LogicalType::INTEGER;
	fn.named_parameters["n_anchors"] = LogicalType::INTEGER;
	fn.named_parameters["batch_size"] = LogicalType::INTEGER;
	fn.named_parameters["seed"] = LogicalType::INTEGER;
	fn.named_parameters["threads"] = LogicalType::INTEGER;
	fn.named_parameters["anchors"] = LogicalType::LIST(LogicalType::VARCHAR);
	fn.named_parameters["global_rotation"] = LogicalType::BOOLEAN;
	loader.RegisterFunction(fn);
}

void RegisterProgressivePcoaFromFeatures(ExtensionLoader &loader) {
	TableFunction fn("progressive_pcoa_from_features", {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                 ProgressivePcoaExecute, ProgressivePcoaFromFeaturesBind, ProgressivePcoaInitGlobal);
	fn.named_parameters["n_dims"] = LogicalType::INTEGER;
	fn.named_parameters["n_anchors"] = LogicalType::INTEGER;
	fn.named_parameters["batch_size"] = LogicalType::INTEGER;
	fn.named_parameters["seed"] = LogicalType::INTEGER;
	fn.named_parameters["threads"] = LogicalType::INTEGER;
	fn.named_parameters["anchors"] = LogicalType::LIST(LogicalType::VARCHAR);
	fn.named_parameters["global_rotation"] = LogicalType::BOOLEAN;
	loader.RegisterFunction(fn);
}

void RegisterProgressivePcoaFromUnifrac(ExtensionLoader &loader) {
	TableFunction fn("progressive_pcoa_from_unifrac", {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                 ProgressivePcoaExecute, ProgressivePcoaFromUnifracBind, ProgressivePcoaInitGlobal);
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
	fn.named_parameters["global_rotation"] = LogicalType::BOOLEAN;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
