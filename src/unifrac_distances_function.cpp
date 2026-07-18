#include "unifrac_table_functions.hpp"

#include <climits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "NewickTree.hpp"
#include "id_column_utils.hpp"
#include "tree_table_reader.hpp"
#include "unifrac_bptree.hpp"
#include "unifrac_condensed.hpp"
#include "unifrac_distance.hpp"
#include "unifrac_function_common.hpp"
#include "unifrac_support_biom.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"

namespace duckdb {
namespace {

using unifrac_internal::AcceptedVariantList;
using unifrac_internal::IsValidVariant;
using unifrac_internal::ReadFeatureTable;
using unifrac_internal::ResolveThreadsParameter;

// Bind data owns the biom + bptree views for the whole scan. Both view types
// have self-wiring move ctors, and this struct is heap-allocated (make_uniq)
// and never moved after construction, so the char** pointers libssu reads stay
// valid across the Bind → Execute boundary. The distance matrix is deliberately
// NOT computed here: unifrac_distances streams the upper triangle one iteration
// at a time (see the global state) rather than materializing all O(n²) rows.
struct UnifracDistancesData : public TableFunctionData {
	miint::unifrac::UnifracSupportBiomView biom_view;
	miint::unifrac::UnifracBptreeView bptree_view;
	std::string variant_fp32;
	bool variance_adjust = false;
	double alpha = 1.0;
	bool bypass_tips = false;
	bool normalize_sample_counts = true;
	uint32_t subsample_depth = 0;
	bool subsample_with_replacement = false;
	int32_t n_subsamples = 1;
	int32_t seed = -1;
	int n_threads = 1;
	// Output type for sample_a/sample_b — mirrors the input sample_id type
	// (BIGINT/UUID) or VARCHAR otherwise. See ResolveSampleIdOutputType.
	LogicalType sample_id_out_type = LogicalType::VARCHAR;

	UnifracDistancesData(miint::unifrac::UnifracSupportBiomView biom, miint::unifrac::UnifracBptreeView bptree)
	    : biom_view(std::move(biom)), bptree_view(std::move(bptree)) {
	}
};

// Holds the current iteration's distance matrix and a cursor into its upper
// triangle. Single-threaded (MaxThreads == 1): DuckDB drives one Execute stream
// which pages the triangle out in STANDARD_VECTOR_SIZE chunks and, when an
// iteration is exhausted, recomputes the next one in place.
struct UnifracDistancesGlobalState : public GlobalTableFunctionState {
	unique_ptr<miint::unifrac::UnifracDistanceMatrix> matrix;
	miint::unifrac::CondensedCursor cursor;
	int32_t iteration = 0;
	bool finished = false;

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

// Compute the distance matrix for `iteration` into the global state and reset
// the cursor to the start of its upper triangle. Subsampling can drop samples
// below subsample_depth, so each iteration's matrix may have a different
// n_samples; the cursor is always sized from the matrix that was produced.
void LoadIteration(const UnifracDistancesData &data, UnifracDistancesGlobalState &gstate, int32_t iteration) {
	// seed_iter mirrors pcoa/permanova: seed + iteration, with -1 meaning
	// "unseeded". The seed + n_subsamples overflow is guarded at bind time.
	const int seed_iter = (data.seed >= 0) ? (data.seed + iteration) : -1;
	// UnifracDistanceMatrix::Compute throws std::runtime_error on libssu errors;
	// tree/feature mismatch is already caught at bind by ValidateTreeCoversFeatures,
	// so anything here is a libssu-internal failure surfaced as InvalidInput.
	gstate.matrix = make_uniq<miint::unifrac::UnifracDistanceMatrix>([&]() {
		try {
			return miint::unifrac::UnifracDistanceMatrix::Compute(
			    data.biom_view, data.bptree_view, data.variant_fp32, data.variance_adjust, data.alpha, data.bypass_tips,
			    data.normalize_sample_counts, data.subsample_depth, data.subsample_with_replacement, seed_iter,
			    data.n_threads);
		} catch (const std::runtime_error &e) {
			throw InvalidInputException("unifrac_distances: %s", e.what());
		}
	}());
	gstate.iteration = iteration;
	gstate.cursor = miint::unifrac::CondensedCursor::Begin(gstate.matrix->n_samples());
}

unique_ptr<FunctionData> UnifracDistancesBind(ClientContext &context, TableFunctionBindInput &input,
                                              vector<LogicalType> &return_types, vector<string> &names) {
	const std::string table_name = input.inputs[0].GetValue<string>();
	const std::string tree_name = input.inputs[1].GetValue<string>();
	if (table_name.empty()) {
		throw BinderException("unifrac_distances: feature-table name must not be empty");
	}
	if (tree_name.empty()) {
		throw BinderException("unifrac_distances: tree name must not be empty");
	}

	// INTEGER-typed named parameters are read as int32_t to match
	// LogicalType::INTEGER (align_common.hpp:75 convention).
	std::string variant = "weighted_normalized";
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
	const int n_threads = ResolveThreadsParameter(context, threads, "unifrac_distances");

	if (!IsValidVariant(variant)) {
		throw BinderException("unifrac_distances: variant '%s' is not valid. Accepted variants: %s", variant,
		                      AcceptedVariantList());
	}
	if (n_subsamples < 1) {
		throw BinderException("unifrac_distances: n_subsamples must be >= 1 (got %d)", n_subsamples);
	}
	if (subsample_depth < 0) {
		throw BinderException("unifrac_distances: subsample_depth must be >= 0 (got %d)", subsample_depth);
	}
	if (n_subsamples > 1 && subsample_depth == 0) {
		throw BinderException("unifrac_distances: n_subsamples > 1 requires subsample_depth > 0 (iterations would "
		                      "otherwise be identical)");
	}
	// Guard the seed_iter = seed + i arithmetic against signed int32 overflow;
	// libssu's ssu_set_random_seed takes a 32-bit int.
	if (seed >= 0 &&
	    static_cast<int64_t>(seed) + static_cast<int64_t>(n_subsamples) - 1 > static_cast<int64_t>(INT_MAX)) {
		throw BinderException("unifrac_distances: seed (%d) + n_subsamples (%d) - 1 exceeds INT_MAX; pick a smaller "
		                      "seed or fewer subsamples",
		                      seed, n_subsamples);
	}

	LogicalType sample_id_type = LogicalType::VARCHAR;
	auto coo_rows = ReadFeatureTable(context, table_name, "unifrac_distances", &sample_id_type);
	if (coo_rows.empty()) {
		throw InvalidInputException("unifrac_distances: feature-table '%s' is empty after dropping NULL/zero rows",
		                            table_name);
	}
	miint::unifrac::UnifracSupportBiomView biom_view = [&]() {
		try {
			return miint::unifrac::UnifracSupportBiomView::FromCoo(std::move(coo_rows));
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("unifrac_distances: %s", e.what());
		}
	}();

	const auto *biom_struct = biom_view.support_biom();
	auto feature_ids = CollectIds(biom_struct->obs_ids, biom_struct->n_obs);

	auto tree_inputs = ReadTreeTable(context, tree_name);
	auto tree = miint::NewickTree::build(tree_inputs);
	try {
		miint::unifrac::ValidateTreeCoversFeatures(tree, feature_ids);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("unifrac_distances: %s", e.what());
	}
	auto bptree_view = miint::unifrac::UnifracBptreeView::FromNewickTree(tree);

	auto data = make_uniq<UnifracDistancesData>(std::move(biom_view), std::move(bptree_view));
	// libssu accepts both bare and _fp32-suffixed variants; append _fp32 so the
	// choice is pinned to fp32 regardless of libssu's default (see pcoa).
	data->variant_fp32 = variant + "_fp32";
	data->variance_adjust = variance_adjust;
	data->alpha = alpha;
	data->bypass_tips = bypass_tips;
	data->normalize_sample_counts = normalize_sample_counts;
	data->subsample_depth = static_cast<uint32_t>(subsample_depth);
	data->subsample_with_replacement = subsample_with_replacement;
	data->n_subsamples = n_subsamples;
	data->seed = seed;
	data->n_threads = n_threads;
	data->sample_id_out_type = unifrac_internal::ResolveSampleIdOutputType(sample_id_type);

	names.emplace_back("iteration");
	return_types.emplace_back(LogicalType::INTEGER);
	names.emplace_back("sample_a");
	return_types.emplace_back(data->sample_id_out_type);
	names.emplace_back("sample_b");
	return_types.emplace_back(data->sample_id_out_type);
	names.emplace_back("distance");
	return_types.emplace_back(LogicalType::DOUBLE);

	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> UnifracDistancesInitGlobal(ClientContext &, TableFunctionInitInput &input) {
	const auto &data = input.bind_data->Cast<UnifracDistancesData>();
	auto gstate = make_uniq<UnifracDistancesGlobalState>();
	LoadIteration(data, *gstate, /*iteration*/ 0);
	return std::move(gstate);
}

void UnifracDistancesExecute(ClientContext &, TableFunctionInput &input, DataChunk &output) {
	const auto &data = input.bind_data->Cast<UnifracDistancesData>();
	auto &gstate = input.global_state->Cast<UnifracDistancesGlobalState>();

	if (gstate.finished) {
		output.SetCardinality(0);
		return;
	}

	auto iter_data = FlatVector::GetData<int32_t>(output.data[0]);
	auto &sample_a_vec = output.data[1];
	auto &sample_b_vec = output.data[2];
	auto dist_data = FlatVector::GetData<double>(output.data[3]);
	const LogicalType &id_type = data.sample_id_out_type;

	idx_t out_n = 0;
	while (out_n < STANDARD_VECTOR_SIZE) {
		// Skip past exhausted iterations, loading the next matrix until one has a
		// pair or we run out. An iteration with fewer than 2 samples (including
		// the base case of a <2-sample table, or subsampling dropping an
		// iteration below 2) has an empty upper triangle, so it contributes no
		// rows. Unlike unifrac_pcoa — which errors when too few samples survive
		// because an ordination is undefined — the condensed form of such a
		// matrix is legitimately empty, so we emit no rows rather than raising.
		while (gstate.cursor.done()) {
			if (gstate.iteration + 1 >= data.n_subsamples) {
				gstate.finished = true;
				break;
			}
			LoadIteration(data, gstate, gstate.iteration + 1);
		}
		if (gstate.finished) {
			break;
		}

		// Hoist the matrix handles out of the per-row loop — they only change
		// when LoadIteration advances to the next iteration's matrix.
		const float *mat = gstate.matrix->matrix();
		const uint32_t n = gstate.matrix->n_samples();
		const auto &ids = gstate.matrix->sample_ids();
		while (out_n < STANDARD_VECTOR_SIZE && !gstate.cursor.done()) {
			const uint32_t i = gstate.cursor.i;
			const uint32_t j = gstate.cursor.j;
			iter_data[out_n] = gstate.iteration;
			// EmitIdCell mirrors the id type onto the output. Its ""/"*"→NULL
			// sentinel branch (SAM-domain) is unreachable here: ReadFeatureTable
			// drops NULL sample_ids and libssu sample ids are never empty.
			EmitIdCell(sample_a_vec, out_n, ids[i], id_type);
			EmitIdCell(sample_b_vec, out_n, ids[j], id_type);
			dist_data[out_n] = static_cast<double>(mat[static_cast<size_t>(i) * n + j]);
			++out_n;
			gstate.cursor.advance();
		}
	}

	output.SetCardinality(out_n);
}

} // namespace

void RegisterUnifracDistances(ExtensionLoader &loader) {
	TableFunction fn("unifrac_distances", {LogicalType::VARCHAR, LogicalType::VARCHAR}, UnifracDistancesExecute,
	                 UnifracDistancesBind, UnifracDistancesInitGlobal);
	fn.named_parameters["variant"] = LogicalType::VARCHAR;
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

} // namespace duckdb
