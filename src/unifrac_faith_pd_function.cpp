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
#include "unifrac_function_common.hpp"
#include "api.hpp"
#include "unifrac_omp_scope.hpp"
#include "unifrac_subsample_bridge.hpp"
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

using unifrac_internal::ReadFeatureTable;
using unifrac_internal::ResolveSampleIdOutputType;
using unifrac_internal::ResolveThreadsParameter;

struct FaithPdRow {
	int32_t iteration;
	std::string sample_id;
	double faith_pd;
};

struct UnifracFaithPdData : public TableFunctionData {
	std::vector<FaithPdRow> rows;
	// Output type for sample_id — mirrors the input sample_id type (BIGINT/UUID)
	// or VARCHAR otherwise. See ResolveSampleIdOutputType.
	LogicalType sample_id_type = LogicalType::VARCHAR;
};

struct UnifracFaithPdGlobalState : public GlobalTableFunctionState {
	std::vector<FaithPdRow> rows;
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

// RAII wrapper around libssu's r_vec*. The faith_pd_inmem path is
// straight-line (no early return on success), but out_rows.push_back may
// throw std::bad_alloc under memory pressure; without RAII the r_vec
// (its values/sample_ids arrays plus the struct) would leak per iteration
// and compound across the n_subsamples loop.
class ResultsVecHandle {
public:
	~ResultsVecHandle() {
		if (ptr_ != nullptr) {
			destroy_results_vec(&ptr_);
		}
	}
	r_vec **out() {
		return &ptr_;
	}
	r_vec *get() const {
		return ptr_;
	}

	ResultsVecHandle() = default;
	ResultsVecHandle(const ResultsVecHandle &) = delete;
	ResultsVecHandle &operator=(const ResultsVecHandle &) = delete;

private:
	r_vec *ptr_ = nullptr;
};

// Run faith_pd_inmem against a support_biom view and emit one row per
// sample. The view may be the original biom (subsample_depth == 0) or a
// bridged subsample (subsample_depth > 0). Note: sample order in the output
// is libssu's order (result->sample_ids), which for our use case is the
// lexicographic order established by UnifracSupportBiomView::FromCoo —
// downstream SQL queries should ORDER BY sample_id regardless rather than
// rely on this implementation detail.
void RunFaithPd(const miint::unifrac::UnifracSupportBiomView &biom_view,
                const miint::unifrac::UnifracBptreeView &bptree_view, int n_threads, int32_t iteration_index,
                std::vector<FaithPdRow> &out_rows) {
	ResultsVecHandle result;
	ComputeStatus status;
	{
		// faith_pd_inmem internally uses libssu's OpenMP-parallel tree
		// traversal; pin its fan-out to n_threads so it honors DuckDB's thread
		// count. A bare pin is enough here — faith_pd_inmem draws no randomness
		// and upstream lists it as safe to call concurrently.
		miint::unifrac::OmpThreadPin omp_pin(n_threads);
		status = faith_pd_inmem(biom_view.support_biom(), bptree_view.support_bptree(), result.out());
	}
	if (status != okay) {
		throw InvalidInputException("unifrac_faith_pd: libssu faith_pd_inmem returned status %d",
		                            static_cast<int>(status));
	}
	const r_vec *r = result.get();
	for (unsigned int i = 0; i < r->n_samples; ++i) {
		FaithPdRow row;
		row.iteration = iteration_index;
		row.sample_id = r->sample_ids[i];
		row.faith_pd = r->values[i];
		out_rows.push_back(std::move(row));
	}
}

unique_ptr<FunctionData> UnifracFaithPdBind(ClientContext &context, TableFunctionBindInput &input,
                                            vector<LogicalType> &return_types, vector<string> &names) {
	const std::string table_name = input.inputs[0].GetValue<string>();
	const std::string tree_name = input.inputs[1].GetValue<string>();
	if (table_name.empty()) {
		throw BinderException("unifrac_faith_pd: feature-table name must not be empty");
	}
	if (tree_name.empty()) {
		throw BinderException("unifrac_faith_pd: tree name must not be empty");
	}

	// INTEGER-typed named params match the project convention (int32_t).
	int32_t subsample_depth = 0;
	bool subsample_with_replacement = false;
	int32_t n_subsamples = 1;
	int32_t seed = -1;
	int32_t threads = 0; // 0 = follow DuckDB's TaskScheduler::NumberOfThreads()
	for (const auto &kv : input.named_parameters) {
		const auto key = StringUtil::Lower(kv.first);
		if (key == "subsample_depth") {
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
	const int n_threads = ResolveThreadsParameter(context, threads, "unifrac_faith_pd");

	if (n_subsamples < 1) {
		throw BinderException("unifrac_faith_pd: n_subsamples must be >= 1 (got %d)", n_subsamples);
	}
	if (subsample_depth < 0) {
		throw BinderException("unifrac_faith_pd: subsample_depth must be >= 0 (got %d)", subsample_depth);
	}
	if (n_subsamples > 1 && subsample_depth == 0) {
		throw BinderException("unifrac_faith_pd: n_subsamples > 1 requires subsample_depth > 0 (iterations would "
		                      "otherwise be identical)");
	}
	if (seed >= 0 &&
	    static_cast<int64_t>(seed) + static_cast<int64_t>(n_subsamples) - 1 > static_cast<int64_t>(INT_MAX)) {
		throw BinderException("unifrac_faith_pd: seed (%d) + n_subsamples (%d) - 1 exceeds INT_MAX; "
		                      "pick a smaller seed or fewer subsamples",
		                      seed, n_subsamples);
	}

	LogicalType sample_id_col_type = LogicalType::VARCHAR;
	auto coo_rows = ReadFeatureTable(context, table_name, "unifrac_faith_pd", &sample_id_col_type);
	if (coo_rows.empty()) {
		throw InvalidInputException("unifrac_faith_pd: feature-table '%s' is empty after dropping NULL/zero rows",
		                            table_name);
	}
	miint::unifrac::UnifracSupportBiomView biom_view = [&]() {
		try {
			return miint::unifrac::UnifracSupportBiomView::FromCoo(std::move(coo_rows));
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("unifrac_faith_pd: %s", e.what());
		}
	}();

	const auto *biom_struct = biom_view.support_biom();
	auto feature_ids = CollectIds(biom_struct->obs_ids, biom_struct->n_obs);
	if (biom_struct->n_samples < 1) {
		throw BinderException("unifrac_faith_pd: feature-table '%s' has no samples", table_name);
	}

	auto tree_inputs = ReadTreeTable(context, tree_name);
	auto tree = miint::NewickTree::build(tree_inputs);
	try {
		miint::unifrac::ValidateTreeCoversFeatures(tree, feature_ids);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("unifrac_faith_pd: %s", e.what());
	}
	auto bptree_view = miint::unifrac::UnifracBptreeView::FromNewickTree(tree);

	auto data = make_uniq<UnifracFaithPdData>();
	data->sample_id_type = ResolveSampleIdOutputType(sample_id_col_type);
	data->rows.reserve(static_cast<size_t>(n_subsamples) * biom_struct->n_samples);

	if (subsample_depth == 0) {
		// No subsampling: n_subsamples > 1 is already rejected above, so this
		// branch always runs exactly once with iteration_index=0.
		RunFaithPd(biom_view, bptree_view, n_threads, /*iteration_index*/ 0, data->rows);
	} else {
		const auto depth = static_cast<uint32_t>(subsample_depth);
		for (int32_t i = 0; i < n_subsamples; ++i) {
			// seed + i overflow prevented by bind-time check above.
			const int seed_iter = (seed >= 0) ? (seed + i) : -1;
			miint::unifrac::UnifracSupportBiomView sub = [&]() {
				try {
					return miint::unifrac::BridgeSubsample(biom_view, depth, subsample_with_replacement, seed_iter);
				} catch (const std::runtime_error &e) {
					throw InvalidInputException("unifrac_faith_pd: %s", e.what());
				} catch (const std::invalid_argument &e) {
					throw InvalidInputException("unifrac_faith_pd: %s", e.what());
				}
			}();
			RunFaithPd(sub, bptree_view, n_threads, i, data->rows);
		}
	}

	names.emplace_back("iteration");
	return_types.emplace_back(LogicalType::INTEGER);
	names.emplace_back("sample_id");
	return_types.emplace_back(data->sample_id_type);
	names.emplace_back("faith_pd");
	return_types.emplace_back(LogicalType::DOUBLE);

	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> UnifracFaithPdInitGlobal(ClientContext &, TableFunctionInitInput &input) {
	auto &data = input.bind_data->CastNoConst<UnifracFaithPdData>();
	auto gstate = make_uniq<UnifracFaithPdGlobalState>();
	gstate->rows = std::move(data.rows);
	gstate->sample_id_type = data.sample_id_type;
	return std::move(gstate);
}

void UnifracFaithPdExecute(ClientContext &, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<UnifracFaithPdGlobalState>();
	const idx_t total = gstate.rows.size();
	if (gstate.cursor >= total) {
		output.SetCardinality(0);
		return;
	}
	const idx_t remaining = total - gstate.cursor;
	const idx_t n = std::min<idx_t>(STANDARD_VECTOR_SIZE, remaining);

	auto iter_data = FlatVector::GetData<int32_t>(output.data[0]);
	auto &sample_id_vec = output.data[1];
	auto faith_pd_data = FlatVector::GetData<double>(output.data[2]);

	for (idx_t i = 0; i < n; ++i) {
		const auto &r = gstate.rows[gstate.cursor + i];
		iter_data[i] = r.iteration;
		// EmitIdCell mirrors the id type; its ""/"*"→NULL sentinel branch is
		// unreachable here (ReadFeatureTable drops NULL sample_ids).
		EmitIdCell(sample_id_vec, i, r.sample_id, gstate.sample_id_type);
		faith_pd_data[i] = r.faith_pd;
	}
	gstate.cursor += n;
	output.SetCardinality(n);
}

} // namespace

void RegisterUnifracFaithPD(ExtensionLoader &loader) {
	TableFunction fn("unifrac_faith_pd", {LogicalType::VARCHAR, LogicalType::VARCHAR}, UnifracFaithPdExecute,
	                 UnifracFaithPdBind, UnifracFaithPdInitGlobal);
	fn.named_parameters["subsample_depth"] = LogicalType::INTEGER;
	fn.named_parameters["subsample_with_replacement"] = LogicalType::BOOLEAN;
	fn.named_parameters["n_subsamples"] = LogicalType::INTEGER;
	fn.named_parameters["seed"] = LogicalType::INTEGER;
	fn.named_parameters["threads"] = LogicalType::INTEGER;
	fn.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
