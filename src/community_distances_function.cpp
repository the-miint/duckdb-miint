#include "community_distances.hpp"

#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "unifrac_function_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace duckdb {

namespace {

using unifrac_internal::ReadFeatureTable;
using unifrac_internal::ResolveSampleIdOutputType;
using unifrac_internal::ResolveThreadsParameter;

struct CommunityDistBindData : public TableFunctionData {
	std::string table_name;
	std::string metric;
	// Output type for sample_a/sample_b, mirrored from the feature table's
	// sample_id column (BIGINT/UUID preserved, else VARCHAR).
	LogicalType sample_id_type = LogicalType::VARCHAR;
	// Threads for the O(n^2 * f) pair loop (resolved from the `threads` named
	// param at bind: 0 = follow DuckDB's thread count; always >= 1 here).
	int n_threads = 1;
};

// Holds the fully-computed condensed distance list. Single-threaded emission
// (MaxThreads == 1): the cursor is a plain counter and the whole result is
// materialized up front in InitGlobal (the sample counts here are small — the
// Kuczynski designs top out at a few thousand samples).
struct CommunityDistGlobalState : public GlobalTableFunctionState {
	std::vector<std::string> sample_ids; // distinct, lexicographically sorted
	std::vector<uint32_t> pair_a;        // index into sample_ids (sample_a)
	std::vector<uint32_t> pair_b;        // index into sample_ids (sample_b)
	std::vector<double> dist;            // condensed distances, aligned to pair_a/pair_b
	LogicalType sample_id_type = LogicalType::VARCHAR;
	idx_t cursor = 0;
	idx_t MaxThreads() const override {
		return 1;
	}
};

unique_ptr<FunctionData> CommunityDistBind(ClientContext &context, TableFunctionBindInput &input,
                                           vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<CommunityDistBindData>();
	data->table_name = input.inputs[0].GetValue<string>();
	data->metric = StringUtil::Lower(input.inputs[1].GetValue<string>());
	if (!miint::IsValidCommunityMetric(data->metric)) {
		throw BinderException("community_distances: unknown metric '%s' (must be one of %s)", data->metric,
		                      miint::CommunityMetricList());
	}

	int32_t threads = 0; // 0 = follow DuckDB's TaskScheduler::NumberOfThreads()
	for (const auto &kv : input.named_parameters) {
		if (StringUtil::Lower(kv.first) == "threads") {
			threads = kv.second.GetValue<int32_t>();
		}
	}
	data->n_threads = ResolveThreadsParameter(context, threads, "community_distances");

	// Mirror the sample_id column's type onto the output ids so BIGINT/UUID
	// results join back to typed metadata without a cast (parity with
	// unifrac_distances). GetTableOrViewColumns throws a clear BinderException if
	// the relation is absent; the full column contract is enforced by
	// ReadFeatureTable in InitGlobal.
	auto cols = GetTableOrViewColumns(context, data->table_name, "feature-table");
	for (idx_t i = 0; i < cols.names.size(); ++i) {
		if (StringUtil::Lower(cols.names[i]) == "sample_id") {
			data->sample_id_type = ResolveSampleIdOutputType(cols.types[i]);
			break;
		}
	}

	names = {"sample_a", "sample_b", "distance"};
	return_types = {data->sample_id_type, data->sample_id_type, LogicalType::DOUBLE};
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> CommunityDistInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<CommunityDistBindData>();
	auto gstate = make_uniq<CommunityDistGlobalState>();
	gstate->sample_id_type = data.sample_id_type;

	// ReadFeatureTable drops NULL/zero/NaN cells (the sparse feature-table
	// contract shared with unifrac_distances). A sample whose every value is 0
	// therefore contributes no rows and is absent from the dictionary below — it
	// is silently dropped, exactly as unifrac_distances drops it, because the
	// sparse COO cannot distinguish an all-zero sample from an absent one. (The
	// pure core's empty-community handling is defensive for dense callers and is
	// not reachable on this path.)
	auto rows = ReadFeatureTable(context, data.table_name, "community_distances");

	// Distinct sample ids (lexicographically sorted — the unifrac dictionary
	// convention, so sample_a < sample_b matches unifrac_distances output) and
	// distinct feature ids (first-seen order; feature order is irrelevant since
	// every metric is a symmetric sum over features).
	std::vector<std::string> sample_ids;
	std::vector<std::string> feature_ids;
	std::unordered_set<std::string> s_seen;
	std::unordered_set<std::string> f_seen;
	for (const auto &r : rows) {
		if (s_seen.insert(r.sample_id).second) {
			sample_ids.push_back(r.sample_id);
		}
		if (f_seen.insert(r.feature_id).second) {
			feature_ids.push_back(r.feature_id);
		}
	}
	std::sort(sample_ids.begin(), sample_ids.end());

	std::unordered_map<std::string, uint32_t> s_index;
	std::unordered_map<std::string, uint32_t> f_index;
	s_index.reserve(sample_ids.size());
	f_index.reserve(feature_ids.size());
	for (uint32_t i = 0; i < sample_ids.size(); ++i) {
		s_index.emplace(sample_ids[i], i);
	}
	for (uint32_t k = 0; k < feature_ids.size(); ++k) {
		f_index.emplace(feature_ids[k], k);
	}

	const uint32_t n = static_cast<uint32_t>(sample_ids.size());
	const uint32_t f = static_cast<uint32_t>(feature_ids.size());

	// Dense sample x feature abundance matrix; duplicate (sample, feature) rows
	// are summed (defensive — a well-formed feature table has one row per cell).
	std::vector<double> matrix(static_cast<size_t>(n) * static_cast<size_t>(f), 0.0);
	for (const auto &r : rows) {
		const uint32_t si = s_index.at(r.sample_id);
		const uint32_t fi = f_index.at(r.feature_id);
		matrix[static_cast<size_t>(si) * f + fi] += r.count;
	}

	// The pure core throws std::invalid_argument on data-dependent degeneracies
	// (e.g. fewer than 2 distinct samples). Surface as InvalidInputException with
	// the core's already-prefixed message.
	std::vector<double> condensed;
	try {
		condensed =
		    miint::CommunityDistancesCondensed(matrix, n, f, data.metric, static_cast<unsigned>(data.n_threads));
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s", e.what());
	}

	// Enumerate the upper-triangle pair indices in the SAME row-major (i<j) order
	// CommunityDistancesCondensed returns, so pair_a/pair_b align with `condensed`.
	gstate->dist = std::move(condensed);
	gstate->pair_a.reserve(gstate->dist.size());
	gstate->pair_b.reserve(gstate->dist.size());
	for (uint32_t i = 0; i + 1 < n; ++i) {
		for (uint32_t j = i + 1; j < n; ++j) {
			gstate->pair_a.push_back(i);
			gstate->pair_b.push_back(j);
		}
	}
	gstate->sample_ids = std::move(sample_ids);
	return std::move(gstate);
}

void CommunityDistExecute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &g = data_p.global_state->Cast<CommunityDistGlobalState>();
	const idx_t total = g.dist.size();
	const idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - g.cursor);

	auto &va = output.data[0];
	auto &vb = output.data[1];
	auto dd = FlatVector::GetData<double>(output.data[2]);

	for (idx_t r = 0; r < count; ++r) {
		const idx_t k = g.cursor + r;
		EmitIdCell(va, r, g.sample_ids[g.pair_a[k]], g.sample_id_type);
		EmitIdCell(vb, r, g.sample_ids[g.pair_b[k]], g.sample_id_type);
		dd[r] = g.dist[k];
	}
	g.cursor += count;
	output.SetCardinality(count);
}

} // namespace

void RegisterCommunityDistances(ExtensionLoader &loader) {
	TableFunction fn("community_distances", {LogicalType::VARCHAR, LogicalType::VARCHAR}, CommunityDistExecute,
	                 CommunityDistBind, CommunityDistInitGlobal);
	// Threads for the internal pair-loop parallelism (0 = follow DuckDB). The
	// distances are computed up front in InitGlobal; row emission stays
	// single-threaded (MaxThreads() == 1), so this only scales the compute.
	fn.named_parameters["threads"] = LogicalType::INTEGER;
	fn.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
