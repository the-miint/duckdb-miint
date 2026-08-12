#include "unifrac_table_functions.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "id_column_utils.hpp"
#include "unifrac_function_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"

namespace duckdb {
namespace {

using unifrac_internal::ReadDistanceTable;

struct PickAnchorsData : public TableFunctionData {
	std::vector<std::string> anchors; // in selection order
	LogicalType sample_id_type = LogicalType::VARCHAR;
};

struct PickAnchorsGlobalState : public GlobalTableFunctionState {
	std::vector<std::string> anchors;
	size_t cursor = 0;
	LogicalType sample_id_type = LogicalType::VARCHAR;
	idx_t MaxThreads() const override {
		return 1;
	}
};

// Greedy farthest-point (k-center / Gonzalez) selection over a dense distance
// matrix, returning `k` sample indices in selection order.
//
// WHY THIS RULE. The progressive PCoA functions default to a seeded RANDOM anchor
// draw, and a random draw's quality varies a lot: measured M² against a full
// ordination ranged 0.087-0.262 across draws at a FIXED anchor count, because a
// draw can leave a region of the space with no nearby anchor. Every batch is
// aligned onto the anchors, so a region no anchor covers is a region whose
// procrustes fit is extrapolated. Farthest-point directly maximizes the minimum
// distance from any sample to its nearest anchor, which is exactly the quantity
// that failure mode is about.
//
// Not to be confused with true max-VOLUME selection (maximizing the determinant of
// the anchor Gram matrix). That needs an embedding to have a Gram matrix at all,
// which on this path would mean ordinating first — circular. Farthest-point is the
// standard cheap surrogate and has a 2-approximation guarantee for the k-center
// objective.
//
// Deterministic, with no seed: rank 0 is the most peripheral sample (largest total
// distance to all others) and every tie — including that one — breaks to the lowest
// sample index, i.e. the lexicographically smallest id, since ReadDistanceTable
// sorts its dictionary. So an anchor set is a reproducible property of the data and
// needs no seed recorded alongside it.
//
// Cost: O(N²) for the row sums (the matrix is already materialized) then O(N·k) for
// the selection, holding one double per sample.
std::vector<uint32_t> GreedyFarthestPoint(const float *matrix, uint32_t n, uint32_t k) {
	// `min_dist[i]` = distance from i to the nearest already-chosen anchor.
	// Chosen samples are set to -1 so they can never be picked twice.
	std::vector<double> min_dist(n, std::numeric_limits<double>::infinity());
	std::vector<uint32_t> chosen;
	chosen.reserve(k);

	uint32_t first = 0;
	{
		double best = -1.0;
		for (uint32_t i = 0; i < n; ++i) {
			double sum = 0.0;
			const float *row = matrix + static_cast<size_t>(i) * n;
			for (uint32_t j = 0; j < n; ++j) {
				sum += row[j];
			}
			if (sum > best) { // strict >: ties keep the lower index
				best = sum;
				first = i;
			}
		}
	}

	const auto take = [&](uint32_t pick) {
		chosen.push_back(pick);
		min_dist[pick] = -1.0;
		const float *row = matrix + static_cast<size_t>(pick) * n;
		for (uint32_t i = 0; i < n; ++i) {
			if (min_dist[i] >= 0.0 && row[i] < min_dist[i]) {
				min_dist[i] = row[i];
			}
		}
	};
	take(first);

	while (chosen.size() < k) {
		uint32_t best_i = 0;
		double best_d = -1.0;
		for (uint32_t i = 0; i < n; ++i) {
			if (min_dist[i] > best_d) { // strict >: ties keep the lower index
				best_d = min_dist[i];
				best_i = i;
			}
		}
		// best_d < 0 would mean every sample is already chosen, which the
		// n_anchors <= n_samples check upstream makes unreachable.
		take(best_i);
	}
	return chosen;
}

unique_ptr<FunctionData> PickAnchorsBind(ClientContext &context, TableFunctionBindInput &input,
                                         vector<LogicalType> &return_types, vector<string> &names) {
	const std::string table_name = input.inputs[0].GetValue<string>();
	if (table_name.empty()) {
		throw BinderException("pick_anchors: distance-table name must not be empty");
	}

	bool has_n = false;
	int32_t n_anchors = 0;
	for (const auto &kv : input.named_parameters) {
		const auto key = StringUtil::Lower(kv.first);
		if (key == "n_anchors") {
			n_anchors = kv.second.GetValue<int32_t>();
			has_n = true;
		}
	}
	if (!has_n) {
		throw BinderException("pick_anchors: n_anchors is required (there is no default anchor count); "
		                      "pass n_anchors := N");
	}
	if (n_anchors < 1) {
		throw BinderException("pick_anchors: n_anchors must be >= 1 (got %d)", n_anchors);
	}

	// ReadDistanceTable carries the fail-loud N² size guard, so a distance table too
	// large for the dense matrix errors with a clear message rather than an OOM kill.
	auto dm = ReadDistanceTable(context, table_name, "pick_anchors");
	if (static_cast<uint32_t>(n_anchors) > dm.n_samples) {
		throw BinderException("pick_anchors: n_anchors (%d) exceeds the %u distinct sample(s) in the distance-table",
		                      n_anchors, dm.n_samples);
	}

	auto data = make_uniq<PickAnchorsData>();
	data->sample_id_type = dm.sample_id_type;
	const auto picked = GreedyFarthestPoint(dm.matrix.data(), dm.n_samples, static_cast<uint32_t>(n_anchors));
	data->anchors.reserve(picked.size());
	for (const auto idx : picked) {
		data->anchors.push_back(dm.sample_ids[idx]);
	}

	names.emplace_back("anchor_rank");
	return_types.emplace_back(LogicalType::INTEGER);
	names.emplace_back("sample_id");
	return_types.emplace_back(data->sample_id_type);
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> PickAnchorsInitGlobal(ClientContext &, TableFunctionInitInput &input) {
	auto &data = input.bind_data->CastNoConst<PickAnchorsData>();
	auto gstate = make_uniq<PickAnchorsGlobalState>();
	gstate->anchors = std::move(data.anchors);
	gstate->sample_id_type = data.sample_id_type;
	return std::move(gstate);
}

void PickAnchorsExecute(ClientContext &, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<PickAnchorsGlobalState>();
	const idx_t total = gstate.anchors.size();
	if (gstate.cursor >= total) {
		output.SetCardinality(0);
		return;
	}
	const idx_t n = std::min<idx_t>(STANDARD_VECTOR_SIZE, total - gstate.cursor);

	auto rank_data = FlatVector::GetData<int32_t>(output.data[0]);
	auto &sample_id_vec = output.data[1];
	for (idx_t i = 0; i < n; ++i) {
		rank_data[i] = static_cast<int32_t>(gstate.cursor + i);
		EmitIdCell(sample_id_vec, i, gstate.anchors[gstate.cursor + i], gstate.sample_id_type);
	}
	gstate.cursor += n;
	output.SetCardinality(n);
}

} // namespace

void RegisterPickAnchors(ExtensionLoader &loader) {
	TableFunction fn("pick_anchors", {LogicalType::VARCHAR}, PickAnchorsExecute, PickAnchorsBind,
	                 PickAnchorsInitGlobal);
	fn.named_parameters["n_anchors"] = LogicalType::INTEGER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
