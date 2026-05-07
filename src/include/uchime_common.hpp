#pragma once

#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"

#include <string>
#include <vector>

namespace miint {

// UCHIME scoring parameters with defaults from Edgar et al. 2011.
struct UchimeParams {
	double minh = 0.28;
	double xn = 8.0;
	double dn = 1.4;
	double mindiv = 0.8;
	int mindiffs = 3;
	double abskew = 2.0; // only used in de novo mode
	// vsearch's internal thread pool for chimera_detect_batch (uchime_ref).
	// 0 = vsearch auto-detect (= all physical cores). Callers normally pass
	// DuckDB's configured thread count so `SET threads=N` is honored. Ignored
	// in de novo mode, which is sequential by construction (one query per call,
	// then index_sequence) — denovo always runs with opt_threads=1.
	int threads = 0;
};

// Full UCHIME result for a single query (mirrors vsearch --uchimeout 18 columns).
struct UchimeResult {
	double score = 0.0;
	std::string query_label;
	std::string parent_a_label;
	std::string parent_b_label;
	std::string closest_parent_label;
	double id_query_model = 0.0;
	double id_query_a = 0.0;
	double id_query_b = 0.0;
	double id_a_b = 0.0;
	double id_query_top = 0.0;
	int left_yes = 0, left_no = 0, left_abstain = 0;
	int right_yes = 0, right_no = 0, right_abstain = 0;
	double divergence = 0.0;
	std::string flag = "N"; // Y, N, or ?
};

} // namespace miint

namespace duckdb {

// Output column names and types for the 18-column uchimeout format.
// Shared between uchime_ref and uchime_denovo.
std::vector<std::string> GetUchimeOutputNames();
std::vector<LogicalType> GetUchimeOutputTypes();

// Write a batch of UchimeResults into a DataChunk starting at column `start_col`.
// When start_col > 0 (per-sample callers), callers must populate columns [0, start_col)
// before/after this call; this function sets the chunk cardinality itself.
// Returns the number of rows written (min of count and remaining results).
idx_t OutputUchimeResults(DataChunk &output, const std::vector<miint::UchimeResult> &results, idx_t offset, idx_t count,
                          idx_t start_col = 0);

} // namespace duckdb
