#pragma once

#include "duckdb.hpp"
#include "duckdb/function/table_function.hpp"

#include <atomic>

namespace duckdb {

// ── Per-sample table-function coordination ──────────────────────────────────
// Shared machinery for table functions that accept a `sample_id` named
// parameter and run their pipeline once per distinct sample value.
//
// Bind time: DiscoverSamples() validates the column, rejects empty names
// and output-column collisions, collects distinct non-NULL values (single
// DISTINCT … NULLS FIRST query + first-row NULL reject), captures the type.
// Exec time: InitPerSampleGlobal() caps thread count; ClaimNextSample() hands
// out sample indices via an atomic counter so each thread processes a
// different sample.
//
// Out of scope: per-sample data loading, result buffering, and the
// non-sample branch. Callers own those; they dispatch on `has_sample_id`
// themselves and use a per-thread Connection in LocalState.

struct PerSampleBindInfo {
	string sample_id_col;
	LogicalType sample_id_type;
	vector<Value> sample_values;
};

// reserved_lowercase_output_names: output column names (already lowercased)
// that the sample column must not match case-insensitively. Pass {} to skip.
// Accepts std::vector so callers with existing std::vector<std::string> fields
// (uchime's `names`) don't need to round-trip through duckdb::vector.
// fn_label: short function name prepended to error messages, e.g. "deblur".
void DiscoverSamples(Connection &conn, const string &source_relation, const string &sample_id_col,
                     const std::vector<std::string> &reserved_lowercase_output_names, const string &fn_label,
                     PerSampleBindInfo &out);

struct PerSampleGlobalState : public GlobalTableFunctionState {
	atomic<idx_t> next_sample_idx {0};
	idx_t max_threads = 1;

	idx_t MaxThreads() const override {
		return max_threads;
	}
};

// max_threads_hint: 0 → use DuckDB scheduler threads (clamped to num_samples).
// N>0 → cap to min(N, num_samples). Pass 1 when the per-sample backend is not
// thread-safe (MAFFT mutex, VsearchChimeraWrapper).
void InitPerSampleGlobal(ClientContext &context, PerSampleGlobalState &gstate, idx_t num_samples,
                         idx_t max_threads_hint = 0);

// Returns false once every sample has been claimed.
bool ClaimNextSample(PerSampleGlobalState &gstate, idx_t num_samples, idx_t &out_idx);

} // namespace duckdb
