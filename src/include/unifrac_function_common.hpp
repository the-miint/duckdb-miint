#pragma once

#include <array>
#include <string>
#include <vector>

#include "duckdb/common/types.hpp"
#include "duckdb/main/client_context.hpp"
#include "unifrac_support_biom.hpp"

namespace duckdb::unifrac_internal {

// libssu's accepted variant names (the bare set, without `_fp32` suffix).
// Callers compose the libssu method string by appending `_fp32`.
//
// Defined inline so PCoA / PERMANOVA / Faith PD can share without forcing
// every translation unit to link a definition.
inline constexpr std::array<const char *, 5> kAcceptedVariants = {
    "unweighted", "weighted_normalized", "weighted_unnormalized", "unweighted_unnormalized", "generalized"};

inline bool IsValidVariant(const std::string &v) {
	for (const auto *name : kAcceptedVariants) {
		if (v == name) {
			return true;
		}
	}
	return false;
}

inline std::string AcceptedVariantList() {
	std::string out;
	for (size_t i = 0; i < kAcceptedVariants.size(); ++i) {
		if (i != 0) {
			out += ", ";
		}
		out += kAcceptedVariants[i];
	}
	return out;
}

// Read the user-named feature-table relation (must expose
// `(sample_id VARCHAR, feature_id VARCHAR, value DOUBLE)`, matching the
// columns produced by read_biom) into long-form COO rows.
//
// `caller_name` is prefixed to error messages so SQL-level binder
// errors identify the offending UDF (e.g., "unifrac_pcoa: ..." vs
// "unifrac_permanova: ...").
//
// Rows with NULL sample_id, feature_id, or value are silently dropped,
// as are zero/NaN values (UnifracSupportBiomView's sparse-storage
// invariant would drop them anyway).
//
// When `sample_id_type` / `feature_id_type` are non-null, the original SQL type of
// that column (before the internal `::VARCHAR` cast) is written through it, so
// callers can mirror BIGINT/UUID identifiers back onto their output columns.
// The ids themselves are always carried as canonical VARCHAR text internally.
// Most callers emit only sample ids and pass nullptr for the feature type;
// `rarefy_feature_table` needs both, because its output IS a feature table and has
// to round-trip both id columns.
std::vector<miint::unifrac::CooRow> ReadFeatureTable(ClientContext &context, const std::string &table_name,
                                                     const std::string &caller_name,
                                                     LogicalType *sample_id_type = nullptr,
                                                     LogicalType *feature_id_type = nullptr);

// A dense, symmetric, zero-diagonal fp32 distance matrix materialized from a
// condensed COO distance relation (`sample_a, sample_b, distance`). This is the
// metric-agnostic input shape shared by the `pcoa` and `permanova` table
// functions: any relation with those three columns can be read — the
// `unifrac_distances` output, a `beta_*` macro result, or a precomputed
// Bray-Curtis/Jaccard/Euclidean table — so ordination and the omnibus test are
// decoupled from UniFrac.
struct DenseDistanceMatrix {
	std::vector<float> matrix;           // n*n row-major, symmetric, zero diagonal (mat[i*n+j])
	std::vector<std::string> sample_ids; // distinct ids from both columns, lexicographically sorted
	uint32_t n_samples = 0;
	// Output type mirrored from sample_a's input column (BIGINT/UUID → same,
	// else VARCHAR); see ResolveSampleIdOutputType.
	LogicalType sample_id_type = LogicalType::VARCHAR;
};

// The distinct, sorted, non-null sample ids of a condensed distance relation,
// plus the resolved output id type. Materializes only the id dictionary (bounded
// by N), never the N×N matrix.
struct DistanceRelationIds {
	std::vector<std::string> sorted_ids;
	LogicalType sample_id_type = LogicalType::VARCHAR;
};

// Enumerate a condensed distance relation's id dictionary: probe the schema,
// reject mismatched sample_a/sample_b id types, and return the distinct non-null
// ids from both columns in lexicographic order. Shared by the progressive PCoA
// functions (which never build a dense matrix) and by ReadDistanceTable, where it
// is the first of two passes — knowing N up front is what lets the second pass
// stream cells straight into the matrix instead of parking every row on the heap.
DistanceRelationIds EnumerateDistanceIds(ClientContext &context, const std::string &table_name,
                                         const std::string &caller_name);

// Read the user-named condensed distance relation into a dense matrix. The
// relation must expose `(sample_a, sample_b, distance)`: sample_a/sample_b of
// any type castable to VARCHAR, distance castable to DOUBLE (mirrors
// ReadFeatureTable's probe-then-cast approach).
//
// Rows with a NULL sample_a/sample_b are skipped entirely (a NULL id has no
// identity). A row with a valid id pair but a NULL/NaN distance still registers
// both ids in the dictionary but contributes no distance ("not provided"), so a
// sample whose every distance is NULL/NaN does NOT silently vanish — it surfaces
// as a completeness error. Distinct ids are collected from BOTH columns and
// sorted lexicographically (the build_dictionary convention shared with
// UnifracSupportBiomView::FromCoo, so sample_a < sample_b holds for VARCHAR ids
// exactly as for unifrac_distances). The fill + validation (n>=2, completeness,
// conflicting duplicates, negative/non-finite, nonzero self-distance) is
// delegated to BuildDenseDistanceMatrix and re-thrown as InvalidInputException
// prefixed with `caller_name`.
//
// sample_a's original SQL type is captured (BIGINT/UUID mirrored, else VARCHAR)
// into DenseDistanceMatrix::sample_id_type. sample_a and sample_b are merged into
// one dictionary emitted under that type, so their resolved output types must
// match (else a BinderException) — sample_a's is mirrored.
DenseDistanceMatrix ReadDistanceTable(ClientContext &context, const std::string &table_name,
                                      const std::string &caller_name);

// Resolve the SQL type that a mirrored sample-id output column should carry.
// BIGINT and UUID are mirrored (so results join back to typed metadata without
// a cast — parity with align_minimap2); every other input type collapses to
// VARCHAR, matching what the ::VARCHAR-cast feature-table reader already
// accepts (so this never rejects a type ReadFeatureTable would have read).
//
// Deliberate divergence from align_minimap2: that function rejects any id
// column outside {VARCHAR, BIGINT, UUID} at bind via IsAllowedIdType. The
// unifrac readers stay intentionally permissive — a DOUBLE/DATE sample_id is
// accepted and simply emitted as VARCHAR, preserving the pre-existing behavior
// of the ::VARCHAR-cast reader rather than tightening it.
inline LogicalType ResolveSampleIdOutputType(const LogicalType &input_type) {
	if (input_type.id() == LogicalTypeId::BIGINT || input_type.id() == LogicalTypeId::UUID) {
		return input_type;
	}
	return LogicalType::VARCHAR;
}

// Resolve the user-supplied `threads` named parameter into a concrete count
// that the libssu / scikit-bio-binaries OpenMP regions will run with.
//
// Convention:
//   * `user_value == 0` (unset or explicit 0) → fall back to DuckDB's
//     TaskScheduler::NumberOfThreads(). We never fall back to OpenMP's
//     default (all cores) — duckdb-miint always pins thread count to
//     DuckDB's understanding of the system.
//   * `user_value > 0`                        → use as-is.
//   * `user_value < 0`                        → BinderException via caller_name.
//
// `caller_name` is the SQL function name used in the error message
// (e.g., "unifrac_pcoa").
//
// Returns a positive int suitable for OmpThreadPin / ComputeCallScope.
int ResolveThreadsParameter(ClientContext &context, int32_t user_value, const std::string &caller_name);

} // namespace duckdb::unifrac_internal
