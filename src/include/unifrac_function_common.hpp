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
// When `sample_id_type` is non-null, the original SQL type of the `sample_id`
// column (before the internal `::VARCHAR` cast) is written through it, so
// callers can mirror BIGINT/UUID identifiers back onto their output columns.
// The ids themselves are always carried as canonical VARCHAR text internally.
std::vector<miint::unifrac::CooRow> ReadFeatureTable(ClientContext &context, const std::string &table_name,
                                                     const std::string &caller_name,
                                                     LogicalType *sample_id_type = nullptr);

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
// Returns a positive int suitable for OmpThreadScope.
int ResolveThreadsParameter(ClientContext &context, int32_t user_value, const std::string &caller_name);

} // namespace duckdb::unifrac_internal
