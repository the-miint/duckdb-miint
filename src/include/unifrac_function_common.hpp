#pragma once

#include <array>
#include <string>
#include <vector>

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
std::vector<miint::unifrac::CooRow> ReadFeatureTable(ClientContext &context, const std::string &table_name,
                                                     const std::string &caller_name);

} // namespace duckdb::unifrac_internal
