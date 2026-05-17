#pragma once
//
// Pure codec for the BIGINT side of id-column round-tripping.
//
// duckdb-miint stores identifier columns as `std::vector<std::string>` in its
// in-memory carriers (SequenceRecordBatch::read_ids, SAMRecordBatch::read_ids,
// references, mate_references), regardless of whether the user's SQL-level
// column is VARCHAR or BIGINT. When the column is BIGINT, the values are
// stringified on ingress and parsed back on egress. These two functions
// encode that contract.
//
// Kept header-light (no DuckDB types) so the unit-test target can link it
// directly. The DuckDB-aware wrappers live in id_column_utils.{hpp,cpp}.

#include <cstdint>
#include <optional>
#include <string>

namespace miint {

// Parse a stringified id back to int64. Locale-independent (uses std::from_chars).
// Returns std::nullopt for empty string and for the SAM-unmapped sentinel "*"
// — both encode NULL in the output column.
// Throws std::invalid_argument on any other unparseable input (including the
// SAM same-as-primary sentinel "=", which is the caller's responsibility to
// resolve to a real reference value before invoking the codec).
// Throws std::out_of_range on overflow.
// Note: zero-prefixed and "-0" inputs are accepted but the round-trip through
// FormatIdFromInt64 will canonicalise them ("01" -> 1 -> "1"; "-0" -> 0 -> "0").
std::optional<int64_t> ParseIdAsInt64(const std::string &s);

// Format an int64 id as a decimal string for embedding in the carrier.
std::string FormatIdFromInt64(int64_t v);

} // namespace miint
