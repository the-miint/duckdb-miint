#pragma once

#include "NewickTree.hpp"
#include "duckdb/main/client_context.hpp"
#include <string>
#include <vector>

namespace duckdb {

// Read a placement table and return vector of Placement structs.
//
// IMPORTANT: Only tables are supported, not views. This is because view execution
// requires running queries during COPY operations, which causes deadlocks in DuckDB.
// View support may be added in the future if DuckDB provides a safe mechanism for
// executing queries during COPY finalize.
//
// Required columns (by name, case-insensitive):
//   - fragment_id (VARCHAR) - NULL values are silently skipped
//   - edge_id (BIGINT, INTEGER, UBIGINT, or UINTEGER) - must not be NULL
//   - like_weight_ratio (DOUBLE or FLOAT) - must not be NULL
//   - distal_length (DOUBLE or FLOAT) - must not be NULL
//   - pendant_length (DOUBLE or FLOAT) - must not be NULL
//
// Throws InvalidInputException if:
//   - Table doesn't exist
//   - Required columns are missing
//   - Table is empty (after skipping NULL fragment_ids)
//   - NULL value in edge_id, like_weight_ratio, distal_length, or pendant_length
//   - UBIGINT edge_id exceeds INT64_MAX
std::vector<miint::Placement> ReadPlacementTable(ClientContext &context, const std::string &table_name);

// Validate that a placement table has the required columns with correct types.
// Called at bind time to provide early error detection.
// Only validates tables (not views) - see note above.
// Throws BinderException if validation fails.
void ValidatePlacementTableSchema(ClientContext &context, const std::string &table_name);

} // namespace duckdb
