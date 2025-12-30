#pragma once

#include "NewickTree.hpp"
#include "duckdb/main/client_context.hpp"
#include <string>
#include <vector>

namespace duckdb {

// Read a placement table or view and return vector of Placement structs.
//
// Both tables and views are supported. The implementation uses ClientContext::Query()
// to execute a SELECT statement, which DuckDB's binder resolves appropriately for
// either tables or views.
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

// Validate that a placement table or view has the required columns with correct types.
// Called at bind time to provide early error detection.
// Supports both tables and views.
// Throws BinderException if validation fails.
void ValidatePlacementTableSchema(ClientContext &context, const std::string &table_name);

} // namespace duckdb
