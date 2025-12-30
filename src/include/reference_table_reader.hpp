#pragma once

#include "duckdb/main/client_context.hpp"
#include <string>
#include <unordered_map>

namespace duckdb {

// Read a reference table and return map of reference name -> length.
// Used by read_alignments and COPY FORMAT SAM/BAM for headerless SAM support.
//
// IMPORTANT: Only tables are supported, not views. This is because view execution
// requires running queries during COPY operations, which causes deadlocks in DuckDB.
// View support may be added in the future if DuckDB provides a safe mechanism for
// executing queries during COPY finalize.
//
// Required columns (by position):
//   - Column 0: reference name (VARCHAR) - must not be NULL
//   - Column 1: reference length (BIGINT, INTEGER, UBIGINT, or UINTEGER) - must not be NULL
//
// Throws InvalidInputException if:
//   - Table doesn't exist (views will fail with this error)
//   - Table has fewer than 2 columns
//   - Column types are incorrect
//   - Table is empty
//   - NULL value in name or length column
//   - Negative length value
//   - UBIGINT length exceeds INT64_MAX
std::unordered_map<std::string, uint64_t> ReadReferenceTable(ClientContext &context, const std::string &table_name);

} // namespace duckdb
