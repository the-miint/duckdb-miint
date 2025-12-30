#pragma once

#include "duckdb/main/client_context.hpp"
#include <string>
#include <unordered_map>

namespace duckdb {

// Read a reference table or view and return map of reference name -> length.
// Used by read_alignments and COPY FORMAT SAM/BAM for headerless SAM support.
//
// Both tables and views are supported. The implementation uses a separate Connection
// to execute a SELECT statement, which DuckDB's binder resolves appropriately for
// either tables or views.
//
// Required columns (by position):
//   - Column 0: reference name (VARCHAR) - must not be NULL
//   - Column 1: reference length (BIGINT, INTEGER, UBIGINT, or UINTEGER) - must not be NULL
//
// Throws InvalidInputException if:
//   - Table/view doesn't exist
//   - Table/view has fewer than 2 columns
//   - Column types are incorrect
//   - Table/view is empty
//   - NULL value in name or length column
//   - Negative length value
//   - UBIGINT length exceeds INT64_MAX
std::unordered_map<std::string, uint64_t> ReadReferenceTable(ClientContext &context, const std::string &table_name);

} // namespace duckdb
