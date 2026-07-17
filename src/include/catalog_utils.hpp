#pragma once

#include "duckdb/common/types.hpp"
#include "duckdb/main/client_context.hpp"
#include <string>
#include <vector>

namespace duckdb {

struct TableOrViewColumns {
	vector<string> names;
	vector<LogicalType> types;
	bool is_physical_table; // true for tables (have rowid), false for views
};

// Look up a table or view by name and return its column names and types.
// entity_type is used in error messages (e.g., "Placement table", "Tree table").
// Throws BinderException if the table/view does not exist.
TableOrViewColumns GetTableOrViewColumns(ClientContext &context, const std::string &table_name,
                                         const std::string &entity_type = "Table");

// Case-insensitive check for whether `columns` contains a column named `col`.
bool HasColumn(const TableOrViewColumns &columns, const std::string &col);

} // namespace duckdb
