#pragma once

#include "NewickTree.hpp"
#include "duckdb/main/client_context.hpp"
#include <string>
#include <vector>

namespace duckdb {

// Validate that a table or view has the required columns for tree data
// Required: node_index (integer), parent_index (integer, nullable)
// Optional: name (VARCHAR), branch_length (DOUBLE/FLOAT), edge_id (integer)
void ValidateTreeTableSchema(ClientContext &context, const std::string &table_name);

// Read tree node data from a table or view and return as NodeInput structs.
// Uses a separate Connection to avoid deadlocking the current context.
std::vector<miint::NodeInput> ReadTreeTable(ClientContext &context, const std::string &table_name);

} // namespace duckdb
