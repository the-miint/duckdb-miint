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

// A built tree plus the caller's original node_index for each dense tree index. Because
// ReadTreeTable sorts by node_index and NewickTree::build preserves input order, dense
// index i corresponds to node_ids[i] — a stable join key for reporting results against the
// caller's node_index rather than the internal dense index.
struct TreeWithNodeIds {
	miint::NewickTree tree;
	std::vector<int64_t> node_ids;
};

// Read a tree table/view, capture its node_index mapping, and build the NewickTree. Shared
// by the comparative-methods table functions (independent contrasts, ancestral states,
// parsimony, ML). Throws InvalidInputException if the tree fails to build.
TreeWithNodeIds BuildTreeAndNodeIds(ClientContext &context, const std::string &table_name);

} // namespace duckdb
