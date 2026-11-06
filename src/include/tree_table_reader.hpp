#pragma once

#include "NewickTree.hpp"
#include "TreeIndex.hpp"
#include "duckdb/main/client_context.hpp"
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
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

// A tree table read as compact arrays, for functions that need the tree's shape and branch
// lengths but only a few of its names.
struct TreeTableArrays {
	miint::DenseTree tree;
	// (dense index, name) for every node whose name is in the requested set, tips and
	// internal nodes alike, in no particular order. A name on several nodes appears once per node.
	std::vector<std::pair<uint32_t, std::string>> named;
};

// Read a tree table or view once, streaming, keeping each node's name only when it is in
// `keep_names`. ReadTreeTable and BuildTreeAndNodeIds hold every row, name included, twice
// over: loading the Greengenes2 2024.09 tree (24.5M nodes) from Parquet that way was
// killed at an 18 GB memory cap, and this way peaked at 2.6 GB resident. Names are read as
// name::VARCHAR, as ReadTreeTable reads them; a missing name or branch_length column reads
// as NULL.
//
// `name_type` is the type Bind saw for the name column. A view can change type without a
// catalog change (one over a file that is rewritten, under a prepared statement), and names
// matched as text would then match the wrong rows or none, so a name column read with
// another type is refused. The check uses the type of the read itself.
//
// Throws std::invalid_argument when the name column's type has changed, or the rows are not
// a numbered tree (no rows, a NULL or repeated node_index, a parent_index that is no row's
// node_index), so the caller can say which function failed; a failed read throws
// InvalidInputException. The tree's shape (one root, no cycles) and branch lengths are left
// to TreeIndex::Build.
TreeTableArrays ReadTreeArrays(ClientContext &context, const std::string &table_name,
                               const std::unordered_set<std::string_view> &keep_names, const LogicalType &name_type);

} // namespace duckdb
