#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// tree_ancestors(tree_table, nodes_table := NULL)
//
// One row per node and ancestor-or-self: the node itself (edges 0, distance 0), its
// parent, and so on to the root, with the edges and the path length between them. The
// distance is NULL when a length on the path is NULL. nodes_table, with a node_index
// column of any integer type, limits the rows to its nodes; without it, every node's.
// The rows stream, so a large tree's are never all held at once. See docs/phylogeny.md.
class TreeAncestorsTableFunction {
public:
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
