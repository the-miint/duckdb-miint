#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// tree_node_depths(tree_table)
//
// Every node's depth (edges from the root) and distance from the root (the sum of the
// branch lengths above it). The distance is NULL when a length above the node is NULL;
// negative lengths are summed like any other. node_index is the caller's. See
// docs/phylogeny.md.
class TreeNodeDepthsTableFunction {
public:
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
