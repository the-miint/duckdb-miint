#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// tree_lca(tree_table, pairs_table, ignore_missing := false)
//
// For each distinct pair (a, b) of tip names in pairs_table, the tips' lowest common
// ancestor (the caller's node_index) and the length of the path between them: the same
// double tree_nearest_tip reports for that pair, or NULL when a length on the path is
// NULL. Tip names may be VARCHAR, BIGINT or UUID, the tree's type in both columns, and
// come back in that type. See docs/phylogeny.md.
class TreeLcaTableFunction {
public:
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
