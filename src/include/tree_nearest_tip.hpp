#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// tree_nearest_tip(tree_table, targets_table, queries_table,
//                  scope := 'path', exclude_self := true, ignore_missing := false)
//
// For each query tip, the nearest target tip by path length (or, with scope := 'clade',
// the nearest within the smallest clade around the query that holds a target), how many
// targets tie at that distance, and the query-target LCA. With a `grp` column in both the
// targets and the queries, a query competes only among the targets of its own group.
// Tip names may be VARCHAR, BIGINT or UUID, one type in all three tables, and come back
// in that type. See docs/phylogeny.md.
class TreeNearestTipTableFunction {
public:
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
