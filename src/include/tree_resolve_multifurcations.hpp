#pragma once

#include "tree_table_reader.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// tree_resolve_multifurcations(tree_table)
//
// Resolve every node with more than two children in a tree (in read_newick
// schema) into a series of bifurcations by inserting unnamed zero-length
// internal connector nodes (a deterministic left-comb). Nodes with two or fewer
// children are unchanged; single-child unifurcations are NOT collapsed (use
// shear_tree for that). Output schema matches read_newick (without filepath), so
// the result round-trips into phylo_independent_contrasts and read_newick, with
// `name` in the tree's id type and the tree's `label` column, when it has one, last.
class TreeResolveMultifurcationsTableFunction {
public:
	// Executes with TreeOutputExecute, from a TreeOutputState
	using Data = TreeOutputData;

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<Identifier> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
