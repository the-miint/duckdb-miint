#pragma once

#include "tree_table_reader.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// shear_tree(tree_table, tips_table, collapse := true, ignore_missing := false)
//
// Subset a tree (in read_newick schema) to the tips named in tips_table.
// - collapse=true  : remove single-descendant ancestors, summing branch
//                    lengths (standard phylogenetic shear); LCA of kept tips
//                    becomes the new root.
// - collapse=false : keep every internal node on the kept root-paths unchanged.
// Missing tips (names not found as tips in the tree) error unless
// ignore_missing:=true. Output schema matches read_newick (without filepath), with
// `name` in the tree's id type and the tree's `label` column, when it has one, last.
class ShearTreeTableFunction {
public:
	// Executes with TreeOutputExecute, from a TreeOutputState
	struct Data : public TreeOutputData {
		std::string tips_table_name;
		bool collapse;
		bool ignore_missing;

		Data(std::string tree_table, TreeNameColumns tree_columns, std::string tips_table, bool collapse,
		     bool ignore_missing);
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<Identifier> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
