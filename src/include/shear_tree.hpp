#pragma once

#include "read_newick.hpp"
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
// ignore_missing:=true. Output schema matches read_newick (without filepath).
class ShearTreeTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string tree_table_name;
		std::string tips_table_name;
		bool collapse;
		bool ignore_missing;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(std::string tree_table, std::string tips_table, bool collapse, bool ignore_missing);
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::vector<ReadNewickTableFunction::NodeRow> rows;
		size_t current_row_idx = 0;

		idx_t MaxThreads() const override {
			return 1;
		}
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
