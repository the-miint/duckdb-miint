#pragma once

#include "tree_table_reader.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class TreeResolvePlacementTableFunction {
public:
	// Executes with TreeOutputExecute, from a TreeOutputState
	struct Data : public TreeOutputData {
		std::string placements_table_name;

		Data(std::string tree_table, TreeNameColumns tree_columns, std::string placements_table);
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<Identifier> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
