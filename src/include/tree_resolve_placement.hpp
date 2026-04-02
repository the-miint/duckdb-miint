#pragma once

#include "read_newick.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class TreeResolvePlacementTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string tree_table_name;
		std::string placements_table_name;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(std::string tree_table, std::string placements_table);
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
