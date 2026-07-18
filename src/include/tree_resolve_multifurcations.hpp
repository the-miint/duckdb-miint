#pragma once

#include "read_newick.hpp"
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
// the result round-trips into phylo_independent_contrasts and read_newick.
class TreeResolveMultifurcationsTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string tree_table_name;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		explicit Data(std::string tree_table);
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
