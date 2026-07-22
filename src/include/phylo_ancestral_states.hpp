#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <string>
#include <vector>

namespace duckdb {

// phylo_ancestral_states(tree_table, traits_table)
//
// Continuous ancestral state reconstruction under a Brownian-motion model over a
// tree (read_newick schema) and one or more numeric per-tip traits. Unlike PIC,
// multifurcating trees are supported (the reconstruction is arity-agnostic).
//
// traits_table is long form (name, trait, value): `name` is type-flexible
// (VARCHAR/UUID/BIGINT, matched to tip labels by canonical text), `trait` is a
// VARCHAR label, `value` is DOUBLE. Every trait must cover exactly the tree's tip
// set (no missing tips, no extra names, no NULLs). Output is long form, one row
// per (internal node, trait): node_index (the input tree's node_index), trait,
// estimate, variance, ci_low, ci_high.
class PhyloAncestralStatesTableFunction {
public:
	struct StateRow {
		int64_t node_index;
		std::string trait;
		double estimate;
		double variance;
		double ci_low;
		double ci_high;
	};

	struct Data : public TableFunctionData {
		std::string tree_table_name;
		std::string traits_table_name;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(std::string tree_table, std::string traits_table);
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::vector<StateRow> rows;
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
