#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <string>
#include <vector>

namespace duckdb {

// phylo_ancestral_parsimony(tree_table, traits_table [, cost_matrix_table])
//
// Discrete ancestral state reconstruction under Sankoff (1975) parsimony over a tree
// (read_newick schema) and one or more categorical per-tip traits. With the default
// unit cost (0 for no change, 1 for any change) this is Fitch (1971) parsimony;
// supply an optional cost_matrix for a general Sankoff cost. Multifurcations AND
// unifurcations are supported (the DP is arity-agnostic); branch lengths are ignored.
//
// traits_table is long form (name, trait, value): `name` is type-flexible
// (VARCHAR/UUID/BIGINT, matched to tip labels by canonical text), `trait` is a VARCHAR
// label, `value` is the categorical state cast to VARCHAR. Every trait must cover
// exactly the tree's tip set (no missing tips, no extra names, no NULLs).
//
// cost_matrix_table (optional) is long form (from_state, to_state, cost): from_state
// and to_state are the state labels (::VARCHAR), cost is a finite non-negative DOUBLE
// (the cost of an edge whose parent end is from_state and child end is to_state).
// When supplied it defines the state alphabet (the union of from_state/to_state); the
// alphabet may include states not observed at any tip (a Sankoff intermediate).
//
// Output is long form, one row per (internal node, trait, state): node_index (the
// input tree's node_index), trait, state, in_mpr (whether the state is in the node's
// most-parsimonious reconstruction set), min_cost (the minimum total tree cost with
// that node fixed to that state). Ties are first-class (multiple in_mpr states).
//
// Designed for SMALL categorical state alphabets (the Sankoff DP is O(n*k^2) time and
// O(n*k) memory in the number of states k). Pointing it at a high-cardinality column
// mistaken for a trait — e.g. a near-unique sample id — makes k approach the tip count
// and is correspondingly slow/memory-hungry.
class PhyloAncestralParsimonyTableFunction {
public:
	struct StateRow {
		int64_t node_index;
		std::string trait;
		std::string state;
		bool in_mpr;
		double min_cost;
	};

	struct Data : public TableFunctionData {
		std::string tree_table_name;
		std::string traits_table_name;
		bool has_cost_matrix = false;
		std::string cost_matrix_table_name;

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

	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
