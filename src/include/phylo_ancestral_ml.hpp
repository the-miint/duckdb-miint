#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <string>
#include <vector>

namespace duckdb {

// phylo_ancestral_ml(tree_table, traits_table, model := 'ER', rate := NULL)
//
// Discrete ancestral state reconstruction by maximum likelihood under the Mk model
// over a tree (read_newick schema) and one or more categorical per-tip traits.
// Felsenstein pruning with per-node rescaling; a uniform root prior; marginal
// posteriors via a two-pass (Yang, Kumar & Nei 1995). Multifurcations and
// unifurcations are supported; branch lengths are used (P(t) depends on t).
//
// Three models are supported: `model := 'ER'` (equal rates, the default), `model := 'SYM'`
// (symmetric rates q_ij = q_ji, k(k-1)/2 free off-diagonal rates), and `model := 'ARD'`
// (all rates different, k(k-1) free off-diagonal rates; non-reversible, so the
// reconstruction depends on root placement). For ER the single rate is fitted by maximum
// likelihood unless `rate := <positive double>` is supplied, in which case that rate is
// used directly. For SYM and ARD the whole rate matrix is always fitted by maximum
// likelihood, so `rate` does not apply (supplying it is an error) and the reported `rate`
// column is NULL. SYM (k > 8) and ARD (k > 6) reject high state counts (the simplex fit is
// unreliable); the ARD fit is best-effort and not guaranteed to be the global optimum (see
// the .cpp / docs/phylogeny.md).
//
// traits_table is long form (name, trait, value): `name` is type-flexible
// (VARCHAR/UUID/BIGINT, matched to tip labels by canonical text), `trait` is a VARCHAR
// label, `value` is the categorical state cast to VARCHAR. The per-trait state
// alphabet is the sorted distinct observed states. Every trait must cover exactly the
// tree's tip set (no missing tips, no extra names, no NULLs); every non-root edge must
// have a finite non-negative branch length and every tip edge a strictly positive one.
//
// Output is long form, one row per (internal node, trait, state): node_index (the
// input tree's node_index), trait, state, probability (marginal posterior, summing to
// 1 per node), rate (fitted or fixed, repeated per trait), log_likelihood (at that
// rate, repeated per trait).
class PhyloAncestralMLTableFunction {
public:
	struct StateRow {
		int64_t node_index;
		std::string trait;
		std::string state;
		double probability;
		double rate;
		bool rate_null; // true for SYM (no single scalar rate) -> emit SQL NULL
		double log_likelihood;
	};

	struct Data : public TableFunctionData {
		std::string tree_table_name;
		std::string traits_table_name;
		std::string model = "ER";
		bool has_fixed_rate = false;
		double fixed_rate = 0.0;

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
