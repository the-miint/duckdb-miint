#include "phylo_ancestral_parsimony.hpp"
#include "tree_table_reader.hpp"
#include "phylo_traits_reader.hpp"
#include "catalog_utils.hpp"
#include "NewickTree.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/parser/keyword_helper.hpp"
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace duckdb {

namespace {

// A resolved substitution cost matrix: a sorted state alphabet plus a k*k row-major
// cost matrix (cost[from*k + to]).
struct CostMatrix {
	std::vector<std::string> alphabet;               // sorted distinct states
	std::unordered_map<std::string, uint32_t> index; // state -> 0..k-1
	std::vector<double> cost;                        // k*k, cost[from*k + to]
};

// Read a long-form (from_state, to_state, cost) cost matrix. The alphabet is the
// sorted union of from_state/to_state. Every off-diagonal ordered pair must be
// present (fail loud on a gap — there is no sensible default for a missing
// transition cost); the diagonal defaults to 0 and must be 0 if given. Costs must be
// finite and non-negative. Uses a separate Connection (docs/internals/
// reading-tables-views.md).
CostMatrix ReadCostMatrix(ClientContext &context, const std::string &table_name) {
	auto conn = MakeReadOnlyHelperConnection(context);
	std::string q = "SELECT from_state::VARCHAR, to_state::VARCHAR, cost::DOUBLE FROM " +
	                KeywordHelper::WriteOptionallyQuoted(table_name);
	auto res = conn.Query(q);
	if (res->HasError()) {
		throw InvalidInputException("Failed to read from cost matrix table '%s': %s", table_name, res->GetError());
	}

	std::map<std::pair<std::string, std::string>, double> entries;
	std::set<std::string> states;
	auto &mat = res->Cast<MaterializedQueryResult>();
	while (true) {
		auto chunk = mat.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}
		UnifiedVectorFormat fd, td, cd;
		chunk->data[0].ToUnifiedFormat(chunk->size(), fd);
		chunk->data[1].ToUnifiedFormat(chunk->size(), td);
		chunk->data[2].ToUnifiedFormat(chunk->size(), cd);
		auto fs = UnifiedVectorFormat::GetData<string_t>(fd);
		auto ts = UnifiedVectorFormat::GetData<string_t>(td);
		auto cs = UnifiedVectorFormat::GetData<double>(cd);
		for (idx_t i = 0; i < chunk->size(); i++) {
			auto fi = fd.sel->get_index(i), ti = td.sel->get_index(i), ci = cd.sel->get_index(i);
			if (!fd.validity.RowIsValid(fi) || !td.validity.RowIsValid(ti)) {
				throw InvalidInputException("NULL from_state/to_state in cost matrix table '%s'", table_name);
			}
			if (!cd.validity.RowIsValid(ci)) {
				throw InvalidInputException("NULL cost in cost matrix table '%s'", table_name);
			}
			std::string f = fs[fi].GetString(), to = ts[ti].GetString();
			double c = cs[ci];
			if (!std::isfinite(c) || c < 0.0) {
				throw InvalidInputException("cost matrix entry (%s -> %s) must be finite and non-negative in '%s'", f,
				                            to, table_name);
			}
			states.insert(f);
			states.insert(to);
			if (!entries.emplace(std::make_pair(f, to), c).second) {
				throw InvalidInputException("duplicate cost matrix entry (%s -> %s) in '%s'", f, to, table_name);
			}
		}
	}
	if (states.empty()) {
		throw InvalidInputException("cost matrix table '%s' is empty", table_name);
	}

	CostMatrix cm;
	cm.alphabet.assign(states.begin(), states.end()); // std::set iterates in sorted order
	for (uint32_t i = 0; i < cm.alphabet.size(); i++) {
		cm.index[cm.alphabet[i]] = i;
	}
	const uint32_t k = static_cast<uint32_t>(cm.alphabet.size());
	cm.cost.assign(static_cast<size_t>(k) * k, 0.0);
	for (uint32_t i = 0; i < k; i++) {
		for (uint32_t j = 0; j < k; j++) {
			auto it = entries.find({cm.alphabet[i], cm.alphabet[j]});
			if (i == j) {
				if (it != entries.end() && it->second != 0.0) {
					throw InvalidInputException("cost matrix diagonal (%s -> %s) must be 0 in '%s'", cm.alphabet[i],
					                            cm.alphabet[j], table_name);
				}
			} else if (it == entries.end()) {
				throw InvalidInputException("cost matrix is missing entry (%s -> %s) in '%s'; every off-diagonal "
				                            "ordered pair over the state alphabet must be specified",
				                            cm.alphabet[i], cm.alphabet[j], table_name);
			} else {
				cm.cost[static_cast<size_t>(i) * k + j] = it->second;
			}
		}
	}
	return cm;
}

} // namespace

PhyloAncestralParsimonyTableFunction::Data::Data(std::string tree_table, std::string traits_table)
    : tree_table_name(std::move(tree_table)), traits_table_name(std::move(traits_table)) {
	names = {"node_index", "trait", "state", "in_mpr", "min_cost"};
	types = {LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::BOOLEAN,
	         LogicalType::DOUBLE};
}

unique_ptr<FunctionData> PhyloAncestralParsimonyTableFunction::Bind(ClientContext &context,
                                                                    TableFunctionBindInput &input,
                                                                    vector<LogicalType> &return_types,
                                                                    vector<std::string> &names) {
	auto tree_table_name = input.inputs[0].ToString();
	auto traits_table_name = input.inputs[1].ToString();

	// Validate schemas at bind time for early error detection.
	ValidateTreeTableSchema(context, tree_table_name);
	ValidateTraitsTableSchema(context, traits_table_name);

	auto data = duckdb::make_uniq<Data>(std::move(tree_table_name), std::move(traits_table_name));

	if (input.inputs.size() >= 3) {
		data->has_cost_matrix = true;
		data->cost_matrix_table_name = input.inputs[2].ToString();
		auto info = GetTableOrViewColumns(context, data->cost_matrix_table_name, "Cost matrix table");
		for (const char *required : {"from_state", "to_state", "cost"}) {
			if (!HasColumn(info, required)) {
				throw BinderException("Cost matrix table '%s' missing required column '%s'",
				                      data->cost_matrix_table_name, required);
			}
		}
	}

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> PhyloAncestralParsimonyTableFunction::InitGlobal(ClientContext &context,
                                                                                      TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = duckdb::make_uniq<GlobalState>();

	auto tree_data = BuildTreeAndNodeIds(context, bind_data.tree_table_name);
	auto &tree = tree_data.tree;
	auto &tree_index_to_node_id = tree_data.node_ids;

	auto traits = ReadDiscreteTraits(context, bind_data.traits_table_name);

	CostMatrix cm;
	if (bind_data.has_cost_matrix) {
		cm = ReadCostMatrix(context, bind_data.cost_matrix_table_name);
	}

	// Build per-trait (tip-state indices, k, cost, alphabet). The unit-cost default
	// builds a per-trait alphabet from that trait's observed states; a supplied cost
	// matrix defines a shared alphabet (possibly wider than the observed states).
	std::vector<std::string> trait_names;
	std::vector<std::unordered_map<std::string, uint32_t>> tip_states_list;
	std::vector<uint32_t> k_list;
	std::vector<std::vector<double>> cost_list;
	std::vector<std::vector<std::string>> alphabet_list; // to map state index -> label on output

	for (auto &entry : traits) {
		const std::string &trait = entry.first;
		const auto &tipmap = entry.second;

		std::vector<std::string> alphabet;
		std::unordered_map<std::string, uint32_t> index;
		std::vector<double> cost;
		if (bind_data.has_cost_matrix) {
			alphabet = cm.alphabet;
			index = cm.index;
			cost = cm.cost;
		} else {
			std::set<std::string> distinct;
			for (const auto &[tip, st] : tipmap) {
				distinct.insert(st);
			}
			alphabet.assign(distinct.begin(), distinct.end()); // std::set -> sorted
			for (uint32_t i = 0; i < alphabet.size(); i++) {
				index[alphabet[i]] = i;
			}
			const uint32_t k = static_cast<uint32_t>(alphabet.size());
			cost.assign(static_cast<size_t>(k) * k, 1.0);
			for (uint32_t i = 0; i < k; i++) {
				cost[static_cast<size_t>(i) * k + i] = 0.0;
			}
		}
		const uint32_t k = static_cast<uint32_t>(alphabet.size());

		std::unordered_map<std::string, uint32_t> tip_states;
		for (const auto &[tip, st] : tipmap) {
			auto it = index.find(st);
			if (it == index.end()) {
				// Only reachable in the cost-matrix case (the unit-cost alphabet is the
				// observed states themselves).
				throw InvalidInputException(
				    "phylo_ancestral_parsimony: trait '%s' tip '%s' has state '%s' not present in the cost matrix",
				    trait, tip, st);
			}
			tip_states[tip] = it->second;
		}

		trait_names.push_back(trait);
		tip_states_list.push_back(std::move(tip_states));
		k_list.push_back(k);
		cost_list.push_back(std::move(cost));
		alphabet_list.push_back(std::move(alphabet));
	}

	std::vector<std::vector<miint::AncestralStateParsimony>> per_trait;
	try {
		per_trait = tree.ancestral_parsimony(tip_states_list, k_list, cost_list, &tree_index_to_node_id);
	} catch (const std::exception &e) {
		// Cold path: re-run per trait to attribute a completeness/state failure to the
		// specific trait (structural errors are trait-independent and surface first).
		for (size_t t = 0; t < tip_states_list.size(); t++) {
			try {
				tree.ancestral_parsimony(tip_states_list[t], k_list[t], cost_list[t], &tree_index_to_node_id);
			} catch (const std::exception &inner) {
				throw InvalidInputException("phylo_ancestral_parsimony failed for trait '%s': %s", trait_names[t],
				                            inner.what());
			}
		}
		throw InvalidInputException("phylo_ancestral_parsimony failed: %s", e.what());
	}

	for (size_t t = 0; t < trait_names.size(); t++) {
		for (const auto &st : per_trait[t]) {
			gstate->rows.push_back(StateRow {tree_index_to_node_id[st.node], trait_names[t], alphabet_list[t][st.state],
			                                 st.in_mpr, st.min_cost});
		}
	}

	return gstate;
}

void PhyloAncestralParsimonyTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p,
                                                   DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	if (gstate.current_row_idx >= gstate.rows.size()) {
		output.SetCardinality(0);
		return;
	}

	size_t count = std::min<size_t>(STANDARD_VECTOR_SIZE, gstate.rows.size() - gstate.current_row_idx);

	auto node_index_data = FlatVector::GetData<int64_t>(output.data[0]);
	auto trait_data = FlatVector::GetData<string_t>(output.data[1]);
	auto state_data = FlatVector::GetData<string_t>(output.data[2]);
	auto in_mpr_data = FlatVector::GetData<bool>(output.data[3]);
	auto min_cost_data = FlatVector::GetData<double>(output.data[4]);

	for (size_t k = 0; k < count; k++) {
		const auto &row = gstate.rows[gstate.current_row_idx + k];
		node_index_data[k] = row.node_index;
		trait_data[k] = StringVector::AddString(output.data[1], row.trait);
		state_data[k] = StringVector::AddString(output.data[2], row.state);
		in_mpr_data[k] = row.in_mpr;
		min_cost_data[k] = row.min_cost;
	}

	gstate.current_row_idx += count;
	output.SetCardinality(count);
}

void PhyloAncestralParsimonyTableFunction::Register(ExtensionLoader &loader) {
	TableFunctionSet set("phylo_ancestral_parsimony");
	set.AddFunction(TableFunction({LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal));
	set.AddFunction(
	    TableFunction({LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal));
	loader.RegisterFunction(set);
}

} // namespace duckdb
