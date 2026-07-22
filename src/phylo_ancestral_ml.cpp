#include "phylo_ancestral_ml.hpp"
#include "tree_table_reader.hpp"
#include "phylo_traits_reader.hpp"
#include "NewickTree.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include <cctype>
#include <cmath>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace duckdb {

namespace {

std::string to_upper(std::string s) {
	for (char &c : s) {
		c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
	}
	return s;
}

} // namespace

PhyloAncestralMLTableFunction::Data::Data(std::string tree_table, std::string traits_table)
    : tree_table_name(std::move(tree_table)), traits_table_name(std::move(traits_table)) {
	names = {"node_index", "trait", "state", "probability", "rate", "log_likelihood"};
	types = {LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::VARCHAR,
	         LogicalType::DOUBLE, LogicalType::DOUBLE,  LogicalType::DOUBLE};
}

unique_ptr<FunctionData> PhyloAncestralMLTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                             vector<LogicalType> &return_types,
                                                             vector<std::string> &names) {
	auto tree_table_name = input.inputs[0].ToString();
	auto traits_table_name = input.inputs[1].ToString();

	ValidateTreeTableSchema(context, tree_table_name);
	ValidateTraitsTableSchema(context, traits_table_name);

	auto data = duckdb::make_uniq<Data>(std::move(tree_table_name), std::move(traits_table_name));

	auto model_it = input.named_parameters.find("model");
	if (model_it != input.named_parameters.end() && !model_it->second.IsNull()) {
		data->model = to_upper(model_it->second.ToString());
	}
	if (data->model != "ER" && data->model != "SYM" && data->model != "ARD") {
		throw BinderException("phylo_ancestral_ml: model must be 'ER' (equal rates), 'SYM' (symmetric rates), or "
		                      "'ARD' (all rates different); got '%s'",
		                      data->model);
	}

	auto rate_it = input.named_parameters.find("rate");
	if (rate_it != input.named_parameters.end() && !rate_it->second.IsNull()) {
		if (data->model != "ER") {
			throw BinderException("phylo_ancestral_ml: the 'rate' parameter fixes the single ER rate and does not "
			                      "apply to model '%s' (which fits a full rate matrix); omit it",
			                      data->model);
		}
		data->has_fixed_rate = true;
		data->fixed_rate = rate_it->second.GetValue<double>();
		if (!(data->fixed_rate > 0.0) || !std::isfinite(data->fixed_rate)) {
			throw BinderException("phylo_ancestral_ml: rate must be a finite strictly positive number; got %f",
			                      data->fixed_rate);
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

unique_ptr<GlobalTableFunctionState> PhyloAncestralMLTableFunction::InitGlobal(ClientContext &context,
                                                                               TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = duckdb::make_uniq<GlobalState>();

	auto tree_data = BuildTreeAndNodeIds(context, bind_data.tree_table_name);
	auto &tree = tree_data.tree;
	auto &tree_index_to_node_id = tree_data.node_ids;

	auto traits = ReadDiscreteTraits(context, bind_data.traits_table_name);

	// Build per-trait (tip-state indices, k, alphabet). The ER alphabet is the sorted
	// distinct observed states for that trait.
	std::vector<std::string> trait_names;
	std::vector<std::unordered_map<std::string, uint32_t>> tip_states_list;
	std::vector<uint32_t> k_list;
	std::vector<std::vector<std::string>> alphabet_list;

	for (auto &entry : traits) {
		const std::string &trait = entry.first;
		const auto &tipmap = entry.second;

		std::set<std::string> distinct;
		for (const auto &[tip, st] : tipmap) {
			distinct.insert(st);
		}
		std::vector<std::string> alphabet(distinct.begin(), distinct.end()); // std::set -> sorted
		std::unordered_map<std::string, uint32_t> index;
		for (uint32_t i = 0; i < alphabet.size(); i++) {
			index[alphabet[i]] = i;
		}
		const uint32_t k = static_cast<uint32_t>(alphabet.size());

		std::unordered_map<std::string, uint32_t> tip_states;
		for (const auto &[tip, st] : tipmap) {
			tip_states[tip] = index.at(st);
		}

		trait_names.push_back(trait);
		tip_states_list.push_back(std::move(tip_states));
		k_list.push_back(k);
		alphabet_list.push_back(std::move(alphabet));
	}

	std::optional<double> fixed_rate;
	if (bind_data.has_fixed_rate) {
		fixed_rate = bind_data.fixed_rate;
	}

	const std::string &model = bind_data.model;
	const bool is_matrix = (model != "ER"); // SYM/ARD fit a rate matrix -> `rate` column is NULL
	auto run_all = [&]() {
		if (model == "SYM") {
			return tree.ancestral_ml_sym(tip_states_list, k_list, &tree_index_to_node_id);
		}
		if (model == "ARD") {
			return tree.ancestral_ml_ard(tip_states_list, k_list, &tree_index_to_node_id);
		}
		return tree.ancestral_ml(tip_states_list, k_list, fixed_rate, &tree_index_to_node_id);
	};
	auto run_one = [&](size_t tt) {
		if (model == "SYM") {
			return tree.ancestral_ml_sym(tip_states_list[tt], k_list[tt], &tree_index_to_node_id);
		}
		if (model == "ARD") {
			return tree.ancestral_ml_ard(tip_states_list[tt], k_list[tt], &tree_index_to_node_id);
		}
		return tree.ancestral_ml(tip_states_list[tt], k_list[tt], fixed_rate, &tree_index_to_node_id);
	};

	std::vector<miint::AncestralMLResult> per_trait;
	try {
		per_trait = run_all();
	} catch (const std::exception &e) {
		// Cold path: re-run per trait to attribute a completeness/branch failure to the
		// specific trait (structural errors are trait-independent and surface first).
		for (size_t t = 0; t < tip_states_list.size(); t++) {
			try {
				run_one(t);
			} catch (const std::exception &inner) {
				throw InvalidInputException("phylo_ancestral_ml failed for trait '%s': %s", trait_names[t],
				                            inner.what());
			}
		}
		throw InvalidInputException("phylo_ancestral_ml failed: %s", e.what());
	}

	for (size_t t = 0; t < trait_names.size(); t++) {
		const auto &res = per_trait[t];
		for (const auto &st : res.states) {
			// SYM/ARD fit a rate matrix, not a scalar -> report the `rate` column as NULL.
			gstate->rows.push_back(StateRow {tree_index_to_node_id[st.node], trait_names[t], alphabet_list[t][st.state],
			                                 st.probability, res.rate, is_matrix, res.log_likelihood});
		}
	}

	return gstate;
}

void PhyloAncestralMLTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	if (gstate.current_row_idx >= gstate.rows.size()) {
		output.SetCardinality(0);
		return;
	}

	size_t count = std::min<size_t>(STANDARD_VECTOR_SIZE, gstate.rows.size() - gstate.current_row_idx);

	auto node_index_data = FlatVector::GetData<int64_t>(output.data[0]);
	auto trait_data = FlatVector::GetData<string_t>(output.data[1]);
	auto state_data = FlatVector::GetData<string_t>(output.data[2]);
	auto prob_data = FlatVector::GetData<double>(output.data[3]);
	auto rate_data = FlatVector::GetData<double>(output.data[4]);
	auto logl_data = FlatVector::GetData<double>(output.data[5]);

	for (size_t k = 0; k < count; k++) {
		const auto &row = gstate.rows[gstate.current_row_idx + k];
		node_index_data[k] = row.node_index;
		trait_data[k] = StringVector::AddString(output.data[1], row.trait);
		state_data[k] = StringVector::AddString(output.data[2], row.state);
		prob_data[k] = row.probability;
		rate_data[k] = row.rate;
		if (row.rate_null) {
			FlatVector::SetNull(output.data[4], k, true);
		}
		logl_data[k] = row.log_likelihood;
	}

	gstate.current_row_idx += count;
	output.SetCardinality(count);
}

TableFunction PhyloAncestralMLTableFunction::GetFunction() {
	TableFunction tf("phylo_ancestral_ml", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal);
	tf.named_parameters["model"] = LogicalType::VARCHAR;
	tf.named_parameters["rate"] = LogicalType::DOUBLE;
	return tf;
}

void PhyloAncestralMLTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
