#include "phylo_independent_contrasts.hpp"
#include "tree_table_reader.hpp"
#include "phylo_traits_reader.hpp"
#include "NewickTree.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace duckdb {

PhyloIndependentContrastsTableFunction::Data::Data(std::string tree_table, std::string traits_table)
    : tree_table_name(std::move(tree_table)), traits_table_name(std::move(traits_table)) {
	names = {"node_index", "trait", "contrast", "ancestral_estimate", "contrast_variance"};
	types = {LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::DOUBLE, LogicalType::DOUBLE, LogicalType::DOUBLE};
}

unique_ptr<FunctionData> PhyloIndependentContrastsTableFunction::Bind(ClientContext &context,
                                                                      TableFunctionBindInput &input,
                                                                      vector<LogicalType> &return_types,
                                                                      vector<std::string> &names) {
	auto tree_table_name = input.inputs[0].ToString();
	auto traits_table_name = input.inputs[1].ToString();

	// Validate schemas at bind time for early error detection.
	ValidateTreeTableSchema(context, tree_table_name);
	ValidateTraitsTableSchema(context, traits_table_name);

	auto data = duckdb::make_uniq<Data>(std::move(tree_table_name), std::move(traits_table_name));

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> PhyloIndependentContrastsTableFunction::InitGlobal(ClientContext &context,
                                                                                        TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = duckdb::make_uniq<GlobalState>();

	auto tree_data = BuildTreeAndNodeIds(context, bind_data.tree_table_name);
	auto &tree = tree_data.tree;
	auto &tree_index_to_node_id = tree_data.node_ids;

	auto traits = ReadContinuousTraits(context, bind_data.traits_table_name);

	// Flatten to parallel name/map vectors (traits is sorted, so output is stable)
	// and compute all traits in one batch: the tree validation and variance
	// precompute happen once, then each trait is a cheap post-order pass.
	std::vector<std::string> trait_names;
	std::vector<std::unordered_map<std::string, double>> trait_maps;
	trait_names.reserve(traits.size());
	trait_maps.reserve(traits.size());
	for (auto &entry : traits) {
		trait_names.push_back(entry.first);
		trait_maps.push_back(std::move(entry.second));
	}

	std::vector<std::vector<miint::IndependentContrast>> per_trait;
	try {
		per_trait = tree.independent_contrasts(trait_maps, &tree_index_to_node_id);
	} catch (const std::exception &e) {
		// Cold path: re-run per trait to attribute a completeness failure to the
		// specific trait for a precise message (structural/branch errors are
		// trait-independent and surface on the first trait).
		for (size_t t = 0; t < trait_maps.size(); t++) {
			try {
				tree.independent_contrasts(trait_maps[t], &tree_index_to_node_id);
			} catch (const std::exception &inner) {
				throw InvalidInputException("phylo_independent_contrasts failed for trait '%s': %s", trait_names[t],
				                            inner.what());
			}
		}
		throw InvalidInputException("phylo_independent_contrasts failed: %s", e.what());
	}

	for (size_t t = 0; t < trait_names.size(); t++) {
		for (const auto &c : per_trait[t]) {
			gstate->rows.push_back(ContrastRow {tree_index_to_node_id[c.node], trait_names[t], c.contrast,
			                                    c.ancestral_estimate, c.contrast_variance});
		}
	}

	return gstate;
}

void PhyloIndependentContrastsTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p,
                                                     DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	if (gstate.current_row_idx >= gstate.rows.size()) {
		output.SetCardinality(0);
		return;
	}

	size_t count = std::min<size_t>(STANDARD_VECTOR_SIZE, gstate.rows.size() - gstate.current_row_idx);

	auto node_index_data = FlatVector::GetData<int64_t>(output.data[0]);
	auto trait_data = FlatVector::GetData<string_t>(output.data[1]);
	auto contrast_data = FlatVector::GetData<double>(output.data[2]);
	auto anc_data = FlatVector::GetData<double>(output.data[3]);
	auto var_data = FlatVector::GetData<double>(output.data[4]);

	for (size_t k = 0; k < count; k++) {
		const auto &row = gstate.rows[gstate.current_row_idx + k];
		node_index_data[k] = row.node_index;
		trait_data[k] = StringVector::AddString(output.data[1], row.trait);
		contrast_data[k] = row.contrast;
		anc_data[k] = row.ancestral_estimate;
		var_data[k] = row.contrast_variance;
	}

	gstate.current_row_idx += count;
	output.SetCardinality(count);
}

TableFunction PhyloIndependentContrastsTableFunction::GetFunction() {
	return TableFunction("phylo_independent_contrasts", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind,
	                     InitGlobal);
}

void PhyloIndependentContrastsTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
