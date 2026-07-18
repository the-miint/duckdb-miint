#include "tree_resolve_multifurcations.hpp"
#include "tree_table_reader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/data_chunk.hpp"

namespace duckdb {

TreeResolveMultifurcationsTableFunction::Data::Data(std::string tree_table) : tree_table_name(std::move(tree_table)) {
	ReadNewickTableFunction::GetSchema(names, types, false);
}

unique_ptr<FunctionData> TreeResolveMultifurcationsTableFunction::Bind(ClientContext &context,
                                                                       TableFunctionBindInput &input,
                                                                       vector<LogicalType> &return_types,
                                                                       vector<std::string> &names) {
	auto tree_table_name = input.inputs[0].ToString();

	// Validate schema at bind time for early error detection.
	ValidateTreeTableSchema(context, tree_table_name);

	auto data = duckdb::make_uniq<Data>(std::move(tree_table_name));

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState>
TreeResolveMultifurcationsTableFunction::InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = duckdb::make_uniq<GlobalState>();

	auto node_inputs = ReadTreeTable(context, bind_data.tree_table_name);

	miint::NewickTree tree;
	try {
		tree = miint::NewickTree::build(node_inputs);
	} catch (const std::exception &e) {
		throw InvalidInputException("Failed to build tree from '%s': %s", bind_data.tree_table_name, e.what());
	}

	// Release the NodeInput vector before resolving to trim transient peak
	// memory (see shear_tree.cpp for the same rationale).
	node_inputs = {};

	miint::NewickTree resolved;
	try {
		resolved = tree.resolve_multifurcations();
	} catch (const std::exception &e) {
		throw InvalidInputException("tree_resolve_multifurcations failed: %s", e.what());
	}

	gstate->rows = ReadNewickTableFunction::TreeToRows(resolved);

	return gstate;
}

void TreeResolveMultifurcationsTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p,
                                                      DataChunk &output) {
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	if (global_state.current_row_idx >= global_state.rows.size()) {
		output.SetCardinality(0);
		return;
	}

	size_t rows_to_output =
	    std::min<size_t>(STANDARD_VECTOR_SIZE, global_state.rows.size() - global_state.current_row_idx);

	ReadNewickTableFunction::EmitNodeRows(global_state.rows, global_state.current_row_idx, rows_to_output, output, 0,
	                                      false, "");

	global_state.current_row_idx += rows_to_output;
	output.SetCardinality(rows_to_output);
}

TableFunction TreeResolveMultifurcationsTableFunction::GetFunction() {
	return TableFunction("tree_resolve_multifurcations", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal);
}

void TreeResolveMultifurcationsTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
