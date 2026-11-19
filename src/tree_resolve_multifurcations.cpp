#include "tree_resolve_multifurcations.hpp"
#include "catalog_utils.hpp"
#include "tree_table_reader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/data_chunk.hpp"

namespace duckdb {

unique_ptr<FunctionData> TreeResolveMultifurcationsTableFunction::Bind(ClientContext &context,
                                                                       TableFunctionBindInput &input,
                                                                       vector<LogicalType> &return_types,
                                                                       vector<Identifier> &names) {
	auto tree_table_name = input.inputs[0].ToString();
	RejectCTERelationName(input, tree_table_name);

	// Validate schema at bind time for early error detection.
	auto tree_columns = BindTreeOutput(context, "tree_resolve_multifurcations", tree_table_name, return_types, names);
	return duckdb::make_uniq<Data>(std::move(tree_table_name), std::move(tree_columns));
}

unique_ptr<GlobalTableFunctionState>
TreeResolveMultifurcationsTableFunction::InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = duckdb::make_uniq<TreeOutputState>();

	auto tree = BuildTree(context, bind_data.tree_table_name, bind_data.tree_columns, gstate->labels);

	miint::NewickTree resolved;
	try {
		resolved = tree.resolve_multifurcations();
	} catch (const std::exception &e) {
		throw InvalidInputException("tree_resolve_multifurcations failed: %s", e.what());
	}

	// resolve_multifurcations keeps input node i at index i and appends the connectors
	gstate->SetRows(resolved);
	return gstate;
}

TableFunction TreeResolveMultifurcationsTableFunction::GetFunction() {
	return TableFunction("tree_resolve_multifurcations", {LogicalType::VARCHAR}, TreeOutputExecute, Bind, InitGlobal);
}

void TreeResolveMultifurcationsTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
