#include "tree_resolve_placement.hpp"
#include "catalog_utils.hpp"
#include "placement_table_reader.hpp"
#include "tree_function_utils.hpp"
#include "tree_table_reader.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

TreeResolvePlacementTableFunction::Data::Data(std::string tree_table, TreeNameColumns tree_columns,
                                              std::string placements_table)
    : TreeOutputData(std::move(tree_table), std::move(tree_columns)),
      placements_table_name(std::move(placements_table)) {
}

unique_ptr<FunctionData> TreeResolvePlacementTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                                 vector<LogicalType> &return_types,
                                                                 vector<Identifier> &names) {
	auto tree_table_name = input.inputs[0].ToString();
	auto placements_table_name = input.inputs[1].ToString();
	RejectCTERelationName(input, tree_table_name);
	RejectCTERelationName(input, placements_table_name);

	// Validate schemas at bind time for early error detection
	auto tree_columns = BindTreeOutput(context, "tree_resolve_placement", tree_table_name, return_types, names);
	ValidatePlacementTableSchema(context, placements_table_name);

	// A fragment's id names the tip it becomes. A VARCHAR tree takes any id type and names the tip
	// with its text; a BIGINT or UUID tree takes its own type, so the tip's name is one.
	if (tree_columns.name_type.id() != LogicalTypeId::VARCHAR) {
		// ValidatePlacementTableSchema, above, has checked that fragment_id exists
		LogicalType fragment_type;
		FindColumnType(GetTableOrViewColumns(context, placements_table_name, "Placement table"), "fragment_id",
		               fragment_type);
		if (fragment_type.id() != tree_columns.name_type.id()) {
			throw BinderException("tree_resolve_placement: the tree's names are %s, but 'fragment_id' is %s; a "
			                      "BIGINT or UUID tree takes fragment ids of its own type, which name its new tips",
			                      tree_columns.name_type.ToString(), fragment_type.ToString());
		}
	}

	return duckdb::make_uniq<Data>(std::move(tree_table_name), std::move(tree_columns),
	                               std::move(placements_table_name));
}

unique_ptr<GlobalTableFunctionState> TreeResolvePlacementTableFunction::InitGlobal(ClientContext &context,
                                                                                   TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = duckdb::make_uniq<TreeOutputState>();

	auto tree = BuildTree(context, bind_data.tree_table_name, bind_data.tree_columns, gstate->labels);
	auto placements = ReadPlacementTable(context, bind_data.placements_table_name);

	// Resolve placements
	try {
		tree.insert_fully_resolved(placements);
	} catch (const std::runtime_error &e) {
		throw InvalidInputException("Failed to resolve placements: %s", e.what());
	}

	// insert_fully_resolved keeps input node i at index i and appends the nodes it creates
	gstate->SetRows(tree);
	return gstate;
}

TableFunction TreeResolvePlacementTableFunction::GetFunction() {
	return TableFunction("tree_resolve_placement", {LogicalType::VARCHAR, LogicalType::VARCHAR}, TreeOutputExecute,
	                     Bind, InitGlobal);
}

void TreeResolvePlacementTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
