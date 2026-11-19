#include "tree_resolve_placement.hpp"
#include "catalog_utils.hpp"
#include "placement_table_reader.hpp"
#include "tree_function_utils.hpp"
#include "tree_table_reader.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

TreeResolvePlacementTableFunction::Data::Data(std::string tree_table, std::string placements_table,
                                              TreeNameColumns tree_columns_)
    : tree_table_name(std::move(tree_table)), placements_table_name(std::move(placements_table)),
      tree_columns(std::move(tree_columns_)) {
	tree_columns.OutputSchema(names, types);
}

unique_ptr<FunctionData> TreeResolvePlacementTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                                 vector<LogicalType> &return_types,
                                                                 vector<Identifier> &names) {
	auto tree_table_name = input.inputs[0].ToString();
	auto placements_table_name = input.inputs[1].ToString();
	RejectCTERelationName(input, tree_table_name);
	RejectCTERelationName(input, placements_table_name);

	// Validate schemas at bind time for early error detection
	ValidateTreeTableSchema(context, tree_table_name);
	ValidatePlacementTableSchema(context, placements_table_name);
	auto tree_columns = BindTreeNameColumns(context, "tree_resolve_placement", tree_table_name);

	// A fragment's id names the tip it becomes. A VARCHAR tree takes any id type and names the tip
	// with its text; a BIGINT or UUID tree takes its own type, so the tip's name is one.
	if (tree_columns.name_type.id() != LogicalTypeId::VARCHAR) {
		LogicalType fragment_type;
		FindColumnType(GetTableOrViewColumns(context, placements_table_name, "Placement table"), "fragment_id",
		               fragment_type);
		if (fragment_type.id() != tree_columns.name_type.id()) {
			throw BinderException("tree_resolve_placement: the tree's names are %s, but 'fragment_id' is %s; a "
			                      "BIGINT or UUID tree takes fragment ids of its own type, which name its new tips",
			                      tree_columns.name_type.ToString(), fragment_type.ToString());
		}
	}

	auto data =
	    duckdb::make_uniq<Data>(std::move(tree_table_name), std::move(placements_table_name), std::move(tree_columns));

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> TreeResolvePlacementTableFunction::InitGlobal(ClientContext &context,
                                                                                   TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = duckdb::make_uniq<GlobalState>();

	// Read tree data and placements
	auto node_inputs = ReadTreeTable(context, bind_data.tree_table_name, bind_data.tree_columns, gstate->labels);
	auto placements = ReadPlacementTable(context, bind_data.placements_table_name);

	// Build tree from node data
	miint::NewickTree tree;
	try {
		tree = miint::NewickTree::build(node_inputs);
	} catch (const std::exception &e) {
		throw InvalidInputException("Failed to build tree from '%s': %s", bind_data.tree_table_name, e.what());
	}

	// Resolve placements
	try {
		tree.insert_fully_resolved(placements);
	} catch (const std::runtime_error &e) {
		throw InvalidInputException("Failed to resolve placements: %s", e.what());
	}

	// Convert to rows
	gstate->rows = ReadNewickTableFunction::TreeToRows(tree);
	// insert_fully_resolved keeps input node i at index i and appends the nodes it creates, which
	// have no label
	if (bind_data.tree_columns.has_label) {
		gstate->labels.resize(tree.num_nodes());
	}

	return gstate;
}

void TreeResolvePlacementTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	if (global_state.current_row_idx >= global_state.rows.size()) {
		output.SetChildCardinality(0);
		return;
	}

	size_t rows_to_output =
	    std::min<size_t>(STANDARD_VECTOR_SIZE, global_state.rows.size() - global_state.current_row_idx);

	ReadNewickTableFunction::EmitNodeRows(global_state.rows, global_state.current_row_idx, rows_to_output, output, 0,
	                                      false, "", bind_data.tree_columns.name_type,
	                                      bind_data.tree_columns.has_label ? &global_state.labels : nullptr);

	global_state.current_row_idx += rows_to_output;
	output.SetChildCardinality(rows_to_output);
}

TableFunction TreeResolvePlacementTableFunction::GetFunction() {
	return TableFunction("tree_resolve_placement", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind,
	                     InitGlobal);
}

void TreeResolvePlacementTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
