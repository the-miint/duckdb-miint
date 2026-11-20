#include "shear_tree.hpp"
#include "tree_table_reader.hpp"
#include "catalog_utils.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/parser/keyword_helper.hpp"
#include <unordered_set>
#include "miint_named_parameter.hpp"

namespace duckdb {

namespace {

// Validate that the tips table/view has the required `name` column.
void ValidateTipsTableSchema(ClientContext &context, const std::string &table_name) {
	auto info = GetTableOrViewColumns(context, table_name, "Tips table");
	if (!HasColumn(info, "name")) {
		throw BinderException("Tips table '%s' missing required column 'name'", table_name);
	}
}

// Read the set of tip names to keep from the tips table/view. Uses a separate
// Connection to avoid re-entering the current context (see
// docs/internals/reading-tables-views.md). NULL names are skipped.
std::unordered_set<std::string> ReadTipNames(ClientContext &context, const std::string &table_name) {
	std::unordered_set<std::string> names;

	auto conn = MakeReadOnlyHelperConnection(context);

	std::string query = "SELECT name::VARCHAR FROM " + SQLIdentifier::ToString(table_name);
	auto query_result = conn.Query(query);
	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from tips table '%s': %s", table_name, query_result->GetError());
	}

	auto &materialized = *query_result;
	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}
		UnifiedVectorFormat name_data;
		chunk->data[0].ToUnifiedFormat(name_data);
		auto name_strs = UnifiedVectorFormat::GetData<string_t>(name_data);
		for (idx_t i = 0; i < chunk->size(); i++) {
			auto idx = name_data.sel->get_index(i);
			if (name_data.validity.RowIsValid(idx)) {
				names.insert(name_strs[idx].GetString());
			}
		}
	}

	return names;
}

} // namespace

ShearTreeTableFunction::Data::Data(std::string tree_table, TreeNameColumns tree_columns, std::string tips_table,
                                   bool collapse_, bool ignore_missing_)
    : TreeOutputData(std::move(tree_table), std::move(tree_columns)), tips_table_name(std::move(tips_table)),
      collapse(collapse_), ignore_missing(ignore_missing_) {
}

unique_ptr<FunctionData> ShearTreeTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                      vector<LogicalType> &return_types, vector<Identifier> &names) {
	auto tree_table_name = input.inputs[0].ToString();
	auto tips_table_name = input.inputs[1].ToString();
	RejectCTERelationName(input, tree_table_name);
	RejectCTERelationName(input, tips_table_name);

	// Validate schemas at bind time for early error detection.
	auto tree_columns = BindTreeOutput(context, "shear_tree", tree_table_name, return_types, names);
	ValidateTipsTableSchema(context, tips_table_name);

	bool collapse = GetBoolParam(input, "collapse", true);
	bool ignore_missing = GetBoolParam(input, "ignore_missing", false);

	return duckdb::make_uniq<Data>(std::move(tree_table_name), std::move(tree_columns), std::move(tips_table_name),
	                               collapse, ignore_missing);
}

unique_ptr<GlobalTableFunctionState> ShearTreeTableFunction::InitGlobal(ClientContext &context,
                                                                        TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = duckdb::make_uniq<TreeOutputState>();

	auto tree = BuildTree(context, "shear_tree", bind_data.tree_table_name, bind_data.tree_columns, gstate->labels);
	auto keep_names = ReadTipNames(context, bind_data.tips_table_name);

	miint::NewickTree sheared;
	std::vector<uint32_t> source;
	try {
		sheared = tree.shear(keep_names, bind_data.collapse, bind_data.ignore_missing, &source);
	} catch (const std::exception &e) {
		throw InvalidInputException("shear_tree failed: %s", e.what());
	}

	// Each kept node keeps its label; a node the collapse removed takes its label with it
	gstate->SetRows(sheared, &source);
	return gstate;
}

TableFunction ShearTreeTableFunction::GetFunction() {
	TableFunction tf("shear_tree", {LogicalType::VARCHAR, LogicalType::VARCHAR}, TreeOutputExecute, Bind, InitGlobal);
	AddNamedParameter(tf, "collapse", LogicalType::BOOLEAN);
	AddNamedParameter(tf, "ignore_missing", LogicalType::BOOLEAN);
	return tf;
}

void ShearTreeTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
