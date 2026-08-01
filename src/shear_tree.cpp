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

	std::string query = "SELECT name::VARCHAR FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);
	auto query_result = conn.Query(query);
	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from tips table '%s': %s", table_name, query_result->GetError());
	}

	auto &materialized = query_result->Cast<MaterializedQueryResult>();
	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}
		UnifiedVectorFormat name_data;
		chunk->data[0].ToUnifiedFormat(chunk->size(), name_data);
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

bool GetBoolParam(TableFunctionBindInput &input, const std::string &key, bool default_value) {
	auto it = input.named_parameters.find(key);
	if (it == input.named_parameters.end() || it->second.IsNull()) {
		return default_value;
	}
	return BooleanValue::Get(it->second);
}

} // namespace

ShearTreeTableFunction::Data::Data(std::string tree_table, std::string tips_table, bool collapse_, bool ignore_missing_)
    : tree_table_name(std::move(tree_table)), tips_table_name(std::move(tips_table)), collapse(collapse_),
      ignore_missing(ignore_missing_) {
	ReadNewickTableFunction::GetSchema(names, types, false);
}

unique_ptr<FunctionData> ShearTreeTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                      vector<LogicalType> &return_types, vector<std::string> &names) {
	auto tree_table_name = input.inputs[0].ToString();
	auto tips_table_name = input.inputs[1].ToString();

	// Validate schemas at bind time for early error detection.
	ValidateTreeTableSchema(context, tree_table_name);
	ValidateTipsTableSchema(context, tips_table_name);

	bool collapse = GetBoolParam(input, "collapse", true);
	bool ignore_missing = GetBoolParam(input, "ignore_missing", false);

	auto data =
	    duckdb::make_uniq<Data>(std::move(tree_table_name), std::move(tips_table_name), collapse, ignore_missing);

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> ShearTreeTableFunction::InitGlobal(ClientContext &context,
                                                                        TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = duckdb::make_uniq<GlobalState>();

	auto node_inputs = ReadTreeTable(context, bind_data.tree_table_name);

	miint::NewickTree tree;
	try {
		tree = miint::NewickTree::build(node_inputs);
	} catch (const std::exception &e) {
		throw InvalidInputException("Failed to build tree from '%s': %s", bind_data.tree_table_name, e.what());
	}

	// The NodeInput vector is not needed once the tree is built; release it
	// before shearing so it doesn't add to peak memory. shear() itself
	// allocates an intermediate NodeInput vector plus the result tree, so for
	// large trees this trims roughly one full copy off the transient peak.
	node_inputs = {};

	auto keep_names = ReadTipNames(context, bind_data.tips_table_name);

	miint::NewickTree sheared;
	try {
		sheared = tree.shear(keep_names, bind_data.collapse, bind_data.ignore_missing);
	} catch (const std::exception &e) {
		throw InvalidInputException("shear_tree failed: %s", e.what());
	}

	gstate->rows = ReadNewickTableFunction::TreeToRows(sheared);

	return gstate;
}

void ShearTreeTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
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

TableFunction ShearTreeTableFunction::GetFunction() {
	TableFunction tf("shear_tree", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal);
	tf.named_parameters["collapse"] = LogicalType::BOOLEAN;
	tf.named_parameters["ignore_missing"] = LogicalType::BOOLEAN;
	return tf;
}

void ShearTreeTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
