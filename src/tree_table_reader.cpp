#include "tree_table_reader.hpp"
#include "catalog_utils.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include <cmath>
#include <limits>

namespace duckdb {

void ValidateTreeTableSchema(ClientContext &context, const std::string &table_name) {
	auto info = GetTableOrViewColumns(context, table_name, "Tree table");
	auto &col_names = info.names;
	auto &col_types = info.types;

	std::unordered_map<string, idx_t> name_to_idx;
	for (idx_t i = 0; i < col_names.size(); i++) {
		name_to_idx[StringUtil::Lower(col_names[i])] = i;
	}

	// Only check that required columns exist — the SQL query casts to correct types
	auto check_exists = [&](const string &col_name) {
		if (name_to_idx.find(col_name) == name_to_idx.end()) {
			throw BinderException("Tree table '%s' missing required column '%s'", table_name, col_name);
		}
	};

	check_exists("node_index");
	check_exists("parent_index");
}

std::vector<miint::NodeInput> ReadTreeTable(ClientContext &context, const std::string &table_name) {
	std::vector<miint::NodeInput> result;

	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	// Cast to canonical types so GetData<> always matches the backing store.
	// DuckDB handles INTEGER->BIGINT, FLOAT->DOUBLE etc. transparently.
	std::string query = "SELECT node_index::BIGINT, parent_index::BIGINT, name::VARCHAR, "
	                    "branch_length::DOUBLE, edge_id::BIGINT FROM " +
	                    KeywordHelper::WriteOptionallyQuoted(table_name);

	auto query_result = conn.Query(query);

	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from tree table '%s': %s", table_name, query_result->GetError());
	}

	auto &materialized = query_result->Cast<MaterializedQueryResult>();

	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}

		auto &node_index_vec = chunk->data[0];
		auto &parent_index_vec = chunk->data[1];
		auto &name_vec = chunk->data[2];
		auto &branch_length_vec = chunk->data[3];
		auto &edge_id_vec = chunk->data[4];

		UnifiedVectorFormat node_data, parent_data, name_data, bl_data, edge_data;
		node_index_vec.ToUnifiedFormat(chunk->size(), node_data);
		parent_index_vec.ToUnifiedFormat(chunk->size(), parent_data);
		name_vec.ToUnifiedFormat(chunk->size(), name_data);
		branch_length_vec.ToUnifiedFormat(chunk->size(), bl_data);
		edge_id_vec.ToUnifiedFormat(chunk->size(), edge_data);

		auto node_indices = UnifiedVectorFormat::GetData<int64_t>(node_data);
		auto parent_indices = UnifiedVectorFormat::GetData<int64_t>(parent_data);
		auto names = UnifiedVectorFormat::GetData<string_t>(name_data);
		auto branch_lengths = UnifiedVectorFormat::GetData<double>(bl_data);
		auto edge_ids = UnifiedVectorFormat::GetData<int64_t>(edge_data);

		for (idx_t i = 0; i < chunk->size(); i++) {
			miint::NodeInput node;

			auto ni_idx = node_data.sel->get_index(i);
			if (!node_data.validity.RowIsValid(ni_idx)) {
				throw InvalidInputException("NULL node_index in tree table '%s'", table_name);
			}
			node.node_id = node_indices[ni_idx];

			auto pi_idx = parent_data.sel->get_index(i);
			if (parent_data.validity.RowIsValid(pi_idx)) {
				node.parent_id = parent_indices[pi_idx];
			} else {
				node.parent_id = std::nullopt; // Root
			}

			auto nm_idx = name_data.sel->get_index(i);
			if (name_data.validity.RowIsValid(nm_idx)) {
				node.name = names[nm_idx].GetString();
			}

			node.branch_length = std::numeric_limits<double>::quiet_NaN();
			auto bl_idx = bl_data.sel->get_index(i);
			if (bl_data.validity.RowIsValid(bl_idx)) {
				node.branch_length = branch_lengths[bl_idx];
			}

			auto ei_idx = edge_data.sel->get_index(i);
			if (edge_data.validity.RowIsValid(ei_idx)) {
				node.edge_id = edge_ids[ei_idx];
			}

			result.push_back(std::move(node));
		}
	}

	if (result.empty()) {
		throw InvalidInputException("Tree table '%s' is empty", table_name);
	}

	return result;
}

} // namespace duckdb
