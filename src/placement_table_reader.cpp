#include "placement_table_reader.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include <limits>

namespace duckdb {

// Helper to get column names and types from either a table or view
static void GetTableOrViewColumns(ClientContext &context, const std::string &table_name,
                                  vector<string> &out_names, vector<LogicalType> &out_types) {
	// Use TABLE_ENTRY lookup which can return either tables or views
	EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY, table_name, QueryErrorContext());
	auto entry = Catalog::GetEntry(context, INVALID_CATALOG, INVALID_SCHEMA, lookup_info, OnEntryNotFound::RETURN_NULL);

	if (!entry) {
		throw BinderException("Placement table or view '%s' does not exist", table_name);
	}

	if (entry->type == CatalogType::TABLE_ENTRY) {
		auto &table = entry->Cast<TableCatalogEntry>();
		auto &columns = table.GetColumns();
		for (idx_t i = 0; i < columns.LogicalColumnCount(); i++) {
			auto &col = columns.GetColumn(LogicalIndex(i));
			out_names.push_back(col.Name());
			out_types.push_back(col.Type());
		}
	} else if (entry->type == CatalogType::VIEW_ENTRY) {
		auto &view = entry->Cast<ViewCatalogEntry>();
		out_names = view.names;
		out_types = view.types;
	} else {
		throw BinderException("'%s' is not a table or view", table_name);
	}
}

void ValidatePlacementTableSchema(ClientContext &context, const std::string &table_name) {
	vector<string> col_names;
	vector<LogicalType> col_types;
	GetTableOrViewColumns(context, table_name, col_names, col_types);

	// Build name-to-index map (case-insensitive)
	std::unordered_map<string, idx_t> name_to_idx;
	for (idx_t i = 0; i < col_names.size(); i++) {
		name_to_idx[StringUtil::Lower(col_names[i])] = i;
	}

	// Check required columns exist and have correct types
	auto check_column = [&](const string &col_name, const vector<LogicalTypeId> &allowed_types,
	                        const string &type_desc) {
		auto it = name_to_idx.find(col_name);
		if (it == name_to_idx.end()) {
			throw BinderException("Placement table '%s' missing required column '%s'", table_name, col_name);
		}
		auto &col_type = col_types[it->second];
		bool valid = false;
		for (auto &allowed : allowed_types) {
			if (col_type.id() == allowed) {
				valid = true;
				break;
			}
		}
		if (!valid) {
			throw BinderException("Column '%s' in table '%s' must be %s", col_name, table_name, type_desc);
		}
	};

	check_column("fragment_id", {LogicalTypeId::VARCHAR}, "VARCHAR");
	check_column("edge_id",
	             {LogicalTypeId::BIGINT, LogicalTypeId::INTEGER, LogicalTypeId::UBIGINT, LogicalTypeId::UINTEGER},
	             "an integer type");
	check_column("like_weight_ratio", {LogicalTypeId::DOUBLE, LogicalTypeId::FLOAT}, "DOUBLE or FLOAT");
	check_column("distal_length", {LogicalTypeId::DOUBLE, LogicalTypeId::FLOAT}, "DOUBLE or FLOAT");
	check_column("pendant_length", {LogicalTypeId::DOUBLE, LogicalTypeId::FLOAT}, "DOUBLE or FLOAT");
}

std::vector<miint::Placement> ReadPlacementTable(ClientContext &context, const std::string &table_name) {
	std::vector<miint::Placement> result;

	// Create a new connection to avoid deadlocking the current context.
	// This is necessary because context.Query() requires locking the context,
	// but during bind/finalize the context is already locked.
	// Using a separate connection allows us to execute queries for both tables and views.
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	// Execute a query to read from the table/view
	// This approach works for both tables and views uniformly
	std::string query = "SELECT fragment_id, edge_id, like_weight_ratio, distal_length, pendant_length FROM " +
	                    KeywordHelper::WriteOptionallyQuoted(table_name);

	auto query_result = conn.Query(query);

	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from placement table '%s': %s", table_name,
		                            query_result->GetError());
	}

	auto &materialized = query_result->Cast<MaterializedQueryResult>();
	auto &result_types = materialized.types;

	// Get column types for proper data extraction
	// Column order: fragment_id(0), edge_id(1), like_weight_ratio(2), distal_length(3), pendant_length(4)
	auto &edge_id_type = result_types[1];
	auto &lwr_type = result_types[2];
	auto &distal_type = result_types[3];
	auto &pendant_type = result_types[4];

	idx_t row_number = 0;

	// Process all chunks from the result
	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}

		auto &fragment_id_vec = chunk->data[0];
		auto &edge_id_vec = chunk->data[1];
		auto &lwr_vec = chunk->data[2];
		auto &distal_vec = chunk->data[3];
		auto &pendant_vec = chunk->data[4];

		// Convert to unified format for proper NULL handling
		UnifiedVectorFormat fragment_data, edge_data, lwr_data, distal_data, pendant_data;
		fragment_id_vec.ToUnifiedFormat(chunk->size(), fragment_data);
		edge_id_vec.ToUnifiedFormat(chunk->size(), edge_data);
		lwr_vec.ToUnifiedFormat(chunk->size(), lwr_data);
		distal_vec.ToUnifiedFormat(chunk->size(), distal_data);
		pendant_vec.ToUnifiedFormat(chunk->size(), pendant_data);

		auto fragment_ids = UnifiedVectorFormat::GetData<string_t>(fragment_data);

		for (idx_t i = 0; i < chunk->size(); i++) {
			row_number++;

			auto frag_idx = fragment_data.sel->get_index(i);
			auto edge_idx = edge_data.sel->get_index(i);
			auto lwr_idx = lwr_data.sel->get_index(i);
			auto distal_idx = distal_data.sel->get_index(i);
			auto pendant_idx = pendant_data.sel->get_index(i);

			// Skip rows with NULL fragment_id
			if (!fragment_data.validity.RowIsValid(frag_idx)) {
				continue;
			}

			// Check for NULL values in required columns
			if (!edge_data.validity.RowIsValid(edge_idx)) {
				throw InvalidInputException("NULL edge_id at row %llu in placement table '%s'", row_number,
				                            table_name);
			}
			if (!lwr_data.validity.RowIsValid(lwr_idx)) {
				throw InvalidInputException("NULL like_weight_ratio at row %llu in placement table '%s'", row_number,
				                            table_name);
			}
			if (!distal_data.validity.RowIsValid(distal_idx)) {
				throw InvalidInputException("NULL distal_length at row %llu in placement table '%s'", row_number,
				                            table_name);
			}
			if (!pendant_data.validity.RowIsValid(pendant_idx)) {
				throw InvalidInputException("NULL pendant_length at row %llu in placement table '%s'", row_number,
				                            table_name);
			}

			miint::Placement p;
			p.fragment_id = fragment_ids[frag_idx].GetString();

			// Get edge_id (handle different integer types)
			switch (edge_id_type.id()) {
			case LogicalTypeId::BIGINT:
				p.edge_id = UnifiedVectorFormat::GetData<int64_t>(edge_data)[edge_idx];
				break;
			case LogicalTypeId::INTEGER:
				p.edge_id = UnifiedVectorFormat::GetData<int32_t>(edge_data)[edge_idx];
				break;
			case LogicalTypeId::UBIGINT: {
				uint64_t uval = UnifiedVectorFormat::GetData<uint64_t>(edge_data)[edge_idx];
				if (uval > static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
					throw InvalidInputException("edge_id %llu at row %llu exceeds INT64_MAX in placement table '%s'",
					                            uval, row_number, table_name);
				}
				p.edge_id = static_cast<int64_t>(uval);
				break;
			}
			case LogicalTypeId::UINTEGER:
				p.edge_id = UnifiedVectorFormat::GetData<uint32_t>(edge_data)[edge_idx];
				break;
			default:
				throw InvalidInputException("Unsupported integer type for edge_id");
			}

			// Get like_weight_ratio
			if (lwr_type.id() == LogicalTypeId::DOUBLE) {
				p.like_weight_ratio = UnifiedVectorFormat::GetData<double>(lwr_data)[lwr_idx];
			} else {
				p.like_weight_ratio = UnifiedVectorFormat::GetData<float>(lwr_data)[lwr_idx];
			}

			// Get distal_length
			if (distal_type.id() == LogicalTypeId::DOUBLE) {
				p.distal_length = UnifiedVectorFormat::GetData<double>(distal_data)[distal_idx];
			} else {
				p.distal_length = UnifiedVectorFormat::GetData<float>(distal_data)[distal_idx];
			}

			// Get pendant_length
			if (pendant_type.id() == LogicalTypeId::DOUBLE) {
				p.pendant_length = UnifiedVectorFormat::GetData<double>(pendant_data)[pendant_idx];
			} else {
				p.pendant_length = UnifiedVectorFormat::GetData<float>(pendant_data)[pendant_idx];
			}

			result.push_back(std::move(p));
		}
	}

	if (result.empty()) {
		throw InvalidInputException("Placement table '%s' is empty", table_name);
	}

	return result;
}

} // namespace duckdb
