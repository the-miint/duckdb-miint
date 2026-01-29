#include "reference_table_reader.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include <limits>

namespace duckdb {

// Helper to validate table/view schema for reference lengths
static void ValidateReferenceTableSchema(ClientContext &context, const std::string &table_name) {
	// Use TABLE_ENTRY lookup which can return either tables or views
	EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY, table_name, QueryErrorContext());
	auto entry = Catalog::GetEntry(context, INVALID_CATALOG, INVALID_SCHEMA, lookup_info, OnEntryNotFound::RETURN_NULL);

	if (!entry) {
		throw InvalidInputException("Reference table or view '%s' does not exist", table_name);
	}

	vector<string> col_names;
	vector<LogicalType> col_types;

	if (entry->type == CatalogType::TABLE_ENTRY) {
		auto &table = entry->Cast<TableCatalogEntry>();
		auto &columns = table.GetColumns();
		for (idx_t i = 0; i < columns.LogicalColumnCount(); i++) {
			auto &col = columns.GetColumn(LogicalIndex(i));
			col_names.push_back(col.Name());
			col_types.push_back(col.Type());
		}
	} else if (entry->type == CatalogType::VIEW_ENTRY) {
		auto &view = entry->Cast<ViewCatalogEntry>();
		col_names = view.names;
		col_types = view.types;
	} else {
		throw InvalidInputException("'%s' is not a table or view", table_name);
	}

	// Validate has at least 2 columns
	if (col_types.size() < 2) {
		throw InvalidInputException("Reference table '%s' must have at least 2 columns (name, length)", table_name);
	}

	// First column should be VARCHAR (reference name)
	if (col_types[0].id() != LogicalTypeId::VARCHAR) {
		throw InvalidInputException("First column of reference table must be VARCHAR (reference name)");
	}

	// Second column should be integer type (length)
	auto length_type = col_types[1].id();
	if (length_type != LogicalTypeId::BIGINT && length_type != LogicalTypeId::INTEGER &&
	    length_type != LogicalTypeId::UBIGINT && length_type != LogicalTypeId::UINTEGER) {
		throw InvalidInputException("Second column of reference table must be an integer type (reference length)");
	}
}

std::unordered_map<std::string, uint64_t> ReadReferenceTable(ClientContext &context, const std::string &table_name) {
	std::unordered_map<std::string, uint64_t> result;

	// Validate schema first
	ValidateReferenceTableSchema(context, table_name);

	// Create a new connection to avoid deadlocking the current context.
	// This is necessary because context.Query() requires locking the context,
	// but during bind/finalize the context is already locked.
	// Using a separate connection allows us to execute queries for both tables and views.
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	// Execute a query to read from the table/view
	// This approach works for both tables and views uniformly
	// We select the first two columns by position (using *)
	std::string query = "SELECT * FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);

	auto query_result = conn.Query(query);

	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from reference table '%s': %s", table_name,
		                            query_result->GetError());
	}

	auto &materialized = query_result->Cast<MaterializedQueryResult>();
	auto &result_types = materialized.types;

	// Get the length column type for proper extraction
	auto &length_type = result_types[1];

	idx_t row_number = 0;

	// Process all chunks from the result
	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}

		auto &name_vector = chunk->data[0];
		auto &length_vector = chunk->data[1];

		// Convert to unified format for proper NULL handling
		UnifiedVectorFormat name_data, length_data;
		name_vector.ToUnifiedFormat(chunk->size(), name_data);
		length_vector.ToUnifiedFormat(chunk->size(), length_data);

		auto name_ptr = UnifiedVectorFormat::GetData<string_t>(name_data);

		for (idx_t i = 0; i < chunk->size(); i++) {
			row_number++;

			auto name_idx = name_data.sel->get_index(i);
			auto length_idx = length_data.sel->get_index(i);

			// Check for NULL values
			if (!name_data.validity.RowIsValid(name_idx)) {
				throw InvalidInputException("NULL reference name at row %llu in table '%s'", row_number, table_name);
			}
			if (!length_data.validity.RowIsValid(length_idx)) {
				throw InvalidInputException("NULL reference length at row %llu in table '%s'", row_number, table_name);
			}

			auto name = name_ptr[name_idx].GetString();

			// Get length value - handle different integer types
			int64_t length;
			switch (length_type.id()) {
			case LogicalTypeId::BIGINT:
				length = UnifiedVectorFormat::GetData<int64_t>(length_data)[length_idx];
				break;
			case LogicalTypeId::INTEGER:
				length = UnifiedVectorFormat::GetData<int32_t>(length_data)[length_idx];
				break;
			case LogicalTypeId::UBIGINT: {
				uint64_t uval = UnifiedVectorFormat::GetData<uint64_t>(length_data)[length_idx];
				if (uval > static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
					throw InvalidInputException("Reference length %llu at row %llu exceeds INT64_MAX in table '%s'",
					                            uval, row_number, table_name);
				}
				length = static_cast<int64_t>(uval);
				break;
			}
			case LogicalTypeId::UINTEGER:
				length = UnifiedVectorFormat::GetData<uint32_t>(length_data)[length_idx];
				break;
			default:
				throw InvalidInputException("Unsupported integer type for reference length");
			}

			if (length < 0) {
				throw InvalidInputException("Reference length must be non-negative for '%s' at row %llu", name,
				                            row_number);
			}

			result.emplace(name, static_cast<uint64_t>(length));
		}
	}

	if (result.empty()) {
		throw InvalidInputException("Reference table '%s' is empty", table_name);
	}

	return result;
}

} // namespace duckdb
