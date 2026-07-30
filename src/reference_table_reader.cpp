#include "reference_table_reader.hpp"
#include "catalog_utils.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include <algorithm>
#include <limits>
#include <unordered_set>

namespace duckdb {

static void ValidateReferenceTableSchema(ClientContext &context, const std::string &table_name) {
	auto info = GetTableOrViewColumns(context, table_name, "Reference table");
	auto &col_types = info.types;

	if (col_types.size() < 2) {
		throw InvalidInputException("Reference table '%s' must have at least 2 columns (name, length)", table_name);
	}

	if (col_types[0].id() != LogicalTypeId::VARCHAR) {
		throw InvalidInputException("First column of reference table must be VARCHAR (reference name)");
	}

	auto length_type = col_types[1].id();
	if (length_type != LogicalTypeId::BIGINT && length_type != LogicalTypeId::INTEGER &&
	    length_type != LogicalTypeId::UBIGINT && length_type != LogicalTypeId::UINTEGER) {
		throw InvalidInputException("Second column of reference table must be an integer type (reference length)");
	}
}

// Shared core: read the table/view into (name, length) pairs in the order the scan produced
// them, collapsing duplicate names to their FIRST occurrence. Both public entry points are thin
// wrappers over this so the validation, error messages and duplicate handling live in one place.
static std::vector<std::pair<std::string, uint64_t>> ReadReferenceRows(ClientContext &context,
                                                                       const std::string &table_name) {
	std::vector<std::pair<std::string, uint64_t>> result;
	std::unordered_set<std::string> seen;

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

			if (seen.insert(name).second) {
				result.emplace_back(name, static_cast<uint64_t>(length));
			}
		}
	}

	if (result.empty()) {
		throw InvalidInputException("Reference table '%s' is empty", table_name);
	}

	return result;
}

std::unordered_map<std::string, uint64_t> ReadReferenceTable(ClientContext &context, const std::string &table_name) {
	auto rows = ReadReferenceRows(context, table_name);

	std::unordered_map<std::string, uint64_t> result;
	result.reserve(rows.size());
	for (auto &row : rows) {
		result.emplace(row.first, row.second);
	}
	return result;
}

std::vector<std::pair<std::string, uint64_t>> ReadReferenceTableSortedByName(ClientContext &context,
                                                                             const std::string &table_name) {
	auto rows = ReadReferenceRows(context, table_name);

	// Byte-wise name order. See the header for why name order is the contract rather than the
	// table's row order, and why this must agree with DuckDB's VARCHAR comparison.
	std::sort(rows.begin(), rows.end(),
	          [](const std::pair<std::string, uint64_t> &a, const std::pair<std::string, uint64_t> &b) {
		          return a.first < b.first;
	          });
	return rows;
}

} // namespace duckdb
