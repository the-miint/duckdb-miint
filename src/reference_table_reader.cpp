#include "reference_table_reader.hpp"
#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include <algorithm>
#include <limits>
#include <unordered_set>

namespace duckdb {

// Check the types the read returned, not the catalog's (docs/internals/reading-tables-views.md)
static void CheckReferenceColumnTypes(const vector<LogicalType> &types, const std::string &table_name) {
	if (types.size() < 2) {
		throw InvalidInputException("Reference table '%s' must have at least 2 columns (name, length)", table_name);
	}

	// A reference name is an id, so it takes the id types
	if (!IsAllowedIdType(types[0])) {
		throw InvalidInputException("First column of reference table '%s' must be %s (reference name)", table_name,
		                            AllowedIdTypeList());
	}

	auto length_type = types[1].id();
	if (length_type != LogicalTypeId::BIGINT && length_type != LogicalTypeId::INTEGER &&
	    length_type != LogicalTypeId::UBIGINT && length_type != LogicalTypeId::UINTEGER) {
		throw InvalidInputException("Second column of reference table '%s' must be an integer type (reference length)",
		                            table_name);
	}
}

// Shared core: read the table/view into (name, length) pairs in the order the scan produced
// them, collapsing duplicate names to their FIRST occurrence. Both public entry points are thin
// wrappers over this so the validation, error messages and duplicate handling live in one place.
// The names are returned as their text; `name_type`, when given, receives the name column's type.
static std::vector<std::pair<std::string, uint64_t>>
ReadReferenceRows(ClientContext &context, const std::string &table_name, optional_ptr<LogicalType> name_type) {
	std::vector<std::pair<std::string, uint64_t>> result;
	std::unordered_set<std::string> seen;

	// Create a new connection to avoid deadlocking the current context.
	// This is necessary because context.Query() requires locking the context,
	// but during bind/finalize the context is already locked.
	// Using a separate connection allows us to execute queries for both tables and views.
	auto conn = MakeReadOnlyHelperConnection(context);

	// Execute a query to read from the table/view
	// This approach works for both tables and views uniformly
	// We select the first two columns by position (using *)
	std::string query = "SELECT * FROM " + SQLIdentifier::ToString(table_name);

	auto query_result = conn.Query(query);

	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from reference table '%s': %s", table_name,
		                            query_result->GetError());
	}

	auto &materialized = *query_result;
	auto &result_types = materialized.GetTypes();
	CheckReferenceColumnTypes(result_types, table_name);
	if (name_type) {
		*name_type = result_types[0];
	}
	std::vector<std::string> names;
	std::vector<bool> name_nulls;

	// Get the length column type for proper extraction
	auto &length_type = result_types[1];

	idx_t row_number = 0;

	// Process all chunks from the result
	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}

		auto &length_vector = chunk->data[1];

		// Convert to unified format for proper NULL handling
		UnifiedVectorFormat length_data;
		length_vector.ToUnifiedFormat(length_data);

		ExtractIdColumnAsStrings(*chunk, 0, result_types[0], names, name_nulls);

		for (idx_t i = 0; i < chunk->size(); i++) {
			row_number++;

			auto length_idx = length_data.sel->get_index(i);

			// Check for NULL values
			if (name_nulls[i]) {
				throw InvalidInputException("NULL reference name at row %llu in table '%s'", row_number, table_name);
			}
			if (!length_data.validity.RowIsValid(length_idx)) {
				throw InvalidInputException("NULL reference length at row %llu in table '%s'", row_number, table_name);
			}

			auto &name = names[i];

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
	auto rows = ReadReferenceRows(context, table_name, nullptr);

	std::unordered_map<std::string, uint64_t> result;
	result.reserve(rows.size());
	for (auto &row : rows) {
		result.emplace(row.first, row.second);
	}
	return result;
}

std::vector<std::pair<std::string, uint64_t>> ReadReferenceTableSortedByName(ClientContext &context,
                                                                             const std::string &table_name) {
	LogicalType name_type;
	auto rows = ReadReferenceRows(context, table_name, name_type);

	// The names' own order: see the header for why name order is the contract rather than the
	// table's row order, and why this must agree with DuckDB's comparison of the name column.
	if (name_type.id() == LogicalTypeId::BIGINT) {
		// Numeric: as text, '10' would sort before '2'
		std::vector<std::pair<int64_t, std::pair<std::string, uint64_t>>> keyed;
		keyed.reserve(rows.size());
		for (auto &row : rows) {
			const auto key = *::miint::ParseIdAsInt64(row.first);
			keyed.emplace_back(key, std::move(row));
		}
		std::sort(keyed.begin(), keyed.end(),
		          [](const std::pair<int64_t, std::pair<std::string, uint64_t>> &a,
		             const std::pair<int64_t, std::pair<std::string, uint64_t>> &b) { return a.first < b.first; });
		for (idx_t i = 0; i < rows.size(); i++) {
			rows[i] = std::move(keyed[i].second);
		}
		return rows;
	}
	// Byte-wise for VARCHAR, and for UUID, which DuckDB orders as its text
	std::sort(rows.begin(), rows.end(),
	          [](const std::pair<std::string, uint64_t> &a, const std::pair<std::string, uint64_t> &b) {
		          return a.first < b.first;
	          });
	return rows;
}

} // namespace duckdb
