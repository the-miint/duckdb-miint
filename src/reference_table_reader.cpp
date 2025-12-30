#include "reference_table_reader.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/storage/data_table.hpp"
#include "duckdb/storage/table/scan_state.hpp"
#include "duckdb/transaction/duck_transaction.hpp"
#include <limits>

namespace duckdb {

std::unordered_map<std::string, uint64_t> ReadReferenceTable(ClientContext &context, const std::string &table_name) {
	std::unordered_map<std::string, uint64_t> result;

	// Get the table catalog entry
	auto catalog_entry = Catalog::GetEntry<TableCatalogEntry>(context, INVALID_CATALOG, INVALID_SCHEMA, table_name,
	                                                          OnEntryNotFound::RETURN_NULL);
	if (!catalog_entry) {
		throw InvalidInputException("Table '%s' does not exist", table_name);
	}

	// Get storage
	auto &storage = catalog_entry->GetStorage();
	auto &columns = catalog_entry->GetColumns();

	// Validate table has at least 2 columns
	if (columns.LogicalColumnCount() < 2) {
		throw InvalidInputException("Reference table '%s' must have at least 2 columns (name, length)", table_name);
	}

	// First column should be VARCHAR (reference name), second should be integer type (length)
	auto &name_type = columns.GetColumn(LogicalIndex(0)).Type();
	auto &length_type = columns.GetColumn(LogicalIndex(1)).Type();

	if (name_type.id() != LogicalTypeId::VARCHAR) {
		throw InvalidInputException("First column of reference table must be VARCHAR (reference name)");
	}

	if (length_type.id() != LogicalTypeId::BIGINT && length_type.id() != LogicalTypeId::INTEGER &&
	    length_type.id() != LogicalTypeId::UBIGINT && length_type.id() != LogicalTypeId::UINTEGER) {
		throw InvalidInputException("Second column of reference table must be an integer type (reference length)");
	}

	// Set up scan (scan only first two columns)
	vector<StorageIndex> column_ids = {StorageIndex(0), StorageIndex(1)};

	auto &transaction = DuckTransaction::Get(context, catalog_entry->catalog);
	TableScanState scan_state;
	scan_state.Initialize(column_ids, &context);
	storage.InitializeScan(context, transaction, scan_state, column_ids);

	// Scan the table
	DataChunk chunk;
	vector<LogicalType> chunk_types = {name_type, length_type};
	chunk.Initialize(Allocator::Get(context), chunk_types);

	idx_t row_number = 0;

	while (true) {
		chunk.Reset();
		storage.Scan(transaction, chunk, scan_state);

		if (chunk.size() == 0) {
			break;
		}

		auto &name_vector = chunk.data[0];
		auto &length_vector = chunk.data[1];

		auto name_data = FlatVector::GetData<string_t>(name_vector);
		auto &name_validity = FlatVector::Validity(name_vector);
		auto &length_validity = FlatVector::Validity(length_vector);

		for (idx_t i = 0; i < chunk.size(); i++) {
			row_number++;

			// Check for NULL values
			if (!name_validity.RowIsValid(i)) {
				throw InvalidInputException("NULL reference name at row %llu in table '%s'", row_number, table_name);
			}
			if (!length_validity.RowIsValid(i)) {
				throw InvalidInputException("NULL reference length at row %llu in table '%s'", row_number, table_name);
			}

			auto name = name_data[i].GetString();

			// Get length value - handle different integer types
			int64_t length;
			switch (length_type.id()) {
			case LogicalTypeId::BIGINT:
				length = FlatVector::GetData<int64_t>(length_vector)[i];
				break;
			case LogicalTypeId::INTEGER:
				length = FlatVector::GetData<int32_t>(length_vector)[i];
				break;
			case LogicalTypeId::UBIGINT: {
				uint64_t uval = FlatVector::GetData<uint64_t>(length_vector)[i];
				if (uval > static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
					throw InvalidInputException(
					    "Reference length %llu at row %llu exceeds INT64_MAX in table '%s'", uval, row_number,
					    table_name);
				}
				length = static_cast<int64_t>(uval);
				break;
			}
			case LogicalTypeId::UINTEGER:
				length = FlatVector::GetData<uint32_t>(length_vector)[i];
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
