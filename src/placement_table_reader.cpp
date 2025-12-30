#include "placement_table_reader.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/storage/data_table.hpp"
#include "duckdb/storage/table/scan_state.hpp"
#include "duckdb/transaction/duck_transaction.hpp"
#include <limits>

namespace duckdb {

void ValidatePlacementTableSchema(ClientContext &context, const std::string &table_name) {
	// Validate table exists at bind time
	auto catalog_entry = Catalog::GetEntry<TableCatalogEntry>(context, INVALID_CATALOG, INVALID_SCHEMA, table_name,
	                                                          OnEntryNotFound::RETURN_NULL);
	if (!catalog_entry) {
		throw BinderException("Placement table '%s' does not exist", table_name);
	}

	auto &columns = catalog_entry->GetColumns();

	// Build name-to-index map
	std::unordered_map<string, idx_t> name_to_idx;
	for (idx_t i = 0; i < columns.LogicalColumnCount(); i++) {
		auto &col = columns.GetColumn(LogicalIndex(i));
		name_to_idx[StringUtil::Lower(col.Name())] = i;
	}

	// Check required columns exist and have correct types
	auto check_column = [&](const string &col_name, const vector<LogicalTypeId> &allowed_types,
	                        const string &type_desc) {
		auto it = name_to_idx.find(col_name);
		if (it == name_to_idx.end()) {
			throw BinderException("Placement table '%s' missing required column '%s'", table_name, col_name);
		}
		auto &col_type = columns.GetColumn(LogicalIndex(it->second)).Type();
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

	// Get the table catalog entry
	auto catalog_entry = Catalog::GetEntry<TableCatalogEntry>(context, INVALID_CATALOG, INVALID_SCHEMA, table_name,
	                                                          OnEntryNotFound::RETURN_NULL);
	if (!catalog_entry) {
		throw InvalidInputException("Placement table '%s' does not exist", table_name);
	}

	auto &storage = catalog_entry->GetStorage();
	auto &columns = catalog_entry->GetColumns();

	// Find column indices by name
	idx_t fragment_id_col = DConstants::INVALID_INDEX;
	idx_t edge_id_col = DConstants::INVALID_INDEX;
	idx_t lwr_col = DConstants::INVALID_INDEX;
	idx_t distal_col = DConstants::INVALID_INDEX;
	idx_t pendant_col = DConstants::INVALID_INDEX;

	for (idx_t i = 0; i < columns.LogicalColumnCount(); i++) {
		auto &col = columns.GetColumn(LogicalIndex(i));
		auto lower_name = StringUtil::Lower(col.Name());
		if (lower_name == "fragment_id") {
			fragment_id_col = i;
		} else if (lower_name == "edge_id") {
			edge_id_col = i;
		} else if (lower_name == "like_weight_ratio") {
			lwr_col = i;
		} else if (lower_name == "distal_length") {
			distal_col = i;
		} else if (lower_name == "pendant_length") {
			pendant_col = i;
		}
	}

	// Validate all required columns found
	if (fragment_id_col == DConstants::INVALID_INDEX) {
		throw InvalidInputException("Placement table '%s' missing required column 'fragment_id'", table_name);
	}
	if (edge_id_col == DConstants::INVALID_INDEX) {
		throw InvalidInputException("Placement table '%s' missing required column 'edge_id'", table_name);
	}
	if (lwr_col == DConstants::INVALID_INDEX) {
		throw InvalidInputException("Placement table '%s' missing required column 'like_weight_ratio'", table_name);
	}
	if (distal_col == DConstants::INVALID_INDEX) {
		throw InvalidInputException("Placement table '%s' missing required column 'distal_length'", table_name);
	}
	if (pendant_col == DConstants::INVALID_INDEX) {
		throw InvalidInputException("Placement table '%s' missing required column 'pendant_length'", table_name);
	}

	// Get column types
	auto &fragment_id_type = columns.GetColumn(LogicalIndex(fragment_id_col)).Type();
	auto &edge_id_type = columns.GetColumn(LogicalIndex(edge_id_col)).Type();
	auto &lwr_type = columns.GetColumn(LogicalIndex(lwr_col)).Type();
	auto &distal_type = columns.GetColumn(LogicalIndex(distal_col)).Type();
	auto &pendant_type = columns.GetColumn(LogicalIndex(pendant_col)).Type();

	// Set up scan for required columns
	vector<StorageIndex> column_ids = {StorageIndex(fragment_id_col), StorageIndex(edge_id_col),
	                                   StorageIndex(lwr_col), StorageIndex(distal_col), StorageIndex(pendant_col)};

	auto &transaction = DuckTransaction::Get(context, catalog_entry->catalog);
	TableScanState scan_state;
	scan_state.Initialize(column_ids, &context);
	storage.InitializeScan(context, transaction, scan_state, column_ids);

	// Set up chunk for results
	DataChunk chunk;
	vector<LogicalType> chunk_types = {fragment_id_type, edge_id_type, lwr_type, distal_type, pendant_type};
	chunk.Initialize(Allocator::Get(context), chunk_types);

	idx_t row_number = 0;

	while (true) {
		chunk.Reset();
		storage.Scan(transaction, chunk, scan_state);

		if (chunk.size() == 0) {
			break;
		}

		auto &fragment_id_vec = chunk.data[0];
		auto &edge_id_vec = chunk.data[1];
		auto &lwr_vec = chunk.data[2];
		auto &distal_vec = chunk.data[3];
		auto &pendant_vec = chunk.data[4];

		auto fragment_ids = FlatVector::GetData<string_t>(fragment_id_vec);
		auto &fragment_validity = FlatVector::Validity(fragment_id_vec);
		auto &edge_validity = FlatVector::Validity(edge_id_vec);
		auto &lwr_validity = FlatVector::Validity(lwr_vec);
		auto &distal_validity = FlatVector::Validity(distal_vec);
		auto &pendant_validity = FlatVector::Validity(pendant_vec);

		for (idx_t i = 0; i < chunk.size(); i++) {
			row_number++;

			// Skip rows with NULL fragment_id
			if (!fragment_validity.RowIsValid(i)) {
				continue;
			}

			// Check for NULL values in required columns
			if (!edge_validity.RowIsValid(i)) {
				throw InvalidInputException("NULL edge_id at row %llu in placement table '%s'", row_number,
				                            table_name);
			}
			if (!lwr_validity.RowIsValid(i)) {
				throw InvalidInputException("NULL like_weight_ratio at row %llu in placement table '%s'", row_number,
				                            table_name);
			}
			if (!distal_validity.RowIsValid(i)) {
				throw InvalidInputException("NULL distal_length at row %llu in placement table '%s'", row_number,
				                            table_name);
			}
			if (!pendant_validity.RowIsValid(i)) {
				throw InvalidInputException("NULL pendant_length at row %llu in placement table '%s'", row_number,
				                            table_name);
			}

			miint::Placement p;
			p.fragment_id = fragment_ids[i].GetString();

			// Get edge_id (handle different integer types)
			switch (edge_id_type.id()) {
			case LogicalTypeId::BIGINT:
				p.edge_id = FlatVector::GetData<int64_t>(edge_id_vec)[i];
				break;
			case LogicalTypeId::INTEGER:
				p.edge_id = FlatVector::GetData<int32_t>(edge_id_vec)[i];
				break;
			case LogicalTypeId::UBIGINT: {
				uint64_t uval = FlatVector::GetData<uint64_t>(edge_id_vec)[i];
				if (uval > static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
					throw InvalidInputException("edge_id %llu at row %llu exceeds INT64_MAX in placement table '%s'",
					                            uval, row_number, table_name);
				}
				p.edge_id = static_cast<int64_t>(uval);
				break;
			}
			case LogicalTypeId::UINTEGER:
				p.edge_id = FlatVector::GetData<uint32_t>(edge_id_vec)[i];
				break;
			default:
				throw InvalidInputException("Unsupported integer type for edge_id");
			}

			// Get like_weight_ratio
			if (lwr_type.id() == LogicalTypeId::DOUBLE) {
				p.like_weight_ratio = FlatVector::GetData<double>(lwr_vec)[i];
			} else {
				p.like_weight_ratio = FlatVector::GetData<float>(lwr_vec)[i];
			}

			// Get distal_length
			if (distal_type.id() == LogicalTypeId::DOUBLE) {
				p.distal_length = FlatVector::GetData<double>(distal_vec)[i];
			} else {
				p.distal_length = FlatVector::GetData<float>(distal_vec)[i];
			}

			// Get pendant_length
			if (pendant_type.id() == LogicalTypeId::DOUBLE) {
				p.pendant_length = FlatVector::GetData<double>(pendant_vec)[i];
			} else {
				p.pendant_length = FlatVector::GetData<float>(pendant_vec)[i];
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
