#include "placement_table_reader.hpp"
#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include <limits>

namespace duckdb {

// The placement columns, in the order ReadPlacementTable selects them, and the types each may
// have. The bind-time check and the read both check them, as docs/internals/reading-tables-views.md
// says to.
struct PlacementColumn {
	const char *name;
	bool (*accepts)(const LogicalType &type);
	const char *description;
};

static bool IsPlacementEdgeIdType(const LogicalType &type) {
	switch (type.id()) {
	case LogicalTypeId::BIGINT:
	case LogicalTypeId::INTEGER:
	case LogicalTypeId::UBIGINT:
	case LogicalTypeId::UINTEGER:
		return true;
	default:
		return false;
	}
}

static bool IsPlacementLengthType(const LogicalType &type) {
	return type.id() == LogicalTypeId::DOUBLE || type.id() == LogicalTypeId::FLOAT;
}

static const PlacementColumn PLACEMENT_COLUMNS[] = {
    // A fragment is a read, so its id takes the id types
    {"fragment_id", IsAllowedIdType, AllowedIdTypeList()},
    {"edge_id", IsPlacementEdgeIdType, "an integer type"},
    {"like_weight_ratio", IsPlacementLengthType, "DOUBLE or FLOAT"},
    {"distal_length", IsPlacementLengthType, "DOUBLE or FLOAT"},
    {"pendant_length", IsPlacementLengthType, "DOUBLE or FLOAT"},
};

void ValidatePlacementTableSchema(ClientContext &context, const std::string &table_name) {
	auto info = GetTableOrViewColumns(context, table_name, "Placement table");
	auto &col_names = info.names;
	auto &col_types = info.types;

	// Build name-to-index map (case-insensitive)
	std::unordered_map<string, idx_t> name_to_idx;
	for (idx_t i = 0; i < col_names.size(); i++) {
		name_to_idx[StringUtil::Lower(col_names[i])] = i;
	}

	// Check required columns exist and have correct types
	for (auto &column : PLACEMENT_COLUMNS) {
		auto it = name_to_idx.find(column.name);
		if (it == name_to_idx.end()) {
			throw BinderException("Placement table '%s' missing required column '%s'", table_name, column.name);
		}
		if (!column.accepts(col_types[it->second])) {
			throw BinderException("Column '%s' in table '%s' must be %s", column.name, table_name, column.description);
		}
	}
}

std::vector<miint::Placement> ReadPlacementTable(ClientContext &context, const std::string &table_name) {
	std::vector<miint::Placement> result;

	// A helper connection, as the caller's context is running its own query
	// (docs/internals/reading-tables-views.md)
	auto conn = MakeReadOnlyHelperConnection(context);

	// Execute a query to read from the table/view
	// This approach works for both tables and views uniformly
	// Pin the row order. NewickTree::insert_fully_resolved deduplicates by
	// fragment_id in a single greedy pass whose comparator is epsilon-tolerant,
	// and therefore not transitive: a ~ b and b ~ c do not imply a ~ c, so which
	// row survives depends on the order they arrive in. Taking whatever the scan
	// yields made the resolved tree vary between runs of the same query.
	//
	// This pins the survivor; it does not make it the highest like_weight_ratio.
	// A chain of rows each within epsilon of the next walks the greedy pass down
	// to the lowest pendant_length in the chain, which can sit outside epsilon of
	// the maximum. test/sql/tree_resolve_placement_order.test carries that case.
	std::string query = "SELECT fragment_id, edge_id, like_weight_ratio, distal_length, pendant_length FROM " +
	                    SQLIdentifier::ToString(table_name) +
	                    " ORDER BY fragment_id, like_weight_ratio DESC, pendant_length, edge_id";

	auto query_result = conn.Query(query);

	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from placement table '%s': %s", table_name,
		                            query_result->GetError());
	}

	auto &materialized = *query_result;
	auto &result_types = materialized.GetTypes();
	D_ASSERT(result_types.size() == sizeof(PLACEMENT_COLUMNS) / sizeof(PLACEMENT_COLUMNS[0]));
	for (idx_t i = 0; i < result_types.size(); i++) {
		auto &column = PLACEMENT_COLUMNS[i];
		if (!column.accepts(result_types[i])) {
			throw InvalidInputException("Column '%s' in table '%s' must be %s", column.name, table_name,
			                            column.description);
		}
	}

	// Get column types for proper data extraction
	// Column order: fragment_id(0), edge_id(1), like_weight_ratio(2), distal_length(3), pendant_length(4)
	auto &fragment_id_type = result_types[0];
	auto &edge_id_type = result_types[1];
	auto &lwr_type = result_types[2];
	auto &distal_type = result_types[3];
	auto &pendant_type = result_types[4];

	idx_t row_number = 0;
	std::vector<std::string> fragment_ids;
	std::vector<bool> fragment_nulls;

	// Process all chunks from the result
	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}

		auto &edge_id_vec = chunk->data[1];
		auto &lwr_vec = chunk->data[2];
		auto &distal_vec = chunk->data[3];
		auto &pendant_vec = chunk->data[4];

		// Convert to unified format for proper NULL handling
		UnifiedVectorFormat edge_data, lwr_data, distal_data, pendant_data;
		edge_id_vec.ToUnifiedFormat(edge_data);
		lwr_vec.ToUnifiedFormat(lwr_data);
		distal_vec.ToUnifiedFormat(distal_data);
		pendant_vec.ToUnifiedFormat(pendant_data);

		ExtractIdColumnAsStrings(*chunk, 0, fragment_id_type, fragment_ids, fragment_nulls);

		for (idx_t i = 0; i < chunk->size(); i++) {
			row_number++;

			auto edge_idx = edge_data.sel->get_index(i);
			auto lwr_idx = lwr_data.sel->get_index(i);
			auto distal_idx = distal_data.sel->get_index(i);
			auto pendant_idx = pendant_data.sel->get_index(i);

			// Skip rows with NULL fragment_id
			if (fragment_nulls[i]) {
				continue;
			}

			// Check for NULL values in required columns
			if (!edge_data.validity.RowIsValid(edge_idx)) {
				throw InvalidInputException("NULL edge_id at row %llu in placement table '%s'", row_number, table_name);
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
			p.fragment_id = std::move(fragment_ids[i]);

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
