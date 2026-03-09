#pragma once

#include "rype.h"

#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/arrow/arrow_wrapper.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/table/arrow.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace duckdb {

// ============================================================================
// Shared Arrow local state for all RYpe table functions.
// Manages per-column ArrowArrayScanState for zero-copy Arrow→DuckDB conversion.
// ============================================================================
struct RypeArrowLocalState : public LocalTableFunctionState {
	std::unordered_map<idx_t, unique_ptr<ArrowArrayScanState>> array_states;
	ClientContext &context;

	explicit RypeArrowLocalState(ClientContext &ctx) : context(ctx) {
	}

	~RypeArrowLocalState() {
		array_states.clear();
	}

	ArrowArrayScanState &GetState(idx_t col_idx) {
		auto it = array_states.find(col_idx);
		if (it == array_states.end()) {
			auto state = make_uniq<ArrowArrayScanState>(context);
			auto &ref = *state;
			array_states.emplace(col_idx, std::move(state));
			return ref;
		}
		return *it->second;
	}

	void ResetStates() {
		for (auto &state : array_states) {
			state.second->Reset();
		}
	}
};

//! Sample average read length from the first 1000 sequences in a table.
//! Returns fallback (default 300) if the query fails or the table is empty.
inline size_t SampleAvgReadLength(Connection &conn, const std::string &table_quoted, size_t fallback = 300) {
	std::string query =
	    "SELECT AVG(LENGTH(sequence1))::BIGINT FROM (SELECT sequence1 FROM " + table_quoted + " LIMIT 1000)";
	auto result = conn.Query(query);
	if (!result->HasError()) {
		auto &materialized = result->Cast<MaterializedQueryResult>();
		auto chunk = materialized.Fetch();
		if (chunk && chunk->size() > 0 && !chunk->data[0].GetValue(0).IsNull()) {
			auto val = chunk->data[0].GetValue(0).GetValue<int64_t>();
			if (val > 0) {
				return static_cast<size_t>(val);
			}
		}
	}
	return fallback;
}

//! Validate that a table/view exists and has the required columns for RYpe functions.
//! Returns true if the optional "sequence2" column is present (used by rype_classify
//! for paired-end reads). All RYpe functions require id_column and "sequence1".
inline bool ValidateSequenceTable(ClientContext &context, const std::string &table_name, const std::string &id_column) {
	EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY, table_name, QueryErrorContext());
	auto entry = Catalog::GetEntry(context, INVALID_CATALOG, INVALID_SCHEMA, lookup_info, OnEntryNotFound::RETURN_NULL);

	if (!entry) {
		throw BinderException("Table or view '%s' does not exist", table_name);
	}

	vector<string> col_names;
	if (entry->type == CatalogType::TABLE_ENTRY) {
		auto &table = entry->Cast<TableCatalogEntry>();
		auto &columns = table.GetColumns();
		for (idx_t i = 0; i < columns.LogicalColumnCount(); i++) {
			col_names.push_back(StringUtil::Lower(columns.GetColumn(LogicalIndex(i)).Name()));
		}
	} else if (entry->type == CatalogType::VIEW_ENTRY) {
		auto &view = entry->Cast<ViewCatalogEntry>();
		view.BindView(context);
		auto col_info = view.GetColumnInfo();
		for (const auto &name : col_info->names) {
			col_names.push_back(StringUtil::Lower(name));
		}
	} else {
		throw BinderException("'%s' is not a table or view", table_name);
	}

	auto id_col_lower = StringUtil::Lower(id_column);
	bool has_id = std::find(col_names.begin(), col_names.end(), id_col_lower) != col_names.end();
	bool has_seq1 = std::find(col_names.begin(), col_names.end(), "sequence1") != col_names.end();
	bool has_seq2 = std::find(col_names.begin(), col_names.end(), "sequence2") != col_names.end();

	if (!has_id) {
		throw BinderException("Table '%s' missing required column '%s'", table_name, id_column);
	}
	if (!has_seq1) {
		throw BinderException("Table '%s' missing required column 'sequence1'", table_name);
	}

	return has_seq2;
}

} // namespace duckdb
