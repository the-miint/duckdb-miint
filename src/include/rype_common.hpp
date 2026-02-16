#pragma once

#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/main/client_context.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace duckdb {

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
		for (const auto &name : view.names) {
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
