#include "catalog_utils.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"

namespace duckdb {

TableOrViewColumns GetTableOrViewColumns(ClientContext &context, const std::string &table_name,
                                         const std::string &entity_type) {
	EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY, table_name, QueryErrorContext());
	auto entry = Catalog::GetEntry(context, INVALID_CATALOG, INVALID_SCHEMA, lookup_info, OnEntryNotFound::RETURN_NULL);

	if (!entry) {
		throw BinderException("%s or view '%s' does not exist", entity_type, table_name);
	}

	TableOrViewColumns result;

	if (entry->type == CatalogType::TABLE_ENTRY) {
		auto &table = entry->Cast<TableCatalogEntry>();
		auto &columns = table.GetColumns();
		for (idx_t i = 0; i < columns.LogicalColumnCount(); i++) {
			auto &col = columns.GetColumn(LogicalIndex(i));
			result.names.push_back(col.Name());
			result.types.push_back(col.Type());
		}
		result.is_physical_table = true;
	} else if (entry->type == CatalogType::VIEW_ENTRY) {
		auto &view = entry->Cast<ViewCatalogEntry>();
		view.BindView(context);
		auto col_info = view.GetColumnInfo();
		result.names = col_info->names;
		result.types = col_info->types;
		result.is_physical_table = false;
	} else {
		throw BinderException("'%s' is not a table or view", table_name);
	}

	return result;
}

bool HasColumn(const TableOrViewColumns &columns, const std::string &col) {
	auto target = StringUtil::Lower(col);
	for (const auto &name : columns.names) {
		if (StringUtil::Lower(name) == target) {
			return true;
		}
	}
	return false;
}

} // namespace duckdb
