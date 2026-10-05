#include "catalog_utils.hpp"
#include "miint_log.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_data.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/parser/keyword_helper.hpp"
#include "duckdb/parser/qualified_name.hpp"
#include "duckdb/planner/binder.hpp"

namespace duckdb {

TableOrViewColumns GetTableOrViewColumns(ClientContext &context, const std::string &table_name,
                                         const std::string &entity_type) {
	EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY, QualifiedName(Identifier(table_name)), QueryErrorContext());
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
			result.names.push_back(col.Name().GetIdentifierName());
			result.types.push_back(col.Type());
		}
		result.is_physical_table = true;
	} else if (entry->type == CatalogType::VIEW_ENTRY) {
		auto &view = entry->Cast<ViewCatalogEntry>();
		view.BindView(context);
		auto col_info = view.GetColumnInfo();
		result.names = IdentifiersToStrings(col_info->names);
		result.types = col_info->types;
		result.is_physical_table = false;
	} else {
		throw BinderException("'%s' is not a table or view", table_name);
	}

	return result;
}

void RejectRelationNameAsLiteral(ClientContext &context, const std::string &function_name, const std::string &literal) {
	// Unqualified lookup first. This honours the session search_path, so a bare name
	// naming a table in a non-default schema is caught.
	EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY, QualifiedName(Identifier(literal)), QueryErrorContext());
	auto entry = Catalog::GetEntry(context, INVALID_CATALOG, INVALID_SCHEMA, lookup_info, OnEntryNotFound::RETURN_NULL);

	if (!entry && literal.find('.') != std::string::npos) {
		// A qualified name ("s.accs", "db.s.accs") is the same mistake and must not slip
		// through. Parse rather than split by hand — QualifiedName handles quoting and
		// the catalog.schema.name form.
		//
		// Note accessions legitimately contain dots ("NC_001416.1"), so this branch is
		// routinely entered for valid input. That is harmless: it resolves schema
		// "NC_001416" / name "1", which does not exist, and we fall through. The guard
		// still keys off catalog residency, never off the string's shape.
		auto qname = QualifiedName::Parse(literal);
		// v2.0: the QualifiedName carries catalog/schema; no separate catalog/schema arguments
		EntryLookupInfo qualified_lookup(CatalogType::TABLE_ENTRY, qname, QueryErrorContext());
		entry = Catalog::GetEntry(context, qualified_lookup, OnEntryNotFound::RETURN_NULL);
	}

	if (!entry) {
		return;
	}

	// TABLE_ENTRY lookup returns tables and views; name whichever it actually is, so
	// the message matches what the user sees in their own schema.
	const char *kind = entry->type == CatalogType::VIEW_ENTRY ? "view" : "table";
	throw InvalidInputException(
	    "%s: '%s' is a %s in the catalog, not an NCBI accession. %s takes literal accessions, e.g. "
	    "%s('NC_001416.1') or %s(['NC_001416.1', 'NC_000913.3']). To use accessions stored in '%s', hoist "
	    "them into a list first — a subquery cannot be a table-function argument: "
	    "SET VARIABLE accs = (SELECT list(<column>) FROM %s); %s(getvariable('accs'));",
	    function_name, literal, kind, function_name, function_name, function_name, literal, literal, function_name);
}

void RejectCTERelationName(TableFunctionBindInput &input, const std::string &relation_name) {
	// Null when DuckDB re-binds a deserialized plan (LogicalGet::Deserialize); no
	// query text, so no CTE, is in scope there.
	if (!input.binder) {
		return;
	}
	// Mirrors Binder::Bind(BaseTableRef), so the name gets the scoping a plain
	// FROM-clause reference would. GetCTEBinding walks parent binders only through
	// REGULAR_BINDERs, which keeps a caller's CTE out of a view body. Inside a CTE's
	// own body (or a recursive CTE's anchor) the lookup returns a CANNOT_BE_REFERENCED
	// self-binding, and there the name still means the catalog relation.
	auto cte = input.binder->GetCTEBinding(BindingAlias(Identifier(relation_name)));
	if (!cte || !cte->CanBeReferenced()) {
		return;
	}
	// The suggested view gets a NEW name: in the shadowing case the CTE body usually
	// reads the same-named table, so reusing the name would make a view that
	// references itself. TEMP is only the common case — massql and some sample_id
	// paths still cannot read TEMP relations (#207).
	auto suggested = KeywordHelper::WriteOptionallyQuoted(relation_name + "_view");
	throw BinderException("%s: '%s' is a common table expression (WITH clause). Relation names are resolved in the "
	                      "catalog, where CTEs do not exist, so it cannot be read by name. Create a view or table "
	                      "under a different name and pass that name instead, e.g. CREATE TEMP VIEW %s AS ...; "
	                      "(massql and some sample_id modes cannot read TEMP relations yet, #207: use CREATE VIEW or "
	                      "CREATE TABLE there).",
	                      input.table_function.GetName().GetIdentifierName(), relation_name, suggested);
}

void InheritTempObjects(ClientContext &context, Connection &conn) {
	// The connection's own temporary_objects — freshly constructed and empty — is
	// released here. Nothing outside ClientData references it: ClientData's
	// constructor only borrows an oid from the DatabaseManager, and lookups read the
	// pointer back out of ClientData dynamically (database_manager.cpp:61), so there
	// is no registry entry or cached oid to invalidate.
	ClientData::Get(*conn.context).temporary_objects = ClientData::Get(context).temporary_objects;
}

Connection MakeReadOnlyHelperConnection(ClientContext &context) {
	Connection conn(DatabaseInstance::GetDatabase(context));
	InheritTempObjects(context, conn);
	return conn;
}

void DropHelperTempRelation(Connection &conn, const std::string &quoted_name) {
	if (quoted_name.empty()) {
		return;
	}
	// Connection::Query converts everything DuckDB throws into an error on the result
	// (ClientContext::Query catches std::exception plus a catch-all), and EmitWarning
	// is itself guarded — but this runs in destructors, where propagating anything at
	// all would terminate, so the whole body is belt-and-braces wrapped.
	try {
		auto result = conn.Query("DROP TABLE IF EXISTS " + quoted_name);
		if (result->HasError()) {
			miint::EmitWarning(*conn.context,
			                   "Failed to drop internal temporary relation %s: %s. It may remain "
			                   "visible in this session until it ends.",
			                   quoted_name, result->GetError());
		}
	} catch (...) { // NOLINT: a destructor must not propagate; see above
	}
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
