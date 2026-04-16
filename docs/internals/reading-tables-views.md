# Reading from Tables and Views in Extensions

When extension code needs to read data from a user-specified table or view (e.g., `PLACEMENTS` parameter in `COPY FORMAT NEWICK`, or `REFERENCE_LENGTHS` in `COPY FORMAT SAM`), use the pattern documented here.

## The Problem
- `Catalog::GetEntry<TableCatalogEntry>` only works for tables, not views
- Views are stored query definitions without physical storage
- `context.Query()` causes deadlocks when called during bind or execution (the context is already locked)

## The Solution: Use a Separate Connection

```cpp
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

// Create a new connection to avoid deadlocking the current context
auto &db = DatabaseInstance::GetDatabase(context);
Connection conn(db);

// Execute a query - works for both tables and views
std::string query = "SELECT col1, col2 FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);
auto result = conn.Query(query);

if (result->HasError()) {
    throw InvalidInputException("Failed to read: %s", result->GetError());
}

// Process the MaterializedQueryResult
auto &materialized = result->Cast<MaterializedQueryResult>();
while (auto chunk = materialized.Fetch()) {
    // Process chunk->data[0], chunk->data[1], etc.
}
```

## Schema Validation for Tables/Views

```cpp
// Use TABLE_ENTRY lookup which returns either tables or views
EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY, table_name, QueryErrorContext());
auto entry = Catalog::GetEntry(context, INVALID_CATALOG, INVALID_SCHEMA, lookup_info, OnEntryNotFound::RETURN_NULL);

if (!entry) {
    throw BinderException("'%s' does not exist", table_name);
}

// Check entry type and get columns
if (entry->type == CatalogType::TABLE_ENTRY) {
    auto &table = entry->Cast<TableCatalogEntry>();
    // Use table.GetColumns()
} else if (entry->type == CatalogType::VIEW_ENTRY) {
    auto &view = entry->Cast<ViewCatalogEntry>();
    // Use view.names and view.types
}
```

## Important Notes
- The separate connection runs in its own transaction, so uncommitted changes from the original context may not be visible
- For most use cases (reading from persistent tables/views), this is fine
- Read placements/references at bind time if possible, storing results in bind data to avoid issues during finalize
- See `src/placement_table_reader.cpp` for a complete example
