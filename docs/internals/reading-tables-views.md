# Reading from Tables and Views in Extensions

When extension code needs to read data from a user-specified table or view (e.g., `PLACEMENTS` parameter in `COPY FORMAT NEWICK`, or `REFERENCE_LENGTHS` in `COPY FORMAT SAM`), use the pattern documented here.

## The Problem
- `Catalog::GetEntry<TableCatalogEntry>` only works for tables, not views
- Views are stored query definitions without physical storage
- `context.Query()` causes deadlocks when called during bind or execution (the context is already locked)

## The Solution: Use a Separate Connection

```cpp
#include "catalog_utils.hpp"

// A helper connection that can resolve the caller's TEMP tables and views.
// Replaces the raw `DatabaseInstance::GetDatabase(context)` + `Connection conn(db)`
// pair — see "TEMP tables and views" below for why the raw form is wrong.
auto conn = MakeReadOnlyHelperConnection(context);

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

## TEMP tables and views

A bare `Connection conn(db)` **cannot see the caller's TEMP objects at all.** DuckDB
constructs `ClientData::temporary_objects` per `ClientContext` (`duckdb/src/main/client_data.cpp`),
so a helper connection starts with its own empty TEMP catalog.

The failure mode is nasty rather than obvious. Bind-time validation
(`GetTableOrViewColumns`) runs `Catalog::GetEntry` against the *caller's* context,
which does see TEMP objects, so bind succeeds; only execution fails, with a message
claiming a relation the user can plainly see does not exist. This was issue #193,
and it is worst under a read-only catalog, where a TEMP table is the only relation
a user can create at all.

Use `MakeReadOnlyHelperConnection(context)` (`src/include/catalog_utils.hpp`), which
shares the caller's `temporary_objects` pointer, or `InheritTempObjects(context, conn)`
when the connection has to be constructed separately (e.g. a `unique_ptr` member of a
global/local state).

**The constraint that comes with it:** an inheriting connection must not create TEMP
objects with fixed names. Anything it creates lands in the *caller's* catalog, so

- parallel workers sharing one fixed name collide there — a silent wrong answer,
  strictly worse than the resolution failure being fixed; and
- the object outlives the helper connection instead of dying with it.

A connection that must create a TEMP object may still inherit if it gives the object a
unique name **and** drops it on every exit path. Two forms, both in `catalog_utils.hpp`:

- `HelperTempRelation` — RAII, for an object created and finished with inside one
  function. `ReadBatchByIds` (`src/sequence_table_reader.cpp`) is the reference caller.
  `WaveDistanceBlockSource::Prefetch` (`src/unifrac_pcoa_function.cpp`) is the second,
  and shows two wrinkles: it re-creates its two staging tables once per wave, so the
  guards are scoped to the wave rather than to the whole run; and it arms them *before*
  the `CREATE`s, because two statements are issued and a failure on the second would
  otherwise leak the first. Arming early is safe — the drop is `DROP ... IF EXISTS`.
- `DropHelperTempRelation` — call it from the destructor of whatever state owns the
  object, when it outlives the creating function. `MaterializeRypeInputTempTable` plus
  the `rype_*` global states are the reference callers.

Both report a failed drop through `miint::EmitWarning` rather than discarding it, because
on an inheriting connection a failed drop leaves an internal relation in the user's
session — it no longer disappears with the connection, so silence would hide a real leak.

Sites that fail either half of that condition therefore keep an isolated catalog and
still cannot read a TEMP relation. Today those all fail the *naming* half: the per-sample
(`sample_id := ...`) paths in `woltka_ogu`, `detect_chimera_uchime` and `sylph_profile`
create a fixed `__<tool>_per_sample` view, and `massql` creates fixed `__massql_base` /
`__massql_ms1` / `__massql_per_sample` / `__massql_src` — note that `__massql_base` is
created on *every* `massql()` call, not only the per-sample path, so `massql` cannot take
a TEMP source at all. Uniquifying those names so they can inherit too is tracked as #207.

`QuerySequenceStream`'s `Connection &` overload exists for exactly those callers — they
pass a connection they created a TEMP view on, so it deliberately does not inherit.

## Helper connections at EXECUTION time

Most callers open a helper connection during `Bind`. It also works from `InitGlobal`
(`DeblurInitGlobal`) and from `Execute` — the progressive PCoA functions run their whole
wave-fill (`CREATE TEMPORARY TABLE`, an `Appender`, and three large joins) from inside
`ProgressivePcoaExecute`, once per wave, so a run can emit rows as it goes instead of
computing everything up front. That means the queries run on a task-scheduler worker
thread rather than the client thread, and it is fine for the same reason the recipe works
at all: the helper has its own `ClientContext`, so there is no re-entrant lock, and a
nested query drives its own execution on the calling thread rather than waiting for a free
worker. Verified by `test/sql/temp_table_resolution_unifrac.test`, which reads a TEMP
distance table through that path and then asserts no `_miint_wave_%` staging table
survived.

Two things to keep in mind if you do this:

- The scan must be single-threaded (`MaxThreads() == 1`) if the helper work is stateful,
  as it is here — otherwise several `Execute` calls issue overlapping nested queries.
- Poll `context.interrupted` around the work and throw `InterruptException`. Nothing else
  will: a blocked `Execute` is not a place DuckDB can cancel on its own, so long work
  behind a helper connection is uninterruptible unless the caller checks.

Failing only the *drop* half is cheap to fix rather than document around, and worth
checking for before writing a caveat: `mzml_peak_pair` already suffixed its two temp
tables from a process-wide counter, so it needed a `HelperTempRelation` guard and
nothing else.

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
- The separate connection runs in its own transaction, so uncommitted changes from the original context are **not** visible. A relation created inside an explicit `BEGIN` and not yet committed fails to resolve, and the failing statement then invalidates that transaction. This applies to ordinary tables as much as to TEMP ones — it is a property of the recipe, not of TEMP, and it is not fixed by inheriting the TEMP catalog.
- The caller's catalog **search path** (`SET search_path`) is not inherited either, so a bare name that resolves for the user through a non-default schema will not resolve on the helper connection.
- **Session variables (`SET VARIABLE`, `getvariable()`) are not visible in the separate connection.** Views that use `getvariable()` will see NULL values. This is a fundamental limitation — each connection has its own session state. Document this in user-facing docs when relevant.
- For most use cases (reading from persistent tables/views), this is fine
- Read placements/references at bind time if possible, storing results in bind data to avoid issues during finalize
- See `src/placement_table_reader.cpp` for a complete example
