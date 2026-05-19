# Architecture (developer reference)

Cross-cutting design patterns and implementation notes for duckdb-miint. This file is for developers working on the extension itself.

**Not here:**
- Per-function API reference → `docs/table-functions.md`, `docs/scalar-functions.md`, `docs/copy-formats.md`, `docs/analysis-functions.md`
- Authoritative list of what's registered → `src/miint_extension.cpp` (`LoadInternal()`)
- How external libraries are embedded → `docs/internals/embedded-tools.md`
- Reading from user-specified tables/views → `docs/internals/reading-tables-views.md`
- Consuming Arrow C Data Interface streams → `docs/internals/arrow-zero-copy.md`
- Adding `sample_id` support to a table function → `docs/internals/per-sample-pattern.md`

## Extension Entry Point

`src/miint_extension.cpp` — `LoadInternal()` registers every table function, scalar function, aggregate, and COPY format. New functionality is registered here.

## Key Design Patterns

### File Reading Pattern
All table functions follow a similar structure:
1. **Data struct**: Stores configuration (paths, parameters, field definitions)
2. **GlobalState struct**: Manages file readers, thread coordination, file iteration
3. **Bind()**: Validates parameters, returns schema
4. **InitGlobal()**: Opens files, creates readers
5. **Execute()**: Reads records in chunks, populates DataChunk output

`MaxThreads()` controls parallelism per function (e.g. `read_fastx` uses 1 for stdin / 8 for files; `read_alignments` uses 4).

### Record Abstraction
Each file format has a record struct (SAMRecord, SequenceRecord, etc.) that:
- Wraps underlying C library objects (HTSlib `bam1_t`, kseq record)
- Provides type-safe field access via enums
- Handles memory management with RAII

### Reference Table Pattern
For headerless SAM files and SAM writing:
- User provides DuckDB table name via parameter
- `src/reference_table_reader.cpp` executes the query and extracts name/length columns
- Validates reference names per SAM spec (no `*`, `=`, tabs, newlines, position patterns)
- Returns `unordered_map<string, uint64_t>` for header construction

### Shared Table Function Utilities
`src/include/table_function_common.hpp` / `.cpp` provides `ParseFilePathsParameter`, `IsStdinPath`, `SetResultVector*`, and similar helpers. Prefer these over re-implementing per-function.

### Dual-path Table Functions (standard + lateral)
Functions like `read_ena_sequences` set both `function` (`Execute`) and `in_out_function` (`ExecuteInOut`) on the same `TableFunction`. DuckDB routes scalar-constant calls to the standard path (parallel, byte-based progress) and correlated-column or subquery calls to the in-out path (one outer row at a time, `LIMIT`-driven short-circuit via `LocalState` destruction). Shared per-run open/read/close logic lives in `miint::PerRunReader` so both paths use the same state machine and retry policy, and shared column-fill logic lives in a `FillOutputFromBatch` helper so both paths emit identical chunk layouts. The `Bind` function detects the in-out path by an empty `input.inputs` (DuckDB's table-in-out dispatcher passes no scalar constants) and sets `Data::deferred_resolution`; the per-query LRU `ENAResolverCache` on `Data` then dedupes metadata lookups across outer rows.

### Filter pushdown via `pushdown_complex_filter`
`read_ena_attributes` is the reference implementation for translating a DuckDB `WHERE` clause into a remote-API query shape. Pattern:

1. **Abstract AST for the extractor.** `src/include/ena_attributes_filter.hpp` defines `miint::ENAFilterNode` (EQ / IN / AND / OR / UNSUPPORTED). `ExtractPushdownPredicates(vector<unique_ptr<ENAFilterNode>> conjuncts) → ENAAttributePushdown` is pure, has no DuckDB dependency, and is exhaustively unit-tested against the AST directly — no binder scaffolding required.
2. **DuckDB → AST translator.** Lives in the table function's `.cpp` (see `read_ena_attributes.cpp`'s anonymous namespace). Handles `BoundComparisonExpression` (COMPARE_EQUAL), `BoundOperatorExpression` (COMPARE_IN), and `BoundConjunctionExpression` (AND / OR). Anything else maps to `ENAFilterNode::MakeUnsupported()`, and the extractor rejects any conjunct list containing an unsupported atom.
3. **Column name resolution.** `BoundColumnRefExpression::binding.column_index` is a local index into `LogicalGet::GetColumnIds()`; `GetPrimaryIndex()` on that entry gives the position in `LogicalGet::names`. (Mirrors `multi_file_list.cpp:68-69`.)
4. **Do not erase filters.** Record the predicate in `Data::pushdown`, but leave the `vector<unique_ptr<Expression>> &filters` unchanged. DuckDB re-applies the filter as a `LogicalFilter` above the scan, so any semantic mismatch between the remote matcher and SQL equality degrades to extra work, never to wrong rows. (See `read_ena_attributes.cpp`'s `PushdownComplexFilter`.)
5. **Dual-path Execute.** `GlobalState` exposes both `FetchNextBatch` (default) and `FetchNextStructuredBatch` (pushdown). `Execute()` branches on `bind_data.pushdown.tags.empty()` at the refill site. Both paths write to the same `current_batch` member so the emission code is shared.
6. **Allowlist-gated pushdown.** The extractor consults `ENASearchFieldRegistry::IsSearchableSampleField` and returns empty (XML fallback) if any referenced tag is unknown — even inside an `IN` list. No per-tag mixing of structured + XML.

Registration: set `tf.pushdown_complex_filter = YourCallback;` in `GetFunction()`. The two filter-pushdown flags are independent:
- `tf.pushdown_complex_filter` (callback) — arbitrary Expression trees, called once during optimization, scope decided by the callback.
- `tf.filter_pushdown = true` (bool) — tells the optimizer to extract simple `col=const` / `LIKE` / `IN` predicates into `get.table_filters`. The scan then receives those via `TableFunctionInitInput::filters` and is expected to consume them at read time (e.g. parquet zonemap pruning). Setting this without implementing consumption does not drop filters for correctness (the optimizer still wraps remaining expressions in a `LogicalFilter` above the scan), but the populated `TableFilterSet` is dead weight. Leave `filter_pushdown = false` unless you actually consume `TableFilter` objects.

### GPL/BSD License Boundary (gpl-boundary subsystem)

Some upstream tools we want to expose are GPL-licensed and cannot be statically linked into miint without changing the extension's BSD licensing. The pattern: link the GPL code into a separate process (`gpl-boundary`), launch it as a child, and communicate via JSON-line control + POSIX shared memory + Arrow IPC bytes.

`src/phylogeny_fasttree.cpp` (FastTree) is the first user. The wire format and process plumbing live in `src/gpl_boundary/` (`process`, `session`, `shm`, `arrow_ipc`); see [`docs/internals/embedded-tools.md`](embedded-tools.md) for the protocol invariants and lifetime rules. Two non-obvious constraints:

1. **SIGPIPE handling is per-thread**, not process-wide. The session blocks SIGPIPE via `pthread_sigmask` and drains pending signals with `sigtimedwait` — never `sigaction(SIGPIPE, SIG_IGN)`, which would mask SIGPIPE for every other thread in the host process.
2. **Output shm is unlinked by miint, not the daemon.** This is the gpl-boundary README's contract; if we ever both unlink, the second unlink races with a shm name reused by another query.

Adding a second tool to the same daemon (or a second daemon for a different GPL library) requires defining a new `output_schema` on the daemon side, then registering a new table function on the miint side that mirrors `phylogeny_fasttree.cpp`. The session/shm/IPC plumbing is generic.

## Testing Strategy

### SQL Tests (`test/sql/`)
Primary test mechanism using DuckDB's test framework.
- Each `.test` file: `statement ok`, `query <types>`, `----` separators
- Coverage: error handling, basic functionality, edge cases, data types, ordering

Typical test-file structure:
1. Error handling (missing files, invalid parameters)
2. Basic functionality (single/multi-file, different input types)
3. Feature-specific tests (include_filepath, compression, etc.)
4. Edge cases (empty files, large values, NULL handling)
5. Data type verification

### C++ Tests (`test/cpp/`)
Unit tests for core components using Catch2 — e.g. `test_SequenceReader.cpp`, `test_SAMReader.cpp`, `test_QualScore.cpp`, `test_IntervalCompressor.cpp`, `test_AlignmentFunctions.cpp`. Use `TEST_CASE()` and `SECTION()` for organization.

## Important Implementation Details

### Thread Safety
- **SAMReader / SequenceReader**: NOT thread-safe for concurrent calls on the same instance
- **Multiple reader instances**: thread-safe (each has independent file handles)
- **GlobalState**: uses a mutex for file iteration coordination
- **compress_intervals**: thread-safe aggregate (each thread has independent state)

### Headerless SAM Support
Implementation in `src/SAMReader.cpp`:
- Constructor accepts `unordered_map<string, uint64_t>` for reference lengths
- Synthetic header created via `sam_hdr_add_line()`
- Reference names validated per SAM spec (no `*`, `=`, special chars, position patterns)
- File position preserved (no reopen) for stdin support
- References present in data but missing from the map → records appear as unmapped with reference `*`

### Quality Score Handling
`src/QualScore.cpp`:
- Offset detection (Phred33 vs Phred64)
- Conversion between offsets
- Validation (scores must fit in valid range)

### Compression Support
All COPY formats:
- Auto-detect from `.gz` extension
- Manual override via the `COMPRESSION` parameter
- Buffered writing for performance

### Stop Position Calculation
SAM alignments compute `stop_position` using HTSlib's `bam_endpos()`:
- Accounts for CIGAR operations (M, D, N, =, X)
- Half-open coordinate: `stop_position = position + reference_length_from_cigar` (exclusive end)
- Example: `10M` at position 2 → `stop_position = 12`, covering bases `[2, 12)`
- Coverage length = `stop_position - position` (no `+1`)
- Critical for interval operations and coverage analysis

### Identifier-column Codec (`id_column_codec` / `id_column_utils`)
Functions that accept user-supplied identifier columns (`read_id`, `reference`, `mate_reference`) may opt in to BIGINT in addition to VARCHAR. The pattern is split across two layers:

- `src/include/id_column_codec.hpp` / `src/id_column_codec.cpp` — pure C++ (`miint::ParseIdAsInt64`, `miint::FormatIdFromInt64`). No DuckDB types, so the `tests` C++ binary can link it directly without pulling in DuckDB.
- `src/include/id_column_utils.hpp` — DuckDB-aware inline header (`duckdb::IsAllowedIdType`, `ExtractIdColumnAsStrings`, `EmitIdColumnFromStrings`). Wraps the codec and handles `Vector` / `DataChunk` / `LogicalType` plumbing.

In-memory carriers (e.g. `SequenceRecordBatch::read_ids`, `SAMRecordBatch::read_ids`/`references`/`mate_references`) always store `std::vector<std::string>` — BIGINT inputs are stringified at the ingress boundary and parsed back at the egress boundary. This keeps the alignment hot loop type-agnostic.

Dispatch is opt-in. Schema validators take an `allow_bigint` flag (default `false`); a caller that doesn't pass `true` keeps the historical VARCHAR-only contract. The captured `LogicalType` is threaded through the bind data so the emit side knows which output type to materialise. Default-constructed `LogicalType` is `INVALID`, so any path that forgets to capture the type fails loud at the helper dispatch rather than silently producing garbage.

Sentinel rules (BIGINT egress): empty string and the SAM `*` sentinel both decode to NULL. The SAM `=` mate-reference sentinel has no BIGINT encoding and is the caller's responsibility to resolve to a real reference value before invoking the codec — `ParseIdAsInt64` throws on `=` to surface that contract. VARCHAR sentinels pass through unchanged.

Currently used by `align_minimap2` + `save_minimap2_index` (both sides), `copy_sam` (write side, decimal stringification only), `align_minimap2_sharded` (query side only — subject side stays VARCHAR because `.mmi` files store subject names as opaque bytes), and `align_bowtie2` / `align_bowtie2_sharded` (the daemon wire schema is always strings — BIGINT support is C++/DuckDB-only at the ingress/egress boundaries). The sharded paths extend to `ReadShardIds` and `ReadBatchByIds` (`sequence_table_reader.cpp`), which thread the captured `LogicalType` so the per-shard `_batch_ids` temp table is declared with the matching type, and to `ValidateReadToShardSchema` (`align_common.hpp`), which takes an optional `expected_read_id_type` parameter and enforces strict equality between `read_to_shard.read_id` and `query_table.read_id`. The bowtie2 daemon-IPC emit path lives in `bt2_daemon::EmitChunkRows` (`align_bowtie2_daemon_common.cpp`); it dispatches the three id columns through small file-local helpers (`emit_id_cell` / `read_arrow_string_owned`) rather than `id_column_utils.hpp` directly — this translation unit also references `duckdb::miint::gpl_boundary`, which collides with `id_column_utils.hpp`'s unqualified `miint::` references (tracked as a separate cleanup). The remaining VARCHAR-only callers — sortmerna, mafft, fasttree, rype, deblur, ena_upload_reads — still pass `allow_bigint=false`. Future extensions that want BIGINT support should reuse these helpers rather than re-implementing dispatch.

## Code Style and Conventions

- **C++ Standard**: C++20 (required for kseq++)
- **Namespaces**: extension code in `duckdb`, readers in `miint`
- **Error Handling**: DuckDB exceptions (`InvalidInputException`, `IOException`, `RuntimeException`)
- **Memory Management**: RAII with smart pointers (HTSlib uses custom deleters)
- **Formatting**: DuckDB style (`.clang-format`); pre-commit hook runs `make format-check`, auto-fix with `conda run -n duckdb-143 make format-fix`

## Adding New Functionality

### Adding a new table function
1. Header in `src/include/` with TableFunction class
2. Implement `Bind()`, `InitGlobal()`, `Execute()` in `src/`
3. Register in `LoadInternal()` in `miint_extension.cpp`
4. Add to `EXTENSION_SOURCES` in `CMakeLists.txt`
5. Create SQL test in `test/sql/`
6. If the function should run per-sample (partitioning an input relation by a column), follow `docs/internals/per-sample-pattern.md`

### Adding a new COPY format
1. Create `copy_<format>.cpp` and header
2. Implement `CopyFunction` with `Bind()`, `InitGlobal()`, `Sink()`, `Finalize()`
3. Use `copy_format_common.cpp` utilities for buffering/compression
4. Register in `LoadInternal()`
5. Add to `EXTENSION_SOURCES`
6. Create comprehensive SQL tests (basic, compression, edge cases)

### Adding scalar/aggregate functions
1. Implement in appropriate source file
2. Create a static `Register()` method
3. Call `Register()` in `LoadInternal()`
4. Add SQL and/or C++ tests
