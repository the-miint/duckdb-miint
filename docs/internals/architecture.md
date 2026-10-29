# Architecture (developer reference)

Cross-cutting design patterns and implementation notes for duckdb-miint. This file is for developers working on the extension itself.

**Not here:**
- Per-function API reference → the task-organized user docs under `docs/` (index: `docs/table_of_contents.md`)
- Authoritative list of what's registered → `src/miint_extension.cpp` (`LoadInternal()`)
- How external libraries are embedded → `docs/internals/embedded-tools.md`
- Reading from user-specified tables/views → `docs/internals/reading-tables-views.md`
- Consuming Arrow C Data Interface streams → `docs/internals/arrow-zero-copy.md`
- Adding `sample_id` support to a table function → `docs/internals/per-sample-pattern.md`

`docs/` is published at https://the-miint.github.io/duckdb-miint/. The site (`the-miint/the-miint.github.io`) copies `docs/` at build time, except `UPDATING.md`, and serves `table_of_contents.md` as the section's index page. Its navigation is a hand-maintained `nav` in its own `mkdocs.yml`, and `mkdocs build --strict` fails when a synced page is missing from that nav or an anchor link is broken, so a new, renamed or deleted page here needs a matching change there. Downstream projects read these pages as miint's contract.

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

`MaxThreads()` controls parallelism per function (e.g. `read_fastx` uses 1 for stdin and one per file, or per pair of files, otherwise). It returns the most threads the scan can use, and DuckDB's `threads` setting bounds it. A file reader therefore returns its file count, that count capped (several readers stop at 8), or a fixed number (`read_biom` returns 4), and does not cap it at `std::thread::hardware_concurrency()`: that can return 0, and it would hold the scan below a `threads` setting larger than the core count.

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
`src/include/table_function_common.hpp` / `.cpp` provides `ExpandGlobPattern` / `ExpandGlobPatternWithInfo` (a VARCHAR path or glob) and `ExpandPathList` (a VARCHAR[] of literal paths), which expand every path as DuckDB does (a `file:` URI to its local path, `~` to home), plus `IsStdinPath`, `SetResultVector*`, and similar helpers. Prefer these over re-implementing per-function: a reader that reads its path argument itself skips the expansion.

A reader whose rows each come from one line of a text file (`read_paf`, `read_gfa`) builds on `src/include/line_file_scan.hpp`. `MakeLineScanFunction` declares the table function with its shared global state and a file-size row estimate; `BindLineScan` resolves the path, glob or list argument and `include_filepath`; `ScanLines` claims files per thread, reads them through DuckDB's file system (gzip, bgzip and remote paths alike), caps a chunk at 16 MiB of lines, and reports a parse error as `<function>: <file>: line N: <defect>`. The reader supplies only a DuckDB-free line parser and the code that writes one row. Its parse step is kept separate from its write step so that only parse errors are labelled with a line.

### Dual-path Table Functions (standard + lateral)
Functions like `read_ena_sequences` set both `function` (`Execute`) and `in_out_function` (`ExecuteInOut`) on the same `TableFunction`. DuckDB routes scalar-constant calls to the standard path (parallel, byte-based progress) and correlated-column or subquery calls to the in-out path (one outer row at a time, `LIMIT`-driven short-circuit via `LocalState` destruction). Shared per-run open/read/close logic lives in `miint::PerRunReader` so both paths use the same state machine and retry policy, and shared column-fill logic lives in a `FillOutputFromBatch` helper so both paths emit identical chunk layouts. The `Bind` function detects the in-out path by an empty `input.inputs` (DuckDB's table-in-out dispatcher passes no scalar constants) and sets `Data::deferred_resolution`; the per-query LRU `ENAResolverCache` on `Data` then dedupes metadata lookups across outer rows.

### Filter pushdown via `pushdown_complex_filter`
A pushed-down path must return exactly the rows of the unpushed scan that pass the filter, value for value: none the scan does not produce, and none of its matching rows missing. Leaving `filters` in place, so that DuckDB re-applies them above the scan, removes rows that fail the filter; it cannot catch a matching row the scan would not produce, or a matching row the remote drops. For example, ENA's Portal API names sample attributes differently from the sample XML that `read_ena_attributes` reads (`country` where the XML says `geographic location (country and/or sea)`), so answering `WHERE tag = 'country'` from the Portal API returns a row the unfiltered scan does not contain. Test a pushed-down path by comparing it with the unpushed scan on data where the two could differ. No table function uses the callback.

To register one, set `tf.pushdown_complex_filter = YourCallback;` in `GetFunction()`. The two filter-pushdown flags are independent:
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
- FASTQ, FASTA, SAM and NEWICK detect gzip output from a `.gz` extension; `COMPRESSION` overrides it.
- BAM and UBAM are always BGZF-compressed.
- BIOM's `COMPRESSION` sets HDF5's internal gzip compression; there is no `.gz` detection.

### Stop Position Calculation
SAM alignments compute `stop_position` using HTSlib's `bam_endpos()`:
- Accounts for CIGAR operations (M, D, N, =, X)
- Half-open coordinate: `stop_position = position + reference_length_from_cigar` (exclusive end)
- Example: `10M` at position 2 → `stop_position = 12`, covering bases `[2, 12)`
- Coverage length = `stop_position - position` (no `+1`)
- Critical for interval operations and coverage analysis

### Coordinates for Readers
A reader converts every interval to miint's 1-based half-open convention when it reads, on a reference and on a read alike, and names the columns for the sequence they are on: `position` / `stop_position` on a reference, `query_start` / `query_stop` on a read. A column that leaves the file's convention in place means every caller has to know the source format to use it correctly — the off-by-one this rule exists to prevent. A column named for neither, such as a bare `start` / `stop`, does not say which sequence it is on. The rule, the per-format conversions and the functions listed as exceptions are in [Coordinate conventions](../reading.md#coordinate-conventions).

### Identifier-column Codec (`id_column_codec` / `id_column_utils`)
Functions that accept user-supplied identifier columns (`read_id`, `reference`, `mate_reference`) may opt in to BIGINT and UUID in addition to VARCHAR; `IsAllowedIdType` in `src/include/id_column_utils.hpp` defines the accepted set. The pattern is split across two layers:

- `src/include/id_column_codec.hpp` / `src/id_column_codec.cpp` — pure C++ (`miint::ParseIdAsInt64`, `miint::FormatIdFromInt64`). No DuckDB types, so the `tests` C++ binary can link it directly without pulling in DuckDB.
- `src/include/id_column_utils.hpp` — DuckDB-aware inline header (`duckdb::IsAllowedIdType`, `ExtractIdColumnAsStrings`, `EmitIdColumnFromStrings`). Wraps the codec and handles `Vector` / `DataChunk` / `LogicalType` plumbing.

In-memory carriers (e.g. `SequenceRecordBatch::read_ids`, `SAMRecordBatch::read_ids`/`references`/`mate_references`) always store `std::vector<std::string>` — BIGINT inputs are stringified at the ingress boundary and parsed back at the egress boundary. This keeps the alignment hot loop type-agnostic.

Dispatch is opt-in. Schema validators take an `allow_bigint` flag (default `false`), which admits UUID as well as BIGINT; a caller that doesn't pass `true` keeps the historical VARCHAR-only contract. The captured `LogicalType` is threaded through the bind data so the emit side knows which output type to materialise. Default-constructed `LogicalType` is `INVALID`, so any path that forgets to capture the type fails loud at the helper dispatch rather than silently producing garbage.

Sentinel rules (BIGINT egress): empty string and the SAM `*` sentinel both decode to NULL. The SAM `=` mate-reference sentinel has no BIGINT encoding and is the caller's responsibility to resolve to a real reference value before invoking the codec — `ParseIdAsInt64` throws on `=` to surface that contract. VARCHAR sentinels pass through unchanged.

Currently used by `align_minimap2` + `save_minimap2_index` (both sides), `copy_sam` (write side, decimal stringification only), `align_minimap2_sharded` (query side only — subject side stays VARCHAR because `.mmi` files store subject names as opaque bytes), `align_bowtie2` / `align_bowtie2_sharded` (the daemon wire schema is always strings — BIGINT support is C++/DuckDB-only at the ingress/egress boundaries), and `align_sortmerna` / `align_sortmerna_rrna` (query side only — subject names come from FASTA files on disk via `ref_paths`, never a user table, so `reference` / `mate_reference` / `ref_name` stay VARCHAR). The sharded paths extend to `ReadShardIds` and `ReadBatchByIds` (`sequence_table_reader.cpp`), which thread the captured `LogicalType` so the per-shard `_batch_ids` temp table is declared with the matching type, and to `ValidateReadToShardSchema` (`align_common.hpp`), which takes an optional `expected_read_id_type` parameter and enforces strict equality between `read_to_shard.read_id` and `query_table.read_id`. The bowtie2 daemon-IPC emit path lives in `bt2_daemon::EmitChunkRows` (`align_bowtie2_daemon_common.cpp`); it dispatches the three id columns through small file-local helpers (`emit_id_cell` / `read_arrow_string_owned`) rather than `id_column_utils.hpp` directly — this translation unit also references `duckdb::miint::gpl_boundary`, which previously collided with `id_column_utils.hpp`'s unqualified `miint::` references; PR B resolved the collision by switching the header to absolute `::miint::` qualifiers. The SortMeRNA emit path lives in two parallel-but-unmerged helpers in `src/include/sortmerna_result_utils.hpp` — `OutputSortMeRNASamBatch` (21-column SAM schema, mirrors `OutputSAMRecordBatch`) and `OutputSortMeRNARRNABatch` (13-column rRNA schema with identity / coverage / e-value). The two are deliberately not merged: their output schemas differ enough that a shared body would either need conditional dispatch on every column or accept dead-code branches in each direction. The remaining VARCHAR-only callers — mafft, fasttree, rype, deblur, ena_upload_reads — still pass `allow_bigint=false`. Future extensions that want BIGINT support should reuse these helpers rather than re-implementing dispatch.

## Code Style and Conventions

- **C++ Standard**: C++20 (required for kseq++)
- **Namespaces**: extension code in `duckdb`, readers in `miint`
- **Error Handling**: DuckDB exceptions (`InvalidInputException`, `IOException`, `RuntimeException`)
- **Memory Management**: RAII with smart pointers (HTSlib uses custom deleters)
- **Formatting**: DuckDB style (`.clang-format`); pre-commit hook runs `make format-check`, auto-fix with `conda run -n duckdb-format make format-fix`

## Adding New Functionality

### Adding a new table function
1. Header in `src/include/` with TableFunction class
2. Implement `Bind()`, `InitGlobal()`, `Execute()` in `src/`
3. Register in `LoadInternal()` in `miint_extension.cpp`
4. Add to `MIINT_EXTENSION_SOURCES` in `CMakeLists.txt`, or in the feature's `cmake/deps` file if it needs an optional library ([build-system.md](build-system.md))
5. Create SQL test in `test/sql/`
6. If the function should run per-sample (partitioning an input relation by a column), follow `docs/internals/per-sample-pattern.md`
7. If an argument names a relation to read, follow `docs/internals/reading-tables-views.md`, including its [CTE guard](reading-tables-views.md#ctes-reject-them-in-bind) and the matching `cte_relation_name_guard` test case

<a name="no-work-in-bind"></a>
### Do the work in `InitGlobal` / `Execute` — never in `Bind`

**`Bind` resolves the schema. It does not compute the answer.** Reading an input
relation, running a compute, building an index: all of that belongs at execution
time. This is easy to get wrong because doing it in `Bind` *works* — the results are
correct, the tests pass, and the cost only shows up in four ways that no test looks
for unless you write one:

1. **A prepared statement returns nothing the second time.** `PREPARE` binds once and
   every `EXECUTE` re-runs that same plan. A `Bind`-computed result therefore has to
   be handed to the scan either by *moving* it out of the bind data — so the second
   `EXECUTE` finds an empty vector and returns **zero rows, silently** — or by
   copying it, paying for the whole result twice. Both `pcoa`/`unifrac_pcoa` and the
   progressive pair shipped with the first bug. Rows a `Bind` reads from a relation
   are also the relation as it was at `PREPARE`, however it changes before `EXECUTE`.
2. **`EXPLAIN` pays full price.** `EXPLAIN` binds without executing, so work in `Bind`
   runs for a query the user never ran. Measured on a 4,000-sample ordination:
   `EXPLAIN` took 0.644 s against 0.638 s to actually run the query. After moving it,
   0.000 s.
3. **The query cannot be cancelled.** Ctrl-C sets the context's interrupt flag
   (`context.IsInterrupted()`), which only helps if something polls it. Nothing
   polls during `Bind`, so a multi-hour run is uninterruptible. Execution-time work
   can poll between units and throw `InterruptException` (see
   `ProgressivePcoaExecute`).
4. **Nothing streams, and the result exists twice.** A `Bind` that computes everything
   must hand over a finished result, which the scan then re-materializes in its own
   row form. Driving the work from the scan instead lets rows go out as they are
   produced and keeps one copy.

**What `Bind` legitimately does**: validate parameters, and resolve the output schema
— including rejecting a mis-shaped input. That last part does not require reading
data. `ProbeFeatureTableIdType` / `ProbeDistanceTableIdType`
(`unifrac_function_common`) are the pattern: a `LIMIT 0` probe proves the columns
exist and cast, and a catalog lookup gives the id type, with no rows read. So a
missing column or a mismatched id type still fails at bind, where a schema error
belongs, while errors about the *data* surface when the query runs.

Consequences to expect when you move work out of `Bind`: some errors change from
`BinderException` to `InvalidInputException` and are raised later (fine — assert on
the message, not the exception class), and per-execution work means an *unseeded*
function redraws on re-execution instead of replaying its first result (also fine,
and more correct).

`pcoa`, `unifrac_pcoa`, `progressive_pcoa_from_distances` and
`progressive_pcoa_from_unifrac` in `src/unifrac_pcoa_function.cpp` are the worked
examples; the progressive pair additionally shows the streaming form, where the scan
steps the computation one wave at a time.

### Adding a new COPY format
1. Create `copy_<format>.cpp` and header
2. Implement `CopyFunction` with `Bind()`, `InitGlobal()`, `Sink()`, `Finalize()`
3. Declare every option in a `copy_options` callback with `DeclareCopyOption` (`src/include/copy_option_parsing.hpp`), and read them in `Bind()` with that header's helpers, ending the if-chain with `ThrowUnreadCopyOption`. Without `copy_options`, DuckDB checks no option: an unknown one reaches `Bind()`, and a bare one arrives with no value.
4. Use `copy_format_common.cpp` utilities for buffering/compression
5. Register in `LoadInternal()`
6. Add to `MIINT_EXTENSION_SOURCES` ([build-system.md](build-system.md))
7. Create comprehensive SQL tests (basic, compression, edge cases)

### Adding scalar/aggregate functions
1. Implement in appropriate source file
2. Create a static `Register()` method
3. Call `Register()` in `LoadInternal()`
4. Add SQL and/or C++ tests
