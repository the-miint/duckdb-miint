# Architecture (developer reference)

Cross-cutting design patterns and implementation notes for duckdb-miint. This file is for developers working on the extension itself.

**Not here:**
- Per-function API reference → `docs/table-functions.md`, `docs/scalar-functions.md`, `docs/copy-formats.md`, `docs/analysis-functions.md`
- Authoritative list of what's registered → `src/miint_extension.cpp` (`LoadInternal()`)
- How external libraries are embedded → `docs/internals/embedded-tools.md`
- Reading from user-specified tables/views → `docs/internals/reading-tables-views.md`
- Consuming Arrow C Data Interface streams → `docs/internals/arrow-zero-copy.md`

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
