# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## ABSOLUTE HARD REQUIREMENTS
- NEVER use `rm` without permission

## Project Overview

**duckdb-miint** is a DuckDB extension for bioinformatics that enables SQL-based analysis of genomic sequence files (FASTA/FASTQ) and alignment files (SAM/BAM). Built on the DuckDB extension template, it integrates HTSlib for BAM/SAM processing and kseq++ for FASTQ/FASTA parsing.

The extension provides:
- Table functions to read bioinformatics file formats as DuckDB tables
- Scalar functions for alignment analysis (flag checking, sequence identity)
- Aggregate functions for genomic interval operations
- Custom COPY formats for writing bioinformatics files

## Priorities

1. Test Driven Development (TDD)
2. Correct code as verified by tests
3. Maintainable code, using Don't Repeat Yourself (DRY) and Keep It Simple Stupid (KISS)
4. Performance

## Build Commands

### Initial Setup (VCPKG for dependencies)
```bash
git clone https://github.com/Microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh
export VCPKG_TOOLCHAIN_PATH=`pwd`/vcpkg/scripts/buildsystems/vcpkg.cmake
```

### Build the extension
```bash
bash build.sh
```

This creates:
- `./build/release/duckdb` - DuckDB shell with extension preloaded
- `./build/release/test/unittest` - DuckDB test runner with extension linked
- `./build/release/extension/miint/miint.duckdb_extension` - Loadable extension binary
- `./build/release/tests` - C++ unit tests (Catch2)

### Testing

If a test produces an **incorrect expected value**: DO NOT change the expected value without permission.

```bash
# Run all SQL tests
bash run_tests.sh

# Run C++ unit tests
./build/release/extension/miint/tests

# Run a single SQL test file
./build/release/test/unittest "test/sql/read_alignments.test"

# Run tests matching a pattern
./build/release/test/unittest "[alignment]"
```

### Development workflow
```bash
# Clean build
make clean

# Build with specific configuration
make release  # or make debug
```

## Architecture Overview

### Extension Entry Point
- **src/miint_extension.cpp**: Main extension registration. The `LoadInternal()` function registers all table functions, scalar functions, and COPY formats. This is where new functionality should be registered.

### Core Components

#### 1. Table Functions (Reading Files)
Table functions allow querying bioinformatics files as SQL tables.

- **read_fastx**: FASTA/FASTQ sequence files
  - Implementation: `src/read_fastx.cpp`, `src/include/read_fastx.hpp`
  - Reader: `src/SequenceReader.cpp` (wraps kseq++)
  - Record: `src/include/SequenceRecord.hpp`
  - Returns: sequence_index, read_id, comment, sequence1, sequence2, qual1, qual2, [filepath]
  - Supports: single/paired-end, multiple files, stdin, quality score offset conversion

- **read_alignments**: SAM/BAM alignment files (previously `read_sam`, which still works as an alias)
  - Implementation: `src/read_alignments.cpp`, `src/include/read_alignments.hpp`
  - Reader: `src/SAMReader.cpp` (wraps HTSlib)
  - Record: `src/include/SAMRecord.hpp`
  - Returns: read_id, flags, reference, position, stop_position, mapq, cigar, mate info, optional tags, [filepath]
  - Supports: headerless SAM (with reference table), multiple files, parallel processing, both SAM and BAM formats

#### 2. Scalar Functions
Individual functions for alignment analysis:

- **Alignment flags**: `src/alignment_flag_functions.cpp`
  - Functions like `alignment_is_paired()`, `alignment_is_unmapped()`, `alignment_is_primary()`
  - Test individual SAM flag bits

- **Alignment analysis**: `src/alignment_functions.cpp`
  - `alignment_seq_identity()`: Calculate sequence identity (gap_compressed, gap_excluded, blast methods)

#### 3. Aggregate Functions

- **compress_intervals**: `src/compress_intervals.cpp`, `src/IntervalCompressor.cpp`
  - Merges overlapping genomic intervals
  - Uses efficient algorithm with periodic compression for large datasets
  - Thread-safe for parallel GROUP BY operations

#### 4. COPY Formats (Writing Files)
Custom formats for exporting query results to bioinformatics files:

- **FORMAT FASTQ**: `src/copy_fastq.cpp`
  - Parameters: QUAL_OFFSET, INCLUDE_COMMENT, ID_AS_SEQUENCE_INDEX, INTERLEAVE, COMPRESSION
  - Supports: single/paired-end, interleaved, split files with {ORIENTATION} placeholder

- **FORMAT FASTA**: `src/copy_fasta.cpp`
  - Parameters: INCLUDE_COMMENT, ID_AS_SEQUENCE_INDEX, INTERLEAVE, COMPRESSION
  - Similar structure to FASTQ without quality scores

- **FORMAT SAM / FORMAT BAM**: `src/copy_sam.cpp`
  - Parameters: INCLUDE_HEADER, REFERENCE_LENGTHS, COMPRESSION (SAM only), COMPRESSION_LEVEL (BAM only)
  - Requires reference_lengths table for header generation
  - BAM format always requires headers (binary format specification)
  - COMPRESSION_LEVEL for BAM: 0-9, default 6 (BGZF compression)
  - Currently writes SEQ/QUAL as `*` (future enhancement opportunity)

- **Common utilities**: `src/copy_format_common.cpp`
  - Shared buffering and compression infrastructure
  - Reference table reading: `src/reference_table_reader.cpp`

### Key Design Patterns

#### File Reading Pattern
All table functions follow similar structure:
1. **Data struct**: Stores configuration (paths, parameters, field definitions)
2. **GlobalState struct**: Manages file readers, thread coordination, file iteration
3. **Bind()**: Validates parameters, returns schema
4. **InitGlobal()**: Opens files, creates readers
5. **Execute()**: Reads records in chunks, populates DataChunk output

MaxThreads() controls parallelism:
- `read_fastx`: 1 for stdin, 8 for files
- `read_alignments`: 4 threads

#### Record Abstraction
Each file format has a record struct (SAMRecord, SequenceRecord) that:
- Wraps underlying C library objects (HTSlib bam1_t, kseq record)
- Provides type-safe field access via enums
- Handles memory management with RAII

#### Reference Table Pattern
For headerless SAM files and SAM writing:
- User provides DuckDB table name via parameter
- `reference_table_reader.cpp` executes query, extracts name/length columns
- Validates reference names per SAM spec (no `*`, `=`, tabs, newlines, position patterns)
- Returns `unordered_map<string, uint64_t>` for header construction

### External Dependencies

- **HTSlib 1.22.1** (`ext/htslib-1.22.1/`): SAM/BAM/CRAM parsing
  - Built as ExternalProject in CMake
  - Configured with zlib and optional zstd support
  - Static library linked into extension

- **kseq++** (`ext/kseq++/`): Modern C++ FASTA/FASTQ parser
  - Header-only library
  - Included directly in compilation

- **VCPKG dependencies**:
  - zlib (required)
  - zstd (optional, for compression)
  - Catch2 (C++ testing framework)

### Testing Strategy

#### SQL Tests (`test/sql/`)
Primary test mechanism using DuckDB's test framework:
- Each `.test` file contains SQL statements and expected outputs
- Format: `statement ok`, `query <types>`, `----` separators
- Coverage: error handling, basic functionality, edge cases, data types, ordering

Test file structure:
1. Error handling tests (missing files, invalid parameters)
2. Basic functionality tests (single/multi-file, different input types)
3. Feature-specific tests (include_filepath, compression, etc.)
4. Edge cases (empty files, large values, NULL handling)
5. Data type verification

#### C++ Tests (`test/cpp/`)
Unit tests for core components using Catch2:
- `test_SequenceReader.cpp`: FASTQ/FASTA parsing
- `test_SAMReader.cpp`: SAM/BAM parsing, headerless support
- `test_QualScore.cpp`: Quality score conversion
- `test_IntervalCompressor.cpp`: Interval merging algorithm
- `test_AlignmentFunctions.cpp`: Sequence identity calculations

Use `TEST_CASE()` and `SECTION()` for test organization.

### Important Implementation Details

#### Thread Safety
- **SAMReader/SequenceReader**: NOT thread-safe for concurrent calls on same instance
- **Multiple reader instances**: Thread-safe (each has independent file handles)
- **GlobalState**: Uses mutex for file iteration coordination
- **compress_intervals**: Thread-safe aggregate (each thread has independent state)

#### Headerless SAM Support
Critical implementation in `SAMReader.cpp`:
- Constructor accepts `unordered_map<string, uint64_t>` for reference lengths
- Synthetic header created via `sam_hdr_add_line()`
- Validation: reference names checked per SAM spec (no `*`, `=`, special chars, position patterns)
- File position preserved (no reopen) for stdin support
- Missing references in data → records appear as unmapped with reference `*`

#### Quality Score Handling
`src/QualScore.cpp` provides:
- Offset detection (Phred33 vs Phred64)
- Conversion between offsets
- Validation (scores must fit in valid range)

#### Compression Support
All COPY formats support compression:
- Auto-detection from `.gz` extension
- Manual override with `COMPRESSION` parameter
- Buffered writing for performance

#### Stop Position Calculation
SAM alignments compute `stop_position` using HTSlib's `bam_endpos()`:
- Accounts for CIGAR operations (M, D, N, =, X)
- 1-based inclusive coordinate
- Critical for interval operations and coverage analysis

#### Reading from Tables and Views in Extensions
When extension code needs to read data from a user-specified table or view (e.g., PLACEMENTS parameter in COPY FORMAT NEWICK, or REFERENCE_LENGTHS in COPY FORMAT SAM), use the following pattern:

**The Problem:**
- `Catalog::GetEntry<TableCatalogEntry>` only works for tables, not views
- Views are stored query definitions without physical storage
- `context.Query()` causes deadlocks when called during bind or execution (context is already locked)

**The Solution: Use a Separate Connection**
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

**Schema Validation for Tables/Views:**
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

**Important Notes:**
- The separate connection runs in its own transaction, so uncommitted changes from the original context may not be visible
- For most use cases (reading from persistent tables/views), this is fine
- Read placements/references at bind time if possible, storing results in bind data to avoid issues during finalize
- See `src/placement_table_reader.cpp` for a complete example

## Code Style and Conventions

- **C++ Standard**: C++20 (required for kseq++)
- **Namespace**: All extension code in `duckdb` namespace, readers in `miint` namespace
- **Error Handling**: Use DuckDB exceptions (`InvalidInputException`, `IOException`, `RuntimeException`)
- **Memory Management**: RAII with smart pointers (HTSlib uses custom deleters)
- **Formatting**: Follow DuckDB style (clang-format config in `.clang-format`)

## Adding New Functionality

### Adding a new table function
1. Create header in `src/include/` with TableFunction class
2. Implement Bind(), InitGlobal(), Execute() in `src/`
3. Register in `LoadInternal()` in `miint_extension.cpp`
4. Add to `EXTENSION_SOURCES` in `CMakeLists.txt`
5. Create SQL test in `test/sql/`

### Adding a new COPY format
1. Create `copy_<format>.cpp` and header
2. Implement CopyFunction with Bind(), InitGlobal(), Sink(), Finalize()
3. Use `copy_format_common.cpp` utilities for buffering/compression
4. Register in `LoadInternal()`
5. Add to `EXTENSION_SOURCES`
6. Create comprehensive SQL tests (basic, compression, edge cases)

### Adding scalar/aggregate functions
1. Implement in appropriate source file
2. Create static Register() method
3. Call Register() in LoadInternal()
4. Add SQL and/or C++ tests

## Common Issues and Solutions

### Build Issues
- **zstd not found**: Optional dependency, build continues without it
- **HTSlib build fails**: Check zlib installation, ensure CFLAGS=-fPIC set
- **VCPKG errors**: Ensure VCPKG_TOOLCHAIN_PATH is set before running cmake

### Runtime Issues
- **"File lacks a header"**: SAM file is headerless, provide `reference_lengths` parameter
- **"Inconsistent headers across files"**: All SAM files must be either header or headerless, not mixed
- **Unknown reference**: Reference in SAM data not in header or reference_lengths table
- **Quality offset errors**: Specify `qual_offset` parameter if auto-detection fails

### Testing Issues
- **Test isolation**: Each test file should be independent, use temp tables for references
- **Data files**: Test data in `data/sam/`, `data/fastq/`, etc. Keep files small
- **Platform differences**: Be aware of path separators, newline conventions
