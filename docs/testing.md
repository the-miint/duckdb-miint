# Testing

The MIINT extension uses three complementary testing approaches to ensure correctness and reliability.

## Table of Contents

- [SQL Logic Tests](#sql-logic-tests)
- [C++ Unit Tests](#c-unit-tests)
- [Shell Script Tests](#shell-script-tests)
- [Test Data](#test-data)

## SQL Logic Tests

The primary testing mechanism for user-facing functionality uses DuckDB's SQL logic test framework. Test files are located in `test/sql/` and use the `.test` extension.

**Running SQL tests:**
```sh
# Run all tests (SQL, C++, and shell)
bash run_tests.sh

# Run all SQL tests
make test

# Run a single test file
./build/release/test/unittest "test/sql/read_alignments.test"

# Run tests matching a pattern
./build/release/test/unittest "[alignment]"
```

**Test file structure:**
SQL logic tests consist of statements and expected outputs. Each test file follows this format:

```sql
# name: test/sql/example.test
# description: Test description
# group: [sql]

require miint

# Test statement that should succeed
statement ok
SELECT * FROM read_fastx('data/fastq/example.fastq');

# Query test with expected output
query I
SELECT COUNT(*) FROM read_fastx('data/fastq/example.fastq');
----
100

# Test error handling
statement error
SELECT * FROM read_fastx('nonexistent.fastq');
----
File not found: nonexistent.fastq
```

**Test organization:**
- `statement ok` - Statement should execute without error
- `statement error` - Statement should raise an error (followed by expected error message)
- `query <types>` - Query should return specific results
  - Types: `I` (INTEGER), `R` (REAL/DOUBLE), `T` (TEXT), etc.
- `----` - Separator between query and expected output

**Example test files:**
- `test/sql/read_alignments.test` - SAM/BAM reading functionality
- `test/sql/read_fastx.test` - FASTA/FASTQ reading functionality
- `test/sql/alignment_functions.test` - Alignment analysis functions
- `test/sql/copy_sam.test` - SAM/BAM writing functionality

**Good practices for SQL tests:**
- Test error cases first (missing files, invalid parameters)
- Test basic functionality with simple inputs
- Test edge cases (empty files, NULL values, boundary conditions)
- Verify data types and column ordering
- Keep test data files small and focused
- Use clear, descriptive comments for each test section

## C++ Unit Tests

Core algorithms and internal utilities are tested using the Catch2 framework. C++ tests are located in `test/cpp/` and compile to `./build/release/extension/miint/tests`.

**Running C++ tests:**
```sh
# Build and run all C++ tests
./build/release/extension/miint/tests

# Run tests matching a specific tag
./build/release/extension/miint/tests "[alignment_functions]"

# Run a specific test case
./build/release/extension/miint/tests "ParseCigar - Basic operations"
```

**Code organization for testing:**

1. **`miint` namespace** - Pure C++ code (algorithms, parsers, utilities)
   - Location: `src/` files, internal headers
   - Tested with: Catch2 C++ unit tests (`test/cpp/`)
   - Examples: CIGAR parsing, MD tag parsing, interval compression

2. **`duckdb` namespace** - DuckDB integration code (table functions, scalar functions)
   - Location: `src/` files (e.g., `alignment_functions.cpp`, `read_fastx.cpp`)
   - Tested with: SQL logic tests (`test/sql/`)
   - Examples: Function registration, DuckDB vector operations, parameter binding

**Example C++ test structure:**
```cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace miint {
    // Implementation code to test
    struct CigarStats { /* ... */ };
    static CigarStats ParseCigar(const std::string &cigar_str);
}

TEST_CASE("ParseCigar - Basic operations", "[alignment_functions]") {
    SECTION("Simple match") {
        auto stats = miint::ParseCigar("10M");
        REQUIRE(stats.matches == 10);
        REQUIRE(stats.alignment_columns == 10);
    }

    SECTION("With insertions and deletions") {
        auto stats = miint::ParseCigar("10M5I3D10M");
        REQUIRE(stats.matches == 20);
        REQUIRE(stats.insertions == 5);
        REQUIRE(stats.deletions == 3);
        REQUIRE(stats.gap_opens == 2);
    }
}
```

**Example test files:**
- `test/cpp/test_AlignmentFunctions.cpp` - CIGAR/MD parsing, identity calculations
- `test/cpp/test_SequenceReader.cpp` - FASTQ/FASTA parsing
- `test/cpp/test_SAMReader.cpp` - SAM/BAM reading logic
- `test/cpp/test_IntervalCompressor.cpp` - Interval merging algorithm

**Good practices for C++ tests:**
- Use `TEST_CASE()` for grouping related tests
- Use `SECTION()` for organizing variations within a test case
- Use descriptive test names that explain what is being tested
- Use `REQUIRE()` for critical assertions, `CHECK()` for non-critical
- Use matchers for floating-point comparisons (`REQUIRE_THAT(value, WithinRel(expected, 0.001))`)
- Test boundary conditions, error cases, and edge cases
- Keep test data minimal and focused on the specific behavior being tested

## Shell Script Tests

Some functionality cannot be tested through SQL logic tests, such as stdin handling. These features are tested using bash scripts located in `test/shell/`.

**Running shell tests:**
```sh
# Run all shell tests
for test in test/shell/*.sh; do bash "$test"; done

# Run a specific shell test
bash test/shell/read_alignments_stdin.sh

# All tests (SQL, C++, and shell)
bash run_tests.sh
```

**Test file structure:**
Shell tests use simple bash scripts with conditional logic and exit codes:

```bash
#!/bin/bash
set -e  # Exit on error

DUCKDB="./build/release/duckdb"
FAILED=0

# Helper function to run a test
run_test() {
    local test_name="$1"
    local expected="$2"
    shift 2
    local cmd="$@"

    echo "Running: $test_name"
    result=$(eval "$cmd" 2>&1 || true)

    if echo "$result" | grep -q "$expected"; then
        echo "  PASS"
    else
        echo "  FAIL"
        FAILED=1
    fi
}

# Test cases
run_test "Test description" \
    "expected output substring" \
    "cat data.txt | $DUCKDB -c 'SELECT * FROM read_alignments(\"/dev/stdin\");'"

# Return status
if [ $FAILED -eq 0 ]; then
    echo "All tests passed!"
    exit 0
else
    echo "Some tests failed!"
    exit 1
fi
```

**Example test files:**
- `test/shell/read_alignments_stdin.sh` - Tests for reading SAM/BAM from stdin

**Good practices for shell tests:**
- Use `set -e` to exit on first error
- Test both success and error cases
- Use `grep -q` to verify expected output substrings
- Always return proper exit codes (0 for success, 1 for failure)
- Use clear, descriptive test names
- Capture both stdout and stderr (`2>&1`)

## Test Data

Test data files are organized in `data/` subdirectories:
- `data/sam/` - SAM/BAM test files
- `data/fastq/` - FASTQ/FASTA test files
- `data/biom/` - BIOM format test files
- `data/newick/` - Newick phylogenetic tree test files

**Important:** Keep test data files small (< 1KB when possible) to maintain fast test execution and minimize repository size.
