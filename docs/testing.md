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

# The same tests, with the SQL files and the C++ tests on 8 parallel workers
MIINT_TEST_JOBS=8 bash run_tests.sh

# Run all SQL tests
make test

# Run a single test file
./build/release/test/unittest "test/sql/read_alignments.test"

# Run tests matching a pattern
./build/release/test/unittest "[alignment]"
```

`MIINT_TEST_JOBS=N` runs the files `make test` runs, each in its own `unittest` process, alongside the C++ test binary (`test/scripts/run_tests_parallel.py`). It prints a line per file as each finishes, then the output of every file that failed, and reports skips per reason as the serial pass does. A job still running after ten minutes is killed and fails the run, so a deadlocked test cannot stall it; each SQL file is one job, and so is the whole C++ test binary. The descriptor-limit rerun and the shell tests stay serial. Its SQL jobs run with `RUST_MIN_STACK=4194304` unless that is already set. Under heavy load on a many-core machine, the workers of rayon's global pool, which rype classifies on, can overflow Rust's default 2 MiB thread stack ([RYpe#28](https://github.com/the-miint/RYpe/issues/28)). The parallel jobs make that load themselves; a serial run meets it only when other heavy work shares the machine. CI leaves `MIINT_TEST_JOBS` unset and runs the serial pass.

`build/release/duckdb`, `test/unittest`, `extension/miint/tests`, `extension/miint/tests_duckdb` and the loadable `miint.duckdb_extension` are separate build targets, so after an incremental `ninja` of one of them the others can be stale. Rebuild every binary you probe with, e.g. `ninja -C build/release duckdb test/unittest extension/miint/tests extension/miint/tests_duckdb`.

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
- A test that reads through a glob and asserts `filepath` values selects `replace(filepath, '\', '/')`: glob expansion returns native separators on Windows. A single path or an explicit path list keeps the caller's string and needs no wrapping.
- Write files under `{TEST_DIR}`, the per-process scratch directory unittest substitutes into a test's SQL. Under `MIINT_TEST_JOBS`, test files and the C++ tests run at the same time in one working directory, so any other path a test writes must be one no other test reads or writes. For the same reason the parallel runner gives each SQL process a spill directory under its `{TEST_DIR}` (`docs/internals/duckdb-engine-notes.md`, "Storage").

**Guards and skipped tests:**

`require-env` guards only what a deployment can lack: an optional embedded library, a runtime binary, a network resource or served data. The one other category is process state that the pass over every test file, serial or parallel, cannot set, such as a lowered descriptor limit, which `run_tests.sh` provides in a rerun. `run_tests.sh` decides which guards are exported. A test of this repository's own C++ code is never guarded: an oracle produced by an external tool is committed as a fixture, and a missing or corrupted fixture fails the test. A checksum gate in `run_tests.sh` would turn a corrupted golden into a skip and leave the run green, so a parity test that pins its golden's checksum checks it with an assertion.

A skipped file is not a passing one, and the suite's skip count is part of its result: a run that skips more files than before covers less, whatever its summary line says. Three skips are easy to miss:

- DuckDB's `unittest` skips a test whose error contains "HTTP" or "Unable to connect" unless its config names other messages. `run_tests.sh` exports `DUCKDB_TEST_CONFIG=test/configs/no_skipped_errors.json`, under which every error fails the test; set it the same way when running `unittest` by hand.
- `require httpfs` skips whenever the build does not include httpfs, and `extension_config.cmake` does not build it on this DuckDB version, so `make test` skips every `require httpfs` file: the HTTP-reader tests, the ENA mock-server tests and the live ENA and NCBI tests. `run_tests.sh` reruns them after `make test`, in what it calls the httpfs pass, with `DUCKDB_TEST_AUTOLOADING=all`, `DUCKDB_TEST_ON_INIT='LOAD httpfs;'` and `HOME=build/release/httpfs_test_home`, into which it installs httpfs once, and passes unittest `--keep-home`, without which DuckDB 2.0's unittest replaces `HOME` with an empty sandbox of its own. Autoloading runs miint from the copy installed in that `HOME`, not the linked one, so the script installs the fresh build there with `FORCE INSTALL` on every run (a file copied over an earlier install does not load, because the install records the build's version). The pass fails when httpfs or miint does not load, when any file in its list matches no test, when a file skips for anything but a `require-env` guard of the kinds above (an optional embedded library, a runtime binary, a network resource) or one on Webin credentials or an opt-in such as `MIINT_BLAST_LIVE_TEST`, or when it outlasts `MIINT_HTTPFS_TIMEOUT_SECONDS` (default 1800), since a stalled transfer otherwise hangs it. Once `run_tests.sh` has set that `HOME` up, one file runs with `HOME=$PWD/build/release/httpfs_test_home DUCKDB_TEST_AUTOLOADING=all DUCKDB_TEST_ON_INIT='LOAD httpfs;' DUCKDB_TEST_CONFIG=$PWD/test/configs/no_skipped_errors.json ./build/release/test/unittest --keep-home test/sql/<file>.test`; after a rebuild, install it there first with `HOME=$PWD/build/release/httpfs_test_home ./build/release/duckdb -unsigned -c "FORCE INSTALL '$PWD/build/release/extension/miint/miint.duckdb_extension'"`. Several files go in one comma-separated argument: separate arguments join into a single filter that matches nothing, and the run reports "No tests ran" with exit status 0.
- The live ENA tests (`require-env ENA_AVAILABLE`, and `test/shell/read_ena_sequences_aspera.sh`) skip unless `MIINT_ENA_LIVE_TEST=1`, which CI sets only in its nightly run. EBI's HTTP 500s, empty answers and md5 mismatches failed 18 of the 175 attempts of CI's test job between 2026-08-01 and 2026-10-22.
- `test/shell/load_miint_concurrent.sh` (many databases in one process loading the extension at once) skips unless `MIINT_LOAD_RACE_TEST=1`, which CI sets only in its nightly run, with 100 processes: it takes minutes, and on a 4-CPU runner the race it guards against hits about one process in 30. With the opt-in set, `run_tests.sh` fails unless `MIINT_TEST_PYTHON` (default `python3`) has the duckdb matching this build.
  - The rest of `run_tests.sh` sends miint's Portal, checklist and Reports requests to the mock `test/scripts/ena_webin_mock.py`, through the base URLs `test/scripts/ena_mock_env.sh` exports. The Portal base covers every Portal request, so a live test run with it would read from the mock; the live tests run after the httpfs pass instead, without those variables and with `ENA_AVAILABLE`, which nothing else sets.
  - On every run, `read_ena_sequences_mock.test` reads runs through the mock, from `data/ena/`, `read_ena_sequences_open_errors_mock.test` checks the warnings for a run whose file cannot be opened, and `read_ena_input_errors.test` checks the ENA read functions' input errors, which need no request. `read_ena_sequences_resume_mock.test` sets DuckDB's `http_retries`, `http_retry_wait_ms` (1 ms) and `http_retry_backoff` itself, and reads the runs that recover before those that never do.
  - `onboarding-tutorials/run_tutorial_tests.sh` starts the same mock unless `MIINT_TUTORIAL_LIVE=1`, which CI also sets only nightly. Live, each tutorial gets three attempts, a minute apart (`MIINT_TUTORIAL_LIVE_ATTEMPTS` changes the count), since EBI also drops FASTQ transfers partway; a pass after a failed attempt shows as a warning on the run. Against the mock a tutorial gets one attempt. The tutorials still fetch a genome from NCBI and install httpfs from DuckDB's extension repository on every run.

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

`tests` links no DuckDB. A C++ test that needs a running database, one that opens an in-process `DuckDB` and drives miint code against it, goes in the second runner instead, `./build/release/extension/miint/tests_duckdb`, which links DuckDB and miint statically. List its file under `tests_duckdb` in `CMakeLists.txt`, not in `TEST_SOURCES`; a file that needs an optional feature, such as HDF5, is added there with `target_sources` under that feature's switch. It is built only where DuckDB builds `unittest`, so not in the WASM builds. `run_tests.sh` runs both.

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

A script that loads httpfs (`read_https.sh`, `read_ena_sequences_aspera.sh`) skips itself when the shell cannot load it, as on any machine with no httpfs under `~/.duckdb`, CI included. `run_tests.sh` runs such a script with the httpfs pass's `HOME=build/release/httpfs_test_home`, where the shell has already loaded httpfs, and picks such scripts by the text `load httpfs`, matched case-insensitively. By hand, run one with `HOME=$PWD/build/release/httpfs_test_home` once `run_tests.sh` has set that `HOME` up.

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
- `test/shell/cron_publish_extension.sh` - Scenario tests for run selection in `scripts/cron-publish-extension.sh` (stubbed `gh`, `mail`, `ssh`, `rsync` and `hostname`; skips without jq, flock and GNU coreutils)

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
