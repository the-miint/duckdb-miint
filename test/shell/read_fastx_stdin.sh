#!/bin/bash
# Tests for read_fastx stdin functionality

set -e

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
        echo "  ✓ PASS"
    else
        echo "  ✗ FAIL"
        echo "  Expected substring: $expected"
        echo "  Got: $result"
        FAILED=1
    fi
}

# Test 1: stdin with /dev/stdin (single-end FASTQ)
run_test "Single-end FASTQ stdin with /dev/stdin" \
    "│            2 │" \
    "cat data/fastq/small_a.fq | $DUCKDB -c \"SELECT COUNT(*) FROM read_fastx('/dev/stdin');\""

# Test 2: stdin with dash (single-end FASTQ)
run_test "Single-end FASTQ stdin with dash" \
    "│            2 │" \
    "cat data/fastq/small_a.fq | $DUCKDB -c \"SELECT COUNT(*) FROM read_fastx('-');\""

# Test 3: stdin with FASTA format
run_test "FASTA stdin" \
    "│            2 │" \
    "cat data/fastq/test.fa | $DUCKDB -c \"SELECT COUNT(*) FROM read_fastx('/dev/stdin');\""

# Test 4: Verify sequence_index is 1-based from stdin
run_test "Verify sequence_index starts at 1 for stdin" \
    "│.*1.*│" \
    "cat data/fastq/small_a.fq | $DUCKDB -c \"SELECT sequence_index FROM read_fastx('/dev/stdin') ORDER BY sequence_index LIMIT 1;\""

# Test 5: Verify all records are read correctly from stdin
run_test "Verify read_id from stdin" \
    "read_a2" \
    "cat data/fastq/small_a.fq | $DUCKDB -c \"SELECT read_id FROM read_fastx('/dev/stdin') ORDER BY sequence_index LIMIT 1 OFFSET 1;\""

# Test 6: Include filepath parameter with stdin (/dev/stdin)
run_test "Include filepath with /dev/stdin" \
    "/dev/stdin" \
    "cat data/fastq/small_a.fq | $DUCKDB -c \"SELECT filepath FROM read_fastx('/dev/stdin', include_filepath=true) LIMIT 1;\""

# Test 7: Include filepath parameter with stdin (dash)
run_test "Include filepath with dash" \
    "/dev/stdin" \
    "cat data/fastq/small_a.fq | $DUCKDB -c \"SELECT filepath FROM read_fastx('-', include_filepath=true) LIMIT 1;\""

# Test 8: Error - stdin in multi-file array with /dev/stdin
run_test "Error: stdin in multi-file array (/dev/stdin)" \
    "stdin" \
    "$DUCKDB -c \"SELECT * FROM read_fastx(['/dev/stdin', 'data/fastq/small_a.fq']);\" 2>&1 || true"

# Test 9: Error - stdin in multi-file array with dash
run_test "Error: stdin in multi-file array (dash)" \
    "stdin" \
    "$DUCKDB -c \"SELECT * FROM read_fastx(['-', 'data/fastq/small_a.fq']);\" 2>&1 || true"

# Test 10: Error - stdin with sequence2 (/dev/stdin)
run_test "Error: stdin with sequence2 (/dev/stdin)" \
    "stdin cannot be used with sequence2" \
    "$DUCKDB -c \"SELECT * FROM read_fastx('/dev/stdin', sequence2='data/fastq/small_a_r2.fq');\" 2>&1 || true"

# Test 11: Error - stdin with sequence2 (dash)
run_test "Error: stdin with sequence2 (dash)" \
    "stdin cannot be used with sequence2" \
    "$DUCKDB -c \"SELECT * FROM read_fastx('-', sequence2='data/fastq/small_a_r2.fq');\" 2>&1 || true"

# Summary
echo ""
if [ $FAILED -eq 0 ]; then
    echo "All tests passed!"
    exit 0
else
    echo "Some tests failed!"
    exit 1
fi
