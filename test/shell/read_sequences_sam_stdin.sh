#!/bin/bash
# read_sequences_sam from stdin: '-' and '/dev/stdin' return what reading the file returns, BAM and SAM alike, with a
# max_batch_bytes small enough to hold back every read but a chunk's first, and stdin is refused in a list of paths.

set -e

# The checks grep for duckbox output (the CLI default). Pin it: DuckDB 2.0's CLI renders markdown in agent mode
# (AI_AGENT, CLAUDECODE, etc. set and stdout not a terminal). There is no -duckbox flag; -box is a different format.
# With several threads a single file is read with an htslib thread pool and stdin without one, so each comparison
# below crosses the two read paths.
DUCKDB="./build/release/duckdb -cmd '.mode duckbox' -cmd 'SET threads=8'"
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

# Every column that varies, in order, from stdin and from the file.
same_as_file() { # same_as_file <stdin path> <file> [named parameters]
    echo "SELECT (SELECT list((sequence_index, read_id, sequence1, qual1) ORDER BY sequence_index) FROM read_sequences_sam('$1'$3))
               = (SELECT list((sequence_index, read_id, sequence1, qual1) ORDER BY sequence_index) FROM read_sequences_sam('$2'))
               AND (SELECT count(*) FROM read_sequences_sam('$2')) = 3 AS same;"
}

run_test "BAM from stdin with dash" \
    "true" \
    "cat data/sam/ubam_ccs_golden.bam | $DUCKDB -c \"$(same_as_file - data/sam/ubam_ccs_golden.bam)\""

run_test "SAM from stdin with /dev/stdin" \
    "true" \
    "cat data/sam/ubam_no_sq.sam | $DUCKDB -c \"$(same_as_file /dev/stdin data/sam/ubam_no_sq.sam)\""

run_test "BAM from stdin, one read per chunk" \
    "true" \
    "cat data/sam/ubam_ccs_golden.bam | $DUCKDB -c \"$(same_as_file - data/sam/ubam_ccs_golden.bam ", max_batch_bytes = '1B'")\""

run_test "Error: stdin in a list of paths" \
    "must be a single file path" \
    "$DUCKDB -c \"SELECT * FROM read_sequences_sam(['-', 'data/sam/ubam_no_sq.sam']);\" 2>&1 || true"

# Summary
echo ""
if [ $FAILED -eq 0 ]; then
    echo "All tests passed!"
    exit 0
else
    echo "Some tests failed!"
    exit 1
fi
