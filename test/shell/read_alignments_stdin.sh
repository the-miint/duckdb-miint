#!/bin/bash
# Tests for read_alignments/read_sam stdin functionality

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

# Test 1: Headerless stdin with reference_lengths using /dev/stdin
run_test "Headerless stdin with /dev/stdin" \
    "│            4 │" \
    "cat data/sam/foo_no_header.sam | $DUCKDB -c \"CREATE TABLE ref AS SELECT 'G1234' AS name, 20 AS length UNION ALL SELECT 'G000144735', 90; SELECT COUNT(*) FROM read_sam('/dev/stdin', reference_lengths='ref');\""

# Test 2: Headerless stdin with reference_lengths using dash
run_test "Headerless stdin with dash" \
    "│            4 │" \
    "cat data/sam/foo_no_header.sam | $DUCKDB -c \"CREATE TABLE ref AS SELECT 'G1234' AS name, 20 AS length UNION ALL SELECT 'G000144735', 90; SELECT COUNT(*) FROM read_sam('-', reference_lengths='ref');\""

# Test 3: Header stdin without reference_lengths
run_test "Header stdin without reference_lengths" \
    "│            4 │" \
    "cat data/sam/foo_has_header.sam | $DUCKDB -c \"SELECT COUNT(*) FROM read_sam('/dev/stdin');\""

# Test 4: Verify all records are read correctly from stdin
run_test "Verify all records read correctly" \
    "foo-3" \
    "cat data/sam/foo_no_header.sam | $DUCKDB -c \"CREATE TABLE ref AS SELECT 'G1234' AS name, 20 AS length UNION ALL SELECT 'G000144735', 90; SELECT read_id FROM read_sam('/dev/stdin', reference_lengths='ref') ORDER BY read_id LIMIT 1 OFFSET 2;\""

# Test 5: Error - stdin in multi-file array with /dev/stdin
run_test "Error: stdin in multi-file array (/dev/stdin)" \
    "Cannot use stdin" \
    "$DUCKDB -c \"SELECT * FROM read_sam(['/dev/stdin', 'data/sam/foo_has_header.sam']);\" 2>&1 || true"

# Test 6: Error - stdin in multi-file array with dash
run_test "Error: stdin in multi-file array (dash)" \
    "Cannot use stdin" \
    "$DUCKDB -c \"SELECT * FROM read_sam(['-', 'data/sam/foo_has_header.sam']);\" 2>&1 || true"

# Test 7: Include filepath parameter with stdin
run_test "Include filepath with stdin" \
    "/dev/stdin" \
    "cat data/sam/foo_no_header.sam | $DUCKDB -c \"CREATE TABLE ref AS SELECT 'G1234' AS name, 20 AS length UNION ALL SELECT 'G000144735', 90; SELECT filepath FROM read_sam('/dev/stdin', reference_lengths='ref', include_filepath=true) LIMIT 1;\""

# Summary
echo ""
if [ $FAILED -eq 0 ]; then
    echo "All tests passed!"
    exit 0
else
    echo "Some tests failed!"
    exit 1
fi
