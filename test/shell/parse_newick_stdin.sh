#!/bin/bash
# Tests for parse_newick stdin functionality

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

# Test 1: stdin with /dev/stdin
run_test "Parse newick from stdin with /dev/stdin" \
    "│            5 │" \
    "cat data/newick/simple.nwk | $DUCKDB -c \"SELECT COUNT(*) FROM parse_newick('/dev/stdin');\""

# Test 2: stdin with dash
run_test "Parse newick from stdin with dash" \
    "│            5 │" \
    "cat data/newick/simple.nwk | $DUCKDB -c \"SELECT COUNT(*) FROM parse_newick('-');\""

# Test 3: Verify tip names from stdin
run_test "Verify tip names from stdin" \
    "A" \
    "cat data/newick/simple.nwk | $DUCKDB -c \"SELECT name FROM parse_newick('/dev/stdin') WHERE is_tip = true ORDER BY name LIMIT 1;\""

# Test 4: Verify tip count from stdin
run_test "Verify tip count from stdin" \
    "│            3 │" \
    "cat data/newick/simple.nwk | $DUCKDB -c \"SELECT COUNT(*) FROM parse_newick('/dev/stdin') WHERE is_tip = true;\""

# Test 5: Include filepath parameter with stdin (/dev/stdin)
run_test "Include filepath with /dev/stdin" \
    "/dev/stdin" \
    "cat data/newick/simple.nwk | $DUCKDB -c \"SELECT filepath FROM parse_newick('/dev/stdin', include_filepath=true) LIMIT 1;\""

# Test 6: Include filepath parameter with stdin (dash)
run_test "Include filepath with dash" \
    "/dev/stdin" \
    "cat data/newick/simple.nwk | $DUCKDB -c \"SELECT filepath FROM parse_newick('-', include_filepath=true) LIMIT 1;\""

# Test 7: Tree with edge IDs from stdin
run_test "Tree with edge IDs from stdin" \
    "│            5 │" \
    "cat data/newick/with_edge_ids.nwk | $DUCKDB -c \"SELECT COUNT(*) FROM parse_newick('/dev/stdin') WHERE edge_id IS NOT NULL;\""

# Test 8: Gzipped input from stdin (use gunzip -c for portability across Linux/macOS)
run_test "Gzipped newick from stdin" \
    "│            5 │" \
    "gunzip -c data/newick/simple.nwk.gz | $DUCKDB -c \"SELECT COUNT(*) FROM parse_newick('/dev/stdin');\""

# Test 9: Error - stdin in multi-file array with /dev/stdin
run_test "Error: stdin in multi-file array (/dev/stdin)" \
    "stdin" \
    "$DUCKDB -c \"SELECT * FROM parse_newick(['/dev/stdin', 'data/newick/simple.nwk']);\" 2>&1 || true"

# Test 10: Error - stdin in multi-file array with dash
run_test "Error: stdin in multi-file array (dash)" \
    "stdin" \
    "$DUCKDB -c \"SELECT * FROM parse_newick(['-', 'data/newick/simple.nwk']);\" 2>&1 || true"

# Summary
echo ""
if [ $FAILED -eq 0 ]; then
    echo "All tests passed!"
    exit 0
else
    echo "Some tests failed!"
    exit 1
fi
