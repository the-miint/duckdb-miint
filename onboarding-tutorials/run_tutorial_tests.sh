#!/bin/bash
# Run tutorial code blocks extracted from markdown files.
#
# Usage:
#   bash onboarding-tutorials/run_tutorial_tests.sh [duckdb_path]
#
# duckdb_path defaults to ./build/release/duckdb
#
# For each tutorial:
#   - beginner.md: concatenate all ```python blocks into one script, run it
#   - intermediate.md: concatenate all ```sql blocks into one session, run via duckdb
#   - advanced.md: run ```bash blocks independently, ```sql blocks as one session
#
# The INSTALL/LOAD miint statements are rewritten to use the local build.

set -euo pipefail

DUCKDB="${1:-./build/release/duckdb}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${TMPDIR:-/tmp}/miint_tutorial_tests"
mkdir -p "$OUTDIR"

PASS=0
FAIL=0

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

extract_blocks() {
    # Extract fenced code blocks of a given language from a markdown file.
    # Usage: extract_blocks <file> <language>
    # Outputs the concatenated content of all matching blocks.
    local file="$1"
    local lang="$2"
    awk -v lang="$lang" '
        /^```'"$lang"'/ { capture=1; next }
        /^```/ { if (capture) { capture=0; print "" } ; next }
        capture { print }
    ' "$file"
}

rewrite_extension_loading() {
    # Replace community install with local extension load.
    # Works for both SQL and Python (con.sql("...")) contexts.
    sed \
        -e "s|INSTALL miint FROM community; ||g" \
        -e "s|INSTALL miint FROM community;||g"
}

# --------------------------------------------------------------------------
# Test: beginner.md (Python)
# --------------------------------------------------------------------------

echo "=== Testing beginner.md ==="

BEGINNER_PY="$OUTDIR/beginner_test.py"
extract_blocks "$SCRIPT_DIR/beginner.md" "python" \
    | rewrite_extension_loading \
    > "$BEGINNER_PY"

PYTHON="${PYTHON:-conda run -n duckdb-151 python3}"
if $PYTHON "$BEGINNER_PY" 2>&1; then
    echo "PASS: beginner.md"
    PASS=$((PASS + 1))
else
    echo "FAIL: beginner.md"
    FAIL=$((FAIL + 1))
fi

# --------------------------------------------------------------------------
# Test: intermediate.md (SQL)
# --------------------------------------------------------------------------

echo ""
echo "=== Testing intermediate.md ==="

INTERMEDIATE_SQL="$OUTDIR/intermediate_test.sql"
extract_blocks "$SCRIPT_DIR/intermediate.md" "sql" \
    | rewrite_extension_loading \
    > "$INTERMEDIATE_SQL"

if "$DUCKDB" < "$INTERMEDIATE_SQL" 2>&1; then
    echo "PASS: intermediate.md"
    PASS=$((PASS + 1))
else
    echo "FAIL: intermediate.md"
    FAIL=$((FAIL + 1))
fi

# --------------------------------------------------------------------------
# Test: advanced.md
# --------------------------------------------------------------------------

echo ""
echo "=== Testing advanced.md (SQL blocks) ==="

ADVANCED_SQL="$OUTDIR/advanced_test.sql"
extract_blocks "$SCRIPT_DIR/advanced.md" "sql" \
    | rewrite_extension_loading \
    > "$ADVANCED_SQL"

if "$DUCKDB" < "$ADVANCED_SQL" 2>&1; then
    echo "PASS: advanced.md (SQL)"
    PASS=$((PASS + 1))
else
    echo "FAIL: advanced.md (SQL)"
    FAIL=$((FAIL + 1))
fi

echo ""
echo "=== Testing advanced.md (bash blocks) ==="

# Extract bash blocks, rewrite duckdb references and extension loading
ADVANCED_BASH="$OUTDIR/advanced_test.sh"
extract_blocks "$SCRIPT_DIR/advanced.md" "bash" \
    | sed "s|duckdb |${DUCKDB} |g" \
    | rewrite_extension_loading \
    > "$ADVANCED_BASH"

if bash "$ADVANCED_BASH" 2>&1; then
    echo "PASS: advanced.md (bash)"
    PASS=$((PASS + 1))
else
    echo "FAIL: advanced.md (bash)"
    FAIL=$((FAIL + 1))
fi

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------

echo ""
echo "==============================="
echo "Tutorial tests: $PASS passed, $FAIL failed"
echo "==============================="

if [ "$FAIL" -gt 0 ]; then
    exit 1
fi
