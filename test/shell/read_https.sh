#!/bin/bash
# Tests for reading files via HTTP URLs
# Requires MIINT_HTTPS_TEST_URL to be set (run_tests.sh sets this automatically)

set -e

DUCKDB="./build/release/duckdb"
FAILED=0

if [ -z "$MIINT_HTTPS_TEST_URL" ]; then
    echo "MIINT_HTTPS_TEST_URL not set, skipping HTTPS tests"
    exit 0
fi

# Check if httpfs is available (use -c mode so DuckDB exits non-zero on failure)
if ! $DUCKDB -c "LOAD httpfs; SELECT 1;" > /dev/null 2>&1; then
    echo "httpfs not available, skipping HTTPS tests"
    exit 0
fi

BASE="$MIINT_HTTPS_TEST_URL"

run_test() {
    local test_name="$1"
    local expected="$2"
    local sql="$3"

    echo "Running: $test_name"

    result=$($DUCKDB -csv -noheader -c "LOAD httpfs; $sql" 2>&1 || true)

    if [ "$result" = "$expected" ]; then
        echo "  PASS"
    else
        echo "  FAIL"
        echo "  Expected: $expected"
        echo "  Got: $result"
        FAILED=1
    fi
}

# Verify a SQL statement produces an error containing the expected substring
run_error_test() {
    local test_name="$1"
    local expected_substring="$2"
    local sql="$3"

    echo "Running: $test_name"

    result=$($DUCKDB -c "LOAD httpfs; $sql" 2>&1 || true)

    if echo "$result" | grep -q "$expected_substring"; then
        echo "  PASS"
    else
        echo "  FAIL"
        echo "  Expected substring: $expected_substring"
        echo "  Got: $result"
        FAILED=1
    fi
}

# --- read_fastx ---

run_test "read_fastx: FASTQ via HTTP" "2" \
    "SELECT COUNT(*) FROM read_fastx('$BASE/fastq/small_a.fq');"

run_test "read_fastx: FASTA via HTTP" "2" \
    "SELECT COUNT(*) FROM read_fastx('$BASE/fastq/test.fa');"

run_test "read_fastx: gzipped FASTQ via HTTP" "2" \
    "SELECT COUNT(*) FROM read_fastx('$BASE/fastq/foo.r1.fastq.gz');"

run_test "read_fastx: paired-end gzipped FASTQ via HTTP" "2,2,2" \
    "SELECT COUNT(*), COUNT(sequence1), COUNT(sequence2) FROM read_fastx('$BASE/fastq/foo.r1.fastq.gz', sequence2='$BASE/fastq/foo.r2.fastq.gz');"

run_error_test "read_fastx: HTTP 404" "HTTP" \
    "SELECT * FROM read_fastx('$BASE/fastq/nonexistent.fq');"

# --- read_alignments ---

run_test "read_alignments: SAM via HTTP" "4" \
    "SELECT COUNT(*) FROM read_alignments('$BASE/sam/foo_has_header.sam');"

run_test "read_alignments: BAM via HTTP" "4" \
    "SELECT COUNT(*) FROM read_alignments('$BASE/sam/foo_has_header.bam');"

run_test "read_alignments: BAM with tags via HTTP" "2" \
    "SELECT COUNT(*) FROM read_alignments('$BASE/sam/foo_with_tags.bam');"

run_test "read_alignments: multiple SAM files via HTTP" "8" \
    "SELECT COUNT(*) FROM read_alignments(['$BASE/sam/foo_has_header.sam', '$BASE/sam/foo_has_header_2.sam']);"

run_error_test "read_alignments: HTTP 404" "Error" \
    "SELECT * FROM read_alignments('$BASE/sam/nonexistent.sam');"

# --- read_sequences_sam ---

run_test "read_sequences_sam: uBAM via HTTP" "3" \
    "SELECT COUNT(*) FROM read_sequences_sam('$BASE/sam/ubam_no_sq.sam');"

run_test "read_sequences_sam: multiple files via HTTP" "5" \
    "SELECT COUNT(*) FROM read_sequences_sam(['$BASE/sam/ubam_no_sq.sam', '$BASE/sam/ubam_no_sq_2.sam']);"

# --- read_newick ---

run_test "read_newick: Newick via HTTP" "5" \
    "SELECT COUNT(*) FROM read_newick('$BASE/newick/simple.nwk');"

run_test "read_newick: gzipped Newick via HTTP" "5" \
    "SELECT COUNT(*) FROM read_newick('$BASE/newick/simple.nwk.gz');"

# --- read_mzml ---

run_test "read_mzml: mzML via HTTP" "3" \
    "SELECT COUNT(*) FROM read_mzml('$BASE/mzml/basic_3spectra.mzML');"

run_test "read_mzml: multiple mzML files via HTTP" "4" \
    "SELECT COUNT(*) FROM read_mzml(['$BASE/mzml/basic_3spectra.mzML', '$BASE/mzml/compressed.mzML']);"

run_error_test "read_mzml: HTTP 404" "HTTP" \
    "SELECT * FROM read_mzml('$BASE/mzml/nonexistent.mzML');"

# --- read_mzxml ---

run_test "read_mzxml: mzXML via HTTP" "3" \
    "SELECT COUNT(*) FROM read_mzxml('$BASE/mzxml/basic_3spectra.mzXML');"

run_test "read_mzxml: multiple mzXML files via HTTP" "4" \
    "SELECT COUNT(*) FROM read_mzxml(['$BASE/mzxml/basic_3spectra.mzXML', '$BASE/mzxml/compressed.mzXML']);"

# --- read_sequences_sff ---

run_test "read_sequences_sff: SFF via HTTP" "2" \
    "SELECT COUNT(*) FROM read_sequences_sff('$BASE/sff/basic_2reads.sff');"

run_test "read_sequences_sff: multiple SFF files via HTTP" "3" \
    "SELECT COUNT(*) FROM read_sequences_sff(['$BASE/sff/basic_2reads.sff', '$BASE/sff/single_read.sff']);"

# --- read_biom (only if HDF5 available) ---

if [ -n "$HDF5_AVAILABLE" ]; then
    run_test "read_biom: BIOM via HTTP" "15" \
        "SELECT COUNT(*) FROM read_biom('$BASE/biom/test.biom');"
else
    echo "Skipping read_biom HTTPS test (HDF5 not available)"
fi

# --- Summary ---

if [ "$FAILED" -ne 0 ]; then
    echo "HTTPS tests FAILED"
    exit 1
fi

echo "All HTTPS tests passed!"
