#!/usr/bin/env bash
# Regenerate golden BLAST output from the native sortmerna CLI for a given
# query + ref. Goldens are hard-coded in test/cpp/test_SortMeRNAAligner.cpp;
# rerun this when test inputs or the embedded sortmerna version changes, and
# copy the BLAST columns into the test literals.
#
# This script is NOT invoked at test time — tests compare the library output
# against the embedded golden constants directly.
#
# Usage:
#   sortmerna_oracle.sh <sortmerna_bin> <ref_fasta> <query_fasta> <out_blast_path>
#
# Optional env:
#   SMR_ORACLE_NUM_ALIGNMENTS  (default: 1)
#   SMR_ORACLE_THREADS         (default: 1)
#
# The workdir is a fresh temp directory and is removed on exit. The oracle
# intentionally uses native (not library) sortmerna so its output is
# independent of the code under test.
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <sortmerna_bin> <ref_fasta> <query_fasta> <out_blast_path>" >&2
    exit 64
fi

SORTMERNA_BIN="$1"
REF_FASTA="$2"
QUERY_FASTA="$3"
OUT_BLAST="$4"
NUM_ALIGNMENTS="${SMR_ORACLE_NUM_ALIGNMENTS:-1}"
THREADS="${SMR_ORACLE_THREADS:-1}"

for f in "$SORTMERNA_BIN" "$REF_FASTA" "$QUERY_FASTA"; do
    if [[ ! -f "$f" ]]; then
        echo "sortmerna_oracle: required file not found: $f" >&2
        exit 66
    fi
done

WORKDIR="$(mktemp -d -t sortmerna-oracle-XXXXXX)"
trap 'rm -rf "$WORKDIR"' EXIT

# Redirect sortmerna's own logs to the workdir so test harness stdout stays clean.
"$SORTMERNA_BIN" \
    --ref "$REF_FASTA" \
    --reads "$QUERY_FASTA" \
    --workdir "$WORKDIR" \
    --aligned "$WORKDIR/aligned" \
    --blast "1 cigar qcov" \
    --num_alignments "$NUM_ALIGNMENTS" \
    --threads "$THREADS" \
    >"$WORKDIR/stdout.log" 2>"$WORKDIR/stderr.log"

if [[ ! -f "$WORKDIR/aligned.blast" ]]; then
    echo "sortmerna_oracle: aligned.blast not produced" >&2
    echo "--- stdout ---" >&2
    cat "$WORKDIR/stdout.log" >&2
    echo "--- stderr ---" >&2
    cat "$WORKDIR/stderr.log" >&2
    exit 70
fi

mkdir -p "$(dirname "$OUT_BLAST")"
cp "$WORKDIR/aligned.blast" "$OUT_BLAST"
