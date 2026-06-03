#!/usr/bin/env bash
#
# Stage real microbiome paired-end sequence data from ENA into a Parquet file for
# the COPY-format write benchmark (run_bench.sh). Uses the extension's own
# `read_ena_sequences` table function (over httpfs) rather than hand-rolled
# downloads, so it goes through the same code path the project ships.
#
# One-time + idempotent: re-running skips staging if data/staged.parquet exists.
#
# Default accession is a paired Illumina WGS metagenome run from study PRJEB11419
# (~7.8M read pairs, ~640MB compressed). Pick a different one with ACCESSION=...,
# or cap the read count with MAX_SEQUENCES for faster iteration.
#
# Env knobs:
#   ACCESSION             ENA run accession (default ERR10677119)
#   MAX_SEQUENCES         cap reads per run, 0 = full run (default 0)
#   MIINT_DUCKDB          duckdb shell with miint preloaded (default build/release/duckdb)
#   MIINT_EXT             if set, `LOAD '<path>'` is issued first (vanilla shell)
#   MIINT_BENCH_DATA_DIR  output dir (default ./data next to this script)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

DUCKDB="${MIINT_DUCKDB:-$REPO_ROOT/build/release/duckdb}"
DATA_DIR="${MIINT_BENCH_DATA_DIR:-$SCRIPT_DIR/data}"
STAGED="$DATA_DIR/staged.parquet"
STAGED_META="$DATA_DIR/staged.meta"
ACCESSION="${ACCESSION:-ERR10677119}"
MAX_SEQUENCES="${MAX_SEQUENCES:-0}"

mkdir -p "$DATA_DIR"

log() { printf '[fetch] %s\n' "$*" >&2; }
die() { printf '[fetch] ERROR: %s\n' "$*" >&2; exit 1; }

[[ -x "$DUCKDB" ]] || die "duckdb shell not found/executable: $DUCKDB (build it, or set MIINT_DUCKDB)"

LOAD_EXT=""
[[ -n "${MIINT_EXT:-}" ]] && LOAD_EXT="LOAD '${MIINT_EXT}';"

if [[ -f "$STAGED" && -f "$STAGED_META" ]]; then
	log "Staged parquet already present: $STAGED"
	sed 's/^/[fetch]   /' "$STAGED_META" >&2
	exit 0
fi

# Optional read cap.
max_clause=""
[[ "$MAX_SEQUENCES" != "0" ]] && max_clause=", max_sequences => $MAX_SEQUENCES"

log "Staging ENA run $ACCESSION -> $STAGED via read_ena_sequences (one-time) ..."
log "(streaming over httpfs; this hits EBI/ENA and may take a minute)"

"$DUCKDB" -batch :memory: 2>&1 <<SQL | sed 's/^/[fetch]   /' >&2
LOAD httpfs;
$LOAD_EXT
COPY (
    SELECT sequence_index, read_id, comment, sequence1, sequence2, qual1, qual2
    FROM read_ena_sequences('$ACCESSION'$max_clause)
) TO '$STAGED' (FORMAT PARQUET);
SQL

[[ -f "$STAGED" ]] || die "staging produced no parquet (see output above)"

log "Computing staged metadata ..."
"$DUCKDB" -batch -noheader -list :memory: > "$STAGED_META" <<SQL || die "metadata query failed"
$LOAD_EXT
SELECT 'accession=$ACCESSION';
SELECT 'rows=' || count(*) FROM read_parquet('$STAGED');
SELECT 'paired_rows=' || count(*) FROM read_parquet('$STAGED') WHERE sequence2 IS NOT NULL;
SELECT 'r1_bytes=' || coalesce(sum(length(sequence1) + length(qual1)), 0) FROM read_parquet('$STAGED');
SELECT 'r2_bytes=' || coalesce(sum(length(sequence2) + length(qual2)), 0) FROM read_parquet('$STAGED') WHERE sequence2 IS NOT NULL;
SQL

log "Staged. Metadata:"
sed 's/^/[fetch]   /' "$STAGED_META" >&2
log "Parquet: $(ls -lh "$STAGED" | awk '{print $5}')  ($STAGED)"
log "Run the benchmark with: bash $SCRIPT_DIR/run_bench.sh"
