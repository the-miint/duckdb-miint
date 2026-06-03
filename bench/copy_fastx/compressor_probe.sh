#!/usr/bin/env bash
#
# Isolate WHY `COPY ... FORMAT fastq` with gzip is slow, and quantify the ceiling
# achievable with compressors we already have (htslib bgzf) vs the current path.
#
# Produces a real uncompressed FASTQ sample from the staged ENA data, then times
# every available compressor at level 6 (DuckDB's default), single- vs
# multi-threaded, plus an isolated zlib mem_level=1-vs-8 test (DuckDB inits gzip
# with mem_level=1) and DuckDB's own COPY on the same sample.
#
# `bgzip` here IS htslib's bgzf CLI (the library we embed). pigz/libdeflate are
# only ceiling references and are NOT proposed as dependencies.
#
# Re-run after changing the writer to confirm the in-process path matches.
#
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DUCKDB="${MIINT_DUCKDB:-$REPO_ROOT/build/release/duckdb}"
PQ="${MIINT_BENCH_DATA_DIR:-$SCRIPT_DIR/data}/staged.parquet"
S="${OUTDIR:-${TMPDIR:-/tmp}/miint_cmp}"
N="${SAMPLE_READS:-2000000}"
mkdir -p "$S"
[[ -f "$PQ" ]] || { echo "stage data first: bash $SCRIPT_DIR/fetch_data.sh" >&2; exit 1; }

export TIMEFORMAT='%R'
run() { local name="$1"; shift; local t; t=$( { time "$@" >/dev/null 2>/dev/null; } 2>&1 ); \
  awk -v n="$name" -v b="$BYTES" -v t="$t" 'BEGIN{printf "%-26s %7.2fs  %8.1f MB/s in\n", n, t, (b/1048576.0)/t}'; }

if [[ ! -f "$S/sample.fq" ]]; then
  echo "Producing ${N}-read uncompressed FASTQ sample ..."
  "$DUCKDB" -batch :memory: \
    "COPY (SELECT read_id, comment, sequence1, qual1 FROM read_parquet('$PQ') LIMIT $N) TO '$S/sample.fq' (FORMAT FASTQ);" >/dev/null
fi
BYTES=$(stat -c %s "$S/sample.fq"); echo "sample: $BYTES bytes (~$((BYTES/1048576)) MB)"

NPROC="$(nproc)"
echo "=== single-threaded (level 6) ==="
command -v gzip            >/dev/null && run "gzip -6"             bash -c "gzip -6 -c '$S/sample.fq'"
command -v libdeflate-gzip >/dev/null && run "libdeflate-gzip -6 (ref)" bash -c "libdeflate-gzip -6 -c '$S/sample.fq'"
command -v bgzip           >/dev/null && run "bgzip -l6 -@1 (htslib)" bash -c "bgzip -l6 -c -@1 '$S/sample.fq'"
command -v pigz            >/dev/null && run "pigz -6 -p1 (ref)"    bash -c "pigz -6 -p1 -c '$S/sample.fq'"
echo "=== multi-threaded (level 6, $NPROC threads) ==="
command -v bgzip           >/dev/null && run "bgzip -l6 -@$NPROC (htslib)" bash -c "bgzip -l6 -c -@$NPROC '$S/sample.fq'"
command -v pigz            >/dev/null && run "pigz -6 -p$NPROC (ref)"  bash -c "pigz -6 -p$NPROC -c '$S/sample.fq'"

echo "=== DuckDB COPY gzip on the same sample (current path) ==="
t=$( { time "$DUCKDB" -batch :memory: \
  "COPY (SELECT read_id, comment, sequence1, qual1 FROM read_parquet('$PQ') LIMIT $N) TO '$S/ddb.fq.gz' (FORMAT FASTQ, COMPRESSION gzip);" >/dev/null 2>&1; } 2>&1 )
awk -v b="$BYTES" -v t="$t" 'BEGIN{printf "DuckDB COPY gzip           %7.2fs  %8.1f MB/s in\n", t, (b/1048576.0)/t}'

echo "=== compressed size / ratio (level 6) ==="
for f in "$S"/*.gz; do [[ -f "$f" ]] && printf '  %-14s %s bytes\n' "$(basename "$f")" "$(stat -c %s "$f")"; done

echo "=== isolate zlib mem_level (DuckDB inits gzip with mem_level=1; default is 8) ==="
head -c $((150*1024*1024)) "$S/sample.fq" > "$S/slice.fq"
python3 - "$S/slice.fq" <<'PY'
import sys, zlib, time
data = open(sys.argv[1], 'rb').read(); mb = len(data)/1048576
for ml in (1, 8):
    c = zlib.compressobj(6, zlib.DEFLATED, -15, ml)
    t0 = time.perf_counter(); out = c.compress(data) + c.flush(); dt = time.perf_counter() - t0
    print(f"  memLevel={ml}: {dt:6.2f}s  {mb/dt:6.1f} MB/s in   out={len(out)} bytes")
PY
