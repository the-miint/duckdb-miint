#!/usr/bin/env bash
#
# Benchmark the COPY ... FORMAT {fastq,fasta} WRITE path, compressed and
# uncompressed, on real microbiome data fetched by fetch_data.sh.
#
# Design goals:
#   * Measure ONLY the COPY statement, not FASTQ decode. We stage the data into a
#     Parquet file once and load it into an in-memory temp table per run, then
#     time just the COPY via DuckDB's `.timer`.
#   * Be re-runnable and append-only: every row is tagged with the git commit and
#     thread count, so a baseline collected today can be compared against a run
#     after the writer is parallelized. Threads=1 vs threads=N showing NO
#     difference today is the expected baseline (the writer ignores parallelism);
#     a future run should show N pulling ahead.
#   * Touch no files destructively. Outputs go to unique paths under OUTDIR and
#     are truncated (not removed) after their size is recorded, to bound disk use
#     without using `rm` (a project hard rule).
#
# Usage:
#   bash run_bench.sh                 # stage (once) + run the default matrix
#   FORMATS="fastq" COMPS="gzip" THREADS="1 8" REPEATS=1 bash run_bench.sh
#
# Env knobs (all optional):
#   MIINT_DUCKDB     duckdb shell with the miint extension (default build/release/duckdb)
#   MIINT_EXT        if set, `LOAD '<path>'` is issued (for a vanilla duckdb shell)
#   MIINT_BENCH_DATA_DIR  data dir written by fetch_data.sh (default ./data)
#   OUTDIR           where COPY outputs are written (default $TMPDIR/miint_copybench_out)
#   RESULTS_CSV      results file (default ./results/bench_<commit>_<ts>.csv)
#   FORMATS          space list subset of {fastq fasta}   (default both)
#   LAYOUTS          space list subset of {single interleave split} (default "single split")
#   COMPS            space list subset of {none gzip}      (default both)
#   ORDERS           space list subset of {reorder preserve} (default both)
#   THREADS          space list of thread counts           (default "1 <nproc>")
#   REPEATS          repeats per cell, fastest is reported  (default 3)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

DUCKDB="${MIINT_DUCKDB:-$REPO_ROOT/build/release/duckdb}"
DATA_DIR="${MIINT_BENCH_DATA_DIR:-$SCRIPT_DIR/data}"
STAGED="$DATA_DIR/staged.parquet"
STAGED_META="$DATA_DIR/staged.meta"
OUTDIR="${OUTDIR:-${TMPDIR:-/tmp}/miint_copybench_out}"

NPROC="$(nproc 2>/dev/null || echo 4)"
FORMATS="${FORMATS:-fastq fasta}"
LAYOUTS="${LAYOUTS:-single split}"
COMPS="${COMPS:-none gzip}"
ORDERS="${ORDERS:-reorder preserve}"
THREADS="${THREADS:-1 $NPROC}"
REPEATS="${REPEATS:-3}"

COMMIT="$(git -C "$REPO_ROOT" rev-parse --short HEAD 2>/dev/null || echo nogit)"
DIRTY=""
git -C "$REPO_ROOT" diff --quiet 2>/dev/null || DIRTY="+dirty"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
RESULTS_CSV="${RESULTS_CSV:-$SCRIPT_DIR/results/bench_${COMMIT}${DIRTY}_${RUN_ID}.csv}"

mkdir -p "$OUTDIR" "$(dirname "$RESULTS_CSV")"

log() { printf '[bench] %s\n' "$*" >&2; }
die() { printf '[bench] ERROR: %s\n' "$*" >&2; exit 1; }

[[ -x "$DUCKDB" ]] || die "duckdb shell not found/executable: $DUCKDB (build it, or set MIINT_DUCKDB)"

# LOAD line for a vanilla shell; empty for the preloaded project shell.
LOAD_SQL=""
[[ -n "${MIINT_EXT:-}" ]] && LOAD_SQL="LOAD '${MIINT_EXT}';"

# Run a SQL script (stdin) through duckdb in-memory; return combined output.
run_sql() { "$DUCKDB" -batch :memory: 2>&1; }

# ---------------------------------------------------------------------------
# Stage once: ENA run -> Parquet (delegated to fetch_data.sh, idempotent).
# ---------------------------------------------------------------------------
stage() {
	if [[ -f "$STAGED" && -f "$STAGED_META" ]]; then
		log "Staged parquet already present: $STAGED"
		return
	fi
	log "No staged parquet yet; running fetch_data.sh (one-time ENA stage) ..."
	bash "$SCRIPT_DIR/fetch_data.sh" || die "fetch_data.sh failed"
	[[ -f "$STAGED" && -f "$STAGED_META" ]] || die "staging did not produce $STAGED"
}

ROWS=0
read_meta() { ROWS="$(awk -F= '/^rows=/{print $2}' "$STAGED_META")"; }

# ---------------------------------------------------------------------------
# Build the SELECT, COPY options, file extension and output spec for a cell.
# ---------------------------------------------------------------------------
make_select() {  # $1=format $2=layout -> echoes SELECT list + FROM
	local fmt="$1" layout="$2"
	case "$layout" in
		single)
			if [[ "$fmt" == fastq ]]; then echo "SELECT read_id, comment, sequence1, qual1 FROM src";
			else echo "SELECT read_id, comment, sequence1 FROM src"; fi ;;
		interleave|split)
			if [[ "$fmt" == fastq ]]; then echo "SELECT read_id, comment, sequence1, sequence2, qual1, qual2 FROM src WHERE sequence2 IS NOT NULL";
			else echo "SELECT read_id, comment, sequence1, sequence2 FROM src WHERE sequence2 IS NOT NULL"; fi ;;
	esac
}

ext_for() {  # $1=format $2=comp -> file extension
	local fmt="$1" comp="$2" base
	[[ "$fmt" == fastq ]] && base="fq" || base="fa"
	[[ "$comp" == gzip ]] && echo "$base.gz" || echo "$base"
}

copy_opts() {  # $1=format $2=layout $3=comp -> COPY options paren content
	local fmt="$1" layout="$2" comp="$3" opts
	[[ "$fmt" == fastq ]] && opts="FORMAT FASTQ" || opts="FORMAT FASTA"
	[[ "$layout" == interleave ]] && opts="$opts, INTERLEAVE true"
	[[ "$comp" == gzip ]] && opts="$opts, COMPRESSION gzip"
	echo "$opts"
}

# Sum the byte size(s) of a run's output, then truncate to reclaim disk (no rm).
measure_output() {  # $1=layout $2=basepath $3=ext -> echoes total bytes
	local layout="$1" base="$2" ext="$3" total=0 f
	local files=()
	if [[ "$layout" == split ]]; then
		files=("$base.R1.$ext" "$base.R2.$ext")
	else
		files=("$base.$ext")
	fi
	for f in "${files[@]}"; do
		[[ -f "$f" ]] || { echo "MISSING:$f" >&2; echo 0; return; }
		total=$(( total + $(stat -c '%s' "$f") ))
		: > "$f"   # truncate, do not delete
	done
	echo "$total"
}

# ---------------------------------------------------------------------------
# Main matrix
# ---------------------------------------------------------------------------
stage
read_meta
log "rows=$ROWS  commit=${COMMIT}${DIRTY}  nproc=$NPROC"
log "matrix: FORMATS=[$FORMATS] LAYOUTS=[$LAYOUTS] COMPS=[$COMPS] ORDERS=[$ORDERS] THREADS=[$THREADS] REPEATS=$REPEATS"
log "results -> $RESULTS_CSV"

printf 'timestamp,commit,nproc,format,layout,compression,order,threads,repeat,rows,output_bytes,wall_s,output_MBps,records_per_s\n' > "$RESULTS_CSV"

counter=0
for fmt in $FORMATS; do
  for layout in $LAYOUTS; do
    sel="$(make_select "$fmt" "$layout")"
    for comp in $COMPS; do
      ext="$(ext_for "$fmt" "$comp")"
      opts="$(copy_opts "$fmt" "$layout" "$comp")"
      for order in $ORDERS; do
        # 'reorder' lets DuckDB drop insertion order (the knob that gates a
        # parallel copy sink); 'preserve' keeps it (the correctness control).
        if [[ "$order" == reorder ]]; then pio="false"; else pio="true"; fi
        for th in $THREADS; do
          for rep in $(seq 1 "$REPEATS"); do
            counter=$((counter+1))
            base="$OUTDIR/r_${fmt}_${layout}_${comp}_${order}_t${th}_${rep}_${counter}"
            if [[ "$layout" == split ]]; then
              outpath="$base.{ORIENTATION}.$ext"
            else
              outpath="$base.$ext"
            fi

            sql="$LOAD_SQL
PRAGMA threads=$th;
SET preserve_insertion_order=$pio;
CREATE TEMP TABLE src AS SELECT * FROM read_parquet('$STAGED');
.timer on
COPY ($sel) TO '$outpath' ($opts);"

            out="$(printf '%s\n' "$sql" | run_sql)" || {
              printf '%s\n' "$out" | sed 's/^/[bench]   /' >&2
              die "COPY failed for $fmt/$layout/$comp/$order/t$th"
            }
            wall="$(printf '%s\n' "$out" | awk '/Run Time/{t=$5} END{print t}')"
            [[ -n "$wall" ]] || { printf '%s\n' "$out" | sed 's/^/[bench]   /' >&2; die "could not parse timer output"; }

            obytes="$(measure_output "$layout" "$base" "$ext")"
            mbps="$(awk -v b="$obytes" -v t="$wall" 'BEGIN{ if(t>0) printf "%.2f", (b/1048576.0)/t; else print "NA" }')"
            rps="$(awk -v r="$ROWS" -v t="$wall" 'BEGIN{ if(t>0) printf "%.0f", r/t; else print "NA" }')"

            printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
              "$(date +%FT%T)" "${COMMIT}${DIRTY}" "$NPROC" "$fmt" "$layout" "$comp" "$order" "$th" "$rep" \
              "$ROWS" "$obytes" "$wall" "$mbps" "$rps" >> "$RESULTS_CSV"
            log "$fmt/$layout/$comp/$order t=$th rep=$rep : ${wall}s  ${mbps} MB/s  out=${obytes}B"
          done
        done
      done
    done
  done
done

# ---------------------------------------------------------------------------
# Summary: fastest repeat per (format,layout,compression,order,threads)
# ---------------------------------------------------------------------------
log "Summary (fastest of $REPEATS repeats per cell):"
"$DUCKDB" -batch -box :memory: 2>/dev/null <<SQL || true
SELECT format, layout, compression, "order", threads,
       min(wall_s) AS best_s,
       max(output_MBps) AS best_MBps,
       max(records_per_s) AS best_rec_s
FROM read_csv('$RESULTS_CSV')
GROUP BY ALL
ORDER BY format, layout, compression, "order", threads;
SQL

log "Done. CSV: $RESULTS_CSV"
log "Compare two CSVs later with, e.g.:"
log "  $DUCKDB -box :memory: \"SELECT * FROM read_csv('$SCRIPT_DIR/results/*.csv')\""
