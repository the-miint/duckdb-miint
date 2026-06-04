# COPY ... FORMAT {fastq,fasta} write-path benchmark

Empirical harness to measure the FASTQ/FASTA **write** (`COPY ... TO`) throughput,
compressed and uncompressed, on real microbiome data — so we can establish a
baseline now and re-run it to verify the gains from any writer changes later.

These scripts only *collect* data. They do not modify the extension.

## Prerequisites

- The miint extension built: `./build/release/duckdb` (the preloaded shell).
  If you only have the loadable extension, point a vanilla shell at it with
  `MIINT_DUCKDB=/path/to/duckdb MIINT_EXT=/path/to/miint.duckdb_extension`.
- Network access to EBI/ENA (the stage step streams over `httpfs`).

## 1. Fetch / stage data (once)

```bash
bash bench/copy_fastx/fetch_data.sh
```

- Uses the extension's own `read_ena_sequences` (over `httpfs`) to stream a real
  paired Illumina WGS metagenome run straight into `data/staged.parquet` — the
  same code path the project ships, no hand-rolled downloads.
- Default run is **ERR10677119** (PRJEB11419, ~7.8M read pairs). Override with
  `ACCESSION=ERRxxxxxx`, or cap reads with `MAX_SEQUENCES=N` for fast iteration.
- Idempotent: re-running skips staging if `data/staged.parquet` exists.

Outputs land in `bench/copy_fastx/data/` (git-ignored): `staged.parquet` plus a
`staged.meta` with the row counts and per-read payload sizes.

## 2. Run the benchmark

```bash
bash bench/copy_fastx/run_bench.sh
```

What it does:
1. **Stages once** (via `fetch_data.sh`): streams the ENA run into
   `data/staged.parquet` so network + FASTQ decode are excluded from write timing.
2. For each matrix cell it starts a fresh in-memory DuckDB, loads the parquet
   into a temp table, sets `PRAGMA threads` and `preserve_insertion_order`,
   then times **only** the `COPY` via `.timer`.
3. Appends one row per repeat to `results/bench_<commit>_<ts>.csv`.

### Matrix (all overridable via env)

| dimension   | default            | why it's here |
|-------------|--------------------|---------------|
| `FORMATS`   | `fastq fasta`      | both share the writer |
| `LAYOUTS`   | `single split`     | `single` = single-end one file; `split` = paired `{ORIENTATION}` R1/R2 (the order-sensitive case) |
| `COMPS`     | `none gzip`        | isolates encoder cost vs gzip cost |
| `ORDERS`    | `reorder preserve` | `reorder` = `preserve_insertion_order=false`, the knob that *gates a parallel copy sink*; `preserve` is the correctness control |
| `THREADS`   | `1 <nproc>`        | 1 vs many |
| `REPEATS`   | `3`                | fastest repeat is reported |

Add `interleave` to `LAYOUTS` for the paired-interleaved path.

### Reading the baseline

On the current code the writer never opts into parallelism, so **`threads=1`
and `threads=N` should produce nearly identical times** — that flat line *is*
the finding. The `none` vs `gzip` gap shows how much of the wall time is gzip
(expected to dominate for `.gz`). After a change, re-run and compare:

```bash
./build/release/duckdb -box :memory: \
  "SELECT format, layout, compression, \"order\", threads,
          min(wall_s) best_s, max(output_MBps) best_MBps
   FROM read_csv('bench/copy_fastx/results/*.csv')
   GROUP BY ALL ORDER BY 1,2,3,4,5"
```

The `commit` column tags every row, so before/after rows are unambiguous even
in the same CSV glob.

## Notes / caveats

- COPY outputs go to `OUTDIR` (default `$TMPDIR/miint_copybench_out`) with unique
  names and are **truncated** (not deleted) after their size is recorded, to
  bound disk use without using `rm`. You can clear that dir yourself.
- `output_MBps` is bytes that actually hit disk / wall. For `gzip` rows that is
  *compressed* throughput; compare `none` vs `gzip` at equal `(format,layout)` to
  see the compression tax.
- This harness covers the FASTX writer only. SAM/BAM is a separate write path
  (needs aligned input) and is analyzed in the parent discussion, not measured
  here.
