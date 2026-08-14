# Reference based alignment methods

MIINT exposes a variety of methods for high throughput alignment. In general, MIINT attempt to provide a similar experience for input and output while honoring the often times comprehensive suite of options these tools provide.

## Table of contents

- [bowtie2](#bowtie2) - Read alignment with bowtie2.
  - [Index creation](#index-creation-with-bowtie2) - Creating an index for use later.
  - [Single index alignment](#single-index-alignment-with-bowtie2) - Aligning to a single index.
  - [Sharded alignment](#sharded-alignment-with-bowtie2) - Aligning to many indices.
- [minimap2](#minimap2) - Read alignment with minimap2.
  - [Index creation](#index-creation-with-minimap2) - Creating an index for use later.
  - [Single index alignment](#single-index-alignment-with-minimap2) - Aligning to a single index.
  - [Sharded alignment](#sharded-alignment-with-minimap2) - Aligning to many indices.
- [SortMeRNA](#sortmerna) - Read alignment with SortMeRNA.
  - [Alignment with SAM output] - Alignment with SAM output like minimap2 and bowtie2.
  - [Alignment with rRNA output] - Aligning with the native rRNA output.

### bowtie2

MIINT supports efficient interaction with bowtie2 through the [GPL-boundary](internals/embedded-tools.md). 

#### Index creation with bowtie2

Build and save a bowtie2 index to disk. This is the bowtie2 analogue of `save_minimap2_index`. The index is produced by the gpl-boundary daemon's bundled `bowtie2-build`, so this function **requires gpl-boundary** (install with `SELECT install_gpl_boundary();`).

**Use case:** Build per-shard bowtie2 indexes once, then align many query sets against them with `align_bowtie2_sharded`. `align_bowtie2_sharded` expects each shard's index at `<shard_directory>/<shard_name>/index.*.bt2`, so build each shard with `output_path = '<shard_directory>/<shard_name>/index'`.

**Parameters:**
- `subject_table` (VARCHAR): Name of table or view containing subject/reference sequences. Must have `read_fastx`-compatible schema. Cannot contain paired-end data (sequence2 must be NULL or absent). The `read_id` column may be `VARCHAR` or `BIGINT`; `BIGINT` ids are stringified before being written into the index (recovered subject names are always `VARCHAR`).
- `output_path` (VARCHAR): Basename **prefix** for the index. bowtie2-build writes multiple files (`<output_path>.1.bt2`, `.2.bt2`, `.3.bt2`, `.4.bt2`, `.rev.1.bt2`, `.rev.2.bt2`). The parent directory is created if it does not exist. Note this differs from `save_minimap2_index`, which writes a single `.mmi` file.
- `threads` (INTEGER, default: 1): Number of threads `bowtie2-build` may use. Must be >= 1.

**Output schema:**
- `success` (BOOLEAN): Always true if function completes successfully
- `index_path` (VARCHAR): The `output_path` prefix the index files were written under
- `num_subjects` (BIGINT): Number of subject sequences indexed

**Behavior:**
- Loads all subject sequences from the table (rejects an empty table and rejects paired subjects)
- Submits them to the gpl-boundary daemon's `bowtie2-build` tool, which writes the `.bt2` files at `output_path`
- Returns a single row with success status and metadata
- Unlike `align_bowtie2`/`align_bowtie2_sharded` (which need gpl-boundary >= 0.4.2 for the alignment-time `memory_mapped` control), index building works with any daemon that advertises `bowtie2-build`.

**Examples:**
```sql
-- Build a single bowtie2 index from references fetched from NCBI
CREATE TABLE refs AS
    SELECT * FROM read_ncbi_fasta(['NC_000913.3', 'NC_002695.2', 'NC_011751.1']);
SELECT * FROM save_bowtie2_index('refs', 'shards/ecoli/index');
-- Returns: true | shards/ecoli/index | 3

-- Build per-shard indexes for align_bowtie2_sharded (one subdir per shard)
CREATE TABLE shard_a AS SELECT * FROM read_fastx('shard_a.fna');
CREATE TABLE shard_b AS SELECT * FROM read_fastx('shard_b.fna');
SELECT * FROM save_bowtie2_index('shard_a', 'shards/shard_a/index');
SELECT * FROM save_bowtie2_index('shard_b', 'shards/shard_b/index');

-- Then align against the shard directory
SELECT * FROM align_bowtie2_sharded('queries',
    shard_directory := 'shards',
    read_to_shard := 'read_to_shard');
```

**Error handling:**
- Error if subject_table does not exist
- Error if subject_table contains paired-end data (sequence2 not NULL)
- Error if subject_table is empty (no sequences to index)
- Error if output_path cannot be written (permissions, disk space, invalid path)
- Error if tables lack required columns (read_id, sequence1)

**Performance notes:**
- Building an index takes time proportional to the total reference size
- Index files are typically 10-30% of the original FASTA size
- Loading a pre-built index is 10-30x faster than building it from sequences
- For workflows that align multiple query sets, build the index once and reuse it
- Index files are portable across systems with the same minimap2 version

**When to use pre-built indexes:**
- **Large reference databases**: WoLr2, RefSeq bacterial genomes, custom marker gene databases
- **Multiple metagenome samples**: When aligning different samples to the same reference database
- **Production pipelines**: Build indexes during setup, use them in production runs
- **Shared infrastructure**: Build indexes once, share them across users/jobs

**When NOT to use pre-built indexes:**
- **Small references**: For small reference sets (<10MB), on-the-fly indexing is fast enough
- **One-time alignments**: If aligning only once, the index-building time won't be recovered

#### Single index alignment with bowtie2

Align query sequences to subject sequences using Bowtie2. This function enables short-read alignment directly within SQL by reading sequences from DuckDB tables/views and returning alignments in the same format as `read_alignments`.

Bowtie2 runs out of process via the `gpl-boundary` daemon (an Apache-licensed process-isolation host shipped separately from this extension). Sequences cross the boundary as Arrow IPC; bowtie2's command line is never user-controlled.

**Requirements:**
- The `gpl-boundary` daemon must be installed. Easiest path: `SELECT install_gpl_boundary();`. The miint extension itself does not link bowtie2.

**Parameters:**
- `query_table` (VARCHAR): Name of table or view containing query sequences. Must have `read_fastx`-compatible schema (read_id, sequence1, optional sequence2/qual1/qual2). The `read_id` column may be `VARCHAR`, `BIGINT`, or `UUID` (see *Identifier-column types* below).
- `subject_table` (VARCHAR): Name of table or view containing subject/reference sequences. Must have `read_fastx`-compatible schema. Cannot contain paired-end data (sequence2 must be NULL or absent). The `read_id` column may be `VARCHAR`, `BIGINT`, or `UUID`.
- `preset` (VARCHAR, optional): Bowtie2 preset for alignment sensitivity
  - `'very-fast'`: Fastest, least sensitive
  - `'fast'`: Fast alignment
  - `'sensitive'`: More sensitive (slower)
  - `'very-sensitive'`: Most sensitive, slowest
- `local` (BOOLEAN, default: false): Use local alignment mode instead of end-to-end
  - `false` (default): End-to-end alignment (entire read must align)
  - `true`: Local alignment (soft-clipping allowed at ends)
- `threads` (INTEGER, default: 1): Number of threads for Bowtie2 alignment (daemon `nthreads`)
- `max_secondary` (INTEGER, default: 1): Maximum alignments to report per query (`-k`)
- `quiet` (BOOLEAN, default: true): Runs Bowtie2 with `--quiet`. Keep the default. miint never surfaces Bowtie2's stderr alignment statistics to SQL, so `quiet := false` has **no user-visible effect** — it only makes the daemon emit per-batch summaries that miint immediately drains and discards, adding overhead for no benefit.
- `memory_mapped` (BOOLEAN, default: true): Memory-map the bowtie2 FM-index (`--mm`) rather than reading it sequentially. The daemon's default is `true` (`--mm` on), which lets concurrent processes share the mmap'd index pages — a win for warm-cache local runs. Set `false` to force a sequential `fread` (faster on a cold network filesystem; this is the default in [`align_bowtie2_sharded`](#sharded-alignment-with-bowtie2)). Requires gpl-boundary ≥ 0.4.2 — older daemons are rejected at query start, since they would otherwise silently ignore this field.

The full bowtie2-align typed parameter set is also exposed (one-to-one with the daemon's `--describe`; integer ranges enforced at bind time):

| Parameter | Type | Notes |
|---|---|---|
| `seed` | INTEGER | Random seed for reproducibility |
| `trim5`, `trim3` | INTEGER | Trim N bases from each end of the read |
| `match_bonus` | INTEGER | `--ma` |
| `mismatch_penalty` | INTEGER | `--mp` MAX (arg1). Pairs with `mismatch_penalty_min` |
| `mismatch_penalty_min` | INTEGER | `--mp` MIN (arg2). Set both equal for a symmetric penalty (e.g. `1,1`). Unset side defaults to bowtie2's compiled-in value (MAX=6, MIN=2); `MIN > MAX` is rejected |
| `n_penalty` | INTEGER | `--np` (ambiguous base penalty) |
| `read_gap_open`, `read_gap_extend` | INTEGER | `--rdg arg1,arg2` |
| `ref_gap_open`, `ref_gap_extend` | INTEGER | `--rfg arg1,arg2` |
| `score_min` | VARCHAR | Min-score function, e.g. `'L,-0.6,-0.6'` |
| `min_insert`, `max_insert` | INTEGER | Paired-end fragment-length bounds (`--minins`, `--maxins`) |
| `mate_orientation` | VARCHAR | One of `'fr'`, `'rf'`, `'ff'` |
| `no_mixed`, `no_discordant`, `dovetail`, `no_contain`, `no_overlap` | BOOLEAN | Paired-end behavior knobs |
| `nofw`, `norc` | BOOLEAN | Suppress forward / reverse-complement alignment |
| `seed_mismatches` | INTEGER | `-N`, must be 0 or 1 |
| `seed_length` | INTEGER | `-L`, must be 1–32 |
| `max_dp_failures` | INTEGER | `-D` |
| `max_seed_rounds` | INTEGER | `-R` |
| `report_all` | BOOLEAN | `-a`, report all alignments |
| `xeq` | BOOLEAN | Use `=`/`X` in CIGAR instead of `M` |
| `rg_id` | VARCHAR | `--rg-id` read group ID |
| `ignore_quals` | BOOLEAN | Treat all qualities as Phred 30 |
| `reorder` | BOOLEAN | Preserve input read order in output |
| `no_exact_upfront` | BOOLEAN | `--no-exact-upfront` (disable the exact end-to-end pass before multiseed) |
| `no_1mm_upfront` | BOOLEAN | `--no-1mm-upfront` (disable the 1-mismatch end-to-end pass before multiseed) |
| `deterministic_seeds` | BOOLEAN | `-d`/`--deterministic-seeds`: disable random subsampling of low-quality seed ranges for fully reproducible seed selection. **Coupling:** requires `report_all := true`, `no_exact_upfront := true`, and `no_1mm_upfront := true`, and is incompatible with `max_secondary`; otherwise the aligner rejects the run. |
| `lowseeds` | VARCHAR | `-l`/`--lowseeds`: discard low-quality seeds whose suffix-array range exceeds a threshold. Passed through verbatim, so a suffix is honored: bare = absolute count, `'%'` = /100, `'m'` = /1000, `'n'` = /1000000 of the reference |

Unknown parameters are rejected at bind time by DuckDB's binder (e.g. `presset := 'fast'` → `Invalid named parameter "presset"`).

**Output schema:**
Returns the same schema as `read_alignments` (21 columns):
- `read_id` (VARCHAR, BIGINT, or UUID — mirrors the query side): Query sequence identifier
- `flags` (USMALLINT): SAM alignment flags
- `reference` (VARCHAR, BIGINT, or UUID — mirrors the subject side): Subject sequence identifier
- `position` (BIGINT): 1-based start position on reference
- `stop_position` (BIGINT): 1-based stop position on reference
- `mapq` (UTINYINT): Mapping quality
- `cigar` (VARCHAR): CIGAR string
- `mate_reference` (VARCHAR, BIGINT, or UUID — mirrors the subject side): Mate reference (for paired-end)
- `mate_position` (BIGINT): Mate position (for paired-end)
- `template_length` (BIGINT): Template length (for paired-end)
- `tag_as` (BIGINT): Alignment score
- `tag_xs` (BIGINT): Suboptimal alignment score
- `tag_ys` (BIGINT): Mate alignment score
- `tag_xn` (BIGINT): Number of ambiguous bases
- `tag_xm` (BIGINT): Number of mismatches
- `tag_xo` (BIGINT): Number of gap opens
- `tag_xg` (BIGINT): Number of gap extensions
- `tag_nm` (BIGINT): Edit distance
- `tag_yt` (VARCHAR): Pair type (UU/CP/DP/UP)
- `tag_md` (VARCHAR): MD tag string
- `tag_sa` (VARCHAR): Supplementary alignment info

**Identifier-column types (`read_id`, `reference`, `mate_reference`):**
- The input columns may be `VARCHAR`, `BIGINT`, or `UUID`. Other numeric types (INTEGER, UBIGINT, HUGEINT, DOUBLE) are rejected at bind time with the message `Column '<name>' in table '<table>' must be VARCHAR, BIGINT, or UUID`.
- The output `read_id` column type mirrors the query side; `reference` and `mate_reference` mirror the subject side. Query and subject id types are independent (mixed BIGINT/VARCHAR is allowed).
- The daemon wire schema is always strings — BIGINT support is C++/DuckDB-only. On egress, the codec parses decimal strings back to `int64_t`; on ingress, `BIGINT` values are stringified by DuckDB's implicit `Value::GetValue<std::string>()` cast before crossing the Arrow IPC boundary.
- For `BIGINT`/`UUID` subjects, the SAM `=` mate-reference sentinel is resolved to the row's `reference` value before being emitted (the literal `=` has no BIGINT/UUID encoding); VARCHAR output preserves `=` verbatim, matching pre-existing behavior.
- For `BIGINT`/`UUID` columns, the SAM `*` unmapped sentinel surfaces as SQL NULL. VARCHAR sentinels pass through verbatim.
- NULL `read_id` rows are rejected at row time with `NULL read_id or sequence1 in query table` — same contract as VARCHAR. The daemon's input schema requires both columns non-null.

**Behavior:**
- Subject sequences are loaded at bind time and indexed via the daemon's `bowtie2-build` tool (index must fit in RAM)
- Query sequences are streamed across the gpl-boundary as Arrow IPC batches
- Supports both single-end and paired-end query sequences (paired-end auto-detected from the query table's `sequence2`/`qual2` columns)
- Bowtie2 is optimized for short reads (Illumina); use `align_minimap2` for long reads

**Examples:**
```sql
-- Create tables with sequence data
CREATE TABLE subjects AS SELECT * FROM read_fastx('references.fasta');
CREATE TABLE queries AS SELECT * FROM read_fastx('reads.fastq');

-- Basic alignment
SELECT * FROM align_bowtie2('queries', 'subjects');

-- Get primary alignments only with high sensitivity
SELECT read_id, reference, position, mapq, cigar
FROM align_bowtie2('queries', 'subjects', preset='very-sensitive', max_secondary=1)
ORDER BY read_id;

-- Use local alignment mode for reads with adapters
SELECT * FROM align_bowtie2('queries', 'subjects', local=true);

-- Multi-threaded alignment
SELECT * FROM align_bowtie2('queries', 'subjects', threads=4);

-- Filter by mapping quality and identity
SELECT read_id, reference, position, alignment_seq_identity(cigar, tag_nm, tag_md) AS identity
FROM align_bowtie2('queries', 'subjects', max_secondary=1)
WHERE mapq >= 30
  AND alignment_seq_identity(cigar, tag_nm, tag_md) > 0.95;

-- Works with views
CREATE VIEW filtered_queries AS
    SELECT * FROM read_fastx('reads.fastq')
    WHERE LENGTH(sequence1) >= 50;

SELECT * FROM align_bowtie2('filtered_queries', 'subjects');

-- Paired-end alignment (table has sequence2 column)
CREATE TABLE paired_queries AS SELECT * FROM read_fastx('R1.fastq', sequence2='R2.fastq');
SELECT * FROM align_bowtie2('paired_queries', 'subjects');

-- Export alignments to SAM format
CREATE TABLE refs AS SELECT read_id AS name, LENGTH(sequence1) AS length FROM subjects;
COPY (
    SELECT * FROM align_bowtie2('queries', 'subjects', max_secondary=1)
) TO 'alignments.sam' (FORMAT SAM, REFERENCE_LENGTHS 'refs');
```

**Error handling:**
- Error if query_table or subject_table does not exist
- Error if subject_table contains paired-end data (sequence2 not NULL)
- Error if tables lack required columns (read_id, sequence1)
- Error if the `gpl-boundary` daemon is not installed (see `SELECT install_gpl_boundary();`)

**Performance notes:**
- Bowtie2 is optimized for short reads (typically <500bp); for long reads, use `align_minimap2`
- The `threads` parameter controls Bowtie2's internal parallelism
- Query sequences are streamed in batches of 1024 to limit memory usage
- Subject sequences must fit in memory for index building

**Comparison with align_minimap2:**
| Feature | align_bowtie2 | align_minimap2 |
|---------|---------------|----------------|
| Best for | Short reads (Illumina) | Long reads (ONT, PacBio) |
| Alignment mode | End-to-end or local | Various presets |
| Index type | FM-index | Minimizer index |
| Paired-end | Interleaved stdin | Native support |
| Per-subject mode | No | Yes |

#### Sharded alignment with bowtie2

Align query sequences against multiple pre-built Bowtie2 index shards. Each shard is a directory containing a Bowtie2 index (with prefix `index`), and a mapping table specifies which reads should be aligned against which shard. This is the Bowtie2 counterpart to `align_minimap2_sharded`, designed for the same large-scale sharded alignment workflows. Bowtie2 runs out of process via the `gpl-boundary` daemon.

The sharded path emits only mapped reads (the daemon is invoked with `--no-unal`); use `align_bowtie2` directly if you need to inspect unaligned records.

**Requirements:**
- The `gpl-boundary` daemon must be installed (see `SELECT install_gpl_boundary();`).

**Parameters:**
- `query_table` (VARCHAR): Name of table or view containing query sequences. Must have `read_fastx`-compatible schema (read_id, sequence1, optional sequence2/qual1/qual2). The `read_id` column may be `VARCHAR`, `BIGINT`, or `UUID` (see *Identifier-column types* below).
- `shard_directory` (VARCHAR, required): Path to directory containing shard subdirectories. Each shard's Bowtie2 index is expected at `<shard_directory>/<shard_name>/index` (i.e., files like `<shard_name>/index.1.bt2`, `<shard_name>/index.rev.1.bt2`, etc.)
- `read_to_shard` (VARCHAR, required): Name of table or view that maps reads to shards. Must have columns:
  - `read_id`: Read identifier. Must match the storage type of `query_table.read_id` exactly — VARCHAR with VARCHAR, BIGINT with BIGINT. Mismatched types are rejected at bind time.
  - `shard_name` (VARCHAR): Name of the shard this read should be aligned against
- `preset` (VARCHAR, optional): Bowtie2 sensitivity preset ('very-fast', 'fast', 'sensitive', 'very-sensitive')
- `local` (BOOLEAN, default: false): Use local alignment mode instead of end-to-end
- `max_secondary` (INTEGER, default: 1): Maximum alignments to report per query (`-k`)
- `max_threads_per_shard` (INTEGER, default: 1, range 1–64): bowtie2's internal `nthreads` *floor* for each per-shard daemon worker, and the lever that trades shard-level concurrency against per-shard threading. The number of shards processed concurrently is `ceil(SET threads / max_threads_per_shard)`. The default of **1** runs one shard per DuckDB worker (each bowtie2 single-threaded), maximizing shard-level concurrency — fastest for the typical many-shard workload on a cold network filesystem, where per-shard index load dominates. **Raise it** for workloads with few, large shards on many cores, where bowtie2's internal `-p` threading matters more than the number of concurrent shards (e.g. 2 shards on a 32-core box would otherwise use only 2 cores). Once every shard has been claimed, a surviving worker grows its `nthreads` to reclaim cores freed by finished workers (up to `SET threads`), so a lone large shard at the end of a run isn't pinned to this value while the rest of the box sits idle. Total CPU stays bounded by `SET threads`; alignments are unaffected (bowtie2 `-p` only partitions reads across threads).
- `memory_mapped` (BOOLEAN, default: **false** in sharded mode): How each shard's bowtie2 FM-index is loaded by the daemon. `false` (the sharded default) reads the index with a sequential `fread`; `true` memory-maps it (bowtie2 `--mm`). Sharded mode defaults to `false` because on a cold network filesystem the sequential read is ~3.7× faster than the random page-faults of `--mm`. (The non-sharded [`align_bowtie2`](#single-index-alignment-with-bowtie2) keeps the daemon's `--mm`-on default, which lets warm-cache local runs share mmap'd index pages across processes.) Requires gpl-boundary ≥ 0.4.2 — older daemons are rejected at query start, since they would otherwise silently ignore this field.
- `prefetch_ahead` (INTEGER, default: auto, range 0–4096): How many upcoming shards' index files to warm into the OS page cache (`POSIX_FADV_WILLNEED`) ahead of the shard-claim frontier, so a worker doesn't stall on a cold network-FS fault when it claims the next shard. The default (`auto`, when unset) warms one full wave ahead — `ceil(SET threads / max_threads_per_shard)` shards; `0` disables prefetch. A pure cache hint: alignment results are identical for any value.
- `include_shard_name` (BOOLEAN, default: false): When true, append a `shard_name` column to the output
- `progress` (BOOLEAN, default: false): Opt-in progress reporting. When true, emit clean, timestamped per-shard lines to **stderr** (`shard i/N 'name': index loaded` → `done - R reads, A alignments (T s)`). Pure side channel — results are byte-identical to the default, which emits nothing, so programmatic callers are unaffected unless they pass `progress := true`. (Distinct from `MIINT_BT2_TELEMETRY`, which emits machine-readable per-batch TSV.)
- `quiet` (BOOLEAN, default: true): Runs Bowtie2 with `--quiet`. Keep the default — miint never surfaces Bowtie2's stderr statistics to SQL, so `quiet := false` has no user-visible effect and only adds overhead (per-batch summaries that miint drains and discards).
- `threads` (INTEGER): Ignored in sharded mode. Use DuckDB's `SET threads=N` to control cross-shard parallelism and `max_threads_per_shard` for per-shard bowtie2 threading. A warning is printed at bind if `threads != 1` is passed directly to this function.

The full bowtie2-align typed parameter set listed under [`align_bowtie2`](#single-index-alignment-with-bowtie2) above is also available here (same names, same bind-time validation): `seed`, `trim5`, `trim3`, `match_bonus`, `mismatch_penalty`, `mismatch_penalty_min`, `n_penalty`, `read_gap_open`, `read_gap_extend`, `ref_gap_open`, `ref_gap_extend`, `score_min`, `min_insert`, `max_insert`, `mate_orientation`, `no_mixed`, `no_discordant`, `dovetail`, `no_contain`, `no_overlap`, `nofw`, `norc`, `seed_mismatches`, `seed_length`, `max_dp_failures`, `max_seed_rounds`, `report_all`, `xeq`, `rg_id`, `ignore_quals`, `reorder`, `no_exact_upfront`, `no_1mm_upfront`, `deterministic_seeds`, `lowseeds`.

The `no_unal` knob is intentionally not exposed: the sharded path always emits only mapped reads (matches the pre-migration `FilterMappedOnly` contract). Use `align_bowtie2` directly if you need to inspect unaligned records.

**Output schema:**
Returns the same 21-column schema as `align_bowtie2` and `read_alignments`.

**Identifier-column types (`read_id`, `reference`, `mate_reference`):**
- The query side (`query_table.read_id` and `read_to_shard.read_id`) may be `VARCHAR`, `BIGINT`, or `UUID`. Both must share the same type — `ValidateReadToShardSchema` enforces strict equality so the underlying JOIN never relies on implicit casts. Other numeric types (INTEGER, UBIGINT, HUGEINT, DOUBLE) are rejected at bind time.
- The output `read_id` column mirrors the query side.
- The output `reference` and `mate_reference` columns are **always VARCHAR** in sharded mode. Sharded alignment always loads prebuilt Bowtie2 indexes, and reference names inside those indexes are opaque bytes — the same contract as `align_minimap2_sharded`. Cast in SQL (`CAST(reference AS BIGINT)`) if you need integer semantics downstream.
- For `BIGINT`/`UUID` `read_id`, NULL input rows that have no entry in `read_to_shard` are silently skipped (the JOIN filters them out before they reach the daemon).

**Behavior:**
- At bind time, reads the `read_to_shard` table to discover shards; per-shard index existence is verified at InitGlobal (not Bind) so the planner doesn't pay filesystem-stat cost on wide shard sets.
- Cross-shard parallelism is driven by `SET threads`: each DuckDB worker thread owns its own gpl-boundary daemon and claims shards atomically. Per-shard bowtie2 internal threading starts at `max_threads_per_shard` and ramps up for surviving shards once all shards are claimed (see `max_threads_per_shard` above) — this keeps the box saturated when a skewed shard-size distribution would otherwise leave a large shard running alone at the end of a run.
- **Memory budget:** with the mm-off default, each concurrently-processed shard holds its FM-index resident in its daemon process (not shared mmap pages). Budget at least **2× the on-disk index size per concurrent shard** (≈2.5–3× for headroom), i.e. roughly `ceil(SET threads / max_threads_per_shard) × 2× <largest index>` of RAM. This memory lives in the gpl-boundary daemon **subprocesses**, so DuckDB's `memory_limit` does *not* bound it. If you hit memory pressure, raise `max_threads_per_shard` to lower the concurrent-shard count. (For reference: WOL3 shard indexes are ~0.75–0.85 GB on disk and ~2.5 GB resident per concurrent shard.)
- For each shard, only the reads assigned to that shard (via the `read_to_shard` mapping) are streamed across the gpl-boundary as Arrow IPC
- A read can appear in multiple shards and will be aligned against each
- Unmapped reads (flag 0x4) are filtered out of results (daemon `--no-unal`)
- Supports both single-end and paired-end query sequences
- Supports views for both `query_table` and `read_to_shard`
- **`query_table` must return the same rows when read more than once.** Each shard opens its own cursor over `query_table`, so a relation that is not stable across repeated reads — a view using `nextval()` / `random()` / `now()`, a view over a table being written concurrently, or a registered single-pass Arrow stream such as a `RecordBatchReader` — would deliver fewer reads to later shards. This function **fails with an error** rather than returning a partial result:

  ```
  align_bowtie2_sharded: shard 'shard_b' delivered 0 of 500 mapped reads. The query
  relation returned different rows when re-read for this shard, ...
  ```

  Materialize the relation first and pass that instead:

  ```sql
  CREATE TEMP TABLE q AS SELECT * FROM my_unstable_view;
  SELECT * FROM align_bowtie2_sharded('q', shard_directory := 'indexes/', read_to_shard := 'read_to_shard');
  ```

  This function deliberately does *not* buffer the reads for you: its per-shard cursors stream with bounded memory, and materializing a large corpus internally would break larger-than-memory input (TEMP tables can only spill to `temp_directory`, never into a persistent database file). A `read_to_shard` that lists reads absent from `query_table` is **not** affected — that stays supported, and does not trigger the error.

**Examples:**
```sql
-- Setup: Pre-build shard indexes (one-time, using bowtie2-build externally)
-- Expected directory structure:
--   indexes/shard_a/index.1.bt2, indexes/shard_a/index.rev.1.bt2, ...
--   indexes/shard_b/index.1.bt2, indexes/shard_b/index.rev.1.bt2, ...

-- Load query sequences
CREATE TABLE queries AS SELECT * FROM read_fastx('reads.fastq');

-- Create read-to-shard mapping
CREATE TABLE read_to_shard AS SELECT * FROM (VALUES
    ('read1', 'shard_a'),
    ('read2', 'shard_b'),
    ('read3', 'shard_a')
) AS t(read_id, shard_name);

-- Align reads against their assigned shards in parallel
SELECT * FROM align_bowtie2_sharded('queries',
    shard_directory := 'indexes/',
    read_to_shard := 'read_to_shard',
    max_secondary := 0);

-- With sensitivity preset and local alignment
SELECT read_id, reference, position, mapq
FROM align_bowtie2_sharded('queries',
    shard_directory := 'indexes/',
    read_to_shard := 'read_to_shard',
    preset := 'very-sensitive',
    local := true,
    max_secondary := 0)
WHERE mapq >= 30;

-- Paired-end alignment
CREATE TABLE paired_queries AS SELECT * FROM read_fastx('R1.fastq', sequence2='R2.fastq');
SELECT * FROM align_bowtie2_sharded('paired_queries',
    shard_directory := 'indexes/',
    read_to_shard := 'read_to_shard',
    max_secondary := 0);
```

**Error handling:**
- Error if `shard_directory` does not exist
- Error if any shard's Bowtie2 index files are missing (e.g., no `index.1.bt2` at `<shard_directory>/<shard_name>/index`)
- Error if `query_table` or `read_to_shard` table/view does not exist
- Error if `read_to_shard` table is missing `read_id` or `shard_name` columns
- Error if `read_to_shard` table contains NULL `shard_name` values
- Error if the `gpl-boundary` daemon is not installed (see `SELECT install_gpl_boundary();`)

**Performance notes:**
- Shards are processed sequentially; per-shard bowtie2 parallelism comes from `max_threads_per_shard`
- The `threads` parameter is ignored (a warning is emitted); use `max_threads_per_shard` instead
- Shards are sorted by read count (largest first) so the first shard's index stays hot in the daemon's worker pool

**Comparison of sharded vs non-sharded alignment functions:**
| Feature | `align_minimap2` / `align_bowtie2` | `align_minimap2_sharded` / `align_bowtie2_sharded` |
|---------|--------------------------------------|------------------------------------------------------|
| Index source | Build on-the-fly or single pre-built index | Multiple pre-built indexes (one per shard) |
| Read routing | All reads against one index | Reads routed to specific shards via mapping table |
| Parallelism | Single aligner thread(s) | One aligner per shard, concurrent |
| Use case | Single reference database | Sharded reference databases (e.g., from prior classification) |

### minimap2

MIINT embeds minimap2 providing a highly efficient interation.

#### Index creation with minimap2

Build and save a minimap2 index to disk for reuse. This provides 10-30x performance improvement when aligning multiple query sets against the same large reference database.

**Use case:** Build indexes once for large reference databases (e.g., WoLr2 phylogenetic markers, RefSeq genomes, custom OGU databases), then use them repeatedly with `align_minimap2(..., index_path='file.mmi')` instead of rebuilding the index each time.

**Parameters:**
- `subject_table` (VARCHAR): Name of table or view containing subject/reference sequences. Must have `read_fastx`-compatible schema. Cannot contain paired-end data (sequence2 must be NULL or absent). The `read_id` column may be `VARCHAR`, `BIGINT`, or `UUID`; `BIGINT` ids are stringified to decimal before being written into the index. When the index is later loaded via `align_minimap2(..., index_path=...)`, the recovered subject names are always `VARCHAR` (see *Identifier-column types* in `align_minimap2`).
- `output_path` (VARCHAR): Path where the index file (.mmi) will be saved
- `preset` (VARCHAR, default: 'sr'): Minimap2 preset (same options as `align_minimap2`)
- `k` (INTEGER, optional): K-mer size (overrides preset default if specified)
- `w` (INTEGER, optional): Minimizer window size (overrides preset default if specified)


**Output schema:**
- `success` (BOOLEAN): Always true if function completes successfully
- `index_path` (VARCHAR): Path to the created index file
- `num_subjects` (BIGINT): Number of subject sequences indexed

**Behavior:**
- Loads all subject sequences from the table
- Builds minimap2 index in memory with specified parameters
- Writes index to disk in .mmi format
- Index file stores k-mer size, window size, and preset configuration
- Returns a single row with success status and metadata

**Examples:**
```sql
-- Create reference table from WoLr2 phylogenetic markers
CREATE TABLE wolr2_refs AS SELECT * FROM read_fastx('WoLr2_db.fna');

-- Build index for short reads (default preset)
SELECT * FROM save_minimap2_index('wolr2_refs', 'wolr2.mmi');

-- Build index with specific preset
SELECT * FROM save_minimap2_index('wolr2_refs', 'wolr2_sr.mmi', preset='sr');
SELECT * FROM save_minimap2_index('wolr2_refs', 'wolr2_ont.mmi', preset='map-ont');

-- Build index with custom k-mer size
SELECT * FROM save_minimap2_index('wolr2_refs', 'wolr2_k15.mmi', k=15);

-- Check the result
SELECT success, index_path, num_subjects
FROM save_minimap2_index('wolr2_refs', 'my_index.mmi', preset='sr');
-- Returns: true | /path/to/my_index.mmi | 10575

-- Use the saved index for alignment (10-30x faster than rebuilding)
CREATE TABLE metagenomic_reads AS SELECT * FROM read_fastx('metagenome.fastq');
SELECT * FROM align_minimap2('metagenomic_reads', index_path='wolr2.mmi', max_secondary=0);

-- Build indexes for different sequencing technologies
CREATE TABLE refseq_genomes AS SELECT * FROM read_fastx('refseq/*.fna');
SELECT * FROM save_minimap2_index('refseq_genomes', 'refseq_sr.mmi', preset='sr');
SELECT * FROM save_minimap2_index('refseq_genomes', 'refseq_ont.mmi', preset='map-ont');
SELECT * FROM save_minimap2_index('refseq_genomes', 'refseq_pb.mmi', preset='map-pb');

-- Then use them as needed
SELECT * FROM align_minimap2('illumina_metagenome', index_path='refseq_sr.mmi');
SELECT * FROM align_minimap2('nanopore_metagenome', index_path='refseq_ont.mmi');
SELECT * FROM align_minimap2('pacbio_metagenome', index_path='refseq_pb.mmi');

-- Build index from bacterial genomes fetched from NCBI
CREATE TABLE ecoli_genomes AS
    SELECT * FROM read_ncbi_fasta(['NC_000913.3', 'NC_002695.2', 'NC_011751.1']);

SELECT * FROM save_minimap2_index('ecoli_genomes', 'ecoli_strains.mmi', preset='sr');

-- Works with views
CREATE VIEW marker_genes AS
    SELECT * FROM read_fastx('markers/*.fna')
    WHERE LENGTH(sequence1) >= 500;

SELECT * FROM save_minimap2_index('marker_genes', 'markers.mmi');
```

#### Single index alignment with minimap2

Align query sequences to subject sequences using minimap2. This function enables sequence alignment directly within SQL by reading sequences from DuckDB tables/views and returning alignments in the same format as `read_alignments`.

**Performance:** For large reference databases (e.g., human genome), use pre-built indexes via `index_path` for 10-30x faster alignment. Build indexes once with `save_minimap2_index`, then reuse them across multiple queries.

**Parameters:**
- `query_table` (VARCHAR): Name of table or view containing query sequences. Must have `read_fastx`-compatible schema (read_id, sequence1, optional sequence2/qual1/qual2). The `read_id` column may be `VARCHAR`, `BIGINT`, or `UUID` (see *Identifier-column types* below).
- `subject_table` (VARCHAR, optional): Name of table or view containing subject/reference sequences. Must have `read_fastx`-compatible schema. Cannot contain paired-end data (sequence2 must be NULL or absent). The `read_id` column may be `VARCHAR`, `BIGINT`, or `UUID`. **Either `subject_table` OR `index_path` must be provided, but not both.**
- `index_path` (VARCHAR, optional): Path to pre-built minimap2 index file (.mmi). Use `save_minimap2_index()` to create index files. **Either `subject_table` OR `index_path` must be provided, but not both.**
- `per_subject_database` (BOOLEAN, default: false): Build separate index for each subject sequence (only valid with `subject_table`, not with `index_path`)
  - `false` (default): Build single index from all subjects, align all queries once (efficient for many subjects)
  - `true`: Build index per subject, align all queries against each (useful for per-genome analysis)
- `preset` (VARCHAR, default: 'sr'): Minimap2 preset
  - `'sr'`: Short reads (Illumina), k=21
  - `'map-ont'`: Oxford Nanopore reads, k=15
  - `'map-pb'`: PacBio reads, k=19
  - Other minimap2 presets also supported
- `max_secondary` (INTEGER, default: 5): Maximum secondary alignments to report per query. Set to 0 for primary alignments only
- `k` (INTEGER, optional): K-mer size (overrides preset default if specified). **Warning:** Ignored when using `index_path` (k-mer size is baked into the pre-built index)
- `w` (INTEGER, optional): Minimizer window size (overrides preset default if specified). **Warning:** Ignored when using `index_path` (window size is baked into the pre-built index)
- `eqx` (BOOLEAN, default: true): Use =/X CIGAR operators instead of M
- `occ_filter` (optional): minimap2's `-f` high-occurrence minimizer filter. Accepts a bare number or the two-value `'INT1,INT2'` string. See *Dense reference sets* below.
- `include_unmapped` (BOOLEAN, default: false): Emit one row per query that produced no alignment, instead of no row at all. See *Unmapped queries* below.
- `min_chain_coverage` (FLOAT, default: 0.0 = disabled, range 0.0–1.0): Skip the expensive dynamic-programming alignment for any query whose best seed **chain** spans less than this fraction of the query. See *Chain-coverage pre-filter* below.
- `debug` (BOOLEAN, default: false): Emit timestamped, thread-tagged diagnostics to **stderr** — index construction (with RSS), per-thread startup, and per-batch alignment counts and timings. Lines from concurrent workers interleave, so redirect and sort by thread id when reading them. A side channel only; results are unchanged.

**Chain-coverage pre-filter (`min_chain_coverage`):**

A speed knob with a **lossy** edge, so it is off by default. Before running DP alignment, each chain's span is measured as `(qe - qs) / qlen`; if no chain reaches the threshold, the query is dropped and produces no alignment at all.

```sql
-- Only attempt full alignment where a chain already spans ≥70% of the query
SELECT * FROM align_minimap2('q', subject_table := 's', min_chain_coverage := 0.7);
```

Points that matter in practice:

- This is **chain-phase** coverage, not post-DP aligned-base coverage. A chain's span includes the gaps between its seeds, so it *overestimates* true coverage — a chain spanning 90% of the query may align far less. Filtering on it is therefore approximate by nature.
- It **discards queries**, it does not merely reorder work. A query filtered here is indistinguishable from one with no alignment, which is exactly the ambiguity `include_unmapped` addresses — set both together if you need to tell "filtered out" from "measured and distant". `include_unmapped` does account for queries dropped by this filter.
- The vendored implementation records **0.70** as the empirical ceiling for zero false negatives on HiFi data. Do not tune close to it: a 51-base query whose last 15 bases diverge has a matching prefix of ~0.70 of its length, but is discarded at a threshold of `0.68` and survives only at `0.65` — the chain span the filter actually measures is around 0.66, below what the prefix fraction suggests. Estimating the right threshold by hand from expected identity will discard real alignments.
- For paired-end input the filter is applied **per segment**, against that segment's own length.
- Values outside 0.0–1.0 are rejected at bind with `min_chain_coverage must be between 0.0 and 1.0`.

> This parameter is **not** a minimap2 command-line option. It is a local addition to the vendored minimap2 (`ext/minimap2`, commit *"Add min_chain_coverage pre-filter to skip DP for low-coverage chains"* on top of upstream 2.30), so it has no `-`-flag equivalent and will not be found in minimap2's own documentation.

**Dense reference sets (`occ_filter`):**

minimap2 defaults to `-f 1000,5000`: it discards minimizers that occur more often than that in the index. That is tuned for genome-scale references with real repeats, and it is **actively wrong for a set of near-identical homologous sequences** — rRNA panels, gene families, dereplicated marker sets — where nearly every minimizer is high-occurrence by construction and gets masked. The symptoms are that a verbatim substring of an indexed sequence may not recover a perfect self match (so reported identity is only a *lower bound*), and that many queries lose their alignment row entirely.

```sql
-- Disable the filter (minimap2's own -f 0)
SELECT * FROM align_minimap2('q', subject_table := 's', occ_filter := 0);

-- Raise the cap to a large explicit value
SELECT * FROM align_minimap2('q', subject_table := 's', occ_filter := 100000);

-- Set both values, as -f INT1,INT2 (mid_occ, max_occ)
SELECT * FROM align_minimap2('q', subject_table := 's', occ_filter := '1000,5000');
```

The semantics are minimap2's, including two that surprise people:

- **A value below 1 is a *fraction*, not a count.** `occ_filter := 0.0002` keeps the top 0.02% of minimizers masked, computed from the index's own distribution. `occ_filter := 0` is therefore "mask the top 0 fraction", i.e. **disabled** — not "a threshold of zero".
- **The preset's *second* value does the real work for short reads.** `sr` sets `mid_occ=1000` *and* `max_occ=5000`, and minimap2 re-chains using `max_occ` when the first pass finds only repetitive seeds. So occurrences between 1000 and 5000 are already rescued, and tightening the filter requires setting both values (`'1,1'`), not just the first.
- `occ_filter := 0` relies on a clamp to `max_mid_occ`, which the `lr:hq`, `map-hifi`, `map-ccs`, `lr:hqae` and `map-iclr` presets narrow to 500. Under those presets `occ_filter := 0` lands at 500 rather than disabling the filter; pass a large explicit value instead.

`align_minimap2_sharded` accepts `occ_filter` too — the threshold is per index, so it applies to each shard independently.

**Unmapped queries (`include_unmapped`):**

By default a query that produces no alignment yields **no row**, so the absence of a row conflates "genuinely distant from every subject" with "the aligner found no seed chain" (repeat masking, short query, low complexity). Consumers then have to reconstruct "not measured" with an anti-join and remember the distinction exists; writing `coalesce(identity, 0) < 0.90` reads "no chain found" as "definitively distant", which inverts the safe default.

With `include_unmapped := true`, every input query is represented:

```sql
SELECT read_id FROM align_minimap2('q', subject_table := 's', include_unmapped := true)
WHERE reference IS NULL;   -- queries that were measured and did not align
```

- `reference`, `position`, `stop_position`, `mapq`, `cigar` and the `tag_*` columns are **NULL** (not the SAM `'*'`/`0` text sentinels), so `IS NULL` is the test. The one exception is `tag_yt`, the pair type, which is known regardless of whether the read aligned and is still `'UU'`/`'UP'`.
- The SAM unmapped flag `0x4` is set, so `flags & 4 != 0` selects the same rows.
- Queries with an empty `sequence1` also get a row — minimap2 cannot be handed a zero-length query, but the query was still submitted and must be accounted for.
- For paired input the accounting is **per segment**: if R1 aligns and R2 does not, you get the R1 alignment plus one unmapped row for R2. Mate flags stay truthful — the unmapped row carries the paired and second-in-pair bits, *not* mate-unmapped, and it carries the mate's `mate_reference`/`mate_position` as SAM requires of a paired record whose mate is mapped.
- These rows can be written straight out: `COPY (SELECT * FROM align_minimap2(…, include_unmapped := true)) TO 'x.bam' (FORMAT BAM, REFERENCE_LENGTHS 'r')` works, because the writer accepts a NULL `reference`/`position`/`mapq`/`cigar` on a record flagged `0x4`. Note the round-trip asymmetry: NULL is written as the SAM `*` sentinel, so reading the file back with `read_alignments` returns `'*'` rather than NULL.
- **Not supported by `align_minimap2_sharded`** (the parameter is rejected), and **not supported with `per_subject_database`** (the combination is rejected at bind). Both re-align each query against a different subject set, so a query that finds no chain in one shard or against one subject routinely maps in another — a synthetic row there would assert "did not align" about a query that did. Correct support needs reconciliation across the whole subject set.

**Output schema:**
Returns the same schema as `read_alignments` (21 columns):
- `read_id` (VARCHAR, BIGINT, or UUID — mirrors `query_table.read_id`): Query sequence identifier
- `flags` (USMALLINT): SAM alignment flags
- `reference` (VARCHAR, BIGINT, or UUID — mirrors `subject_table.read_id`; always VARCHAR when using `index_path`): Subject sequence identifier
- `position` (BIGINT): 1-based start position on reference
- `stop_position` (BIGINT): 1-based stop position on reference
- `mapq` (UTINYINT): Mapping quality
- `cigar` (VARCHAR): CIGAR string (with =/X by default)
- `mate_reference` (VARCHAR, BIGINT, or UUID — same type as `reference`): Mate reference (for paired-end)
- `mate_position` (BIGINT): Mate position (for paired-end)
- `template_length` (BIGINT): Template length (for paired-end)
- `tag_as` (BIGINT): Alignment score
- `tag_xs` (BIGINT): Suboptimal alignment score
- `tag_ys` (BIGINT): Mate alignment score
- `tag_xn` (BIGINT): Number of ambiguous bases
- `tag_xm` (BIGINT): Number of mismatches
- `tag_xo` (BIGINT): Number of gap opens
- `tag_xg` (BIGINT): Number of gap extensions
- `tag_nm` (BIGINT): Edit distance
- `tag_yt` (VARCHAR): Pair type (UU/CP/DP/UP)
- `tag_md` (VARCHAR): MD tag string
- `tag_sa` (VARCHAR): Supplementary alignment info

**Identifier-column types (`read_id`, `reference`, `mate_reference`):**
- The input columns may be `VARCHAR`, `BIGINT`, or `UUID`. Other numeric types (INTEGER, UBIGINT, HUGEINT, DOUBLE) are rejected at bind time with the message `Column '<name>' in table '<table>' must be VARCHAR, BIGINT, or UUID`.
- The output `read_id` column type mirrors the query side. The output `reference` and `mate_reference` types mirror the subject side, **except when using `index_path`** — subject names stored in a `.mmi` file are opaque bytes, so both default to `VARCHAR` regardless of what the source table looked like when the index was built.
- `BIGINT` ids are stringified to decimal on the wire; `UUID` ids are stringified to their canonical 36-char lowercase form. Both are parsed back to the original storage type on output (a `UUID` round-trips case-insensitively — uppercase input comes back canonical lowercase).
- For non-VARCHAR subjects (`BIGINT`/`UUID`), the SAM `=` mate-reference sentinel (which has no `BIGINT`/`UUID` encoding) is resolved to the primary reference value before being emitted; downstream consumers see the resolved id, never `'='`.
- For `BIGINT`/`UUID` columns, NULL input values and the SAM `'*'` unmapped sentinel are both surfaced as SQL NULL in the output. VARCHAR sentinels (`'*'`, `'='`) pass through verbatim, preserving pre-existing behaviour.

**Behavior:**
- Subject sequences are loaded into memory at bind time (must fit in RAM)
- Query sequences are processed in batches for memory efficiency
- Supports both single-end and paired-end query sequences
- Uses minimap2's in-memory indexing for fast alignment
- Secondary alignments are controlled by `max_secondary` parameter

**Examples:**
```sql
-- Create tables with sequence data
CREATE TABLE subjects AS SELECT * FROM read_fastx('references.fasta');
CREATE TABLE queries AS SELECT * FROM read_fastx('reads.fastq');

-- Basic alignment using subject_table (builds index on-the-fly)
SELECT * FROM align_minimap2('queries', subject_table='subjects');

-- Get only primary alignments
SELECT read_id, reference, position, mapq, cigar
FROM align_minimap2('queries', subject_table='subjects', max_secondary=0)
ORDER BY read_id;

-- Use Oxford Nanopore preset for long reads
SELECT * FROM align_minimap2('queries', subject_table='subjects', preset='map-ont');

-- === Using Pre-built Indexes (Recommended for large references) ===

-- Step 1: Build and save index once (see save_minimap2_index for details)
SELECT * FROM save_minimap2_index('subjects', 'references.mmi', preset='sr');

-- Step 2: Use pre-built index for fast alignment (10-30x faster for large references)
SELECT * FROM align_minimap2('queries', index_path='references.mmi', max_secondary=0);

-- Align different query sets against the same index
CREATE TABLE queries_batch1 AS SELECT * FROM read_fastx('batch1.fastq');
CREATE TABLE queries_batch2 AS SELECT * FROM read_fastx('batch2.fastq');

SELECT * FROM align_minimap2('queries_batch1', index_path='references.mmi');
SELECT * FROM align_minimap2('queries_batch2', index_path='references.mmi');

-- Pre-built index with different preset
SELECT * FROM save_minimap2_index('subjects', 'references_ont.mmi', preset='map-ont');
SELECT * FROM align_minimap2('long_reads', index_path='references_ont.mmi');

-- === Advanced Usage ===

-- Filter by mapping quality and identity
SELECT read_id, reference, position, alignment_seq_identity(cigar, tag_nm, tag_md) AS identity
FROM align_minimap2('queries', subject_table='subjects', max_secondary=0)
WHERE mapq >= 30
  AND alignment_seq_identity(cigar, tag_nm, tag_md) > 0.95;

-- Works with views
CREATE VIEW filtered_queries AS
    SELECT * FROM read_fastx('reads.fastq')
    WHERE LENGTH(sequence1) >= 50;

SELECT * FROM align_minimap2('filtered_queries', subject_table='subjects', max_secondary=0);

-- Align all queries against each subject separately (per-subject mode)
SELECT reference, COUNT(*) AS aligned_reads
FROM align_minimap2('queries', subject_table='subjects', per_subject_database=true, max_secondary=0)
GROUP BY reference;

-- Export alignments to SAM format
CREATE TABLE refs AS SELECT read_id AS name, LENGTH(sequence1) AS length FROM subjects;
COPY (
    SELECT * FROM align_minimap2('queries', subject_table='subjects', max_secondary=0)
) TO 'alignments.sam' (FORMAT SAM, REFERENCE_LENGTHS 'refs');

-- Calculate coverage per reference
-- position/stop_position are 1-based half-open, so the span is (stop - position) with no + 1.
SELECT reference,
       compress_intervals(position, stop_position) AS coverage_regions,
       SUM(stop_position - position) AS total_aligned_bases
FROM align_minimap2('queries', subject_table='subjects', max_secondary=0)
GROUP BY reference;

-- Paired-end alignment
CREATE TABLE paired_queries AS SELECT * FROM read_fastx('R1.fastq', sequence2='R2.fastq');
SELECT * FROM align_minimap2('paired_queries', subject_table='subjects', max_secondary=0);
```

**Error handling:**
- Error if query_table does not exist
- Error if neither `subject_table` nor `index_path` is provided
- Error if both `subject_table` and `index_path` are provided
- Error if `index_path` file does not exist or is not a valid minimap2 index
- Error if subject_table contains paired-end data (sequence2 not NULL)
- Error if tables lack required columns (read_id, sequence1)
- Error if preset is unknown to minimap2
- Error if `per_subject_database=true` is used with `index_path` (incompatible modes)
- Warning if `k` or `w` parameters are specified with `index_path` (these values are ignored)

**Performance notes:**
- **Pre-built indexes provide 10-30x faster alignment** for large reference databases
- Use `save_minimap2_index()` to build indexes once, then reuse them across multiple query sets
- Index files (.mmi) store k-mer size and window size, so `k` and `w` parameters are ignored when using `index_path`
- For large reference sets, the default mode (single index) is most efficient
- The `per_subject_database=true` mode rebuilds the index for each subject, which is slower but useful for specific analyses
- Query sequences are streamed in batches of 1024 to limit memory usage
- **`query_table` is read exactly once**, in both the default and `per_subject_database` modes. This makes single-pass relations safe: a registered Arrow relation — including a `RecordBatchReader` streamed from an external source such as Arrow Flight — can be passed by name with no intermediate file. Earlier versions paged `per_subject_database` through the relation with repeated `LIMIT/OFFSET` queries, which silently dropped reads for any relation that did not return identical rows on every read.
- Secondary alignments can significantly increase output size; use `max_secondary=0` for primary-only results

**Limitations:**
- Subject sequences must fit in memory (loaded at bind time for indexing when using `subject_table`)
- No support for reading sequences directly from files (use tables/views from `read_fastx`)

#### Sharded alignment with minimap2

Align query sequences against multiple pre-built minimap2 index shards in parallel. Each shard is a separate `.mmi` index file, and a mapping table specifies which reads should be aligned against which shard. This is designed for large-scale metagenomic workflows where the reference database is split across multiple shards and reads have been pre-assigned to shards (e.g., by a prior classification step).

**Parameters:**
- `query_table` (VARCHAR): Name of table or view containing query sequences. Must have `read_fastx`-compatible schema (read_id, sequence1, optional sequence2/qual1/qual2). The `read_id` column may be `VARCHAR`, `BIGINT`, or `UUID` (see *Identifier-column types* below).
- `shard_directory` (VARCHAR, required): Path to directory containing pre-built minimap2 index files. Each shard's index is expected at `<shard_directory>/<shard_name>.mmi`
- `read_to_shard` (VARCHAR, required): Name of table or view that maps reads to shards. Must have columns:
  - `read_id`: Read identifier. Must match the storage type of `query_table.read_id` exactly — VARCHAR with VARCHAR, BIGINT with BIGINT. Mismatched types are rejected at bind time.
  - `shard_name` (VARCHAR): Name of the shard this read should be aligned against
- `preset` (VARCHAR, default: 'sr'): Minimap2 preset ('sr', 'map-ont', 'map-pb', etc.)
- `max_secondary` (INTEGER, default: 5): Maximum secondary alignments per query. Set to 0 for primary only
- `eqx` (BOOLEAN, default: true): Use =/X CIGAR operators instead of M
- `progress` (BOOLEAN, default: false): Opt-in progress reporting. When true, emit clean, timestamped per-shard lines to **stderr** (`shard i/N 'name': index loaded, R reads` → `done - R reads, A alignments (T s)`). Pure side channel — results are byte-identical to the default, which emits nothing, so programmatic callers are unaffected unless they pass `progress := true`.
- `occ_filter` (optional): minimap2's `-f` high-occurrence minimizer filter, as for `align_minimap2` — see *Dense reference sets* above. The threshold is per index, so it applies to each shard independently.
- `min_chain_coverage` (FLOAT, default: 0.0 = disabled, range 0.0–1.0): Chain-coverage pre-filter, as for `align_minimap2` — see *Chain-coverage pre-filter* above.
- `debug` (BOOLEAN, default: false): Emit per-shard diagnostic detail to **stderr**, including memory checkpoints. A side channel only.

**`include_unmapped` is not available here.** A query that finds no seed chain in one shard routinely aligns in another, so a per-shard unmapped row would assert "did not align" about a query that did. The parameter is deliberately unregistered, so passing it is a binder error rather than a source of wrong rows; doing it correctly requires reconciling results across every shard a read was assigned to.

**Output schema:**
Returns the same 21-column schema as `align_minimap2` and `read_alignments`.

**Identifier-column types (`read_id`, `reference`, `mate_reference`):**
- The query side (`query_table.read_id` and `read_to_shard.read_id`) may be `VARCHAR`, `BIGINT`, or `UUID`. Both must share the same type — `ValidateReadToShardSchema` enforces strict equality so the underlying JOIN never relies on implicit casts. Other numeric types (INTEGER, UBIGINT, HUGEINT, DOUBLE) are rejected at bind time with the message `Column 'read_id' in table '<table>' must be VARCHAR, BIGINT, or UUID`.
- The output `read_id` column mirrors the query side.
- The output `reference` and `mate_reference` columns are **always VARCHAR** in sharded mode. Sharded alignment always loads prebuilt `.mmi` indexes, and subject names inside an `.mmi` are opaque bytes — the same contract as `align_minimap2(index_path=...)`. If the index was built from a `BIGINT` subject table via `save_minimap2_index`, the values come back as their decimal string form (e.g. `'2001'`); cast in SQL if you need integer semantics downstream.
- For `BIGINT`/`UUID` `read_id`, NULL input rows that have no entry in `read_to_shard` are silently skipped (the pipeline never sees them). The SAM `'='` mate-reference sentinel is irrelevant here because the subject side is VARCHAR — `'='` passes through verbatim, matching the existing VARCHAR contract.

**Behavior:**
- At bind time, reads the `read_to_shard` table to discover shards and validate that each `<shard_name>.mmi` file exists in `shard_directory`
- Shards are processed in parallel (one DuckDB thread per shard), each loading its `.mmi` index independently
- For each shard, only the reads assigned to that shard (via the `read_to_shard` mapping) are queried
- A read can appear in multiple shards (mapped to multiple shard_name values) and will be aligned against each
- Unmapped reads (flag 0x4) are automatically filtered out of results
- Supports both single-end and paired-end query sequences
- Supports views for both `query_table` and `read_to_shard`
- **`query_table` is read exactly once**, at the start of the scan, into a per-call TEMP table keyed by shard. Earlier versions re-read it once per shard, which silently dropped reads whenever the relation was not stable across repeated reads (a view using `nextval()` / `random()` / `now()`, a view over a concurrently-written table, or a registered single-pass Arrow stream). Reading once removes that class of failure, and is also faster on multi-shard runs — measured 39% on an 8-shard, 200k-read workload — because the relation is scanned once instead of N times.
- Because that snapshot is a TEMP table, a very large query set needs a usable `temp_directory`: TEMP data can only be offloaded there, never into a persistent database file. Single-shard runs skip the snapshot entirely (nothing is re-read, so there is nothing to guard against).

**Examples:**
```sql
-- Setup: Pre-build shard indexes (one-time)
-- Assume reference database is split into shards, each saved as a .mmi file
SELECT * FROM save_minimap2_index('shard_a_refs', 'indexes/shard_a.mmi', preset='sr');
SELECT * FROM save_minimap2_index('shard_b_refs', 'indexes/shard_b.mmi', preset='sr');

-- Load query sequences
CREATE TABLE queries AS SELECT * FROM read_fastx('reads.fastq');

-- Create read-to-shard mapping (from prior classification, e.g., RYpe)
CREATE TABLE read_to_shard AS SELECT * FROM (VALUES
    ('read1', 'shard_a'),
    ('read2', 'shard_b'),
    ('read3', 'shard_a')
) AS t(read_id, shard_name);

-- Align reads against their assigned shards in parallel
SELECT * FROM align_minimap2_sharded('queries',
    shard_directory := 'indexes/',
    read_to_shard := 'read_to_shard',
    max_secondary := 0);

-- Filter by mapping quality
SELECT read_id, reference, position, mapq
FROM align_minimap2_sharded('queries',
    shard_directory := 'indexes/',
    read_to_shard := 'read_to_shard',
    max_secondary := 0)
WHERE mapq >= 30;

-- Works with views
CREATE VIEW filtered_queries AS
    SELECT * FROM read_fastx('reads.fastq') WHERE LENGTH(sequence1) >= 50;
CREATE VIEW shard_mapping AS
    SELECT * FROM read_to_shard WHERE shard_name != 'excluded_shard';

SELECT * FROM align_minimap2_sharded('filtered_queries',
    shard_directory := 'indexes/',
    read_to_shard := 'shard_mapping',
    max_secondary := 0);

-- Use with long-read preset
SELECT * FROM align_minimap2_sharded('nanopore_reads',
    shard_directory := 'indexes/',
    read_to_shard := 'read_to_shard',
    preset := 'map-ont',
    max_secondary := 0);
```

**Error handling:**
- Error if `shard_directory` does not exist
- Error if any `<shard_name>.mmi` file referenced in `read_to_shard` does not exist in `shard_directory`
- Error if a `.mmi` file is not a valid minimap2 index
- Error if `query_table` or `read_to_shard` table/view does not exist
- Error if `read_to_shard` table is missing `read_id` or `shard_name` columns
- Error if `read_to_shard` table contains NULL `shard_name` values

**Performance notes:**
- Parallelism is one DuckDB thread per shard; control with `SET threads=N`
- Shards are sorted by read count (largest first) for better load balancing
- Pre-built indexes avoid redundant index building across runs
- `k` and `w` parameters are ignored (baked into the pre-built index); a warning is printed if specified

### SortMeRNA

#### Alignment with SAM output

rRNA filtering / alignment against one or more rRNA reference databases using [SortMeRNA](https://github.com/biocore/sortmerna) (Kopylova et al. 2012, Bioinformatics 28:3211-3217). Embedded as a statically linked library (SortMeRNA 4.4.0 fork; LGPL-3.0-or-later).

Emits the standard 21-column SAM schema shared with `align_minimap2` / `align_bowtie2`. For a schema preserving SortMeRNA's native identity / coverage / e-value / edit-distance fields, use [`align_sortmerna_rrna`](#alignment-with-rrna-output) below.

**Parameters:**
- `query_table` (VARCHAR, positional): Name of a table or view with columns `read_id` (VARCHAR, BIGINT, or UUID — see *Identifier-column types* below), `sequence1` (VARCHAR), and optionally `sequence2` (VARCHAR) when `paired := true`.
- `ref_paths` (VARCHAR[], required): List of FASTA paths for the reference database(s). The index is built once per query in-memory — re-using references across queries rebuilds the index each time.
- `num_threads` (INTEGER, default = DuckDB's thread count): Number of threads SortMeRNA's internal pool uses. The DuckDB function runs on a single DuckDB thread; parallelism is inside SortMeRNA.
- `match`, `mismatch`, `gap_open`, `gap_ext`, `score_N` (INTEGER): SW scoring. Defaults `2 / -3 / 5 / 2 / 0`.
- `evalue` (DOUBLE, default 1.0): E-value threshold.
- `seed_win_len` (UINTEGER, default 18): Seed window length.
- `num_alignments` (UINTEGER, default 1): Max alignments reported per read.
- `best` (BOOLEAN, default true): Keep only the best-scoring alignment per read.
- `paired` (BOOLEAN, default false): Paired-end mode. The query table must have a `sequence2` column; bind fails otherwise.
- `forward_only`, `reverse_only`, `full_search` (BOOLEAN): Strand-search controls.

**Output schema:** Same 21 columns as [`align_minimap2`](#single-index-alignment-with-minimap2). SortMeRNA-specific notes:
- `mapq` is always `255` — SortMeRNA does not compute mapping quality; 255 is the SAM convention for "unavailable".
- `tag_as` carries the raw Smith-Waterman score; `tag_nm` carries edit distance. Both are NULL for unaligned rows.
- `tag_xs`, `tag_ys`, `tag_xn`, `tag_xm`, `tag_xo`, `tag_xg`, `tag_yt`, `tag_md`, `tag_sa` are always NULL — SortMeRNA does not produce these.
- `stop_position` = `ref_end + 1` (1-based half-open), matching the convention used by `align_minimap2` and `align_bowtie2`.
- Paired-end: `flags & 0x2` (proper pair) is set when both mates aligned, regardless of reference. This is weaker than SAM's standard "concordant orientation within insert size" meaning — SortMeRNA is an rRNA classifier with no notion of insert size or orientation. When both mates aligned but to different references, `mate_reference` reports the actual partner reference name (rather than `=`), so cross-reference pairs remain distinguishable.

**Identifier-column types (`read_id`):**
- The input `read_id` column may be `VARCHAR`, `BIGINT`, or `UUID`. Other numeric types (INTEGER, UBIGINT, HUGEINT, DOUBLE) are rejected at bind time with the message `Column 'read_id' in table '<table>' must be VARCHAR, BIGINT, or UUID`.
- The output `read_id` column type mirrors the query side. `reference` and `mate_reference` are always `VARCHAR` — SortMeRNA references come from FASTA files on disk (`ref_paths`), never a user-provided table.
- Because the subject side is always VARCHAR, the SAM `=` mate-reference sentinel is preserved verbatim for paired-end output when both mates map to the same reference (no resolution needed — VARCHAR carries `=` directly, unlike the BIGINT-subject path in `align_minimap2` / `align_bowtie2`).
- For `BIGINT`/`UUID` `read_id`, rows with NULL `read_id` are skipped at ingress (matching the existing `ProcessSingleChunk` row-skip convention used by VARCHAR query tables).

**Caveats:**
- **No minimum-score filter:** the embedded library returns every positive Smith-Waterman hit. The `sortmerna` CLI applies an internal minimum-score threshold that our streaming API path bypasses; to reproduce CLI output row-for-row, filter on `score` (`tag_as`) or on `e_value` in SQL. Be aware that e-values are not directly comparable between the library and the CLI (see below), so a library-side e-value threshold is not guaranteed to reproduce the exact CLI row set.
- **E-value divergence from the CLI:** the library computes per-query Karlin-Altschul e-values; the CLI sums / corrects across the full database. Identity, coverage, score, CIGAR, ref_start, ref_end, and edit_distance are bit-identical between the two.
- **Process-wide serialization:** SortMeRNA's `g_run_mutex` serializes calls process-wide. Concurrent `align_sortmerna` / `align_sortmerna_rrna` queries will block on each other at the library level. The DuckDB function deliberately runs on a single DuckDB thread (`MaxThreads() == 1`) to avoid queueing behind the mutex twice.

**Examples:**

```sql
CREATE TABLE reads AS SELECT read_id, sequence1 FROM read_fastx('metaT.fastq.gz');

-- Keep only reads that align to any SILVA reference with a reasonable e-value.
CREATE TABLE rrna_reads AS
  SELECT read_id, flags, reference, position, cigar, tag_as AS score
  FROM align_sortmerna('reads', ref_paths := ['silva-bac-16s.fasta', 'silva-arc-16s.fasta'])
  WHERE (flags & 0x4) = 0;  -- drop unaligned
```

#### Alignment with rRNA output

Same aligner as `align_sortmerna` but emits SortMeRNA's native output schema, which preserves identity / coverage / e-value / edit-distance (SAM cannot carry all of these as first-class columns).

**Parameters:** identical to `align_sortmerna` above.

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `read_id` | VARCHAR, BIGINT, or UUID | Query read identifier from the input table (mirrors the query column type — see *Identifier-column types* below) |
| `aligned` | INTEGER | `1` if the read aligned, `0` otherwise |
| `strand` | INTEGER | `1` for forward, `0` for reverse-complement. `0` for unaligned rows has no meaning. |
| `ref_name` | VARCHAR | Reference sequence ID (empty for unaligned rows) |
| `ref_start` | INTEGER | 1-based inclusive start position on reference (`0` for unaligned) |
| `ref_end` | INTEGER | 1-based inclusive end position on reference (`0` for unaligned) |
| `cigar` | VARCHAR | CIGAR string, SSW convention (uses `M/I/D/S`) |
| `score` | INTEGER | Raw Smith-Waterman score (not bitscore) |
| `e_value` | DOUBLE | Per-query Karlin-Altschul e-value |
| `identity` | DOUBLE | Percent identity (0..100) at full precision |
| `coverage` | DOUBLE | Query coverage (0..100) at full precision |
| `edit_distance` | INTEGER | Number of mismatches + indels |
| `segment_idx` | INTEGER | `0` for single-end or the forward mate; `1` for the reverse mate in paired-end output |

Paired-end mode produces two rows per input row with `segment_idx` 0 (forward) and 1 (reverse), even when one or both mates failed to align.

**Identifier-column types (`read_id`):**
- The input `read_id` column may be `VARCHAR`, `BIGINT`, or `UUID`. Other numeric types (INTEGER, UBIGINT, HUGEINT, DOUBLE) are rejected at bind time with the message `Column 'read_id' in table '<table>' must be VARCHAR, BIGINT, or UUID`.
- The output `read_id` column type mirrors the query side. `ref_name` is always `VARCHAR` — it carries the free-form FASTA header text from `ref_paths`, not an identifier column.
- For `BIGINT`/`UUID` `read_id`, rows with NULL `read_id` are skipped at ingress (matching the existing `ProcessSingleChunk` row-skip convention used by VARCHAR query tables).

**Example:**

```sql
-- Assign reads to GG references, filter, aggregate.
SELECT ref_name, COUNT(*) AS hits
FROM align_sortmerna_rrna('reads',
       ref_paths := ['gg_13_8.fasta'])
WHERE aligned = 1 AND e_value <= 1e-5 AND coverage >= 80.0
GROUP BY ref_name ORDER BY hits DESC;
```
