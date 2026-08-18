# Taxonomic classification

Assign sequencing reads to references or taxa by their minimizer/k-mer content against a prebuilt index, without performing a full alignment. Use this when you want fast, robust bucket (reference/taxon) assignments for reads — for example, screening reads as host vs. microbe, or binning reads by which reference set they best match.

These functions are powered by [RYpe](https://github.com/biocore/rype), a minimizer-based sequence classification tool that uses RY-space encoding (purine/pyrimidine) for robust sequence matching. They require a RYpe index directory (`.ryxdi`), which contains a Parquet-based inverted index built from reference sequences.

Most of these functions take a `sequence_table` parameter that references a DuckDB table or view containing sequence data. The table must have a `sequence1` column and an identifier column (default `read_id`). The exception is `rype_index_create`, which *builds* a `.ryxdi` index from a chunked reference table. Sequence tables are typically produced by the readers in [reading](reading.md). For turning classification hits into per-sample abundance tables, see [profiling](profiling.md). When your buckets are NCBI taxids, resolve them to ranks, names, and full lineages with the [NCBI taxonomy](insdc_ncbi.md#taxonomy) functions.

## Table of Contents

- [`rype_index_create`](#rype_index_create) - Build a `.ryxdi` index from a chunked sequence table
- [`rype_classify`](#rype_classify) - Classify sequences against an index
- [`rype_log_ratio`](#rype_log_ratio) - Log-ratio classification between two single-bucket indices
- [`rype_extract_minimizer_set`](#rype_extract_minimizer_set) - Extract deduplicated minimizer hash sets
- [`rype_extract_strand_minimizers`](#rype_extract_strand_minimizers) - Extract minimizers with positions

### rype_index_create

`rype_index_create(chunk_table, output_path, [mapping_table], [k=64], [w=50], [salt=6148914691236517205], [orient=true], [max_memory=0])`

Build a RYpe `.ryxdi` index directly from a DuckDB table of chunked reference sequences, without going through the `rype` CLI. The references are supplied in an at-rest *chunked* layout (one row per fixed-size block of a sequence) — the same shape microbiome reference data is stored in for efficient columnar storage. Chunks belonging to a feature are reassembled, in order, before minimizers are extracted, so minimizers spanning chunk boundaries are computed correctly.

**Parameters:**
- `chunk_table` (VARCHAR): Name of a DuckDB table or view of chunked sequences. Must have columns:
  - `feature_idx` (BIGINT): identifier of the reference sequence a chunk belongs to
  - `chunk_index` (INTEGER): 0-based block order within the feature
  - `chunk_data` (VARCHAR or BLOB): the sequence bytes for the block
- `output_path` (VARCHAR): Path of the `.ryxdi` index directory to create
- `mapping_table` (VARCHAR, optional): Name of a table or view mapping each feature to a bucket. Must have columns `feature_idx` (BIGINT) and `bucket_name` (VARCHAR), with each `feature_idx` appearing at most once. If omitted, every feature is placed into a single bucket named `unnamed-bucket`
- `k` (INTEGER, optional, default 64): K-mer size. Must be 16, 32, or 64 (constrained by RY-space 1-bit encoding fitting in uint64)
- `w` (INTEGER, optional, default 50): Minimizer window size. Must be >= 1
- `salt` (UBIGINT, optional, default 6148914691236517205): Hash salt. An index can only be classified against with matching `k`, `w`, and `salt`
- `orient` (BOOLEAN, optional, default true): Orient sequences within each bucket for better overlap before extraction
- `max_memory` (BIGINT, optional, default 0): Approximate memory budget in bytes for the build; 0 auto-detects available memory
- `feed_window_features` (BIGINT, optional, default 0): Advanced. Number of features per internal feed window (see *Behavior*); 0 auto-sizes each window from a memory budget. Set a small value to force more, smaller windows. Rarely needed.

**Output schema:** a single status row.
- `output_path` (VARCHAR): The path of the index that was created (echoes the input)
- `k` (INTEGER): The k-mer size used
- `w` (INTEGER): The window size used
- `status` (VARCHAR): `ok` on success

**Input order:** the `chunk_table` may be in **any** physical order — chunks are sorted as part of the build. Within each feature, `chunk_index` must still be 0-based and gap-free; a missing or duplicate `chunk_index` (a genuinely malformed feature) fails the build with a clear RYpe error (e.g. `expected chunk_index N but got M`, or `first chunk_index must be 0`).

**Behavior:**
- Memory is bounded regardless of corpus size. The build never sorts the whole corpus at once — a single `ORDER BY` over many 64 KB chunks would buffer everything and run out of memory, because DuckDB cannot spill large variable-length sort payloads. Instead the features are partitioned into windows (contiguous `feature_idx` ranges, auto-sized to a per-window memory budget, or set explicitly via `feed_window_features`), each window is read and sorted independently (multi-threaded, spillable), and the windows are spliced into one stream so chunks still reach RYpe fully ordered. The small bucket mapping is read into memory.
- Inputs are validated at bind time: a missing table, a missing required column, `k` not in {16, 32, 64}, or `w < 1` all raise an error before any build work begins.
- A duplicate `feature_idx` in `mapping_table` is rejected.
- The index directory is written **non-atomically**: if the build fails partway, a partial, unusable directory may remain at `output_path` and should be discarded. Build to a temporary path and move it into place if atomicity is required.
- The resulting `.ryxdi` is usable by `rype_classify`, `rype_log_ratio`, and the other functions here.

**Examples:**
```sql
-- Chunked reference table (e.g. genomes stored as 64 KB blocks)
--   chunks(feature_idx BIGINT, chunk_index INTEGER, chunk_data VARCHAR)
-- and a feature -> bucket mapping
--   refmap(feature_idx BIGINT, bucket_name VARCHAR)
SELECT * FROM rype_index_create('chunks', 'bacteria.ryxdi', mapping_table := 'refmap');

-- Build with a larger k-mer and custom window
SELECT * FROM rype_index_create('chunks', 'bacteria.ryxdi', mapping_table := 'refmap', k := 64, w := 50);

-- No mapping: all features go into a single 'unnamed-bucket'
SELECT * FROM rype_index_create('chunks', 'wholeset.ryxdi');

-- Build, then immediately classify reads against the new index
SELECT * FROM rype_index_create('chunks', 'bacteria.ryxdi', mapping_table := 'refmap');
CREATE TABLE reads AS SELECT * FROM read_fastx('reads.fastq');
SELECT * FROM rype_classify('bacteria.ryxdi', 'reads');
```

**Error conditions:**
- `chunk_table` (or `mapping_table`) does not exist
- `chunk_table` missing a required column (`feature_idx`, `chunk_index`, `chunk_data`)
- `k` not in {16, 32, 64}, or `w < 1`
- Duplicate `feature_idx` in `mapping_table`
- Malformed feature ordering (missing or duplicate `chunk_index`, or first `chunk_index` not 0)

### rype_classify

`rype_classify(index_path, sequence_table, [id_column='read_id'], [threshold=0.1], [negative_index=path], [max_memory=0], [debug=false])`

Classify sequences against a RYpe index, returning bucket assignments with confidence scores.

**Parameters:**
- `index_path` (VARCHAR): Path to a RYpe index directory (`.ryxdi`)
- `sequence_table` (VARCHAR): Name of a DuckDB table or view containing sequences. Must have columns: identifier column + `sequence1` + optional `sequence2`
- `id_column` (VARCHAR, optional, default `'read_id'`): Name of the identifier column in the sequence table
- `threshold` (DOUBLE, optional, default 0.1): Minimum score threshold (0.0-1.0). Only matches with score >= threshold are returned
- `negative_index` (VARCHAR, optional): Path to a second RYpe index used as a negative filter
- `max_memory` (BIGINT, optional, default 0): Byte budget for RYpe's internal batch sizing (how many reads it accumulates before each pass over the index). 0 derives one — see [Memory budget](#memory-budget) below
- `debug` (BOOLEAN, optional, default false): Report the chosen batch size and the memory estimates behind it through `miint_warnings()` — see [Memory budget](#memory-budget)

**Output schema:**
- `read_id` (VARCHAR): Sequence identifier from the input table
- `bucket_id` (UINTEGER): Numeric bucket identifier
- `bucket_name` (VARCHAR): Human-readable bucket name from the index
- `score` (DOUBLE): Classification confidence score (0.0-1.0)

**Behavior:**
- A sequence can match multiple buckets (one row per match above threshold)
- Sequences with no matches above the threshold produce no output rows
- Paired-end handling follows the **contents** of `sequence2`, not merely the presence of the column. This matters because `read_fastx` always emits a `sequence2` column, so single-end reads loaded the obvious way carry an all-NULL one; treating that as paired-end would halve the batch size and double the number of index passes. The check samples the first chunk of the single pass over your relation rather than scanning it separately — the relation is read exactly once, which it must be, since a view or a registered Arrow relation need not return the same rows twice. A relation that only becomes paired well after its first few thousand rows is therefore sized as single-end; set `max_memory` if that under-budgets
- Works with both tables and views

**Examples:**
```sql
-- Load sequences and classify
CREATE TABLE seqs AS SELECT * FROM read_fastx('reads.fastq');
SELECT * FROM rype_classify('my_index.ryxdi', 'seqs');

-- Classify with stricter threshold
SELECT * FROM rype_classify('my_index.ryxdi', 'seqs', threshold := 0.5);

-- Use a custom identifier column
CREATE TABLE custom_seqs AS SELECT read_id AS sample_id, sequence1 FROM read_fastx('reads.fastq');
SELECT * FROM rype_classify('my_index.ryxdi', 'custom_seqs', id_column := 'sample_id');

-- Count hits per bucket
SELECT bucket_name, COUNT(*) as hits
FROM rype_classify('my_index.ryxdi', 'seqs')
GROUP BY bucket_name
ORDER BY hits DESC;

-- Use negative index to filter host sequences
SELECT * FROM rype_classify('microbe.ryxdi', 'seqs', negative_index := 'host.ryxdi');

-- Classify paired-end reads (table must have sequence2 column)
CREATE TABLE paired AS SELECT * FROM read_fastx('R1.fastq', sequence2='R2.fastq');
SELECT * FROM rype_classify('my_index.ryxdi', 'paired');
```

**Error conditions:**
- `index_path` does not point to a valid `.ryxdi` index directory
- `sequence_table` does not exist, or is missing the identifier column or `sequence1`
- `threshold` outside the 0.0-1.0 range
- `negative_index` does not point to a valid index

### rype_log_ratio

`rype_log_ratio(numerator_path, denominator_path, sequence_table, [id_column='read_id'], [skip_threshold=0.5], [max_memory=0], [debug=false])`

Compute the log-ratio of classification scores between two single-bucket RYpe indices. For each input sequence, this returns `log10(numerator_score / denominator_score)`, indicating which index the sequence is more similar to.

**Parameters:**
- `numerator_path` (VARCHAR): Path to a single-bucket RYpe index directory (`.ryxdi`) used as the numerator
- `denominator_path` (VARCHAR): Path to a single-bucket RYpe index directory (`.ryxdi`) used as the denominator
- `sequence_table` (VARCHAR): Name of a DuckDB table or view containing sequences. Must have columns: identifier column + `sequence1` + optional `sequence2`
- `id_column` (VARCHAR, optional, default `'read_id'`): Name of the identifier column in the sequence table
- `skip_threshold` (DOUBLE, optional, default 0.5): Numerator score threshold for fast-path. Reads with numerator score >= this value skip denominator classification and get `+inf` log_ratio immediately. Set to 0 or negative to disable the fast-path
- `max_memory` (BIGINT, optional, default 0): Byte budget for RYpe's internal batch sizing (how many reads it accumulates before each pass over the index). 0 derives one — see [Memory budget](#memory-budget) below
- `debug` (BOOLEAN, optional, default false): Report the chosen batch size and the memory estimates behind it through `miint_warnings()` — see [Memory budget](#memory-budget)

**Output schema:**
- `read_id` (VARCHAR): Sequence identifier from the input table
- `log_ratio` (DOUBLE): `log10(numerator_score / denominator_score)`. Positive values indicate the sequence matches the numerator index more strongly; negative values indicate the denominator. Special values: `+inf` (only numerator matches), `-inf` (only denominator matches), `NaN` (neither index matches)
- `fast_path` (INTEGER): 1 if the read was classified via the fast-path (skipped denominator), 0 otherwise

**Behavior:**
- Both indices must be single-bucket indices (multi-bucket indices are rejected)
- Returns exactly one row per input sequence
- Paired-end handling follows the **contents** of `sequence2`, not merely the presence of the column. This matters because `read_fastx` always emits a `sequence2` column, so single-end reads loaded the obvious way carry an all-NULL one; treating that as paired-end would halve the batch size and double the number of index passes. The check samples the first chunk of the single pass over your relation rather than scanning it separately — the relation is read exactly once, which it must be, since a view or a registered Arrow relation need not return the same rows twice. A relation that only becomes paired well after its first few thousand rows is therefore sized as single-end; set `max_memory` if that under-budgets
- Works with both tables and views
- Swapping numerator and denominator negates the log_ratio (symmetry property)

**Examples:**
```sql
-- Load sequences and compute log-ratio between two indices
CREATE TABLE seqs AS SELECT * FROM read_fastx('reads.fastq');
SELECT * FROM rype_log_ratio('host.ryxdi', 'microbe.ryxdi', 'seqs');

-- Disable fast-path for exact classification against both indices
SELECT * FROM rype_log_ratio('host.ryxdi', 'microbe.ryxdi', 'seqs', skip_threshold := 0.0);

-- Filter for sequences that strongly match the numerator
SELECT read_id, log_ratio
FROM rype_log_ratio('host.ryxdi', 'microbe.ryxdi', 'seqs', skip_threshold := 0.0)
WHERE log_ratio > 1.0
ORDER BY log_ratio DESC;

-- Classify reads as host vs microbe based on log-ratio sign
SELECT read_id,
       CASE WHEN isinf(log_ratio) AND log_ratio > 0 THEN 'host_only'
            WHEN isinf(log_ratio) AND log_ratio < 0 THEN 'microbe_only'
            WHEN isnan(log_ratio) THEN 'unclassified'
            WHEN log_ratio > 0 THEN 'host'
            WHEN log_ratio < 0 THEN 'microbe'
            ELSE 'ambiguous' END as classification
FROM rype_log_ratio('host.ryxdi', 'microbe.ryxdi', 'seqs', skip_threshold := 0.0);

-- Use with paired-end reads
CREATE TABLE paired AS SELECT * FROM read_fastx('R1.fastq', sequence2='R2.fastq');
SELECT * FROM rype_log_ratio('host.ryxdi', 'microbe.ryxdi', 'paired');

-- Use a custom identifier column
SELECT * FROM rype_log_ratio('host.ryxdi', 'microbe.ryxdi', 'seqs', id_column := 'sample_id');
```

**Error conditions:**
- Either `numerator_path` or `denominator_path` is not a valid `.ryxdi` index directory
- Either index is not a single-bucket index (multi-bucket indices are rejected)
- `sequence_table` does not exist, or is missing the identifier column or `sequence1`

### rype_extract_minimizer_set

`rype_extract_minimizer_set(sequence_table, k, w, [salt=6148914691236517205], [id_column='read_id'], [max_memory=0], [debug=false])`

Extract deduplicated minimizer hash sets from sequences for both forward and reverse complement strands.

**Parameters:**
- `sequence_table` (VARCHAR): Name of a DuckDB table or view with identifier column + `sequence1`
- `k` (BIGINT): K-mer size. Must be 16, 32, or 64 (constrained by RY-space 1-bit encoding fitting in uint64)
- `w` (BIGINT): Window size for minimizer selection. Must be > 0
- `salt` (UBIGINT, optional, default 6148914691236517205): Hash salt for reproducible but varied minimizer selection
- `id_column` (VARCHAR, optional, default `'read_id'`): Name of the identifier column
- `max_memory` (BIGINT, optional, default 0): Byte budget for RYpe's internal batch sizing (how many reads it accumulates before each pass over the index). 0 derives one — see [Memory budget](#memory-budget) below
- `debug` (BOOLEAN, optional, default false): Report the chosen batch size and the memory estimates behind it through `miint_warnings()` — see [Memory budget](#memory-budget)

**Output schema:**
- `read_id` (VARCHAR): Sequence identifier from the input table
- `fwd_set` (UBIGINT[]): Sorted, deduplicated minimizer hashes from the forward strand
- `rc_set` (UBIGINT[]): Sorted, deduplicated minimizer hashes from the reverse complement strand

**Behavior:**
- Returns one row per input sequence
- Sequences shorter than k produce empty lists
- Hash sets are sorted and deduplicated (set semantics)
- Different salt values produce different hash sets for the same input

**Examples:**
```sql
-- Extract minimizer sets with k=32, w=10
CREATE TABLE seqs AS SELECT * FROM read_fastx('reads.fastq');
SELECT * FROM rype_extract_minimizer_set('seqs', 32, 10);

-- Use k=16 for shorter k-mers
SELECT read_id, len(fwd_set) as num_minimizers
FROM rype_extract_minimizer_set('seqs', 16, 5);

-- Extract with a salt value
SELECT * FROM rype_extract_minimizer_set('seqs', 32, 10, salt := 42);

-- Compare minimizer overlap between two sequences
WITH mins AS (
    SELECT read_id, fwd_set
    FROM rype_extract_minimizer_set('seqs', 32, 10)
)
SELECT a.read_id, b.read_id,
       len(list_intersect(a.fwd_set, b.fwd_set)) as shared_minimizers
FROM mins a, mins b
WHERE a.read_id < b.read_id;
```

**Error conditions:**
- `sequence_table` does not exist, or is missing the identifier column or `sequence1`
- `k` not in {16, 32, 64}
- `w` not greater than 0

### rype_extract_strand_minimizers

`rype_extract_strand_minimizers(sequence_table, k, w, [salt=6148914691236517205], [id_column='read_id'], [max_memory=0], [debug=false])`

Extract minimizer hashes with their positions for both forward and reverse complement strands. Unlike `rype_extract_minimizer_set`, this preserves positional information and may contain duplicate hashes at different positions.

**Parameters:**
- `sequence_table` (VARCHAR): Name of a DuckDB table or view with identifier column + `sequence1`
- `k` (BIGINT): K-mer size. Must be 16, 32, or 64
- `w` (BIGINT): Window size for minimizer selection. Must be > 0
- `salt` (UBIGINT, optional, default 6148914691236517205): Hash salt for reproducible but varied minimizer selection
- `id_column` (VARCHAR, optional, default `'read_id'`): Name of the identifier column
- `max_memory` (BIGINT, optional, default 0): Byte budget for RYpe's internal batch sizing (how many reads it accumulates before each pass over the index). 0 derives one — see [Memory budget](#memory-budget) below
- `debug` (BOOLEAN, optional, default false): Report the chosen batch size and the memory estimates behind it through `miint_warnings()` — see [Memory budget](#memory-budget)

**Output schema:**
- `read_id` (VARCHAR): Sequence identifier from the input table
- `fwd_hashes` (UBIGINT[]): Minimizer hashes from the forward strand
- `fwd_positions` (UBIGINT[]): Positions corresponding to each forward hash
- `rc_hashes` (UBIGINT[]): Minimizer hashes from the reverse complement strand
- `rc_positions` (UBIGINT[]): Positions corresponding to each reverse complement hash

**Behavior:**
- Returns one row per input sequence
- Sequences shorter than k produce empty lists
- Hash and position arrays always have the same length per strand (i.e., `len(fwd_hashes) = len(fwd_positions)`)
- Positions are 0-based offsets into the sequence

**Examples:**
```sql
-- Extract strand minimizers with positions
CREATE TABLE seqs AS SELECT * FROM read_fastx('reads.fastq');
SELECT * FROM rype_extract_strand_minimizers('seqs', 32, 10);

-- Verify hash/position alignment
SELECT read_id,
       len(fwd_hashes) = len(fwd_positions) as fwd_aligned,
       len(rc_hashes) = len(rc_positions) as rc_aligned
FROM rype_extract_strand_minimizers('seqs', 32, 10);

-- Examine minimizer positions along a sequence
SELECT read_id, unnest(fwd_hashes) as hash, unnest(fwd_positions) as pos
FROM rype_extract_strand_minimizers('seqs', 32, 10)
ORDER BY read_id, pos;
```

**Error conditions:**
- `sequence_table` does not exist, or is missing the identifier column or `sequence1`
- `k` not in {16, 32, 64}
- `w` not greater than 0

---

## Memory budget

`rype_classify`, `rype_log_ratio`, `rype_extract_minimizer_set` and
`rype_extract_strand_minimizers` all accept `max_memory` (BIGINT, bytes,
default 0), which bounds how many reads RYpe accumulates before each pass over
the index. That accumulation is the dominant memory term on a large corpus, so
this is the knob to reach for if a classify is being killed.

**What 0 does.** RYpe's own auto-detection resolves to the cgroups or SLURM
limit for the whole process, which inside DuckDB double-counts whatever DuckDB
is already holding — both would size themselves against the same allocation
without either knowing about the other. Instead, miint subtracts what DuckDB is
actually holding when the function starts, plus a tenth of the allocation as
room for it to grow into.

It deliberately does not subtract DuckDB's `memory_limit`. That is a ceiling on
the buffer pool rather than a reservation, and it defaults to about 80% of RAM,
so subtracting it would cut the budget roughly fourfold on an unconfigured host.
Since batch size scales with the budget and every batch costs a full pass over
the index, that would trade an occasional out-of-memory kill for a reliable
severalfold slowdown. If you need the two budgets to be provably disjoint rather
than estimated, say so with `max_memory` — that is what the parameter is for.

**Bounding the Arrow batch.** Separately from `max_memory`, the sequence bytes
handed to RYpe at any one moment are capped by the
`miint_rype_arrow_batch_bytes` setting (bytes, default 256 MiB). This is what
stops peak memory from scaling with the corpus: without it a batch is sized by a
row count, so its byte size is unbounded.

```sql
SET miint_rype_arrow_batch_bytes = 67108864;  -- 64 MiB
SET miint_rype_arrow_batch_bytes = 0;         -- no ceiling (one unbounded batch)
```

Lowering it does not change how often RYpe passes over the index — that is
`max_memory`'s job — so it is close to free. Raising it is rarely useful; the
default is already large enough that per-batch fixed costs are negligible.

**When to set it explicitly.** Pass a value when you need a hard guarantee
rather than a heuristic — a container sized close to the job, a shared node, or
a scheduler that kills on RSS. A stated `max_memory` is used verbatim.

```sql
-- Bound RYpe to 8 GB regardless of what is detected
SELECT * FROM rype_classify('host.ryxdi', 'reads', max_memory := 8000000000);

-- Pair it with DuckDB's own limit when the two must fit one allocation
SET memory_limit = '16GB';
SELECT * FROM rype_classify('host.ryxdi', 'reads', max_memory := 8000000000);
```

Lowering `max_memory` lowers peak memory but increases the number of index
passes, and index loading dominates runtime — so it trades wall clock for
footprint, not the other way round.

Two notes on what `max_memory` does **not** cover:

- It does not bound the Arrow batches miint builds to feed RYpe. Those are
  capped separately at a fixed 256 MiB of sequence payload, independent of
  corpus size and of this setting.
- DuckDB's `memory_limit` does not bound RYpe either, in the other direction:
  RYpe allocates outside DuckDB's buffer manager, so lowering `memory_limit`
  alone will not reduce peak RSS for these functions.

### Seeing what was chosen

Batch sizing is invisible by default, which makes a slow classify hard to
diagnose: index loading is roughly 99.9% of the work, so the number of batches
is very nearly a direct multiplier on runtime. Pass `debug := true` to have the
chosen batch size and the estimates behind it recorded in `miint_warnings()`:

```sql
SELECT count(*) FROM rype_classify('host.ryxdi', 'reads', debug := true);

SELECT message FROM miint_warnings();
-- rype_classify debug: classification batch = 1354645 reads; memory budget
-- 21.30 GB, estimated 2.87 GB per batch and 3.12 GB peak; avg read length 150,
-- is_paired 0
```

Read it as: RYpe will classify 1,354,645 reads per pass over the index, and
expects to need about 3.12 GB at peak against the 21.30 GB it was given. If the
batch size is far below your read count, the corpus is being classified in
several passes — raise `max_memory` if the memory is genuinely available, or
accept the extra passes if it is not.

`debug := true` also reports on `rype_log_ratio`,
`rype_extract_minimizer_set` and `rype_extract_strand_minimizers`. The
extraction functions report a per-read cost rather than shard-aware estimates,
because they do not load an index at all.

A `max_memory` too small for even the minimum batch is not an error: RYpe falls
back to a 1000-read floor. `debug := true` is how you notice that has happened.
