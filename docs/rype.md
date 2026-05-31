# RYpe Classification and Extraction Functions

[RYpe](https://github.com/biocore/rype) is a minimizer-based sequence classification tool that uses RY-space encoding (purine/pyrimidine) for robust sequence matching. These functions require a RYpe index directory (`.ryxdi`), which contains a Parquet-based inverted index built from reference sequences.

Most RYpe functions take a `sequence_table` parameter that references a DuckDB table or view containing sequence data. The table must have a `sequence1` column and an identifier column (default `read_id`). The exception is `rype_index_create`, which *builds* a `.ryxdi` index from a chunked reference table.

## Table of Contents

- [`rype_index_create`](#rype_index_createchunk_table-output_path-mapping_table-k64-w50-salt6148914691236517205-orienttrue-max_memory0) - Build a `.ryxdi` index from a chunked sequence table
- [`rype_classify`](#rype_classifyindex_path-sequence_table-id_columnread_id-threshold01-negative_indexpath) - Classify sequences against an index
- [`rype_log_ratio`](#rype_log_rationumerator_path-denominator_path-sequence_table-id_columnread_id-skip_threshold05) - Log-ratio classification between two single-bucket indices
- [`rype_extract_minimizer_set`](#rype_extract_minimizer_setsequence_table-k-w-salt6148914691236517205-id_columnread_id) - Extract deduplicated minimizer hash sets
- [`rype_extract_strand_minimizers`](#rype_extract_strand_minimizerssequence_table-k-w-salt6148914691236517205-id_columnread_id) - Extract minimizers with positions

## `rype_index_create(chunk_table, output_path, [mapping_table], [k=64], [w=50], [salt=6148914691236517205], [orient=true], [max_memory=0])`

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

**Output schema:** a single status row.
- `output_path` (VARCHAR): The path of the index that was created (echoes the input)
- `k` (INTEGER): The k-mer size used
- `w` (INTEGER): The window size used
- `status` (VARCHAR): `ok` on success

**Behavior:**
- A feature's chunks must be contiguous and form an ascending, 0-based, gap-free `chunk_index` sequence; the function reads them ordered by `(feature_idx, chunk_index)`.
- The large sequence/chunk data is streamed (never fully materialized); the small bucket mapping is read into memory.
- Inputs are validated at bind time: a missing table, a missing required column, `k` not in {16, 32, 64}, or `w < 1` all raise an error before any build work begins.
- A duplicate `feature_idx` in `mapping_table` is rejected.
- The index directory is written **non-atomically**: if the build fails partway, a partial, unusable directory may remain at `output_path` and should be discarded. Build to a temporary path and move it into place if atomicity is required.
- The resulting `.ryxdi` is usable by `rype_classify`, `rype_log_ratio`, and the other RYpe functions.

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

## `rype_classify(index_path, sequence_table, [id_column='read_id'], [threshold=0.1], [negative_index=path])`

Classify sequences against a RYpe index, returning bucket assignments with confidence scores.

**Parameters:**
- `index_path` (VARCHAR): Path to a RYpe index directory (`.ryxdi`)
- `sequence_table` (VARCHAR): Name of a DuckDB table or view containing sequences. Must have columns: identifier column + `sequence1` + optional `sequence2`
- `id_column` (VARCHAR, optional, default `'read_id'`): Name of the identifier column in the sequence table
- `threshold` (DOUBLE, optional, default 0.1): Minimum score threshold (0.0-1.0). Only matches with score >= threshold are returned
- `negative_index` (VARCHAR, optional): Path to a second RYpe index used as a negative filter

**Output schema:**
- `read_id` (VARCHAR): Sequence identifier from the input table
- `bucket_id` (UINTEGER): Numeric bucket identifier
- `bucket_name` (VARCHAR): Human-readable bucket name from the index
- `score` (DOUBLE): Classification confidence score (0.0-1.0)

**Behavior:**
- A sequence can match multiple buckets (one row per match above threshold)
- Sequences with no matches above the threshold produce no output rows
- If the sequence table has a `sequence2` column, paired-end classification is used
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

## `rype_log_ratio(numerator_path, denominator_path, sequence_table, [id_column='read_id'], [skip_threshold=0.5])`

Compute the log-ratio of classification scores between two single-bucket RYpe indices. For each input sequence, this returns `log10(numerator_score / denominator_score)`, indicating which index the sequence is more similar to.

**Parameters:**
- `numerator_path` (VARCHAR): Path to a single-bucket RYpe index directory (`.ryxdi`) used as the numerator
- `denominator_path` (VARCHAR): Path to a single-bucket RYpe index directory (`.ryxdi`) used as the denominator
- `sequence_table` (VARCHAR): Name of a DuckDB table or view containing sequences. Must have columns: identifier column + `sequence1` + optional `sequence2`
- `id_column` (VARCHAR, optional, default `'read_id'`): Name of the identifier column in the sequence table
- `skip_threshold` (DOUBLE, optional, default 0.5): Numerator score threshold for fast-path. Reads with numerator score >= this value skip denominator classification and get `+inf` log_ratio immediately. Set to 0 or negative to disable the fast-path

**Output schema:**
- `read_id` (VARCHAR): Sequence identifier from the input table
- `log_ratio` (DOUBLE): `log10(numerator_score / denominator_score)`. Positive values indicate the sequence matches the numerator index more strongly; negative values indicate the denominator. Special values: `+inf` (only numerator matches), `-inf` (only denominator matches), `NaN` (neither index matches)
- `fast_path` (INTEGER): 1 if the read was classified via the fast-path (skipped denominator), 0 otherwise

**Behavior:**
- Both indices must be single-bucket indices (multi-bucket indices are rejected)
- Returns exactly one row per input sequence
- If the sequence table has a `sequence2` column, paired-end classification is used
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

## `rype_extract_minimizer_set(sequence_table, k, w, [salt=6148914691236517205], [id_column='read_id'])`

Extract deduplicated minimizer hash sets from sequences for both forward and reverse complement strands.

**Parameters:**
- `sequence_table` (VARCHAR): Name of a DuckDB table or view with identifier column + `sequence1`
- `k` (BIGINT): K-mer size. Must be 16, 32, or 64 (constrained by RY-space 1-bit encoding fitting in uint64)
- `w` (BIGINT): Window size for minimizer selection. Must be > 0
- `salt` (UBIGINT, optional, default 6148914691236517205): Hash salt for reproducible but varied minimizer selection
- `id_column` (VARCHAR, optional, default `'read_id'`): Name of the identifier column

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

## `rype_extract_strand_minimizers(sequence_table, k, w, [salt=6148914691236517205], [id_column='read_id'])`

Extract minimizer hashes with their positions for both forward and reverse complement strands. Unlike `rype_extract_minimizer_set`, this preserves positional information and may contain duplicate hashes at different positions.

**Parameters:**
- `sequence_table` (VARCHAR): Name of a DuckDB table or view with identifier column + `sequence1`
- `k` (BIGINT): K-mer size. Must be 16, 32, or 64
- `w` (BIGINT): Window size for minimizer selection. Must be > 0
- `salt` (UBIGINT, optional, default 6148914691236517205): Hash salt for reproducible but varied minimizer selection
- `id_column` (VARCHAR, optional, default `'read_id'`): Name of the identifier column

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
