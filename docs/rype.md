# RYpe Classification and Extraction Functions

[RYpe](https://github.com/biocore/rype) is a minimizer-based sequence classification tool that uses RY-space encoding (purine/pyrimidine) for robust sequence matching. These functions require a RYpe index directory (`.ryxdi`), which contains a Parquet-based inverted index built from reference sequences.

All RYpe functions take a `sequence_table` parameter that references a DuckDB table or view containing sequence data. The table must have a `sequence1` column and an identifier column (default `read_id`).

## Table of Contents

- [`rype_classify`](#rype_classifyindex_path-sequence_table-id_columnread_id-threshold01-negative_indexpath) - Classify sequences against an index
- [`rype_log_ratio`](#rype_log_rationumerator_path-denominator_path-sequence_table-id_columnread_id-skip_threshold05) - Log-ratio classification between two single-bucket indices
- [`rype_extract_minimizer_set`](#rype_extract_minimizer_setsequence_table-k-w-salt6148914691236517205-id_columnread_id) - Extract deduplicated minimizer hash sets
- [`rype_extract_strand_minimizers`](#rype_extract_strand_minimizerssequence_table-k-w-salt6148914691236517205-id_columnread_id) - Extract minimizers with positions

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
