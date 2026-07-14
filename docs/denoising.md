# Denoising

Turn raw amplicon reads into exact, error-corrected sequence variants (sub-OTUs). Use this when you have dereplicated amplicon sequences and want to strip out the spurious variants introduced by sequencing error before you build a feature table.

Quality trimming typically happens upstream (see [quality control](qc.md)); the denoised sequences and their corrected abundances feed downstream into feature tables and profiles (see [profiling](profiling.md)).

## Table of Contents

- [Deblur](#deblur) - Deblur amplicon sequence denoising into sub-OTUs.

### Deblur

Deblur amplicon sequence denoising (Amir et al. 2017, mSystems 2:e00191-16). A greedy deconvolution algorithm that removes sequencing errors from amplicon data by iteratively subtracting expected error-derived reads from less-abundant sequences. Sequences whose corrected abundance rounds to zero are removed as errors; the remainder are denoised "sub-OTUs" (sOTUs).

Designed as a composable SQL building block. Dereplication is native SQL (`GROUP BY`), alignment is `align_mafft()`, and `deblur()` does the core denoising. See the full workflow example below.

**Function signature**:

`deblur(input_table, [sample_id='col'], [options])`

**Parameters:**
- `input_table` (VARCHAR): Name of a table or view containing pre-aligned, pre-dereplicated sequences. By default must have `read_id` (VARCHAR), `sequence1` (VARCHAR), and `abundance` (integer type) columns; use `id_col`/`sequence_col`/`count_col` to override. All sequences in the sequence column must have the same aligned length and the same unaligned length (number of non-gap characters).
- `sample_id` (VARCHAR, optional): Name of a column in `input_table` to partition by. Deblur is applied independently per sample, and the sample column is prepended to the output. Unlike the two uchime functions, deblur's backend is re-entrant so samples run across DuckDB worker threads in parallel (bounded by `min(num_threads, num_samples)`).
- `id_col` (VARCHAR, default `'read_id'`): Name of the read identifier column in `input_table`.
- `sequence_col` (VARCHAR, default `'sequence1'`): Name of the sequence column. Set to `'aligned_sequence'` to chain `align_mafft(...)` directly into this function.
- `count_col` (VARCHAR, default `'abundance'`): Name of the per-sequence count column.
- `mean_error` (DOUBLE, default `0.005`): Per-base Illumina error rate. **This is the primary tuning knob.** The default 0.005 reflects MiSeq/HiSeq circa 2015. For modern NovaSeq or stitched reads (~250nt), use 0.001-0.002. Lowering `mean_error` makes denoising more conservative (fewer sequences removed). Must be > 0 and < 1.
- `error_profile` (LIST(DOUBLE), optional): Override the default 12-element error probability profile. Each element represents the fraction of reads from a true sequence that land at exactly that Hamming distance. Default: `[1, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, 0.005, 0.001, 0.001, 0.001, 0.0005]`. All values must be non-negative.
- `indel_prob` (DOUBLE, default `0.01`): Multiplicative penalty applied to corrections involving indels. A value of 0 disables indel-based corrections entirely.
- `indel_max` (INTEGER, default `3`): Maximum number of indels before a sequence is protected from correction (treated as a real variant, not an error).

**Error model:** The error profile is normalized by sequence length: `mod_factor = (1 - mean_error)^unaligned_length`. This is the probability a read has zero errors. The profile is divided by `mod_factor`, so longer sequences or higher error rates produce larger corrections. The profile shape was empirically derived from Illumina data and is deliberately lower than a binomial prediction (accounts for error collision). Use `error_profile` to provide custom calibration from mock community data.

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `read_id` | VARCHAR | Sequence identifier from input |
| `sequence` | VARCHAR | Denoised sequence (gaps stripped) |
| `abundance` | BIGINT | Corrected abundance (banker's rounding) |

**Behavior:**
- Sequences are uppercased internally (matching the Python deblur reference implementation)
- Output is ordered by corrected abundance descending
- Single-threaded (inherently sequential — each sequence's correction depends on all previous corrections)
- Empty input tables return zero rows (not an error)
- Uses banker's rounding (round half to even) for the final abundance, matching Python 3's `round()`

**Examples:**

Full workflow example:

```sql
-- 1. Read and trim to fixed length
CREATE TABLE trimmed AS
  SELECT read_id, substr(sequence1, 1, 150) AS sequence1
  FROM read_fastx('sample.fq')
  WHERE length(sequence1) >= 150;

-- 2. Dereplicate (SQL GROUP BY replaces vsearch --derep_fulllength)
CREATE TABLE dereplicated AS
  SELECT MIN(read_id) AS read_id, sequence1, COUNT(*) AS abundance
  FROM trimmed
  GROUP BY sequence1
  HAVING COUNT(*) >= 2;

-- 3. Align (required when sequences may have indels relative to each other)
-- For same-length amplicons without indels, this step can be skipped.
CREATE VIEW aligned AS
  SELECT a.read_id, a.aligned_sequence, d.abundance
  FROM align_mafft('dereplicated') a
  JOIN dereplicated d ON a.read_id = d.read_id;

-- 4. Deblur — point sequence_col at the MAFFT output column directly.
CREATE TABLE denoised AS
  SELECT * FROM deblur('aligned', sequence_col := 'aligned_sequence');

-- 5. Optional: de novo chimera removal — overrides match deblur's output
-- columns ('sequence' and 'abundance') so no aliasing subquery is needed.
SELECT * FROM detect_chimera_uchime_denovo('denoised',
                                            sequence_col := 'sequence',
                                            count_col := 'abundance');
```

Minimal example (pre-aligned sequences of equal length):

```sql
CREATE TABLE seqs(read_id VARCHAR, sequence1 VARCHAR, abundance BIGINT);
INSERT INTO seqs VALUES
  ('true_seq', 'ACGTACGTACGTACGT', 1000),
  ('error_seq', 'ACGTACGTACGTACGA', 3);

SELECT * FROM deblur('seqs');
-- Only true_seq survives (error_seq is explained by sequencing error)
```

Tuning for modern platforms:

```sql
-- NovaSeq with stitched 250nt reads
SELECT * FROM deblur('aligned_seqs', mean_error := 0.002);

-- Custom error profile from mock community calibration
SELECT * FROM deblur('aligned_seqs', error_profile := [1, 0.04, 0.01, 0.005]);
```

**Error conditions:**
- Error if table does not exist or lacks the resolved id/sequence/count columns
- Error if any of `id_col`/`sequence_col`/`count_col` is the empty string
- Error if `abundance` column is not an integer type
- Error if `mean_error` is not in the open interval (0, 1)
- Error if `indel_max` is negative
- Error if `error_profile` contains negative values or is empty
- Error if sequences have different aligned lengths or different unaligned lengths
