# Multiple sequence alignment

MIINT provides capability for multiple sequence alignment.

## Table of contents

- [MAFFT](#multiple-sequence-alignment-with-mafft) - Multiple sequence alignment with MAFFT (PartTree).
- [abPOA](#multiple-sequence-alignment-with-abpoa) - Multiple sequence alignment with abPOA. 
- [abPOA Consensus Sequence](#consensus-sequence-generation-with-abpoa) - Generate consensus sequences using abPOA. 

### Multiple sequence alignment with MAFFT

Multiple sequence alignment using MAFFT's PartTree algorithm. Reads all sequences from a sequence table/view, aligns them, and returns the aligned sequences with gap characters inserted.

MAFFT is embedded as a statically linked C library (no external binary required). The PartTree algorithm uses O(N log N) guide tree construction with k-tuple distances, making it suitable for large datasets (tested up to 5,000+ sequences).

**Parameters:**
- `table_name` (VARCHAR): Name of a table or view containing sequences to align. Must have `read_id` (VARCHAR) and `sequence1` (VARCHAR) columns. Minimum 2 sequences, each at least 6 characters. DNA and protein sequences are auto-detected. Paired-end tables (those with a `sequence2` column) are rejected at bind.
- `sample_id` (VARCHAR, optional): Name of a column in `table_name` to partition by. When provided, `align_mafft` runs one MSA per distinct sample value and prepends the sample column to the output. The per-sample ≥2-sequences and ≥6-char validations apply within each sample; the whole query aborts if any sample violates them. `sequence_index` is per-sample (0..n-1 within each sample). Join back to the input on `(<sample_id>, read_id)`.

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `sequence_index` | BIGINT | 0-based position matching input order |
| `read_id` | VARCHAR | Sequence identifier from the input `read_id` column |
| `aligned_sequence` | VARCHAR | Aligned sequence with `-` gap characters |
| `original_length` | INTEGER | Length of the input sequence (no gaps) |
| `aligned_length` | INTEGER | Length of aligned sequences (same for all rows) |

**Algorithm:** Equivalent to `mafft --quiet --preservecase --parttree`. Original case is preserved (MAFFT internally lowercases DNA; the wrapper restores the original characters after alignment).

**Thread safety:** Alignment uses a process-wide mutex (one alignment at a time). The function is safe to call from concurrent queries but they will serialize.

**Examples:**

```sql
-- Basic multiple sequence alignment
CREATE TABLE seqs AS SELECT read_id, sequence1 FROM read_fastx('sequences.fasta');
SELECT * FROM align_mafft('seqs');

-- Align, then analyze gap content
CREATE TABLE ref_16s AS SELECT read_id, sequence1 FROM read_fastx('16s_sequences.fasta');
SELECT read_id,
       original_length,
       aligned_length,
       aligned_length - original_length AS gaps_inserted
FROM align_mafft('ref_16s')
ORDER BY gaps_inserted DESC;

-- Filter before alignment — operate entirely in SQL, no temp files
CREATE TABLE filtered AS
  SELECT read_id, sequence1 FROM read_fastx('large_dataset.fasta')
  WHERE length(sequence1) >= 100;
SELECT * FROM align_mafft('filtered');

-- Align per-sample: one MSA per distinct sample value
CREATE VIEW cohort AS
  SELECT 'S1' AS sample, * FROM read_fastx('sample1.fasta')
  UNION ALL
  SELECT 'S2' AS sample, * FROM read_fastx('sample2.fasta');
SELECT * FROM align_mafft('cohort', sample_id := 'sample') ORDER BY sample, sequence_index;

-- Compute pairwise identity from aligned sequences
CREATE TABLE seqs2 AS SELECT read_id, sequence1 FROM read_fastx('seqs.fasta');
WITH aligned AS (
  SELECT read_id, aligned_sequence FROM align_mafft('seqs2')
)
SELECT a.read_id AS seq1, b.read_id AS seq2,
       sum(CASE WHEN a.c = b.c AND a.c != '-' THEN 1 ELSE 0 END)::FLOAT
       / nullif(sum(CASE WHEN a.c != '-' OR b.c != '-' THEN 1 ELSE 0 END), 0) AS identity
FROM aligned a, aligned b,
     unnest(string_split(a.aligned_sequence, '')) WITH ORDINALITY AS a(c, pos),
     unnest(string_split(b.aligned_sequence, '')) WITH ORDINALITY AS b(c, bpos)
WHERE a.read_id < b.read_id AND a.pos = b.bpos
GROUP BY a.read_id, b.read_id;
```

**Performance:** ~2x faster than native `mafft --parttree` for 1,000-5,000 sequences (no shell script overhead or temp file I/O). For 36 sequences, ~15x faster.

**Limitations:**
- All sequences are materialized in memory (required by the MSA algorithm)
- Sequences must be at least 6 characters (MAFFT internal requirement)
- Single-threaded alignment (process-wide mutex)

### Multiple sequence alignment with abPOA 

Multiple sequence alignment using [abPOA](https://github.com/yangao07/abPOA) (adaptive banded Partial Order Alignment). Reads all sequences from a sequence table/view, aligns them via a partial order graph, and returns the aligned sequences with gap characters inserted.

abPOA is embedded as a statically linked C library (no external binary required). It supports adaptive banding for 3-10x speedup over full-DP POA, minimizer-based progressive guide tree construction, and convex gap penalties.

**Parameters:**
- `table_name` (VARCHAR): Name of a table or view containing sequences to align. Must have `read_id` (VARCHAR) and `sequence1` (VARCHAR) columns. Minimum 2 sequences. DNA only. Paired-end tables (those with a `sequence2` column) are rejected at bind.
- `sample_id` (VARCHAR, optional): Name of a column in `table_name` to partition by. When provided, `align_abpoa` runs one MSA per distinct sample value and prepends the sample column to the output. `sequence_index` is per-sample (0..n-1 within each sample). Unlike `align_mafft`, per-sample alignments run in parallel (abPOA is thread-safe).
- `match` (INTEGER, default 2): Match score.
- `mismatch` (INTEGER, default 4): Mismatch penalty.
- `gap_open1` (INTEGER, default 4): First-layer gap opening penalty (affine/convex).
- `gap_open2` (INTEGER, default 24): Second-layer gap opening penalty (convex gap mode).
- `gap_ext1` (INTEGER, default 2): First-layer gap extension penalty.
- `gap_ext2` (INTEGER, default 1): Second-layer gap extension penalty.
- `align_mode` (VARCHAR, default `'global'`): Alignment mode — `'global'`, `'local'`, or `'extension'`.
- `progressive` (BOOLEAN, default `true`): Use minimizer-based guide tree for progressive alignment order. Improves quality by aligning most-similar sequences first.
- `disable_seeding` (BOOLEAN, default `false`): Disable minimizer seeding (forces full-length alignment).
- `amb_strand` (BOOLEAN, default `false`): Try reverse complement if alignment score is too low.
- `k` (INTEGER, default 19): Minimizer k-mer size.
- `w` (INTEGER, default 10): Minimizer window size.
- `min_w` (INTEGER, default 500): Minimum window size for POA between anchors.
- `bandwidth` (INTEGER, default 10): Adaptive banding base width.
- `bandwidth_frac` (FLOAT, default 0.01): Adaptive banding fraction of query length.

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `sequence_index` | BIGINT | 0-based position matching input order |
| `read_id` | VARCHAR | Sequence identifier from the input `read_id` column |
| `aligned_sequence` | VARCHAR | Aligned sequence with `-` gap characters |
| `original_length` | INTEGER | Length of the input sequence (no gaps) |
| `aligned_length` | INTEGER | Length of aligned sequences (same for all rows) |

**Thread safety:** Each alignment creates its own abPOA instance. No process-wide mutex — per-sample alignments run in true parallel.

**Examples:**

```sql
-- Basic multiple sequence alignment
CREATE TABLE seqs AS SELECT read_id, sequence1 FROM read_fastx('amplicons.fasta');
SELECT * FROM align_abpoa('seqs');

-- With progressive guide tree (default) and custom scoring
SELECT * FROM align_abpoa('seqs', match := 1, mismatch := 3, progressive := true);

-- Per-sample alignment (parallel across samples)
CREATE VIEW cohort AS
  SELECT 'S1' AS sample, * FROM read_fastx('sample1.fasta')
  UNION ALL
  SELECT 'S2' AS sample, * FROM read_fastx('sample2.fasta');
SELECT * FROM align_abpoa('cohort', sample_id := 'sample') ORDER BY sample, sequence_index;
```

**Limitations:**
- All sequences are materialized in memory (required by the POA algorithm)
- DNA only (nucleotide mode)

### Consensus sequence generation with abPOA 

Consensus sequence generation using [abPOA](https://github.com/yangao07/abPOA). Reads all sequences from a sequence table/view, builds a partial order alignment graph, and extracts one or more consensus sequences using heaviest bundling or majority voting.

Particularly useful for amplicon denoising, strain-level variant calling, and long-read error correction where multiple consensus sequences from heterogeneous samples are needed.

**Parameters:**

Same alignment parameters as `align_abpoa`, plus:
- `max_num_cons` (INTEGER, default 1): Maximum number of consensus sequences to generate. When > 1, abPOA clusters reads and produces one consensus per cluster (up to `max_num_cons`).
- `min_freq` (FLOAT, default 0.25): Minimum cluster frequency for multi-consensus mode. Clusters below this threshold are merged.
- `algorithm` (VARCHAR, default `'heaviest_bundling'`): Consensus algorithm — `'heaviest_bundling'` (highest weight path through the graph) or `'majority_voting'` (most frequent base per position).

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `consensus_id` | INTEGER | 0-based identifier for each consensus sequence |
| `consensus_sequence` | VARCHAR | The consensus nucleotide sequence |
| `consensus_length` | INTEGER | Length of the consensus sequence |
| `num_reads` | INTEGER | Number of input reads assigned to this consensus cluster |

**Examples:**

```sql
-- Single consensus from amplicon reads
CREATE TABLE amplicons AS SELECT read_id, sequence1 FROM read_fastx('amplicons.fasta');
SELECT * FROM consensus_abpoa('amplicons');

-- Multi-consensus for strain detection
SELECT * FROM consensus_abpoa('amplicons', max_num_cons := 3, min_freq := 0.1);

-- Per-sample consensus (parallel)
SELECT * FROM consensus_abpoa('cohort', sample_id := 'sample', max_num_cons := 2)
  ORDER BY sample, consensus_id;

-- Majority voting instead of heaviest bundling
SELECT * FROM consensus_abpoa('amplicons', algorithm := 'majority_voting');
```
