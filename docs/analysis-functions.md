# Analysis Functions

Functions for higher-level genomic analysis, sequence manipulation, and pairwise alignment.

## Table of Contents

- [`woltka_ogu`](#woltka_ogurelation-sequence_id_field-sample_id) - OGU counts (global or per-sample)
- [`sequence_dna_reverse_complement` / `sequence_rna_reverse_complement`](#sequence_dna_reverse_complementsequence-and-sequence_rna_reverse_complementsequence) - Reverse complement
- [`sequence_dna_as_regexp` / `sequence_rna_as_regexp`](#sequence_dna_as_regexpsequence-and-sequence_rna_as_regexpsequence) - IUPAC to regex
- [`compress_intervals`](#compress_intervalsstart-stop) - Merge overlapping intervals
- [`compute_coverage_depth`](#compute_coverage_depthposition-stop_position-cigar-reference_length-mode) - Per-position depth aggregate
- [`compute_msa_consensus`](#compute_msa_consensusaligned_seq-qual) - Q-aware MSA column consensus with HP post-correction
- [`genome_coverage`](#genome_coveragealignments-subject_total_length-subject_genome_id) - Compute genome coverage from alignments
- [Pairwise Alignment Functions](#pairwise-alignment-functions) - WFA2-based pairwise alignment
- [`formula`](#formulaformula_string) - Chemical formula to monoisotopic mass
- [`massql`](#massqlquery-source) - MassQL query language for mass spectrometry
- [Utility Functions](#utility-functions) - `miint_version()` and others

## `woltka_ogu(relation, sequence_id_field [, sample_id])`

Compute [Woltka](https://github.com/qiyunzhu/woltka) OGU (Operational Genomic Unit) counts over SAM-like alignment data. Implements Woltka's classification algorithm, assigning reads to taxonomic units while fractionally distributing multi-mapped reads. When the optional `sample_id` named parameter is supplied, the aggregation runs in parallel across distinct sample values — one DuckDB query per sample on a dedicated per-thread connection — which bounds memory to a single sample's footprint.

**Parameters:**
- `relation` (VARCHAR): Name of a table or view containing SAM-like alignment data (must be a resolvable catalog name; pass as a string literal)
- `sequence_id_field` (VARCHAR): Name of the column holding sequence identifiers — typically `read_id`, or a numeric index column for better hash performance
- `sample_id` (VARCHAR, named, optional): Name of the column holding sample identifiers. When supplied, adds that column to the output and runs the aggregation per distinct sample value in parallel.

**Required columns in relation:**
- Column named by `sequence_id_field`: read/sequence identifier
- `reference` (VARCHAR): reference sequence name (becomes `feature_id`)
- `flags` (USMALLINT): SAM alignment flags
- When `sample_id` is supplied: the named column (any comparable type) — NULLs are rejected at bind time

**Returns:**
- When `sample_id` is omitted: `(feature_id VARCHAR, value DOUBLE)`
- When `sample_id` is supplied: `(<sample_id_column> <its_type>, feature_id VARCHAR, value DOUBLE)` — the first column's name matches the value you passed to `sample_id`.

**Algorithm:**
1. Orients reads using alignment flags (forward/reverse via `alignment_is_read1`).
2. For each read orientation, divides 1 by the number of unique features aligned to.
3. Aggregates fractional counts per feature (and per sample, when `sample_id` is used).

**Correctness assumption for per-sample mode:** read IDs are unique across samples. The per-sample subset then yields the same distribution as the global aggregation.

**Examples:**
```sql
-- Global aggregation across the whole relation
SELECT * FROM woltka_ogu('my_alignments', 'read_id');

-- Per-sample aggregation — one aggregation per distinct sample value, in parallel
SELECT * FROM woltka_ogu(
    'my_alignments',
    'read_id',
    sample_id := 'sample_id'
);

-- Filter high-quality alignments via a view, then classify
CREATE OR REPLACE VIEW primary_alignments AS
    SELECT * FROM read_alignments('alignments.bam')
    WHERE alignment_is_primary(flags) AND mapq >= 20;

SELECT * FROM woltka_ogu('primary_alignments', 'read_id');

-- Multi-sample: union sources under a sample_id column, then classify per sample
CREATE OR REPLACE VIEW all_samples AS
    SELECT *, 'sample1' AS sample_id FROM read_alignments('sample1.bam')
    UNION ALL
    SELECT *, 'sample2' AS sample_id FROM read_alignments('sample2.bam')
    UNION ALL
    SELECT *, 'sample3' AS sample_id FROM read_alignments('sample3.bam');

SELECT * FROM woltka_ogu('all_samples', 'read_id', sample_id := 'sample_id')
ORDER BY sample_id, feature_id;

-- Export per-sample results to BIOM for downstream analysis
COPY (
    SELECT * FROM woltka_ogu('my_alignments', 'read_id', sample_id := 'sample_id')
) TO 'ogu_table.biom' (FORMAT BIOM, COMPRESSION 'gzip');

-- Export global results to BIOM (add a sample_id column)
COPY (
    SELECT feature_id, 'MySample' AS sample_id, value
    FROM woltka_ogu('my_alignments', 'read_id')
) TO 'ogu_single.biom' (FORMAT BIOM);
```

**Notes:**
- Arguments `relation` and `sequence_id_field` (and `sample_id` if used) must be quoted string literals; they are catalog/column names resolved at bind time.
- Multi-mapped reads are fractionally assigned: each mapping receives weight 1/N where N is the number of unique references the read maps to, within the same orientation.
- Read orientation (forward/reverse) is considered separately via the `alignment_is_read1()` flag — paired-end reads are handled automatically.
- For better performance on large datasets, add a numeric index column and pass it as `sequence_id_field` instead of `read_id`.
- Output row order is non-deterministic when `sample_id` is used (parallel per-sample execution). Use an explicit `ORDER BY` if stable ordering is required.

## `sequence_dna_reverse_complement(sequence)` and `sequence_rna_reverse_complement(sequence)`

Calculate the reverse complement of DNA or RNA sequences. Supports full IUPAC nucleotide ambiguity codes and preserves case.

**Parameters:**
- `sequence` (VARCHAR): DNA or RNA sequence string

**Returns:** VARCHAR - The reverse complement of the input sequence

**Behavior:**
- Reverses the sequence order (5' to 3' becomes 3' to 5')
- Complements each base according to Watson-Crick pairing rules
- Preserves uppercase/lowercase in the input
- Supports gap characters (`.` and `-`) which map to themselves
- Strict molecular type validation: DNA function rejects U bases, RNA function rejects T bases

**Supported bases:**
- **DNA**: A<->T, C<->G, plus IUPAC codes (R<->Y, S<->S, W<->W, K<->M, B<->V, D<->H, N<->N)
- **RNA**: A<->U, C<->G, plus IUPAC codes (R<->Y, S<->S, W<->W, K<->M, B<->V, D<->H, N<->N)

**Examples:**
```sql
-- Basic DNA reverse complement
SELECT sequence_dna_reverse_complement('ATCG');
-- Returns: CGAT

-- Basic RNA reverse complement
SELECT sequence_rna_reverse_complement('AUCG');
-- Returns: CGAU

-- Works with IUPAC ambiguity codes
SELECT sequence_dna_reverse_complement('ACGTMRWSYKVHDBN.-');
-- Returns: -.NVHDBMRSWYKACGT

-- Case is preserved
SELECT sequence_dna_reverse_complement('AcGt');
-- Returns: aCgT

-- Process sequences from FASTQ files
SELECT read_id,
       sequence1,
       sequence_dna_reverse_complement(sequence1) AS rev_comp
FROM read_fastx('sequences.fastq');

-- Find palindromic sequences (equal to their reverse complement)
SELECT read_id, sequence1
FROM read_fastx('sequences.fastq')
WHERE sequence1 = sequence_dna_reverse_complement(sequence1);

-- Error: DNA function rejects U bases
SELECT sequence_dna_reverse_complement('AUCG');
-- Error: Invalid DNA base 'U'

-- Error: RNA function rejects T bases
SELECT sequence_rna_reverse_complement('ATCG');
-- Error: Invalid RNA base 'T'
```

**IUPAC Ambiguity Code Reference:**
- R = A or G (purine)
- Y = C or T/U (pyrimidine)
- S = G or C (strong)
- W = A or T/U (weak)
- K = G or T/U (keto)
- M = A or C (amino)
- B = not A (C, G, T/U)
- D = not C (A, G, T/U)
- H = not G (A, C, T/U)
- V = not T/U (A, C, G)
- N = any base

## `sequence_dna_as_regexp(sequence)` and `sequence_rna_as_regexp(sequence)`

Convert DNA or RNA sequences with IUPAC ambiguity codes to regular expression patterns. Useful for pattern matching with degenerate primers and probes.

**Parameters:**
- `sequence` (VARCHAR): DNA or RNA sequence string with IUPAC codes

**Returns:** VARCHAR - Regular expression pattern

**Behavior:**
- Unambiguous bases (A, C, G, T/U) remain unchanged
- Ambiguous IUPAC codes expand to character classes (e.g., R -> `[AG]`, N -> `[ACGT]`)
- Preserves uppercase/lowercase in the output
- Gap characters (`-` and `.`) convert to `.` (regex wildcard matching any character)
- Strict molecular type validation: DNA function rejects U bases, RNA function rejects T bases

**Expansion rules:**
- **Unambiguous bases**: A, C, G, T (DNA) or U (RNA) -> no brackets
- **Two-base codes**: R -> `[AG]`, Y -> `[CT]` or `[CU]`, S -> `[CG]`, W -> `[AT]` or `[AU]`, K -> `[GT]` or `[GU]`, M -> `[AC]`
- **Three-base codes**: B -> `[CGT]` or `[CGU]`, D -> `[AGT]` or `[AGU]`, H -> `[ACT]` or `[ACU]`, V -> `[ACG]`
- **Any base**: N -> `[ACGT]` or `[ACGU]`
- **Gap characters**: `-` or `.` -> `.` (matches any character)

**Examples:**
```sql
-- Basic DNA sequence (unambiguous)
SELECT sequence_dna_as_regexp('ATCG');
-- Returns: ATCG

-- Degenerate primer with ambiguous positions
SELECT sequence_dna_as_regexp('ATNGG');
-- Returns: AT[ACGT]GG

-- Multiple IUPAC codes
SELECT sequence_dna_as_regexp('RYMKSW');
-- Returns: [AG][CT][AC][GT][CG][AT]

-- Case is preserved
SELECT sequence_dna_as_regexp('AcGtRy');
-- Returns: AcGt[AG][ct]

-- RNA sequence with ambiguity
SELECT sequence_rna_as_regexp('AUNGG');
-- Returns: AU[ACGU]GG

-- Use with DuckDB's regexp_matches for pattern searching
SELECT read_id, sequence1
FROM read_fastx('sequences.fastq')
WHERE regexp_matches(sequence1, sequence_dna_as_regexp('ATNGG'));

-- Find sequences matching a degenerate probe
CREATE TABLE probes AS
  SELECT 'probe1' AS name, 'GCRAA' AS sequence
  UNION ALL SELECT 'probe2', 'ATNGG';

SELECT p.name, f.read_id, f.sequence1
FROM read_fastx('sequences.fastq') f
CROSS JOIN probes p
WHERE regexp_matches(f.sequence1, sequence_dna_as_regexp(p.sequence));

-- Count reads matching a consensus pattern
SELECT COUNT(*) AS matching_reads
FROM read_fastx('sequences.fastq')
WHERE regexp_matches(sequence1, sequence_dna_as_regexp('ACGTNNNNACGT'));

-- Error: DNA function rejects U bases
SELECT sequence_dna_as_regexp('AUNGG');
-- Error: Invalid DNA base 'U'

-- Error: RNA function rejects T bases
SELECT sequence_rna_as_regexp('ATNGG');
-- Error: Invalid RNA base 'T'
```

**Use cases:**
- **Degenerate primer matching**: Search for sequences matching primers with ambiguous positions
- **Motif finding**: Identify consensus sequences with variable positions
- **Probe design validation**: Check which sequences match a degenerate probe
- **Quality control**: Filter reads matching specific sequence patterns
- **Pattern-based classification**: Group sequences by motif presence

**Note:** Gap characters become `.` (regex wildcard), which matches any single character. This is useful for representing unknown or variable positions in alignments.

## `compress_intervals(start, stop)`

Aggregate function that merges overlapping genomic intervals into a minimal set of non-overlapping intervals. Useful for computing coverage regions, reducing redundant intervals, and analyzing read depth.

**Parameters:**
- `start` (BIGINT): Start position of interval
- `stop` (BIGINT): Stop position of interval

**Returns:** `LIST<STRUCT(start BIGINT, stop BIGINT)>` - Array of merged intervals, sorted by start position

**Behavior:**
- Intervals that overlap or touch (stop1 = start2) are merged into a single interval
- Automatically compresses state when accumulating >1M intervals (prevents memory issues)
- Thread-safe: works correctly with parallel GROUP BY operations
- Returns NULL for empty groups

**Examples:**
```sql
-- Calculate coverage regions per reference from SAM/BAM alignments
SELECT reference,
       compress_intervals(position, stop_position) AS coverage
FROM read_alignments('alignments.bam')
GROUP BY reference;

-- Count number of distinct coverage regions per reference
SELECT reference,
       LEN(compress_intervals(position, stop_position)) AS num_regions
FROM read_alignments('alignments.bam')
GROUP BY reference;

-- Calculate total covered bases per reference
SELECT reference,
       SUM(c.stop - c.start) AS total_coverage
FROM (
  SELECT reference,
         UNNEST(compress_intervals(position, stop_position)) AS c
  FROM read_alignments('alignments.bam')
  GROUP BY reference
);

-- Find gaps in coverage (regions with no reads)
WITH coverage AS (
  SELECT reference,
         UNNEST(compress_intervals(position, stop_position)) AS interval
  FROM read_alignments('alignments.bam')
  GROUP BY reference
)
SELECT reference,
       interval.stop AS gap_start,
       LEAD(interval.start) OVER (PARTITION BY reference ORDER BY interval.start) AS gap_end
FROM coverage
WHERE LEAD(interval.start) OVER (PARTITION BY reference ORDER BY interval.start) IS NOT NULL;

-- Merge overlapping genomic features from different sources
CREATE TABLE features AS
  SELECT 'gene1' AS name, 'genome1' AS ref, 1000 AS start, 2000 AS stop
  UNION ALL SELECT 'gene2', 'genome1', 1500, 2500
  UNION ALL SELECT 'gene3', 'genome1', 3000, 4000;

SELECT ref,
       compress_intervals(start, stop) AS merged_regions
FROM features
GROUP BY ref;
```

**Performance Notes:**
- Automatic periodic compression at 1M intervals prevents memory bloat with large datasets
- Multi-threaded aggregation: each thread maintains its own state, merged at finalization
- Algorithm: sorts intervals by start position, then single-pass merge (O(n log n))

## `compute_coverage_depth(position, stop_position, cigar, reference_length, mode)`

Aggregate function that computes per-position depth of coverage across a reference sequence. Follows `samtools depth` semantics: two modes for handling deletion (D) operations, and N (ref-skip) operations are always excluded.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `position` | BIGINT | 1-based inclusive start position |
| `stop_position` | BIGINT | 1-based exclusive end position (half-open, matches `read_alignments` output) |
| `cigar` | VARCHAR | CIGAR string |
| `reference_length` | BIGINT | Determines output list size (must be positive, max 2,000,000,000) |
| `mode` | VARCHAR | `'include_deletions'` or `'exclude_deletions'` (must be a constant) |

**Returns:** `LIST(UINTEGER)` — element `i` is the depth at 1-based position `i+1`.

**Mode Semantics:**
- `'include_deletions'`: M/=/X/D count as coverage, N excluded (equivalent to `samtools depth -J`)
- `'exclude_deletions'`: only M/=/X count, D and N excluded (equivalent to `samtools depth` default)

**Behavior:**
- NULL rows are ignored
- Empty groups (all NULLs) return NULL
- All rows in a group must have the same `reference_length`
- Mode parameter must be a bind-time constant string (not a column reference)

**SQL Examples:**

```sql
-- Per-position depth for a single reference
SELECT compute_coverage_depth(position, stop_position, cigar, 1000, 'exclude_deletions')
FROM read_alignments('alignments.bam')
WHERE reference = 'chr1';

-- Per-reference depth with GROUP BY
SELECT reference,
       compute_coverage_depth(position, stop_position, cigar,
           ref_lengths.length, 'include_deletions') AS depths
FROM read_alignments('alignments.bam') AS a
JOIN ref_lengths ON a.reference = ref_lengths.name
GROUP BY reference, ref_lengths.length;

-- Mean depth via UNNEST
SELECT AVG(depth) AS mean_depth
FROM (
  SELECT UNNEST(compute_coverage_depth(position, stop_position, cigar, 1000, 'exclude_deletions')) AS depth
  FROM read_alignments('alignments.bam')
  WHERE reference = 'chr1'
);

-- Positions with depth >= 10
SELECT position, depth
FROM (
  SELECT ROW_NUMBER() OVER () AS position,
         UNNEST(compute_coverage_depth(position, stop_position, cigar, 1000, 'exclude_deletions')) AS depth
  FROM read_alignments('alignments.bam')
  WHERE reference = 'chr1'
)
WHERE depth >= 10;
```

**Performance Notes:**
- Memory: allocates `reference_length * 4` bytes per group (e.g., ~1GB for a 250M-position human chromosome)
- Multi-threaded aggregation: each thread maintains independent state, merged at finalization
- Fast path: reads with only M/=/X operations (no N or excluded D) skip CIGAR walking
- Maximum reference_length is 2,000,000,000; use `compress_intervals` for larger references

## `compute_msa_consensus(aligned_seq, qual)`

Aggregate function that computes a quality-aware column consensus from a
multiple sequence alignment (MSA). Designed as the consensus-building
step for long-read amplicon UMI pipelines (Karst-protocol consensus on
Revio HiFi), replacing Racon polishing for HiFi-grade per-base quality.
General enough to be useful anywhere you have a per-bin MSA + per-read
quality.

Algorithm (full design rationale: HP-indel-dominated residual error model
on modern HiFi makes column voting alone insufficient):

1. **Per-column 5-state log-likelihood vote** over `{A, C, G, T, '-'}`.
   Each base observation uses `p_err = 10^(-q/10)`; gap observations use
   a fixed `p_err = 0.05`. The argmax base is emitted; columns where gap
   wins are suppressed from the output.
2. **Posterior Q** derived from `log(sum_{k≠best} exp(ll[k]) / sum_all)`,
   clamped to `[0, 60]` and emitted as UTINYINT alongside the base.
3. **HP-aware post-correction**: each homopolymer run (length ≥ 2) in the
   column-voted consensus is replaced by `median(per-read HP length at the
   corresponding ungapped locus)`. Critical: MAFFT places HP gaps
   inconsistently and naive column voting biases HP length short.
4. **Bin of size 1** bypasses the MSA pipeline entirely and emits the
   lone read's gap-stripped sequence + unchanged qual (the MSA is
   degenerate at n=1).

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `aligned_seq` | VARCHAR | MSA row for one read (gaps encoded as `-`), as produced by [`align_mafft`](table-functions.md#align_maffttable_name-sample_idcol) |
| `qual` | LIST(UTINYINT) | The read's original *ungapped* per-base Phred qual (no ASCII offset; matches `read_fastx`) |

The deliberate type asymmetry means callers join the MAFFT row back to
the original quality list by `read_id`; Q never enters MAFFT, so it
never has to be reassembled across gap insertions.

**Returns:** STRUCT with fields:

| Field | Type | Description |
|-------|------|-------------|
| `seq` | VARCHAR | Consensus sequence (no gaps; HP-corrected) |
| `qual` | LIST(UTINYINT) | Posterior per-base Q list (parallel to `seq`, capped at 60) |

**Behavior:**
- NULL rows are ignored (`IgnoreNull() = true`). Groups where every row
  is NULL return NULL.
- All rows within a group must share the same MSA width; throws on a
  width mismatch.
- Each row's `qual` length must equal its sequence's ungapped count;
  throws on mismatch.
- Tie-break is deterministic: log-likelihood difference < 1e-9 falls back
  to per-candidate Q-sum, then alphabetical (`A < C < G < T < -`).
- Lower-median is used for even HP-length samples (a length some read
  actually observed).
- Eager materialization: all `(aligned_seq, qual)` rows are held in the
  aggregate state. Memory is `O(rows × MSA_width)` per group, which is
  comfortable for typical UMI bin sizes (5–30 reads × ~5 kb amplicons).

**Example:**
```sql
-- Per-bin consensus via the standard MAFFT → consensus chain.
-- See test/sql/umi_pipeline_integration.test for an end-to-end example.

-- 1. Per-bin MSA
CREATE TABLE mafft_out AS
SELECT * FROM align_mafft('amplicons_with_bin', sample_id := 'bin_id');

-- 2. Join MAFFT's gapped seq back to the original qual list
CREATE VIEW consensus_input AS
SELECT m.bin_id, m.read_id, m.aligned_sequence, a.qual
FROM mafft_out m
JOIN read_amplicons a USING (read_id);

-- 3. Aggregate per bin
SELECT bin_id, (compute_msa_consensus(aligned_sequence, qual)).seq AS consensus
FROM consensus_input
GROUP BY bin_id;
```

**See also:** [`extract_linked_amplicon`](scalar-functions.md#extract_linked_ampliconseq-qual-anchor5-anchor3-min_len-max_len-error_rate), [`match_short_barcodes`](table-functions.md#match_short_barcodesquery_table-ref_table-max_nmn-report_alltrue), [`compute_pileup`](table-functions.md#compute_pileupalignments_table-reference_table), [`align_mafft`](table-functions.md#align_maffttable_name-sample_idcol) — together these four primitives compose into a long-read UMI consensus pipeline.

## `genome_coverage(alignments, subject_total_length, subject_genome_id)`

Table macro that computes genome coverage from alignment data. It compresses overlapping alignment intervals per reference contig using `compress_intervals`, maps contigs to genomes, sums covered bases, and joins with total genome lengths to compute the proportion covered.

**Parameters:**

All three parameters are unquoted table/view names (not string literals):

- `alignments`: A relation with columns `reference` (VARCHAR), `position` (BIGINT), `stop_position` (BIGINT)
- `subject_total_length`: A relation with columns `genome_id` (VARCHAR), `total_length` (BIGINT)
- `subject_genome_id`: A relation with columns `contig_id` (VARCHAR), `genome_id` (VARCHAR)

**Returns:**
| Column | Type | Description |
|--------|------|-------------|
| `genome_id` | VARCHAR | Genome identifier |
| `covered` | BIGINT | Total number of bases covered |
| `proportion_covered` | DOUBLE | Fraction of genome covered (`covered / total_length`) |

**Behavior:**
- Overlapping alignments on the same contig are merged before counting (via `compress_intervals`)
- Multiple contigs mapping to the same genome have their coverage summed
- Contigs in `alignments` that are not present in `subject_genome_id` are excluded from output
- Uses half-open coordinates consistent with `read_alignments` output

**Examples:**
```sql
-- Setup: alignment data and genome metadata
CREATE TABLE alignments AS
  SELECT reference, position, stop_position
  FROM read_alignments('alignments.bam');

CREATE TABLE genome_lengths (genome_id VARCHAR, total_length BIGINT);
INSERT INTO genome_lengths VALUES ('genomeA', 100000), ('genomeB', 200000);

CREATE TABLE contig_to_genome (contig_id VARCHAR, genome_id VARCHAR);
INSERT INTO contig_to_genome VALUES
  ('contig1', 'genomeA'), ('contig2', 'genomeA'),
  ('contig3', 'genomeB');

-- Compute genome coverage
SELECT * FROM genome_coverage(alignments, genome_lengths, contig_to_genome);

-- Filter to genomes with >50% coverage
SELECT * FROM genome_coverage(alignments, genome_lengths, contig_to_genome)
WHERE proportion_covered > 0.5;
```

## Pairwise Alignment Functions

Gap-affine pairwise sequence alignment powered by [WFA2-lib](https://github.com/smarco/WFA2-lib) (Wavefront Alignment Algorithm). Three functions at increasing detail levels:

- `align_pairwise_score` -- alignment score only (fastest)
- `align_pairwise_cigar` -- score + extended CIGAR string
- `align_pairwise_full` -- score + CIGAR + aligned sequences with gap characters

Each function has a 2-argument form (uses defaults) and a 6-argument form (explicit penalties):

```sql
-- 2-arg: uses default penalties (mismatch=4, gap_open=6, gap_extend=2)
SELECT align_pairwise_score(query, subject);

-- 6-arg: custom penalties
SELECT align_pairwise_score(query, subject, 'wfa2', 2, 6, 2);
```

**Parameters (6-arg form):**
- `query` (VARCHAR): Query sequence
- `subject` (VARCHAR): Subject sequence
- `method` (VARCHAR): Alignment method, currently only `'wfa2'`
- `mismatch` (INTEGER): Mismatch penalty (must be > 0)
- `gap_open` (INTEGER): Gap opening penalty (must be >= 0)
- `gap_extend` (INTEGER): Gap extension penalty (must be > 0)

Penalty parameters must be constant values, not column references.

**Returns:**
- `align_pairwise_score` -> `INTEGER` -- alignment score (0 = identical, higher = more divergent)
- `align_pairwise_cigar` -> `STRUCT(score INTEGER, cigar VARCHAR)` -- score and extended CIGAR (`=`/`X` ops, not `M`)
- `align_pairwise_full` -> `STRUCT(score INTEGER, cigar VARCHAR, query_aligned VARCHAR, subject_aligned VARCHAR)` -- score, CIGAR, and aligned sequences with `-` gap characters

NULL inputs produce NULL output. Alignment failure (e.g., excessive divergence) also produces NULL.

**Examples:**
```sql
-- Score identical sequences
SELECT align_pairwise_score('ACGT', 'ACGT');
-- 0

-- Score with single mismatch (default mismatch penalty = 4)
SELECT align_pairwise_score('ACGT', 'ACAT');
-- 4

-- Get CIGAR string
SELECT (align_pairwise_cigar('ACGT', 'ACAT')).cigar;
-- 2=1X1=

-- Get full alignment with gap characters
SELECT (align_pairwise_full('ACGT', 'AGT')).query_aligned,
       (align_pairwise_full('ACGT', 'AGT')).subject_aligned;
-- Shows aligned sequences with '-' for gaps

-- Use with table data
SELECT name,
       align_pairwise_score(query_seq, ref_seq) AS score,
       (align_pairwise_cigar(query_seq, ref_seq)).cigar AS cigar
FROM sequence_pairs;

-- Custom penalties for more sensitive alignment
SELECT (align_pairwise_full(seq_a, seq_b, 'wfa2', 2, 4, 1)).query_aligned
FROM amplicon_pairs;
```

## Utility Functions

## `formula(formula_string)`

Compute the monoisotopic mass of a chemical formula.

**Parameters:**
- `formula_string` (VARCHAR): Chemical formula (e.g., `'C6H12O6'`, `'Fe'`, `'H2O'`)

**Returns:** DOUBLE - Monoisotopic mass in Daltons.

**Examples:**
```sql
SELECT formula('H2O');       -- 18.010565
SELECT formula('C6H12O6');   -- 180.063388
SELECT formula('Fe');        -- 55.934936
```

## `massql(query, source)`

Execute a [MassQL](https://github.com/mwang87/MassQueryLanguage) query against mass spectrometry data. Parses the MassQL query string, transpiles it to SQL, and executes it against the source table or file.

**Parameters:**
- `query` (VARCHAR): MassQL query string
- `source` (VARCHAR): Table name or file path (`.mzML` files are detected automatically)

**Returns:** Table with columns depending on the aggregation function used.

See [Mass Spectrometry & MassQL](massql.md) for full syntax reference and examples.

**Quick examples:**
```sql
-- Find MS2 spectra with a product ion at m/z 167.0857
SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857', 'my_table');

-- Find MS1 spectra with a peak at m/z 200 within 5 ppm
SELECT * FROM massql('QUERY scannum(MS1DATA) WHERE MS1MZ=200:TOLERANCEPPM=5', 'sample.mzML');

-- Iron isotope pattern matching
SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=X AND MS2PROD=2*(X-formula(Fe))', 'my_table');
```

## Utility Functions

### `miint_version()`

Returns the MIINT extension version string.

**Returns:** VARCHAR - Version string derived from the git tag at build time (e.g., `v0.1.0`, or `v0.1.0-3-gabcdef1` if built from a commit after a tag).

**Example:**
```sql
SELECT miint_version();
```
