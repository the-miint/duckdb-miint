# Analysis Functions

Functions for higher-level genomic analysis, sequence manipulation, and pairwise alignment.

## Table of Contents

- [`woltka_ogu_per_sample`](#woltka_ogu_per_samplerelation-sample_id_field-sequence_id_field) - Multi-sample OGU counts
- [`woltka_ogu`](#woltka_ogurelation-sequence_id_field) - Single-sample OGU counts
- [`sequence_dna_reverse_complement` / `sequence_rna_reverse_complement`](#sequence_dna_reverse_complementsequence-and-sequence_rna_reverse_complementsequence) - Reverse complement
- [`sequence_dna_as_regexp` / `sequence_rna_as_regexp`](#sequence_dna_as_regexpsequence-and-sequence_rna_as_regexpsequence) - IUPAC to regex
- [`compress_intervals`](#compress_intervalsstart-stop) - Merge overlapping intervals
- [`genome_coverage`](#genome_coveragealignments-subject_total_length-subject_genome_id) - Compute genome coverage from alignments
- [Pairwise Alignment Functions](#pairwise-alignment-functions) - WFA2-based pairwise alignment
- [`formula`](#formulaformula_string) - Chemical formula to monoisotopic mass
- [`massql`](#massqlquery-source) - MassQL query language for mass spectrometry
- [Utility Functions](#utility-functions) - `miint_version()` and others

## `woltka_ogu_per_sample(relation, sample_id_field, sequence_id_field)`

Compute [Woltka](https://github.com/qiyunzhu/woltka) OGU (Operational Genomic Unit) counts over SAM-like alignment data for multiple samples. This function implements Woltka's classification algorithm, which assigns reads to taxonomic units while accounting for multi-mapped reads.

**IMPORTANT**: Function parameters should NOT be quoted. Pass table and column names directly as identifiers, not as string literals.

**Parameters:**
- `relation`: A table, view, or subquery containing SAM-like alignment data
- `sample_id_field`: Column name containing sample identifiers
- `sequence_id_field`: Column name containing sequence identifiers (can be `read_id` or a numeric index for better performance)

**Required columns in relation:**
- Column specified by `sample_id_field`: Sample identifier
- Column specified by `sequence_id_field`: Read/sequence identifier
- `reference` (VARCHAR): Reference sequence name (feature ID)
- `flags` (USMALLINT): SAM alignment flags

**Returns:**
- `sample_id`: Sample identifier
- `feature_id`: Reference/feature identifier
- `value`: OGU count (fractional, accounts for multi-mapping)

**Algorithm:**
1. Orients reads using alignment flags (forward/reverse)
2. For each read orientation, divides 1 by the number of unique features aligned to
3. Aggregates fractional counts per sample and feature

**Examples:**
```sql
-- Basic usage: count OGUs per sample from alignment table
SELECT * FROM woltka_ogu_per_sample(
    my_alignments,
    sample_id,
    read_id
);

-- Using with a filtered view
CREATE VIEW high_quality_alignments AS
    SELECT *, 'sample1' AS sample_id
    FROM read_alignments('alignments.bam')
    WHERE mapq >= 30 AND alignment_is_primary(flags);

SELECT * FROM woltka_ogu_per_sample(
    high_quality_alignments,
    sample_id,
    read_id
);

-- Process multiple samples with UNION
CREATE VIEW all_samples AS
    SELECT *, 'sample1' AS sample_id FROM read_alignments('sample1.bam')
    UNION ALL
    SELECT *, 'sample2' AS sample_id FROM read_alignments('sample2.bam')
    UNION ALL
    SELECT *, 'sample3' AS sample_id FROM read_alignments('sample3.bam');

SELECT * FROM woltka_ogu_per_sample(
    all_samples,
    sample_id,
    read_id
) ORDER BY sample_id, feature_id;

-- Export to BIOM format for downstream analysis
COPY (
    SELECT * FROM woltka_ogu_per_sample(my_alignments, sample_id, read_id)
) TO 'ogu_table.biom' (FORMAT BIOM, COMPRESSION 'gzip');

-- WRONG: Do not quote parameters (this will cause an error)
-- SELECT * FROM woltka_ogu_per_sample('my_alignments', 'sample_id', 'read_id');
```

## `woltka_ogu(relation, sequence_id_field)`

Compute Woltka OGU counts over SAM-like alignment data for a single sample or aggregated across all samples.

**IMPORTANT**: Function parameters should NOT be quoted. Pass table and column names directly as identifiers, not as string literals.

**Parameters:**
- `relation`: A table, view, or subquery containing SAM-like alignment data
- `sequence_id_field`: Column name containing sequence identifiers (can be `read_id` or a numeric index for better performance)

**Required columns in relation:**
- Column specified by `sequence_id_field`: Read/sequence identifier
- `reference` (VARCHAR): Reference sequence name (feature ID)
- `flags` (USMALLINT): SAM alignment flags

**Returns:**
- `feature_id`: Reference/feature identifier
- `value`: OGU count (fractional, accounts for multi-mapping)

**Examples:**
```sql
-- Basic usage: count OGUs from alignment table
SELECT * FROM woltka_ogu(
    my_alignments,
    read_id
);

-- Direct query from read_alignments
SELECT * FROM woltka_ogu(
    (SELECT * FROM read_alignments('alignments.bam')),
    read_id
) ORDER BY value DESC;

-- Filter high-quality primary alignments
CREATE VIEW primary_alignments AS
    SELECT * FROM read_alignments('alignments.bam')
    WHERE alignment_is_primary(flags) AND mapq >= 20;

SELECT * FROM woltka_ogu(primary_alignments, read_id);

-- Combine with filtering for specific references
CREATE VIEW bacterial_alignments AS
    SELECT * FROM read_alignments('metagenome.bam')
    WHERE reference LIKE 'bacteria_%';

SELECT feature_id, value
FROM woltka_ogu(bacterial_alignments, read_id)
WHERE value > 10
ORDER BY value DESC;

-- Export to BIOM format (add sample_id column)
COPY (
    SELECT feature_id, 'MySample' AS sample_id, value
    FROM woltka_ogu(my_alignments, read_id)
) TO 'ogu_single.biom' (FORMAT BIOM);
```

**Notes:**
- Multi-mapped reads (reads aligning to multiple references) are fractionally assigned: each mapping receives weight 1/N where N is the number of unique references
- Read orientation (forward/reverse) is considered separately using SAM flags
- For better performance with large datasets, consider adding a numeric index column and using it as `sequence_id_field` instead of `read_id`
- This function handles paired-end data by distinguishing R1 and R2 reads via the `alignment_is_read1()` flag
- **Implementation note:** Implemented as a DuckDB macro (table-returning expression), so parameters are not quoted

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
