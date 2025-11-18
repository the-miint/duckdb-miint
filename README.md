# MIINT: MIcrobiome INTelligence.

This extension, Miint, allow you to ... <extension_goal>.

## Installing
To install MIINT the extension binary, DuckDB should be launched with the
`allow_unsigned_extensions` option set to true. How to set this will depend on the client you're using. Some examples:

CLI:
```shell
duckdb -unsigned
```

Python:
```python
con = duckdb.connect(':memory:', config={'allow_unsigned_extensions' : 'true'})
```

NodeJS:
```js
db = new duckdb.Database(':memory:', {"allow_unsigned_extensions": "true"});
```

After running these steps, you can install and load your extension using the regular INSTALL/LOAD commands in DuckDB:
```sql
LOAD '/path/to/miint.duckdb_extension';
```

## Building
### Managing dependencies
DuckDB extensions uses VCPKG for dependency management. Enabling VCPKG is very simple: follow the [installation instructions](https://vcpkg.io/en/getting-started) or just run the following:
```shell
git clone https://github.com/Microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh
export VCPKG_TOOLCHAIN_PATH=`pwd`/vcpkg/scripts/buildsystems/vcpkg.cmake
```

### Build system
We use [Ninja](https://ninja-build.org) for quick builds. The easiest install is to use prebuild [release](https://github.com/ninja-build/ninja/releases) binaries.

### Build steps
Now to build the extension, run:
```sh
GEN=ninja make
```

The main binaries that will be built are:

```sh
./build/release/duckdb
./build/release/test/unittest
./build/release/extension/miint/miint.duckdb_extension
```

- `duckdb` is the binary for the duckdb shell with the extension code automatically loaded.
- `unittest` is the test runner of duckdb. Again, the extension is already linked into the binary.
- `miint.duckdb_extension` is the loadable binary as it would be distributed.

## Running the extension
To run the extension code, simply start the shell with `./build/release/duckdb`.

## Functions

### `read_alignments(filename, [reference_lengths='table_name'], [include_filepath=false])`
Read SAM/BAM alignment files.

**Note:** `read_sam` is still supported as a backward-compatible alias.

**Parameters:**
- `filename` (VARCHAR or VARCHAR[]): Path to SAM/BAM file(s)
- `reference_lengths` (VARCHAR, optional): Table name containing reference sequences for headerless SAM files. Table must have at least 2 columns: first column = reference name (VARCHAR), second column = reference length (INTEGER/BIGINT). Column names don't matter.
- `include_filepath` (BOOLEAN, optional, default false): Add filepath column to output

**Output schema includes:**
- `position` (BIGINT): 1-based start position
- `stop_position` (BIGINT): 1-based stop position (computed from CIGAR using `bam_endpos`)
- `cigar` (VARCHAR): CIGAR string
- Plus other standard SAM fields and optional tags

**Examples:**
```sql
-- Read SAM file with header
SELECT * FROM read_alignments('alignments.sam');

-- Read BAM file
SELECT * FROM read_alignments('alignments.bam');

-- Read headerless SAM file (requires reference table)
CREATE TABLE my_refs AS
  SELECT 'chr1' AS name, 248956422 AS length
  UNION ALL SELECT 'chr2', 242193529;
SELECT * FROM read_alignments('headerless.sam', reference_lengths='my_refs');

-- Read multiple files with filepath tracking
SELECT * FROM read_alignments(['file1.sam', 'file2.bam'], include_filepath=true);

-- Backward compatible: read_sam still works
SELECT * FROM read_sam('alignments.sam');
```

### `read_fastx(filename)`
Read FASTA/FASTQ sequence files.

### `read_biom(filename, [include_filepath=false])`
Read BIOM (Biological Observation Matrix) format files.

**Parameters:**
- `filename` (VARCHAR or VARCHAR[]): Path to BIOM file(s)
- `include_filepath` (BOOLEAN, optional, default false): Add filepath column to output

**Output schema:**
- `sample_id` (VARCHAR): Sample identifier
- `feature_id` (VARCHAR): Feature/OGU/OTU/ASV identifier
- `value` (DOUBLE): Abundance/count value
- `filepath` (VARCHAR, optional): File path when include_filepath=true

**Behavior:**
- Reads BIOM format v2.1 files (HDF5-based)
- Returns data in sparse COO (coordinate) format: one row per non-zero (sample, feature, value) entry
- Supports parallel processing for faster reads of large files
- Zero values are not returned (sparse representation)
- Supports reading multiple files which are concatenated in the output

**Examples:**
```sql
-- Read single BIOM file
SELECT * FROM read_biom('ogu_table.biom');

-- Read multiple BIOM files
SELECT * FROM read_biom(['sample1.biom', 'sample2.biom', 'sample3.biom']);

-- Include source filepath for each record
SELECT * FROM read_biom('ogu_table.biom', include_filepath=true);

-- Aggregate counts per sample
SELECT sample_id, SUM(value) as total_count
FROM read_biom('ogu_table.biom')
GROUP BY sample_id
ORDER BY total_count DESC;

-- Aggregate counts per feature
SELECT feature_id, SUM(value) as total_abundance
FROM read_biom('ogu_table.biom')
GROUP BY feature_id
ORDER BY total_abundance DESC
LIMIT 10;

-- Filter for specific samples
SELECT feature_id, value
FROM read_biom('ogu_table.biom')
WHERE sample_id IN ('Sample1', 'Sample2', 'Sample3')
ORDER BY value DESC;

-- Count unique features per sample
SELECT sample_id, COUNT(DISTINCT feature_id) as n_features
FROM read_biom('ogu_table.biom')
GROUP BY sample_id;

-- Join with metadata table
CREATE TABLE sample_metadata AS
SELECT 'Sample1' as sample_id, 'Control' as group
UNION ALL SELECT 'Sample2', 'Treatment'
UNION ALL SELECT 'Sample3', 'Treatment';

SELECT sm.group, b.feature_id, SUM(b.value) as group_total
FROM read_biom('ogu_table.biom') b
JOIN sample_metadata sm ON b.sample_id = sm.sample_id
GROUP BY sm.group, b.feature_id
ORDER BY sm.group, group_total DESC;

-- Round-trip: read BIOM, filter, write back to BIOM
COPY (
    SELECT sample_id, feature_id, value
    FROM read_biom('input.biom')
    WHERE value >= 10.0  -- Filter low abundance features
) TO 'filtered.biom' (FORMAT BIOM, COMPRESSION 'gzip');

-- Convert BIOM to CSV
COPY (
    SELECT * FROM read_biom('ogu_table.biom')
    ORDER BY sample_id, feature_id
) TO 'ogu_table.csv' (HEADER, DELIMITER ',');

-- Read multiple files and track source
SELECT filepath, sample_id, COUNT(*) as n_features, SUM(value) as total_count
FROM read_biom(['batch1.biom', 'batch2.biom', 'batch3.biom'], include_filepath=true)
GROUP BY filepath, sample_id
ORDER BY filepath, sample_id;
```

**Compatibility:**
- Reads BIOM files created by:
  - QIIME2 (FeatureTable exports)
  - biom-format Python package
  - Woltka
  - This extension's `COPY ... FORMAT BIOM`
- Follows BIOM format specification v2.1
- Only supports HDF5-based BIOM files (not JSON format)

**Performance:**
- Efficiently handles large sparse matrices
- Parallel processing enabled by default
- Only non-zero values are read/returned, optimizing memory usage

### SAM Flag Functions

Test individual SAM flag bits. Each function takes a `USMALLINT` (the flags column from `read_sam`) and returns a `BOOLEAN`.

**Primary names:**
- `alignment_is_paired(flags)` - Read is paired (0x1)
- `alignment_is_proper_pair(flags)` - Read is properly paired (0x2)
- `alignment_is_unmapped(flags)` - Read is unmapped (0x4)
- `alignment_is_mate_unmapped(flags)` - Mate is unmapped (0x8)
- `alignment_is_reverse(flags)` - Read is reverse strand (0x10)
- `alignment_is_mate_reverse(flags)` - Mate is reverse strand (0x20)
- `alignment_is_read1(flags)` - Read is first in pair (0x40)
- `alignment_is_read2(flags)` - Read is second in pair (0x80)
- `alignment_is_secondary(flags)` - Secondary alignment (0x100)
- `alignment_is_primary(flags)` - NOT alignment_is_secondary
- `alignment_is_qc_failed(flags)` - QC failure (0x200)
- `alignment_is_duplicate(flags)` - PCR/optical duplicate (0x400)
- `alignment_is_supplementary(flags)` - Supplementary alignment (0x800)

**HTSlib-compatible aliases:**
`is_paired`, `is_proper_pair`, `is_unmapped`, `is_munmap`, `is_reverse`, `is_mreverse`, `is_read1`, `is_read2`, `is_secondary`, `is_qcfail`, `is_dup`, `is_supplementary`

**Example:**
```sql
SELECT read_id, flags
FROM read_alignments('alignments.sam')
WHERE alignment_is_paired(flags)
  AND NOT alignment_is_unmapped(flags);
```

### `alignment_seq_identity(cigar, nm, md, type)`

Calculate sequence identity between read and reference using three different methods. They are derived from Heng Li's [blog post](https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity).

*IMPORTANT*: these calculations assume use of '='/'X' (e.g., `--xeq` with bowtie2, `-eqx` with minimap2, etc).

**Parameters:**
- `cigar` (VARCHAR): CIGAR string from alignment.
- `nm` (BIGINT): NM tag value (edit distance)
- `md` (VARCHAR): MD tag value (mismatch positions)
- `type` (VARCHAR, default='gap_compressed'): Identity calculation method

**Identity types:**

1. **`'gap_excluded'`**: Ignore gaps, only consider match/mismatch positions
   - Formula: `#matches / (#matches + #mismatches)`
   - Requires: CIGAR + MD tag
   - Use case: Genetic divergence between species

2. **`'blast'`**: Traditional BLAST-style identity
   - Formula: `#matches / alignment_columns`
   - Requires: CIGAR + NM tag
   - Use case: General similarity measurement
   - Note: Large indels significantly lower identity

3. **`'gap_compressed'`** (default): Count consecutive gaps as single events
   - Formula: `(m - n + g) / (m + o)` where m=M_columns, n=NM, g=gap_bases, o=gap_opens
   - Equivalent to: `1 - (n - g + o) / (m + o)` from the blog post
   - Requires: CIGAR + NM tag
   - Use case: Filtering alignments (recommended)
   - Note: More robust to structural variations

**Returns:** DOUBLE between 0.0 and 1.0, or NULL for unmapped reads or missing required tags

**Example:**
```sql
-- Calculate gap-compressed identity (default)
SELECT read_id, alignment_seq_identity(cigar, tag_nm, tag_md, 'gap_compressed') AS identity
FROM read_alignments('alignments.sam')
WHERE tag_nm IS NOT NULL;

-- Filter high-quality alignments
SELECT COUNT(*)
FROM read_alignments('alignments.bam')
WHERE alignment_seq_identity(cigar, tag_nm, tag_md) > 0.95;

-- Compare different identity methods
SELECT read_id,
  alignment_seq_identity(cigar, tag_nm, tag_md, 'gap_excluded') AS gap_excl,
  alignment_seq_identity(cigar, tag_nm, tag_md, 'blast') AS blast,
  alignment_seq_identity(cigar, tag_nm, tag_md, 'gap_compressed') AS gap_comp
FROM read_alignments('alignments.sam')
WHERE tag_nm IS NOT NULL AND tag_md IS NOT NULL;
```

**Reference:** [On the definition of sequence identity](https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity) by Heng Li

### `alignment_query_length(cigar, [include_hard_clips=true])`

Calculate the total query length from a CIGAR string. This is useful for understanding read lengths and query coverage.

**Parameters:**
- `cigar` (VARCHAR): CIGAR string from alignment
- `include_hard_clips` (BOOLEAN, default=true): Whether to include hard-clipped bases in the length

**Returns:** BIGINT - Total query length, or 0 for NULL/unmapped reads

**Behavior:**
- When `include_hard_clips=true`: Returns M + I + S + = + X + H (all query-consuming operations)
- When `include_hard_clips=false`: Returns M + I + S + = + X (matches HTSlib's `bam_cigar2qlen`)
- Soft clips (S) are always included (they're present in the sequence field)
- Hard clips (H) are only included when parameter is true
- Deletions (D) and reference skips (N) don't consume query, so they're never counted

**Examples:**
```sql
-- Get query length including hard clips (default)
SELECT read_id, alignment_query_length(cigar) AS query_len
FROM read_alignments('alignments.sam');

-- Get query length excluding hard clips (matches bam_cigar2qlen)
SELECT read_id, alignment_query_length(cigar, false) AS query_len
FROM read_alignments('alignments.sam');

-- Compare lengths with and without hard clips
SELECT read_id, cigar,
  alignment_query_length(cigar, true) AS len_with_hard,
  alignment_query_length(cigar, false) AS len_without_hard
FROM read_alignments('alignments.sam')
WHERE cigar LIKE '%H%';

-- Calculate average query length per reference
SELECT reference, AVG(alignment_query_length(cigar)) AS avg_query_len
FROM read_alignments('alignments.bam')
WHERE NOT alignment_is_unmapped(flags)
GROUP BY reference;
```

**Note:** When `include_hard_clips=false`, this function's output matches HTSlib's `bam_cigar2qlen` behavior, which counts M, I, S, =, and X operations.

### `alignment_query_coverage(cigar, [type='aligned'])`

Calculate the proportion of query bases covered by the reference alignment. This helps assess how much of a read actually aligns versus being clipped.

**Parameters:**
- `cigar` (VARCHAR): CIGAR string from alignment
- `type` (VARCHAR, default='aligned'): Coverage calculation method

**Coverage types:**

1. **`'aligned'`** (default): Bases that align to the reference
   - Formula: `(M + = + X) / (M + I + S + = + X + H)`
   - Only counts bases that match/mismatch the reference
   - Insertions and clips reduce coverage

2. **`'mapped'`**: Bases that are mapped (not clipped)
   - Formula: `(M + I + = + X) / (M + I + S + = + X + H)`
   - Counts insertions as "mapped" even though they don't align
   - Only clips reduce coverage

**Returns:** DOUBLE between 0.0 and 1.0, or NULL for NULL CIGAR

**Behavior:**
- Query length denominator always includes hard clips
- Returns 0.0 for reads with only clipping operations (no alignment)
- Deletions (D) and reference skips (N) don't affect coverage (they don't consume query)

**Examples:**
```sql
-- Get aligned coverage (default)
SELECT read_id, alignment_query_coverage(cigar) AS aligned_cov
FROM read_alignments('alignments.sam');

-- Get mapped coverage (includes insertions)
SELECT read_id, alignment_query_coverage(cigar, 'mapped') AS mapped_cov
FROM read_alignments('alignments.sam');

-- Compare aligned vs mapped coverage
SELECT read_id, cigar,
  alignment_query_coverage(cigar, 'aligned') AS aligned_cov,
  alignment_query_coverage(cigar, 'mapped') AS mapped_cov
FROM read_alignments('alignments.sam')
WHERE cigar LIKE '%I%';  -- Reads with insertions show the difference

-- Filter reads with high query coverage
SELECT COUNT(*)
FROM read_alignments('alignments.bam')
WHERE alignment_query_coverage(cigar, 'aligned') > 0.9;

-- Find heavily clipped reads
SELECT read_id, cigar, alignment_query_coverage(cigar) AS coverage
FROM read_alignments('alignments.sam')
WHERE alignment_query_coverage(cigar) < 0.5
ORDER BY coverage;

-- Calculate average coverage per reference
SELECT reference,
  AVG(alignment_query_coverage(cigar, 'aligned')) AS avg_aligned_cov,
  AVG(alignment_query_coverage(cigar, 'mapped')) AS avg_mapped_cov
FROM read_alignments('alignments.bam')
WHERE NOT alignment_is_unmapped(flags)
GROUP BY reference;
```

**Use cases:**
- **Aligned coverage**: Assess quality of alignment (how much of read actually matches reference)
- **Mapped coverage**: Identify heavily clipped reads (adapters, chimeras, low-quality ends)
- Filter reads based on alignment quality thresholds
- QC metrics for sequencing runs

## Analysis Macros

### `woltka_ogu_per_sample(relation, sample_id_field, sequence_id_field)`

Compute [Woltka](https://github.com/qiyunzhu/woltka) OGU (Operational Genomic Unit) counts over SAM-like alignment data for multiple samples. This macro implements Woltka's classification algorithm, which assigns reads to taxonomic units while accounting for multi-mapped reads.

**IMPORTANT**: Macro parameters should NOT be quoted. Pass table and column names directly as identifiers, not as string literals.

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

### `woltka_ogu(relation, sequence_id_field)`

Compute Woltka OGU counts over SAM-like alignment data for a single sample or aggregated across all samples.

**IMPORTANT**: Macro parameters should NOT be quoted. Pass table and column names directly as identifiers, not as string literals.

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
- The macro handles paired-end data by distinguishing R1 and R2 reads via the `alignment_is_read1()` flag

### `sequence_dna_reverse_complement(sequence)` and `sequence_rna_reverse_complement(sequence)`

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
- **DNA**: A↔T, C↔G, plus IUPAC codes (R↔Y, S↔S, W↔W, K↔M, B↔V, D↔H, N↔N)
- **RNA**: A↔U, C↔G, plus IUPAC codes (R↔Y, S↔S, W↔W, K↔M, B↔V, D↔H, N↔N)

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

### `compress_intervals(start, stop)`

Aggregate function that merges overlapping genomic intervals into a minimal set of non-overlapping intervals. Useful for computing coverage regions, reducing redundant intervals, and analyzing read depth.

**Parameters:**
- `start` (BIGINT): Start position of interval
- `stop` (BIGINT): Stop position of interval

**Returns:** `LIST<STRUCT(start BIGINT, stop BIGINT)>` - Array of merged intervals, sorted by start position

**Behavior:**
- Intervals that overlap or touch (stop₁ = start₂) are merged into a single interval
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

## COPY Formats

DuckDB miint provides custom COPY formats for writing bioinformatics file formats.

### `COPY ... TO '...' (FORMAT FASTQ)`

Write query results to FASTQ format files. Requires `read_id`, `sequence1`, and `qual1` columns from `read_fastx` output.

**Required columns:**
- `read_id` (VARCHAR): Sequence identifier
- `sequence1` (VARCHAR): DNA/RNA sequence
- `qual1` (BLOB): Quality scores as raw bytes

**Optional columns:**
- `comment` (VARCHAR): Comment line (only included if `INCLUDE_COMMENT=true`)
- `sequence_index` (BIGINT): Used as identifier if `ID_AS_SEQUENCE_INDEX=true`
- `sequence2` (VARCHAR): Second read for paired-end data
- `qual2` (BLOB): Quality scores for second read

**Parameters:**
- `QUAL_OFFSET` (default: 33): Quality score encoding offset (33 or 64)
- `INCLUDE_COMMENT` (default: false): Include comment field in output
- `ID_AS_SEQUENCE_INDEX` (default: false): Use `sequence_index` as identifier instead of `read_id`
- `INTERLEAVE` (default: false): Write paired reads interleaved in single file
- `COMPRESSION` (default: auto): Enable gzip compression (auto-detected from `.gz` extension)

**Examples:**
```sql
-- Basic single-end FASTQ output
COPY (SELECT * FROM read_fastx('input.fastq'))
TO 'output.fastq' (FORMAT FASTQ);

-- Paired-end interleaved output
COPY (SELECT * FROM read_fastx('R1.fastq', 'R2.fastq'))
TO 'output.fastq' (FORMAT FASTQ, INTERLEAVE true);

-- Paired-end split files (use {ORIENTATION} placeholder)
COPY (SELECT * FROM read_fastx('R1.fastq', 'R2.fastq'))
TO 'output_{ORIENTATION}.fastq' (FORMAT FASTQ);

-- Compressed output with custom quality offset
COPY (SELECT * FROM read_fastx('input.fastq'))
TO 'output.fastq.gz' (FORMAT FASTQ, QUAL_OFFSET 33, COMPRESSION gzip);

-- Use sequence index as identifier
COPY (SELECT * FROM read_fastx('input.fastq'))
TO 'output.fastq' (FORMAT FASTQ, ID_AS_SEQUENCE_INDEX true);
```

### `COPY ... TO '...' (FORMAT FASTA)`

Write query results to FASTA format files. Requires `read_id` and `sequence1` columns from `read_fastx` output.

**Required columns:**
- `read_id` (VARCHAR): Sequence identifier
- `sequence1` (VARCHAR): DNA/RNA/protein sequence

**Optional columns:**
- `comment` (VARCHAR): Comment line (only included if `INCLUDE_COMMENT=true`)
- `sequence_index` (BIGINT): Used as identifier if `ID_AS_SEQUENCE_INDEX=true`
- `sequence2` (VARCHAR): Second read for paired-end data

**Parameters:**
- `INCLUDE_COMMENT` (default: false): Include comment field in output
- `ID_AS_SEQUENCE_INDEX` (default: false): Use `sequence_index` as identifier instead of `read_id`
- `INTERLEAVE` (default: false): Write paired reads interleaved in single file
- `COMPRESSION` (default: auto): Enable gzip compression (auto-detected from `.gz` extension)

**Examples:**
```sql
-- Basic FASTA output
COPY (SELECT * FROM read_fastx('input.fasta'))
TO 'output.fasta' (FORMAT FASTA);

-- Paired-end split files
COPY (SELECT * FROM read_fastx('R1.fasta', 'R2.fasta'))
TO 'output_{ORIENTATION}.fasta' (FORMAT FASTA);

-- Compressed with comments
COPY (SELECT * FROM read_fastx('input.fasta'))
TO 'output.fasta.gz' (FORMAT FASTA, INCLUDE_COMMENT true);
```

### `COPY ... TO '...' (FORMAT SAM)` and `COPY ... TO '...' (FORMAT BAM)`

Write query results to SAM or BAM format files. Requires all mandatory SAM columns from `read_alignments` output.

**Required columns:**
- `read_id` (VARCHAR): Query template name
- `flags` (USMALLINT): Bitwise flags
- `reference` (VARCHAR): Reference sequence name
- `position` (BIGINT): 1-based leftmost mapping position
- `mapq` (UTINYINT): Mapping quality
- `cigar` (VARCHAR): CIGAR string
- `mate_reference` (VARCHAR): Reference name of mate/next read
- `mate_position` (BIGINT): Position of mate/next read
- `template_length` (BIGINT): Observed template length

**Optional columns:**
- `tag_as`, `tag_xs`, `tag_ys`, `tag_xn`, `tag_xm`, `tag_xo`, `tag_xg`, `tag_nm` (BIGINT): Optional integer tags
- `tag_yt`, `tag_md`, `tag_sa` (VARCHAR): Optional string tags

**Parameters:**
- `INCLUDE_HEADER` (default: true): Include header with reference sequences
  - **Note:** BAM format requires `INCLUDE_HEADER=true` (headers are mandatory in BAM files)
- `REFERENCE_LENGTHS` (VARCHAR, required if INCLUDE_HEADER=true): Table name containing reference sequences. Table must have at least 2 columns: first column = reference name (VARCHAR), second column = reference length (INTEGER/BIGINT). Column names don't matter.
- `COMPRESSION` (default: auto, SAM only): Enable gzip compression (auto-detected from `.gz` extension)
- `COMPRESSION_LEVEL` (BAM only): BGZF compression level 0-9 (default: 6). Higher = better compression, slower speed.

**SAM Format Examples:**
```sql
-- Create reference table (recommended for reuse, especially with large reference sets)
CREATE TABLE ref_table AS
  SELECT 'genome1' AS name, 248956422 AS length
  UNION ALL SELECT 'genome2', 242193529;

-- Basic SAM output with header
COPY (SELECT * FROM read_alignments('input.sam'))
TO 'output.sam' (FORMAT SAM, REFERENCE_LENGTHS 'ref_table');

-- Headerless SAM output (no reference lengths needed)
COPY (SELECT * FROM read_alignments('input.sam'))
TO 'output.sam' (FORMAT SAM, INCLUDE_HEADER false);

-- Compressed SAM output with header
COPY (SELECT * FROM read_alignments('input.sam'))
TO 'output.sam.gz' (FORMAT SAM, COMPRESSION gzip, REFERENCE_LENGTHS 'ref_table');

-- Filter and write high-quality alignments (headerless)
COPY (
  SELECT * FROM read_alignments('input.sam')
  WHERE mapq >= 30 AND NOT alignment_is_unmapped(flags)
) TO 'filtered.sam' (FORMAT SAM, INCLUDE_HEADER false);
```

**BAM Format Examples:**
```sql
-- Basic BAM output (always includes header)
COPY (SELECT * FROM read_alignments('input.sam'))
TO 'output.bam' (FORMAT BAM, REFERENCE_LENGTHS 'ref_table');

-- BAM with maximum compression
COPY (SELECT * FROM read_alignments('input.bam'))
TO 'compressed.bam' (FORMAT BAM, COMPRESSION_LEVEL 9, REFERENCE_LENGTHS 'ref_table');

-- BAM with no compression (fastest)
COPY (SELECT * FROM read_alignments('input.bam'))
TO 'uncompressed.bam' (FORMAT BAM, COMPRESSION_LEVEL 0, REFERENCE_LENGTHS 'ref_table');

-- Convert SAM to BAM
COPY (SELECT * FROM read_alignments('input.sam'))
TO 'output.bam' (FORMAT BAM, REFERENCE_LENGTHS 'ref_table');

-- Convert BAM to SAM
COPY (SELECT * FROM read_alignments('input.bam'))
TO 'output.sam' (FORMAT SAM, REFERENCE_LENGTHS 'ref_table');

-- Filter alignments and write to BAM
COPY (
  SELECT * FROM read_alignments('input.bam')
  WHERE mapq >= 30 AND alignment_is_primary(flags)
) TO 'filtered.bam' (FORMAT BAM, REFERENCE_LENGTHS 'ref_table');
```

**Notes:**
- SEQ and QUAL fields are always written as `*` in current implementation
- Reference lengths must be provided explicitly when writing headers - they cannot be inferred from the data
- All optional tags present in the input are preserved in the output
- BAM files always require headers (binary format specification)

### `COPY ... TO '...' (FORMAT BIOM)`

Write query results to BIOM (Biological Observation Matrix) format files. BIOM is an HDF5-based format commonly used for representing OGU/OTU/ASV tables in microbiome analyses.

**Required columns:**
- `feature_id` (VARCHAR): Feature/OGU/OTU/ASV identifier
- `sample_id` (VARCHAR): Sample identifier
- `value` (DOUBLE): Abundance/count value

**Parameters:**
- `COMPRESSION` (default: none): HDF5 internal compression algorithm ('gzip', 'lzf', or 'none')
- `ID` (default: auto-generated): Custom table identifier for BIOM metadata
- `GENERATED_BY` (default: 'DuckDB-MIINT'): Tool/version string for provenance tracking

**Behavior:**
- **Automatic deduplication**: Duplicate (feature_id, sample_id) pairs are automatically summed
- **Sparse optimization**: Zero values are automatically removed from output
- **Ordering**: Feature and sample IDs appear in order of first occurrence in the input data
- **NULL handling**: NULL values in any required column cause an error

**Examples:**
```sql
-- Basic BIOM output from Woltka results
COPY (
    SELECT * FROM woltka_ogu_per_sample(my_alignments, sample_id, read_id)
) TO 'ogu_table.biom' (FORMAT BIOM);

-- With HDF5 gzip compression (recommended)
COPY (
    SELECT * FROM woltka_ogu_per_sample(my_alignments, sample_id, read_id)
) TO 'ogu_table.biom' (FORMAT BIOM, COMPRESSION 'gzip');

-- With custom metadata
COPY (
    SELECT * FROM woltka_ogu_per_sample(my_alignments, sample_id, read_id)
) TO 'ogu_table.biom' (FORMAT BIOM,
                       COMPRESSION 'gzip',
                       ID 'MyStudy_16S',
                       GENERATED_BY 'DuckDB-MIINT v1.0 + Woltka algorithm');

-- Export filtered high-quality alignments
CREATE VIEW high_qual AS
    SELECT *, 'sample1' AS sample_id
    FROM read_alignments('sample1.bam')
    WHERE mapq >= 30 AND alignment_is_primary(flags);

COPY (
    SELECT * FROM woltka_ogu_per_sample(high_qual, sample_id, read_id)
) TO 'high_qual.biom' (FORMAT BIOM, COMPRESSION 'gzip');

-- Convert any compatible table to BIOM format
CREATE TABLE feature_counts AS
SELECT * FROM (VALUES
    ('OGU_001', 'Sample_A', 45.0),
    ('OGU_001', 'Sample_B', 32.0),
    ('OGU_002', 'Sample_A', 18.0),
    ('OGU_002', 'Sample_B', 27.0)
) AS t(feature_id, sample_id, value);

COPY (SELECT * FROM feature_counts)
TO 'counts.biom' (FORMAT BIOM, COMPRESSION 'gzip');

-- Process multiple BAM files and export
CREATE VIEW all_samples AS
    SELECT *, 'sample1' AS sample_id FROM read_alignments('sample1.bam')
    UNION ALL
    SELECT *, 'sample2' AS sample_id FROM read_alignments('sample2.bam')
    UNION ALL
    SELECT *, 'sample3' AS sample_id FROM read_alignments('sample3.bam');

COPY (
    SELECT * FROM woltka_ogu_per_sample(all_samples, sample_id, read_id)
) TO 'all_samples.biom' (FORMAT BIOM, COMPRESSION 'gzip');

-- Export single-sample results (add sample_id column)
COPY (
    SELECT feature_id, 'MySample' AS sample_id, value
    FROM woltka_ogu(my_alignments, read_id)
) TO 'single_sample.biom' (FORMAT BIOM);
```

**Deduplication behavior:**
```sql
-- Duplicates are automatically summed
CREATE TABLE with_dups AS
SELECT * FROM (VALUES
    ('F1', 'S1', 10.0),
    ('F1', 'S1', 5.0),   -- Same feature + sample
    ('F2', 'S1', 8.0)
) AS t(feature_id, sample_id, value);

COPY (SELECT * FROM with_dups)
TO 'merged.biom' (FORMAT BIOM);

-- Result: F1,S1,15.0 and F2,S1,8.0 (duplicates merged)

-- Zeros are removed
CREATE TABLE with_zeros AS
SELECT * FROM (VALUES
    ('F1', 'S1', 10.0),
    ('F1', 'S2', 0.0),   -- Zero value
    ('F2', 'S1', 8.0)
) AS t(feature_id, sample_id, value);

COPY (SELECT * FROM with_zeros)
TO 'sparse.biom' (FORMAT BIOM);

-- Result: Only F1,S1,10.0 and F2,S1,8.0 (zero removed)
```

**Compression options:**
- `'gzip'`: Good compression, widely compatible (recommended)
- `'lzf'`: Faster compression/decompression, less space savings
- `'none'`: No compression, fastest but largest files

**Compatibility:**
- Output files are compatible with:
  - QIIME2 (qiime tools import --type FeatureTable[Frequency])
  - phyloseq (import_biom)
  - biom-format Python package
  - The `read_biom()` function in this extension

**Notes:**
- BIOM format uses HDF5, which provides efficient storage for sparse matrices
- The format automatically handles very large feature/sample counts
- Feature and sample metadata columns are not currently supported (data only)
- The output follows BIOM format specification v2.1

## Running the tests
Different tests can be created for DuckDB extensions. The primary way of testing DuckDB extensions should be the SQL tests in `./test/sql`. These SQL tests can be run using:
```sh
make test
```

---

This repository is based on https://github.com/duckdb/extension-template, check it out if you want to build and ship your own DuckDB extension.
