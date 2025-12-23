# MIINT: MIcrobiome INTelligence.

The MIINT [DuckDB](https://duckdb.org) extension enables DuckDB to interoperate with file formats and operations central to microbiome research.

## Quick Start

MIINT brings the power of SQL to microbiome data analysis. Here's a complete workflow analyzing metagenomic alignments:

```sql
-- Load the extension
-- IMPORTANT: see installing, it is at the moment necessary to allow "unsigned" extensions _before_ loading the extension.
-- This is because We currently are not available in the DuckDB community repository
-- but will be in the future.
LOAD '/path/to/miint.duckdb_extension';

-- Read alignment files and combine samples
CREATE VIEW all_alignments AS
    SELECT *, 'sample1' AS sample_id FROM read_alignments('sample1.bam')
    UNION ALL
    SELECT *, 'sample2' AS sample_id FROM read_alignments('sample2.bam')
    UNION ALL
    SELECT *, 'sample3' AS sample_id FROM read_alignments('sample3.bam');

-- Filter high-quality primary alignments
CREATE VIEW high_quality AS
    SELECT * FROM all_alignments
    WHERE alignment_is_primary(flags)
      AND mapq >= 30
      AND alignment_query_coverage(cigar) > 0.9;

-- Compute OGU counts using Woltka algorithm
CREATE VIEW ogu_counts AS
    SELECT * FROM woltka_ogu_per_sample(high_quality, sample_id, read_id);

-- Export to BIOM format for downstream analysis (QIIME2, phyloseq, etc.)
COPY (SELECT * FROM ogu_counts)
TO 'ogu_table.biom' (FORMAT BIOM, COMPRESSION 'gzip');

-- Or analyze directly in SQL
SELECT sample_id,
       COUNT(DISTINCT feature_id) AS richness,
       SUM(value) AS total_reads
FROM ogu_counts
GROUP BY sample_id;
```

**Key capabilities:**
- **Read bioinformatics formats as tables**: FASTQ/FASTA, SAM/BAM, BIOM
- **Analyze with SQL**: Filter, aggregate, join with metadata tables
- **Write back to standard formats**: Export results to FASTQ, SAM/BAM, BIOM
- **Powerful functions**: Sequence identity, coverage, flag checking, Woltka classification
- **Performance**: Parallel processing, efficient compression, streaming I/O

## Table of Contents

- [Quick Start](#quick-start)
- [Installing](#installing)
- [Building](#building)
  - [Managing dependencies](#managing-dependencies)
  - [Build system](#build-system)
  - [Build steps](#build-steps)
- [Running the extension](#running-the-extension)
- [Functions](#functions)
  - [read_alignments](#read_alignmentsfilename-reference_lengthstable_name-include_filepathfalse)
  - [read_fastx](#read_fastxfilename-sequence2filename-include_filepathfalse-qual_offset33)
  - [read_biom](#read_biomfilename-include_filepathfalse)
  - [read_gff](#read_gffpath)
  - [SAM Flag Functions](#sam-flag-functions)
  - [alignment_seq_identity](#alignment_seq_identitycigar-nm-md-type)
  - [alignment_query_length](#alignment_query_lengthcigar-include_hard_clipstrue)
  - [alignment_query_coverage](#alignment_query_coveragecigar-typealigned)
- [Analysis Functions](#analysis-functions)
  - [woltka_ogu_per_sample](#woltka_ogu_per_samplerelation-sample_id_field-sequence_id_field)
  - [woltka_ogu](#woltka_ogurelation-sequence_id_field)
  - [sequence_dna_reverse_complement / sequence_rna_reverse_complement](#sequence_dna_reverse_complementsequence-and-sequence_rna_reverse_complementsequence)
  - [sequence_dna_as_regexp / sequence_rna_as_regexp](#sequence_dna_as_regexpsequence-and-sequence_rna_as_regexpsequence)
  - [compress_intervals](#compress_intervalsstart-stop)
- [COPY Formats](#copy-formats)
  - [FORMAT FASTQ](#copy--to--format-fastq)
  - [FORMAT FASTA](#copy--to--format-fasta)
  - [FORMAT SAM / FORMAT BAM](#copy--to--format-sam-and-copy--to--format-bam)
  - [FORMAT BIOM](#copy--to--format-biom)
- [Running the tests](#running-the-tests)
  - [SQL Logic Tests](#sql-logic-tests)
  - [C++ Unit Tests](#c-unit-tests)
  - [Test Data](#test-data)

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

### `read_alignments(filename, [reference_lengths='table_name'], [include_filepath=false], [include_seq_qual=false])`
Read SAM/BAM alignment files.

**Note:** `read_sam` is still supported as a backward-compatible alias.

**Parameters:**
- `filename` (VARCHAR or VARCHAR[]): Path to SAM/BAM file(s), or `-` / `/dev/stdin` for standard input
- `reference_lengths` (VARCHAR, optional): Table name containing reference sequences for headerless SAM files. Table must have at least 2 columns: first column = reference name (VARCHAR), second column = reference length (INTEGER/BIGINT). Column names don't matter.
- `include_filepath` (BOOLEAN, optional, default false): Add filepath column to output
- `include_seq_qual` (BOOLEAN, optional, default false): Add sequence and quality score columns to output. When enabled, primary alignments (non-secondary, non-supplementary) and unmapped reads must have SEQ/QUAL data or an error will be raised.

**Output schema includes:**
- `position` (BIGINT): 1-based start position
- `stop_position` (BIGINT): 1-based stop position (computed from CIGAR using `bam_endpos`)
- `cigar` (VARCHAR): CIGAR string
- `sequence` (VARCHAR, optional): Read sequence from SEQ field (when include_seq_qual=true)
- `qual` (UTINYINT[], optional): Quality scores as array of integers 0-93 (when include_seq_qual=true)
- Plus other standard SAM fields and optional tags

**Behavior:**
- Supports both SAM (text) and BAM (binary) formats
- Auto-detects file format from content
- Supports gzip-compressed SAM files
- Supports stdin input using `-` or `/dev/stdin` (single file only, not in arrays)
- Supports parallel processing (4 threads for files, single-threaded for stdin)
- For headerless SAM files, provide reference information via `reference_lengths` parameter

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

-- Include sequence and quality scores
SELECT read_id, sequence, qual
FROM read_alignments('alignments.sam', include_seq_qual=true);

-- Analyze quality scores
SELECT read_id, sequence, list_avg(qual) as avg_qual
FROM read_alignments('alignments.bam', include_seq_qual=true)
WHERE list_avg(qual) >= 30;

-- Filter by sequence length
SELECT read_id, len(sequence) as seq_len
FROM read_alignments('alignments.sam', include_seq_qual=true)
WHERE len(sequence) >= 100;

-- Read from stdin with header
SELECT * FROM read_alignments('/dev/stdin');

-- Read headerless SAM from stdin (requires reference table)
CREATE TABLE refs AS SELECT 'chr1' AS name, 248956422 AS length;
SELECT * FROM read_alignments('-', reference_lengths='refs');

-- Backward compatible: read_sam still works
SELECT * FROM read_sam('alignments.sam');
```

**Stdin limitations:**
- Cannot be used in multi-file arrays (e.g., `['/dev/stdin', 'file.sam']` will error)
- Data with headers works without `reference_lengths`
- Headerless data requires `reference_lengths` parameter
- User must know whether their stdin data contains headers

### `read_fastx(filename, [sequence2=filename], [include_filepath=false], [qual_offset=33])`
Read FASTA/FASTQ sequence files.

**Parameters:**
- `filename` (VARCHAR or VARCHAR[]): Path to FASTA/FASTQ file(s) for single-end reads, or R1 files for paired-end
- `sequence2` (VARCHAR or VARCHAR[], optional): Path to R2 file(s) for paired-end reads. Must have same number of files as `filename`
- `include_filepath` (BOOLEAN, optional, default false): Add filepath column to output
- `qual_offset` (INTEGER, optional, default 33): Quality score offset (33 for Phred+33, 64 for Phred+64)

**Output schema:**
- `sequence_index` (BIGINT): 1-based sequential index per file (resets to 1 for each file when reading multiple files)
- `read_id` (VARCHAR): Sequence identifier (without '@' or '>' prefix)
- `comment` (VARCHAR, nullable): Comment line after identifier
- `sequence1` (VARCHAR): DNA/RNA sequence (R1 for paired-end)
- `sequence2` (VARCHAR, nullable): Second sequence for paired-end reads
- `qual1` (UINT8[], nullable): Quality scores as array of integers (NULL for FASTA)
- `qual2` (UINT8[], nullable): Quality scores for R2 (NULL for FASTA or single-end)
- `filepath` (VARCHAR, optional): File path when include_filepath=true

**Behavior:**
- Auto-detects FASTA (.fasta, .fa, .fna) vs FASTQ (.fastq, .fq) by file extension
- Supports gzip-compressed files (.gz extension)
- Supports stdin input using `-` or `/dev/stdin` (single file only, no paired-end)
- Quality scores converted to integers using specified offset (Phred+33 or Phred+64)
- Supports parallel processing (8 threads for files, 1 thread for stdin)
- For paired-end data, reads are matched by position in files (not by ID)

**Examples:**
```sql
-- Read single-end FASTQ file
SELECT * FROM read_fastx('reads.fastq');

-- Read single-end FASTA file
SELECT * FROM read_fastx('sequences.fasta');

-- Read gzip-compressed FASTQ
SELECT * FROM read_fastx('reads.fastq.gz');

-- Read paired-end FASTQ files
SELECT * FROM read_fastx('R1.fastq', sequence2='R2.fastq');

-- Read multiple single-end files
SELECT * FROM read_fastx(['sample1.fastq', 'sample2.fastq', 'sample3.fastq']);

-- Read multiple paired-end files
SELECT * FROM read_fastx(
    ['sample1_R1.fastq', 'sample2_R1.fastq'],
    sequence2=['sample1_R2.fastq', 'sample2_R2.fastq']
);

-- Include source filepath for tracking (recommended for multiple files)
-- Note: sequence_index resets to 1 for each file
SELECT * FROM read_fastx(['batch1.fastq', 'batch2.fastq'], include_filepath=true)
ORDER BY filepath, sequence_index;

-- Read from stdin
SELECT * FROM read_fastx('-');

-- Specify quality offset for older Illumina data (Phred+64)
SELECT * FROM read_fastx('old_illumina.fastq', qual_offset=64);

-- Get basic statistics
SELECT COUNT(*) as num_reads,
       AVG(LENGTH(sequence1)) as avg_length,
       MIN(LENGTH(sequence1)) as min_length,
       MAX(LENGTH(sequence1)) as max_length
FROM read_fastx('reads.fastq');

-- Filter by sequence length
SELECT * FROM read_fastx('reads.fastq')
WHERE LENGTH(sequence1) >= 100 AND LENGTH(sequence1) <= 150;

-- Count reads per file
SELECT filepath, COUNT(*) as read_count
FROM read_fastx(['file1.fastq', 'file2.fastq', 'file3.fastq'], include_filepath=true)
GROUP BY filepath;

-- Extract sequence IDs
SELECT read_id FROM read_fastx('reads.fastq')
ORDER BY sequence_index;

-- Quality control: check average quality scores
SELECT read_id,
       CAST(AVG(q) AS INTEGER) as avg_qual
FROM (
    SELECT read_id, UNNEST(qual1) as q
    FROM read_fastx('reads.fastq')
)
GROUP BY read_id
HAVING AVG(q) >= 30;
```

**Performance:**
- Multithreaded, but not yet optimized.
- Streaming I/O minimizes memory usage
- Quality scores stored as efficient UINT8 arrays

**Notes:**
- Read IDs must match between R1 and R2 files for paired-end data (not validated, matched by position)
- For FASTA files, `qual1` and `qual2` are NULL
- The `sequence_index` resets to 1 for each file. To distinguish sequences from different files, use `include_filepath=true` and order by `filepath, sequence_index`
- Comment field is NULL if no comment is present in the sequence header

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

### `read_gff(path)`

Read GFF3 (General Feature Format) annotation files. GFF is a standard format for genomic feature annotations including genes, transcripts, exons, and other biological features.

**Parameters:**
- `path` (VARCHAR): Path to GFF/GFF3 file

**Output schema:**
- `seqid` (VARCHAR): Sequence/chromosome identifier
- `source` (VARCHAR): Annotation source (e.g., 'NCBI', 'Ensembl')
- `type` (VARCHAR): Feature type (e.g., 'gene', 'mRNA', 'exon', 'CDS')
- `position` (INTEGER): 1-based start position
- `stop_position` (INTEGER): 1-based end position (inclusive)
- `score` (DOUBLE, nullable): Feature score (NULL if '.')
- `strand` (VARCHAR, nullable): Strand ('+', '-', or NULL if '.')
- `phase` (INTEGER, nullable): CDS phase (0, 1, 2, or NULL if '.')
- `attributes` (MAP(VARCHAR, VARCHAR)): Parsed key-value attributes

**Behavior:**
- Reads tab-delimited GFF3 format files
- Automatically filters comment lines starting with '##'
- Converts '.' to NULL for score, strand, and phase fields
- Parses the attributes column (semicolon-separated key=value pairs) into a SQL MAP
- Supports standard GFF3 specification

**Examples:**
```sql
-- Read a GFF file
SELECT * FROM read_gff('annotations.gff');

-- Extract genes only
SELECT seqid, position, stop_position, attributes['ID'] AS gene_id, attributes['Name'] AS gene_name
FROM read_gff('annotations.gff')
WHERE type = 'gene';

-- Find protein-coding genes
SELECT seqid, position, stop_position, attributes
FROM read_gff('annotations.gff')
WHERE type = 'gene' AND attributes['biotype'] = 'protein_coding';

-- Get feature counts by type
SELECT type, COUNT(*) AS count
FROM read_gff('annotations.gff')
GROUP BY type
ORDER BY count DESC;

-- Find genes on the plus strand
SELECT attributes['ID'] AS gene_id, attributes['Name'] AS gene_name
FROM read_gff('annotations.gff')
WHERE type = 'gene' AND strand = '+';

-- Calculate feature lengths
SELECT type,
       AVG(stop_position - position + 1) AS avg_length,
       MIN(stop_position - position + 1) AS min_length,
       MAX(stop_position - position + 1) AS max_length
FROM read_gff('annotations.gff')
GROUP BY type;

-- Extract CDS features with phase information
SELECT seqid, position, stop_position, phase, attributes['Parent'] AS parent_transcript
FROM read_gff('annotations.gff')
WHERE type = 'CDS' AND phase IS NOT NULL;

-- Find overlapping features between two positions
SELECT type, position, stop_position, attributes['ID'] AS feature_id
FROM read_gff('annotations.gff')
WHERE seqid = 'chr1'
  AND position <= 5000
  AND stop_position >= 1000
ORDER BY position;

-- Access nested attributes
SELECT seqid,
       type,
       attributes['ID'] AS id,
       attributes['Parent'] AS parent,
       attributes['Name'] AS name
FROM read_gff('annotations.gff')
WHERE type = 'exon';

-- Join genes with their transcripts
SELECT g.attributes['ID'] AS gene_id,
       g.attributes['Name'] AS gene_name,
       t.attributes['ID'] AS transcript_id
FROM read_gff('annotations.gff') g
JOIN read_gff('annotations.gff') t
  ON t.attributes['Parent'] = g.attributes['ID']
WHERE g.type = 'gene' AND t.type = 'mRNA';
```

**Attribute Parsing:**
The `attributes` column is automatically parsed from GFF format (semicolon-separated key=value pairs) into a DuckDB MAP:
- Input: `ID=gene1;Name=TEST1;biotype=protein_coding`
- Output: `{'ID': 'gene1', 'Name': 'TEST1', 'biotype': 'protein_coding'}`
- Access values using bracket notation: `attributes['ID']`

**GFF3 Format Notes:**
- Coordinates are 1-based and inclusive (both start and end)
- Strand: '+' (forward), '-' (reverse), '.' (unknown/not applicable)
- Phase: Indicates position within codon (0, 1, or 2) for CDS features
- Comment lines starting with '##' are filtered out (metadata/directives)
- The 9th column (attributes) must contain at least an ID for most feature types

**Use Cases:**
- **Gene annotation analysis**: Extract genes, transcripts, exons from genome annotations
- **Feature filtering**: Select specific feature types or genomic regions
- **Structural analysis**: Calculate feature lengths, count features, analyze distributions
- **Hierarchical queries**: Join parent-child relationships (gene → transcript → exon)
- **Integration**: Combine annotations with alignment data or other genomic datasets

**Compatibility:**
- Supports GFF3 format specification
- Compatible with annotations from:
  - NCBI RefSeq
  - Ensembl
  - GENCODE
  - Other standard GFF3 producers

**Implementation note:** Implemented as a DuckDB macro using `read_csv` with GFF-specific parsing.

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
- `alignment_is_primary(flags)` - Primary alignment (neither secondary nor supplementary)
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

## Analysis Functions

### `woltka_ogu_per_sample(relation, sample_id_field, sequence_id_field)`

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

### `woltka_ogu(relation, sequence_id_field)`

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

### `sequence_dna_as_regexp(sequence)` and `sequence_rna_as_regexp(sequence)`

Convert DNA or RNA sequences with IUPAC ambiguity codes to regular expression patterns. Useful for pattern matching with degenerate primers and probes.

**Parameters:**
- `sequence` (VARCHAR): DNA or RNA sequence string with IUPAC codes

**Returns:** VARCHAR - Regular expression pattern

**Behavior:**
- Unambiguous bases (A, C, G, T/U) remain unchanged
- Ambiguous IUPAC codes expand to character classes (e.g., R → `[AG]`, N → `[ACGT]`)
- Preserves uppercase/lowercase in the output
- Gap characters (`-` and `.`) convert to `.` (regex wildcard matching any character)
- Strict molecular type validation: DNA function rejects U bases, RNA function rejects T bases

**Expansion rules:**
- **Unambiguous bases**: A, C, G, T (DNA) or U (RNA) → no brackets
- **Two-base codes**: R → `[AG]`, Y → `[CT]` or `[CU]`, S → `[CG]`, W → `[AT]` or `[AU]`, K → `[GT]` or `[GU]`, M → `[AC]`
- **Three-base codes**: B → `[CGT]` or `[CGU]`, D → `[AGT]` or `[AGU]`, H → `[ACT]` or `[ACU]`, V → `[ACG]`
- **Any base**: N → `[ACGT]` or `[ACGU]`
- **Gap characters**: `-` or `.` → `.` (matches any character)

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

The MIINT extension uses three complementary testing approaches to ensure correctness and reliability.

### SQL Logic Tests

The primary testing mechanism for user-facing functionality uses DuckDB's SQL logic test framework. Test files are located in `test/sql/` and use the `.test` extension.

**Running SQL tests:**
```sh
# Run all SQL tests
make test

# Run a single test file
./build/release/test/unittest "test/sql/read_alignments.test"

# Run tests matching a pattern
./build/release/test/unittest "[alignment]"
```

**Test file structure:**
SQL logic tests consist of statements and expected outputs. Each test file follows this format:

```sql
# name: test/sql/example.test
# description: Test description
# group: [sql]

require miint

# Test statement that should succeed
statement ok
SELECT * FROM read_fastx('data/fastq/example.fastq');

# Query test with expected output
query I
SELECT COUNT(*) FROM read_fastx('data/fastq/example.fastq');
----
100

# Test error handling
statement error
SELECT * FROM read_fastx('nonexistent.fastq');
----
File not found: nonexistent.fastq
```

**Test organization:**
- `statement ok` - Statement should execute without error
- `statement error` - Statement should raise an error (followed by expected error message)
- `query <types>` - Query should return specific results
  - Types: `I` (INTEGER), `R` (REAL/DOUBLE), `T` (TEXT), etc.
- `----` - Separator between query and expected output

**Example test files:**
- `test/sql/read_alignments.test` - SAM/BAM reading functionality
- `test/sql/read_fastx.test` - FASTA/FASTQ reading functionality
- `test/sql/alignment_functions.test` - Alignment analysis functions
- `test/sql/copy_sam.test` - SAM/BAM writing functionality

**Good practices for SQL tests:**
- Test error cases first (missing files, invalid parameters)
- Test basic functionality with simple inputs
- Test edge cases (empty files, NULL values, boundary conditions)
- Verify data types and column ordering
- Keep test data files small and focused
- Use clear, descriptive comments for each test section

### C++ Unit Tests

Core algorithms and internal utilities are tested using the Catch2 framework. C++ tests are located in `test/cpp/` and compile to `./build/release/extension/miint/tests`.

**Running C++ tests:**
```sh
# Build and run all C++ tests
./build/release/extension/miint/tests

# Run tests matching a specific tag
./build/release/extension/miint/tests "[alignment_functions]"

# Run a specific test case
./build/release/extension/miint/tests "ParseCigar - Basic operations"
```

**Code organization for testing:**

1. **`miint` namespace** - Pure C++ code (algorithms, parsers, utilities)
   - Location: `src/` files, internal headers
   - Tested with: Catch2 C++ unit tests (`test/cpp/`)
   - Examples: CIGAR parsing, MD tag parsing, interval compression

2. **`duckdb` namespace** - DuckDB integration code (table functions, scalar functions)
   - Location: `src/` files (e.g., `alignment_functions.cpp`, `read_fastx.cpp`)
   - Tested with: SQL logic tests (`test/sql/`)
   - Examples: Function registration, DuckDB vector operations, parameter binding

**Example C++ test structure:**
```cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace miint {
    // Implementation code to test
    struct CigarStats { /* ... */ };
    static CigarStats ParseCigar(const std::string &cigar_str);
}

TEST_CASE("ParseCigar - Basic operations", "[alignment_functions]") {
    SECTION("Simple match") {
        auto stats = miint::ParseCigar("10M");
        REQUIRE(stats.matches == 10);
        REQUIRE(stats.alignment_columns == 10);
    }

    SECTION("With insertions and deletions") {
        auto stats = miint::ParseCigar("10M5I3D10M");
        REQUIRE(stats.matches == 20);
        REQUIRE(stats.insertions == 5);
        REQUIRE(stats.deletions == 3);
        REQUIRE(stats.gap_opens == 2);
    }
}
```

**Example test files:**
- `test/cpp/test_AlignmentFunctions.cpp` - CIGAR/MD parsing, identity calculations
- `test/cpp/test_SequenceReader.cpp` - FASTQ/FASTA parsing
- `test/cpp/test_SAMReader.cpp` - SAM/BAM reading logic
- `test/cpp/test_IntervalCompressor.cpp` - Interval merging algorithm

**Good practices for C++ tests:**
- Use `TEST_CASE()` for grouping related tests
- Use `SECTION()` for organizing variations within a test case
- Use descriptive test names that explain what is being tested
- Use `REQUIRE()` for critical assertions, `CHECK()` for non-critical
- Use matchers for floating-point comparisons (`REQUIRE_THAT(value, WithinRel(expected, 0.001))`)
- Test boundary conditions, error cases, and edge cases
- Keep test data minimal and focused on the specific behavior being tested

### Shell Script Tests

Some functionality cannot be tested through SQL logic tests, such as stdin handling. These features are tested using bash scripts located in `test/shell/`.

**Running shell tests:**
```sh
# Run all shell tests
for test in test/shell/*.sh; do bash "$test"; done

# Run a specific shell test
bash test/shell/read_alignments_stdin.sh

# All tests (SQL, C++, and shell)
bash run_tests.sh
```

**Test file structure:**
Shell tests use simple bash scripts with conditional logic and exit codes:

```bash
#!/bin/bash
set -e  # Exit on error

DUCKDB="./build/release/duckdb"
FAILED=0

# Helper function to run a test
run_test() {
    local test_name="$1"
    local expected="$2"
    shift 2
    local cmd="$@"

    echo "Running: $test_name"
    result=$(eval "$cmd" 2>&1 || true)

    if echo "$result" | grep -q "$expected"; then
        echo "  ✓ PASS"
    else
        echo "  ✗ FAIL"
        FAILED=1
    fi
}

# Test cases
run_test "Test description" \
    "expected output substring" \
    "cat data.txt | $DUCKDB -c 'SELECT * FROM read_alignments(\"/dev/stdin\");'"

# Return status
if [ $FAILED -eq 0 ]; then
    echo "All tests passed!"
    exit 0
else
    echo "Some tests failed!"
    exit 1
fi
```

**Example test files:**
- `test/shell/read_alignments_stdin.sh` - Tests for reading SAM/BAM from stdin

**Good practices for shell tests:**
- Use `set -e` to exit on first error
- Test both success and error cases
- Use `grep -q` to verify expected output substrings
- Always return proper exit codes (0 for success, 1 for failure)
- Use clear, descriptive test names
- Capture both stdout and stderr (`2>&1`)

### Test Data

Test data files are organized in `data/` subdirectories:
- `data/sam/` - SAM/BAM test files
- `data/fastq/` - FASTQ/FASTA test files
- `data/biom/` - BIOM format test files

**Important:** Keep test data files small (< 1KB when possible) to maintain fast test execution and minimize repository size.

---

This repository is based on https://github.com/duckdb/extension-template, check it out if you want to build and ship your own DuckDB extension.
