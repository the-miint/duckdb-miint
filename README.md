# MIINT: MIcrobiome INTelligence.

The MIINT [DuckDB](https://duckdb.org) extension enables DuckDB to interoperate with file formats and operations central to microbiome research.

## Quick Start

MIINT brings the power of SQL to microbiome data analysis. Here's a complete workflow analyzing metagenomic alignments:

```sql
-- Install and load the extension
INSTALL miint FROM community;
LOAD miint;

-- Read alignment files using glob pattern (files sorted alphabetically)
-- Use include_filepath to track source file as sample identifier
CREATE VIEW all_alignments AS
    SELECT *, regexp_extract(filepath, '.*/(.*)\.bam', 1) AS sample_id
    FROM read_alignments('samples/*.bam', include_filepath=true);

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
- **Read bioinformatics formats as tables**: FASTQ/FASTA, SAM/BAM, SFF, BIOM, Newick trees
- **Analyze with SQL**: Filter, aggregate, join with metadata tables
- **Write back to standard formats**: Export results to FASTQ, SAM/BAM, BIOM, Newick
- **Powerful functions**: Sequence identity, coverage, flag checking, Woltka classification, RYpe sequence classification
- **Performance**: Parallel processing, efficient compression, streaming I/O

## Table of Contents

- [Quick Start](#quick-start)
- [Installing](#installing)
  - [Installing from a local build](#installing-from-a-local-build)
- [Building](#building)
  - [Managing dependencies](#managing-dependencies)
  - [Build system](#build-system)
  - [Build steps](#build-steps)
- [Running the extension](#running-the-extension)
- [Functions](#functions)
  - [read_alignments](#read_alignmentsfilename-reference_lengthstable_name-include_filepathfalse-include_seq_qualfalse)
  - [read_fastx](#read_fastxfilename-sequence2filename-include_filepathfalse-qual_offset33)
  - [read_sequences_sff](#read_sequences_sfffilename-include_filepathfalse-trimtrue)
  - [read_biom](#read_biomfilename-include_filepathfalse)
  - [read_gff](#read_gffpath)
  - [read_ncbi](#read_ncbiaccession-api_key)
  - [read_ncbi_fasta](#read_ncbi_fastaaccession-api_key-include_filepathfalse)
  - [read_ncbi_annotation](#read_ncbi_annotationaccession-api_key-include_filepathfalse)
  - [read_jplace](#read_jplacepath)
  - [read_newick](#read_newickfilename-include_filepathfalse)
  - [align_minimap2](#align_minimap2query_table-subject_tablenull-index_pathnull-options)
  - [save_minimap2_index](#save_minimap2_indexsubject_table-output_path-options)
  - [align_minimap2_sharded](#align_minimap2_shardedquery_table-shard_directory-read_to_shard-options)
  - [align_bowtie2](#align_bowtie2query_table-subject_table-options)
  - [align_bowtie2_sharded](#align_bowtie2_shardedquery_table-shard_directory-read_to_shard-options)
  - [SAM Flag Functions](#sam-flag-functions)
  - [alignment_seq_identity](#alignment_seq_identitycigar-nm-md-type)
  - [alignment_query_length](#alignment_query_lengthcigar-include_hard_clipstrue)
  - [alignment_query_coverage](#alignment_query_coveragecigar-typealigned)
  - [RYpe Classification and Extraction Functions](#rype-classification-and-extraction-functions)
- [Analysis Functions](#analysis-functions)
  - [woltka_ogu_per_sample](#woltka_ogu_per_samplerelation-sample_id_field-sequence_id_field)
  - [woltka_ogu](#woltka_ogurelation-sequence_id_field)
  - [sequence_dna_reverse_complement / sequence_rna_reverse_complement](#sequence_dna_reverse_complementsequence-and-sequence_rna_reverse_complementsequence)
  - [sequence_dna_as_regexp / sequence_rna_as_regexp](#sequence_dna_as_regexpsequence-and-sequence_rna_as_regexpsequence)
  - [compress_intervals](#compress_intervalsstart-stop)
  - [genome_coverage](#genome_coveragealignments-subject_total_length-subject_genome_id)
  - [Pairwise Alignment Functions](#pairwise-alignment-functions)
- [Utility Functions](#utility-functions)
  - [miint_version](#miint_version)
- [COPY Formats](#copy-formats)
  - [FORMAT FASTQ](#copy--to--format-fastq)
  - [FORMAT FASTA](#copy--to--format-fasta)
  - [FORMAT SAM / FORMAT BAM](#copy--to--format-sam-and-copy--to--format-bam)
  - [FORMAT BIOM](#copy--to--format-biom)
  - [FORMAT NEWICK](#copy--to--format-newick)
- [Running the tests](#running-the-tests)
  - [SQL Logic Tests](#sql-logic-tests)
  - [C++ Unit Tests](#c-unit-tests)
  - [Shell Script Tests](#shell-script-tests)
  - [Test Data](#test-data)

## Installing

MIINT is available in the [DuckDB Community Extensions](https://community-extensions.duckdb.org/) repository:

```sql
INSTALL miint FROM community;
LOAD miint;
```

You can verify the installed version with:
```sql
SELECT miint_version();
```

### Installing from a local build

If you are building MIINT from source (see [Building](#building)), launch DuckDB with the `allow_unsigned_extensions` option:

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

Then load the extension binary directly:
```sql
LOAD '/path/to/miint.duckdb_extension';
```

## Building
### Managing dependencies

**Rust toolchain (required):** The RYpe sequence classification library is written in Rust. Install Rust via [rustup](https://rustup.rs/):
```shell
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

**VCPKG (required):** DuckDB extensions use VCPKG for dependency management. Enabling VCPKG is very simple: follow the [installation instructions](https://vcpkg.io/en/getting-started) or just run the following:
```shell
git clone https://github.com/Microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh
export VCPKG_TOOLCHAIN_PATH=`pwd`/vcpkg/scripts/buildsystems/vcpkg.cmake
```

### Build system
We use [Ninja](https://ninja-build.org) for quick builds. The easiest install is to use prebuilt [release](https://github.com/ninja-build/ninja/releases) binaries.

### Build steps
Now to build the extension, run:
```sh
bash build.sh
```

This is equivalent to:
```sh
export VCPKG_TOOLCHAIN_PATH=`pwd`/vcpkg/scripts/buildsystems/vcpkg.cmake
GEN=ninja make
```

The main binaries that will be built are:

```sh
./build/release/duckdb
./build/release/test/unittest
./build/release/extension/miint/miint.duckdb_extension
./build/release/extension/miint/tests
```

- `duckdb` is the binary for the DuckDB shell with the extension code automatically loaded.
- `unittest` is the DuckDB test runner with the extension linked in (for SQL logic tests).
- `miint.duckdb_extension` is the loadable binary as it would be distributed.
- `tests` is the Catch2 C++ unit test runner.

## Running the extension
To run the extension code, simply start the shell with `./build/release/duckdb`.

## Functions

### `read_alignments(filename, [reference_lengths='table_name'], [include_filepath=false], [include_seq_qual=false])`
Read SAM/BAM alignment files.

**Note:** `read_sam` is still supported as a backward-compatible alias.

**Parameters:**
- `filename` (VARCHAR or VARCHAR[]): Path to SAM/BAM file(s), glob pattern (e.g., `'data/*.bam'`), or `-` / `/dev/stdin` for standard input
  - **Glob patterns**: When a single VARCHAR contains glob characters (`*`, `?`, `[`), files are expanded and sorted alphabetically
  - **Arrays**: VARCHAR[] elements are treated as literal paths (no glob expansion)
- `reference_lengths` (VARCHAR, optional): Table or view name containing reference sequences for headerless SAM files. Must have at least 2 columns: first column = reference name (VARCHAR), second column = reference length (INTEGER/BIGINT). Column names don't matter. Views are fully supported and can include computed columns.
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

-- Read all BAM files matching a glob pattern (sorted alphabetically)
SELECT * FROM read_alignments('samples/*.bam', include_filepath=true);

-- Glob pattern matching specific files
SELECT COUNT(*) FROM read_alignments('data/sample_*.bam');

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
- `filename` (VARCHAR or VARCHAR[]): Path to FASTA/FASTQ file(s), glob pattern (e.g., `'data/*.fastq'`), or R1 files for paired-end
  - **Glob patterns**: When a single VARCHAR contains glob characters (`*`, `?`, `[`), files are expanded and sorted alphabetically
  - **Arrays**: VARCHAR[] elements are treated as literal paths (no glob expansion)
- `sequence2` (VARCHAR or VARCHAR[], optional): Path to R2 file(s) for paired-end reads. Must have same number of files as `filename`
  - **Paired-end with globs**: When `filename` is a glob pattern, `sequence2` must also be a glob pattern. Both are expanded and sorted independently, then paired by position. The expanded file counts must match.
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

-- Read all FASTQ files matching a glob pattern (sorted alphabetically)
SELECT * FROM read_fastx('samples/*.fastq', include_filepath=true);

-- Paired-end with glob patterns (both must be globs, paired by sorted order)
SELECT * FROM read_fastx(
    'samples/*_R1.fastq',
    sequence2='samples/*_R2.fastq'
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

### `read_sequences_sff(filename, [include_filepath=false], [trim=true])`
Read SFF (Standard Flowgram Format) files from 454/Roche sequencing platforms.

**Parameters:**
- `filename` (VARCHAR or VARCHAR[]): Path to SFF file(s) or glob pattern (e.g., `'data/*.sff'`)
  - **Glob patterns**: When a single VARCHAR contains glob characters (`*`, `?`, `[`), files are expanded and sorted alphabetically
  - **Arrays**: VARCHAR[] elements are treated as literal paths (no glob expansion)
  - **Stdin not supported**: SFF is a binary format that requires seeking, so `-` and `/dev/stdin` are not allowed
- `include_filepath` (BOOLEAN, optional, default false): Add filepath column to output
- `trim` (BOOLEAN, optional, default true): Apply quality and adapter clip positions from the SFF header to trim sequences and quality scores

**Output schema:**
- `sequence_index` (BIGINT): 1-based sequential index per file (resets to 1 for each file)
- `read_id` (VARCHAR): Read name from the SFF read header
- `comment` (VARCHAR, nullable): Always NULL (SFF format has no comment field)
- `sequence1` (VARCHAR): DNA sequence (trimmed by default)
- `sequence2` (VARCHAR, nullable): Always NULL (SFF is single-end only)
- `qual1` (UINT8[], nullable): Quality scores as array of integers (trimmed by default)
- `qual2` (UINT8[], nullable): Always NULL (SFF is single-end only)
- `filepath` (VARCHAR, optional): File path when include_filepath=true

**Behavior:**
- Schema matches `read_fastx` for `UNION ALL` compatibility
- When `trim=true` (default), quality and adapter clip positions from the SFF header are applied. Overlapping clips produce an empty sequence
- When `trim=false`, the full untrimmed sequence and quality scores are returned
- SFF index blocks (if present) are automatically skipped
- Supports parallel processing (up to 8 threads, one per file)

**Examples:**
```sql
-- Read a single SFF file
SELECT * FROM read_sequences_sff('reads.sff');

-- Read without trimming
SELECT * FROM read_sequences_sff('reads.sff', trim=false);

-- Read multiple SFF files with glob
SELECT * FROM read_sequences_sff('data/*.sff', include_filepath=true);

-- Read multiple SFF files via array
SELECT * FROM read_sequences_sff(['batch1.sff', 'batch2.sff']);

-- Combine SFF and FASTQ data (schemas match)
SELECT sequence_index, read_id, comment, sequence1, sequence2, qual1, qual2
FROM read_sequences_sff('legacy_454.sff')
UNION ALL
SELECT * FROM read_fastx('modern_illumina.fastq');

-- Filter by average quality
SELECT read_id, list_avg(qual1)::INTEGER as avg_qual
FROM read_sequences_sff('reads.sff')
WHERE list_avg(qual1) >= 30;
```

**Notes:**
- The `comment`, `sequence2`, and `qual2` columns are always NULL for SFF files
- The `sequence_index` resets to 1 for each file. Use `include_filepath=true` to distinguish sequences from different files
- Stdin is not supported because SFF is a binary format that requires file seeking

### `read_biom(filename, [include_filepath=false])`
Read BIOM (Biological Observation Matrix) format files.

**Parameters:**
- `filename` (VARCHAR or VARCHAR[]): Path to BIOM file(s) or glob pattern (e.g., `'data/*.biom'`)
  - **Glob patterns**: When a single VARCHAR contains glob characters (`*`, `?`, `[`), files are expanded and sorted alphabetically
  - **Arrays**: VARCHAR[] elements are treated as literal paths (no glob expansion)
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

-- Read all BIOM files matching a glob pattern (sorted alphabetically)
SELECT * FROM read_biom('results/*.biom', include_filepath=true);

-- Glob pattern for batch processing
SELECT sample_id, SUM(value) as total_count
FROM read_biom('batches/batch_*.biom')
GROUP BY sample_id;

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

### `read_ncbi(accession, [api_key])`

Fetch GenBank metadata from NCBI by accession number. This function queries NCBI's E-utilities API to retrieve sequence metadata without downloading the full sequence.

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to NCBI servers (eutils.ncbi.nlm.nih.gov)

**Parameters:**
- `accession` (VARCHAR or VARCHAR[]): NCBI accession number(s) (e.g., 'NC_001416.1', 'NM_001101.5')
- `api_key` (VARCHAR, optional): NCBI API key for higher rate limits (10 req/s vs 3 req/s)

**Output schema:**
- `accession` (VARCHAR): Accession with version (e.g., 'NC_001416.1')
- `version` (INTEGER): Version number extracted from accession
- `description` (VARCHAR): Sequence description/definition
- `organism` (VARCHAR): Source organism name
- `taxonomy_id` (BIGINT): NCBI taxonomy ID
- `length` (BIGINT): Sequence length in base pairs
- `molecule_type` (VARCHAR): Molecule type (e.g., 'DNA', 'RNA', 'protein')
- `update_date` (DATE): Last modification date

**Behavior:**
- Fetches GenBank XML format from E-utilities
- Automatically detects accession type (RefSeq NC_/NM_/NP_, GenBank, etc.)
- Rate-limited to respect NCBI guidelines (3 req/s without key, 10 req/s with key)
- Retry logic for transient failures (429, 500, 502, 503) with exponential backoff

**Examples:**
```sql
-- Get metadata for a single accession
SELECT * FROM read_ncbi('NC_001416.1');

-- Get metadata for multiple accessions
SELECT * FROM read_ncbi(['NC_001416.1', 'NC_001422.1']);

-- Use API key for higher rate limits
SELECT * FROM read_ncbi('NC_001416.1', api_key='your_ncbi_api_key');

-- Get organism and length for a set of sequences
SELECT accession, organism, length
FROM read_ncbi(['NC_001416.1', 'NC_001422.1', 'NC_000913.3']);

-- Filter by molecule type
SELECT accession, description
FROM read_ncbi(['NC_001416.1', 'NM_001101.5'])
WHERE molecule_type = 'DNA';
```

**Notes:**
- Empty accessions will raise an error
- Invalid accessions will raise an error with the specific accession that failed
- For bulk downloads, consider using the API key to avoid rate limiting

### `read_ncbi_fasta(accession, [api_key], [include_filepath=false])`

Fetch FASTA sequences from NCBI by accession number. Returns data in the same schema as `read_fastx`, making it easy to combine NCBI sequences with local files.

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to NCBI servers (eutils.ncbi.nlm.nih.gov)

**Parameters:**
- `accession` (VARCHAR or VARCHAR[]): NCBI accession number(s) (e.g., 'NC_001416.1')
- `api_key` (VARCHAR, optional): NCBI API key for higher rate limits
- `include_filepath` (BOOLEAN, optional, default false): Add filepath column showing NCBI URL

**Output schema (matches `read_fastx`):**
- `sequence_index` (BIGINT): 0-based sequential index
- `read_id` (VARCHAR): Accession number extracted from FASTA header
- `comment` (VARCHAR, nullable): Description from FASTA header
- `sequence1` (VARCHAR): DNA/RNA sequence
- `sequence2` (VARCHAR): Always NULL (single sequences only)
- `qual1` (UINT8[]): Always NULL (FASTA has no quality scores)
- `qual2` (UINT8[]): Always NULL
- `filepath` (VARCHAR, optional): NCBI E-utilities URL when include_filepath=true

**Behavior:**
- Fetches FASTA format from E-utilities
- Parses pipe-delimited FASTA headers (e.g., `gi|123|ref|NC_001416.1|description`)
- Extracts accession as `read_id`, remainder as `comment`
- Rate-limited with retry logic for robustness

**Examples:**
```sql
-- Fetch a single sequence
SELECT * FROM read_ncbi_fasta('NC_001416.1');

-- Fetch multiple sequences
SELECT read_id, LENGTH(sequence1) AS seq_length
FROM read_ncbi_fasta(['NC_001416.1', 'NC_001422.1']);

-- Combine NCBI sequences with local files
SELECT read_id, sequence1 FROM read_ncbi_fasta('NC_001416.1')
UNION ALL
SELECT read_id, sequence1 FROM read_fastx('local_sequences.fasta');

-- Export NCBI sequence to local FASTA file
COPY (SELECT * FROM read_ncbi_fasta('NC_001416.1'))
TO 'lambda_phage.fasta' (FORMAT FASTA);

-- Get sequence with source URL tracking
SELECT read_id, LENGTH(sequence1) AS length, filepath
FROM read_ncbi_fasta('NC_001416.1', include_filepath=true);

-- Use as reference for alignment
CREATE TABLE reference AS SELECT * FROM read_ncbi_fasta('NC_001416.1');
CREATE TABLE reads AS SELECT * FROM read_fastx('reads.fastq');
SELECT * FROM align_minimap2('reads', 'reference');
```

**Notes:**
- Large sequences (e.g., complete chromosomes) may take time to download
- Output is compatible with all functions expecting `read_fastx` schema
- Empty accessions will raise an error

### `read_ncbi_annotation(accession, [api_key], [include_filepath=false])`

Fetch feature annotations from NCBI by accession number. Returns data in the same schema as `read_gff`, making it easy to analyze NCBI annotations with SQL.

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to NCBI servers (eutils.ncbi.nlm.nih.gov)

**Parameters:**
- `accession` (VARCHAR or VARCHAR[]): NCBI accession number(s) (e.g., 'NC_001416.1')
- `api_key` (VARCHAR, optional): NCBI API key for higher rate limits
- `include_filepath` (BOOLEAN, optional, default false): Add filepath column showing NCBI URL

**Output schema (matches `read_gff`):**
- `seqid` (VARCHAR): Sequence/chromosome identifier (the accession)
- `source` (VARCHAR): Annotation source ('RefSeq', 'GenBank', or 'NCBI')
- `type` (VARCHAR): Feature type (e.g., 'gene', 'CDS', 'mRNA')
- `position` (INTEGER): 1-based start position
- `stop_position` (INTEGER): 1-based end position
- `score` (DOUBLE, nullable): Feature score (typically NULL for NCBI data)
- `strand` (VARCHAR): Strand ('+' or '-')
- `phase` (INTEGER, nullable): CDS phase (0, 1, 2) derived from codon_start qualifier
- `attributes` (MAP(VARCHAR, VARCHAR)): Key-value attributes (gene, product, locus_tag, etc.)
- `filepath` (VARCHAR, optional): NCBI E-utilities URL when include_filepath=true

**Behavior:**
- Fetches INSDC feature table format from E-utilities
- Converts to GFF3-compatible schema
- Detects source automatically from accession prefix (RefSeq for NC_/NM_, etc.)
- Handles complement strand (reversed positions)
- Parses codon_start qualifier to set correct CDS phase
- Warns on complex locations (join, complement) - uses outer bounds

**Examples:**
```sql
-- Get all annotations for Lambda phage
SELECT * FROM read_ncbi_annotation('NC_001416.1');

-- Get only genes
SELECT seqid, position, stop_position, attributes['gene'] AS gene_name
FROM read_ncbi_annotation('NC_001416.1')
WHERE type = 'gene';

-- Get CDS features with their products
SELECT position, stop_position, strand, phase,
       attributes['gene'] AS gene,
       attributes['product'] AS product
FROM read_ncbi_annotation('NC_001416.1')
WHERE type = 'CDS';

-- Count feature types
SELECT type, COUNT(*) AS count
FROM read_ncbi_annotation('NC_001416.1')
GROUP BY type
ORDER BY count DESC;

-- Combine annotations from multiple genomes
SELECT seqid, type, COUNT(*) AS feature_count
FROM read_ncbi_annotation(['NC_001416.1', 'NC_001422.1'])
GROUP BY seqid, type
ORDER BY seqid, feature_count DESC;

-- Find genes on the minus strand
SELECT attributes['gene'] AS gene_name, position, stop_position
FROM read_ncbi_annotation('NC_001416.1')
WHERE type = 'gene' AND strand = '-';

-- Export to local GFF-like format
COPY (
    SELECT seqid, source, type, position, stop_position,
           NULL AS score, strand, phase
    FROM read_ncbi_annotation('NC_001416.1')
) TO 'lambda_features.tsv' (DELIMITER '\t', HEADER);
```

**Notes:**
- Output is compatible with queries designed for `read_gff`
- Complex feature locations (join, complement) emit a warning and use outer bounds only
- The `phase` column is set from the `codon_start` qualifier for CDS features
- Source is detected from accession prefix: NC_/NM_/NP_ → 'RefSeq', others → 'GenBank' or 'NCBI'

### `read_jplace(path)`

Read jplace phylogenetic placement files. The jplace format stores query sequence placements onto a reference phylogenetic tree, as defined in [Matsen et al. 2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009).

**Parameters:**
- `path` (VARCHAR): Path to jplace file(s), supports glob patterns (e.g., `'data/*.jplace'`)

**Output schema:**
- `fragment` (VARCHAR): Fragment/sequence name (from `nm` or `n` field)
- `edge_num` (INTEGER): Edge number in the reference tree where placement occurs
- `likelihood` (DOUBLE): Log likelihood of the placement
- `like_weight_ratio` (DOUBLE): Likelihood weight ratio (proportion of total likelihood)
- `distal_length` (DOUBLE): Distance from distal end of edge to placement point
- `pendant_length` (DOUBLE): Pendant branch length (branch to the placed sequence)
- `filepath` (VARCHAR): Source file path

**Behavior:**
- Returns only the best placement (first in `p` array) for each fragment
- Supports both `nm` (named multiplicities) and `n` (names) formats
- Supports glob patterns for reading multiple files
- Requires the json extension (automatically loaded)

**Examples:**
```sql
-- Read a single jplace file
SELECT * FROM read_jplace('placements.jplace');

-- Read multiple jplace files with glob pattern
SELECT * FROM read_jplace('results/*.jplace');

-- Filter placements by likelihood weight ratio
SELECT fragment, edge_num, like_weight_ratio
FROM read_jplace('placements.jplace')
WHERE like_weight_ratio > 0.5
ORDER BY like_weight_ratio DESC;

-- Count placements per edge
SELECT edge_num, COUNT(*) AS num_placements
FROM read_jplace('placements.jplace')
GROUP BY edge_num
ORDER BY num_placements DESC;

-- Aggregate placements from multiple files
SELECT filepath, COUNT(*) AS num_fragments
FROM read_jplace('batch_*.jplace')
GROUP BY filepath;

-- Find high-confidence placements
SELECT fragment, edge_num, likelihood, like_weight_ratio
FROM read_jplace('placements.jplace')
WHERE like_weight_ratio >= 0.9;

-- Join with edge metadata (if available)
CREATE TABLE edge_taxa AS
SELECT 0 AS edge_num, 'Bacteria' AS taxon
UNION ALL SELECT 1, 'Archaea'
UNION ALL SELECT 2, 'Eukarya';

SELECT p.fragment, e.taxon, p.like_weight_ratio
FROM read_jplace('placements.jplace') p
JOIN edge_taxa e ON p.edge_num = e.edge_num;
```

**jplace Format Notes:**
- Standard format for phylogenetic placement tools (pplacer, EPA-ng, etc.)
- The `fields` array in the file defines column order: typically `[edge_num, likelihood, like_weight_ratio, distal_length, pendant_length]`
- Multiple placements per fragment are supported in the format, but this function returns only the best (first) placement
- The `nm` field contains `[[name, multiplicity], ...]` pairs; `n` field contains simple name arrays

**Use Cases:**
- **Taxonomic profiling**: Analyze where metagenomic reads place on a reference tree
- **Community composition**: Aggregate placements to quantify taxa
- **Quality filtering**: Filter by likelihood weight ratio to retain confident placements
- **Multi-sample analysis**: Process multiple jplace files with glob patterns

**Compatibility:**
- Reads jplace files created by:
  - pplacer
  - EPA-ng
  - SEPP
  - Other tools following the jplace specification

**Implementation note:** Implemented as a DuckDB macro using `read_json` with JSON path extraction.

### `read_newick(filename, [include_filepath=false])`

Read Newick phylogenetic tree files and return a table with one row per node.

**Parameters:**
- `filename` (VARCHAR or VARCHAR[]): Path to Newick file(s), glob pattern (e.g., `'data/*.nwk'`), or `-` / `/dev/stdin` for standard input
  - **Glob patterns**: When a single VARCHAR contains glob characters (`*`, `?`, `[`), files are expanded and sorted alphabetically
  - **Arrays**: VARCHAR[] elements are treated as literal paths (no glob expansion)
- `include_filepath` (BOOLEAN, optional, default false): Add filepath column to output

**Output schema:**
- `node_index` (BIGINT): 0-based index of node in tree (internal representation)
- `name` (VARCHAR): Node label (empty string for unlabeled internal nodes)
- `branch_length` (DOUBLE, nullable): Branch length (NULL if not specified in file)
- `edge_id` (BIGINT, nullable): Edge identifier from jplace format `{n}` syntax (NULL if not specified)
- `parent_index` (BIGINT, nullable): Parent node's node_index (NULL for root)
- `is_tip` (BOOLEAN): Whether node is a tip/leaf (has no children)
- `filepath` (VARCHAR, optional): File path when include_filepath=true

**Behavior:**
- Parses standard Newick format including:
  - Node names (quoted or unquoted)
  - Branch lengths (`:0.123`)
  - Edge identifiers for jplace compatibility (`{0}`, `{1}`, etc.)
- Supports gzip-compressed files (auto-detected from `.gz` extension)
- Supports stdin input using `-` or `/dev/stdin` (single file only)
- Returns exactly one row per node in the tree
- Root node has `parent_index = NULL`

**Examples:**
```sql
-- Read a single Newick file
SELECT * FROM read_newick('tree.nwk');

-- Get tip names only
SELECT name FROM read_newick('tree.nwk')
WHERE is_tip = true;

-- Count nodes in tree
SELECT COUNT(*) AS total_nodes,
       COUNT(*) FILTER (WHERE is_tip) AS tips,
       COUNT(*) FILTER (WHERE NOT is_tip) AS internal_nodes
FROM read_newick('tree.nwk');

-- Read tree with edge IDs (jplace format)
SELECT name, edge_id, branch_length
FROM read_newick('reference.nwk')
WHERE edge_id IS NOT NULL;

-- Read multiple trees with glob pattern
SELECT filepath, COUNT(*) AS num_nodes
FROM read_newick('trees/*.nwk', include_filepath=true)
GROUP BY filepath;

-- Read gzip-compressed tree
SELECT * FROM read_newick('tree.nwk.gz');

-- Read from stdin
SELECT * FROM read_newick('/dev/stdin');

-- Find the root node
SELECT * FROM read_newick('tree.nwk')
WHERE parent_index IS NULL;

-- Calculate tree depth (distance from root)
WITH RECURSIVE node_depth AS (
    SELECT node_index, name, parent_index, 0 AS depth
    FROM read_newick('tree.nwk')
    WHERE parent_index IS NULL
    UNION ALL
    SELECT t.node_index, t.name, t.parent_index, nd.depth + 1
    FROM read_newick('tree.nwk') t
    JOIN node_depth nd ON t.parent_index = nd.node_index
)
SELECT name, depth FROM node_depth
WHERE is_tip ORDER BY depth DESC;
```

**Stdin limitations:**
- Cannot be used in multi-file arrays (e.g., `['/dev/stdin', 'file.nwk']` will error)
- User must know the format of their stdin data

**Roundtrip with COPY FORMAT NEWICK:**
Trees can be read with `read_newick`, modified via SQL, and written back to Newick format:
```sql
-- Read, filter to subtree, write back
COPY (
    SELECT node_index, name, branch_length, edge_id, parent_index
    FROM read_newick('input.nwk')
    -- Your filtering/modification logic here
) TO 'output.nwk' (FORMAT NEWICK);
```

**Newick Format Notes:**
- Standard format for representing phylogenetic trees as text
- Parentheses denote tree structure: `(A,B)` means A and B share a parent
- Colons precede branch lengths: `A:0.1` means tip A with branch length 0.1
- Semicolons terminate the tree: `(A,B);`
- Edge IDs in braces are an extension for jplace: `A:0.1{0}` has edge_id 0

### `align_minimap2(query_table, [subject_table=NULL], [index_path=NULL], [options])`

Align query sequences to subject sequences using minimap2. This function enables sequence alignment directly within SQL by reading sequences from DuckDB tables/views and returning alignments in the same format as `read_alignments`.

**Performance:** For large reference databases (e.g., human genome), use pre-built indexes via `index_path` for 10-30x faster alignment. Build indexes once with `save_minimap2_index`, then reuse them across multiple queries.

**Parameters:**
- `query_table` (VARCHAR): Name of table or view containing query sequences. Must have `read_fastx`-compatible schema (read_id, sequence1, optional sequence2/qual1/qual2)
- `subject_table` (VARCHAR, optional): Name of table or view containing subject/reference sequences. Must have `read_fastx`-compatible schema. Cannot contain paired-end data (sequence2 must be NULL or absent). **Either `subject_table` OR `index_path` must be provided, but not both.**
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

**Output schema:**
Returns the same schema as `read_alignments` (21 columns):
- `read_id` (VARCHAR): Query sequence identifier
- `flags` (USMALLINT): SAM alignment flags
- `reference` (VARCHAR): Subject sequence identifier
- `position` (BIGINT): 1-based start position on reference
- `stop_position` (BIGINT): 1-based stop position on reference
- `mapq` (UTINYINT): Mapping quality
- `cigar` (VARCHAR): CIGAR string (with =/X by default)
- `mate_reference` (VARCHAR): Mate reference (for paired-end)
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
SELECT reference,
       compress_intervals(position, stop_position) AS coverage_regions,
       SUM(stop_position - position + 1) AS total_aligned_bases
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
- Secondary alignments can significantly increase output size; use `max_secondary=0` for primary-only results

**Limitations:**
- Subject sequences must fit in memory (loaded at bind time for indexing when using `subject_table`)
- No support for reading sequences directly from files (use tables/views from `read_fastx`)

### `save_minimap2_index(subject_table, output_path, [options])`

Build and save a minimap2 index to disk for reuse. This provides 10-30x performance improvement when aligning multiple query sets against the same large reference database.

**Use case:** Build indexes once for large reference databases (e.g., WoLr2 phylogenetic markers, RefSeq genomes, custom OGU databases), then use them repeatedly with `align_minimap2(..., index_path='file.mmi')` instead of rebuilding the index each time.

**Parameters:**
- `subject_table` (VARCHAR): Name of table or view containing subject/reference sequences. Must have `read_fastx`-compatible schema. Cannot contain paired-end data (sequence2 must be NULL or absent)
- `output_path` (VARCHAR): Path where the index file (.mmi) will be saved
- `preset` (VARCHAR, default: 'sr'): Minimap2 preset (same options as `align_minimap2`)
- `k` (INTEGER, optional): K-mer size (overrides preset default if specified)
- `w` (INTEGER, optional): Minimizer window size (overrides preset default if specified)
- `eqx` (BOOLEAN, default: true): Use =/X CIGAR operators instead of M

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

### `align_minimap2_sharded(query_table, shard_directory, read_to_shard, [options])`

Align query sequences against multiple pre-built minimap2 index shards in parallel. Each shard is a separate `.mmi` index file, and a mapping table specifies which reads should be aligned against which shard. This is designed for large-scale metagenomic workflows where the reference database is split across multiple shards and reads have been pre-assigned to shards (e.g., by a prior classification step).

**Parameters:**
- `query_table` (VARCHAR): Name of table or view containing query sequences. Must have `read_fastx`-compatible schema (read_id, sequence1, optional sequence2/qual1/qual2)
- `shard_directory` (VARCHAR, required): Path to directory containing pre-built minimap2 index files. Each shard's index is expected at `<shard_directory>/<shard_name>.mmi`
- `read_to_shard` (VARCHAR, required): Name of table or view that maps reads to shards. Must have columns:
  - `read_id` (VARCHAR): Read identifier (must match read_id in query_table)
  - `shard_name` (VARCHAR): Name of the shard this read should be aligned against
- `preset` (VARCHAR, default: 'sr'): Minimap2 preset ('sr', 'map-ont', 'map-pb', etc.)
- `max_secondary` (INTEGER, default: 5): Maximum secondary alignments per query. Set to 0 for primary only
- `eqx` (BOOLEAN, default: true): Use =/X CIGAR operators instead of M

**Output schema:**
Returns the same 21-column schema as `align_minimap2` and `read_alignments`.

**Behavior:**
- At bind time, reads the `read_to_shard` table to discover shards and validate that each `<shard_name>.mmi` file exists in `shard_directory`
- Shards are processed in parallel (one DuckDB thread per shard), each loading its `.mmi` index independently
- For each shard, only the reads assigned to that shard (via the `read_to_shard` mapping) are queried
- A read can appear in multiple shards (mapped to multiple shard_name values) and will be aligned against each
- Unmapped reads (flag 0x4) are automatically filtered out of results
- Supports both single-end and paired-end query sequences
- Supports views for both `query_table` and `read_to_shard`

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

### `align_bowtie2(query_table, subject_table, [options])`

Align query sequences to subject sequences using Bowtie2. This function enables short-read alignment directly within SQL by reading sequences from DuckDB tables/views and returning alignments in the same format as `read_alignments`.

**Requirements:**
- Bowtie2 must be installed and available in PATH (`bowtie2` and `bowtie2-build` commands)

**Parameters:**
- `query_table` (VARCHAR): Name of table or view containing query sequences. Must have `read_fastx`-compatible schema (read_id, sequence1, optional sequence2/qual1/qual2)
- `subject_table` (VARCHAR): Name of table or view containing subject/reference sequences. Must have `read_fastx`-compatible schema. Cannot contain paired-end data (sequence2 must be NULL or absent)
- `preset` (VARCHAR, optional): Bowtie2 preset for alignment sensitivity
  - `'very-fast'`: Fastest, least sensitive
  - `'fast'`: Fast alignment
  - `'sensitive'`: More sensitive (slower)
  - `'very-sensitive'`: Most sensitive, slowest
- `local` (BOOLEAN, default: false): Use local alignment mode instead of end-to-end
  - `false` (default): End-to-end alignment (entire read must align)
  - `true`: Local alignment (soft-clipping allowed at ends)
- `threads` (INTEGER, default: 1): Number of threads for Bowtie2 alignment (-p parameter)
- `max_secondary` (INTEGER, default: 1): Maximum alignments to report per query (-k parameter)
- `extra_args` (VARCHAR, optional): Additional Bowtie2 command-line arguments (space-separated)
- `quiet` (BOOLEAN, default: true): Suppress Bowtie2 stderr output (alignment statistics)

**Output schema:**
Returns the same schema as `read_alignments` (21 columns):
- `read_id` (VARCHAR): Query sequence identifier
- `flags` (USMALLINT): SAM alignment flags
- `reference` (VARCHAR): Subject sequence identifier
- `position` (BIGINT): 1-based start position on reference
- `stop_position` (BIGINT): 1-based stop position on reference
- `mapq` (UTINYINT): Mapping quality
- `cigar` (VARCHAR): CIGAR string
- `mate_reference` (VARCHAR): Mate reference (for paired-end)
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

**Behavior:**
- Subject sequences are loaded into memory at bind time and indexed (must fit in RAM)
- Query sequences are processed in batches for memory efficiency
- Supports both single-end and paired-end query sequences (paired-end uses interleaved format internally)
- Uses Bowtie2's subprocess interface for alignment
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

-- Pass additional Bowtie2 arguments
SELECT * FROM align_bowtie2('queries', 'subjects', extra_args='--no-unal --rdg 5,3');
```

**Error handling:**
- Error if query_table or subject_table does not exist
- Error if subject_table contains paired-end data (sequence2 not NULL)
- Error if tables lack required columns (read_id, sequence1)
- Error if bowtie2 or bowtie2-build is not found in PATH

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

### `align_bowtie2_sharded(query_table, shard_directory, read_to_shard, [options])`

Align query sequences against multiple pre-built Bowtie2 index shards in parallel. Each shard is a directory containing a Bowtie2 index (with prefix `index`), and a mapping table specifies which reads should be aligned against which shard. This is the Bowtie2 counterpart to `align_minimap2_sharded`, designed for the same large-scale sharded alignment workflows.

**Requirements:**
- Bowtie2 must be installed and available in PATH (`bowtie2` command)

**Parameters:**
- `query_table` (VARCHAR): Name of table or view containing query sequences. Must have `read_fastx`-compatible schema (read_id, sequence1, optional sequence2/qual1/qual2)
- `shard_directory` (VARCHAR, required): Path to directory containing shard subdirectories. Each shard's Bowtie2 index is expected at `<shard_directory>/<shard_name>/index` (i.e., files like `<shard_name>/index.1.bt2`, `<shard_name>/index.rev.1.bt2`, etc.)
- `read_to_shard` (VARCHAR, required): Name of table or view that maps reads to shards. Must have columns:
  - `read_id` (VARCHAR): Read identifier (must match read_id in query_table)
  - `shard_name` (VARCHAR): Name of the shard this read should be aligned against
- `preset` (VARCHAR, optional): Bowtie2 sensitivity preset ('very-fast', 'fast', 'sensitive', 'very-sensitive')
- `local` (BOOLEAN, default: false): Use local alignment mode instead of end-to-end
- `max_secondary` (INTEGER, default: 1): Maximum alignments to report per query (-k parameter)
- `extra_args` (VARCHAR, optional): Additional Bowtie2 command-line arguments (space-separated)
- `quiet` (BOOLEAN, default: true): Suppress Bowtie2 stderr output
- `threads` (INTEGER): Ignored in sharded mode. Parallelism comes from running multiple single-threaded Bowtie2 processes (one per shard). A warning is printed if set to a value other than 1

**Output schema:**
Returns the same 21-column schema as `align_bowtie2` and `read_alignments`.

**Behavior:**
- At bind time, reads the `read_to_shard` table to discover shards and validate that each Bowtie2 index exists at `<shard_directory>/<shard_name>/index`
- Shards are processed in parallel, with each DuckDB thread running a separate single-threaded Bowtie2 process
- For each shard, only the reads assigned to that shard (via the `read_to_shard` mapping) are queried
- A read can appear in multiple shards and will be aligned against each
- Unmapped reads (flag 0x4) are automatically filtered out of results
- Supports both single-end and paired-end query sequences
- Supports views for both `query_table` and `read_to_shard`

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

-- Pass additional Bowtie2 arguments
SELECT * FROM align_bowtie2_sharded('queries',
    shard_directory := 'indexes/',
    read_to_shard := 'read_to_shard',
    extra_args := '--no-unal --rdg 5,3');
```

**Error handling:**
- Error if `shard_directory` does not exist
- Error if any shard's Bowtie2 index files are missing (e.g., no `index.1.bt2` at `<shard_directory>/<shard_name>/index`)
- Error if `query_table` or `read_to_shard` table/view does not exist
- Error if `read_to_shard` table is missing `read_id` or `shard_name` columns
- Error if `read_to_shard` table contains NULL `shard_name` values
- Error if `bowtie2` is not found in PATH

**Performance notes:**
- Each shard runs a single-threaded Bowtie2 process; parallelism comes from running multiple shards concurrently
- Control shard parallelism with `SET threads=N` in DuckDB
- The `threads` parameter is ignored (always 1 per shard) to avoid CPU oversubscription
- Shards are sorted by read count (largest first) for better load balancing

**Comparison of sharded vs non-sharded alignment functions:**
| Feature | `align_minimap2` / `align_bowtie2` | `align_minimap2_sharded` / `align_bowtie2_sharded` |
|---------|--------------------------------------|------------------------------------------------------|
| Index source | Build on-the-fly or single pre-built index | Multiple pre-built indexes (one per shard) |
| Read routing | All reads against one index | Reads routed to specific shards via mapping table |
| Parallelism | Single aligner thread(s) | One aligner per shard, concurrent |
| Use case | Single reference database | Sharded reference databases (e.g., from prior classification) |

### SAM Flag Functions

Test individual SAM flag bits. Each function takes a `USMALLINT` (the flags column from `read_alignments`) and returns a `BOOLEAN`.

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

### RYpe Classification and Extraction Functions

[RYpe](https://github.com/biocore/rype) is a minimizer-based sequence classification tool that uses RY-space encoding (purine/pyrimidine) for robust sequence matching. These functions require a RYpe index directory (`.ryxdi`), which contains a Parquet-based inverted index built from reference sequences.

All RYpe functions take a `sequence_table` parameter that references a DuckDB table or view containing sequence data. The table must have a `sequence1` column and an identifier column (default `read_id`).

#### `rype_classify(index_path, sequence_table, [id_column='read_id'], [threshold=0.1], [negative_index=path])`

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

#### `rype_extract_minimizer_set(sequence_table, k, w, [salt=6148914691236517205], [id_column='read_id'])`

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

#### `rype_extract_strand_minimizers(sequence_table, k, w, [salt=6148914691236517205], [id_column='read_id'])`

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

### `genome_coverage(alignments, subject_total_length, subject_genome_id)`

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

### Pairwise Alignment Functions

Gap-affine pairwise sequence alignment powered by [WFA2-lib](https://github.com/smarco/WFA2-lib) (Wavefront Alignment Algorithm). Three functions at increasing detail levels:

- `align_pairwise_score` — alignment score only (fastest)
- `align_pairwise_cigar` — score + extended CIGAR string
- `align_pairwise_full` — score + CIGAR + aligned sequences with gap characters

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
- `align_pairwise_score` → `INTEGER` — alignment score (0 = identical, higher = more divergent)
- `align_pairwise_cigar` → `STRUCT(score INTEGER, cigar VARCHAR)` — score and extended CIGAR (`=`/`X` ops, not `M`)
- `align_pairwise_full` → `STRUCT(score INTEGER, cigar VARCHAR, query_aligned VARCHAR, subject_aligned VARCHAR)` — score, CIGAR, and aligned sequences with `-` gap characters

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

### `miint_version()`

Returns the MIINT extension version string.

**Returns:** VARCHAR - Version string derived from the git tag at build time (e.g., `v0.1.0`, or `v0.1.0-3-gabcdef1` if built from a commit after a tag).

**Example:**
```sql
SELECT miint_version();
```

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
- `REFERENCE_LENGTHS` (VARCHAR, required if INCLUDE_HEADER=true): Table or view name containing reference sequences. Must have at least 2 columns: first column = reference name (VARCHAR), second column = reference length (INTEGER/BIGINT). Column names don't matter. Views are fully supported and can include computed columns.
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

### `COPY ... TO '...' (FORMAT NEWICK)`

Write query results to Newick phylogenetic tree format. Reconstructs a tree from tabular node data and serializes to standard Newick format.

**Required columns:**
- `node_index` (BIGINT): Unique identifier for each node
- `parent_index` (BIGINT, nullable): Parent node's node_index (NULL for root)

**Optional columns:**
- `name` (VARCHAR): Node label
- `branch_length` (DOUBLE): Branch length
- `edge_id` (BIGINT): Edge identifier for jplace format

**Parameters:**
- `EDGE_IDS` (BOOLEAN, default: auto): Include edge identifiers `{n}` in output
  - Default: true if `edge_id` column exists, false otherwise
  - Set to `false` to explicitly exclude edge IDs
- `COMPRESSION` (default: auto): Enable gzip compression
  - `'gzip'`: Enable gzip compression
  - `'none'`: No compression
  - Auto-detected from `.gz` extension
- `PLACEMENTS` (VARCHAR): Name of a table containing phylogenetic placement data to insert
  - Inserts fragment sequences into the tree at their placement locations
  - Requires tree to have `edge_id` column with valid edge identifiers
  - Required table columns: `fragment_id` (VARCHAR), `edge_id` (BIGINT/INTEGER), `like_weight_ratio` (DOUBLE), `distal_length` (DOUBLE), `pendant_length` (DOUBLE)
  - Duplicate fragment_ids are deduplicated (keeps highest like_weight_ratio)
  - Multiple placements on same edge are handled correctly

**Behavior:**
- Reconstructs tree structure from parent-child relationships
- Validates tree structure:
  - Exactly one root (node with NULL parent_index)
  - All parent references valid
  - No cycles
  - All nodes connected
- Serializes to standard Newick format with semicolon terminator
- Extra columns (e.g., `is_tip`, `filepath`) are ignored

**Examples:**
```sql
-- Basic roundtrip: read and write tree
COPY (SELECT * FROM read_newick('input.nwk'))
TO 'output.nwk' (FORMAT NEWICK);

-- Write compressed tree
COPY (SELECT * FROM read_newick('input.nwk'))
TO 'output.nwk.gz' (FORMAT NEWICK);

-- Explicit compression parameter
COPY (SELECT * FROM read_newick('input.nwk'))
TO 'output.nwk' (FORMAT NEWICK, COMPRESSION 'gzip');

-- Include edge IDs (for jplace compatibility)
COPY (SELECT * FROM read_newick('reference.nwk'))
TO 'with_edges.nwk' (FORMAT NEWICK, EDGE_IDS true);

-- Exclude edge IDs even if present in input
COPY (SELECT * FROM read_newick('reference.nwk'))
TO 'no_edges.nwk' (FORMAT NEWICK, EDGE_IDS false);

-- Modify tree before writing (e.g., scale branch lengths)
COPY (
    SELECT node_index, name, branch_length * 2.0 AS branch_length,
           edge_id, parent_index
    FROM read_newick('input.nwk')
) TO 'scaled.nwk' (FORMAT NEWICK);

-- Create tree from scratch
CREATE TABLE my_tree AS
SELECT * FROM (VALUES
    (0, NULL::BIGINT, '', 0.0),      -- root
    (1, 0, 'A', 0.1),                 -- tip A
    (2, 0, 'B', 0.2)                  -- tip B
) AS t(node_index, parent_index, name, branch_length);

COPY (SELECT * FROM my_tree)
TO 'new_tree.nwk' (FORMAT NEWICK);
-- Produces: (A:0.1,B:0.2);

-- Filter to subtree (advanced - requires careful node_index management)
-- Note: Ensure parent references remain valid after filtering

-- Insert phylogenetic placements into a reference tree
-- First, create a placement table (e.g., from read_jplace output)
CREATE TABLE placements AS
SELECT * FROM (VALUES
    ('fragment_1', 0::BIGINT, 0.95::DOUBLE, 0.05::DOUBLE, 0.001::DOUBLE),
    ('fragment_2', 1::BIGINT, 0.80::DOUBLE, 0.10::DOUBLE, 0.002::DOUBLE)
) AS t(fragment_id, edge_id, like_weight_ratio, distal_length, pendant_length);

-- Write tree with placements inserted
COPY (SELECT * FROM read_newick('reference.nwk'))
TO 'with_placements.nwk' (FORMAT NEWICK, PLACEMENTS 'placements');

-- Combine with read_jplace for seamless jplace to newick workflow
CREATE TABLE jplace_placements AS
SELECT fragment_id, edge_id, like_weight_ratio, distal_length, pendant_length
FROM read_jplace('placements.jplace');

COPY (SELECT * FROM read_newick('reference.nwk'))
TO 'resolved_tree.nwk' (FORMAT NEWICK, PLACEMENTS 'jplace_placements');
```

**Validation errors:**
- "No data to write" - Empty result set
- "no root" - All nodes have a parent_index
- "multiple roots" - More than one node has NULL parent_index
- "invalid parent reference" - parent_index references non-existent node
- "cycle detected" - Circular parent-child relationship
- "disconnected nodes" - Nodes not reachable from root

**Placement-specific errors (when using PLACEMENTS):**
- "Placement table does not exist" - The specified table doesn't exist
- "Placement table is empty" - The table has no rows
- "missing required column" - Missing one of: fragment_id, edge_id, like_weight_ratio, distal_length, pendant_length
- "Tree has no edge_id values" - Tree data lacks edge identifiers
- "Unknown edge_id" - Placement references an edge_id not present in tree
- "distal_length is negative" - distal_length must be ≥ 0
- "pendant_length is negative" - pendant_length must be ≥ 0
- "distal_length exceeds edge length" - distal_length is larger than the edge's branch_length

**Roundtrip guarantee:**
Trees written with `COPY FORMAT NEWICK` can be read back with `read_newick` and will produce equivalent structure:
```sql
-- Write tree
COPY (SELECT * FROM read_newick('original.nwk'))
TO 'copy.nwk' (FORMAT NEWICK);

-- Read back - should have same structure
SELECT * FROM read_newick('copy.nwk');
```

**Newick Format Output:**
- Node names are included without quotes (unless containing special characters)
- Branch lengths are included after colon (`:0.123`)
- Edge IDs are included in braces when EDGE_IDS=true (`{0}`)
- Tree is terminated with semicolon (`;`)

## Running the tests

The MIINT extension uses three complementary testing approaches to ensure correctness and reliability.

### SQL Logic Tests

The primary testing mechanism for user-facing functionality uses DuckDB's SQL logic test framework. Test files are located in `test/sql/` and use the `.test` extension.

**Running SQL tests:**
```sh
# Run all tests (SQL, C++, and shell)
bash run_tests.sh

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
- `data/newick/` - Newick phylogenetic tree test files

**Important:** Keep test data files small (< 1KB when possible) to maintain fast test execution and minimize repository size.

---

This repository is based on https://github.com/duckdb/extension-template, check it out if you want to build and ship your own DuckDB extension.
