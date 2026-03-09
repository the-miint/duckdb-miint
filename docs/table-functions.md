# Table Functions

Table functions allow querying bioinformatics files as SQL tables.

## Table of Contents

- [`read_alignments`](#read_alignmentsfilename-reference_lengthstable_name-include_filepathfalse-include_seq_qualfalse) - SAM/BAM alignment files
- [`read_fastx`](#read_fastxfilename-sequence2filename-include_filepathfalse-qual_offset33) - FASTA/FASTQ sequence files
- [`read_sequences_sff`](#read_sequences_sfffilename-include_filepathfalse-trimtrue) - SFF (454/Roche) sequence files
- [`read_biom`](#read_biomfilename-include_filepathfalse) - BIOM observation matrix files
- [`read_gff`](#read_gffpath) - GFF3 genome annotation files
- [`read_ncbi`](#read_ncbiaccession-api_key) - NCBI accession metadata
- [`read_ncbi_fasta`](#read_ncbi_fastaaccession-api_key-include_filepathfalse) - NCBI FASTA sequences
- [`read_ncbi_annotation`](#read_ncbi_annotationaccession-api_key-include_filepathfalse) - NCBI genome annotations
- [`read_jplace`](#read_jplacepath) - Phylogenetic placement files
- [`read_newick`](#read_newickfilename-include_filepathfalse) - Newick phylogenetic trees
- [`align_minimap2`](#align_minimap2query_table-subject_tablenull-index_pathnull-options) - Minimap2 alignment
- [`save_minimap2_index`](#save_minimap2_indexsubject_table-output_path-options) - Save minimap2 index
- [`align_minimap2_sharded`](#align_minimap2_shardedquery_table-shard_directory-read_to_shard-options) - Sharded minimap2 alignment
- [`align_bowtie2`](#align_bowtie2query_table-subject_table-options) - Bowtie2 alignment
- [`align_bowtie2_sharded`](#align_bowtie2_shardedquery_table-shard_directory-read_to_shard-options) - Sharded bowtie2 alignment

## `read_alignments(filename, [reference_lengths='table_name'], [include_filepath=false], [include_seq_qual=false])`
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

## `read_fastx(filename, [sequence2=filename], [include_filepath=false], [qual_offset=33])`
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

## `read_sequences_sff(filename, [include_filepath=false], [trim=true])`
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

## `read_biom(filename, [include_filepath=false])`
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

## `read_gff(path)`

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
- **Hierarchical queries**: Join parent-child relationships (gene -> transcript -> exon)
- **Integration**: Combine annotations with alignment data or other genomic datasets

**Compatibility:**
- Supports GFF3 format specification
- Compatible with annotations from:
  - NCBI RefSeq
  - Ensembl
  - GENCODE
  - Other standard GFF3 producers

**Implementation note:** Implemented as a DuckDB macro using `read_csv` with GFF-specific parsing.

## `read_ncbi(accession, [api_key])`

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

## `read_ncbi_fasta(accession, [api_key], [include_filepath=false])`

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

## `read_ncbi_annotation(accession, [api_key], [include_filepath=false])`

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
- Source is detected from accession prefix: NC_/NM_/NP_ -> 'RefSeq', others -> 'GenBank' or 'NCBI'

## `read_jplace(path)`

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

## `read_newick(filename, [include_filepath=false])`

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

## `align_minimap2(query_table, [subject_table=NULL], [index_path=NULL], [options])`

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

## `save_minimap2_index(subject_table, output_path, [options])`

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

## `align_minimap2_sharded(query_table, shard_directory, read_to_shard, [options])`

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

## `align_bowtie2(query_table, subject_table, [options])`

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

## `align_bowtie2_sharded(query_table, shard_directory, read_to_shard, [options])`

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
