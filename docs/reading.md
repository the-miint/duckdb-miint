# Reading files

MIINT supports reading FASTA, FASTQ, SAM, BAM, SFF, BIOM, mzML, mzXML, GFF, jplace and Newick.  

## Table of contents

- [FASTA / FASTQ](#fasta-and-fastq) - FASTA and FASTQ sequence files.
- [SAM / BAM](#sam-and-bam) - SAM and BAM alignment files.
- [SFF](#sff) - SFF sequence files.
- [mzML / mzXML](#mzml-and-mzxml) - mzML and mzXML files.
- [BIOM](#biom) - BIOM-Format v2.1 files.
- [jplace](#jplace) - jplace formatted files.
- [Newick](#newick) - Phylogenies in Newick.
- [GFFv3](#gff) - GFF formatted annotation data.

### FASTA and FASTQ

Reading FASTA and FASTQ sequence files occurs through the `read_fastx` table function.

**Basic example:**

```sql
SELECT * FROM read_fastx('some_file.fasta');
```

**Parameters:**
- `filename` (VARCHAR or VARCHAR[]): Path to FASTA/FASTQ file(s), glob pattern (e.g., `'data/*.fastq'`), or R1 files for paired-end
  - **Glob patterns**: When a single VARCHAR contains glob characters (`*`, `?`, `[`), files are expanded and sorted alphabetically
  - **Arrays**: VARCHAR[] elements are treated as literal paths (no glob expansion)
- `sequence2` (VARCHAR or VARCHAR[], optional): Path to R2 file(s) for paired-end reads. Must have same number of files as `filename`
  - **Paired-end with globs**: When `filename` is a glob pattern, `sequence2` must also be a glob pattern. Both are expanded and sorted independently, then paired by position. The expanded file counts must match.
- `include_filepath` (BOOLEAN, optional, default false): Add filepath column to output
- `qual_offset` (INTEGER, optional, default 33): Quality score offset (33 for Phred+33, 64 for Phred+64)
- `max_batch_bytes` (VARCHAR, optional, default `'512MiB'`): Soft cap on the uncompressed sequence+quality bytes buffered per output chunk, given as a formatted byte size (e.g. `'256MiB'`, `'2GB'`). The reader stops adding records to a chunk once their combined bytes reach this budget, which bounds memory when individual records are very large (e.g. assembled genomes at multiple MB each). Short-read data is unaffected — a chunk fills to the standard vector size long before the cap. Must be greater than 0.

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
- Streaming I/O minimizes memory usage
- Quality scores stored as efficient UINT8 arrays

**Notes:**
- Read IDs must match between R1 and R2 files for paired-end data (not validated, matched by position)
- For FASTA files, `qual1` and `qual2` are NULL
- The `sequence_index` resets to 1 for each file. To distinguish sequences from different files, use `include_filepath=true` and order by `filepath, sequence_index`
- Comment field is NULL if no comment is present in the sequence header

### SAM and BAM

Reading SAM/BAM alignment files occurs through the `read_alignments` table function. Header information is optional. SEQ and QUAL fields are by default not included. 

Sequence data can be read from SAM/BAM into a `read_fastx` compatible schema using `read_sequences_sam`. 

**Basic example:**

```sql
SELECT * FROM read_alignments('some_alignments.sam')
```

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

### SFF

Read SFF (Standard Flowgram Format) files from 454/Roche sequencing platforms.

**Basic example:**

```sql
SELECT * FROM read_sequences_sff('some_file.sff');
```

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

### mzML and mzXML

Read mzML and mzXML mass spectrometry files and chromatograms. Returns one row per spectrum with metadata and peak arrays. 

**Basic example:**

```
SELECT * FROM read_mzml('some_file.mzml');
SELECT * FROM read_mzxml('some_file.mzxml');
SELECT * FROM read_mzml_chromatograms('some_file.mzml');
```

**Parameters:**

Note the parameters for `read_mzml` and `read_mzxml` are the same.

- `filename` (VARCHAR or VARCHAR[]): Path to mzML file(s) or glob pattern
- `include_filepath` (BOOLEAN, default false): Add a `filepath` column

**Output columns (27):**

Note the output schema for `read_mzml` and `read_mzxml` are the same.

| Column | Type | Description |
|--------|------|-------------|
| `spectrum_index` | INTEGER | 0-based spectrum position |
| `spectrum_id` | VARCHAR | Spectrum identifier from file |
| `scan_number` | INTEGER | Scan number (extracted from spectrum_id, nullable) |
| `ms_level` | INTEGER | MS level (1, 2, ...) |
| `retention_time` | DOUBLE | Retention time in minutes (nullable) |
| `spectrum_type` | VARCHAR | "centroid" or "profile" (nullable) |
| `polarity` | VARCHAR | "positive" or "negative" (nullable) |
| `base_peak_mz` | DOUBLE | m/z of base peak (nullable) |
| `base_peak_intensity` | DOUBLE | Intensity of base peak (nullable) |
| `total_ion_current` | DOUBLE | Total ion current (nullable) |
| `lowest_mz` | DOUBLE | Lowest observed m/z (nullable) |
| `highest_mz` | DOUBLE | Highest observed m/z (nullable) |
| `default_array_length` | INTEGER | Number of peaks |
| `precursor_mz` | DOUBLE | Precursor m/z for MS2+ (nullable) |
| `precursor_charge` | INTEGER | Precursor charge state (nullable) |
| `precursor_intensity` | DOUBLE | Precursor intensity (nullable) |
| `isolation_window_target` | DOUBLE | Isolation window center (nullable) |
| `isolation_window_lower` | DOUBLE | Isolation window lower offset (nullable) |
| `isolation_window_upper` | DOUBLE | Isolation window upper offset (nullable) |
| `activation_method` | VARCHAR | "CID", "HCD", "ETD" (nullable) |
| `collision_energy` | DOUBLE | Collision energy (nullable) |
| `mz_array` | DOUBLE[] | Array of m/z values |
| `intensity_array` | DOUBLE[] | Array of intensity values |
| `filter_string` | VARCHAR | Instrument filter string (nullable) |
| `scan_window_lower` | DOUBLE | Scan window lower bound (nullable) |
| `scan_window_upper` | DOUBLE | Scan window upper bound (nullable) |
| `ms1_scan_index` | INTEGER | spectrum_index of parent MS1 scan (nullable) |

Chromatogram output schema:

**Output columns (7):**

| Column | Type | Description |
|--------|------|-------------|
| `chromatogram_index` | INTEGER | 0-based chromatogram position |
| `chromatogram_id` | VARCHAR | Chromatogram identifier |
| `chromatogram_type` | VARCHAR | "TIC", "BPC", "SRM", "SIC" (nullable) |
| `precursor_mz` | DOUBLE | Precursor m/z (nullable) |
| `product_mz` | DOUBLE | Product m/z (nullable) |
| `time_array` | DOUBLE[] | Array of time values |
| `intensity_array` | DOUBLE[] | Array of intensity values |

**Examples:**

The examples below are interchangable with `read_mzxml`. 

```sql
-- Count spectra by MS level
SELECT ms_level, COUNT(*) FROM read_mzml('sample.mzML') GROUP BY ms_level;

-- Get MS2 spectra with precursor info
SELECT scan_number, precursor_mz, precursor_charge, activation_method
FROM read_mzml('sample.mzML')
WHERE ms_level = 2;

-- Unnest peak arrays for peak-level analysis
SELECT s.spectrum_index, s.ms_level, unnest(s.mz_array) AS mz, unnest(s.intensity_array) AS intensity
FROM read_mzml('sample.mzML') s;

-- Read multiple files
SELECT * FROM read_mzml('data/mzml/*.mzML', include_filepath=true);
```

**mzML Notes:**
- Stdin is not supported (mzML requires file seeking)
- Supports zlib-compressed and uncompressed binary arrays, 32-bit and 64-bit precision
- Also reads chromatograms via `read_mzml_chromatograms`
- See [Mass Spectrometry & MassQL](massql.md) for higher-level analysis macros

**mzXML Notes:**
- Handles big-endian interleaved binary data (mzXML spec)
- Supports nested `<scan>` elements (child scans emitted before parents)
- `scan_number` always populated from the `<scan num>` attribute
- Supports zlib compression, 32-bit and 64-bit precision

### BIOM

Read BIOM (Biological Observation Matrix) format files.

**Basic exmaple:**

```sql
SELECT * FROM read_biom('some_file.biom');
```

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

### jplace

Read jplace phylogenetic placement files. The jplace format stores query sequence placements onto a reference phylogenetic tree, as defined in [Matsen et al. 2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009). We support reading the placement data, and separately, reading the phylogeny stored in `jplace`.
These methods are differentiated as the output schemas are not the same.

**Basic example:**

```sql
SELECT * FROM read_jplace('some_data.jplace');
SELECT * FROM read_jplace_newick('some_data.jplace');
```

**Parameters:**
- `path` (VARCHAR): Path to jplace file(s), supports glob patterns (e.g., `'data/*.jplace'`)

**Output schema:**

For `read_jplace`, the output schema is:

- `fragment` (VARCHAR): Fragment/sequence name (from `nm` or `n` field)
- `edge_num` (INTEGER): Edge number in the reference tree where placement occurs
- `likelihood` (DOUBLE): Log likelihood of the placement
- `like_weight_ratio` (DOUBLE): Likelihood weight ratio (proportion of total likelihood)
- `distal_length` (DOUBLE): Distance from distal end of edge to placement point
- `pendant_length` (DOUBLE): Pendant branch length (branch to the placed sequence)
- `filepath` (VARCHAR): Source file path

For `read_jplace_newick`, the output schema is consistent with `read_newick`:

- `node_index` (BIGINT): 0-based index of node in tree
- `name` (VARCHAR): Node label (empty string for unlabeled internal nodes)
- `branch_length` (DOUBLE, nullable): Branch length
- `edge_id` (BIGINT, nullable): Edge identifier from `{n}` syntax
- `parent_index` (BIGINT, nullable): Parent node's node_index (NULL for root)
- `is_tip` (BOOLEAN): Whether node is a tip/leaf
- `filepath` (VARCHAR, optional): File path when include_filepath=true

**Behavior:**

`read_jplace` specific items:

- Returns only the best placement (first in `p` array) for each fragment
- Supports both `nm` (named multiplicities) and `n` (names) formats

`read_jplace_newick` specific items:

- Reads the `"tree"` field from each jplace JSON file and parses it as Newick
- Schema is UNION ALL-compatible with `read_newick`

General items:

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

-- Extract the reference tree from a jplace file
SELECT * FROM read_jplace_newick('placements.jplace');

-- Get tip names from the reference tree
SELECT name FROM read_jplace_newick('placements.jplace')
WHERE is_tip = true;

-- Combine with placement data for a full workflow
CREATE TABLE ref_tree AS
SELECT * FROM read_jplace_newick('placements.jplace');

CREATE TABLE placements AS
SELECT fragment AS fragment_id, edge_num::BIGINT AS edge_id,
       like_weight_ratio, distal_length, pendant_length
FROM read_jplace('placements.jplace');

-- Now use tree_resolve_placement to insert placements into the tree
SELECT * FROM tree_resolve_placement('ref_tree', 'placements');

-- Compare trees across multiple jplace files
SELECT filepath, COUNT(*) AS num_nodes,
       COUNT(*) FILTER (WHERE is_tip) AS num_tips
FROM read_jplace_newick('results/*.jplace', include_filepath=true)
GROUP BY filepath;
```

**jplace Format Notes:**
- Standard format for phylogenetic placement tools (pplacer, EPA-ng, etc.)
- The `fields` array in the file defines column order: typically `[edge_num, likelihood, like_weight_ratio, distal_length, pendant_length]`
- Multiple placements per fragment are supported in the format, but this function returns only the best (first) placement
- The `nm` field contains `[[name, multiplicity], ...]` pairs; `n` field contains simple name arrays

**Compatibility:**
- Reads jplace files created by:
  - pplacer
  - EPA-ng
  - SEPP
  - Other tools following the jplace specification

**Implementation note:** `read_jplace` is implemented as a DuckDB macro using `read_json` with JSON path extraction.

### Newick

Read Newick phylogenetic tree files and return a table with one row per node.

**Basic example:**

```sql
SELECT * FROM read_newick('some_file.newick');
```

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

**Newick Format Notes:**
- Standard format for representing phylogenetic trees as text
- Parentheses denote tree structure: `(A,B)` means A and B share a parent
- Colons precede branch lengths: `A:0.1` means tip A with branch length 0.1
- Semicolons terminate the tree: `(A,B);`
- Edge IDs in braces are an extension for jplace: `A:0.1{0}` has edge_id 0

### GFF

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
- Percent-decodes (`%XX`) attribute keys and values per the GFF3 spec (e.g. `%3D` → `=`, `%2C` → `,`, `%3B` → `;`, `%26` → `&`); the value is everything after the first `=`, so a raw `=` inside a value is preserved rather than truncated
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

#### `parse_gff_attributes(kvp_string)`

The parsing above is done by `parse_gff_attributes`, which `read_gff` applies to column 9. It is also callable directly on any GFF-style attribute string.

**Parameters:**
- `kvp_string` (VARCHAR): A semicolon-separated `key=value` string (a GFF3 column-9 value).

**Returns:** `MAP(VARCHAR, VARCHAR)` — the parsed key/value pairs.

```sql
SELECT parse_gff_attributes('ID=gene1;Name=TEST1;biotype=protein_coding');
-- {'ID': 'gene1', 'Name': 'TEST1', 'biotype': 'protein_coding'}
```
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
