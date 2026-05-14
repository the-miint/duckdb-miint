# Table Functions

Table functions allow querying bioinformatics files as SQL tables.

## Table of Contents

- [`read_alignments`](#read_alignmentsfilename-reference_lengthstable_name-include_filepathfalse-include_seq_qualfalse) - SAM/BAM alignment files
- [`read_fastx`](#read_fastxfilename-sequence2filename-include_filepathfalse-qual_offset33) - FASTA/FASTQ sequence files
- [`read_sequences_sff`](#read_sequences_sfffilename-include_filepathfalse-trimtrue) - SFF (454/Roche) sequence files
- [`read_mzml`](#read_mzmlfilename-include_filepathfalse) - mzML mass spectrometry files
- [`read_mzxml`](#read_mzxmlfilename-include_filepathfalse) - mzXML mass spectrometry files
- [`read_mzml_chromatograms`](#read_mzml_chromatogramsfilename-include_filepathfalse) - mzML chromatogram data
- [`read_biom`](#read_biomfilename-include_filepathfalse) - BIOM observation matrix files
- [`read_gff`](#read_gffpath) - GFF3 genome annotation files
- [`read_ncbi`](#read_ncbiaccession-api_key) - NCBI accession metadata
- [`read_ncbi_fasta`](#read_ncbi_fastaaccession-api_key-include_filepathfalse) - NCBI FASTA sequences
- [`read_ncbi_annotation`](#read_ncbi_annotationaccession-api_key-include_filepathfalse) - NCBI genome annotations
- [`read_ena`](#read_enaaccession-resultread_run-fields) - EBI/ENA metadata queries
- [`read_ena_attributes`](#read_ena_attributesaccession) - EBI/ENA custom sample attributes
- [`ena_searchable_fields`](#ena_searchable_fieldsresult_type) - Enumerate ENA structured-search fields per result type
- [`read_ena_sequences`](#read_ena_sequencesaccession-include_filepathfalse-qual_offset33-max_sequences0) - Stream FASTA/FASTQ from EBI/ENA with accession columns
- [`read_jplace`](#read_jplacepath) - Phylogenetic placement files
- [`read_jplace_newick`](#read_jplace_newickpath-include_filepathfalse) - Newick tree from jplace files
- [`read_newick`](#read_newickfilename-include_filepathfalse) - Newick phylogenetic trees
- [`tree_resolve_placement`](#tree_resolve_placementtree_table-placements_table) - Resolve phylogenetic placements into a tree
- [`phylogeny_fasttree`](#phylogeny_fasttreetable_name-options) - Approximately-maximum-likelihood phylogenetic tree (FastTree via gpl-boundary daemon)
- [`align_minimap2`](#align_minimap2query_table-subject_tablenull-index_pathnull-options) - Minimap2 alignment
- [`save_minimap2_index`](#save_minimap2_indexsubject_table-output_path-options) - Save minimap2 index
- [`align_minimap2_sharded`](#align_minimap2_shardedquery_table-shard_directory-read_to_shard-options) - Sharded minimap2 alignment
- [`align_bowtie2`](#align_bowtie2query_table-subject_table-options) - Bowtie2 alignment
- [`align_bowtie2_sharded`](#align_bowtie2_shardedquery_table-shard_directory-read_to_shard-options) - Sharded bowtie2 alignment
- [`align_mafft`](#align_maffttable_name) - MAFFT multiple sequence alignment (PartTree)
- [`align_sortmerna`](#align_sortmernaquery_table-ref_pathspaths-options) - SortMeRNA rRNA alignment (SAM schema)
- [`align_sortmerna_rrna`](#align_sortmerna_rrnaquery_table-ref_pathspaths-options) - SortMeRNA rRNA alignment (rRNA schema with identity/coverage/e-value)
- [`detect_chimera_uchime`](#detect_chimera_uchimequery_table-dbrefs_table-options) - Reference-based chimera detection (UCHIME)
- [`detect_chimera_uchime_denovo`](#detect_chimera_uchime_denovoinput_table-options) - De novo chimera detection (UCHIME)
- [`search_sequences_vsearch`](#search_sequences_vsearchquery_table-dbref_table-idthreshold-options) - Global sequence search
- [`cluster_sequences_vsearch`](#cluster_sequences_vsearchinput_table-idthreshold-options) - Greedy sequence clustering
- [`deblur`](#deblurinput_table-options) - Deblur amplicon sequence denoising
- [`alignment_slice`](#alignment_slicetable_name-start-stop-include_deletionsfalse) - Slice alignments to a genomic region
- [`unifrac_pcoa`](#unifrac_pcoaobservations-tree-options) - UniFrac distance + Principal Coordinates Analysis
- [`unifrac_permanova`](#unifrac_permanovaobservations-tree-metadata-options) - UniFrac distance + PERMANOVA pseudo-F + p-value
- [`unifrac_faith_pd`](#unifrac_faith_pdobservations-tree-options) - Faith's phylogenetic diversity per sample
- [`miint_warnings`](#miint_warnings) - Query miint's operational warnings as a table

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

## `alignment_slice(table_name, start, stop, [include_deletions=false])`

Slice alignment data from a table or view to a genomic region. Each alignment is trimmed to the specified region `[start, stop)`, with trimmed portions represented as hard clips (H) in the output CIGAR.

**Parameters:**
- `table_name` (VARCHAR): Name of a table or view containing alignment data
- `start` (BIGINT): Region start position (1-based, inclusive)
- `stop` (BIGINT): Region stop position (1-based, exclusive)
- `include_deletions` (BOOLEAN, default `false`): Whether deletion (D) operations count as overlap evidence. When `false`, reads whose only overlap with the region is through deletions are excluded.

**Coordinates:** 1-based half-open `[start, stop)`, consistent with `position` and `stop_position` from `read_alignments`.

**Required input columns:** `cigar` (VARCHAR), `position` (BIGINT), `stop_position` (BIGINT)

**Output schema:** Same columns as found in the input table (from the recognized alignment column set), with adjusted values for overlapping reads.

**Behavior:**
- Reads that don't overlap the region are excluded
- Soft clips (S) do not count as overlap evidence
- Trimmed portions of the CIGAR become hard clips (H)
- Tags (`tag_as` through `tag_sa`) are set to NULL when trimming occurs
- `template_length` is set to NULL when trimming occurs
- `mapq` and mate fields are preserved
- If the input table has a `reference` column, all rows must have the same reference (single-region slicing)
- Rows with NULL `cigar`, `position`, or `stop_position` are skipped

**Examples:**

```sql
-- Load alignments and filter to a single reference
CREATE VIEW chr1_alns AS
  SELECT * FROM read_alignments('sample.bam')
  WHERE reference = 'chr1';

-- Slice to a region of interest
SELECT * FROM alignment_slice('chr1_alns', 1000, 2000);

-- Include reads that overlap only via deletions
SELECT * FROM alignment_slice('chr1_alns', 1000, 2000, include_deletions := true);

-- Composable: compute coverage depth on sliced alignments
SELECT compute_coverage_depth(position, stop_position, cigar, 1000, 'exclude_deletions')
FROM alignment_slice('chr1_alns', 1000, 2000);
```

**Notes:**
- The function reads from the input table via a separate connection, so uncommitted changes in the current transaction may not be visible
- Multi-region slicing (different regions per reference) is not yet supported; use separate queries per region
- `alignment_seq_identity` with the `'cigar'` method works on sliced output because it reads identity directly from `=`/`X` CIGAR ops without needing tags. Other methods (`gap_compressed`, `blast`, `gap_excluded`) require NM or MD tags which are NULLed after trimming.

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

## `read_mzml(filename, [include_filepath=false])`

Read mzML mass spectrometry files. Returns one row per spectrum with metadata and peak arrays.

**Parameters:**
- `filename` (VARCHAR or VARCHAR[]): Path to mzML file(s) or glob pattern
- `include_filepath` (BOOLEAN, default false): Add a `filepath` column

**Output columns (27):**

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

**Examples:**
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

**Notes:**
- Stdin is not supported (mzML requires file seeking)
- Supports zlib-compressed and uncompressed binary arrays, 32-bit and 64-bit precision
- Also reads chromatograms via `read_mzml_chromatograms`
- See [Mass Spectrometry & MassQL](massql.md) for higher-level analysis macros

## `read_mzxml(filename, [include_filepath=false])`

Read mzXML mass spectrometry files. Returns the **same 27-column schema** as `read_mzml`, enabling `UNION ALL` across formats.

**Parameters:** Same as `read_mzml`.

**Output columns:** Identical to `read_mzml` (see above).

**Examples:**
```sql
-- Read mzXML file
SELECT * FROM read_mzxml('sample.mzXML');

-- Combine mzML and mzXML files
SELECT * FROM read_mzml('batch1.mzML')
UNION ALL
SELECT * FROM read_mzxml('batch2.mzXML');
```

**Notes:**
- Handles big-endian interleaved binary data (mzXML spec)
- Supports nested `<scan>` elements (child scans emitted before parents)
- `scan_number` always populated from the `<scan num>` attribute
- Supports zlib compression, 32-bit and 64-bit precision

## `read_mzml_chromatograms(filename, [include_filepath=false])`

Read chromatogram data from mzML files (TIC, BPC, SRM, SIC).

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

## `read_ena(accession, [result='read_run'], [fields])`

Query metadata from the EBI European Nucleotide Archive (ENA) Portal API. Returns run, sample, or study metadata as a table of VARCHAR columns.

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to EBI servers (www.ebi.ac.uk)

**Parameters:**
- `accession` (VARCHAR or VARCHAR[]): ENA/SRA accession(s). Supports study (PRJNA\*, PRJEB\*, ERP\*, SRP\*), sample (SAMN\*, SAME\*), run (SRR\*, ERR\*), and experiment (SRX\*, ERX\*) accessions.
- `result` (VARCHAR, optional, default 'read_run'): ENA result type. One of: `read_run`, `sample`, `study`. When the accession type doesn't match the result type (e.g., run accession with `result='study'`), the accession is automatically resolved.
- `fields` (VARCHAR, optional): Comma-separated list of ENA field names to return. If omitted, a sensible default set is used per result type.

**Output schema:**
- All columns are VARCHAR, named according to the requested fields
- Default `read_run` fields include: `run_accession`, `experiment_accession`, `sample_accession`, `study_accession`, `fastq_ftp`, `fastq_bytes`, `fastq_md5`, `library_strategy`, `library_layout`, `instrument_model`, `read_count`, `base_count`, and more
- Default `sample` fields include: `sample_accession`, `scientific_name`, `tax_id`, `collection_date`, `country`, `lat`, `lon`, `depth`, and more
- Default `study` fields include: `study_accession`, `study_title`, `study_description`, `center_name`, and more

**Behavior:**
- Queries the ENA Portal API (`/search` endpoint) with `format=tsv`
- Accession type is auto-detected from prefix and mapped to the appropriate query parameter
- Cross-type resolution: e.g., a run accession with `result='study'` first resolves to the study accession, then queries the study
- Rate-limited to ~3 requests/second with retry on 429/500/502/503

**Examples:**
```sql
-- Get all run metadata for a single run
SELECT * FROM read_ena('ERR1074767');

-- Get specific fields
SELECT run_accession, library_layout, read_count
FROM read_ena('ERR1074767', fields='run_accession,library_layout,read_count');

-- Get sample metadata (auto-resolves run -> sample)
SELECT * FROM read_ena('ERR1074767', result='sample');

-- Get study info
SELECT study_title FROM read_ena('PRJEB11419', result='study');

-- Get all runs for a BioProject
CREATE TABLE project_runs AS SELECT * FROM read_ena('PRJNA555783');
```

## `read_ena_attributes(accession)`

Fetch custom sample attributes from EBI/ENA via the Browser XML API. Returns all submitter-defined key-value attributes that are not available through the Portal API's fixed schema.

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to EBI servers (www.ebi.ac.uk)

**Parameters:**
- `accession` (VARCHAR or VARCHAR[]): ENA/SRA accession(s). Supports study, sample, run, and experiment accessions. Non-sample accessions are automatically resolved to their associated sample accession(s).

**Output schema:**
- `sample_accession` (VARCHAR): BioSample accession (SAMN\*/SAME\*)
- `tag` (VARCHAR): Attribute name (e.g., 'collection date', 'geographic location', custom fields)
- `value` (VARCHAR): Attribute value

**Behavior:**
- Resolves non-sample accessions to sample accessions via the Portal API
- Fetches XML from the Browser API in batches of 50 samples
- Parses `<SAMPLE_ATTRIBUTE>` `<TAG>`/`<VALUE>` pairs
- Returns ALL attributes including custom/submitter-defined ones (e.g., primer sequences, custom identifiers) that are not available via `read_ena`

**Predicate pushdown:**
When the `WHERE` clause references `tag` (and optionally `value`) with equality-only operators, and every referenced tag is a known ENA sample-search field, the scan switches from per-sample XML fetches to a single `/search?result=sample` TSV request per batch. This converts an O(N)-per-sample XML scan into one HTTP call per 200 samples and can cut minutes off large studies (e.g., a 33 000-sample study drops from a ~3.7-minute rate-limit floor to a few seconds for `WHERE tag='host_body_site'`).

Pushdown triggers when:
- The filter is a conjunction (AND) of `tag = 'X'`, `tag IN ('X','Y',...)`, or `tag = 'X' AND value = 'Y'`
- Every referenced tag passes `ena_searchable_fields('sample')` (see the curated allowlist under the hood; uppercase/mixed-case names are accepted)

Pushdown is **declined** (falls back to the XML path, preserving correctness) when:
- Any referenced tag is not in the searchable-field allowlist (including any single unknown tag inside an `IN` list)
- Any `OR`, `LIKE`, `!=`, `NOT IN`, etc. appears anywhere in the filter tree
- `value` is constrained without a single pinned `tag` (ambiguous)
- Any predicate is on `sample_accession` or another non-`tag`/`value` column

DuckDB always re-applies the original filter above the scan, so the output of the pushdown path is guaranteed to be a subset of what the XML path would have returned — a pushdown mistake degrades to extra work, never to wrong rows.

**Examples:**
```sql
-- Get all attributes for a run's sample
SELECT * FROM read_ena_attributes('ERR1074767');

-- Get attributes for a specific sample
SELECT * FROM read_ena_attributes('SAMEA3610311');

-- Find primer information
SELECT tag, value FROM read_ena_attributes('ERR1074767')
WHERE tag LIKE '%primer%' OR tag LIKE '%Primer%';

-- Pivot attributes to wide format for a study
SELECT sample_accession,
       MAX(CASE WHEN tag = 'collection date' THEN value END) AS collection_date,
       MAX(CASE WHEN tag = 'geographic location (country and/or sea)' THEN value END) AS country
FROM read_ena_attributes('PRJEB11419')
GROUP BY sample_accession;

-- Pushdown-enabled: uses /search endpoint (single TSV per batch) instead of
-- per-sample XML. Fast on large studies.
SELECT sample_accession
FROM read_ena_attributes('PRJEB11419')
WHERE tag='host_body_site' AND value='UBERON:feces';
```

## `ena_searchable_fields(result_type)`

Enumerate the fields that ENA's Portal API `/search?result=<result_type>` endpoint accepts as structured filters. Use this to discover what field names are available for a given result type (sample, read_run, study, experiment, etc.) before constructing a query — either a direct Portal API `/search` URL or a `WHERE tag='X'` filter on `read_ena_attributes` that can be pushed down to the structured search.

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to EBI servers (www.ebi.ac.uk)

**Parameters:**
- `result_type` (VARCHAR, required): An ENA Portal API result type (e.g., `'sample'`, `'read_run'`, `'study'`, `'experiment'`)

**Output schema:**
- `field_name` (VARCHAR): Canonical field identifier accepted by ENA's `/search?fields=...` query parameter
- `type` (VARCHAR): ENA's declared type for the field (e.g., `text`, `number`, `date`, `controlled value`, `taxonomy`)
- `description` (VARCHAR, nullable): ENA's human-readable description, when available

**Behavior:**
- Issues exactly one HTTP call to `/returnFields?result=<result_type>&format=tsv` on first use and caches the parsed result for the remainder of the scan
- Validates `result_type` as alphanumeric + `_-.` to prevent URL injection

**Examples:**
```sql
-- All searchable fields for the sample result type
SELECT field_name, type FROM ena_searchable_fields('sample') ORDER BY field_name;

-- Check whether a specific field exists before using it in a filter
SELECT COUNT(*) > 0 AS exists
FROM ena_searchable_fields('sample')
WHERE field_name = 'host_body_site';
```

## `read_ena_sequences(accession, [include_filepath=false], [qual_offset=33], [max_sequences=0])`

Stream FASTA/FASTQ sequence data from EBI/ENA with run, sample, and experiment accession columns. Returns data in the same schema as `read_fastx` plus accession metadata columns, enabling direct association of sequence data with project metadata. Supports mixed FASTA/FASTQ paired-end data (unlike `read_fastx` which rejects format mismatches).

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to EBI servers (ftp.sra.ebi.ac.uk via HTTPS)

**Parameters:**
- `accession` (VARCHAR or VARCHAR[]): ENA/SRA accession(s). Supports study (bulk download all runs), sample, run, and experiment accessions.
- `include_filepath` (BOOLEAN, optional, default false): Add filepath column with the HTTPS download URL(s). For paired-end runs, URLs are semicolon-separated.
- `qual_offset` (BIGINT, optional, default 33): Quality score offset (33 for Phred+33/Sanger, 64 for Phred+64/Illumina 1.3+)
- `max_sequences` (BIGINT, optional, default 0): If `> 0`, stop emitting from each run after this many sequences. `0` (or NULL / absent) means unlimited. For paired-end runs the cap counts **pairs** (one output row per pair), not underlying FASTQ records — `max_sequences=N` yields at most N rows and corresponds to 2N downloaded reads. When downloading via Aspera the cap tears down the `ascp` transfer early, saving real bandwidth. For SFF runs the cap applies but the full file is downloaded before any record is parsed; a loud warning is printed in that case.
- `trim_sff` (BOOLEAN, optional, default true): For SFF runs, apply the quality and adapter clip positions from the SFF header to trim sequences and quality scores. Ignored for FASTQ runs. Named `trim_sff` rather than `trim` because `TRIM` is a SQL function keyword and this function is dual-path (supports both scalar and lateral invocation), which together prevent DuckDB's binder from accepting `trim=...`.

Because `read_ena_sequences` supports lateral / correlated invocation, named parameters must be passed with arrow syntax (`name => value`), not `name = value`. For example: `read_ena_sequences('X', prefer_format => 'sff', trim_sff => false)`.

**Output schema:**
- `sequence_index` (BIGINT): 1-based sequence index (per run)
- `read_id` (VARCHAR): FASTQ read identifier
- `comment` (VARCHAR): FASTQ comment line (nullable)
- `sequence1` (VARCHAR): Forward read sequence
- `sequence2` (VARCHAR): Reverse read sequence (NULL for single-end)
- `qual1` (UTINYINT[]): Forward quality scores
- `qual2` (UTINYINT[]): Reverse quality scores (NULL for single-end)
- `run_accession` (VARCHAR): SRR/ERR run accession
- `sample_accession` (VARCHAR): SAMN/SAME BioSample accession
- `experiment_accession` (VARCHAR): SRX/ERX experiment accession
- `filepath` (VARCHAR, optional): HTTPS download URL(s)

**Behavior:**
- At bind time, queries ENA Portal API to discover runs, FASTQ URLs, and library layout (paired vs single-end)
- Constructs HTTPS URLs from ENA's FTP paths
- Streams FASTQ data through the same `SequenceReader`/`DuckDBSeqStream` infrastructure as `read_fastx`
- Paired-end detection is automatic from the `library_layout` field
- Supports parallel streaming across multiple runs (up to 8 concurrent)

**Examples:**
```sql
-- Stream reads from a single run
SELECT * FROM read_ena_sequences('ERR1074767') LIMIT 100;

-- Stream all reads for a BioProject directly to Parquet
COPY (SELECT * FROM read_ena_sequences('PRJNA555783'))
  TO 'project_sequences.parquet' (FORMAT PARQUET);

-- Join sequence data with metadata
CREATE TABLE runs AS SELECT * FROM read_ena('PRJNA555783');
CREATE TABLE seqs AS SELECT * FROM read_ena_sequences('PRJNA555783');

SELECT s.run_accession, r.library_strategy, COUNT(*) AS read_count
FROM seqs s
JOIN runs r USING (run_accession)
GROUP BY ALL;

-- Stream with filepath for provenance tracking
SELECT run_accession, read_id, filepath
FROM read_ena_sequences('ERR1074767', include_filepath=true) LIMIT 5;
```

**Notes:**
- For large projects, the FASTQ download can take significant time and bandwidth
- The `run_accession` column enables easy JOIN back to metadata from `read_ena`
- Quality scores are converted to numeric values using the specified offset (default Phred+33)
- **Failure handling**: on a transient open failure, the run is retried once. On
  a mid-stream failure (connection dropped after some rows have been emitted),
  the run is NOT retried — re-reading from scratch would emit the same reads
  again with duplicate `sequence_index` values downstream. Instead, the run is
  recorded as skipped, a loud warning reports the run accession and the
  number of reads emitted before the failure, and the scan continues with
  other runs. If you see such a warning, re-run the query (the metadata
  lookup is cached) or use a smaller per-run selection to recover the
  truncated data.

**Lateral invocation (correlated arguments):**

Accessions can come from a correlated column via `LATERAL`. Each outer row opens its
own run; a `LIMIT` inside the lateral short-circuits that run's download (the reader
and any HTTP connection are torn down as soon as the limit is reached). Metadata
lookups share an in-memory LRU cache scoped to the query, so repeated outer-row
values and within-query batches avoid redundant ENA Portal API calls.

```sql
-- Find runs containing a probe sequence, downloading only until each match is seen.
SELECT r.run_accession
FROM read_ena('PRJEB11419') AS r,
     LATERAL (SELECT 1 FROM read_ena_sequences(r.run_accession)
              WHERE sequence1.contains('AACGTAGGTCACAAGCGTTGTCCGGA')
              LIMIT 1);
```

Limitations in lateral mode:
- `download_method='aspera'` is not supported (use HTTP; the lateral use case is
  short-circuit-driven, not throughput-driven).
- `table_scan_progress` reports `-1.0` (indeterminate) because the total work is
  driven by the outer side and not known at bind time.

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

## `read_jplace_newick(path, [include_filepath=false])`

Extract the reference Newick tree from jplace phylogenetic placement files. Returns the tree in the same schema as `read_newick`, with one row per node.

**Parameters:**
- `path` (VARCHAR or VARCHAR[]): Path to jplace file(s), supports glob patterns (e.g., `'data/*.jplace'`)
- `include_filepath` (BOOLEAN, optional, default false): Add filepath column to output

**Output schema:** Same as [`read_newick`](#read_newickfilename-include_filepathfalse):
- `node_index` (BIGINT): 0-based index of node in tree
- `name` (VARCHAR): Node label (empty string for unlabeled internal nodes)
- `branch_length` (DOUBLE, nullable): Branch length
- `edge_id` (BIGINT, nullable): Edge identifier from `{n}` syntax
- `parent_index` (BIGINT, nullable): Parent node's node_index (NULL for root)
- `is_tip` (BOOLEAN): Whether node is a tip/leaf
- `filepath` (VARCHAR, optional): File path when include_filepath=true

**Behavior:**
- Reads the `"tree"` field from each jplace JSON file and parses it as Newick
- Schema is UNION ALL-compatible with `read_newick`
- Supports glob patterns for reading trees from multiple jplace files
- Supports gzip-compressed jplace files

**Examples:**
```sql
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

## `tree_resolve_placement(tree_table, placements_table)`

Resolve phylogenetic placements into a reference tree, returning a fully resolved tree with placed fragments as new tips. This exposes the `insert_fully_resolved` algorithm as a SQL-accessible table function.

**Parameters:**
- `tree_table` (VARCHAR): Name of a table or view containing tree data in `read_newick` schema (requires `node_index` and `parent_index` columns; `name`, `branch_length`, `edge_id` are optional)
- `placements_table` (VARCHAR): Name of a table or view containing placement data (requires `fragment_id`, `edge_id`, `like_weight_ratio`, `distal_length`, `pendant_length` columns)

**Output schema:** Same as [`read_newick`](#read_newickfilename-include_filepathfalse) (without filepath):
- `node_index` (BIGINT): 0-based index of node in resolved tree
- `name` (VARCHAR): Node label (placed fragments use their fragment_id as name)
- `branch_length` (DOUBLE, nullable): Branch length
- `edge_id` (BIGINT, nullable): Edge identifier (NULL for newly created nodes)
- `parent_index` (BIGINT, nullable): Parent node's node_index (NULL for root)
- `is_tip` (BOOLEAN): Whether node is a tip/leaf

**Behavior:**
- Reads tree data from the tree table and builds a `NewickTree`
- Reads placements and inserts each fragment as a new tip on the specified edge
- Each placement creates 2 new nodes: an internal node (splitting the edge) and a fragment tip
- Deduplicates placements by `fragment_id` (keeps highest `like_weight_ratio`, then lowest `pendant_length`)
- Multiple placements on the same edge are sorted by `distal_length` and inserted as a chain
- Preserves original tip-to-tip distances in the tree
- Works with both tables and views for either parameter
- Schema is UNION ALL-compatible with `read_newick`

**Examples:**
```sql
-- Basic workflow: load tree, load placements, resolve
CREATE TABLE ref_tree AS
SELECT * FROM read_newick('reference.nwk');

CREATE TABLE placements AS
SELECT * FROM (VALUES
    ('seq1', 0::BIGINT, 0.95::DOUBLE, 0.05::DOUBLE, 0.001::DOUBLE),
    ('seq2', 1::BIGINT, 0.80::DOUBLE, 0.10::DOUBLE, 0.002::DOUBLE)
) AS t(fragment_id, edge_id, like_weight_ratio, distal_length, pendant_length);

SELECT * FROM tree_resolve_placement('ref_tree', 'placements');

-- Full jplace workflow: extract tree and placements from jplace file
CREATE TABLE jplace_tree AS
SELECT * FROM read_jplace_newick('results.jplace');

CREATE TABLE jplace_placements AS
SELECT fragment AS fragment_id, edge_num::BIGINT AS edge_id,
       like_weight_ratio, distal_length, pendant_length
FROM read_jplace('results.jplace');

-- Resolve and inspect the result
SELECT name FROM tree_resolve_placement('jplace_tree', 'jplace_placements')
WHERE is_tip = true
ORDER BY name;

-- Write resolved tree to Newick format
COPY (
    SELECT node_index, name, branch_length, edge_id, parent_index
    FROM tree_resolve_placement('ref_tree', 'placements')
) TO 'resolved.nwk' (FORMAT NEWICK);

-- Use a view for filtered placements
CREATE VIEW confident_placements AS
SELECT * FROM jplace_placements WHERE like_weight_ratio > 0.8;

SELECT COUNT(*) FROM tree_resolve_placement('jplace_tree', 'confident_placements')
WHERE is_tip = true;
```

**Error conditions:**
- Tree table or placements table does not exist
- Tree table missing required `node_index` or `parent_index` columns
- Placements table missing required columns (`fragment_id`, `edge_id`, `like_weight_ratio`, `distal_length`, `pendant_length`)
- Placement references an `edge_id` not present in the tree
- `distal_length` exceeds the edge's branch length
- Negative `distal_length` or `pendant_length`

## `phylogeny_fasttree(table_name, [options])`

Build an approximately-maximum-likelihood phylogenetic tree from a multiple sequence alignment using [FastTree 2](http://www.microbesonline.org/fasttree/) (Price, Dehal, and Arkin 2010, Mol. Biol. Evol. 27:1641-1650). Returns one row per node (tips and internal nodes) with parent links suitable for joining back to the input or rendering to Newick via `COPY ... (FORMAT 'newick')`.

**Architecture:** FastTree is statically linked into a separate process (`gpl-boundary`) that miint launches per query. The daemon is GPL-licensed; miint stays BSD-clean by communicating with it over POSIX shared memory and Arrow IPC. See [`docs/internals/embedded-tools.md`](internals/embedded-tools.md) for the protocol details and process-lifecycle invariants.

**Requirements:**
- The `gpl-boundary` binary must be on `PATH`. Use [`phylogeny_fasttree_available()`](scalar-functions.md#phylogeny_fasttree_available) to check at runtime.

**Parameters:**

- `table_name` (VARCHAR): Name of a table or view containing the input alignment. Must have columns `name` (VARCHAR) and `sequence` (VARCHAR), one row per tip. Minimum 3 sequences (FastTree cannot build a meaningful tree below that). All sequences should be the same length (a true MSA — pad with `-` if needed; see `align_mafft`).

The remaining parameters mirror the FastTree CLI flags. Defaults are the daemon's FastTree defaults; omitted parameters are not emitted into the wire config.

*Sequence type and model:*

- `seq_type` (VARCHAR, default `'auto'`): One of `'auto'`, `'nucleotide'`, `'protein'`. Drives which substitution models are valid.
- `model` (VARCHAR, default `'auto'`): Substitution model. One of `'auto'`, `'jtt'`, `'lg'`, `'wag'` (protein), `'jc'`, `'gtr'` (nucleotide).
- `gtrrates` (DOUBLE[], optional): GTR rate parameters (6 entries: AC, AG, AT, CG, CT, GT, all ≥ 0). Only valid with `model='gtr'`.
- `gtrfreq` (DOUBLE[], optional): GTR base frequencies (4 entries: A, C, G, T, all ≥ 0). Only valid with `model='gtr'`.
- `gamma` (BOOLEAN, optional): Use a discretized gamma distribution for among-site rate variation.
- `cat` (BIGINT, optional, ≥ 1): Number of rate categories for CAT approximation.

*Tree search:*

- `nni` (BIGINT, optional, ≥ 0): Number of NNI rounds. `0` disables NNI.
- `spr` (BIGINT, optional, ≥ 0): Number of SPR rounds (subtree prune-regraft). `0` disables SPR.
- `mlnni` (BIGINT, optional, ≥ 0): Number of ML NNI rounds.
- `mlacc` (BIGINT, optional, ≥ 1): ML accuracy iterations per branch.
- `noml` (BOOLEAN, optional): Skip the ML refinement phase entirely. Mutually exclusive with `mlnni > 0`.
- `slow` (BOOLEAN, optional): Use slower / more accurate neighbor-joining (disables top hits).
- `bionj` (BOOLEAN, optional): Use BIONJ joins. Mutually exclusive with `nj`.
- `nj` (BOOLEAN, optional): Use vanilla NJ joins. Mutually exclusive with `bionj`.

*Top-hits heuristics:*

- `top` (BOOLEAN, optional): Force top-hits heuristics on.
- `notop` (BOOLEAN, optional): Force top-hits heuristics off.
- `topm` (DOUBLE, optional, > 0): Top-hits multiplier (typical: `1.0`-`2.0`). Only meaningful when top hits are enabled (rejected if `notop=true` or `slow=true`).

*Bootstrap / support:*

- `bootstrap` (BIGINT, optional, ≥ 0): Number of resamples for SH-test branch support. `0` disables. Mutually exclusive with `nosupport=true`.
- `nosupport` (BOOLEAN, optional): Suppress branch support output. Mutually exclusive with `bootstrap > 0`.

*Pseudocounts (low-confidence trees):*

- `pseudo` (BOOLEAN, optional): Enable pseudocounts.
- `pseudo_weight` (DOUBLE, optional, ≥ 0): Pseudocount weight. Requires `pseudo=true`.

*Misc:*

- `seed` (BIGINT, optional): RNG seed. Pin for reproducibility (combine with `threads=1`).
- `threads` (BIGINT, optional, ≥ 1): Thread count. **Output is non-deterministic across thread counts** — FastTree's parallel NJ/NNI/SPR sections produce floating-point variation. Force `threads=1` for reproducibility.
- `verbose` (BOOLEAN, optional): Forward FastTree's diagnostic output to the daemon log.
- `quote` (BOOLEAN, optional): Quote tip names in the rendered Newick (when going through `COPY FORMAT 'newick'`).

**Cross-parameter validation rules** (raised at SQL bind time, before any process spawn):

| Constraint | Error |
|---|---|
| `bootstrap > 0` + `nosupport=true` | "bootstrap=N cannot be combined with nosupport=true" |
| `pseudo_weight` without `pseudo=true` | "pseudo_weight requires pseudo=true" |
| `mlnni > 0` + `noml=true` | "mlnni=N cannot be combined with noml=true" |
| Protein model + `seq_type='nucleotide'` | "model='X' is a protein model and cannot be combined with seq_type='nucleotide'" |
| Nucleotide model + `seq_type='protein'` | "model='X' is a nucleotide model and cannot be combined with seq_type='protein'" |
| `gtrrates` / `gtrfreq` without `model='gtr'` | "gtrrates/gtrfreq are only valid with model='gtr'" |
| `bionj=true` + `nj=true` | "bionj=true and nj=true are mutually exclusive" |
| `top=true` + `notop=true` | "top=true and notop=true are mutually exclusive" |
| `slow=true` + `top=true` | "slow=true disables top hits; cannot also set top=true" |
| `topm` + (`notop=true` or `slow=true`) | "topm is only meaningful when top hits are enabled" |
| `fastest=true` | Currently rejected — the embedded library API is a strict subset of the CLI's `-fastest`. Will be supported once the upstream surface catches up. |

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `node_index` | BIGINT | 0-based row index, also the node identifier within this tree |
| `parent_index` | BIGINT (nullable) | `node_index` of this node's parent; `NULL` for the root |
| `edge_id` | BIGINT (nullable) | Edge identifier (matches the `{N}` placement notation, when used with `read_jplace`) |
| `branch_length` | DOUBLE (nullable) | Length of the edge from this node to its parent (NULL on the root) |
| `support` | DOUBLE (nullable) | SH-like branch support in `[0,1]` when `bootstrap > 0`; NULL otherwise |
| `is_tip` | BOOLEAN | True if this is a leaf (tip) node |
| `name` | VARCHAR (nullable) | Tip name (matches the input `name` column) for tips; NULL or empty for internal nodes |
| `n_children` | BIGINT | Number of direct children (0 for tips, 2 for binary internals, 3 for the unrooted root) |

The schema is `UNION ALL`-compatible with `read_newick` if you select `node_index, parent_index, name, branch_length, edge_id, is_tip` (i.e., drop `support` and `n_children`).

**Behavior:**

- One `gpl-boundary` daemon per query call (no cross-query connection caching yet).
- Sequences are read from `table_name` via a fresh DuckDB connection (avoids re-entrant `context.Query()` deadlock — see [`docs/internals/reading-tables-views.md`](internals/reading-tables-views.md)).
- The full alignment is materialized in shared memory and submitted as a single Arrow IPC batch (FastTree is not a streaming algorithm).
- The daemon is sent SIGTERM in the function's destructor; SIGKILL after a 30 × 10ms grace period if it doesn't exit cleanly.
- FastTree's text output uses `%.9f` for branch lengths; miint's intermediate Arrow representation is full IEEE-754 double precision. Round-trip parity to the CLI binary is verified in `test/sql/phylogeny_fasttree_parity.test`.

**Examples:**

```sql
-- Basic workflow: align with MAFFT, then build a tree
CREATE TABLE seqs AS SELECT read_id AS name, sequence1 AS sequence FROM read_fastx('16s.fasta');

-- Default tree (auto-detected seq_type, FastTree defaults for everything else)
SELECT * FROM phylogeny_fasttree('seqs');

-- Reproducible run: pin seed and force threads=1
SELECT * FROM phylogeny_fasttree('seqs', seed := 42, threads := 1);

-- Nucleotide GTR with branch support from 100 bootstraps
SELECT * FROM phylogeny_fasttree(
    'seqs',
    seq_type := 'nucleotide',
    model := 'gtr',
    bootstrap := 100,
    seed := 42
);

-- GTR with explicit rate parameters
SELECT * FROM phylogeny_fasttree(
    'seqs',
    seq_type := 'nucleotide',
    model := 'gtr',
    gtrrates := [1.0, 2.5, 1.0, 1.0, 2.5, 1.0],
    gtrfreq  := [0.25, 0.25, 0.25, 0.25]
);

-- Render the tree to a Newick file
COPY (
    SELECT node_index, parent_index, name, branch_length, edge_id
    FROM phylogeny_fasttree('seqs', seed := 42, threads := 1)
) TO 'tree.nwk' (FORMAT 'newick');

-- Tip-only summary
SELECT name, branch_length AS terminal_branch_length
FROM phylogeny_fasttree('seqs', seed := 42)
WHERE is_tip
ORDER BY branch_length DESC;

-- Pipeline: align → tree → resolve placements (round-trip via read_newick)
COPY (
    SELECT node_index, parent_index, name, branch_length, edge_id
    FROM phylogeny_fasttree('seqs', seed := 42, threads := 1)
) TO 'ref.nwk' (FORMAT 'newick');

CREATE TABLE ref_tree AS SELECT * FROM read_newick('ref.nwk');
SELECT * FROM tree_resolve_placement('ref_tree', 'placements');
```

**Error conditions:**

- `gpl-boundary` not on PATH (use `phylogeny_fasttree_available()` to gate)
- Input table does not exist or is missing `name` / `sequence` columns
- Fewer than 3 input sequences
- Any of the cross-parameter rules above
- Daemon process exits non-zero or returns a non-OK protocol response (exit message includes the daemon's stderr)

**Reproducibility:** With `threads=1` and a fixed `seed`, the tree is bit-deterministic across runs of the same `gpl-boundary` build. With `threads>1`, FastTree's parallel sections produce floating-point variation in branch lengths and (rarely) topology — pin `threads=1` for deterministic output.

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

Bowtie2 runs out of process via the `gpl-boundary` daemon (an Apache-licensed process-isolation host shipped separately from this extension). Sequences cross the boundary as Arrow IPC; bowtie2's command line is never user-controlled.

**Requirements:**
- The `gpl-boundary` daemon must be installed. Easiest path: `SELECT install_gpl_boundary();`. The miint extension itself does not link bowtie2.

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
- `threads` (INTEGER, default: 1): Number of threads for Bowtie2 alignment (daemon `nthreads`)
- `max_secondary` (INTEGER, default: 1): Maximum alignments to report per query (`-k`)
- `quiet` (BOOLEAN, default: true): Suppress Bowtie2 stderr output (alignment statistics)

The full bowtie2-align typed parameter set is also exposed (one-to-one with the daemon's `--describe`; integer ranges enforced at bind time):

| Parameter | Type | Notes |
|---|---|---|
| `seed` | INTEGER | Random seed for reproducibility |
| `trim5`, `trim3` | INTEGER | Trim N bases from each end of the read |
| `match_bonus` | INTEGER | `--ma` |
| `mismatch_penalty` | INTEGER | `--mp` |
| `n_penalty` | INTEGER | `--np` (ambiguous base penalty) |
| `read_gap_open`, `read_gap_extend` | INTEGER | `--rdg arg1,arg2` |
| `ref_gap_open`, `ref_gap_extend` | INTEGER | `--rfg arg1,arg2` |
| `score_min` | VARCHAR | Min-score function, e.g. `'L,-0.6,-0.6'` |
| `min_insert`, `max_insert` | INTEGER | Paired-end fragment-length bounds (`--minins`, `--maxins`) |
| `mate_orientation` | VARCHAR | One of `'fr'`, `'rf'`, `'ff'` |
| `no_mixed`, `no_discordant`, `dovetail`, `no_contain`, `no_overlap` | BOOLEAN | Paired-end behavior knobs |
| `nofw`, `norc` | BOOLEAN | Suppress forward / reverse-complement alignment |
| `seed_mismatches` | INTEGER | `-N`, must be 0 or 1 |
| `seed_length` | INTEGER | `-L`, must be 1–32 |
| `max_dp_failures` | INTEGER | `-D` |
| `max_seed_rounds` | INTEGER | `-R` |
| `report_all` | BOOLEAN | `-a`, report all alignments |
| `xeq` | BOOLEAN | Use `=`/`X` in CIGAR instead of `M` |
| `rg_id` | VARCHAR | `--rg-id` read group ID |
| `ignore_quals` | BOOLEAN | Treat all qualities as Phred 30 |
| `reorder` | BOOLEAN | Preserve input read order in output |

Unknown parameters are rejected at bind time by DuckDB's binder (e.g. `presset := 'fast'` → `Invalid named parameter "presset"`).

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
- Subject sequences are loaded at bind time and indexed via the daemon's `bowtie2-build` tool (index must fit in RAM)
- Query sequences are streamed across the gpl-boundary as Arrow IPC batches
- Supports both single-end and paired-end query sequences (paired-end auto-detected from the query table's `sequence2`/`qual2` columns)
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
```

**Error handling:**
- Error if query_table or subject_table does not exist
- Error if subject_table contains paired-end data (sequence2 not NULL)
- Error if tables lack required columns (read_id, sequence1)
- Error if the `gpl-boundary` daemon is not installed (see `SELECT install_gpl_boundary();`)

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

Align query sequences against multiple pre-built Bowtie2 index shards. Each shard is a directory containing a Bowtie2 index (with prefix `index`), and a mapping table specifies which reads should be aligned against which shard. This is the Bowtie2 counterpart to `align_minimap2_sharded`, designed for the same large-scale sharded alignment workflows. Bowtie2 runs out of process via the `gpl-boundary` daemon.

The sharded path emits only mapped reads (the daemon is invoked with `--no-unal`); use `align_bowtie2` directly if you need to inspect unaligned records.

**Requirements:**
- The `gpl-boundary` daemon must be installed (see `SELECT install_gpl_boundary();`).

**Parameters:**
- `query_table` (VARCHAR): Name of table or view containing query sequences. Must have `read_fastx`-compatible schema (read_id, sequence1, optional sequence2/qual1/qual2)
- `shard_directory` (VARCHAR, required): Path to directory containing shard subdirectories. Each shard's Bowtie2 index is expected at `<shard_directory>/<shard_name>/index` (i.e., files like `<shard_name>/index.1.bt2`, `<shard_name>/index.rev.1.bt2`, etc.)
- `read_to_shard` (VARCHAR, required): Name of table or view that maps reads to shards. Must have columns:
  - `read_id` (VARCHAR): Read identifier (must match read_id in query_table)
  - `shard_name` (VARCHAR): Name of the shard this read should be aligned against
- `preset` (VARCHAR, optional): Bowtie2 sensitivity preset ('very-fast', 'fast', 'sensitive', 'very-sensitive')
- `local` (BOOLEAN, default: false): Use local alignment mode instead of end-to-end
- `max_secondary` (INTEGER, default: 1): Maximum alignments to report per query (`-k`)
- `max_threads_per_shard` (INTEGER, default: 4, range 1–64): bowtie2's internal `nthreads` for each per-shard daemon worker
- `include_shard_name` (BOOLEAN, default: false): When true, append a `shard_name` column to the output
- `quiet` (BOOLEAN, default: true): Suppress Bowtie2 stderr output
- `threads` (INTEGER): Ignored in sharded mode. Use DuckDB's `SET threads=N` to control cross-shard parallelism and `max_threads_per_shard` for per-shard bowtie2 threading. A warning is printed at bind if `threads != 1` is passed directly to this function.

The full bowtie2-align typed parameter set listed under [`align_bowtie2`](#align_bowtie2query_table-subject_table-options) above is also available here (same names, same bind-time validation): `seed`, `trim5`, `trim3`, `match_bonus`, `mismatch_penalty`, `n_penalty`, `read_gap_open`, `read_gap_extend`, `ref_gap_open`, `ref_gap_extend`, `score_min`, `min_insert`, `max_insert`, `mate_orientation`, `no_mixed`, `no_discordant`, `dovetail`, `no_contain`, `no_overlap`, `nofw`, `norc`, `seed_mismatches`, `seed_length`, `max_dp_failures`, `max_seed_rounds`, `report_all`, `xeq`, `rg_id`, `ignore_quals`, `reorder`.

The `no_unal` knob is intentionally not exposed: the sharded path always emits only mapped reads (matches the pre-migration `FilterMappedOnly` contract). Use `align_bowtie2` directly if you need to inspect unaligned records.

**Output schema:**
Returns the same 21-column schema as `align_bowtie2` and `read_alignments`.

**Behavior:**
- At bind time, reads the `read_to_shard` table to discover shards; per-shard index existence is verified at InitGlobal (not Bind) so the planner doesn't pay filesystem-stat cost on wide shard sets.
- Cross-shard parallelism is driven by `SET threads`: each DuckDB worker thread owns its own gpl-boundary daemon and claims shards atomically. Per-shard bowtie2 internal threading is driven by `max_threads_per_shard`.
- For each shard, only the reads assigned to that shard (via the `read_to_shard` mapping) are streamed across the gpl-boundary as Arrow IPC
- A read can appear in multiple shards and will be aligned against each
- Unmapped reads (flag 0x4) are filtered out of results (daemon `--no-unal`)
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
```

**Error handling:**
- Error if `shard_directory` does not exist
- Error if any shard's Bowtie2 index files are missing (e.g., no `index.1.bt2` at `<shard_directory>/<shard_name>/index`)
- Error if `query_table` or `read_to_shard` table/view does not exist
- Error if `read_to_shard` table is missing `read_id` or `shard_name` columns
- Error if `read_to_shard` table contains NULL `shard_name` values
- Error if the `gpl-boundary` daemon is not installed (see `SELECT install_gpl_boundary();`)

**Performance notes:**
- Shards are processed sequentially; per-shard bowtie2 parallelism comes from `max_threads_per_shard`
- The `threads` parameter is ignored (a warning is emitted); use `max_threads_per_shard` instead
- Shards are sorted by read count (largest first) so the first shard's index stays hot in the daemon's worker pool

**Comparison of sharded vs non-sharded alignment functions:**
| Feature | `align_minimap2` / `align_bowtie2` | `align_minimap2_sharded` / `align_bowtie2_sharded` |
|---------|--------------------------------------|------------------------------------------------------|
| Index source | Build on-the-fly or single pre-built index | Multiple pre-built indexes (one per shard) |
| Read routing | All reads against one index | Reads routed to specific shards via mapping table |
| Parallelism | Single aligner thread(s) | One aligner per shard, concurrent |
| Use case | Single reference database | Sharded reference databases (e.g., from prior classification) |

## `align_mafft(table_name, [sample_id='col'])`

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

---

## `align_sortmerna(query_table, ref_paths=paths, [options])`

rRNA filtering / alignment against one or more rRNA reference databases using [SortMeRNA](https://github.com/biocore/sortmerna) (Kopylova et al. 2012, Bioinformatics 28:3211-3217). Embedded as a statically linked library (SortMeRNA 4.4.0 fork; LGPL-3.0-or-later).

Emits the standard 21-column SAM schema shared with `align_minimap2` / `align_bowtie2`. For a schema preserving SortMeRNA's native identity / coverage / e-value / edit-distance fields, use [`align_sortmerna_rrna`](#align_sortmerna_rrnaquery_table-ref_pathspaths-options) below.

**Parameters:**
- `query_table` (VARCHAR, positional): Name of a table or view with columns `read_id` (VARCHAR), `sequence1` (VARCHAR), and optionally `sequence2` (VARCHAR) when `paired := true`.
- `ref_paths` (VARCHAR[], required): List of FASTA paths for the reference database(s). The index is built once per query in-memory — re-using references across queries rebuilds the index each time.
- `num_threads` (INTEGER, default = DuckDB's thread count): Number of threads SortMeRNA's internal pool uses. The DuckDB function runs on a single DuckDB thread; parallelism is inside SortMeRNA.
- `match`, `mismatch`, `gap_open`, `gap_ext`, `score_N` (INTEGER): SW scoring. Defaults `2 / -3 / 5 / 2 / 0`.
- `evalue` (DOUBLE, default 1.0): E-value threshold.
- `seed_win_len` (UINTEGER, default 18): Seed window length.
- `num_alignments` (UINTEGER, default 1): Max alignments reported per read.
- `best` (BOOLEAN, default true): Keep only the best-scoring alignment per read.
- `paired` (BOOLEAN, default false): Paired-end mode. The query table must have a `sequence2` column; bind fails otherwise.
- `forward_only`, `reverse_only`, `full_search` (BOOLEAN): Strand-search controls.

**Output schema:** Same 21 columns as [`align_minimap2`](#align_minimap2query_table-subject_tablenull-index_pathnull-options). SortMeRNA-specific notes:
- `mapq` is always `255` — SortMeRNA does not compute mapping quality; 255 is the SAM convention for "unavailable".
- `tag_as` carries the raw Smith-Waterman score; `tag_nm` carries edit distance. Both are NULL for unaligned rows.
- `tag_xs`, `tag_ys`, `tag_xn`, `tag_xm`, `tag_xo`, `tag_xg`, `tag_yt`, `tag_md`, `tag_sa` are always NULL — SortMeRNA does not produce these.
- `stop_position` = `ref_end + 1` (1-based half-open), matching the convention used by `align_minimap2` and `align_bowtie2`.
- Paired-end: `flags & 0x2` (proper pair) is set when both mates aligned, regardless of reference. This is weaker than SAM's standard "concordant orientation within insert size" meaning — SortMeRNA is an rRNA classifier with no notion of insert size or orientation. When both mates aligned but to different references, `mate_reference` reports the actual partner reference name (rather than `=`), so cross-reference pairs remain distinguishable.

**Caveats:**
- **No minimum-score filter:** the embedded library returns every positive Smith-Waterman hit. The `sortmerna` CLI applies an internal minimum-score threshold that our streaming API path bypasses; to reproduce CLI output row-for-row, filter on `score` (`tag_as`) or on `e_value` in SQL. Be aware that e-values are not directly comparable between the library and the CLI (see below), so a library-side e-value threshold is not guaranteed to reproduce the exact CLI row set.
- **E-value divergence from the CLI:** the library computes per-query Karlin-Altschul e-values; the CLI sums / corrects across the full database. Identity, coverage, score, CIGAR, ref_start, ref_end, and edit_distance are bit-identical between the two.
- **Process-wide serialization:** SortMeRNA's `g_run_mutex` serializes calls process-wide. Concurrent `align_sortmerna` / `align_sortmerna_rrna` queries will block on each other at the library level. The DuckDB function deliberately runs on a single DuckDB thread (`MaxThreads() == 1`) to avoid queueing behind the mutex twice.

**Examples:**

```sql
CREATE TABLE reads AS SELECT read_id, sequence1 FROM read_fastx('metaT.fastq.gz');

-- Keep only reads that align to any SILVA reference with a reasonable e-value.
CREATE TABLE rrna_reads AS
  SELECT read_id, flags, reference, position, cigar, tag_as AS score
  FROM align_sortmerna('reads', ref_paths := ['silva-bac-16s.fasta', 'silva-arc-16s.fasta'])
  WHERE (flags & 0x4) = 0;  -- drop unaligned
```

---

## `align_sortmerna_rrna(query_table, ref_paths=paths, [options])`

Same aligner as `align_sortmerna` but emits SortMeRNA's native output schema, which preserves identity / coverage / e-value / edit-distance (SAM cannot carry all of these as first-class columns).

**Parameters:** identical to `align_sortmerna` above.

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `read_id` | VARCHAR | Query read identifier from the input table |
| `aligned` | INTEGER | `1` if the read aligned, `0` otherwise |
| `strand` | INTEGER | `1` for forward, `0` for reverse-complement. `0` for unaligned rows has no meaning. |
| `ref_name` | VARCHAR | Reference sequence ID (empty for unaligned rows) |
| `ref_start` | INTEGER | 1-based inclusive start position on reference (`0` for unaligned) |
| `ref_end` | INTEGER | 1-based inclusive end position on reference (`0` for unaligned) |
| `cigar` | VARCHAR | CIGAR string, SSW convention (uses `M/I/D/S`) |
| `score` | INTEGER | Raw Smith-Waterman score (not bitscore) |
| `e_value` | DOUBLE | Per-query Karlin-Altschul e-value |
| `identity` | DOUBLE | Percent identity (0..100) at full precision |
| `coverage` | DOUBLE | Query coverage (0..100) at full precision |
| `edit_distance` | INTEGER | Number of mismatches + indels |
| `segment_idx` | INTEGER | `0` for single-end or the forward mate; `1` for the reverse mate in paired-end output |

Paired-end mode produces two rows per input row with `segment_idx` 0 (forward) and 1 (reverse), even when one or both mates failed to align.

**Example:**

```sql
-- Assign reads to GG references, filter, aggregate.
SELECT ref_name, COUNT(*) AS hits
FROM align_sortmerna_rrna('reads',
       ref_paths := ['gg_13_8.fasta'])
WHERE aligned = 1 AND e_value <= 1e-5 AND coverage >= 80.0
GROUP BY ref_name ORDER BY hits DESC;
```

---

## `detect_chimera_uchime(query_table, db='refs_table', [sample_id='col'], [options])`

Reference-based chimera detection using the UCHIME algorithm (Edgar et al. 2011, Bioinformatics 27:2194-2200), powered by the [vsearch](https://github.com/torognes/vsearch) library (Rognes et al. 2016, PeerJ 4:e2584). Detects chimeric sequences by comparing queries against a trusted chimera-free reference database.

**Parameters:**
- `query_table` (VARCHAR): Name of a table or view containing query sequences. Must have `read_id` (VARCHAR) and `sequence1` (VARCHAR) columns.
- `db` (VARCHAR, required): Name of a table or view containing reference sequences. Same schema requirements as `query_table`.
- `sample_id` (VARCHAR, optional): Name of a column in `query_table` to partition by. When provided, queries are scored per-sample against the (shared, load-once) reference database, and the sample column is prepended to the output. Execution is serialized (the vsearch wrapper is not thread-safe across concurrent calls).
- `minh` (DOUBLE, default 0.28): Minimum h-score to flag as chimeric. Range [0, 1].
- `xn` (DOUBLE, default 8.0): Weight of "no" votes in h-score computation. Must be >= 1.0.
- `dn` (DOUBLE, default 1.4): Pseudo-count prior on "no" votes. Must be >= 0.
- `mindiv` (DOUBLE, default 0.8): Minimum divergence (percentage points) from closest parent. Must be >= 0.
- `mindiffs` (INTEGER, default 3): Minimum number of diffs in each segment (left and right). Must be >= 1.
- `threads` (INTEGER, optional): Number of threads vsearch uses for the internal `chimera_detect_batch` parallel scan. Defaults to DuckDB's configured thread count (`SET threads=N`) at bind time; pass an explicit value to override. Must be 1–1024 (matching vsearch's CLI ceiling).

**Output schema (18 columns, compatible with vsearch `--uchimeout`):**

| Column | Type | Description |
|--------|------|-------------|
| `score` | DOUBLE | Chimera h-score (higher = more likely chimeric) |
| `read_id` | VARCHAR | Query sequence identifier |
| `parent_a_id` | VARCHAR | Parent A identifier (NULL if non-chimeric) |
| `parent_b_id` | VARCHAR | Parent B identifier (NULL if non-chimeric) |
| `closest_parent_id` | VARCHAR | Closest parent to query (NULL if non-chimeric) |
| `id_query_model` | DOUBLE | Query-to-chimeric-model identity % |
| `id_query_a` | DOUBLE | Query-to-parent-A identity % |
| `id_query_b` | DOUBLE | Query-to-parent-B identity % |
| `id_a_b` | DOUBLE | Parent-A-to-parent-B identity % |
| `id_query_top` | DOUBLE | Query-to-closest-parent identity % |
| `left_yes` | INTEGER | Left segment yes votes |
| `left_no` | INTEGER | Left segment no votes |
| `left_abstain` | INTEGER | Left segment abstain votes |
| `right_yes` | INTEGER | Right segment yes votes |
| `right_no` | INTEGER | Right segment no votes |
| `right_abstain` | INTEGER | Right segment abstain votes |
| `divergence` | DOUBLE | Model divergence (id_query_model - id_query_top) |
| `flag` | VARCHAR | Classification: `Y` (chimera), `N` (non-chimera), `?` (borderline) |

**Example:**
```sql
-- Load sequences from FASTA files into tables
CREATE TABLE refs AS SELECT read_id, sequence1 FROM read_fastx('gold.fasta');
CREATE TABLE queries AS SELECT read_id, sequence1 FROM read_fastx('amplicons.fasta');

-- Detect chimeras
SELECT * FROM detect_chimera_uchime('queries', db:='refs');

-- Filter chimeric sequences
CREATE TABLE clean_seqs AS
SELECT q.* FROM queries q
JOIN detect_chimera_uchime('queries', db:='refs') u ON q.read_id = u.read_id
WHERE u.flag = 'N';

-- Count chimeras
SELECT flag, count(*) FROM detect_chimera_uchime('queries', db:='refs') GROUP BY flag;
```

**Behavior:**
- One output row per query sequence (all queries are reported, not just chimeras)
- Non-chimeric sequences (flag=`N`) have NULL for parent and identity columns, and 0 for vote columns
- `id_a_b` (parent A vs parent B identity) is only computed for chimeric (`Y`) and borderline (`?`) results. Non-chimeric results report `id_a_b=0.0`. This avoids an extra pairwise alignment per query that is not needed for classification. Note: vsearch computes `id_a_b` for all queries unconditionally.
- Multi-threaded: vsearch's internal `chimera_detect_batch` parallelizes across the `threads` parameter (defaults to DuckDB's configured thread count; override per-call with `threads:=N`)
- The reference database is fully materialized in memory at init time
- Tables and views are both supported for query and reference inputs

**Error handling:**
- Error if `db` parameter is missing
- Error if query or reference table does not exist or lacks `read_id`/`sequence1` columns
- Error if reference table is empty
- Error if query table contains NULL `read_id` values
- Error if scoring parameters are out of valid range

**Algorithm:**
1. For each query, partition into 4 chunks and search the reference DB using an 8-mer index for candidate parents (up to 16)
2. Align query to each candidate using WFA2 global alignment
3. Select the 2 best parents via smoothed identity (32bp sliding window)
4. Build a 3-way star alignment and classify each column (match-A, match-B, no-vote, abstain)
5. Sweep all breakpoints left-to-right, computing h-score at each position
6. Classify based on h-score, divergence, and minimum diff thresholds

---

## `detect_chimera_uchime_denovo(input_table, [sample_id='col'], [options])`

De novo chimera detection using the UCHIME algorithm, powered by the [vsearch](https://github.com/torognes/vsearch) library (Rognes et al. 2016, PeerJ 4:e2584). Detects chimeric sequences without a reference database by using abundance information: more abundant sequences are assumed to be non-chimeric and serve as parents for less abundant sequences.

**Parameters:**
- `input_table` (VARCHAR): Name of a table or view containing sequences with abundance. By default must have `read_id` (VARCHAR), `sequence1` (VARCHAR), and `size` (integer type) columns; use `id_col`/`sequence_col`/`count_col` to override.
- `sample_id` (VARCHAR, optional): Name of a column in `input_table` to partition by. Each sample gets its own k-mer index and bootstrap; a read_id that appears in multiple samples is therefore scored independently. The sample column is prepended to the output. Execution is serialized per the vsearch wrapper's thread-safety constraints.
- `id_col` (VARCHAR, default `'read_id'`): Name of the read identifier column in `input_table`.
- `sequence_col` (VARCHAR, default `'sequence1'`): Name of the sequence column.
- `count_col` (VARCHAR, default `'size'`): Name of the per-sequence count column. Set to `'abundance'` to chain `deblur(...)` directly into this function.
- `abskew` (DOUBLE, default 2.0): Abundance skew. Candidate parents must have abundance >= abskew * query abundance. Must be >= 1.0.
- `minh`, `xn`, `dn`, `mindiv`, `mindiffs`: Same as `detect_chimera_uchime`. (No `threads` parameter — de novo detection is sequential by construction; vsearch is run with `opt_threads=1`.)

**Output schema:** Same 18 columns as `detect_chimera_uchime`.

**Example:**
```sql
-- Load sequences and compute abundance
CREATE TABLE seqs AS
SELECT read_id, sequence1, count(*) AS size
FROM read_fastx('amplicons.fasta')
GROUP BY read_id, sequence1;

-- De novo chimera detection
SELECT * FROM detect_chimera_uchime_denovo('seqs');

-- Filter out chimeras
SELECT s.* FROM seqs s
JOIN detect_chimera_uchime_denovo('seqs') u ON s.read_id = u.read_id
WHERE u.flag != 'Y';
```

**Behavior:**
- Sequences are processed in decreasing abundance order (highest first)
- The first two sequences (most abundant) are unconditionally treated as non-chimeric to seed the reference database
- Non-chimeric and borderline sequences are added to the reference DB incrementally
- Chimeric sequences are NOT added to the reference DB
- Single-threaded (inherently sequential — each result depends on previous classifications)
- One output row per input sequence

**Error handling:**
- Error if table does not exist or lacks the resolved id/sequence/count columns
- Error if the count column is not an integer type
- Error if table is empty
- Error if scoring parameters are out of valid range
- Error if any of `id_col`/`sequence_col`/`count_col` is the empty string

---

## `search_sequences_vsearch(query_table, db='ref_table', id=threshold, [options])`

Global pairwise sequence search, powered by the [vsearch](https://github.com/torognes/vsearch) library (Rognes et al. 2016, PeerJ 4:e2584). Finds the best matching sequences in a reference database for each query sequence using SIMD-optimized Needleman-Wunsch alignment with k-mer candidate filtering.

**Parameters:**
- `query_table` (VARCHAR): Name of a table or view containing query sequences. Must have `read_id` (VARCHAR) and `sequence1` (VARCHAR) columns.
- `db` (VARCHAR, required): Name of a table or view containing reference sequences. Same schema requirements as `query_table`.
- `id` (DOUBLE, required): Minimum identity threshold (0.0-1.0). No silent default — must be specified explicitly.
- `maxaccepts` (INTEGER, default 1): Maximum number of accepted hits per query. Must be >= 1.
- `maxrejects` (INTEGER, default 32): Maximum rejected targets before stopping search. Must be >= 1.
- `threads` (INTEGER, optional): Number of threads vsearch uses for its internal `search_batch` parallel scan. Defaults to DuckDB's configured thread count (`SET threads=N`) at bind time; pass an explicit value to override. Must be 1–1024 (matching vsearch's CLI ceiling).

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `read_id` | VARCHAR | Query sequence identifier |
| `target_id` | VARCHAR | Reference sequence identifier |
| `identity` | DOUBLE | Percent identity (0-100) |
| `matches` | INTEGER | Number of matching columns |
| `mismatches` | INTEGER | Number of mismatching columns |
| `gaps` | INTEGER | Number of gap columns |
| `alignment_length` | INTEGER | Total alignment length |
| `query_length` | INTEGER | Query sequence length |
| `target_length` | INTEGER | Target sequence length |
| `accepted` | BOOLEAN | True if hit passes identity threshold |

**Example:**
```sql
CREATE TABLE refs AS SELECT read_id, sequence1 FROM read_fastx('database.fasta');
CREATE TABLE queries AS SELECT read_id, sequence1 FROM read_fastx('queries.fasta');

-- Search at 97% identity
SELECT * FROM search_sequences_vsearch('queries', db:='refs', id:=0.97);

-- Top 3 hits per query at 90% identity
SELECT * FROM search_sequences_vsearch('queries', db:='refs', id:=0.90, maxaccepts:=3);

-- Count queries with hits
SELECT count(DISTINCT read_id) FROM search_sequences_vsearch('queries', db:='refs', id:=0.97);
```

**Behavior:**
- Each query produces 0 to `maxaccepts` output rows
- Results include both accepted (above threshold) and weak (near-miss) hits; use the `accepted` column to filter
- Plus-strand only (no reverse complement search)
- Multi-threaded: vsearch's internal `search_batch` parallelizes across the `threads` parameter (defaults to DuckDB's configured thread count; override per-call with `threads:=N`)
- Reference database is fully materialized in memory at init time
- RNA sequences (U) are automatically converted to DNA (T)

**Error handling:**
- Error if `db` or `id` parameter is missing
- Error if `id` is not between 0.0 and 1.0
- Error if query or reference table does not exist or lacks required columns
- Error if reference table is empty

---

## `cluster_sequences_vsearch(input_table, id=threshold, [options])`

Greedy sequence clustering, powered by the [vsearch](https://github.com/torognes/vsearch) library (Rognes et al. 2016, PeerJ 4:e2584). Clusters sequences by iterating in input order: each sequence is compared against existing centroids, and either joins the best matching cluster (if above the identity threshold) or becomes a new centroid.

**Parameters:**
- `input_table` (VARCHAR): Name of a table or view containing sequences. Must have `read_id` (VARCHAR) and `sequence1` (VARCHAR) columns.
- `id` (DOUBLE, required): Minimum identity threshold (0.0-1.0). No silent default — must be specified explicitly.
- `strand` (VARCHAR, default `'plus'`): `'plus'` for plus-strand only, `'both'` to also search reverse complements.
- `threads` (INTEGER, optional): Number of threads vsearch uses for its internal `cluster_assign_batch` parallel scan. Defaults to DuckDB's configured thread count (`SET threads=N`) at bind time; pass an explicit value to override. Must be 1–1024 (matching vsearch's CLI ceiling).

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `read_id` | VARCHAR | Input sequence identifier |
| `is_centroid` | BOOLEAN | True if this sequence started a new cluster |
| `cluster_id` | INTEGER | Cluster number (0-based) |
| `centroid_id` | VARCHAR | Identifier of the cluster's centroid |
| `identity` | DOUBLE | Percent identity to centroid (100.0 if centroid) |
| `cigar` | VARCHAR | CIGAR alignment to centroid (empty if centroid) |
| `cigar_truncated` | BOOLEAN | True if CIGAR was truncated (>4096 chars) |

**Sort order is the caller's responsibility.** The function clusters sequences in the order they appear in the input table. For `cluster_fast`-equivalent behavior (longest first), sort by length descending. For `cluster_size`-equivalent behavior (most abundant first), sort by abundance descending:

```sql
-- Load and sort by length descending (like vsearch --cluster_fast)
CREATE TABLE sorted_seqs AS
  SELECT * FROM read_fastx('sequences.fasta')
  ORDER BY length(sequence1) DESC;
SELECT * FROM cluster_sequences_vsearch('sorted_seqs', id:=0.97);

-- Cluster by abundance (like vsearch --cluster_size)
CREATE TABLE by_abundance AS
  SELECT read_id, sequence1, count(*) AS size
  FROM read_fastx('amplicons.fasta')
  GROUP BY read_id, sequence1
  ORDER BY size DESC;
SELECT * FROM cluster_sequences_vsearch('by_abundance', id:=0.97);

-- Count clusters
SELECT count(*) FROM cluster_sequences_vsearch('sorted_seqs', id:=0.97) WHERE is_centroid;

-- Get cluster sizes
SELECT centroid_id, count(*) AS size
FROM cluster_sequences_vsearch('sorted_seqs', id:=0.97)
GROUP BY centroid_id ORDER BY size DESC;
```

**Behavior:**
- One output row per input sequence
- Single-threaded (inherently sequential — each centroid must be indexed before the next sequence is processed)
- RNA sequences (U) are automatically converted to DNA (T)
- All results are materialized before returning (vsearch session mutex held for the duration)

**Error handling:**
- Error if `id` parameter is missing
- Error if `id` is not between 0.0 and 1.0
- Error if `strand` is not `'plus'` or `'both'`
- Error if table does not exist or lacks required columns
- Error if table is empty or contains NULL read_ids, NULL sequences, or empty sequences
- Error if any read_id exceeds 1023 characters

## `deblur(input_table, [sample_id='col'], [options])`

Deblur amplicon sequence denoising (Amir et al. 2017, mSystems 2:e00191-16). A greedy deconvolution algorithm that removes sequencing errors from amplicon data by iteratively subtracting expected error-derived reads from less-abundant sequences. Sequences whose corrected abundance rounds to zero are removed as errors; the remainder are denoised "sub-OTUs" (sOTUs).

Designed as a composable SQL building block. Dereplication is native SQL (`GROUP BY`), alignment is `align_mafft()`, and `deblur()` does the core denoising. See the full workflow example below.

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

**Full workflow example:**

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

**Minimal example (pre-aligned sequences of equal length):**

```sql
CREATE TABLE seqs(read_id VARCHAR, sequence1 VARCHAR, abundance BIGINT);
INSERT INTO seqs VALUES
  ('true_seq', 'ACGTACGTACGTACGT', 1000),
  ('error_seq', 'ACGTACGTACGTACGA', 3);

SELECT * FROM deblur('seqs');
-- Only true_seq survives (error_seq is explained by sequencing error)
```

**Tuning for modern platforms:**

```sql
-- NovaSeq with stitched 250nt reads
SELECT * FROM deblur('aligned_seqs', mean_error := 0.002);

-- Custom error profile from mock community calibration
SELECT * FROM deblur('aligned_seqs', error_profile := [1, 0.04, 0.01, 0.005]);
```

**Behavior:**
- Sequences are uppercased internally (matching the Python deblur reference implementation)
- Output is ordered by corrected abundance descending
- Single-threaded (inherently sequential — each sequence's correction depends on all previous corrections)
- Empty input tables return zero rows (not an error)
- Uses banker's rounding (round half to even) for the final abundance, matching Python 3's `round()`

**Error handling:**
- Error if table does not exist or lacks the resolved id/sequence/count columns
- Error if any of `id_col`/`sequence_col`/`count_col` is the empty string
- Error if `abundance` column is not an integer type
- Error if `mean_error` is not in the open interval (0, 1)
- Error if `indel_max` is negative
- Error if `error_profile` contains negative values or is empty
- Error if sequences have different aligned lengths or different unaligned lengths

---

## `unifrac_pcoa(observations, tree, [options])`

The three `unifrac_*` functions share enough context (input schemas, accepted variant strings, subsampling semantics) that they have their own reference at **[`docs/unifrac.md`](unifrac.md)**. The entries here are short signatures for discoverability; the linked doc is authoritative.

Compute UniFrac distances and reduce to PCoA coordinates via randomized FSVD (fp32). Operates on a long-form feature-table relation (matching the [`read_biom`](#read_biomfilename-include_filepathfalse) schema) and a tree relation (matching the [`read_newick`](#read_newickfilename-include_filepathfalse) schema).

```sql
SELECT * FROM unifrac_pcoa('observations', 'tree',
    variant := 'weighted_normalized', n_dims := 3, seed := 42);
```

Output: `(iteration INTEGER, sample_id VARCHAR, axis INTEGER, coordinate DOUBLE, eigenvalue DOUBLE, proportion_explained DOUBLE)`. The `eigenvalue` and `proportion_explained` columns repeat across samples within each `(iteration, axis)` pair — extract per-axis summaries with `SELECT DISTINCT`.

Full parameter reference, replication semantics, fp32 reproducibility note, and worked examples: see **[`docs/unifrac.md`](unifrac.md#unifrac_pcoa)**.

---

## `unifrac_permanova(observations, tree, metadata, [options])`

PERMANOVA pseudo-F + p-value on a UniFrac distance matrix, against a wide-form metadata relation (one row per sample, one column per variable; `sample_id` column is required).

```sql
SELECT * FROM unifrac_permanova('observations', 'tree', 'metadata',
    variables := ['body_site', 'treatment'], n_permutations := 999, seed := 42);
```

Output: `(iteration INTEGER, variable VARCHAR, n_groups INTEGER, f_stat DOUBLE, p_value DOUBLE, n_permutations INTEGER)`. Identical groupings under the same seed produce identical `f_stat` (Rule-7 invariant — value-driven, order-stable factorization).

Full parameter reference, variant strings, and worked examples: see **[`docs/unifrac.md`](unifrac.md#unifrac_permanova)**.

---

## `unifrac_faith_pd(observations, tree, [options])`

Faith's phylogenetic diversity per sample (sum of branch lengths on the spanning subtree). Optionally rarefies via `subsample_depth` and `n_subsamples` to produce a multi-iteration distribution.

```sql
SELECT * FROM unifrac_faith_pd('observations', 'tree',
    subsample_depth := 3, n_subsamples := 100, seed := 42);
```

Output: `(iteration INTEGER, sample_id VARCHAR, faith_pd DOUBLE)`. Deterministic — same seed produces byte-identical output (no fp32 tolerance needed).

Full parameter reference and worked examples: see **[`docs/unifrac.md`](unifrac.md#unifrac_faith_pd)**.

---

## `miint_warnings()`

Returns miint's operational warnings as a queryable table. Populated by user-facing warning sites — skipped accessions in `read_ena_sequences`, the SFF `max_sequences` caveat, the `threads` parameter being ignored in `align_bowtie2_sharded`, mid-stream download failures, etc. Every entry is also printed to stderr (today's behavior), so interactive users see no change; pipeline and `COPY TO` users now have a way to inspect warnings after the fact.

**Returns:**
| Column | Type | Description |
|---|---|---|
| `timestamp` | `TIMESTAMP WITH TIME ZONE` | When the warning was emitted |
| `message` | `VARCHAR` | Human-readable warning text |

**Example:**
```sql
-- Scan a study, then see what got skipped
SELECT COUNT(*) FROM read_ena_sequences('PRJEB1234');
SELECT timestamp, message FROM miint_warnings();
```

**Behavior:**
- The log is process-scoped and in-memory by default; entries accumulate across queries within a single DuckDB session.
- Warnings are captured regardless of the `enable_logging` / `logging_level` settings — miint writes directly to DuckDB's global log sink so `miint_warnings()` works without any setup.
- Under the hood, entries have type `'MiintWarning'` and live alongside DuckDB's own logs in `duckdb_logs()`; the macro just filters on type.

**Interactions with DuckDB logging settings:**
- `SET logging_storage='stderr'` — entries go straight to stderr instead of the in-memory sink, so `miint_warnings()` stops returning rows in that mode (the stderr output is the workaround). Default storage is in-memory; leave it alone unless you have reason to change it.
- `SET warnings_as_errors=true` — miint suppresses the log-storage write for that query so a skip-warning doesn't abort the retry or the partial-data path it was designed to allow. The stderr message still fires, so nothing is silently dropped, but `miint_warnings()` will not see the row for that call.
