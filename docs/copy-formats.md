# COPY Formats

DuckDB miint provides custom COPY formats for writing bioinformatics file formats.

## Table of Contents

- [FORMAT FASTQ](#copy-to-format-fastq) - Write FASTQ sequence files
- [FORMAT FASTA](#copy-to-format-fasta) - Write FASTA sequence files
- [FORMAT SAM / FORMAT BAM](#copy-to-format-sam-and-copy-to-format-bam) - Write SAM/BAM alignment files
- [FORMAT BIOM](#copy-to-format-biom) - Write BIOM observation matrix files
- [FORMAT NEWICK](#copy-to-format-newick) - Write Newick phylogenetic trees

## `COPY ... TO '...' (FORMAT FASTQ)`

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

## `COPY ... TO '...' (FORMAT FASTA)`

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

## `COPY ... TO '...' (FORMAT SAM)` and `COPY ... TO '...' (FORMAT BAM)`

Write query results to SAM or BAM format files. Requires all mandatory SAM columns from `read_alignments` output.

**Required columns:**
- `read_id` (VARCHAR or BIGINT): Query template name
- `flags` (USMALLINT): Bitwise flags
- `reference` (VARCHAR or BIGINT): Reference sequence name
- `position` (BIGINT): 1-based leftmost mapping position
- `mapq` (UTINYINT): Mapping quality
- `cigar` (VARCHAR): CIGAR string
- `mate_reference` (VARCHAR or BIGINT): Reference name of mate/next read
- `mate_position` (BIGINT): Position of mate/next read
- `template_length` (BIGINT): Observed template length

**Identifier-column types (`read_id`, `reference`, `mate_reference`):**
- Each column is independently `VARCHAR` or `BIGINT`. The three need not match. Other numeric types are rejected at bind time with the message `Column '<name>' must be VARCHAR or BIGINT`.
- `BIGINT` values are stringified to decimal on write — the on-disk SAM/BAM record always carries text. NULL `BIGINT` values are written as the SAM `*` sentinel.
- Round-trip through `read_alignments` / `read_sam` returns these columns as `VARCHAR` regardless of how they were written; the BIGINT type is not preserved on disk.
- HTSlib normalises `RNEXT` to `=` when it equals `RNAME` at the reference-tid level. This applies to both VARCHAR and BIGINT writes — a `BIGINT` mate_reference that matches `reference` still appears as `=` when the file is read back.

**Optional columns:**
- `tag_as`, `tag_xs`, `tag_ys`, `tag_xn`, `tag_xm`, `tag_xo`, `tag_xg`, `tag_nm` (BIGINT): Optional integer tags
- `tag_yt`, `tag_md`, `tag_sa` (VARCHAR): Optional string tags

**Parameters:**
- `INCLUDE_HEADER` (default: true): Include header with reference sequences
  - **Note:** BAM format requires `INCLUDE_HEADER=true` (headers are mandatory in BAM files)
- `REFERENCE_LENGTHS` (VARCHAR, required if INCLUDE_HEADER=true): Table or view name containing reference sequences. Must have at least 2 columns: first column = reference name (VARCHAR), second column = reference length (INTEGER/BIGINT). Column names don't matter. Views are fully supported and can include computed columns.
- `SEQUENCE_DATA` (VARCHAR, optional): Table or view name containing original read sequences from `read_fastx`. When provided, writes actual SEQ and QUAL fields into the output instead of `*`. See [Sequence Data](#sequence-data) below.
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

-- Write BAM with sequence data (SEQ and QUAL populated from FASTQ)
CREATE TABLE sequences AS SELECT * FROM read_fastx('reads_R1.fastq', 'reads_R2.fastq');

COPY (SELECT * FROM read_alignments('input.bam'))
TO 'with_seq.bam' (FORMAT BAM, REFERENCE_LENGTHS 'ref_table', SEQUENCE_DATA 'sequences');
```

### Sequence Data

The `SEQUENCE_DATA` parameter allows writing actual SEQ and QUAL fields into SAM/BAM output by looking up original read sequences from a table or view. Without this parameter, SEQ and QUAL are written as `*`.

This design avoids duplicating immutable sequence data across alignment records. The original FASTQ sequences are stored once and joined at write time.

**Required columns in the SEQUENCE_DATA table:**
- `read_id` (VARCHAR): Read identifier (must match alignment read_id values)
- `sequence1` (VARCHAR): Forward-strand DNA sequence for read 1
- `qual1` (LIST(UTINYINT)): Phred quality scores for read 1 (raw values 0-93, no ASCII offset)

**Optional columns (for paired-end):**
- `sequence2` (VARCHAR): Forward-strand DNA sequence for read 2
- `qual2` (LIST(UTINYINT)): Phred quality scores for read 2

These columns match the output schema of `read_fastx`, so the typical workflow is:

```sql
-- Load sequences
CREATE TABLE sequences AS SELECT * FROM read_fastx('R1.fastq', 'R2.fastq');

-- Load alignments
CREATE TABLE alignments AS SELECT * FROM read_alignments('aligned.bam');

-- Write with sequences
COPY (SELECT * FROM alignments)
TO 'output.bam' (FORMAT BAM, REFERENCE_LENGTHS 'ref_table', SEQUENCE_DATA 'sequences');
```

**Behavior:**
- **Paired-end reads**: Uses SAM flag bit 0x80 (second in pair) to select `sequence2`/`qual2` vs `sequence1`/`qual1`
- **Reverse strand**: When flag bit 0x10 is set, the sequence is reverse-complemented and quality scores are reversed
- **Hard clipping**: CIGAR H operations trim the original sequence accordingly (leading H trims from start, trailing H from end)
- **CIGAR validation**: After hard clipping, sequence length must match the sum of query-consuming CIGAR operations (M, I, S, =, X)
- **NULL quality**: If `qual1`/`qual2` is NULL, sequence is written but quality is marked as missing
- **Missing reads**: An error is thrown if any alignment's `read_id` is not found in the SEQUENCE_DATA table

**Notes:**
- The SEQUENCE_DATA table is loaded entirely into memory at initialization. For very large datasets, ensure sufficient RAM is available.
- Reference lengths must be provided explicitly when writing headers - they cannot be inferred from the data
- All optional tags present in the input are preserved in the output
- BAM files always require headers (binary format specification)

## `COPY ... TO '...' (FORMAT BIOM)`

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
    SELECT * FROM woltka_ogu('my_alignments', 'read_id', sample_id := 'sample_id')
) TO 'ogu_table.biom' (FORMAT BIOM);

-- With HDF5 gzip compression (recommended)
COPY (
    SELECT * FROM woltka_ogu('my_alignments', 'read_id', sample_id := 'sample_id')
) TO 'ogu_table.biom' (FORMAT BIOM, COMPRESSION 'gzip');

-- With custom metadata
COPY (
    SELECT * FROM woltka_ogu('my_alignments', 'read_id', sample_id := 'sample_id')
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
    SELECT * FROM woltka_ogu('high_qual', 'read_id', sample_id := 'sample_id')
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
    SELECT * FROM woltka_ogu('all_samples', 'read_id', sample_id := 'sample_id')
) TO 'all_samples.biom' (FORMAT BIOM, COMPRESSION 'gzip');

-- Export single-sample results (add sample_id column)
COPY (
    SELECT feature_id, 'MySample' AS sample_id, value
    FROM woltka_ogu('my_alignments', 'read_id')
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

## `COPY ... TO '...' (FORMAT NEWICK)`

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
- "distal_length is negative" - distal_length must be >= 0
- "pendant_length is negative" - pendant_length must be >= 0
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
