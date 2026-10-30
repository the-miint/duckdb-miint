# Reading files

MIINT supports reading FASTA, FASTQ, SAM, BAM, PAF, GFA, SFF, BIOM, mzML, mzXML, GFF, jplace and Newick, and the output directory of the myloasm assembler.  

A reader expands each path it is given as DuckDB's own readers do: a `file:` URI names the local path it spells (`file:///p`, `file://localhost/p` and `file:/p` name `/p`; on Windows, `file:///C:/p` names `C:/p`), and a leading `~` names the directory of DuckDB's `home_directory` setting (by default your home directory). This holds for a single path, for each path in a list (a list's paths are literal paths, never globs), and for the files a glob matches. Index and database paths given to aligners and classifiers are passed to their libraries as given, so give those as plain paths.

## Table of contents

- [FASTA / FASTQ](#fasta-and-fastq) - FASTA and FASTQ sequence files.
- [SAM / BAM](#sam-and-bam) - SAM and BAM alignment files.
  - [`read_alignment_header(path)`](#read_alignment_headerpath) - A header's `@SQ` lines as `(tid, reference, length)`.
- [PAF](#paf) - PAF pairwise alignment files, as written by minimap2 and other long-read tools.
- [GFA](#gfa) - GFA 1 assembly graphs.
  - [`read_gfa_segments(filename)`](#read_gfa_segmentsfilename) - A graph's segments as typed columns.
  - [`read_gfa_links(filename)`](#read_gfa_linksfilename) - A graph's links as typed columns.
- [myloasm](#myloasm) - A myloasm assembly's output directory as tables.
  - [`read_myloasm_contigs(dir)`](#read_myloasm_contigsdir) - Contigs, with dereplication status, lengths and depths.
  - [`read_myloasm_layout(dir)`](#read_myloasm_layoutdir) - The read slices each contig is built from.
  - [`read_myloasm_mappings(dir)`](#read_myloasm_mappingsdir) - Every read's mappings to the contigs.
  - [`read_myloasm_links(dir)`](#read_myloasm_linksdir) - The contig graph's edges.
  - [`read_myloasm_splits(dir)`](#read_myloasm_splitsdir) - Where chimeric reads were split.
  - [`read_myloasm_run(dir)`](#read_myloasm_rundir) - Each invocation's version, command, k and c.
- [SFF](#sff) - SFF sequence files.
- [mzML / mzXML](#mzml-and-mzxml) - mzML and mzXML files.
- [BIOM](#biom) - BIOM-Format v2.1 files.
- [jplace](#jplace) - jplace formatted files.
- [Newick](#newick) - Phylogenies in Newick.
- [GFFv3](#gff) - GFF formatted annotation data.
- [Coordinate conventions](#coordinate-conventions) - How positions on references and on reads are reported, and how each format is converted.

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
- `allow_empty` (BOOLEAN, optional, default false): Read an empty file as no records instead of failing, for tools that write a file whether or not they have records for it. A pair with exactly one empty mate still fails: the other mate's records would have no mates.
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
- Detects FASTA (`>` headers) or FASTQ (`@` headers) from the file's content, whatever its extension. Both files of a pair must be the same format: a FASTA file paired with a FASTQ file fails the query
- Reads gzip-compressed files. A local file is recognized as gzip by its content, whatever its name; a remote file (`s3://`, `https://`) only by a `.gz` extension
- Supports stdin input using `-` or `/dev/stdin` (single file only, no paired-end)
- Quality scores converted to integers using specified offset (Phred+33 or Phred+64)
- FASTQ records, those whose header starts with `@`, must be four lines each: header, sequence, `+`, quality, with the quality as long as the sequence. Blank lines, including lines of nothing but whitespace, may come before the first record and between records, and are ignored. Wrapped (multi-line) FASTQ is rejected; FASTA sequences may span any number of lines, and text before the first FASTA header (such as `;` comments) is skipped
- A malformed record, a file that ends inside a FASTQ record, a truncated or corrupt gzip file, or NUL bytes where a record should start (such as the zero-filled block a crash can leave) fail the query rather than returning the records before it. The error names the file, and either the malformed record or how many records were read before the read failed; for damaged gzip data that count can stop short of the damage, since the data is inflated in blocks
- Reads files in parallel, one thread per file (per pair of files for paired input), up to DuckDB's `threads` setting; stdin is read on one thread
- An `http://` or `https://` file whose connection drops while it opens or partway through is retried by DuckDB's HTTP settings, `http_retries`, `http_retry_wait_ms` and `http_retry_backoff`. A retried read asks for the same byte range again, not the whole file. An HTTP error status DuckDB does not treat as transient, such as a 404, or a file whose ETag changed since it opened, fails at once. The defaults give up within a few seconds, which a long download will outlast; for one, raise them, for example `SET http_retries = 8; SET http_retry_wait_ms = 500; SET http_retry_backoff = 2;`. DuckDB also limits retries per host, so once many requests to one host have failed, later drops to that host can fail at once, even on a file that would have recovered
- Each local regular file can be decompressed and parsed on a thread of its own, ahead of the rows being produced. This is decided once for the whole scan: it happens only if every file in the scan, each taking one of DuckDB's threads plus one more thread (two for a pair), fits in `threads` at once. With `SET threads = 8` that is up to four single-end files or two pairs. The count covers every file even when DuckDB reads them one at a time, as it does for a sink that keeps insertion order (a `CREATE TABLE ... AS` under the default `preserve_insertion_order`). Each read-ahead stops once it holds 4 MiB of records, or `max_batch_bytes` if that is smaller, passing it by about a record (by up to sixteen where record lengths jump). Remote files and pipes (stdin from a pipe included) are read without it, as is everything in DuckDB-Wasm
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
- `reference_lengths` (VARCHAR, optional): Table or view name containing reference sequences for headerless SAM files. Must have at least 2 columns: first column = reference name (VARCHAR, BIGINT or UUID, matched as its text (BIGINT in decimal, UUID in the lowercase hyphenated form)), second column = reference length (INTEGER/BIGINT). Column names don't matter. Views are fully supported and can include computed columns.
- `include_filepath` (BOOLEAN, optional, default false): Add filepath column to output
- `include_seq_qual` (BOOLEAN, optional, default false): Add sequence and quality score columns to output. When enabled, primary alignments (non-secondary, non-supplementary) and unmapped reads must have SEQ/QUAL data or an error will be raised.

**Output schema includes:**
- `position` (BIGINT): 1-based start position
- `stop_position` (BIGINT): 1-based stop position (computed from CIGAR using `bam_endpos`)
- `cigar` (VARCHAR): CIGAR string; `''` for a record whose CIGAR field is `*`, never NULL
- `sequence` (VARCHAR, optional): Read sequence from SEQ field (when include_seq_qual=true)
- `qual` (UTINYINT[], optional): Quality scores as array of integers 0-93 (when include_seq_qual=true)
- Plus other standard SAM fields and optional tags

**Behavior:**
- Supports both SAM (text) and BAM (binary) formats
- Auto-detects file format from content
- Supports gzip-compressed SAM files
- Supports stdin input using `-` or `/dev/stdin` (single file only, not in arrays)
- Reads files in parallel, one thread per file, up to DuckDB's `threads` setting; stdin is read on one thread
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

**Malformed records fail the scan.** A record htslib cannot parse — a malformed field, a truncated record, a reference id outside the header's range — raises `Failed to read SAM/BAM record: the file is truncated or malformed`, with the line number for text SAM. Earlier versions dropped such a record and carried on, so a corrupt file read back as a silently short table. If you need to salvage what is readable from a damaged file, repair it first (e.g. `samtools view`), rather than relying on the reader to skip past the damage.

One case is *not* detectable and matches `samtools` exactly: a BGZF file truncated at a block boundary looks like a cleanly-ended file to htslib, so it reads back as a short table with no error. `samtools quickcheck <file>` is the tool for that — it verifies the BGZF EOF marker and exits non-zero when it is missing.

#### `read_alignment_header(path)`

Read a SAM/BAM **header's `@SQ` lines** as rows of `(tid, reference, length)`. This is the only way to
see a reference's `tid` — its index in the `@SQ` list — from SQL.

```sql
SELECT tid, reference, length FROM read_alignment_header('reads.bam') ORDER BY tid;
```

**Output schema:**
- `tid` (INTEGER): the reference's index in the `@SQ` list — this *is* the record's `refID` on disk
- `reference` (VARCHAR): reference name (`SN:`)
- `length` (BIGINT): reference length (`LN:`). BIGINT, not INTEGER: `LN` can exceed 2^31

**Why it exists.** A BAM is coordinate-sorted by **tid**, not by reference name, and no other function
exposes tid. `read_alignments` / `read_sam` decode each record's `reference` through the file's own
header, and name↔tid is a bijection fixed when the header is built — so reading names back yields
ascending *names* whatever the `@SQ` order was. That makes tid ordering unobservable through the record
readers, which is why the `@SQ` ordering contract in [Writing](writing.md#sam-and-bam) could not be
checked for BAM from SQL before this function existed.

**Behavior:**
- Reads **only** the header. `sam_read1` is never called, so a file with intact `@SQ` lines and a
  corrupt *record* is read successfully here even though `read_alignments` fails the scan on it.
- A header with **no `@SQ` lines returns zero rows, not an error** — unaligned BAM/SAM legitimately has
  no references, including the uBAM files `COPY … (FORMAT UBAM)` writes.
- Works on both BAM and text SAM.
- **Local paths only.** `https://`, `s3://` and friends are rejected at bind with a clear error; use
  `read_alignments` if you need remote streaming, or download the file first. Single path only — no
  glob, no list.

**Verifying a coordinate-sorted BAM:**
```sql
-- @SQ is emitted name-sorted, so tid order must equal name order.
SELECT COUNT(*) = 0 AS header_is_name_sorted FROM (
  SELECT tid, ROW_NUMBER() OVER (ORDER BY reference) - 1 AS expected_tid
  FROM read_alignment_header('sorted.bam')
) WHERE tid != expected_tid;
```

### PAF

Read PAF (Pairwise mApping Format) files, the tab-separated alignment format minimap2 and many long-read tools write, through the `read_paf` table function. Each row is one PAF record.

**Basic example:**

```sql
SELECT * FROM read_paf('alignments.paf');
```

**Parameters:**
- `filename` (VARCHAR or VARCHAR[]): Path to PAF file(s), or a glob pattern (e.g., `'data/*.paf.gz'`)
  - **Glob patterns**: When a single VARCHAR contains glob characters (`*`, `?`, `[`), files are expanded and sorted alphabetically
  - **Arrays**: VARCHAR[] elements are treated as literal paths (no glob expansion)
- `include_filepath` (BOOLEAN, optional, default false): Add a `filepath` column to the output

**Output schema:**
- `read_id` (VARCHAR): Query sequence name (PAF column 1)
- `query_length` (BIGINT): Query sequence length
- `query_start` (BIGINT): Start of the aligned block on the query, 1-based; PAF's 0-based start `+ 1`
- `query_stop` (BIGINT): End of the aligned block on the query, exclusive; PAF's end `+ 1`
- `strand` (VARCHAR): `+` or `-`; `*` for an unmapped record
- `reference` (VARCHAR): Target sequence name; `*` for an unmapped record
- `reference_length` (BIGINT): Target sequence length
- `position` (BIGINT): Start of the aligned block on the target, 1-based; PAF's 0-based start `+ 1`
- `stop_position` (BIGINT): End of the aligned block on the target, exclusive; PAF's end `+ 1`
- `residue_matches` (BIGINT): Number of matching bases (PAF column 10)
- `block_length` (BIGINT): Alignment block length (PAF column 11)
- `mapq` (UTINYINT): Mapping quality; 255 means missing
- minimap2's tags, one typed column each, NULL when the record does not carry it. These are the tags its manual lists for PAF, except `SA`, which minimap2 writes only in SAM; the names follow `read_alignments`' `tag_<name>` columns:

  | Column | Type | Tag |
  |---|---|---|
  | `tag_tp` | VARCHAR | `tp:A`, alignment type: `P` primary, `S` secondary, `I` / `i` inversion |
  | `tag_cm` | BIGINT | `cm:i`, minimizers on the chain |
  | `tag_s1` | BIGINT | `s1:i`, chaining score |
  | `tag_s2` | BIGINT | `s2:i`, chaining score of the best secondary chain |
  | `tag_nm` | BIGINT | `NM:i`, mismatches and gaps in the alignment |
  | `tag_md` | VARCHAR | `MD:Z`, for regenerating the reference sequence of the alignment |
  | `tag_as` | BIGINT | `AS:i`, DP alignment score |
  | `tag_ms` | BIGINT | `ms:i`, DP score of the best-scoring segment |
  | `tag_nn` | BIGINT | `nn:i`, ambiguous bases in the alignment |
  | `tag_ts` | VARCHAR | `ts:A`, transcript strand (splice mode) |
  | `tag_cg` | VARCHAR | `cg:Z`, CIGAR |
  | `tag_cs` | VARCHAR | `cs:Z`, difference string |
  | `tag_dv` | DOUBLE | `dv:f`, approximate per-base divergence |
  | `tag_de` | DOUBLE | `de:f`, gap-compressed per-base divergence |
  | `tag_rl` | BIGINT | `rl:i`, length of query regions with repetitive seeds |
  | `tag_zd` | BIGINT | `zd:i`, Z-drop break: bit 1 left, bit 2 right |
- `tags` (MAP(VARCHAR, VARCHAR)): Every other optional field after column 12, as `tag → value`
- `filepath` (VARCHAR, optional): File path when `include_filepath=true`

**Behavior:**
- Both intervals move to miint's 1-based half-open convention ([Coordinate conventions](#coordinate-conventions)): the target interval matches `read_alignments`, so `stop_position - position` is the aligned target length, and `query_stop - query_start` is the aligned query length. The query interval is on the read as sequenced, even on the `-` strand, as PAF writes it (see [Reads](#reads)).
- An unmapped record (`minimap2 --paf-no-hit`) reports `strand` and `reference` as `*` and `position` / `stop_position` as 0, as `read_alignments` does for an unplaced unmapped SAM record (RNAME `*`, POS 0). It has no aligned block on the query either, so `query_start` / `query_stop` are 0 too.
- minimap2's tags are typed as they are read: `NM:i:3` becomes `tag_nm = 3`. Any other tag goes into `tags` as text, `zz:i:7` as `tags['zz'] = '7'`, its type letter checked against its value and dropped; a single map has one value type, so cast these when you need numbers. A float (`f`) may also be `nan`, `inf` or `infinity`, in any case and with an optional sign, as htslib (and so `read_alignments`) reads one; a float too large for a DOUBLE reads as an infinity, as htslib's `strtod` gives. Untyped `key:value` fields, which some tools write instead of SAM-style tags, are accepted as well. An untyped field with one of minimap2's names is read as its column's type.
- gzip-compressed files (`.gz`) are decompressed transparently, including multi-member (bgzip) files. Remote paths (`s3://`, `https://`) are read through DuckDB's file system.
- One thread reads each file; several files are read in parallel.
- An empty file returns no rows.
- stdin is not supported.

**Malformed records fail the scan.** The error names the file, the line and the defect:
- fewer than 12 fields, or a blank line
- an empty query or target name
- a field that is not a whole non-negative integer, or a mapping quality above 255
- an interval that runs backwards or past its sequence's length
- a strand other than `+` or `-`, except `*` on an unmapped record
- `*` as only one of strand and target name, or an unmapped record with a non-zero query or target start or end
- a query or target end too large to take the `+ 1`
- a malformed tag, or a tag that appears twice on one line
- one of minimap2's tags with a type letter other than the one in the table, or a value that does not fit its column (an integer beyond 64 bits)

**Examples:**
```sql
-- Primary alignments of reads with at most 2% divergence
SELECT read_id, reference, position, stop_position
FROM read_paf('reads.paf.gz')
WHERE tag_tp = 'P' AND tag_de <= 0.02;

-- Fraction of each read covered by its alignments, pooling fragments on the read
SELECT read_id, query_length,
       list_sum([iv.stop - iv.start FOR iv IN compress_intervals(query_start, query_stop)]) / query_length
           AS fraction_aligned
FROM read_paf('reads.paf')
WHERE reference <> '*'
GROUP BY read_id, query_length;

-- Unaligned reads
SELECT read_id FROM read_paf('reads.paf') WHERE reference = '*';
```

### GFA

Read GFA 1 assembly graphs (GFA 1.0–1.2), the format most assemblers use for their contig graphs. `read_gfa` returns every record as it is written; `read_gfa_segments` and `read_gfa_links` give typed views of the two record types every graph has.

**Basic example:**

```sql
SELECT * FROM read_gfa('assembly_graph.gfa');
```

#### `read_gfa(filename)`

**Parameters:**
- `filename` (VARCHAR or VARCHAR[]): Path to GFA file(s), or a glob pattern (e.g., `'graphs/*.gfa'`)
  - **Glob patterns**: When a single VARCHAR contains glob characters (`*`, `?`, `[`), files are expanded and sorted alphabetically
  - **Arrays**: VARCHAR[] elements are treated as literal paths (no glob expansion)
- `include_filepath` (BOOLEAN, optional, default false): Add a `filepath` column to the output

**Output schema:**
- `line_number` (BIGINT): 1-based line of the record in its file. Comment lines are counted, so this matches what an editor shows; it is also the record order, which matters for records such as myloasm's read placements
- `record_type` (VARCHAR): The first field: `H`, `S`, `L`, `C`, `P`, `W`, `J`, or any other printable text as written
- `fields` (VARCHAR[]): The record type's mandatory positional fields, verbatim (table below). For a record type GFA 1 does not define, every field after the type
- `tags` (MAP(VARCHAR, VARCHAR), nullable): Optional `TAG:TYPE:VALUE` fields as `tag → value`. Empty for a known record with no tags; NULL for an unknown record type, whose fields are never split into positional fields and tags
- `filepath` (VARCHAR, optional): File path when `include_filepath=true`

| Record type | Positional fields |
|---|---|
| `H` header | none |
| `S` segment | name, sequence (`*` when not stored) |
| `L` link | from, from orientation, to, to orientation, overlap (CIGAR or `*`) |
| `C` containment | container, container orientation, contained, contained orientation, position, overlap |
| `P` path | path name, segment names, overlaps |
| `W` walk | sample, haplotype index, sequence name, start, end, walk |
| `J` jump | from, from orientation, to, to orientation, distance (or `*`) |

**Behavior:**
- Lines starting with `#` are comments and produce no row.
- `fields` are the file's own text. No coordinate conversion is applied: a `W` walk's start and end, and a `C` containment's position, stay 0-based as GFA writes them. See [Coordinate conventions](#coordinate-conventions).
- Tags must be typed (`TAG:TYPE:VALUE`); the type letter is validated and dropped. `read_gfa` keeps every record type, and each type has its own tags, so they stay in one text map: cast values when you need numbers, `tags['LN']::BIGINT`. `read_gfa_segments` and `read_gfa_links` type the specification's tags as columns. A float (`f`) may also be `nan`, `inf` or `infinity`, in any case and with an optional sign, as htslib reads one. A tag name may be any identifier (letters, digits and `_`, not starting with a digit), not only the specification's two characters: assemblers write longer ones, such as myloasm's `MEDIAN_DP1:f:33.0`.
- A header that declares a version other than GFA 1 (`VN:Z:1` or `VN:Z:1.x`) is rejected rather than read as GFA 1. GFA 2 (`VN:Z:2.0`) redefines the record types.
- gzip-compressed files (`.gz`) are decompressed transparently, including multi-member (bgzip) files. Remote paths are read through DuckDB's file system.
- One thread reads each file; several files are read in parallel.
- An empty file returns no rows. stdin is not supported.

**Malformed records fail the scan.** The error names the file, the line and the defect: a blank line, an empty record type, a record type that is not printable text (as when a gzipped file is named without `.gz`), a known record type with fewer positional fields than it needs, an `L`, `C` or `J` orientation other than `+` or `-`, a malformed or untyped tag, or a tag that appears twice on one line.

#### `read_gfa_segments(filename)`

The `S` records of a GFA file.

- `segment_id` (VARCHAR): Segment name
- `sequence` (VARCHAR, nullable): Segment sequence; NULL when the file stores `*`
- `tag_ln`, `tag_rc`, `tag_fc`, `tag_kc` (BIGINT, nullable): The specification's `LN` (segment length), `RC` (read count), `FC` (fragment count) and `KC` (k-mer count) tags
- `tag_sh` (VARCHAR, nullable): `SH`, the sequence's SHA-256 checksum, in hex
- `tag_ur` (VARCHAR, nullable): `UR`, a URI or path for the sequence
- `tags` (MAP(VARCHAR, VARCHAR)): The segment's other tags, e.g. an assembler's depth tags

#### `read_gfa_links(filename)`

The `L` records of a GFA file.

- `from_segment` (VARCHAR), `from_orientation` (VARCHAR): The link's source segment and its orientation, `+` or `-`
- `to_segment` (VARCHAR), `to_orientation` (VARCHAR): The link's target segment and its orientation
- `overlap` (VARCHAR, nullable): The overlap CIGAR; NULL when the file stores `*`
- `tag_mq`, `tag_nm`, `tag_rc`, `tag_fc`, `tag_kc` (BIGINT, nullable): The specification's `MQ` (mapping quality), `NM` (mismatches and gaps), `RC` (read count), `FC` (fragment count) and `KC` (k-mer count) tags
- `tag_id` (VARCHAR, nullable): `ID`, the link's identifier
- `tags` (MAP(VARCHAR, VARCHAR)): The link's other tags

Both views read through `read_gfa`, so a malformed record of any type fails either one, as it fails `read_gfa`. A typed tag column is NULL when the record does not carry the tag. Read with a value that is not a 64-bit integer (a file may give `LN` a type other than `i`), it is an error naming the record. Both views accept the same paths as `read_gfa` (a path, a glob, or a list); for a `filepath` column, query `read_gfa` directly. The functions named `_miint_gfa_*` in `duckdb_functions()` are internal helpers.

**Examples:**
```sql
-- Segment lengths, from the sequence when it is stored and from the LN tag otherwise
SELECT segment_id, coalesce(length(sequence), tag_ln) AS length
FROM read_gfa_segments('assembly_graph.gfa');

-- Every segment linked to one segment, in either direction
SELECT DISTINCT CASE WHEN from_segment = 'contig_1' THEN to_segment ELSE from_segment END AS neighbour
FROM read_gfa_links('assembly_graph.gfa')
WHERE 'contig_1' IN (from_segment, to_segment);

-- Which record types a graph contains
SELECT record_type, count(*) AS n FROM read_gfa('assembly_graph.gfa') GROUP BY ALL;

-- Paths as lists of oriented segments
SELECT fields[1] AS path_name, string_split(fields[2], ',') AS steps
FROM read_gfa('assembly_graph.gfa')
WHERE record_type = 'P';
```

### myloasm

Read one [myloasm](https://github.com/bluenote-1577/myloasm) output directory (its `-o`) as tables: which reads each contig is built from, where every read maps, how contigs connect, and where chimeric reads were split. Each function takes the directory and reads only the files it needs:

| Function | Reads | One row per |
|---|---|---|
| [`read_myloasm_contigs(dir)`](#read_myloasm_contigsdir) | `assembly_primary.fa`, `alternate_assemblies/assembly_alternate.fa`, `alternate_assemblies/duplicated_contigs.fa`, `final_contig_graph.gfa` | contig |
| [`read_myloasm_layout(dir)`](#read_myloasm_layoutdir) | `final_contig_graph.gfa` (`a` lines), breakpoint logs | read slice in a contig's backbone |
| [`read_myloasm_mappings(dir)`](#read_myloasm_mappingsdir) | `3-mapping/map_to_unitigs.paf.gz`, breakpoint logs | alignment of a read to a contig |
| [`read_myloasm_links(dir)`](#read_myloasm_linksdir) | `final_contig_graph.gfa` (`L` lines) | contig graph edge |
| [`read_myloasm_splits(dir)`](#read_myloasm_splitsdir) | `0-cleaning_and_unitigs/read_coverages.txt.gz`, `0-cleaning_and_unitigs/remap_temp/read_coverages.txt.gz` | chimera breakpoint |
| [`read_myloasm_run(dir)`](#read_myloasm_rundir) | `myloasm*.log` | myloasm invocation |

The formats were verified against myloasm 0.6.0's and 0.7.0's source and their real output; another version may write them differently. The two versions differ only in the GFA's `L` lines (see `read_myloasm_links`). What the values mean (the depth tiers, thresholds and split conditions below) comes from myloasm's source. No function returns sequences. Contig sequences, alternates included, are in the FASTA files: join `read_fastx` on `fasta_id`.

myloasm tells strains apart with **SNPmers**: k-mers that span a variable site. Two reads, or a read and a contig, that share many SNPmers and contradict few come from the same strain.

**Coordinates.** Every interval is 1-based half-open, as everywhere in miint ([Coordinate conventions](#coordinate-conventions)). Contig intervals (`position`, `stop_position`) are on the **unpolished** contig: the coordinates of the GFA and the PAF. Polishing shifts positions slightly, and trims the overlapping ends of a circular contig, whose unpolished length then exceeds its polished length by the overlap. Use the polished FASTA for which contig a read belongs to, not for positions in it. Read intervals (`query_start`, `query_stop`) are on the **original read as sequenced**. myloasm writes every interval 0-based half-open, and its read intervals on something other than the read, so these functions convert them:
- For a `-` alignment, myloasm writes the query interval on the reverse complement of the read. It is flipped back.
- myloasm splits reads it judges chimeric and renames each piece `<read>+split<i>`. A second cleaning round can split a piece again, giving `<read>+split<i>+split<j>`. Coordinates are written relative to the piece. `read_id` is the read's own name, `split_suffix` keeps the suffix, `piece_offset` is the number of read bases before the piece (0 for the first), and `query_start` / `query_stop` are moved onto the read. A `+split0` suffix does not by itself mean the read was split. In each cleaning round, myloasm checks for chimeras only the reads no other read contains (its "outer" reads); each of those gets a suffix, and one with no breakpoint becomes a single full-length `+split0` piece.

See [Coordinate conventions](#coordinate-conventions).

**Behavior:**
- Every member file a function reads is required. A missing one fails with the reader's file error, which names it.
- A directory path containing `*`, `?` or `[` is refused: the readers would expand it as a glob, so `run[1]` would read a directory named `run1`.
- Values that are not in myloasm's format fail the scan, naming the record. This covers headers, `a` lines, tags, breakpoint lines and PAF names, as well as files that contradict each other: a FASTA contig with no `S` record, or a split piece no breakpoint creates. The checks fire whichever columns are selected.
- myloasm writes its three FASTA files even when it has no contigs for one, so an empty one is read as no contigs.
- myloasm drops reads under 1 kb and reads it finds repetitive before assembling. They appear in no output.
- A read can map to several contigs (strain variants, repeats), so `read_myloasm_mappings` does not partition the reads.
- A read whose own name ends in `+split<number>` cannot be told apart from a myloasm piece name.
- The functions named `_miint_myloasm_*` in `duckdb_functions()` are internal helpers.

#### `read_myloasm_contigs(dir)`

Every contig in the GFA, with its FASTA header's values. myloasm writes the GFA before its last filter, so a contig dropped after that is in the GFA and in no FASTA: `filtered`.

- `contig_id` (VARCHAR): The contig's name, e.g. `u491ctg`
- `dereplication_status` (VARCHAR): `primary`, `alternate`, `duplicated` (which FASTA the contig is in) or `filtered` (none)
- `fasta_id` (VARCHAR, nullable): The contig's FASTA record name, which is `read_fastx`'s `read_id`
- `polished_length` (BIGINT, nullable): Length of the polished FASTA sequence
- `unpolished_length` (BIGINT): The GFA's `LN`, the length every coordinate in these tables refers to
- `circularity` (VARCHAR, nullable): `yes`, `possibly` or `no`
- `depth_identity_99`, `depth_identity_99_75`, `depth_identity_100` (DOUBLE): The contig's depth at three SNPmer-identity tiers: the median, over its reads, of each read's minimum depth counting only reads at ≥ 99%, ≥ 99.75% and 100% identity to it. These are the GFA's `DP`, `DP2` and `DP3` and the header's `depth-<d1>-<d2>-<d3>`.
- `median_depth_identity_100` (DOUBLE): The median of the reads' median depths at the 100% tier: the GFA's `MEDIAN_DP1`, although the GFA's read placements use `DP1` for the 99% tier
- `kmer_multiplicity` (DOUBLE, nullable) and `kmer_multiplicity_flag` (VARCHAR, nullable): The header's `mult=` and `duplicated-<no|possibly|yes>`: how repetitive the contig's own k-mers are (below 1.10 is `no`, below 1.50 `possibly`). **This is not the dereplication status**: a contig in `duplicated_contigs.fa` can say `duplicated-no`.
- `secondary_to_contig_id` (VARCHAR, nullable): For an alternate or duplicated contig, the contig myloasm compared it to and kept in its place. That contig can itself be an alternate.
- `secondary_aligned_percent`, `secondary_ani_percent` (DOUBLE, nullable): The aligned fraction and ANI of that comparison, in percent

#### `read_myloasm_layout(dir)`

The read slices each unpolished contig is built from, in order: the contig is these slices laid end to end. Only backbone reads are listed, a small fraction of the reads; `read_myloasm_mappings` has every read.

- `contig_id` (VARCHAR), `layout_order` (BIGINT): The contig and the slice's 1-based place in it
- `position`, `stop_position` (BIGINT): The slice on the unpolished contig. Consecutive slices abut, and the last stops at `unpolished_length + 1`. A slice can be empty when its neighbours cover the read entirely.
- `read_id` (VARCHAR), `split_suffix` (VARCHAR, nullable), `piece_offset` (BIGINT), `piece_length` (BIGINT): The read, and the piece of it this is
- `query_start`, `query_stop` (BIGINT): The slice on the original read
- `strand` (VARCHAR): `-` when the read runs reverse-complemented along the contig
- `read_depth_identity_99`, `read_depth_identity_99_75`, `read_depth_identity_100` (DOUBLE): The read's minimum depth at each identity tier
- `next_overlap_length`, `next_shared_snpmers`, `next_differing_snpmers` (BIGINT, nullable): The overlap with the next read in the layout. NULL on a contig's last read, where myloasm writes 0.

#### `read_myloasm_mappings(dir)`

myloasm's final mapping of every read, contained ones included, to the unpolished contigs. It is the input to polishing.

- `read_id` (VARCHAR), `split_suffix` (VARCHAR, nullable), `piece_offset` (BIGINT), `piece_length` (BIGINT): The read, and the piece of it that was mapped
- `query_start`, `query_stop` (BIGINT): The aligned block on the original read
- `strand` (VARCHAR): `+` or `-`
- `contig_id` (VARCHAR), `position`, `stop_position` (BIGINT): The aligned block on the unpolished contig
- `query_anchor_span`, `contig_anchor_span` (BIGINT): PAF columns 10 and 11, which myloasm fills with the distance from the first anchor's start to the last anchor's end, not residue matches and block length
- `shared_minimizers`, `shared_snpmers`, `differing_snpmers` (BIGINT): The strain-level evidence: minimizers the read shares with the contig, and SNPmers it shares and contradicts

myloasm's own SNPmer identity estimate for a mapping is `1 / (1 + differing_snpmers / shared_minimizers / c)`, with `c` from `read_myloasm_run`'s `compression`.

#### `read_myloasm_links(dir)`

The contig graph's edges, from the GFA's `L` lines.

- `from_contig_id`, `from_orientation`, `to_contig_id`, `to_orientation` (VARCHAR): The edge
- `overlap_length` (BIGINT): Overlap between the two contig ends, in bases
- `shared_snpmers`, `differing_snpmers`, `shared_minimizers` (BIGINT): The overlap's evidence
- `snpmer_identity_percent` (DOUBLE, nullable): myloasm's identity estimate for the overlap; NULL when the ends share no minimizers
- `from_read_id`, `from_split_suffix`, `to_read_id`, `to_split_suffix` (VARCHAR): The reads at the two ends that meet. Leaving a contig on `+` leaves from its last read and on `-` from its first; entering on `+` enters at its first read and on `-` at its last.

myloasm 0.6.0 writes `L` lines without these tags, so for its links every column after the edge is NULL. The overlap length is unknown there, not 0: the GFA's overlap field is a `0M` placeholder in every version. A link that has some of the tags but not all of them is an error, as is a tagless link whose overlap is not `0M` and a file with both tagged and tagless links.

#### `read_myloasm_splits(dir)`

Every breakpoint myloasm split a read at, from both cleaning rounds' logs. Only reads with a breakpoint are listed.

- `read_id` (VARCHAR): The read
- `parent_split_suffix` (VARCHAR, nullable): The piece that was split, when round 2 split a round-1 piece; NULL when the read itself was split
- `cleaning_round` (INTEGER): 1 or 2
- `breakpoint_order` (BIGINT): 1-based order along the piece. Piece `+split<i>` starts at breakpoint `i`'s `query_stop` and ends at breakpoint `i + 1`'s `query_start`, or at the end of the parent. `+split0` starts with its parent. myloasm drops pieces of 1000 bases or fewer but still counts them, so the numbering can skip.
- `query_start`, `query_stop` (BIGINT): The bases myloasm discarded around the breakpoint, on the original read
- `split_condition` (INTEGER): myloasm's code for why it split there, e.g. `1` for a low-depth stretch between well-covered ones, and negative codes for a read end (`-1`/`-4` unmapped, `-2`/`-3` low depth)

#### `read_myloasm_run(dir)`

One row per `myloasm*.log`, one per invocation; a resumed run has several.

- `log_file` (VARCHAR): The log's path
- `started_at` (TIMESTAMP): When the invocation began, in the log's local time
- `myloasm_version` (VARCHAR), `command` (VARCHAR): myloasm's version and command line
- `kmer_size` (INTEGER), `compression` (INTEGER): The run's k-mer size k, and its compression ratio c: myloasm keeps 1 in c k-mers as minimizers. myloasm logs its command line, not these, so they are derived by its own rules, checked against real runs. Without `--nano-r9`, k is `-k`/`--kmer-size` (default 21) and c is the hidden `--c` (default 11). `--nano-r9` sets k = 17 and c = 7, overriding `-k` and `--c`. `-c`/`--compression` then sets c in every case. The log's "Setting -c to N" line prints c before `-c` is applied. Arguments containing spaces cannot be told apart in the logged command line.
- `hifi` (BOOLEAN): Whether `--hifi` was given

**Examples:**
```sql
-- Reads a group of contigs was assembled from, gated on strain-level agreement
SELECT DISTINCT read_id
FROM read_myloasm_mappings('myloasm_out')
WHERE contig_id IN ('u491ctg', 'u826ctg') AND differing_snpmers = 0;

-- myloasm's own SNPmer identity estimate for every mapping
SELECT read_id, contig_id,
       1 / (1 + differing_snpmers / shared_minimizers
                / (SELECT compression FROM read_myloasm_run('myloasm_out') ORDER BY started_at DESC LIMIT 1))
           AS snpmer_identity
FROM read_myloasm_mappings('myloasm_out');

-- A contig's neighbours in the assembly graph
SELECT DISTINCT CASE WHEN from_contig_id = 'u491ctg' THEN to_contig_id ELSE from_contig_id END AS neighbour
FROM read_myloasm_links('myloasm_out')
WHERE 'u491ctg' IN (from_contig_id, to_contig_id);

-- Reads that map to no contig, e.g. to recruit against another sample's contigs
SELECT r.read_id
FROM read_fastx('reads.fastq.gz') r
ANTI JOIN read_myloasm_mappings('myloasm_out') m ON m.read_id = r.read_id;

-- The bases of each read that align, cut from the read itself
SELECT m.read_id, m.contig_id, substr(r.sequence1, m.query_start, m.query_stop - m.query_start) AS aligned
FROM read_myloasm_mappings('myloasm_out') m
JOIN read_fastx('reads.fastq.gz') r USING (read_id);

-- Primary contigs with their sequences
SELECT c.contig_id, c.depth_identity_100, f.sequence1
FROM read_myloasm_contigs('myloasm_out') c
JOIN read_fastx('myloasm_out/assembly_primary.fa') f ON f.read_id = c.fasta_id;
```

**Assembling pooled reads with myloasm.** myloasm 0.7.0 sizes its Bloom filter by applying a ×1.5 gzip factor to its running total once per input file, so the estimate compounds with the number of files. By its sizing code (`seq_parse.rs`, `estimate_bf_size`), 50 GB of gzipped reads is sized at 37.5 GB as one file, about 76 GB as four, and the 200 GB cap as ten, and that size is allocated. When co-assembling reads gathered with these tables, write them to one FASTQ or set `-b` explicitly.

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
- Reads up to 4 files at a time, one thread per file
- Zero values are not returned (sparse representation)
- Supports reading multiple files which are concatenated in the output
- A local file is checked when the query is planned: a missing file, a file HDF5 cannot open and a file that is not BIOM are refused then. A file HDF5 cannot open, being unreadable or damaged, is refused with HDF5's reason (`read_biom: cannot open '<file>': H5Fopen failed: truncated file ...`). A remote file is checked when it is read, and a read that fails names the file and why (`read_biom: cannot read '<file>': ...`)

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
- `stop_position` (INTEGER): 1-based **half-open** end position — i.e. GFF3's `end` **+ 1**. See [Coordinate conventions](#coordinate-conventions) below.
- `score` (DOUBLE, nullable): Feature score (NULL if '.')
- `strand` (VARCHAR, nullable): Strand ('+', '-', or NULL if '.')
- `phase` (INTEGER, nullable): CDS phase (0, 1, 2, or NULL if '.')
- `attributes` (MAP(VARCHAR, VARCHAR)): Parsed key-value attributes

**Behavior:**
- Reads tab-delimited GFF3 format files
- Skips comment and directive lines (those starting with '#'); a '#' later in a line is kept as data
- Every line, comment lines included, must be UTF-8 and have at most nine tab-separated fields; a line that breaks either fails the scan. GFF3 has no quoting, so '"' is an ordinary character
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

-- Calculate feature lengths (half-open: no + 1)
SELECT type,
       AVG(stop_position - position) AS avg_length,
       MIN(stop_position - position) AS min_length,
       MAX(stop_position - position) AS max_length
FROM read_gff('annotations.gff')
GROUP BY type;

-- Extract CDS features with phase information
SELECT seqid, position, stop_position, phase, attributes['Parent'] AS parent_transcript
FROM read_gff('annotations.gff')
WHERE type = 'CDS' AND phase IS NOT NULL;

-- Find features overlapping the window [1000, 5000].
-- stop_position is half-open, so the end test is strict (>): an interval whose
-- stop_position equals 1000 ends just before the window and does not overlap it.
SELECT type, position, stop_position, attributes['ID'] AS feature_id
FROM read_gff('annotations.gff')
WHERE seqid = 'chr1'
  AND position <= 5000
  AND stop_position > 1000
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
- Coordinates **in the file** are 1-based and inclusive (both start and end). miint normalizes them on read — see [Coordinate conventions](#coordinate-conventions).
- Strand: '+' (forward), '-' (reverse), '.' (unknown/not applicable)
- Phase: Indicates position within codon (0, 1, or 2) for CDS features
- A line starting with '#' is a comment or a directive (such as '##gff-version'). GFF3 has no end-of-line comments, and a seqid must escape '#', so a '#' later in a line, such as in an attribute value, is data.
- The 9th column (attributes) may be empty; a feature line must carry the 8 mandatory fields
- **The GFF3 sequence section is honoured**: content after a `##FASTA` directive is not returned as feature rows. prokka and bakta always append the genome, so this is the common case. A feature line carrying only *some* of the 8 mandatory fields is a malformed GFF line and raises an error rather than being silently dropped, as does a feature line with an empty seqid, a start or end that is not an integer, or a score or phase that is neither a number nor `.`, whichever columns the query selects or filters on. (A query that stops early, such as one with `LIMIT`, may finish before reaching it.)

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

### Coordinate conventions

**Every miint interval is 1-based half-open**, on a reference and on a read alike, so the same base of a sequence has the same number wherever it appears and intervals from different functions compose without adjustment. A length is `stop - start`, never `+ 1`, and the bases of an interval are `substr(sequence, start, stop - start)`. A reader converts the file's own convention when it reads; the [exceptions](#exceptions) are listed below.

| Sequence | Columns | Interval |
|---|---|---|
| A reference, contig or annotated sequence | `position`, `stop_position` | `[position, stop_position)` |
| A read, in the orientation it was sequenced | `query_start`, `query_stop` | `[query_start, query_stop)` |

#### References

`position` / `stop_position` are used by `read_alignments`, `align_minimap2` / `align_bowtie2` and friends, `alignment_slice`, `compute_coverage_depth`, `genome_coverage` / `genome_coverage_per_sample`, `read_gff`, `read_ncbi_annotation`, `read_paf`, `read_myloasm_layout` and `read_myloasm_mappings`.

| Source convention | Formats | How miint converts it |
|---|---|---|
| 1-based closed `[start, end]` | GFF3, NCBI annotation | `position = start`, `stop_position = end + 1` |
| 0-based half-open `[start, end)` | PAF target columns | `position = start + 1`, `stop_position = end + 1` |
| 1-based `POS` plus the CIGAR's reference span | SAM / BAM | `position = POS`, `stop_position = POS + span` |

The value a file contained is always recoverable: GFF's closed `end` is `stop_position - 1`, and PAF's 0-based `start` is `position - 1`. An unplaced unmapped record, a SAM record with RNAME `*` and `POS` 0 or a PAF record with target `*`, reports `position = 0` and `stop_position = 0`.

> **Changed behaviour.** `read_gff` and `read_ncbi_annotation` previously emitted the source format's closed `end` directly under the same `stop_position` column name. Queries that consumed it now see a value one larger. Anything computing a length as `stop_position - position + 1` should drop the `+ 1`; anything that compensated with its own `+ 1` before feeding another miint function should remove that compensation. The motivation was that the mismatch was silent: feeding `read_gff` output to `genome_coverage` (which computes `SUM(stop - start)`) type-checked, ran, raised nothing, and under-counted every feature by exactly one base.

#### Reads

Read intervals say which bases of a read an alignment covers. `read_paf` reports them as `query_start` / `query_stop`, and the [myloasm](#myloasm) functions report them on the original read.

| Source convention | Formats | How miint converts it |
|---|---|---|
| 0-based half-open `[start, end)` on the read as sequenced | PAF query columns | `query_start = start + 1`, `query_stop = end + 1` |
| 0-based half-open on a read piece, on its reverse complement for `-` | myloasm's PAF, GFA `a` lines and breakpoint logs | mirrored back for `-`, moved onto the original read, then `+ 1` |

Read intervals are always on the read **as sequenced**, whatever strand the alignment is on. A reverse-strand alignment covering the read's last 3000 of 4500 bases is `[1501, 4501)`, not `[1, 3001)`. PAF already writes forward-strand query coordinates; myloasm does not, and its functions flip its `-` intervals. The aligned bases of a read are `substr(sequence1, query_start, query_stop - query_start)` on `read_fastx`'s `sequence1`, the read as sequenced. `read_alignments`' `sequence` is SAM's SEQ instead: reverse-complemented on a reverse-strand record, and missing any hard-clipped bases.

An unmapped PAF record has no aligned block, so it reports `query_start = 0` and `query_stop = 0`, like its `position` / `stop_position`.

#### Exceptions

- `cigar_query_intervals` and `extract_linked_amplicon` return **0-based** half-open read intervals as `start` / `stop`. Add 1 to both ends to compare them with `read_paf` or the myloasm functions. [#330](https://github.com/the-miint/duckdb-miint/issues/330) is to move them to the rule, which will add 1 to each value they return and may rename the fields, so that a query written for the old values fails instead of shifting.
- `blast` reports BLAST's own 1-based **closed** `query_start` / `query_end` and `subject_start` / `subject_end`; the subject pair is reversed for a minus-strand hit. `align_sortmerna_rrna` reports 1-based **closed** `ref_start` / `ref_end`. [#330](https://github.com/the-miint/duckdb-miint/issues/330) covers these too.
- `rype_extract_strand_minimizers` reports each minimizer's position as a 0-based offset into the sequence ([Classification](classification.md)).
- `read_gfa` returns a file's fields as written, so a `W` walk's start and end and a `C` containment's position stay 0-based. See [GFA](#gfa).
- Bin indices from `bin_of` / `bin_start` / `interval_bins` ([Genome bins](alignment_analysis.md#genome-bins)) are 0-based array positions, not coordinates.
- `absquant_orf_copies` takes 1-based **closed** ORF coordinates, the convention woltka's `coords.txt` uses; see [Absolute quantification](absolute_quantification.md).
