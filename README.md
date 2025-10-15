# Miint

This repository is based on https://github.com/duckdb/extension-template, check it out if you want to build and ship your own DuckDB extension.

---

This extension, Miint, allow you to ... <extension_goal>.


## Building
### Managing dependencies
DuckDB extensions uses VCPKG for dependency management. Enabling VCPKG is very simple: follow the [installation instructions](https://vcpkg.io/en/getting-started) or just run the following:
```shell
git clone https://github.com/Microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh
export VCPKG_TOOLCHAIN_PATH=`pwd`/vcpkg/scripts/buildsystems/vcpkg.cmake
```
Note: VCPKG is only required for extensions that want to rely on it for dependency management. If you want to develop an extension without dependencies, or want to do your own dependency management, just skip this step. Note that the example extension uses VCPKG to build with a dependency for instructive purposes, so when skipping this step the build may not work without removing the dependency.

### Build steps
Now to build the extension, run:
```sh
make
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

### `read_sam(filename)`
Read SAM/BAM alignment files.

**Output schema includes:**
- `position` (BIGINT): 1-based start position
- `stop_position` (BIGINT): 1-based stop position (computed from CIGAR using `bam_endpos`)
- `cigar` (VARCHAR): CIGAR string
- Plus other standard SAM fields and optional tags

### `read_fastx(filename)`
Read FASTA/FASTQ sequence files.

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
- `alignment_is_qc_failed(flags)` - QC failure (0x200)
- `alignment_is_duplicate(flags)` - PCR/optical duplicate (0x400)
- `alignment_is_supplementary(flags)` - Supplementary alignment (0x800)

**HTSlib-compatible aliases:**
`is_paired`, `is_proper_pair`, `is_unmapped`, `is_munmap`, `is_reverse`, `is_mreverse`, `is_read1`, `is_read2`, `is_secondary`, `is_qcfail`, `is_dup`, `is_supplementary`

**Example:**
```sql
SELECT read_id, flags
FROM read_sam('alignments.sam')
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
FROM read_sam('alignments.sam')
WHERE tag_nm IS NOT NULL;

-- Filter high-quality alignments
SELECT COUNT(*)
FROM read_sam('alignments.sam')
WHERE alignment_seq_identity(cigar, tag_nm, tag_md) > 0.95;

-- Compare different identity methods
SELECT read_id,
  alignment_seq_identity(cigar, tag_nm, tag_md, 'gap_excluded') AS gap_excl,
  alignment_seq_identity(cigar, tag_nm, tag_md, 'blast') AS blast,
  alignment_seq_identity(cigar, tag_nm, tag_md, 'gap_compressed') AS gap_comp
FROM read_sam('alignments.sam')
WHERE tag_nm IS NOT NULL AND tag_md IS NOT NULL;
```

**Reference:** [On the definition of sequence identity](https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity) by Heng Li

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

### `COPY ... TO '...' (FORMAT SAM)`

Write query results to SAM format files. Requires all mandatory SAM columns from `read_sam` output.

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
- `INCLUDE_HEADER` (default: true): Include @SQ header lines with reference sequences
- `REFERENCE_LENGTHS` (required if INCLUDE_HEADER=true): MAP of reference names to their lengths (e.g., `MAP{'chr1': 248956422, 'chr2': 242193529}`)
- `COMPRESSION` (default: auto): Enable gzip compression (auto-detected from `.gz` extension)

**Examples:**
```sql
-- Basic SAM output with header (requires reference lengths)
COPY (SELECT * FROM read_sam('input.sam'))
TO 'output.sam' (FORMAT SAM, REFERENCE_LENGTHS MAP{'chr1': 248956422, 'chr2': 242193529});

-- Using a variable for reference lengths (recommended for reuse)
SET VARIABLE my_refs = MAP{'chr1': 248956422, 'chr2': 242193529};
COPY (SELECT * FROM read_sam('input.sam'))
TO 'output.sam' (FORMAT SAM, REFERENCE_LENGTHS getvariable('my_refs'));

-- Constructing reference lengths from a query/table
CREATE TABLE refs AS SELECT 'chr1' as name, 248956422 as length 
                     UNION ALL SELECT 'chr2', 242193529;
SET VARIABLE ref_map = (SELECT MAP(LIST(name), LIST(length)) FROM refs);
COPY (SELECT * FROM read_sam('input.sam'))
TO 'output.sam' (FORMAT SAM, REFERENCE_LENGTHS getvariable('ref_map'));

-- Headerless SAM output (no reference lengths needed)
COPY (SELECT * FROM read_sam('input.sam'))
TO 'output.sam' (FORMAT SAM, INCLUDE_HEADER false);

-- Compressed SAM output with header
COPY (SELECT * FROM read_sam('input.sam'))
TO 'output.sam.gz' (FORMAT SAM, COMPRESSION gzip, REFERENCE_LENGTHS MAP{'chr1': 248956422});

-- Filter and write high-quality alignments (headerless)
COPY (
  SELECT * FROM read_sam('input.sam')
  WHERE mapq >= 30 AND NOT alignment_is_unmapped(flags)
) TO 'filtered.sam' (FORMAT SAM, INCLUDE_HEADER false);
```

**Notes:**
- SEQ and QUAL fields are always written as `*` in current implementation
- Reference lengths must be provided explicitly when writing headers - they cannot be inferred from the data
- All optional tags present in the input are preserved in the output

## Running the tests
Different tests can be created for DuckDB extensions. The primary way of testing DuckDB extensions should be the SQL tests in `./test/sql`. These SQL tests can be run using:
```sh
make test
```

### Installing the deployed binaries
To install your extension binaries from S3, you will need to do two things. Firstly, DuckDB should be launched with the
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

Secondly, you will need to set the repository endpoint in DuckDB to the HTTP url of your bucket + version of the extension
you want to install. To do this run the following SQL query in DuckDB:
```sql
SET custom_extension_repository='bucket.s3.eu-west-1.amazonaws.com/<your_extension_name>/latest';
```
Note that the `/latest` path will allow you to install the latest extension version available for your current version of
DuckDB. To specify a specific version, you can pass the version instead.

After running these steps, you can install and load your extension using the regular INSTALL/LOAD commands in DuckDB:
```sql
INSTALL miint
LOAD miint
```
