# Writing files

MIINT provides custom COPY formats for writing bioinformatics file formats. The writers extend the native SQL methods for writing files through the `COPY ... TO '...' (FORMAT csv)` style. A quick example is below, where we read a FASTQ file and write it as FASTA:

```sql
COPY (SELECT * FROM 'some_file.fastq')
TO 'some_other_file.fasta' (FORMAT FASTA);
```

**FASTQ, FASTA, SAM, BAM and UBAM records come out in no fixed order unless the query has an `ORDER BY`.** Their records are a set (a SAM or BAM file is sorted only when its records are), so these formats write on all of DuckDB's threads by default, and the order can change from run to run. This is not what DuckDB's own formats do: CSV and Parquet keep the order rows come in. To write records in an order, sort them with `ORDER BY`; to keep the order they come in without sorting, such as the input file's own order, add `PRESERVE_ORDER true`. Each format's "Threads, order and memory" has the details.

Options follow DuckDB's rules for COPY options in every format here, as they do for CSV and Parquet:

- A boolean option given alone means true: `(FORMAT FASTQ, INTERLEAVE)` is `INTERLEAVE true`.
- An option the format does not have is refused, and the error may suggest options with similar names. The names of DuckDB's own COPY options, such as `OVERWRITE`, are not refused as unknown; DuckDB handles them, and refuses those a format cannot support, such as `FILE_SIZE_BYTES`, or `ROW_GROUPS_PER_FILE`, which in the record formats above would drop an `ORDER BY`.
- Other values are cast to the option's type, so `COMPRESSION_LEVEL '3'` is level 3 and `COMPRESSION_LEVEL 1.6` is level 2. An option with no value or a NULL value is refused.

`COMPRESSION` takes the same values in the text formats, FASTQ, FASTA, NEWICK and SAM, in any case: `'gzip'` or `'gz'`, and `'none'` or `'uncompressed'`. Without it, a path ending in `.gz` is written as gzip. zstd (`'zstd'`, or a `.zst` path) is refused. BIOM's `COMPRESSION` takes the same values but sets HDF5's compression inside the file, and BAM's takes only gzip; see [BIOM](#biom) and [SAM / BAM](#sam-and-bam).

FASTQ, FASTA and NEWICK write a path with a URI scheme other than `file://` (`s3://`, `https://`, `az://`, ...) through DuckDB's file systems, which load the extension that serves the scheme. For a scheme no file system serves, DuckDB opens the path as a local file, as it does for its own formats. SAM, BAM, UBAM and BIOM write local files only and refuse a remote target. In every format a `file:` URI is the local file it names, as DuckDB reads it (`file:///p`, `file://localhost/p` and `file:/p` are `/p`; on Windows, `file:///C:/p` is `C:/p`), and gzip output to it takes `COMPRESSION_LEVEL` as any local path does.

## Table of Contents

- [FASTQ](#fastq) - Write FASTQ sequence files
- [FASTA](#fasta) - Write FASTA sequence files
- [SAM / BAM](#sam-and-bam) - Write SAM/BAM alignment files
- [UBAM](#ubam) - Write unaligned reads BAM (uBAM) files
- [BIOM](#biom) - Write BIOM observation matrix files
- [NEWICK](#newick) - Write Newick phylogenetic trees
- [GFF](#gff) - Write GFF3 annotations, and Prodigal's GFF from `predict_genes`, as lines for `FORMAT CSV`

### FASTQ 

Write query results to FASTQ format files. Requires `read_id`, `sequence1`, and `qual1` columns from `read_fastx` output.

```sql
COPY (SELECT * FROM 'input.fastq')
TO 'output.fastq' (FORMAT FASTQ);
```

The writer treats the NULL mate columns from a single-end read as single-end input. Set `INCLUDE_COMMENT true` to preserve FASTQ header comments on a round trip.

**Required columns:**

- `read_id` (VARCHAR, BIGINT, or UUID): Sequence identifier. `BIGINT` is stringified to decimal and `UUID` to its canonical 36-char lowercase form; the record name on disk is always text. Round-trip via `read_fastx` returns `read_id` as `VARCHAR` (the `BIGINT`/`UUID` type is not preserved on disk — recover it with e.g. `read_id::BIGINT`).
- `sequence1` (VARCHAR): DNA/RNA sequence
- `qual1` (`UTINYINT[]`): Decoded Phred scores for R1

**Optional columns:**

- `comment` (VARCHAR): Comment line (only included if `INCLUDE_COMMENT=true`). Note this defaults to **false**, so a `read_fastx` → `COPY` round trip drops the comment unless you pass `INCLUDE_COMMENT true`. For Illumina data the comment carries the read number, filter flag and index (`1:N:0:ATCACG`), so pass it when a faithful round trip matters.
- `sequence_index` (BIGINT): Used as identifier if `ID_AS_SEQUENCE_INDEX=true`
- `sequence2` (VARCHAR): Second read for paired-end data
- `qual2` (`UTINYINT[]`): Decoded Phred scores for R2

**COPY options**

| Option | Default | Purpose |
| --- | --- | --- |
| `QUAL_OFFSET` | `33` | Encode quality scores with offset `33` or `64`. |
| `INCLUDE_COMMENT` | `false` | Include the comment from the input. |
| `ID_AS_SEQUENCE_INDEX` | `false` | Use `sequence_index` as the identifier. |
| `INTERLEAVE` | `false` | Write paired reads into one file. Incompatible with `{ORIENTATION}`. |
| `COMPRESSION` | Auto | Infer gzip from the `.gz` extension or set explicitly. |
| `COMPRESSION_LEVEL` | `1` | Set the local gzip compression level from `0` to `9`. |
| `ROW_GROUP_SIZE` | `64` | Rows in each batch DuckDB hands the writer. See "Threads, order and memory" below. |
| `ROW_GROUP_SIZE_BYTES` | `'8MB'` | Bytes in each batch, as DuckDB measures the rows. See "Threads, order and memory" below. |

**Compression**

Local gzip output uses BGZF (block gzip via HTSlib), which standard gzip tools can read. It uses DuckDB's thread count, does not change the record order (see "Threads, order and memory" below), supports block indexing, and runs about 10× faster than single-threaded gzip on multi-core hosts. Lower compression levels are faster and produce larger files. Remote targets use single-threaded gzip. `COMPRESSION_LEVEL` is available only for local gzip output.

**Single-end vs paired-end output:**

You do not need to project away `sequence2`/`qual2` or pre-declare whether the data is paired. `read_fastx` always emits `sequence2`/`qual2` (NULL for single-end), and the writer decides per record from whether those values are NULL:

- **`{ORIENTATION}` in the path** selects **split** output: R1 → the `R1` file, R2 → the `R2` file. If the data turns out to be entirely single-end, only the R1 file is written (no empty R2 file is created).
- **No `{ORIENTATION}`** selects a **single file**. Single-end data is written as-is. Paired-end data is interleaved only if `INTERLEAVE=true`.
- A record is **paired** when both `sequence2` and `qual2` are non-NULL, and **single-end** when both are NULL. A single COPY must be all one or all the other.

The writer fails loud (rather than guessing or silently dropping data) when:

- paired data is sent to a single-file path without `INTERLEAVE=true` (add `{ORIENTATION}` to split, or `INTERLEAVE=true` to interleave),
- `{ORIENTATION}` is combined with `INTERLEAVE=true` (contradictory),
- `sequence2` and `qual2` disagree on NULL-ness within a record, or
- the input mixes single-end and paired-end records.

An error about a record's data names its `read_id`: a NULL `sequence1`, `qual1` or (with `ID_AS_SEQUENCE_INDEX`) `sequence_index`; a NULL quality score; a score that would encode past `~` (126); a quality length that does not match the sequence; `sequence2` and `qual2` disagreeing on NULL; a paired record with nowhere to put R2.

**Threads, order and memory:**

- **No record order, by default.** Unless the query has an `ORDER BY`, the query and the writer use all of DuckDB's threads, and the records come out in an order that can change from run to run (DuckDB's own formats keep the order rows come in). Mates stay together: an interleaved file has each R2 right after its R1, and split R1/R2 files stay in step.
- **Write a pair's files in one COPY**, with `{ORIENTATION}` or `INTERLEAVE`. Two COPYs, one for R1 and one for R2, write each file in its own order, and the files stay in step only if both have the same `ORDER BY`, or both say `PRESERVE_ORDER true`.
- **An `ORDER BY` keeps its order when it sorts the rows written**: one at the end of the query, or at the end of a subquery or CTE that the query selects from, with nothing after it but `WHERE`, `LIMIT` and the select list (none of them holding a subquery). The writer runs on one thread; the sort still uses all threads. An `ORDER BY` whose rows feed a join, an aggregate, or an `IN`, `EXISTS` or scalar subquery sorts nothing the COPY writes, and the COPY runs on all threads.
- **The file is the same whatever `threads` is set to only when no two rows tie on the `ORDER BY`.** Rows that tie come out of DuckDB's sort in either order, so add a column that tells them apart, such as `read_id`.
- **`PRESERVE_ORDER true` keeps the order the rows come in**, without sorting: the input file's own order, for instance, or a table's. Records are written in the query's order, on one thread. The file is the same whatever `threads` is set to only when that order is fixed, as a file's or a table's is: an aggregate hands on its groups in an order that can change with `threads` and from run to run. The work that feeds the COPY directly runs on that thread too; work behind a sort, an aggregate or a join's build side still uses all threads.
- **`PRESERVE_ORDER false` drops an `ORDER BY`**, and the records come out unsorted. `SET preserve_insertion_order` makes no difference to these writers: an `ORDER BY` that sorts the rows written reaches the COPY, and DuckDB keeps its order whatever the setting.
- `PARTITION_BY` and `PER_THREAD_OUTPUT` always use all threads, and name their files `data_0.fastq` and so on (`data_0.fastq.gz` with gzip).
- **Batches.** DuckDB hands the writer its rows in batches. Unless the COPY sets them, `ROW_GROUP_SIZE 64` and `ROW_GROUP_SIZE_BYTES '8MB'` are added; a value given replaces its default, and the other default stays. A batch ends at whichever it reaches first: `ROW_GROUP_SIZE` rows (a pair is one row), or `ROW_GROUP_SIZE_BYTES`, measured as DuckDB's in-memory size of the rows rather than the text written. A batch can pass either: DuckDB adds rows a chunk at a time and cuts a chunk only at the byte cap, and keeps the row that crosses it. `FORMAT _miint_fastq`, the writer behind `FORMAT FASTQ`, adds no defaults. How a batch is written: [architecture notes](internals/architecture.md#record-writers-on-duckdbs-batch-api).
- **Memory.** The writer's text buffers sit outside DuckDB's `memory_limit`; the rows themselves are held by DuckDB and counted against it.

**More examples**

```sql
-- Paired-end interleaved output (single file)
COPY (SELECT * FROM read_fastx('R1.fastq', sequence2 = 'R2.fastq'))
TO 'output.fastq' (FORMAT FASTQ, INTERLEAVE true);

-- Paired-end split files (the {ORIENTATION} placeholder selects split output)
COPY (SELECT * FROM read_fastx('R1.fastq', sequence2 = 'R2.fastq'))
TO 'output_{ORIENTATION}.fastq' (FORMAT FASTQ);

-- Single-end data with {ORIENTATION} simply writes the R1 file (output_R1.fastq.gz)
COPY (SELECT * FROM read_fastx('input.fastq'))
TO 'output_{ORIENTATION}.fastq.gz' (FORMAT FASTQ, COMPRESSION gzip);

-- Compressed output with custom quality offset
COPY (SELECT * FROM read_fastx('input.fastq'))
TO 'output.fastq.gz' (FORMAT FASTQ, QUAL_OFFSET 33, COMPRESSION gzip);

-- Use sequence index as identifier
COPY (SELECT * FROM read_fastx('input.fastq'))
TO 'output.fastq' (FORMAT FASTQ, ID_AS_SEQUENCE_INDEX true);

-- Filter on all threads; the records come out in no fixed order, each R2 file in step with its R1
COPY (SELECT * FROM read_fastx('R1.fastq.gz', sequence2 = 'R2.fastq.gz') WHERE length(sequence1) >= 100)
TO 'filtered_{ORIENTATION}.fastq.gz' (FORMAT FASTQ);

-- Records sorted by name
COPY (SELECT * FROM read_fastx('input.fastq') ORDER BY read_id)
TO 'sorted.fastq' (FORMAT FASTQ);

-- Records in the input files' own order, without sorting
COPY (SELECT * FROM read_fastx('run_*.fastq.gz'))
TO 'all_runs.fastq.gz' (FORMAT FASTQ, PRESERVE_ORDER true);
```

<!-- BEGIN REFERENCE: copy_fasta -->
### FASTA

Write query results to FASTA format files. Requires `read_id` and `sequence1` columns from `read_fastx` output.

```sql
COPY (SELECT * FROM 'input.fasta')
TO 'output.fasta' (FORMAT FASTA);
```

**Required columns:**

- `read_id` (VARCHAR, BIGINT, or UUID): Sequence identifier. `BIGINT` is stringified to decimal and `UUID` to its canonical 36-char lowercase form; the record name on disk is always text. Round-trip via `read_fastx` returns `read_id` as `VARCHAR` (the `BIGINT`/`UUID` type is not preserved on disk — recover it with e.g. `read_id::BIGINT`).
- `sequence1` (VARCHAR): DNA/RNA/protein sequence

**Optional columns:**

| Column | Type | Description |
| --- | --- | --- |
| `comment` | VARCHAR | Comment line (only included if `INCLUDE_COMMENT=true`) |
| `sequence_index` | BIGINT | Used as identifier if `ID_AS_SEQUENCE_INDEX=true` |
| `sequence2` | VARCHAR | Second read for paired-end data |

**COPY options**

| Option | Default | Purpose |
| --- | --- | --- |
| `INCLUDE_COMMENT` | `false` | Include the comment from the input. |
| `ID_AS_SEQUENCE_INDEX` | `false` | Use `sequence_index` as the identifier. |
| `INTERLEAVE` | `false` | Write paired reads into one file. Incompatible with `{ORIENTATION}`. |
| `COMPRESSION` | Auto | Infer gzip from the `.gz` extension or set explicitly. |
| `COMPRESSION_LEVEL` | `1` | Set the local gzip compression level from `0` to `9`. |
| `ROW_GROUP_SIZE` | `64` | Rows in each batch DuckDB hands the writer. See "Threads, order and memory" below. |
| `ROW_GROUP_SIZE_BYTES` | `'8MB'` | Bytes in each batch, as DuckDB measures the rows. See "Threads, order and memory" below. |
| `LINE_WIDTH` | One line per sequence | Wrap sequences at this many characters. The minimum is `1`. |

**Compression**

Local gzip output uses BGZF (block gzip via HTSlib), which standard gzip tools can read. It uses DuckDB's thread count, does not change the record order (see "Threads, order and memory" below), supports block indexing, and runs about 10× faster than single-threaded gzip on multi-core hosts. Lower compression levels are faster and produce larger files. Remote targets use single-threaded gzip. `COMPRESSION_LEVEL` is available only for local gzip output.

`LINE_WIDTH` wraps the sequence only. The last line may be shorter, and an exact multiple of the width adds no empty line.

**Single-end vs paired-end output:**

FASTA paired-end output works the same way as FASTQ (a record is paired when `sequence2` is non-NULL), without the quality columns:

- **`{ORIENTATION}` in the path** writes split R1/R2 files; an entirely single-end dataset writes only the R1 file.
- **No `{ORIENTATION}`** writes a single file; paired data is interleaved only if `INTERLEAVE=true`.
- The writer fails loud when paired data is sent to a single-file path without `INTERLEAVE=true`, when `{ORIENTATION}` is combined with `INTERLEAVE=true`, or when the input mixes single-end and paired-end records.

An error about a record's data (a NULL `sequence1`, or with `ID_AS_SEQUENCE_INDEX` a NULL `sequence_index`; a paired record with nowhere to put R2) names its `read_id`.

**Threads, order and memory:** as for [FASTQ](#fastq), with files named `data_0.fasta` (`data_0.fasta.gz` with gzip), and `FORMAT _miint_fasta` as the writer behind `FORMAT FASTA`.

**More examples**

```sql
-- Paired-end split files (the {ORIENTATION} placeholder selects split output)
COPY (SELECT * FROM read_fastx('R1.fasta', sequence2 = 'R2.fasta'))
TO 'output_{ORIENTATION}.fasta' (FORMAT FASTA);

-- Compressed with comments
COPY (SELECT * FROM read_fastx('input.fasta'))
TO 'output.fasta.gz' (FORMAT FASTA, INCLUDE_COMMENT true);

-- Lines of 60
COPY (SELECT * FROM read_fastx('input.fasta'))
TO 'output.fasta' (FORMAT FASTA, LINE_WIDTH 60);

-- Genomes on all threads, in no fixed order
COPY (SELECT read_id, sequence1 FROM genomes)
TO 'genomes.fasta.gz' (FORMAT FASTA);

-- Genomes sorted by name
COPY (SELECT read_id, sequence1 FROM genomes ORDER BY read_id)
TO 'genomes_sorted.fasta.gz' (FORMAT FASTA);
```

Prodigal's `-a` (proteins) and `-d` (gene nucleotides) FASTA, from `predict_genes`, are written with this format and the [`prodigal_proteins` and `prodigal_genes`](gene_prediction.md#writing-prodigals-fasta) macros.
<!-- END REFERENCE: copy_fasta -->

###  SAM and BAM

Write query results to SAM or BAM format files. Requires all mandatory SAM columns from `read_alignments` output.

**Required columns:**

| Column | Type | Description |
| --- | --- | --- |
| `read_id` | VARCHAR, BIGINT, or UUID | Query template name |
| `flags` | USMALLINT | Bitwise flags |
| `reference` | VARCHAR, BIGINT, or UUID | Reference sequence name |
| `position` | BIGINT | 1-based leftmost mapping position |
| `mapq` | UTINYINT | Mapping quality |
| `cigar` | VARCHAR | CIGAR string |
| `mate_reference` | VARCHAR, BIGINT, or UUID | Reference name of mate/next read |
| `mate_position` | BIGINT | Position of mate/next read |
| `template_length` | BIGINT | Observed template length |

**Identifier-column types (`read_id`, `reference`, `mate_reference`):**

- Each column is independently `VARCHAR`, `BIGINT`, or `UUID`. The three need not match. Other numeric types are rejected at bind time with the message `Column '<name>' must be VARCHAR, BIGINT, or UUID`.
- `BIGINT` values are stringified to decimal on write; `UUID` values are stringified to their canonical 36-char lowercase form. The on-disk SAM/BAM record always carries text.
- A NULL in any of these three columns is written as the SAM `*` sentinel, whatever the column's type. So a single-end record can carry `NULL` (or `'*'`) for `mate_reference` interchangeably.
- Round-trip through `read_alignments` / `read_sam` returns these columns as `VARCHAR` regardless of how they were written; the `BIGINT`/`UUID` type is not preserved on disk.
- HTSlib normalises `RNEXT` to `=` when it equals `RNAME` at the reference-tid level. This applies to all id types — a `BIGINT`/`UUID` mate_reference that matches `reference` still appears as `=` when the file is read back.

**NULL handling in the required columns:** a NULL is rejected in `flags`, `mate_position` and `template_length` (`Column '<name>' is NULL for read: <id>`), because it has no SAM spelling and no defensible default. `read_id`, `reference` and `mate_reference` accept NULL and write the `*` sentinel, as above. `position`, `mapq` and `cigar` depend on the unmapped flag — see below. (`stop_position` is not consulted by the writer at all.)

**`cigar`, `position`, `mapq`:** whether a record may omit these depends on its unmapped flag (`0x4`):

- **Unmapped** (`0x4` set): none of them is expected. `cigar` may be `'*'`, `''` or NULL and is written as `*`; `position` and `mapq` may be NULL and are written as 0. Nothing can be lost, since there is no alignment to describe. This is what lets `align_minimap2(…, include_unmapped := true)` output — which carries NULL in exactly these columns — be written directly.
- **Mapped** (`0x4` clear): a CIGAR is required, and `position`/`mapq` must be present. A missing CIGAR (`'*'`, `''` or NULL) is rejected with `Mapped read '<id>' has no CIGAR (flags=<n>)`.

The mapped case is refused rather than written because htslib cannot represent it: samtools reads such a record back and silently rewrites it as unmapped (`[W::sam_parse1] mapped query must have a CIGAR; treated as unmapped`, `FLAG` 0 → 4). A CIGAR dropped by a bad join would otherwise become a mapped record with no alignment, and the COPY would report success. If you genuinely mean "did not align", set the `0x4` flag.

`read_alignments` / `read_sam` return `''` for a record whose on-disk CIGAR is `*`, so both spellings occur naturally and `COPY (SELECT * FROM read_alignments(...))` round-trips unmapped reads unchanged. One consequence worth knowing: a **mapped** record with no CIGAR — which older BAM files can contain, since `bam_read1` accepts `n_cigar == 0` without the rewrite `sam_parse1` applies — was previously copied through and is now rejected. Set `0x4` on such records, or filter them.

BAM stores POS, PNEXT and TLEN in 32 bits, so a record past that range (or with a CIGAR of more than 65,535 operations spanning 2^28 bases or more) cannot be written as BAM and is refused with an error naming the column, its value and the read, such as `Column 'position' is 2147483649 for read 'r5', past 2147483648, the largest a BAM record holds (SAM text holds it)`; SAM text writes it. A `read_id` that is empty (htslib would write it as `*`) or longer than 254 bytes is refused in SAM and BAM alike, as for [UBAM](#ubam). A `SEQUENCE_DATA` read whose qualities do not match its bases is refused as FORMAT FASTQ words it, with the table named.

A data error fails the COPY. The writer raises it before it writes any of the batch the record is in, but earlier batches may be written by then; DuckDB then removes the file the COPY created (each of them, with `PARTITION_BY` or `PER_THREAD_OUTPUT`), and a file that was already there keeps its old contents. A failed write, as on a full disk, does the same, with an error naming the format and the file, such as `COPY FORMAT BAM: cannot write 'out.bam': sam_close failed: File too large` (the end of a file is written when it is closed). Over an existing file, DuckDB writes to a temporary file beside it, `tmp_out.bam`, and the error names that.

**Optional columns:**

- `tag_as`, `tag_xs`, `tag_ys`, `tag_xn`, `tag_xm`, `tag_xo`, `tag_xg`, `tag_nm` (BIGINT): Optional integer tags
- `tag_yt`, `tag_md`, `tag_sa` (VARCHAR): Optional string tags

A tag column of another type is refused before anything is written (`Column 'tag_as' must be BIGINT`). A NULL writes no tag. An integer tag is a 32-bit signed value in SAM and BAM alike, so a value outside that range is an error naming the read (`Column 'tag_as' is 2147483648 for read 'r1', outside the 32-bit signed range of a SAM integer tag`), not a value cut to 32 bits.

SAM and BAM are written to local files only; a remote target (`s3://...`) is refused.

**Parameters:**

- `INCLUDE_HEADER` (default: true): Include header with reference sequences
    - **Note:** BAM format requires `INCLUDE_HEADER=true` (headers are mandatory in BAM files)
- `REFERENCE_LENGTHS` (VARCHAR, required if INCLUDE_HEADER=true): Table or view name containing reference sequences. Must have at least 2 columns: first column = reference name (VARCHAR, BIGINT or UUID, as `reference` takes, matched as its text (BIGINT in decimal, UUID in the lowercase hyphenated form)), second column = reference length (INTEGER/BIGINT). Column names don't matter. Views are fully supported and can include computed columns. Duplicate reference names are collapsed, keeping the first occurrence.

**`@SQ` order, and how to write a coordinate-sorted BAM:**

`@SQ` lines are emitted **sorted by reference name**, in the name column's own order (numeric for BIGINT names), regardless of the order of rows in `REFERENCE_LENGTHS`. This matters because a reference's index in the `@SQ` list *is* its TID, and a BAM is coordinate-sorted by **TID**, not by reference name — so the header order determines what "sorted" means.

Because the order is name-sorted, sorting your records the obvious way produces a genuinely coordinate-sorted file, as long as the name column and `reference` have the same type (a VARCHAR `reference` sorts `'10'` before `'2'`, a BIGINT name column puts 2 first) and VARCHAR names are compared with the default binary collation:

```sql
COPY (SELECT * FROM aln ORDER BY reference, position)
TO 'sorted.bam' (FORMAT BAM, REFERENCE_LENGTHS 'refs');
-- samtools index sorted.bam   -> succeeds
```

The `ORDER BY` is what sorts the file: without it the records come out in no fixed order, and `PRESERVE_ORDER false` drops it (see "Threads, order and memory" below).

Two consequences worth knowing:

- The order does **not** depend on how `REFERENCE_LENGTHS` was built (ASC, DESC or shuffled all give the same header), and it does not depend on `PRAGMA threads`. This is deliberate: DuckDB does not guarantee the row order of a parallel table scan, so the table's own row order could not be a reproducible contract.
- The header is readable back from SQL with `read_alignment_header('file.bam')`, which exposes `@SQ` as `(tid, reference, length)` — see [Reading](reading.md#read_alignment_headerpath). That is the only way to observe `tid`, and therefore the only way to check this ordering contract on a BAM.
- **Karyotype order is not expressible**: `'chr10'` sorts before `'chr2'`, so a human-genome BAM will be ordered `chr1, chr10, chr11, …, chr2` rather than `chr1, chr2, …, chr10`. The file is still internally consistent and indexable; it simply won't match a reference-order convention that isn't lexicographic. Sorting your records by `reference` keeps them aligned with whatever the header says.

With `INCLUDE_HEADER=false` no `@SQ` lines are written at all, so the ordering above does not apply; each record's RNAME is written verbatim from its `reference` value. A `REFERENCE_LENGTHS` given with it is still read and checked (an empty table, a NULL name or length, or a negative length fails the COPY), but changes nothing in the file. (This mode is SAM-only — BAM requires a header.)

> **`REFERENCE_LENGTHS` must be exhaustive when a header is written.** Supplying the table asserts it is the complete reference set, so a record whose `reference` (or `mate_reference`) is absent from it is rejected: `Column 'reference' is 'ctg000002' for read 'r2', which is not in the REFERENCE_LENGTHS table`. The header is serialized before the records are, so it cannot gain `@SQ` lines afterwards — a missing reference has no representation in the file being written, and earlier versions lost such records silently (SAM read the name back as `*`, BAM dropped the record). Derive the table from the same relation you are writing, e.g. `SELECT DISTINCT reference … FROM aln`. With `INCLUDE_HEADER=false` there are no `@SQ` lines to be exhaustive about, so any reference is accepted and written verbatim.

- `SEQUENCE_DATA` (VARCHAR, optional): Table or view name containing original read sequences from `read_fastx`. When provided, writes actual SEQ and QUAL fields into the output instead of `*`. See [Sequence Data](#sequence-data) below.
- `COMPRESSION` (default: auto): Enable gzip compression for SAM (auto-detected from `.gz` extension); values as for the other text formats, [above](#writing-files). BAM is always BGZF-compressed: it accepts `'gzip'` or `'gz'`, in any case, and refuses any other value, which it would ignore; `COMPRESSION_LEVEL 0` writes its BGZF blocks uncompressed.
- `COMPRESSION_LEVEL` (BAM, and gzip SAM): BGZF compression level 0-9 (default: 6). Higher = better compression, slower speed. Refused for uncompressed SAM, as for FASTQ, FASTA and NEWICK.
- `ROW_GROUP_SIZE`, `ROW_GROUP_SIZE_BYTES` (defaults: 64 and `'8MB'`): how DuckDB batches the records it hands the writer; see "Threads, order and memory" below.

A COPY reads `REFERENCE_LENGTHS` and `SEQUENCE_DATA` once, when it runs, so a COPY that writes several files gives every file the same `@SQ` lines and the same sequences.

**Threads, order and memory:**

- **No record order, by default.** Unless an `ORDER BY` sorts the rows written (as for [FASTQ](#fastq)), records are formatted (SAM) or built (BAM) on all of DuckDB's threads, and follow the header in an order that can change from run to run; a pair's records need not stay next to each other. A tool that needs a pair's records together, reading name-grouped input (`samtools fixmate`, for one), needs `ORDER BY read_id`.
- **A sorted file needs an `ORDER BY`.** `ORDER BY reference, position` writes a coordinate-sorted file ([above](#sam-and-bam)), on one thread (the sort still uses all threads). Alignments at the same position tie, so the file is the same whatever `threads` is set to only with columns that tell them apart, such as `ORDER BY reference, position, read_id, flags`.
- **`PRESERVE_ORDER true` keeps the order the rows come in**, without sorting: converting a coordinate-sorted BAM keeps it sorted only with this (or an `ORDER BY`). Records are written in the query's order, on one thread, and, as for [FASTQ](#fastq), the file is the same whatever `threads` is set to only when that order is fixed; work behind a sort, an aggregate or a join's build side still uses all threads.
- **`PRESERVE_ORDER false` drops an `ORDER BY`**, and a BAM meant to be sorted comes out unsorted. `SET preserve_insertion_order` makes no difference here either, as for [FASTQ](#fastq).
- BAM and gzip SAM compress their blocks on all threads whatever the order, which does not change the bytes.
- `PARTITION_BY` and `PER_THREAD_OUTPUT` always use all threads, and name their files `data_0.sam` (`data_0.sam.gz` with gzip) or `data_0.bam`.
- **Batches.** As for [FASTQ](#fastq): `ROW_GROUP_SIZE 64` and `ROW_GROUP_SIZE_BYTES '8MB'` are added unless the COPY sets them, and `FORMAT _miint_sam` and `FORMAT _miint_bam`, the writers behind `FORMAT SAM` and `FORMAT BAM`, add no defaults. `ROW_GROUP_SIZE_BYTES` measures the rows the query hands the writer, which do not hold the reads `SEQUENCE_DATA` adds, so a batch of long reads can make far more than 8 MB of records.
- **Memory.** The records of a batch, once built, sit outside DuckDB's `memory_limit` until written: up to 64 MiB a batch; the rest of a larger batch is built as it is written. The rows themselves are held by DuckDB and counted against it. `SEQUENCE_DATA` is held in memory whole ([below](#sequence-data)).

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

-- Convert SAM to BAM (records in no fixed order)
COPY (SELECT * FROM read_alignments('input.sam'))
TO 'output.bam' (FORMAT BAM, REFERENCE_LENGTHS 'ref_table');

-- Convert BAM to SAM, keeping the input's record order (a sorted input stays sorted)
COPY (SELECT * FROM read_alignments('input.bam'))
TO 'output.sam' (FORMAT SAM, REFERENCE_LENGTHS 'ref_table', PRESERVE_ORDER true);

-- Filter alignments and write to BAM
COPY (
  SELECT * FROM read_alignments('input.bam')
  WHERE mapq >= 30 AND alignment_is_primary(flags)
) TO 'filtered.bam' (FORMAT BAM, REFERENCE_LENGTHS 'ref_table');

-- A coordinate-sorted BAM, ready for samtools index
COPY (SELECT * FROM read_alignments('input.bam') WHERE mapq >= 30 ORDER BY reference, position)
TO 'sorted.bam' (FORMAT BAM, REFERENCE_LENGTHS 'ref_table');

-- Write BAM with sequence data (SEQ and QUAL populated from FASTQ)
CREATE TABLE sequences AS SELECT * FROM read_fastx('reads_R1.fastq', sequence2 = 'reads_R2.fastq');

COPY (SELECT * FROM read_alignments('input.bam'))
TO 'with_seq.bam' (FORMAT BAM, REFERENCE_LENGTHS 'ref_table', SEQUENCE_DATA 'sequences');
```

#### Sequence Data

The `SEQUENCE_DATA` parameter allows writing actual SEQ and QUAL fields into SAM/BAM output by looking up original read sequences from a table or view. Without this parameter, SEQ and QUAL are written as `*`.

This design avoids duplicating immutable sequence data across alignment records. The original FASTQ sequences are stored once and joined at write time.

**Required columns in the SEQUENCE_DATA table:**

- `read_id` (VARCHAR, BIGINT or UUID): Read identifier, matched as its text (BIGINT in decimal, UUID in the lowercase hyphenated form) against the alignment's read_id values, so a BIGINT or UUID read_id needs no cast
- `sequence1` (VARCHAR): Forward-strand DNA sequence for read 1
- `qual1` (LIST(UTINYINT)): Phred quality scores for read 1 (raw values 0-93, no ASCII offset)

**Optional columns (for paired-end):**

- `sequence2` (VARCHAR): Forward-strand DNA sequence for read 2
- `qual2` (LIST(UTINYINT)): Phred quality scores for read 2

These columns match the output schema of `read_fastx`, so the typical workflow is:

```sql
-- Load sequences
CREATE TABLE sequences AS SELECT * FROM read_fastx('R1.fastq', sequence2 = 'R2.fastq');

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

### UBAM

Write query results to an **unaligned reads BAM** (uBAM): reads with SEQ/QUAL and no reference alignment. This is the BAM counterpart of `FORMAT FASTQ` — it takes the same `read_id`/`sequence1`/`qual1` columns — and exists for tools that require a BAM (rather than FASTQ) input for unaligned reads, e.g. PacBio's `lima`, which selects CCS-vs-CLR demultiplexing from the input **format**.

Every record is written unmapped (`FLAG=4`), with no `@SQ` header and no reference table. `read_id` is the record name (QNAME), written unchanged: the writer makes no assumption about its contents, but refuses an empty one and one longer than 254 bytes, which a BAM record cannot hold. Output is always BAM (block-gzip).

**Required columns:**

| Column | Type | Description |
| --- | --- | --- |
| `read_id` | VARCHAR, BIGINT, or UUID | Record name (QNAME), written unchanged; empty or longer than 254 bytes is an error. `BIGINT`/`UUID` are stringified as for FASTQ. |
| `sequence1` | VARCHAR | Read sequence. |
| `qual1` | UTINYINT[] | Phred quality scores, raw values 0–93 (no ASCII offset) — exactly what `read_fastx` / `read_sequences_sam` emit. Must be the same length as `sequence1`; a score past 93 is refused. |

Single-end only (no `sequence2`/`qual2`). Any additional projected columns are ignored unless referenced by `TAGS`.

UBAM is written to local files only; a remote target is refused.

**Parameters:**

- `READ_GROUP` (STRUCT, optional): Emits an `@RG` header line and a per-record `RG:Z` tag. An `ID` field is required; any other fields (`PL`, `DS`, `SM`, `LB`, `PU`, `CN`, …) are passed through. Field names are case-insensitive (upper-cased to the SAM codes). Example: `READ_GROUP {ID: 'qiita', PL: 'PACBIO', DS: 'READTYPE=CCS'}`. `DS:READTYPE=CCS` is the field `lima` keys on to treat the input as CCS.
- `TAGS` (STRUCT, optional): Maps a 2-character BAM aux tag to a projected **integer** column, written as a signed 32-bit (`i`) tag per record. A NULL cell omits the tag for that record. Example: `TAGS {zm: 'zmw'}` writes `zm:i:<zmw>` (PacBio ZMW hole number). Non-integer columns are rejected at bind. Note: struct keys are lower-cased unless quoted, so an uppercase tag needs a quoted key (e.g. `TAGS {'NM': 'edit_distance'}`), and two tags whose names differ only in case, such as `XU` and `xu`, cannot both be given (DuckDB refuses a struct with both keys); the `RG` tag is reserved for `READ_GROUP`.
- `COMPRESSION_LEVEL` (0–9, optional): BAM compression level (default: HTSlib's level 6).
- `ROW_GROUP_SIZE`, `ROW_GROUP_SIZE_BYTES` (defaults: 64 and `'8MB'`): how DuckDB batches the reads it hands the writer; see "Threads, order and memory" below.

`REFERENCE_LENGTHS` is **not** accepted — a uBAM's header has no `@SQ` lines, as its reads are unaligned. Use `FORMAT BAM` for aligned records.

An error about a record's data names its `read_id`, as FORMAT FASTQ's do: a NULL `sequence1` or `qual1` (a NULL `read_id` has no name to give); a NULL quality score; a quality length that does not match the sequence; an empty `sequence1` (a record without bases reads as one whose sequence was left out); a score past 93, the highest SAM text can show (a first score of 255 would mark the record as having none); a `read_id` longer than 254 bytes, the longest a BAM record holds (`Column 'read_id' is 300 bytes for read '…', past 254, the longest a BAM record holds`), or empty, which a BAM record cannot hold (htslib would write `*`); and a `TAGS` value outside the 32-bit signed range of a SAM integer tag, from any integer type `TAGS` takes (`Column 'zmw' (tag 'zm') is 2147483648 for read 'm1/7/ccs', outside the 32-bit signed range of a SAM integer tag`). A data error or a failed write fails the COPY as for [SAM and BAM](#sam-and-bam): DuckDB removes the file the COPY created, and a file that was already there keeps its old contents. A failed write names the format and the file, such as `COPY FORMAT UBAM: cannot write 'out.bam': sam_close failed: File too large`.

**Threads, order and memory:**

- **No record order, by default.** Unless an `ORDER BY` sorts the rows written (as for [FASTQ](#fastq)), records are built on all of DuckDB's threads and come out in an order that can change from run to run.
- **Keep a ZMW's reads together when a tool reads them as a group**, as `ccs` does with a subreads BAM: sort by the ZMW column (`ORDER BY zmw`), or, when the input is already grouped, keep its order with `PRESERVE_ORDER true`. A ZMW's reads tie on `ORDER BY zmw`, so the file is the same whatever `threads` is set to only with a column that tells them apart, such as `ORDER BY zmw, read_id`.
- **`PRESERVE_ORDER true` keeps the order the rows come in**, without sorting, on one thread, and, as for [FASTQ](#fastq), the file is the same whatever `threads` is set to only when that order is fixed; work behind a sort, an aggregate or a join's build side still uses all threads.
- **`PRESERVE_ORDER false` drops an `ORDER BY`.** `SET preserve_insertion_order` makes no difference here either, as for [FASTQ](#fastq).
- BGZF blocks are compressed on all threads whatever the order, which does not change the bytes.
- `PARTITION_BY` and `PER_THREAD_OUTPUT` always use all threads, and name their files `data_0.bam`; `EXPORT DATABASE … (FORMAT UBAM)` names each table's file `<table>.bam`. A `.bam` path with no `FORMAT` is written as aligned BAM (`FORMAT BAM`), never as uBAM.
- **Batches.** As for [FASTQ](#fastq): `ROW_GROUP_SIZE 64` and `ROW_GROUP_SIZE_BYTES '8MB'` are added unless the COPY sets them, and `FORMAT _miint_ubam`, the writer behind `FORMAT UBAM`, adds no defaults. `ROW_GROUP_SIZE_BYTES` measures the rows, about 2 bytes a base (the sequence and its qualities); a record holds about 1.5.
- **Memory.** As for [SAM and BAM](#sam-and-bam): the records of a batch, once built, sit outside DuckDB's `memory_limit` until written, up to 64 MiB a batch; the rest of a larger batch is built as it is written.

**Examples:**
```sql
-- CCS reads from the lake -> a uBAM lima will accept. zmw is a projected
-- integer column carrying the ZMW hole number for the zm tag.
COPY (SELECT read_id, sequence1, qual1, zmw FROM reads)
TO 'out.bam' (
  FORMAT UBAM,
  READ_GROUP {ID: 'qiita', PL: 'PACBIO', DS: 'READTYPE=CCS'},
  TAGS {zm: 'zmw'}
);

-- Minimal uBAM (no @RG, no tags), records in no fixed order.
COPY (SELECT read_id, sequence1, qual1 FROM read_fastx('input.fastq'))
TO 'reads.bam' (FORMAT UBAM);

-- Each ZMW's reads together, the same file at any thread count (read_id is
-- unique within a ZMW).
COPY (SELECT read_id, sequence1, qual1, zmw FROM reads ORDER BY zmw, read_id)
TO 'by_zmw.bam' (FORMAT UBAM, READ_GROUP {ID: 'qiita', PL: 'PACBIO'}, TAGS {zm: 'zmw'});
```

Round-trip via `read_sequences_sam('out.bam')` returns `read_id`/`sequence1`/`qual1` (with `read_id` as VARCHAR, as for the FASTQ writer).

### BIOM

Write query results to BIOM (Biological Observation Matrix) format files. BIOM is an HDF5-based format commonly used for representing OGU/OTU/ASV tables in microbiome analyses.

**Required columns:**

| Column | Type | Description |
| --- | --- | --- |
| `feature_id` | VARCHAR, BIGINT or UUID | Feature/OGU/OTU/ASV identifier |
| `sample_id` | VARCHAR, BIGINT or UUID | Sample identifier |
| `value` | DOUBLE | Abundance/count value |

An id is written as its text, as `read_id` is elsewhere: BIGINT in decimal, UUID in the lowercase hyphenated form. `read_biom` returns ids as VARCHAR. A column of another type is refused when the COPY is planned, with an error naming the column. An untyped decimal literal such as `45.0` is DECIMAL, so cast a literal `value` to DOUBLE.

BIOM is written to local files only; a remote target is refused. `PARTITION_BY` and `PER_THREAD_OUTPUT` write `.biom` files, each a complete table of its own rows: one per thread, and one or more per partition (DuckDB starts another file for a partition when more partitions are open than `partitioned_write_max_open_files`), so a (feature, sample) pair can be split across files. DuckDB leaves the partition columns out of the rows it writes unless `WRITE_PARTITION_COLUMNS true`, so partitioning by `sample_id` or `feature_id` needs it.

A target file that already exists is refused, and so is one made between `PREPARE` and `EXECUTE`; the file is left as it was. The error is `Cannot overwrite existing file`, or `H5Fcreate failed: File exists` when DuckDB does not bind the COPY again at `EXECUTE`. `USE_TMP_FILE true`, which has DuckDB rename a temporary file over the target, is refused when the COPY runs. With `PARTITION_BY` or `PER_THREAD_OUTPUT` and `OVERWRITE_OR_IGNORE`, a file of the same name already in the directory is replaced, as that option asks. A write that fails partway, on a full disk for example, fails the COPY with an IO Error naming the file, the HDF5 call that failed and the system's reason, such as `No space left on device`, and DuckDB removes the file the failed COPY created.

**Parameters:**

- `COMPRESSION` (default: 'gzip'): HDF5's compression of the datasets inside the file, `'gzip'` or `'gz'` (deflate level 4, h5py's default), or `'none'` or `'uncompressed'`, in any case. Any other value, zstd included, is refused. The `.biom` file itself is always HDF5, never gzip, whatever its name.
- `ID` (default: 'No Table ID'): Custom table identifier for BIOM metadata
- `GENERATED_BY` (default: 'miint'): Tool/version string for provenance tracking

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
    ('OGU_001', 'Sample_A', 45.0::DOUBLE),
    ('OGU_001', 'Sample_B', 32.0::DOUBLE),
    ('OGU_002', 'Sample_A', 18.0::DOUBLE),
    ('OGU_002', 'Sample_B', 27.0::DOUBLE)
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
    ('F1', 'S1', 10.0::DOUBLE),
    ('F1', 'S1', 5.0::DOUBLE),   -- Same feature + sample
    ('F2', 'S1', 8.0::DOUBLE)
) AS t(feature_id, sample_id, value);

COPY (SELECT * FROM with_dups)
TO 'merged.biom' (FORMAT BIOM);

-- Result: F1,S1,15.0 and F2,S1,8.0 (duplicates merged)

-- Zeros are removed
CREATE TABLE with_zeros AS
SELECT * FROM (VALUES
    ('F1', 'S1', 10.0::DOUBLE),
    ('F1', 'S2', 0.0::DOUBLE),   -- Zero value
    ('F2', 'S1', 8.0::DOUBLE)
) AS t(feature_id, sample_id, value);

COPY (SELECT * FROM with_zeros)
TO 'sparse.biom' (FORMAT BIOM);

-- Result: Only F1,S1,10.0 and F2,S1,8.0 (zero removed)
```

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

## Newick

Write query results to Newick phylogenetic tree format. Reconstructs a tree from tabular node data and serializes to standard Newick format.

**Required columns:**

| Column | Type | Description |
| --- | --- | --- |
| `node_index` | BIGINT | Unique identifier for each node |
| `parent_index` | BIGINT, nullable | Parent node's node_index (NULL for root) |

**Optional columns:**

| Column | Type | Description |
| --- | --- | --- |
| `name` | VARCHAR, BIGINT or UUID | Node name, an id such as a tip's feature id, written as its text |
| `label` | VARCHAR | Text label, written where the node has no `name` |
| `branch_length` | DOUBLE | Branch length |
| `edge_id` | BIGINT | Edge identifier for jplace format |

Each column takes only the types listed. A column of another type, such as the INTEGER and DECIMAL that untyped literals like `1` and `0.1` are, is refused when the COPY is planned; cast it in the query (`node_index::BIGINT`). `name` takes the types a `read_id` does. An untyped `NULL` column is cast to the column's type, so `NULL AS name` writes no names.

A node's Newick label is its `name`, or its `label` where the name is NULL or `''`, so a tree with BIGINT or UUID tip ids and text labels on its internal nodes (see "Tree tables" under [Newick](reading.md#newick)) is written whole. A `label` of `''` is no label. Newick holds one label per node, so a node with both a non-empty `name` and a non-empty `label` is an error naming its `node_index`. `read_newick` reads every label back into `name`.

**Parameters:**

- `EDGE_IDS` (BOOLEAN, default: auto): Include edge identifiers `{n}` in output
    - Default: true if `edge_id` column exists, false otherwise
    - Set to `false` to explicitly exclude edge IDs
- `COMPRESSION` (default: auto): Enable gzip compression, auto-detected from the `.gz` extension; values as for the other text formats, [above](#writing-files)
- `COMPRESSION_LEVEL` (0-9, default: 6): BGZF compression level for local `.gz` output, as for FASTQ and FASTA. Refused for uncompressed output and for remote targets, whose gzip writer takes no level.
- `PLACEMENTS` (VARCHAR): Name of a table containing phylogenetic placement data to insert
    - Inserts fragment sequences into the tree at their placement locations
    - Requires tree to have `edge_id` column with valid edge identifiers
    - Required table columns: `fragment_id` (VARCHAR, BIGINT or UUID, as a `read_id`; the placed tip is named with its text), `edge_id` (BIGINT, INTEGER, UBIGINT or UINTEGER), `like_weight_ratio`, `distal_length` and `pendant_length` (each DOUBLE or FLOAT). A column of another type is refused, when the COPY is planned and again when the placements are read
    - Duplicate fragment_ids are deduplicated as in [`tree_resolve_placement`](phylogeny.md#resolve-placements), which does not always keep the placement with the highest `like_weight_ratio`
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
- Writes one tree to one file: `PER_THREAD_OUTPUT`, `PARTITION_BY` and `ORDER_BY` (which DuckDB writes as a directory of files) are refused, as each file would hold only part of the tree. The error comes when the COPY runs, not when it is planned (`EXPLAIN` and `PREPARE` do not raise it), and a partitioned or ordered COPY of no rows succeeds and writes no file

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

-- Create tree from scratch (the casts give each column its type)
CREATE TABLE my_tree AS
SELECT * FROM (VALUES
    (0::BIGINT, NULL::BIGINT, '', NULL::DOUBLE),  -- root
    (1, 0, 'A', 0.1),                              -- tip A
    (2, 0, 'B', 0.2)                               -- tip B
) AS t(node_index, parent_index, name, branch_length);

COPY (SELECT * FROM my_tree)
TO 'new_tree.nwk' (FORMAT NEWICK);
-- Produces: (A:0.1,B:0.2);

-- BIGINT tip ids, with a text label on the internal node
CREATE TABLE labeled_tree AS
SELECT * FROM (VALUES
    (0::BIGINT, NULL::BIGINT, NULL::BIGINT, 'g__Escherichia'),  -- root
    (1, 0, 101, NULL),                                          -- tip 101
    (2, 0, 102, NULL)                                           -- tip 102
) AS t(node_index, parent_index, name, label);

COPY (SELECT * FROM labeled_tree)
TO 'labeled_tree.nwk' (FORMAT NEWICK);
-- Produces: (101,102)g__Escherichia;

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
- "node N has both a name ('...') and a label ('...')" - A node with a non-empty `name` and a non-empty `label`

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

## GFF

`COPY ... TO 'out.gff' (FORMAT GFF)` writes GFF, a line per row. A path ending in `.gff` is FORMAT GFF without saying so. The rows are made into lines in SQL and DuckDB's CSV writer writes them (as DuckDB's own `FORMAT JSON` does), so `COMPRESSION` (and gzip for a path ending in `.gz`), `PARTITION_BY`, `FILE_EXTENSION`, `OVERWRITE`, `OVERWRITE_OR_IGNORE`, `FILENAME_PATTERN`, `USE_TMP_FILE` and `RETURN_FILES` are the CSV writer's. Options that would break the file (`APPEND`, `PER_THREAD_OUTPUT`, `FILE_SIZE_BYTES`, `WRITE_PARTITION_COLUMNS`, `PRESERVE_ORDER`, and CSV's own) are refused, and so is `PARTITION_PATH`. The lines come out in the query's order, and the count COPY returns is the lines written, the version and header lines included.

**Parameters:**

- `DIALECT` (VARCHAR, default `'gff3'`): `'gff3'` or `'prodigal'`.
- `SOURCE` (VARCHAR, default `'Prodigal_v2.6.5'`, `'prodigal'` only): column 2.
- `PARTITION_BY (genome_id)` (`'prodigal'` only): a file per genome.

- **DIALECT `'gff3'`**, from **[`read_gff`](reading.md#gff)'s columns**, `seqid`, `source`, `type` (VARCHAR), `position`, `stop_position` (integers), `score` (a number), `strand` (VARCHAR), `phase` (an integer) and `attributes` (`MAP(VARCHAR, VARCHAR)`, or VARCHAR written as given): GFF3, a `##gff-version 3` line, then a line per row. What `read_gff` reads, this writes and `read_gff` reads back the same. A missing column or a column of another type is refused as the query is bound, even when every value is NULL.
    - `stop_position` is half-open, as `read_gff` returns it, so the GFF `end` written is `stop_position - 1`. `position` must be at least 1 and `stop_position` greater than `position`.
    - A NULL `score`, `strand`, `phase` or `attributes` is `.`. `strand` is `+`, `-`, `.` or `?`; `phase` 0, 1 or 2.
    - A MAP's keys and values are percent-encoded where GFF3 reserves the character (`;`, `=`, `&`, `,`, `%`, tab, line breaks and other control characters), and `read_gff` decodes them. A key with a NULL value is written alone, without `=`, as `read_gff` reads such a token (an attributes column ending in `;`, or `.`, reads back the same). Text `attributes` is written as given and may not hold a control character (a tab or line break, say).
    - A value is one value: `read_gff` reads `Parent=a,b` as the text `a,b` (it does not split GFF3's comma-separated lists), so this writes it back as `Parent=a%2Cb`, which a reader that splits lists takes as one parent. To write a list, give `attributes` as text.
    - `seqid`, `source` and `type` are written as given (`read_gff` does not decode them), so a NULL or empty one, a control character, or a `seqid` starting with `#` (which would make the line a comment) is refused.
- **DIALECT `'prodigal'`**, from **[`predict_genes`](gene_prediction.md)' columns**: Prodigal's own GFF, as `prodigal -f gff` writes it (but for column 6 in metagenomic mode, below), with a `# Sequence Data` and a `# Model Data` line before each contig's genes. CheckM reads the translation table from the Model Data line. Its columns' types are checked as the query is bound (`read_id` may be any type). Without DIALECT, predict_genes' columns (with `gene_index` and without read_gff's) are refused with a pointer to it.
    - The rows must be in contig order, as `predict_genes` returns them: the header lines are written each time the contig (`genome_id` and `read_id`) changes. Without a `genome_id` column, the rows are one genome's.
    - A contig without genes (its `predict_genes` row has a NULL `gene_index`) gets its header lines and no feature, as in the CLI's file, and is counted in `seqnum`, which numbers the contigs from 1. A genome's own row (`training_failed`, or `no_genes` with a NULL `read_id`) is skipped.
    - `seqhdr` is the `read_id`. A value a line needs that is NULL, a `read_id` with a control character or starting with `#`, and a `strand` other than `+` or `-` are refused.
    - Column 6 is `score` (column 9's `score=`) to one decimal, which is what the CLI writes there in single-genome mode. In metagenomic mode the CLI writes the gene's score (`cscore + sscore`) computed a second time under the chosen model, which can differ from `score=` and which `predict_genes` does not return, so there column 6 can differ from the CLI's. CheckM does not read column 6 (`checkm/prodigal.py`, 1.2.5 and 1.3.0).
    - With `PARTITION_BY (genome_id)`, each genome is a file of its own, with its own version line and `seqnum` from 1.

```sql
-- Round trip
CREATE TABLE features AS SELECT * FROM read_gff('in.gff');
COPY features TO 'out.gff.gz' (FORMAT GFF);

-- Prodigal's GFF of every genome's genes in one file, contigs in file order (predict_genes' own order)
CREATE TABLE genes AS SELECT * FROM predict_genes('contigs', genome_id := 'filepath');
COPY genes TO 'genes.gff' (FORMAT GFF, DIALECT 'prodigal');

-- A file per genome, each as the Prodigal CLI would write it for that genome alone
COPY genes TO 'gff_by_bin' (FORMAT GFF, DIALECT 'prodigal', PARTITION_BY (genome_id));
```
