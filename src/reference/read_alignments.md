### SAM and BAM

<!-- manual:description -->
Read SAM and BAM alignment records. Use `include_seq_qual` to return sequence and quality scores.

<!-- manual:example -->
```sql
SELECT * FROM 'alignments.sam';
SELECT * FROM 'alignments.bam';
```

DuckDB uses `read_alignments` for these [file suffixes](#reading-a-file-by-its-name). Use the function to pass options, read stdin, or open a file with a different suffix.

#### Read sequence data or track source files

<!-- manual:example -->
```sql
SELECT read_id, sequence, qual FROM read_alignments('alignments.sam', include_seq_qual = true);
```

```sql
SELECT * FROM read_alignments('samples/*.bam', include_filepath = true);
SELECT * FROM read_alignments(['file1.sam', 'file2.bam'], include_filepath = true);
```

Use [`read_sequences_sam`](#read_sequences_samfilename) when you want sequence rows compatible with `read_fastx`. Use [`read_alignment_header`](#read_alignment_headerpath) to inspect the reference dictionary.

#### SAM and BAM parameters

| Parameter | Type | Default | Purpose |
| --- | --- | --- | --- |
| `filename` | `VARCHAR` or `VARCHAR[]` | Required | Path, glob, list of literal paths, or `-` / `/dev/stdin` for standard input. |
| `reference_lengths` | `VARCHAR` | Omitted | Name of a table or view supplying reference names and lengths for headerless SAM. |
| `include_filepath` | `BOOLEAN` | `false` | Include the source path. |
| `include_seq_qual` | `BOOLEAN` | `false` | Include the SAM sequence and Phred scores. |

#### SAM and BAM output columns

| Column | Type | Meaning |
| --- | --- | --- |
| `read_id` | `VARCHAR` | Query name (QNAME). |
| `flags` | `USMALLINT` | SAM bit flags (FLAG). |
| `reference` | `VARCHAR` | Reference name (RNAME), or `*` if unavailable. |
| `position` | `BIGINT` | 1-based start (POS), or 0 if unavailable. |
| `stop_position` | `BIGINT` | 1-based exclusive end from the CIGAR, or 0 for an unmapped record. See [Coordinate conventions](#coordinate-conventions). |
| `mapq` | `UTINYINT` | Mapping quality (MAPQ). |
| `cigar` | `VARCHAR` | CIGAR string, or an empty string when the SAM field is `*`. |
| `mate_reference` | `VARCHAR` | Mate reference (RNEXT), `=` for the same reference, or `*` if unavailable. |
| `mate_position` | `BIGINT` | 1-based mate start (PNEXT), or 0 if unavailable. |
| `template_length` | `BIGINT` | Signed template length (TLEN). |
| `tag_as`, `tag_xs`, `tag_ys`, `tag_xn`, `tag_xm`, `tag_xo`, `tag_xg`, `tag_nm` | `BIGINT` | The corresponding integer SAM tags, or `NULL` when absent. |
| `tag_yt`, `tag_md`, `tag_sa` | `VARCHAR` | The corresponding text SAM tags, or `NULL` when absent. |
| `sequence` | `VARCHAR` | SEQ, returned with `include_seq_qual = true`. |
| `qual` | `UTINYINT[]` | Phred scores (0–93), returned with `include_seq_qual = true`. Missing qualities produce an empty list. |
| `filepath` | `VARCHAR` | Source path, returned with `include_filepath = true`. |

#### Paths, compression, and parallel reads

The reader detects SAM or BAM from the file contents and supports gzip-compressed SAM. A glob's matches are sorted alphabetically. List entries are literal paths. Files can be read in parallel, with one thread per file up to DuckDB's `threads` setting. `read_sam` is a compatibility alias for `read_alignments`.

#### Headerless SAM and standard input

For headerless SAM, pass a reference table or view. The first column contains reference names as `VARCHAR`, `BIGINT`, or `UUID`. Matching uses the text representation, including decimal integers and lowercase hyphenated UUIDs. The second column contains lengths as `INTEGER` or `BIGINT`. Column names are arbitrary and views may compute the columns. When reading multiple local files, all inputs must agree on whether they contain reference dictionary entries (`@SQ`).

```sql
CREATE TABLE my_refs AS
  SELECT 'chr1' AS name, 248956422 AS length
  UNION ALL SELECT 'chr2', 242193529;
SELECT * FROM read_alignments('headerless.sam', reference_lengths = 'my_refs');
```

Read stdin with `-` or `/dev/stdin` as the only input. It uses one thread. Apply the same header and reference-table rules to stdin.

```sql
SELECT * FROM read_alignments('/dev/stdin');
```

Use [`read_sequences_sam`](#read_sequences_samfilename) to read sequences from an unaligned BAM with no reference dictionary.

#### Missing sequence and malformed records

With `include_seq_qual = true`, primary and unmapped records must carry SEQ. On reverse-strand alignments, SEQ is the sequence as stored in SAM, reverse-complemented relative to the original read. Hard-clipped bases are absent.

A record htslib cannot parse, such as a malformed field, truncated record, or out-of-range reference ID, raises `Failed to read SAM/BAM record: the file is truncated or malformed`. The error includes the line number for text SAM. A BGZF file truncated at a block boundary can look like a clean end of file to htslib. `samtools quickcheck <file>` can detect a missing BGZF EOF marker.

<!-- manual:url /reading/#sam-and-bam -->
