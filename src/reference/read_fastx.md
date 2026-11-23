### FASTA and FASTQ

<!-- manual:url /reading/#fasta-and-fastq -->
<!-- manual:description -->
Read FASTA and FASTQ files, including paired-end reads. Pass a path, a glob, or a list of literal paths.

After [loading miint](#reading-a-file-by-its-name), read a file directly by name. DuckDB calls `read_fastx` for these suffixes.

<!-- manual:example -->
```sql
SELECT * FROM 'sequences.fasta';
SELECT * FROM 'reads.fastq.gz';
```

Use the [file suffix table](#reading-a-file-by-its-name) to see which names DuckDB recognizes. Call `read_fastx` when you need to pass options, read from stdin, or use a path with an unrecognized suffix.

#### Read paired FASTQ or several files

<!-- manual:example -->
```sql
SELECT * FROM read_fastx('R1.fastq', sequence2 = 'R2.fastq');
```

```sql
SELECT * FROM 'samples/*.fastq';

SELECT * FROM read_fastx('samples/*.fastq', include_filepath = true)
ORDER BY filepath, sequence_index;
```

Reads are paired by position. Each pair must have matching read IDs after removing mate suffixes such as `/1` and `/2`.

#### Parameters

| Parameter | Type | Default | Purpose |
| --- | --- | --- | --- |
| `filename` | `VARCHAR` or `VARCHAR[]` | Required | Path, glob, or list of paths for a single-end input or R1. |
| `sequence2` | `VARCHAR` or `VARCHAR[]` | Omitted | R2 path, glob, or list of paths. |
| `include_filepath` | `BOOLEAN` | `false` | Include the source path. |
| `qual_offset` | `INTEGER` | `33` | Decode quality scores using Phred+33 or Phred+64. |
| `allow_empty` | `BOOLEAN` | `false` | Accept empty files as zero records. A pair with exactly one empty mate fails. |
| `max_batch_bytes` | `VARCHAR` | `'64MiB'` | Set the soft limit on the bytes of the records produced at a time. See [Memory and parallel reads](#memory-and-parallel-reads). |

#### Output columns

| Column | Type | Meaning |
| --- | --- | --- |
| `sequence_index` | `BIGINT` | 1-based record index, starting over for each file. |
| `read_id` | `VARCHAR` | Sequence ID without the leading `@` or `>`. |
| `comment` | `VARCHAR` | Comment after the ID, or `NULL` if absent. |
| `sequence1` | `VARCHAR` | Sequence, or R1 for paired reads. |
| `sequence2` | `VARCHAR` | R2 sequence, or `NULL` for single-end reads. |
| `qual1` | `UINT8[]` | Integer quality scores for R1, or `NULL` for FASTA. A FASTQ record with no bases has an empty list. |
| `qual2` | `UINT8[]` | Integer quality scores for R2, or `NULL` for FASTA or single-end reads. A FASTQ record with no bases has an empty list. |
| `filepath` | `VARCHAR` | Source path, included only with `include_filepath = true`. |

#### Paths and file pairing

A single path containing `*`, `?`, or `[` is expanded and sorted alphabetically. List elements are literal paths. For paired-end reads, give the same number of files on each side. If R1 is a glob, R2 must also be a glob. Each glob is expanded and sorted independently, then files are paired by position.

```sql
SELECT * FROM read_fastx(['sample1.fastq', 'sample2.fastq']);

SELECT * FROM read_fastx(
    ['sample1_R1.fastq', 'sample2_R1.fastq'],
    sequence2 = ['sample1_R2.fastq', 'sample2_R2.fastq']
);

SELECT * FROM read_fastx('samples/*_R1.fastq', sequence2 = 'samples/*_R2.fastq');
```

#### Formats, compression, and invalid input

Detects FASTA (`>` headers) or FASTQ (`@` headers) from the file's first header, whatever its extension. Every record of a FASTA file has a `NULL` `qual1`, and every record of a FASTQ file has a list, even a record with no bases. A file that mixes the two formats fails the query, naming its first record of the other format. Both files of a pair must use the same format. A FASTA file paired with a FASTQ file fails the query.

Reads gzip-compressed files. The reader detects gzip from content for local files and from the `.gz` suffix for remote files (`s3://`, `https://`).

Supports stdin input using `-` or `/dev/stdin` (single file only, no paired-end).

Each FASTQ record takes four lines: header, sequence, `+`, and quality. The quality and sequence must have the same length. Blank or whitespace-only lines before the first record and between records are ignored. FASTA sequences may span multiple lines, and text before the first FASTA header is skipped. Wrapped FASTQ records fail.

A malformed record, a file that ends inside a FASTQ record, a truncated or corrupt gzip file, or NUL bytes where a record should start (such as the zero-filled block a crash can leave) fail the query rather than returning the records before it. The error names the file, and either the malformed record or how many records were read before the read failed; for damaged gzip data that count can stop short of the damage, since the data is inflated in blocks.

#### Memory and parallel reads

The reader streams records from each file rather than loading it whole. Rows come out in groups of up to 2,048 records. `max_batch_bytes` is a positive formatted byte size such as `'256MiB'` or `'2GB'`. A group closes early once the uncompressed bytes of its read IDs, comments, sequences and qualities reach it. This bounds memory for large records such as assembled genomes. With short reads, a group reaches 2,048 records long before the cap.

The reader uses one thread per file or file pair, up to DuckDB's `threads` setting. It reads stdin on one thread.

Under the default `preserve_insertion_order`, rows come out file by file, in the order the files are given (a glob's sorted by name), and each file's records in order. `CREATE TABLE ... AS`, `INSERT`, and `COPY` to CSV or Parquet keep that order and still read the files in parallel. `COPY ... (FORMAT FASTQ)` and `(FORMAT FASTA)` do not keep it unless asked. By default they read and write on all threads, in no fixed order. With `PRESERVE_ORDER true` they keep it, writing on one thread and reading on it too. See [Writing](writing.md).

Reading several files in parallel with insertion order kept costs memory. Each file is read whole by one thread, and the rows of later files wait in memory until the files before them are done. For a few large files, that is most of the input. A consumer that keeps the order, such as `CREATE TABLE ... AS`, holds whole groups until the ones before them arrive, so its memory also follows `max_batch_bytes`. DuckDB bounds the wait by holding back the later files once it passes a share of `memory_limit`, and the read then slows toward one thread's speed. The [engine notes](internals/duckdb-engine-notes.md#output-order-costs-threads-and-memory) have the measurements. Where record order does not matter, `SET preserve_insertion_order = false` (or `PRESERVE_ORDER false` on a `COPY`) avoids the wait, as a FASTQ or FASTA `COPY` with no `ORDER BY` does by default. A single file is read on one thread, so nothing waits.

Each local regular file can be decompressed and parsed on a thread of its own, ahead of the rows being produced. This is decided once for the whole scan: it happens only if every file in the scan, each taking one of DuckDB's threads plus one more thread (two for a pair), fits in `threads` at once. With `SET threads = 8` that is up to four single-end files or two pairs. The count covers every file even when DuckDB reads them one at a time, as it does for `COPY ... (FORMAT FASTQ)` with insertion order kept. Each read-ahead stops once it holds 4 MiB of records, or `max_batch_bytes` if that is smaller, passing it by about a record (by up to sixteen where record lengths jump). Remote files and pipes (stdin from a pipe included) are read without it, as is everything in DuckDB-Wasm.

#### Remote reads and retries

An `http://` or `https://` file whose connection drops while it opens or partway through is retried by DuckDB's HTTP settings, `http_retries`, `http_retry_wait_ms` and `http_retry_backoff`. A retried read asks for the same byte range again, not the whole file. An HTTP error status DuckDB does not treat as transient, such as a 404, or a file whose ETag changed since it opened, fails at once. DuckDB also limits retries per host. After many failures on one host, a later drop can fail immediately.

The default retry settings give up within a few seconds. Raise them for a long download if needed.

```sql
SET http_retries = 8;
SET http_retry_wait_ms = 500;
SET http_retry_backoff = 2;
```

#### More SQL examples

```sql
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
