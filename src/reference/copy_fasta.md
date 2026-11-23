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
