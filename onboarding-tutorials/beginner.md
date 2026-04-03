# Beginner: Exploring Microbiome Sequences with Python

This tutorial introduces miint through Python and Jupyter notebooks using the
[DuckDB Python client](https://duckdb.org/docs/current/clients/python/overview).
No SQL or command-line experience is required &mdash; every example returns a
[pandas](https://pandas.pydata.org/) DataFrame you can inspect, plot, or save.

By the end you will be able to read a FASTQ file directly from the internet,
explore its contents, and save the results as a Parquet file for downstream
analysis.

> **Prerequisites:** Python 3.9+, a Jupyter environment (JupyterLab, VS Code,
> Google Colab, etc.), and the `duckdb` and `pandas` packages.
> Install them with `pip install duckdb pandas matplotlib`.

## Setup

Create a new notebook and run this cell once to install the extensions and
open a connection:

```python
import duckdb

con = duckdb.connect()

# httpfs lets DuckDB read files over HTTPS
# https://duckdb.org/docs/current/core_extensions/httpfs/overview
con.sql("INSTALL httpfs; LOAD httpfs;")

# miint adds bioinformatics functions
con.sql("INSTALL miint FROM community; LOAD miint;")
```

The extensions are cached locally after the first install, so subsequent
sessions only need the `LOAD` statements.

## A quick introduction to SQL

Every example in this tutorial uses a small SQL query inside `con.sql(...)`.
SQL (Structured Query Language) reads like an English sentence once you know a
few keywords:

| Keyword | What it does | Plain English |
|---|---|---|
| [`SELECT`](https://duckdb.org/docs/current/sql/statements/select) | Choose which columns to return | "Give me these columns..." |
| `FROM` | Specify the data source (a table, file, or function) | "...from this data..." |
| [`WHERE`](https://duckdb.org/docs/current/sql/query_syntax/where) | Filter rows by a condition | "...but only rows where..." |
| [`LIMIT`](https://duckdb.org/docs/current/sql/query_syntax/limit) | Return at most *n* rows | "...and stop after *n* rows" |
| [`ORDER BY`](https://duckdb.org/docs/current/sql/query_syntax/orderby) | Sort the results | "...sorted by this column" |
| [`GROUP BY`](https://duckdb.org/docs/current/sql/query_syntax/groupby) | Group rows to compute summaries (count, average, ...) | "...one row per group" |
| `AS` | Rename a column in the output | "...and call it this" |

A minimal query looks like:

```text
SELECT column1, column2
FROM some_data
WHERE column1 > 10
LIMIT 5
```

You don't need to memorize anything &mdash; each example below explains the
SQL it uses. For a broader introduction, see DuckDB's
[SQL Introduction](https://duckdb.org/docs/current/sql/introduction).

## Reading a FASTQ file

The
[`read_fastx`](../docs/table-functions.md#read_fastxfilename-sequence2filename-include_filepathfalse-qual_offset33)
table function reads FASTA and FASTQ files. It works with local paths, glob
patterns, and HTTPS URLs. Let's read an amplicon sequencing run from the
[American Gut Project](https://americangut.org/) (McDonald et al., 2018):

```python
FASTQ_URL = (
    "https://ftp.sra.ebi.ac.uk/vol1/run/ERR107/"
    "ERR1074767/10317.000004216.fastq.gz"
)

df = con.sql(f"""
    SELECT read_id,
           comment,
           sequence1
    FROM read_fastx('{FASTQ_URL}')
    LIMIT 5
""").df()

print(df)
```

Let's break this query down:

- **`SELECT read_id, comment, sequence1`** &mdash; we want three columns from
  each sequencing read.
- **`FROM read_fastx(...)`** &mdash; the data comes from a FASTQ file. `read_fastx`
  is a miint function that opens the file and presents each read as a row.
- **`LIMIT 5`** &mdash; only return the first five reads (useful for a quick peek).

The columns are:

| Column | Description |
|---|---|
| `read_id` | Unique identifier for the read |
| `comment` | Header metadata (here, barcode information) |
| `sequence1` | The nucleotide sequence (A, C, G, T) |

There are additional columns such as `qual1` (per-base quality scores) and
`sequence_index` (0-based read number) &mdash; see the
[`read_fastx` reference](../docs/table-functions.md#read_fastxfilename-sequence2filename-include_filepathfalse-qual_offset33)
for the full schema.

## How many reads are there?

SQL has built-in
[aggregate functions](https://duckdb.org/docs/current/sql/functions/aggregates)
for summarizing data. `count(*)` counts rows, and `min(...)` / `max(...)` find
the smallest and largest values:

```python
counts = con.sql(f"""
    SELECT count(*) AS n_reads,
           min(len(sequence1)) AS min_length,
           max(len(sequence1)) AS max_length
    FROM read_fastx('{FASTQ_URL}')
""").df()

print(counts)
```

Here `len(sequence1)` computes the length of each sequence, and the `AS`
keyword gives each result column a readable name. There is no `LIMIT` because
the summary is already a single row.

All 20,939 reads are exactly 151 bp &mdash; typical for Illumina amplicon
sequencing of the 16S rRNA V4 region.

## Quality scores

Every base call in a FASTQ file has an associated
[Phred quality score](https://en.wikipedia.org/wiki/Phred_quality_score). A
score of 30 means a 1-in-1,000 chance the base call is wrong; 40 means
1-in-10,000.

`read_fastx` returns quality scores as a list of integers in the `qual1`
column. We can compute the mean quality per read and bring the result into
pandas:

```python
quality = con.sql(f"""
    SELECT read_id,
           round(list_aggregate(qual1, 'avg'), 1)
               AS mean_quality
    FROM read_fastx('{FASTQ_URL}')
""").df()

print(quality.describe())
```

[`list_aggregate(qual1, 'avg')`](https://duckdb.org/docs/current/sql/functions/list)
takes the list of per-base scores for each read and computes their average.
[`round(..., 1)`](https://duckdb.org/docs/current/sql/functions/numeric)
rounds to one decimal place.

### Plotting the distribution

```python
import os, tempfile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

outdir = tempfile.gettempdir()

fig, ax = plt.subplots(figsize=(8, 4))
quality["mean_quality"].hist(bins=30, edgecolor="black", ax=ax)
ax.set_xlabel("Mean Phred quality score")
ax.set_ylabel("Number of reads")
ax.set_title("Per-read quality distribution (ERR1074767)")
quality_hist = os.path.join(outdir, "miint_quality_hist.png")
fig.savefig(quality_hist, dpi=150, bbox_inches="tight")
print(f"Saved to {quality_hist}")
```

Most reads cluster around Phred 37&ndash;38, indicating high-quality data.

## GC content

GC content (the fraction of bases that are G or C) varies across microbial
taxa and can reveal composition at a glance. Let's compute it for every read:

```python
gc = con.sql(f"""
    SELECT read_id,
           round(
               len(replace(replace(upper(sequence1), 'A', ''), 'T', ''))
               * 100.0 / len(sequence1),
           1) AS gc_pct
    FROM read_fastx('{FASTQ_URL}')
""").df()

fig, ax = plt.subplots(figsize=(8, 4))
gc["gc_pct"].hist(bins=40, edgecolor="black", ax=ax)
ax.set_xlabel("GC content (%)")
ax.set_ylabel("Number of reads")
ax.set_title("GC content distribution (ERR1074767)")
gc_hist = os.path.join(outdir, "miint_gc_hist.png")
fig.savefig(gc_hist, dpi=150, bbox_inches="tight")
print(f"Saved to {gc_hist}")
```

The trick here is `replace(replace(upper(sequence1), 'A', ''), 'T', '')` which
removes all A and T characters, leaving only G and C. The length of what
remains, divided by the original length, gives the GC fraction.

A broad GC distribution is expected in an amplicon dataset &mdash; different
microbial taxa have different GC content, and the sample contains a mixture of
organisms.

## Saving results to Parquet

[Parquet](https://parquet.apache.org/) is a columnar file format widely used in
data science. It is fast to read, compact, and supported by pandas, Polars,
Spark, and many other tools. Saving your sequences as Parquet makes them easy
to work with outside of miint:

```python
parquet_path = os.path.join(outdir, "miint_sequences.parquet")

con.sql(f"""
    COPY (
        SELECT sequence_index,
               read_id,
               sequence1,
               round(list_aggregate(qual1, 'avg'), 1)
                   AS mean_quality
        FROM read_fastx('{FASTQ_URL}')
    ) TO '{parquet_path}' (FORMAT parquet)
""")

print(f"Wrote {parquet_path}")
```

[`COPY (...) TO`](https://duckdb.org/docs/current/sql/statements/copy) takes
any query and writes its output to a file. The `FORMAT parquet` option selects
the Parquet format (you could also use `FORMAT csv`).

You can read the file back with pandas to confirm it worked:

```python
import pandas as pd

df = pd.read_parquet(parquet_path)
print(f"{len(df)} reads saved")
print(df.head())
```

## What's next?

You've learned how to:
- Read a remote FASTQ file into a pandas DataFrame
- Compute per-read quality and GC content
- Save processed data to Parquet

In the [intermediate tutorial](intermediate.md), we'll align these reads
against the *E. coli* K-12 genome using
[minimap2](https://github.com/lh3/minimap2) (Li, 2018) and learn how to tell
real alignments from noise.
