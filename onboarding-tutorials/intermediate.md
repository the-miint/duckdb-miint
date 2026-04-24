# Intermediate: Alignment and Quality Filtering with the DuckDB CLI

This tutorial uses the
[DuckDB command-line interface](https://duckdb.org/docs/current/clients/cli/overview)
to align reads against a reference genome, assess alignment quality, and count
taxa. You will learn
[common table expressions](https://duckdb.org/docs/current/sql/query_syntax/with),
[views](https://duckdb.org/docs/current/sql/statements/create_view),
[macros](https://duckdb.org/docs/current/sql/statements/create_macro), and
&mdash; most importantly &mdash; why query coverage filtering matters.

> **Prerequisites:** The [beginner tutorial](beginner.md). Install the DuckDB
> command-line interface:
>
> ```bash
> pip install duckdb-cli               # pip
> uv pip install duckdb-cli            # or uv
> conda install -c conda-forge duckdb-cli  # or conda
> ```
>
> See the
> [DuckDB installation page](https://duckdb.org/docs/current/installation/)
> for other options (Homebrew, standalone binaries, etc.).

## Setup

Open a terminal and start the DuckDB shell by running `duckdb`. You'll see an
interactive prompt where you can type SQL queries. Load the extensions we need:

```sql
INSTALL httpfs;
LOAD httpfs;
INSTALL miint FROM community;
LOAD miint;
```

See the DuckDB documentation on
[installing and loading extensions](https://duckdb.org/docs/current/sql/statements/load_and_install)
for details.

We'll use the same American Gut Project FASTQ from the beginner tutorial, plus
the *Escherichia coli* K-12 MG1655 reference genome (Blattner et al., 1997)
fetched from NCBI.

## Step 1: Load the data

### Fetch the reference genome

[`read_ncbi_fasta`](../docs/table-functions.md#read_ncbi_fastaaccession-api_key-include_filepathfalse)
downloads a FASTA record from NCBI by accession and returns it as a table.
We'll save it into a local table with
[`CREATE TABLE ... AS`](https://duckdb.org/docs/current/sql/statements/create_table)
so we only download it once:

```sql
CREATE TABLE ref_genome AS
    SELECT read_id, sequence1
    FROM read_ncbi_fasta('U00096.3');
```

[`U00096.3`](https://www.ncbi.nlm.nih.gov/nuccore/U00096.3) is the *E. coli*
K-12 MG1655 complete genome (4,641,652 bp).

### Load the reads

```sql
CREATE TABLE reads AS
    SELECT *
    FROM read_fastx(
        'https://ftp.sra.ebi.ac.uk/vol1/run/ERR107/ERR1074767/10317.000004216.fastq.gz'
    );
```

This loads all 20,939 reads into memory. For larger files you could query the
file directly, but loading into a table avoids re-downloading on every query.

## Step 2: Align reads to the reference

[`align_minimap2`](../docs/table-functions.md#align_minimap2query_table-subject_tablenull-index_pathnull-options)
aligns reads against a reference using
[minimap2](https://github.com/lh3/minimap2) (Li, 2018). It reads sequences
from DuckDB tables and returns standard SAM alignment fields:

```sql
CREATE TABLE alignments AS
    SELECT *
    FROM align_minimap2(
        'reads',
        subject_table='ref_genome',
        preset='sr',
        max_secondary=5
    );
```

- `preset='sr'` &mdash; short-read mode, appropriate for 151 bp Illumina data.
- `max_secondary=5` &mdash; report up to 5 secondary alignments per read.
  *E. coli* has seven copies of the 16S rRNA operon, so a single read can
  align to multiple locations.

How many alignments did we get?

```sql
SELECT count(*) AS n_alignments FROM alignments;
```

About 199 alignments from 20,939 reads. Most reads in this diverse microbiome
sample come from organisms other than *E. coli*, so a low mapping rate is
expected.

## Step 3: Primary vs. secondary alignments

Because we allowed secondary alignments (`max_secondary=5`), a single read can
appear in the results multiple times. Let's see the breakdown using miint's
[SAM flag functions](../docs/scalar-functions.md#sam-flag-functions):

```sql
SELECT alignment_is_primary(flags) AS is_primary,
       count(*) AS n
FROM alignments
GROUP BY is_primary;
```

| is_primary | n |
|---|---|
| true | 34 |
| false | 165 |

Only 34 of the 199 rows are primary alignments &mdash; the rest are secondary
hits to other copies of the rRNA operon.

For single-end data like ours, filtering to **primary alignments** removes
duplicate counting. For paired-end data, you would also want to filter to
**proper pairs** using
[`alignment_is_proper_pair(flags)`](../docs/scalar-functions.md#sam-flag-functions),
which ensures both mates mapped in the expected orientation and distance. miint
provides flag-checking functions for all standard
[SAM flags](../docs/scalar-functions.md#sam-flag-functions).

## Step 4: Examine alignment quality

Even among primary alignments, not all are trustworthy. Two key metrics:

- [`alignment_query_coverage(cigar)`](../docs/scalar-functions.md#alignment_query_coveragecigar-typealigned)
  &mdash; what fraction of the read actually aligned? A value of 1.0 means the
  entire read matched; 0.6 means 40% of the read was soft-clipped (ignored by
  the aligner).

- [`alignment_seq_identity(cigar, nm, md, type)`](../docs/scalar-functions.md#alignment_seq_identitycigar-nm-md-type)
  &mdash; of the bases that *did* align, what fraction are identical to the
  reference?

```sql
SELECT read_id,
       round(alignment_query_coverage(cigar), 2) AS query_cov,
       round(alignment_seq_identity(cigar, tag_nm, tag_md, 'blast'), 3)
           AS seq_identity,
       cigar
FROM alignments
WHERE alignment_is_primary(flags)
ORDER BY alignment_query_coverage(cigar) ASC
LIMIT 10;
```

[`WHERE`](https://duckdb.org/docs/current/sql/query_syntax/where) filters to
primary alignments,
[`ORDER BY`](https://duckdb.org/docs/current/sql/query_syntax/orderby) sorts
by query coverage (worst first), and
[`LIMIT`](https://duckdb.org/docs/current/sql/query_syntax/limit) caps the
output at 10 rows.

Notice the pattern: many alignments have ~90% sequence identity but only
60&ndash;70% query coverage. These reads have large soft clips (`S` operations
in the CIGAR string), meaning the aligner could only match a portion of the
read. The matched portion looks decent, but the overall alignment is poor.

## Step 5: The query coverage lesson

This is the single most important quality-control concept when working with
alignments: **sequence identity alone can be misleading**.

A read with 90% identity and 60% query coverage means only 60% of the read
aligned, and of that portion, 90% of bases matched. The other 40% of the read
didn't match the reference at all.

Let's use a
[Common Table Expression](https://duckdb.org/docs/current/sql/query_syntax/with)
(CTE) to see the full picture. A CTE is introduced with `WITH` and creates a
temporary named result you can query like a table:

```sql
WITH quality AS (
    SELECT
        read_id,
        alignment_query_coverage(cigar) AS qcov,
        alignment_seq_identity(cigar, tag_nm, tag_md, 'blast') AS identity
    FROM alignments
    WHERE alignment_is_primary(flags)
)
SELECT
    CASE
        WHEN qcov >= 0.9 THEN 'high (>= 90%)'
        WHEN qcov >= 0.8 THEN 'medium (80-90%)'
        ELSE 'low (< 80%)'
    END AS coverage_bin,
    count(*) AS n_alignments,
    round(min(identity), 3) AS min_identity,
    round(max(identity), 3) AS max_identity,
    round(avg(identity), 3) AS mean_identity
FROM quality
GROUP BY coverage_bin
ORDER BY coverage_bin;
```

The
[`CASE`](https://duckdb.org/docs/current/sql/expressions/case) expression
assigns each alignment to a bin based on its query coverage, and
[`GROUP BY`](https://duckdb.org/docs/current/sql/query_syntax/groupby) computes
summary statistics per bin. You should see something like:

| coverage_bin | n_alignments | min_identity | max_identity | mean_identity |
|---|---|---|---|---|
| high (>= 90%) | 11 | 0.914 | 1.000 | 0.963 |
| low (< 80%) | 23 | 0.890 | 0.932 | 0.907 |

**23 out of 34 primary alignments** (68%) have low query coverage. If you
filtered only on sequence identity (say, >= 85%), you would keep all of them.
Those are reads from other organisms that happen to share partial 16S rRNA
homology with *E. coli* &mdash; not genuine *E. coli* reads.

## Step 6: Filter and count with `woltka_ogu`

Now let's apply a proper filter and count with
[`woltka_ogu`](../docs/analysis-functions.md#woltka_ogurelation-sequence_id_field-sample_id),
a built-in table macro that implements the Woltka OGU (Operational Genomic
Unit) counting method (Zhu et al., 2022). `woltka_ogu` is primarily designed
for shotgun metagenomic data, where reads span the entire genome; however, its
handling of multi-mapped reads (splitting counts proportionally across
references) is equally useful here, and the workflow translates directly to
metagenomic analyses.

First, create a
[view](https://duckdb.org/docs/current/sql/statements/create_view) with our
quality filter. A view is a saved query &mdash; it doesn't copy data, it
re-runs the filter each time you reference it:

```sql
CREATE VIEW filtered_alignments AS
    SELECT *
    FROM alignments
    WHERE alignment_is_primary(flags)
      AND alignment_query_coverage(cigar) >= 0.8;
```

Now compare unfiltered vs. filtered counts using
[`UNION ALL`](https://duckdb.org/docs/current/sql/query_syntax/setops):

```sql
SELECT 'unfiltered' AS source, *
    FROM woltka_ogu('alignments', 'read_id')
UNION ALL
SELECT 'filtered' AS source, *
    FROM woltka_ogu('filtered_alignments', 'read_id');
```

| source | feature_id | value |
|---|---|---|
| unfiltered | U00096.3 | 34.0 |
| filtered | U00096.3 | 11.0 |

Without the query coverage filter, you would report 3x more *E. coli* than is
actually present. Those extra counts come from reads that share partial 16S
sequence with *E. coli* but belong to other organisms.

## Step 7: Tightening the filter

All 11 reads that pass the query coverage filter have 100% query coverage
&mdash; the entire read aligned. But their sequence identity ranges from 91%
to 100%. How many survive at stricter identity thresholds?

```sql
CREATE VIEW strict_97 AS
    SELECT *
    FROM alignments
    WHERE alignment_is_primary(flags)
      AND alignment_query_coverage(cigar) >= 0.99
      AND alignment_seq_identity(cigar, tag_nm, tag_md, 'blast') >= 0.97;

CREATE VIEW strict_99 AS
    SELECT *
    FROM alignments
    WHERE alignment_is_primary(flags)
      AND alignment_query_coverage(cigar) >= 0.99
      AND alignment_seq_identity(cigar, tag_nm, tag_md, 'blast') >= 0.99;

SELECT 'qcov >= 99%, id >= 97%' AS threshold, *
    FROM woltka_ogu('strict_97', 'read_id')
UNION ALL
SELECT 'qcov >= 99%, id >= 99%' AS threshold, *
    FROM woltka_ogu('strict_99', 'read_id');
```

| threshold | feature_id | value |
|---|---|---|
| qcov >= 99%, id >= 97% | U00096.3 | 5.0 |
| qcov >= 99%, id >= 99% | U00096.3 | 4.0 |

At 97% identity, five reads confidently match *E. coli*. At 99%, four remain
&mdash; including three that are perfect matches (`151=` in the CIGAR). The
reads between 91&ndash;97% identity likely come from closely related
Enterobacteriaceae that share near-identical 16S rRNA with *E. coli*.

The right threshold depends on your question. For a conservative "is *E. coli*
present?" answer, 99% identity with full query coverage is hard to argue with.

## Step 8: Create a reusable macro

If you find yourself applying the same filter repeatedly, you can wrap it in a
DuckDB [macro](https://duckdb.org/docs/current/sql/statements/create_macro). A
table macro is a parameterized query you can call like a table:

```sql
CREATE OR REPLACE MACRO quality_filtered(
    min_query_coverage := 0.99,
    min_seq_identity := 0.97
) AS TABLE
    SELECT *
    FROM alignments
    WHERE alignment_is_primary(flags)
      AND alignment_query_coverage(cigar) >= min_query_coverage
      AND alignment_seq_identity(cigar, tag_nm, tag_md, 'blast')
          >= min_seq_identity;
```

`woltka_ogu` resolves its `relation` argument through the catalog, so wrap the
macro call in a view before counting:

```sql
CREATE OR REPLACE VIEW strict_filtered AS SELECT * FROM quality_filtered();
SELECT * FROM woltka_ogu('strict_filtered', 'read_id');
```

Or adjust the thresholds:

```sql
CREATE OR REPLACE VIEW strictest_filtered AS
    SELECT * FROM quality_filtered(min_seq_identity := 0.99);
SELECT * FROM woltka_ogu('strictest_filtered', 'read_id');
```

## Step 9: Genome coverage

How much of the *E. coli* genome did our reads actually cover? The
[`genome_coverage`](../docs/analysis-functions.md#genome_coveragealignments-subject_total_length-subject_genome_id)
macro computes this by merging overlapping alignment intervals using
[`compress_intervals`](../docs/analysis-functions.md#compress_intervalsstart-stop):

```sql
CREATE TABLE ecoli_genome_id AS
    SELECT 'U00096.3' AS contig_id,
           'ecoli_k12' AS genome_id;

CREATE TABLE ecoli_total_length AS
    SELECT 'ecoli_k12' AS genome_id,
           len(sequence1) AS total_length
    FROM ref_genome;

SELECT genome_id,
       covered,
       round(proportion_covered * 100, 4) AS pct_covered
FROM genome_coverage(
    filtered_alignments,
    ecoli_total_length,
    ecoli_genome_id
);
```

| genome_id | covered | pct_covered |
|---|---|---|
| ecoli_k12 | 906 | 0.0195 |

About 906 bases covered (0.02% of the genome). This is expected: 16S V4
amplicons target a ~250 bp region, and *E. coli* has seven rRNA operons.
Coverage is concentrated in those operons, not spread across the genome.

## What's next?

You've learned how to:
- Align reads against a reference genome with minimap2
- Filter to primary alignments (and, for paired-end data, proper pairs)
- Evaluate alignments using query coverage and sequence identity
- Separate genuine alignments from partial-homology noise
- Choose identity thresholds appropriate to your biological question
- Count taxa with `woltka_ogu` and wrap filtering logic in a macro
- Measure genome coverage

In the [advanced tutorial](advanced.md), we'll run these analyses from the
command line using Bash, batch-process multiple samples, and write complex
analytical SQL.
