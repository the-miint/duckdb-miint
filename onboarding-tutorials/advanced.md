# Advanced: Bash Pipelines and Complex SQL

This tutorial shows how to use DuckDB from the command line for batch
processing, complex analytical queries, and integration with shell workflows.

> **Prerequisites:** The [intermediate tutorial](intermediate.md). Comfort with
> Bash scripting and SQL.

## Using DuckDB from Bash

The
[DuckDB CLI](https://duckdb.org/docs/current/clients/cli/overview) accepts SQL
directly via the `-c` flag, making it easy to embed queries in shell scripts
and one-liners:

```bash
duckdb -c "
INSTALL httpfs; LOAD httpfs;
INSTALL miint FROM community; LOAD miint;
SELECT count(*) AS n_reads
FROM read_fastx(
    'https://ftp.sra.ebi.ac.uk/vol1/run/ERR107/ERR1074767/10317.000004216.fastq.gz'
);
"
```

This is DuckDB's equivalent of `awk` or `csvkit` &mdash; a full analytical
engine you can call from any script.

## Batch processing multiple samples

A common workflow is to process many FASTQ files identically. Let's simulate
this by splitting our dataset into three synthetic "samples" and then
processing each one in a loop.

### Create the samples

```bash
OUTDIR="${TMPDIR:-/tmp}"

duckdb -c "
INSTALL httpfs; LOAD httpfs;
INSTALL miint FROM community; LOAD miint;

COPY (
    SELECT * FROM read_fastx(
        'https://ftp.sra.ebi.ac.uk/vol1/run/ERR107/ERR1074767/10317.000004216.fastq.gz'
    ) WHERE sequence_index % 3 = 0
) TO '${OUTDIR}/sample_a.parquet' (FORMAT parquet);

COPY (
    SELECT * FROM read_fastx(
        'https://ftp.sra.ebi.ac.uk/vol1/run/ERR107/ERR1074767/10317.000004216.fastq.gz'
    ) WHERE sequence_index % 3 = 1
) TO '${OUTDIR}/sample_b.parquet' (FORMAT parquet);

COPY (
    SELECT * FROM read_fastx(
        'https://ftp.sra.ebi.ac.uk/vol1/run/ERR107/ERR1074767/10317.000004216.fastq.gz'
    ) WHERE sequence_index % 3 = 2
) TO '${OUTDIR}/sample_c.parquet' (FORMAT parquet);
"
```

We use
[`COPY ... TO`](https://duckdb.org/docs/current/sql/statements/copy) with
[Parquet](https://parquet.apache.org/) format so each sample is a self-contained
file that DuckDB can read back efficiently.

### Process each sample in a loop

```bash
OUTDIR="${TMPDIR:-/tmp}"

for sample in "${OUTDIR}"/sample_a.parquet "${OUTDIR}"/sample_b.parquet "${OUTDIR}"/sample_c.parquet; do
    name=$(basename "$sample" .parquet)
    duckdb -c "
    LOAD miint;
    SELECT '${name}' AS sample,
           count(*) AS n_reads,
           round(avg(len(sequence1)), 0) AS avg_len,
           round(avg(list_aggregate(qual1, 'avg')), 1) AS avg_quality
    FROM read_parquet('${sample}');
    "
done
```

Each iteration spawns a fresh DuckDB process, reads one Parquet file, and
prints summary statistics. In a real pipeline you might align each sample to a
reference and write per-sample OGU tables.

### Align and count per sample

Building on the [intermediate tutorial](intermediate.md), here is a complete
per-sample pipeline:

```bash
OUTDIR="${TMPDIR:-/tmp}"

# Fetch the reference once and save as Parquet
duckdb -c "
INSTALL httpfs; LOAD httpfs;
INSTALL miint FROM community; LOAD miint;
COPY (
    SELECT read_id, sequence1
    FROM read_ncbi_fasta('U00096.3')
) TO '${OUTDIR}/ecoli_ref.parquet' (FORMAT parquet);
"

for sample in "${OUTDIR}"/sample_a.parquet "${OUTDIR}"/sample_b.parquet "${OUTDIR}"/sample_c.parquet; do
    name=$(basename "$sample" .parquet)
    duckdb -c "
    INSTALL miint FROM community; LOAD miint;
    CREATE TABLE ref_genome AS SELECT * FROM read_parquet('${OUTDIR}/ecoli_ref.parquet');
    CREATE TABLE reads AS SELECT * FROM read_parquet('${sample}');
    CREATE TABLE alns AS
        SELECT * FROM align_minimap2('reads', subject_table='ref_genome', preset='sr');
    CREATE VIEW filtered AS
        SELECT * FROM alns
        WHERE alignment_is_primary(flags)
          AND cigar_query_coverage(cigar) >= 0.99
          AND alignment_seq_identity(cigar, tag_nm, tag_md, 'blast') >= 0.97;
    SELECT '${name}' AS sample, * FROM woltka_ogu('filtered', 'read_id');
    "
done
```

Notice two patterns here:
- The reference genome is saved as Parquet once and re-read from disk in each
  iteration, avoiding repeated NCBI downloads.
- Each iteration runs in a fresh in-memory DuckDB process &mdash; no persistent
  state to manage.

## Complex analytical SQL

The remaining examples in this tutorial are meant to be run interactively in
the DuckDB shell (`duckdb`), rather than inlined in Bash scripts. Start a new
shell session and paste the SQL blocks below.

### Verifying alignments against gene annotations

In the intermediate tutorial we saw that 16S V4 reads cover only a tiny
fraction of the *E. coli* genome. But *where* exactly do they land, and do
they match the known rRNA operons?

[`read_ncbi_annotation`](../docs/table-functions.md#read_ncbi_annotationaccession-api_key-include_filepathfalse)
fetches gene annotations from NCBI. Let's set up the full pipeline and
retrieve the 16S rRNA gene coordinates:

```sql
INSTALL httpfs;
LOAD httpfs;
INSTALL miint FROM community;
LOAD miint;

CREATE TABLE ref_genome AS
    SELECT read_id, sequence1
    FROM read_ncbi_fasta('U00096.3');

CREATE TABLE reads AS
    SELECT *
    FROM read_fastx(
        'https://ftp.sra.ebi.ac.uk/vol1/run/ERR107/ERR1074767/10317.000004216.fastq.gz'
    );

CREATE TABLE alignments AS
    SELECT *
    FROM align_minimap2(
        'reads',
        subject_table='ref_genome',
        preset='sr',
        max_secondary=5
    );

CREATE TABLE ncbi_annotations AS
    SELECT * FROM read_ncbi_annotation('U00096.3');

CREATE TABLE rrna_16s AS
    SELECT position AS gene_start,
           stop_position AS gene_end,
           strand
    FROM ncbi_annotations
    WHERE type = 'rRNA'
      AND attributes['product'] LIKE '%16S%';

SELECT * FROM rrna_16s
ORDER BY gene_start;
```

| gene_start | gene_end | strand |
|---|---|---|
| 223771 | 225312 | + |
| 2729616 | 2731157 | - |
| 3427221 | 3428762 | - |
| 3941808 | 3943349 | + |
| 4035531 | 4037072 | + |
| 4166659 | 4168200 | + |
| 4208147 | 4209688 | + |

Seven 16S rRNA operons, as expected.

### Crossing alignments with gene coordinates

Now, merge the alignment positions using
[`compress_intervals`](../docs/analysis-functions.md#compress_intervalsstart-stop)
and
[`JOIN`](https://duckdb.org/docs/current/sql/query_syntax/from) them against
the gene annotations to see which operons our reads hit:

We filter to primary alignments with high
[query coverage](../docs/scalar-functions.md#cigar_query_coveragecigar-typealigned)
(&ge; 99% of the read aligned) to exclude partial soft-clipped matches that
inflate counts &mdash; see the
[intermediate tutorial](intermediate.md#step-5-the-query-coverage-lesson) for
a detailed explanation of why this matters.

```sql
WITH good_alns AS (
    SELECT reference, position, stop_position
    FROM alignments
    WHERE alignment_is_primary(flags)
      AND cigar_query_coverage(cigar) >= 0.99
),
compressed AS (
    SELECT reference,
           unnest(compress_intervals(position, stop_position)) AS ci
    FROM good_alns
    GROUP BY reference
)
SELECT r.gene_start,
       r.gene_end,
       r.strand AS gene_strand,
       ci.start AS read_start,
       ci.stop AS read_end,
       ci.stop - ci.start AS read_length
FROM compressed c
JOIN rrna_16s r
  ON ci.start >= r.gene_start
 AND ci.stop <= r.gene_end
ORDER BY r.gene_start;
```

| gene_start | gene_end | gene_strand | read_start | read_end | read_length |
|---|---|---|---|---|---|
| 2729616 | 2731157 | - | 2730474 | 2730625 | 151 |
| 3427221 | 3428762 | - | 3428079 | 3428230 | 151 |
| 3941808 | 3943349 | + | 3942341 | 3942492 | 151 |
| 4166659 | 4168200 | + | 4167192 | 4167343 | 151 |
| 4208147 | 4209688 | + | 4208680 | 4208831 | 151 |

Five out of seven operons have primary alignments, each with exactly one
151 bp region. The V4 amplicon reads land at a consistent offset within each
operon &mdash; precisely what you would expect from a targeted amplicon
protocol.

### Primary vs. secondary: seeing the full picture

The two operons without primary hits (*rrsA* at 223,771 and *rrsE* at
4,035,531) aren't missing from the data &mdash; they show up as secondary
alignments. We can see this by looking at all alignments for a single read:

```sql
SELECT alignment_is_primary(a.flags) AS is_primary,
       a.position,
       a.stop_position,
       a.cigar,
       r.gene_start AS in_16s_operon
FROM alignments a
LEFT JOIN rrna_16s r
  ON a.position >= r.gene_start
 AND a.stop_position <= r.gene_end
WHERE a.read_id = '10317.000004216_6975'
ORDER BY alignment_is_primary(a.flags) DESC, a.position;
```

| is_primary | position | stop_position | cigar | in_16s_operon |
|---|---|---|---|---|
| true | 4167192 | 4167343 | 151= | 4166659 |
| false | 2730474 | 2730625 | 151= | 2729616 |
| false | 3428079 | 3428230 | 151= | 3427221 |
| false | 3942341 | 3942492 | 151= | 3941808 |
| false | 4036064 | 4036215 | 151= | 4035531 |
| false | 4208680 | 4208831 | 151= | 4208147 |

One primary alignment plus five secondaries &mdash; all perfect matches
(`151=`), all falling within annotated 16S operons. That accounts for six of
the seven operons. The seventh (*rrsA*) was not reported because we set
`max_secondary=5`.

minimap2 must choose one location for the primary alignment, but the choice
among identical targets is arbitrary. We can confirm the V4 region is identical
across all seven operons:

```sql
CREATE TABLE v4_regions AS
    SELECT r.gene_start,
           r.strand,
           CASE
               WHEN r.strand = '+' THEN
                   substring(g.sequence1, r.gene_start + 533, 151)
               ELSE
                   sequence_dna_reverse_complement(
                       substring(g.sequence1, r.gene_end - 533 - 151 + 1, 151)
                   )
           END AS v4_sequence
    FROM rrna_16s r
    CROSS JOIN ref_genome g;

SELECT a.gene_start,
       a.strand,
       a.v4_sequence = b.v4_sequence AS identical_to_ref
FROM v4_regions a
CROSS JOIN v4_regions b
WHERE b.gene_start = 3941808
ORDER BY a.gene_start;
```

| gene_start | strand | identical_to_ref |
|---|---|---|
| 223771 | + | true |
| 2729616 | - | true |
| 3427221 | - | true |
| 3941808 | + | true |
| 4035531 | + | true |
| 4166659 | + | true |
| 4208147 | + | true |

All seven operons have **identical** V4 sequences. minimap2 cannot distinguish
which operon a read came from, so it picks one as primary and reports the rest
as secondary.

### An unexpected alignment outside 16S

Our `compress_intervals` query found six alignment regions, but the `JOIN`
matched only five to annotated 16S genes. Where is the sixth?

```sql
WITH good_alns AS (
    SELECT reference, position, stop_position
    FROM alignments
    WHERE alignment_is_primary(flags)
      AND cigar_query_coverage(cigar) >= 0.99
),
compressed AS (
    SELECT reference,
           unnest(compress_intervals(position, stop_position)) AS ci
    FROM good_alns
    GROUP BY reference
),
unmatched AS (
    SELECT c.ci
    FROM compressed c
    LEFT JOIN rrna_16s r
      ON ci.start >= r.gene_start
     AND ci.stop <= r.gene_end
    WHERE r.gene_start IS NULL
)
SELECT ci.start AS read_start,
       ci.stop AS read_end,
       a.type,
       a.position AS gene_start,
       a.stop_position AS gene_end,
       a.attributes['product'] AS product
FROM unmatched
JOIN ncbi_annotations a
  ON ci.start >= a.position
 AND ci.stop <= a.stop_position
WHERE a.type = 'CDS'
ORDER BY a.position;
```

| read_start | read_end | type | gene_start | gene_end | product |
|---|---|---|---|---|---|
| 1565200 | 1565351 | CDS | 1563334 | 1565733 | oxygen-sensing c-di-GMP phosphodiesterase DosP |

A read aligned with 100% identity and 100% query coverage to a region inside
the *dosP* gene, which encodes an oxygen-sensing phosphodiesterase &mdash; not
a 16S rRNA gene at all.

Which read produced this alignment? We can find it by joining back to the
alignments table on position:

```sql
SELECT read_id, position, cigar
FROM alignments
WHERE alignment_is_primary(flags)
  AND position = 1565200;
```

| read_id | position | cigar |
|---|---|---|
| 10317.000004216_6599 | 1565200 | 151= |

Unlike the 16S reads, which have multiple secondary alignments across operons,
this read has **no secondary alignments** &mdash; the DosP location is its
only hit in the entire genome:

```sql
SELECT alignment_is_primary(flags) AS is_primary,
       position,
       cigar,
       round(alignment_seq_identity(cigar, tag_nm, tag_md, 'blast'), 3)
           AS identity
FROM alignments
WHERE read_id = '10317.000004216_6599'
ORDER BY position;
```

| is_primary | position | cigar | identity |
|---|---|---|---|
| true | 1565200 | 151= | 1.0 |

And indeed, the read sequence does not appear in any of the seven 16S operons:

```sql
WITH dosp_seq AS (
    SELECT substring(sequence1, 1565200, 151) AS seq,
           sequence_dna_reverse_complement(
               substring(sequence1, 1565200, 151)
           ) AS revcomp
    FROM ref_genome
),
operons AS (
    SELECT r.gene_start,
           substring(g.sequence1, r.gene_start, r.gene_end - r.gene_start)
               AS gene_seq
    FROM rrna_16s r
    CROSS JOIN ref_genome g
)
SELECT o.gene_start,
       position(d.seq IN o.gene_seq) AS fwd_match,
       position(d.revcomp IN o.gene_seq) AS rev_match
FROM operons o, dosp_seq d;
```

| gene_start | fwd_match | rev_match |
|---|---|---|
| 223771 | 0 | 0 |
| 2729616 | 0 | 0 |
| 3427221 | 0 | 0 |
| 3941808 | 0 | 0 |
| 4035531 | 0 | 0 |
| 4166659 | 0 | 0 |
| 4208147 | 0 | 0 |

Zero matches in either orientation. To confirm there is no broader homology, we
can excise a 16S-length (1,541 bp) region around the DosP alignment &mdash;
using the same V4 offset of 533 bp observed in the real operons &mdash; and
align it against an annotated 16S gene using
[WFA2](../docs/analysis-functions.md#pairwise-alignment-functions) (Marco-Sola
et al., 2023):

```sql
SELECT
    align_pairwise_full(
        substring(sequence1, 1565200 - 533, 1541),
        substring(sequence1, 3941808, 1541)
    ).score AS wfa2_penalty
FROM ref_genome;
```

| wfa2_penalty |
|---|
| 3572 |

A WFA2 penalty of 3,572 (where 0 is a perfect match) confirms the *dosP*
region has no meaningful homology to the 16S rRNA gene. The 151 bp match is a
coincidental sequence identity &mdash; likely an off-target amplification
product that was picked up during PCR.

## What else can miint do?

These tutorials focused on amplicon sequencing against a single genome, but
miint supports a much broader range of bioinformatics workflows:

- **File formats** &mdash;
  [read](../docs/table-functions.md) FASTA, FASTQ, SAM, BAM, SFF, BIOM,
  Newick, GFF, mzML, mzXML, and jplace files;
  [write](../docs/copy-formats.md) FASTA, FASTQ, SAM, BAM, BIOM, and
  Newick.
- **Alignment** &mdash; align with
  [minimap2](../docs/table-functions.md#align_minimap2query_table-subject_tablenull-index_pathnull-options),
  [Bowtie2](../docs/table-functions.md#align_bowtie2query_table-subject_table-options),
  [MAFFT](../docs/table-functions.md#align_maffttable_name) (multiple sequence
  alignment), or compute
  [pairwise alignment scores](../docs/analysis-functions.md#pairwise-alignment-functions)
  with WFA2.
- **Quality control** &mdash;
  [detect chimeras](../docs/table-functions.md#detect_chimera_uchimequery_table-dbrefs_table-options)
  with UCHIME,
  [mask low-complexity regions](../docs/scalar-functions.md#mask_dustsequence-hardmaskfalse)
  with DUST,
  [merge overlapping paired-end reads](../docs/scalar-functions.md#merge_pairs_vsearchfwd_seq-fwd_qual-rev_seq-rev_qual-options).
- **Sequence analysis** &mdash;
  [reverse complement](../docs/analysis-functions.md#sequence_dna_reverse_complementsequence-and-sequence_rna_reverse_complementsequence),
  [degenerate base matching](../docs/analysis-functions.md#sequence_dna_as_regexpsequence-and-sequence_rna_as_regexpsequence),
  [sequence search](../docs/table-functions.md#search_sequences_vsearchquery_table-dbref_table-idthreshold-options)
  and
  [clustering](../docs/table-functions.md#cluster_sequences_vsearchinput_table-idthreshold-options).
- **Community ecology** &mdash;
  [OGU counting](../docs/analysis-functions.md#woltka_ogurelation-sequence_id_field-sample_id)
  per sample,
  [genome coverage](../docs/analysis-functions.md#genome_coveragealignments-subject_total_length-subject_genome_id),
  reading phylogenetic
  [placement data](../docs/table-functions.md#read_jplacepath) and
  [resolving placements onto trees](../docs/table-functions.md#tree_resolve_placementtree_table-placements_table).
- **Mass spectrometry** &mdash;
  [MassQL queries](../docs/massql.md),
  read [mzML](../docs/table-functions.md#read_mzmlfilename-include_filepathfalse)
  and
  [mzXML](../docs/table-functions.md#read_mzxmlfilename-include_filepathfalse)
  files, isotope pattern matching.
- **NCBI integration** &mdash;
  fetch
  [genomes](../docs/table-functions.md#read_ncbi_fastaaccession-api_key-include_filepathfalse),
  [metadata](../docs/table-functions.md#read_ncbiaccession-api_key), and
  [annotations](../docs/table-functions.md#read_ncbi_annotationaccession-api_key-include_filepathfalse)
  directly from NCBI by accession.

See the full [table functions](../docs/table-functions.md),
[scalar functions](../docs/scalar-functions.md),
[analysis functions](../docs/analysis-functions.md), and
[COPY formats](../docs/copy-formats.md) reference documentation.

## Extending miint

If the built-in functions don't cover your use case, you can extend miint
yourself. The extension is written in C++20 and follows DuckDB's extension
framework.

### Anatomy of a scalar function

Every scalar function in miint follows the same pattern:

1. **Define the function** with input types and return type.
2. **Implement the logic** in a C++ lambda or function.
3. **Register it** in `LoadInternal()` in
   `src/miint_extension.cpp`.

For example, `cigar_query_coverage` is defined in
`src/alignment_functions.cpp` and registered via a static `Register()` method.

### Getting started

1. Fork the [duckdb-miint](https://github.com/the-miint/duckdb-miint)
   repository.
2. Read the [DuckDB extension template](../docs/README.md) documentation.
3. Look at existing scalar functions in `src/` for patterns to follow.
4. Add tests in `test/sql/` and `test/cpp/` &mdash; see the
   [testing guide](../docs/testing.md).

DuckDB also provides an
[extension template](https://github.com/duckdb/extension-template) for
starting from scratch.

## Summary

Across these three tutorials, you've gone from reading a FASTQ file in a
Jupyter notebook to running batch alignment pipelines from the command line.
The core workflow &mdash; read sequences, align, filter by query coverage and
identity, count &mdash; applies whether you're checking for a single organism
or analyzing thousands of samples.
