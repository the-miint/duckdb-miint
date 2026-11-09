# GTDB fixtures

Real GTDB release files cut to a few rows, for `read_gtdb_metadata` and `read_gtdb_qc_failed`. They sit at the
paths GTDB serves them from, under `releases/`, so a test can read a release's directory or point the readers at
this tree in place of `https://data.gtdb.ecogenomic.org/releases`.

GTDB's data are released under the [Creative Commons Attribution-ShareAlike 4.0 International
licence](https://creativecommons.org/licenses/by-sa/4.0/) (CC BY-SA 4.0). These files are adapted from them, by
removing rows, and are shared under the same licence. Source: the Genome Taxonomy Database,
<https://gtdb.ecogenomic.org/>.

## How they were cut

On 2026-11-09 the first bytes of each file were fetched with an HTTP range request, which is enough for GTDB's
header line and a few hundred rows. Each fixture keeps the header line and the rows listed below, byte for byte
and in the file's order. The metadata files are gzipped again with no name or timestamp in the gzip header
(Python's `gzip.compress(data, mtime=0)`); `qc_failed` is plain text, as GTDB serves it.

## Files

Under `releases/`:

| File | Columns | Rows |
|---|---|---|
| `release232/232.0/bac120_metadata_r232.tsv.gz` | 113 | `RS_GCF_035810195.1`, `RS_GCF_001697265.1`, `RS_GCF_030769525.1`, `RS_GCF_026891715.1` |
| `release232/232.0/ar53_metadata_r232.tsv.gz` | 113 | `RS_GCF_000485535.1`, `GB_GCA_041569275.1` |
| `release232/232.0/auxillary_files/qc_failed_r232.tsv` | 26 | the first 3 |
| `release214/214.1/bac120_metadata_r214.tsv.gz` | 110 | `RS_GCF_000381545.1`, `RS_GCF_001206855.1` |
| `release214/214.1/ar53_metadata_r214.tsv.gz` | 110 | `RS_GCF_000246985.2` |
| `release214/214.1/auxillary_files/qc_failed_r214.tsv` | 20 | the first 2 |

Why these rows:
- **Release 232's bacteria** are the three genomes of *Bacillus subtilis* (NCBI taxid 1423) in the fetched rows,
  of which `RS_GCF_001697265.1` is the type strain (`gtdb_type_designation_ncbi_taxa` is
  `type strain of species`). They show the link from an NCBI species to its type strain's genome. The fourth,
  `RS_GCF_026891715.1`, has an `ncbi_submitter` GTDB writes CSV-quoted:
  `"National Research Center ""Kurchatov Institute"""`.
- **Release 232's archaea** are a type strain and a genome that is not type material.
- **Release 214.1** is the last release before CheckM2's columns, so its files have 110 columns to 232's 113.
- **Every metadata row** holds GTDB's missing-value marker `none`: 54 times in all in release 232's bacterial file.
