# Getting Started with miint

[miint](https://github.com/the-miint/duckdb-miint) is a
[DuckDB](https://duckdb.org/) extension for bioinformatics. It lets you read,
analyze, and write common bioinformatics file formats using SQL &mdash; or from
Python, R, and the command line. If you have never used SQL before, don't
worry: the [beginner tutorial](beginner.md) walks you through everything using
Python and Jupyter notebooks.

## What you'll learn

These tutorials follow a single question through three levels of complexity:

> **Does my sample contain *Escherichia coli*, and how confidently can I tell?**

Starting from a raw FASTQ file, you will read and explore the sequences, align
them against the *E. coli* K-12 reference genome, and learn how to separate
real signal from alignment noise. Each tutorial uses the same public dataset
&mdash; an amplicon sequencing run from the
[American Gut Project](https://americangut.org/) (McDonald et al., 2018), part
of the [Earth Microbiome Project](https://earthmicrobiome.org/) (Thompson
et al., 2017) &mdash; and builds on the ideas introduced in the previous one.

| Tutorial | Interface | You'll learn to... |
|---|---|---|
| [Beginner](beginner.md) | Python / Jupyter | Read sequences, explore quality, export to Parquet |
| [Intermediate](intermediate.md) | DuckDB CLI | Align reads, filter by quality metrics, count taxa |
| [Advanced](advanced.md) | Bash + SQL | Batch-process samples, compute coverage depth, extend miint |

## Installation

### Install DuckDB

DuckDB is available for many languages and platforms &mdash; see the
[DuckDB installation page](https://duckdb.org/docs/current/installation/) for
the full list. These tutorials use the
[Python client](https://duckdb.org/docs/current/clients/python/overview) and
the [command-line interface](https://duckdb.org/docs/current/clients/cli/overview);
installation instructions for each are below.

**Python library** (for the [beginner tutorial](beginner.md)):

```bash
pip install duckdb pyarrow pandas matplotlib           # pip
uv pip install duckdb pyarrow pandas matplotlib        # or uv
conda install -c conda-forge python-duckdb pyarrow pandas matplotlib  # or conda
```

**Command-line interface** (for the
[intermediate](intermediate.md) and [advanced](advanced.md) tutorials):

```bash
pip install duckdb-cli               # pip
uv pip install duckdb-cli            # or uv
conda install -c conda-forge duckdb-cli  # or conda
```

### Install miint

miint is published in the
[DuckDB Community Extensions](https://community-extensions.duckdb.org/) repository.
No compilation required &mdash; DuckDB will download and verify the extension
automatically:

```sql
INSTALL miint FROM community;
LOAD miint;
```

For building from source or other installation methods, see the
[installation guide](../docs/installation.md).

## The dataset

All tutorials use the same publicly available FASTQ file from the American Gut
Project:

| Field | Value |
|---|---|
| Run accession | [ERR1074767](https://www.ebi.ac.uk/ena/browser/view/ERR1074767) |
| Project | [American Gut Project](https://americangut.org/) / [Earth Microbiome Project](https://earthmicrobiome.org/) |
| Region | 16S rRNA gene, V4 hypervariable region |
| Read count | 20,939 |
| Read length | 151 bp |
| Platform | Illumina |

The file is hosted by the European Nucleotide Archive and can be read directly
over HTTPS &mdash; no download step needed.

## References and citations

These tutorials rely on several open-source tools and public datasets. Please
cite them if you use miint in published work.

- **DuckDB** &mdash; Raasveldt, M. & M&uuml;hleisen, H. (2019). DuckDB: an
  embeddable analytical database. *SIGMOD* 2019.
  [duckdb.org](https://duckdb.org/)

- **minimap2** &mdash; Li, H. (2018). Minimap2: pairwise alignment for
  nucleotide sequences. *Bioinformatics*, 34(18), 3094&ndash;3100.
  [doi:10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)

- **Woltka** &mdash; Zhu, Q. et al. (2022). Phylogeny-aware analysis of
  metagenome community ecology based on matched reference genomes while
  leveraging phylogeny. *mSystems*, 7(2).
  [doi:10.1128/msystems.00167-22](https://doi.org/10.1128/msystems.00167-22)

- **HTSlib** &mdash; Bonfield, J.K. et al. (2021). HTSlib: C library for
  reading/writing high-throughput sequencing data. *GigaScience*, 10(2).
  [doi:10.1093/gigascience/giab007](https://doi.org/10.1093/gigascience/giab007)

- **WFA2-lib** &mdash; Marco-Sola, S. et al. (2023). Optimal gap-affine
  alignment in O(s) space. *Bioinformatics*, 39(2).
  [doi:10.1093/bioinformatics/btad074](https://doi.org/10.1093/bioinformatics/btad074)

- **American Gut Project** &mdash; McDonald, D. et al. (2018). American Gut:
  an open platform for citizen science microbiome research. *mSystems*, 3(3).
  [doi:10.1128/mSystems.00031-18](https://doi.org/10.1128/mSystems.00031-18)

- **Earth Microbiome Project** &mdash; Thompson, L.R. et al. (2017). A
  communal catalogue reveals Earth's multiscale microbial diversity. *Nature*,
  551(7681), 457&ndash;463.
  [doi:10.1038/nature24621](https://doi.org/10.1038/nature24621)

- ***Escherichia coli* K-12 MG1655 genome** &mdash; Blattner, F.R. et al.
  (1997). The complete genome sequence of *Escherichia coli* K-12. *Science*,
  277(5331), 1453&ndash;1462.
  [doi:10.1126/science.277.5331.1453](https://doi.org/10.1126/science.277.5331.1453)
