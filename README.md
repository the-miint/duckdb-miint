# MIINT: MIcrobiome INTelligence.

The MIINT [DuckDB](https://duckdb.org) extension enables DuckDB to interoperate with file formats and operations central to microbiome research.

## Quick Start

```sql
-- Install and load the extension
INSTALL miint FROM community;
LOAD miint;

-- Read alignment files using glob pattern (files sorted alphabetically)
-- Use include_filepath to track source file as sample identifier
CREATE VIEW all_alignments AS
    SELECT *, regexp_extract(filepath, '.*/(.*)\.bam', 1) AS sample_id
    FROM read_alignments('samples/*.bam', include_filepath=true);

-- Filter high-quality primary alignments
CREATE VIEW high_quality AS
    SELECT * FROM all_alignments
    WHERE alignment_is_primary(flags)
      AND mapq >= 30
      AND alignment_query_coverage(cigar) > 0.9;

-- Compute OGU counts using Woltka algorithm
CREATE VIEW ogu_counts AS
    SELECT * FROM woltka_ogu_per_sample(high_quality, sample_id, read_id);

-- Export to BIOM format for downstream analysis (QIIME2, phyloseq, etc.)
COPY (SELECT * FROM ogu_counts)
TO 'ogu_table.biom' (FORMAT BIOM, COMPRESSION 'gzip');

-- Or analyze directly in SQL
SELECT sample_id,
       COUNT(DISTINCT feature_id) AS richness,
       SUM(value) AS total_reads
FROM ogu_counts
GROUP BY sample_id;
```

## Key Capabilities

- **Read bioinformatics formats as tables**: FASTQ/FASTA, SAM/BAM, SFF, BIOM, GFF, Newick trees, jplace, mzML, mzXML
- **Mass spectrometry analysis**: MassQL query language, `read_mzml`/`read_mzxml` with `UNION ALL`-compatible schemas, 30+ helper macros for peak matching, isotope patterns, and cross-level queries
- **Align sequences in SQL**: minimap2 and Bowtie2 integration with sharded parallel alignment
- **Classify sequences**: RYpe minimizer-based sequence classification
- **Analyze with SQL**: Filter, aggregate, join with metadata tables
- **Write back to standard formats**: Export results to FASTQ, FASTA, SAM/BAM, BIOM, Newick
- **Powerful functions**: Sequence identity, coverage, flag checking, Woltka classification, pairwise alignment, chemical formula parsing
- **Performance**: Parallel processing, efficient compression, streaming I/O

## Documentation

| Document | Description |
|----------|-------------|
| [Installation & Building](docs/installation.md) | Installing from community extensions, building from source, dependencies |
| [Table Functions](docs/table-functions.md) | `read_alignments`, `read_fastx`, `read_mzml`, `read_mzxml`, `read_biom`, `align_minimap2`, `align_bowtie2`, and more |
| [Mass Spectrometry & MassQL](docs/massql.md) | MassQL query language, `read_mzml`/`read_mzxml`, helper macros, formula functions |
| [Scalar Functions](docs/scalar-functions.md) | SAM flag functions, sequence identity, query length, query coverage |
| [RYpe Functions](docs/rype.md) | RYpe sequence classification and minimizer extraction |
| [Analysis Functions](docs/analysis-functions.md) | Woltka OGU, reverse complement, IUPAC regexp, interval compression, pairwise alignment, `formula()`, `miint_version()` |
| [COPY Formats](docs/copy-formats.md) | Writing FASTQ, FASTA, SAM, BAM, BIOM, and Newick files |
| [Testing](docs/testing.md) | SQL logic tests, C++ unit tests, shell tests, test data |

## Installing

```sql
INSTALL miint FROM community;
LOAD miint;
```

See [Installation & Building](docs/installation.md) for building from source.

## Python CLI

A lightweight Python CLI wraps the extension for common workflows (format conversion, alignment, feature table generation):

```sh
# Install from a local checkout (requires duckdb Python package)
cd python && pip install -e .

# Convert FASTQ to parquet using a local extension build
miint --extension-path ./build/release/extension/miint/miint.duckdb_extension \
    convert sequence -1 reads.fastq.gz -o reads.parquet

# Align and compute OGU feature table
miint align minimap2 -1 reads.parquet -d reference.fasta.gz -o alignments.parquet
miint transform woltka-ogu -i alignments.parquet -o feature_table.parquet
```

See `miint --help` for all commands. Without `--extension-path`, the CLI installs miint from community extensions automatically.

## Running the Tests

```sh
# All tests (SQL, C++, and shell)
bash run_tests.sh

# SQL tests only
make test

# C++ unit tests only
./build/release/extension/miint/tests
```

See [Testing](docs/testing.md) for details.

---

This repository is based on https://github.com/duckdb/extension-template, check it out if you want to build and ship your own DuckDB extension.
