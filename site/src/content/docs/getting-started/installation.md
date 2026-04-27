---
title: Installation
description: Building duckdb-miint from source and loading it into DuckDB.
---

## Table of Contents

- [Installing](#installing)
- [Python CLI](#python-cli)
- [Building](#building)
- [Running the extension](#running-the-extension)

## Installing

MIINT is available in the [DuckDB Community Extensions](https://community-extensions.duckdb.org/) repository:

```sql
INSTALL miint FROM community;
LOAD miint;
```

You can verify the installed version with:
```sql
SELECT miint_version();
```

### Installing from a local build

If you are building MIINT from source (see [Building](#building)), launch DuckDB with the `allow_unsigned_extensions` option:

CLI:
```shell
duckdb -unsigned
```

Python:
```python
con = duckdb.connect(':memory:', config={'allow_unsigned_extensions' : 'true'})
```

NodeJS:
```js
db = new duckdb.Database(':memory:', {"allow_unsigned_extensions": "true"});
```

Then load the extension binary directly:
```sql
LOAD '/path/to/miint.duckdb_extension';
```

## Python CLI

A lightweight Python CLI (`miint`) wraps the DuckDB extension for common bioinformatics workflows: format conversion, sequence alignment, and feature table generation. All computation happens in SQL — Python handles argument parsing and query construction.

### Installing from a local checkout

```shell
cd python
pip install -e .
```

This requires the `duckdb` Python package (installed automatically as a dependency). The CLI is **not** published to PyPI at this time.

### Usage

```shell
# Use a local extension build
miint --extension-path ./build/release/extension/miint/miint.duckdb_extension \
    convert sequence -1 reads.fastq.gz -o reads.parquet

# Without --extension-path, the CLI installs miint from community extensions
miint convert alignment -i alignments.bam -o alignments.parquet
```

### Available commands

```
miint convert sequence    Convert FASTQ/FASTA to parquet
miint convert alignment   Convert SAM/BAM to parquet
miint convert biom        Convert BIOM to parquet
miint convert parquet     Convert parquet to FASTQ/FASTA/SAM/BAM/BIOM

miint transform genome-coverage   Compute genome coverage from alignments
miint transform woltka-ogu        Generate OGU feature table (with optional filters)

miint align minimap2              Align sequences with minimap2
miint align minimap2-sharded      Sharded minimap2 alignment with RYpe classification
```

Run `miint --help` or `miint <command> --help` for detailed usage.

## Building

### Managing dependencies

**Rust toolchain (required):** The RYpe sequence classification library is written in Rust. Install Rust via [rustup](https://rustup.rs/):
```shell
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

**VCPKG (required):** DuckDB extensions use VCPKG for dependency management. Enabling VCPKG is very simple: follow the [installation instructions](https://vcpkg.io/en/getting-started) or just run the following:
```shell
git clone https://github.com/Microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh
export VCPKG_TOOLCHAIN_PATH=`pwd`/vcpkg/scripts/buildsystems/vcpkg.cmake
```

### Build system
We use [Ninja](https://ninja-build.org) for quick builds. The easiest install is to use prebuilt [release](https://github.com/ninja-build/ninja/releases) binaries.

### Build steps
Now to build the extension, run:
```sh
bash build.sh
```

This is equivalent to:
```sh
export VCPKG_TOOLCHAIN_PATH=`pwd`/vcpkg/scripts/buildsystems/vcpkg.cmake
GEN=ninja make
```

The main binaries that will be built are:

```sh
./build/release/duckdb
./build/release/test/unittest
./build/release/extension/miint/miint.duckdb_extension
./build/release/extension/miint/tests
```

- `duckdb` is the binary for the DuckDB shell with the extension code automatically loaded.
- `unittest` is the DuckDB test runner with the extension linked in (for SQL logic tests).
- `miint.duckdb_extension` is the loadable binary as it would be distributed.
- `tests` is the Catch2 C++ unit test runner.

## Running the extension
To run the extension code, simply start the shell with `./build/release/duckdb`.
