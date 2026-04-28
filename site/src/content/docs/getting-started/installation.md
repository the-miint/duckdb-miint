---
title: Installation
description: Installing the released MIINT extension into DuckDB; also covers building from source.
---

## Table of Contents

- [Installing the released extension](#installing-the-released-extension)
- [Python CLI](#python-cli)
- [Building from source](#building-from-source)
- [Running a from-source build](#running-a-from-source-build)

## Installing the released extension

MIINT is published in the [DuckDB Community Extensions](https://community-extensions.duckdb.org/)
repository. From any standard DuckDB shell:

```sql
INSTALL httpfs; INSTALL miint FROM community;
LOAD httpfs; LOAD miint;
```

`httpfs` is loaded explicitly because some MIINT code paths fetch remote
resources and `httpfs` does not always autoload. Verify the install:

```sql
SELECT miint_version();
```

That's all most users need. The rest of this page covers building from
source for development.

## Python CLI

A lightweight Python CLI (`miint`) wraps the DuckDB extension for common
bioinformatics workflows: format conversion, sequence alignment, and feature
table generation. All computation happens in SQL — Python handles argument
parsing and query construction.

### Installing from a local checkout

```shell
cd python
pip install -e .
```

This requires the `duckdb` Python package (installed automatically as a
dependency). The CLI is **not** published to PyPI at this time.

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

## Building from source

This is for contributors and for environments that need a custom build.
End users should prefer the released-extension flow above.

### Managing dependencies

**Rust toolchain (required):** The RYpe sequence classification library is
written in Rust. Install Rust via [rustup](https://rustup.rs/):
```shell
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

**VCPKG (required):** DuckDB extensions use VCPKG for dependency
management. Enabling VCPKG is very simple: follow the
[installation instructions](https://vcpkg.io/en/getting-started) or just
run the following:
```shell
git clone https://github.com/Microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh
export VCPKG_TOOLCHAIN_PATH=`pwd`/vcpkg/scripts/buildsystems/vcpkg.cmake
```

### Build system
We use [Ninja](https://ninja-build.org) for quick builds. The easiest
install is to use prebuilt [release](https://github.com/ninja-build/ninja/releases)
binaries.

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

- `duckdb` is the binary for the DuckDB shell with the extension code
  automatically loaded — no further `LOAD` is needed.
- `unittest` is the DuckDB test runner with the extension linked in (for
  SQL logic tests).
- `miint.duckdb_extension` is the loadable binary as it would be
  distributed.
- `tests` is the Catch2 C++ unit test runner.

## Running a from-source build

The simplest path is to use the bundled binary that already has miint
linked in:

```sh
./build/release/duckdb
```

No `LOAD` and no `-unsigned` flag are needed — miint is already
available.

If you want to load the freshly-built `.duckdb_extension` into a stock
DuckDB shell of the matching version (e.g. for testing the loadable
artifact), launch with `-unsigned` since the local build is unsigned:

```sh
duckdb -unsigned -c "LOAD '$(pwd)/build/release/extension/miint/miint.duckdb_extension';"
```

The stock-shell + LOAD path is sometimes finicky depending on how DuckDB
was packaged on your system; if you hit signature errors that `-unsigned`
doesn't resolve, use `./build/release/duckdb` instead.
