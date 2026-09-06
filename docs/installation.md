# Installation and Building

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

### System prerequisites

These are the system-level tools the build invokes directly. CI installs the same set via `duckdb/extension-ci-tools/.github/workflows/_extension_distribution.yml` (with `extra_toolchains: 'rust;omp'`); see `.github/workflows/MainDistributionPipeline.yml` for the upstream reference.

**macOS (Homebrew):**
```shell
brew install ninja pkg-config autoconf automake libtool autoconf-archive libomp
```
- `ninja` — build driver (`bash build.sh` runs `GEN=ninja make`).
- `pkg-config` — required by vcpkg's `openssl` port (transitively pulled in by `curl` / `rocksdb`).
- `autoconf`, `automake`, `libtool`, `autoconf-archive` — vsearch ships `autogen.sh`; without these the configure step fails with `autoreconf: command not found`.
- `libomp` — required when `MIINT_ENABLE_UNIFRAC=ON` (the default). Apple Clang doesn't ship libomp, and the CMake check looks for Homebrew's `libomp.a` at `/opt/homebrew/opt/libomp` (Apple Silicon) or `/usr/local/opt/libomp` (Intel). Pass `-DMIINT_ENABLE_UNIFRAC=OFF` to build without UniFrac and skip this.

**Linux (Debian/Ubuntu):**
```shell
sudo apt-get install -y ninja-build pkg-config autoconf automake libtool autoconf-archive
```
OpenMP comes from `libgomp` shipped with gcc, so no separate install is needed for UniFrac.

**Windows:** see CI for the reference setup. The `windows_amd64` MSVC and `windows_amd64_rtools` builds are currently excluded from CI (see `exclude_archs` in `MainDistributionPipeline.yml`); `windows_amd64_mingw` is the supported path.

### Managing dependencies

**Rust toolchain (required):** The RYpe sequence classification library, and the embedded sylph profiler, are written in Rust. Install Rust via [rustup](https://rustup.rs/):
```shell
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

**VCPKG (required):** DuckDB extensions use VCPKG for dependency management. Enabling VCPKG is very simple: follow the [installation instructions](https://vcpkg.io/en/getting-started) or just run the following:
```shell
git clone https://github.com/Microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh
export VCPKG_TOOLCHAIN_PATH=`pwd`/vcpkg/scripts/buildsystems/vcpkg.cmake
```

**Submodules (required):** Several embedded libraries (MAFFT, vsearch, sylph, scikit-bio-binaries, unifrac-binaries, etc.) are git submodules. After cloning:
```shell
git submodule update --init --recursive
```

### Build system
We use [Ninja](https://ninja-build.org) for quick builds. The easiest install is to use prebuilt [release](https://github.com/ninja-build/ninja/releases) binaries, or `brew install ninja` / `apt-get install ninja-build` as listed above.

### Optional feature flags

Pass these to CMake via `EXT_FLAGS` (e.g. `EXT_FLAGS="-DMIINT_ENABLE_UNIFRAC=OFF" bash build.sh`) to opt out of optional subsystems:

| Flag | Default | What it controls |
|------|---------|------------------|
| `MIINT_ENABLE_CURL` | `ON` | libcurl-based streaming uploads (auto-disabled on macOS + WASM) |
| `MIINT_ENABLE_HDF5` | `ON` | HDF5/BIOM reader and writer |
| `MIINT_ENABLE_MAFFT` | `ON` | MAFFT multiple sequence alignment |
| `MIINT_ENABLE_SORTMERNA` | `ON` | SortMeRNA rRNA alignment (requires RocksDB) |
| `MIINT_ENABLE_SYLPH` | `ON` | sylph FracMinHash relative-abundance profiling (Rust) |
| `MIINT_ENABLE_GPL_BOUNDARY` | `ON` | gpl-boundary subsystem (process pipes + Arrow IPC over POSIX shm) |
| `MIINT_ENABLE_UNIFRAC` | `ON` | UniFrac (PCoA / PERMANOVA / Faith PD) — requires libomp on macOS |
| `MIINT_ENABLE_KREPP` | `ON` | krepp phylogenetic placement (`place_krepp`) — requires the parallel-hashmap and boost-math headers from vcpkg |

WASM (Emscripten) builds automatically disable HDF5, SortMeRNA, sylph, gpl-boundary, krepp, vsearch, and libcurl. UniFrac and MAFFT remain enabled — UniFrac builds against a dedicated single-threaded WASM target (`libssu_wasm.a` / `libskbb_wasm.a`, Eigen-backed, OpenMP stubbed out).

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
