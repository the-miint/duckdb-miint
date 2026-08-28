# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## ABSOLUTE HARD REQUIREMENTS
- NEVER use `rm` without permission

## Project Overview

**duckdb-miint** is a DuckDB extension that brings columnar analytics to microbiome research. It enables SQL-based analysis of microbiome sequence (FASTA/FASTQ), alignment (SAM/BAM/MSA), phylogeny (Newick), and mass spectrometry (mzML, mzXML) data. It provides powerful, and correct, mechanisms to manipulate and store these data.

The extension provides:
- Functions to read and write common formats.
- Functions to perform important operations on data such as query coverage filtering.
- Embedded libraries to expose critical capabilities such as high throughput alignment.

## Priorities

1. Red/green/refactor Test Driven Development (TDD)
2. Verifiably correct code
3. Maintainable code, using Don't Repeat Yourself (DRY) and Keep It Simple Stupid (KISS)
4. Performance

## Rules

These rules apply to every task in this project unless explicitly overridden.

### Rule 1 — Think Before Coding
Bias: caution over speed on non-trivial work.
State assumptions explicitly. Ask rather than guess.
Push back when a simpler approach exists. Stop when confused.

### Rule 2 — Simplicity First
Minimum code that solves the problem. Nothing speculative.
No abstractions for single-use code.

### Rule 3 — Surgical Changes
Touch only what you must. Don't improve adjacent code.
Match existing style. Don't refactor what isn't broken.

### Rule 4 — Goal-Driven Execution
Define success criteria. Loop until verified.

### Rule 5 — Surface conflicts, don't average them
If two patterns contradict, pick one (more recent / more tested).
Explain why. Flag the other for cleanup.

### Rule 6 — Read before you write
Before adding code, read exports, immediate callers, shared utilities.
If unsure why existing code is structured a certain way, ask.

### Rule 7 — Tests verify intent, not just behavior
Tests must encode WHY behavior matters, not just WHAT it does.
A test that can't fail when business logic changes is wrong.

### Rule 8 — Checkpoint after every significant step
Summarize what was done, what's verified, what's left.
Don't continue from a state you can't describe back.

### Rule 9 — Match the codebase's conventions, even if you disagree
Conformance > taste inside the codebase.
If you think a convention is harmful, surface it. Don't fork silently.

### Rule 10 — Fail loud
"Completed" is wrong if anything was skipped silently.
Don't silently skip tests you caused to be skipped (commented out, xfail'd, missed) — existing `require-env` skips for optional deps/oracles are legitimate and expected.
Default to surfacing uncertainty, not hiding it.

## Build Commands

### Initial Setup (VCPKG for dependencies)
```bash
git clone https://github.com/Microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh
export VCPKG_TOOLCHAIN_PATH=`pwd`/vcpkg/scripts/buildsystems/vcpkg.cmake
```

### Build the extension
```bash
bash build.sh
```

This creates:
- `./build/release/duckdb` — DuckDB shell with extension preloaded
- `./build/release/test/unittest` — DuckDB test runner with extension linked
- `./build/release/extension/miint/miint.duckdb_extension` — Loadable extension binary
- `./build/release/tests` — C++ unit tests (Catch2)

### Development workflow
```bash
make clean            # clean build
bash build.sh         # always use build.sh, not make directly
```

### Formatting
Pre-commit hook (and CI's "Code Quality Check / Format Check" job) runs `make format-check`. CI pins specific versions — match them locally or newer `black` will reformat files differently than CI expects and the check will fail on files you never touched.

Required versions (from `extension-ci-tools/.github/workflows/_extension_code_quality.yml`):
- `black==24.*`
- `clang_format==11.0.1`
- `cmake-format` (any recent version)

One-time setup of a dedicated conda env:
```bash
conda create -n duckdb-format -c conda-forge python=3.12 pip -y
conda activate duckdb-format
pip install 'black==24.*' 'clang_format==11.0.1' cmake-format
```

Then to check / auto-fix:
```bash
conda activate duckdb-format
make format-check                  # verify (CI runs this)
make format-fix                    # auto-fix differences
```

Note: black ships breaking style changes between major versions — pinning to `24.*` is required, not optional.

## Testing

If a test produces an **incorrect expected value**: DO NOT change the expected value without permission.

```bash
bash run_tests.sh                                                 # all SQL tests
./build/release/extension/miint/tests                             # C++ unit tests
./build/release/test/unittest "test/sql/read_alignments.test"     # single SQL test
./build/release/test/unittest "[alignment]"                       # tests matching a pattern
```

### Python CLI smoke tests
```bash
cd python && conda run -n duckdb-144 pip install -e . && cd ..
EXT=./build/release/extension/miint/miint.duckdb_extension
conda run -n duckdb-144 miint --extension-path $EXT convert sequence \
  -1 data/fastq/small_a.fq -o /tmp/test.parquet
```
Full manual smoke test procedure: `localdocs/cli-smoke-tests.md`.

### Optional-dependency guards in SQL tests
SQL tests use `require-env VAR_NAME` (after `require miint`) to skip gracefully when a dependency isn't present. If the env var isn't set, the test file is skipped, not failed.

`run_tests.sh` is the authoritative source for what gets detected and exported. It covers three kinds of guards:
- **Runtime binaries** detected via `command -v` (e.g., `BOWTIE2_AVAILABLE`, `ASPERA_AVAILABLE`, `MASSQL_PYTHON_AVAILABLE`, `GPL_BOUNDARY_AVAILABLE` for the FastTree process-isolation host, `FASTTREE_AVAILABLE` for the bioconda parity oracle binary)
- **Compile-time features** detected by querying the built extension (e.g., `HDF5_AVAILABLE`, `VSEARCH_AVAILABLE`, `MAFFT_AVAILABLE` — each checks for a registered function or library entry)
- **Downloaded / served test data** (e.g., `MIINT_HTTPS_TEST_URL`, `MASSQL_TEST_DATA`, `MASSQL_GNPS_DATA`, `MZXML_REAL_DATA`)
- **SHA-pinned parity oracles** (e.g., `SORTMERNA_REAL_ORACLE`, `MIINT_FASTTREE_TINY_PARITY_OK`, `MIINT_FASTTREE_MODERATE_PARITY_OK`) — exported only when the oracle file's SHA matches the recorded sidecar. To regenerate stale oracles: `MIINT_SORTMERNA_REAL_DATA=1` or `MIINT_FASTTREE_REGENERATE=1` and rerun `run_tests.sh`.
- **Process state the main pass cannot provide** (currently `MIINT_LOW_FD_TEST`, guarding `test/sql/read_fastx_fd_release.test`) — not an availability check. `make test` has no way to set `RLIMIT_NOFILE`, so the file opts out of that pass and `run_tests.sh` reruns it in a subshell with the limit lowered. Reach for this only when the test needs process state, never to work around missing data.

When adding a new guard: detect in `run_tests.sh`, add `require-env` to the test file(s), and leave availability-check tests (that only verify the scalar/feature-flag query itself) unguarded so they always run.

## Architecture and Patterns (detailed docs)

The entry point is `src/miint_extension.cpp` — `LoadInternal()` registers every table function, scalar function, aggregate, and COPY format. That file is the authoritative catalog.

Developer-facing deep dives live under `docs/internals/`:

- **[`docs/internals/architecture.md`](docs/internals/architecture.md)** — design patterns (file reading, record abstraction, reference table), code style, testing strategy, cross-cutting impl details (thread safety, headerless SAM, stop-position math, quality scores, compression), and how to add new table/COPY/scalar/aggregate functions — including **[why the work goes in `InitGlobal`/`Execute` and never in `Bind`](docs/internals/architecture.md#no-work-in-bind)**, which is easy to get wrong because it appears to work.
- **[`docs/internals/embedded-tools.md`](docs/internals/embedded-tools.md)** — how every external library/tool is embedded: static libraries from source (HTSlib, minimap2, WFA2, vsearch, MAFFT, rype), header-only (kseq++), vcpkg/system (zlib, zstd, expat, HDF5, Catch2), and runtime binaries (bowtie2, Aspera). Platform-specific gotchas and feature flags.
- **[`docs/internals/reading-tables-views.md`](docs/internals/reading-tables-views.md)** — the separate-connection recipe for reading user-specified tables/views from extension code (avoids the `context.Query()` deadlock). Covers both data reads and schema validation.
- **[`docs/internals/arrow-zero-copy.md`](docs/internals/arrow-zero-copy.md)** — zero-copy Arrow C Data Interface → DuckDB Vector conversion, with lifetime management and reference implementations.
- **[`docs/internals/per-sample-pattern.md`](docs/internals/per-sample-pattern.md)** — the `sample_id` named parameter that partitions an input relation and runs a function's pipeline once per distinct sample, prepending a sample column. Covers `src/include/per_sample_table_function.hpp` (bind-time discovery, exec-time atomic claim) and which functions use it.
- **[`docs/internals/duckdb-engine-notes.md`](docs/internals/duckdb-engine-notes.md)** — non-obvious DuckDB engine behavior that affects us: storage version 64 as the default (and the compression it gates), non-flat input vectors, `string_t` inlining limits, single-pass optimizer passes, zonemap vs. ART index reality, and jemalloc knobs.

User-facing API reference (parameters, return types, examples) is organized by task under `docs/` — start from the index at [`docs/table_of_contents.md`](docs/table_of_contents.md) (e.g. `docs/reading.md`, `docs/writing.md`, `docs/alignment_*.md`, `docs/qc.md`, `docs/diversity.md`, `docs/insdc_ena.md`, `docs/utilities.md`).

