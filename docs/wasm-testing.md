# WASM Extension Testing

This document describes how to build and verify the miint WASM extension
locally **before** publishing to community extensions. There are four levels
of verification, each catching different classes of bugs.

## Prerequisites

- **emscripten** (emsdk) installed and activated (`emcc` in PATH)
- **VCPKG_TOOLCHAIN_PATH** set, or `vcpkg/` symlinked in repo root
- **rustup**. The build installs the pinned toolchain (1.86.0), its
  `wasm32-unknown-emscripten` target and, for `wasm_eh`, its `rust-src`
  component; configure stops without rustup
- **Node.js** for Levels 2 and 3; Level 4 needs **Node.js 22.12 or newer** and
  **npm**
- **Google Chrome** or Chromium, and network access to cdn.jsdelivr.net, for
  Level 4

## Quick Start

```bash
# Full pipeline: build extension, verify imports, load test, browser check
./scripts/build_wasm.sh
./scripts/test_wasm_load.sh
./scripts/test_wasm_browser.sh
```

## Level 1: Build Succeeds

Catches: linker errors, missing symbols, `-fPIC` issues, Rust PIC problems.

```bash
./scripts/build_wasm.sh --build-only
```

This runs `make wasm_eh` which builds the extension as a WASM side module.
The output is `build/wasm_eh/extension/miint/miint.duckdb_extension.wasm`.

Use `--clean` for a full rebuild (necessary after changing CMake flags or
updating submodules):

```bash
./scripts/build_wasm.sh --clean --build-only
```

### What it verifies

- All C libraries (htslib, minimap2, WFA2, vsearch) compile for WASM
- Rust rype library compiles with PIC for WASM
- The emcc post-build linking step succeeds (all static libraries linked
  into the side module via `DUCKDB_EXTENSION_MIINT_LINKED_LIBS`)
- A `.wasm` file is produced

The Catch2 `tests` target (`add_executable(tests …)` in `CMakeLists.txt`) is not gated on
Emscripten, so `test/cpp/` is held to the same rule as `src/`: a POSIX-only call outside a
`MIINT_ENABLE_*` or `__EMSCRIPTEN__` guard breaks the WASM build. Pull-request CI does not
build `wasm_eh`; it builds on tag pushes, so such a break first appears at release unless
this level is run locally. The second C++ target, `tests_duckdb`, is the exception: DuckDB
builds no `unittest` for WASM, and `tests_duckdb` is built only where it does, so its cases
never run in a WASM build.

### What it does NOT catch

- Unresolved symbols that become WASM imports (deferred to load time)
- ABI mismatches between extension and DuckDB runtime
- Runtime errors in extension code

## Level 2: Import Verification

Catches: unresolved library symbols that would fail at load time (like the
`gzclose` bug).

```bash
./scripts/build_wasm.sh    # builds + verifies (default)
```

This uses Node.js to compile the `.wasm` file and inspect its imports/exports.
It checks that no symbols from our static libraries (zlib, htslib, minimap2,
WFA2, expat, rype) leak as unresolved imports.

### How PIC side modules work

In a PIC WASM side module (`-sSIDE_MODULE=2`), library symbols appear as:
- `GOT.func.X` / `GOT.mem.X` **imports** (PIC mechanism)
- Direct **exports** of `X` (the actual symbol definition)

A GOT import with a matching export is correct — the emscripten runtime
resolves it at load time by pointing the GOT entry at the module's own export.
A GOT import **without** a matching export means the library wasn't linked.

The verification script checks:
1. No direct `env.*` imports of library symbols (= not linked at all)
2. No `GOT.*` imports without matching exports (= linked but PIC broken)
3. `miint_duckdb_cpp_init` export exists (extension entry point)
4. In a build with wasm-native exceptions (`wasm_eh`), no imports of
   Emscripten's JS exception-handling helpers (`invoke_*`,
   `__cxa_find_matching_catch_*`, `__resumeException`). The build records its
   exception model as `MIINT_WASM_NATIVE_EH` in `CMakeCache.txt`; a build that
   does not record it fails the check rather than skipping it

### Exceptions

The [Rust glue](internals/embedded-tools.md#rust-crates-via-the-miint_rust_glue-umbrella)
uses the C++ side's exception model, which on `wasm_eh` means wasm-native
exceptions, and check 4 fails a build that still imports the JS helpers. Its
"Emscripten exceptions" entry says why, and what that does not fix.

## Level 3: Headless Load Test

Catches: ABI mismatches against the `duckdb/` submodule's core, and extension
initialization failures.

```bash
./scripts/test_wasm_load.sh
```

It:

1. Builds DuckDB WASM (`libduckdb_static.a`) from the `duckdb/` submodule
   using emscripten — guaranteeing the same ABI as the extension
2. Compiles `scripts/test_wasm_extension.c` as a WASM `MAIN_MODULE` linked
   against DuckDB, with the miint extension embedded via `--preload-file`
3. Runs the resulting `test_wasm_harness.js` in Node.js
4. Verifies: DuckDB opens, json extension loads, miint extension loads,
   `miint_version()` returns successfully, and UniFrac's `faith_pd` computes

Its main module and JS glue are built here, not taken from the duckdb-wasm
release users run, and it runs only the queries in
`scripts/test_wasm_extension.c`. For how the extension behaves in that release,
use Level 4.

### First run

The first run builds DuckDB WASM from source (~5 minutes). Subsequent runs
reuse the cached build:

```bash
# First time: builds DuckDB + harness + runs test
./scripts/test_wasm_load.sh

# After code changes to extension only (DuckDB core unchanged):
./scripts/build_wasm.sh --build-only    # rebuild extension
./scripts/test_wasm_load.sh --harness-only  # recompile harness + run test

# Full clean rebuild:
./scripts/test_wasm_load.sh --clean
```

### Why build DuckDB here rather than use the npm package?

A C++ extension loads only into a DuckDB WASM of the same core version, built
with a compatible emsdk. A core mismatch stops at `LOAD` with DuckDB's own error
that the extension was built for another version; builds from emsdk versions
that generate different function signatures fail with
`LinkError: imported function does not match the expected type`. Building the
core from the submodule matches both by construction. An npm
`@duckdb/duckdb-wasm` release matches when it reports the same core and was
built with a compatible emsdk; the release the playground loads is one, since
the extension built here with emsdk 3.1.71 loads into it. Level 4 uses that
release, pinned in `scripts/wasm_browser/index.html`, and stops at boot when its
core differs.

## Level 4: Browser Check

Catches: Rust code that fails or computes a different answer in the runtime the
playground uses.

```bash
./scripts/test_wasm_browser.sh            # the build/wasm_eh extension
./scripts/test_wasm_browser.sh path/to/miint.duckdb_extension.wasm
```

It boots duckdb-wasm in headless Chrome the way the playground does (the
`@duckdb/duckdb-wasm` release pinned in `scripts/wasm_browser/index.html`,
`selectBundle`, the eh bundle, `allowUnsignedExtensions`), serves the extension
with no cross-origin isolation headers, like GitHub Pages, `LOAD`s it, and runs
`scripts/wasm_browser/steps.json`. Each step with an `expect` must return that
value. The run passes only when every step ran and passed; a mismatch, an
error, a step left unrun because the database died, duckdb-wasm picking a
bundle other than eh, or the pin running a DuckDB core other than the
workflow's `duckdb_version` fails it. Timeouts: in the browser, 2 min for boot
and for `LOAD`, 3 min per step and 20 s for the alive check after it; in native
mode, 15 min for the whole run and 256 MiB of shell output.

**This is the check for Rust code paths.** Level 3 links a main module and JS
glue built from the `duckdb/` submodule; this check runs the duckdb-wasm release
users get. Rust code depends on that release's JS (see "Exceptions" under
Level 2), so judge it here.

### Adding steps

A step has a unique `name` and its `sql`, and runs in order on one connection.
It checks one of three things:

- `expect`: the first cell of its first row equals this. Make the step return
  one VARCHAR cell (an `md5` over an ordered `string_agg` for anything bigger),
  so both runtimes print it the same way.
- `expect_error`: it fails with an error containing this text.
- neither (setup): it succeeds.

After every step, in both modes, the database must still answer
`SELECT 'alive'`; a step after which it does not fails, and the run stops there.
The first committed step, `self_check_expected_error`, raises an error on
purpose, so every run goes through that path.

Take `expect` from a native build, which runs the same steps:

```bash
./scripts/test_wasm_browser.sh --native   # build/release/duckdb
```

The rype steps pass whether the Rust glue is built with JS or wasm-native
exceptions, so check 4 in Level 2 is what fails a build with JS exceptions.

## CI Integration

Two CI jobs in `.github/workflows/MainDistributionPipeline.yml`:

- **`verify-wasm`** (Level 2) — Downloads the extension `.wasm` artifacts
  from the distribution build and checks for unresolved library imports.
  Runs in ~30 seconds.

- **`test-wasm-load`** (Level 3) — Builds the extension and DuckDB WASM
  from source (matching emsdk 3.1.71 + Rust 1.86.0 + vcpkg), compiles the
  test harness, and runs the headless load test. This is a full independent
  build that catches ABI mismatches. Runs in ~10-15 minutes.

## Build Artifacts

These files are generated by the scripts and should **not** be committed:

- `scripts/test_wasm_harness.js` — emscripten-generated JS loader
- `scripts/test_wasm_harness.wasm` — compiled test harness WASM binary
- `scripts/test_wasm_harness.data` — preloaded extension data
- `scripts/wasm_browser/node_modules/` — puppeteer-core, installed by `test_wasm_browser.sh`
- `build/wasm_eh/` — WASM extension build directory
- `build/wasm_shell/` — DuckDB WASM build directory

## Troubleshooting

### `bad export type for 'X': undefined`

A library symbol wasn't linked into the WASM side module. Check that
`DUCKDB_EXTENSION_MIINT_LINKED_LIBS` in CMakeLists.txt includes the library.

### `relocation R_WASM_MEMORY_ADDR_SLEB cannot be used against symbol`

The library wasn't compiled with `-fPIC`. WASM side modules require PIC for
all linked objects. Check the EMSCRIPTEN blocks in CMakeLists.txt — each
external library needs `-fPIC` in its compile flags.

### `imported function does not match the expected type`

The extension and the DuckDB WASM runtime were built with emsdk versions that
generate different function signatures. A different core version does not get
this far; see the next entry.

### `The file was built specifically for DuckDB version ...`

DuckDB refuses the `LOAD` because the extension was built for another core
version. In Level 3, `build/wasm_shell/` is left over from an older DuckDB;
rerun with `--clean`. Level 4 stops before `LOAD`, at boot, naming the pin's file
to bump.

### `function signature mismatch` from a Rust function

A Rust function reached through a path Emscripten's JS mishandles: see
"Emscripten exceptions" in
[embedded-tools.md](internals/embedded-tools.md#rust-crates-via-the-miint_rust_glue-umbrella).

### Extension loads but `miint_version()` fails

The extension loaded but initialization threw an error. Common cause: a
dependency extension (like `json`) isn't available. The test harness
pre-loads the json extension to avoid this.
