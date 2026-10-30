# WASM Extension Testing

This document describes how to build and verify the miint WASM extension
locally **before** publishing to community extensions. There are four levels
of verification, each catching different classes of bugs.

## Prerequisites

- **emscripten** (emsdk) 6.0.9, the version extension-ci-tools uses for DuckDB
  v2.0, installed and activated (`emcc` in PATH)
- **VCPKG_TOOLCHAIN_PATH** set, or `vcpkg/` symlinked in repo root
- **rustup**. The build installs the pinned toolchain (1.96.0), its
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
   using emscripten — guaranteeing the same ABI as the extension — and
   configures it with the DuckDB version stamped in the extension, so `LOAD`
   accepts it. A core cached in `build/wasm_shell/` keeps the version it was
   built with, and the script stops when that differs from the extension's or
   was not recorded; `--clean` rebuilds it
2. Compiles `scripts/test_wasm_extension.c` as a WASM `MAIN_MODULE` linked
   against DuckDB, with the miint extension embedded via `--preload-file`, and
   registers DuckDB's statically linked extensions through the loader DuckDB's
   `scripts/generate_static_extension_loader.py` writes. Which extensions, and
   which undefined symbols fail the link, are in the script's step 2
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
built with a compatible emsdk. The version must match exactly: a dev build of
the same commit reports `vX.Y.Z-devN`, and DuckDB refuses to `LOAD` an
extension stamped `vX.Y.Z` into it. Level 4 stops at boot, naming the
mismatch, while its pin in `scripts/wasm_browser/index.html` reports another
version. The v2.0.0 runtime at
[shell.duckdb.org](https://shell.duckdb.org) is not an npm release: it serves
its own build with its own worker protocol, which the duckdb-wasm npm API does
not drive.

## Level 4: Browser Check

Catches: Rust code that fails or computes a different answer in the runtime the
playground uses.

```bash
./scripts/test_wasm_browser.sh            # the build/wasm_eh extension
./scripts/test_wasm_browser.sh path/to/miint.duckdb_extension.wasm
```

It boots duckdb-wasm in headless Chrome the way the playground does (the
`@duckdb/duckdb-wasm` release pinned in `scripts/wasm_browser/index.html`,
`selectBundle`, the base bundle, which runs `wasm_eh`,
`allowUnsignedExtensions`), serves the extension with no cross-origin isolation
headers, like GitHub Pages, `LOAD`s it, and runs
`scripts/wasm_browser/steps.json`. Each step with an `expect` must return that
value. The run passes only when every step ran and passed; a mismatch, an
error, a step left unrun because the database died, duckdb-wasm picking a
bundle other than base, or the pin running a DuckDB core other than the
workflow's `duckdb_version` fails it. Timeouts: in the browser, 2 min for boot
and for `LOAD`, 3 min per step and 20 s for the alive check after it; in native
mode, 15 min for the whole run and 256 MiB of shell output.

**This is the check for Rust code paths.** Level 3 links a main module and JS
glue built from the `duckdb/` submodule; this check runs the duckdb-wasm release
users get. Rust calls cross into that release's JS (exception handling, and
function pointers resolved by its loader; see "Exceptions" under Level 2), so
judge them here.

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

On push to a tag, the distribution pipeline builds `wasm_eh` (Level 1); pull
requests and branch pushes do not. The jobs for Levels 2 and 3 (`verify-wasm`,
`test-wasm-load`) and the WASM-compatible native build are commented out in
`.github/workflows/MainDistributionPipeline.yml`, so run Levels 2–4 locally
before a release.

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
the library's `cmake/deps/<name>.cmake` appends it to `MIINT_WASM_LIBS`, which
becomes `DUCKDB_EXTENSION_MIINT_LINKED_LIBS`.

### `relocation R_WASM_MEMORY_ADDR_SLEB cannot be used against symbol`

The library wasn't compiled with `-fPIC`. WASM side modules require PIC for
all linked objects. Check the EMSCRIPTEN blocks in `cmake/deps/` — each
external library needs `-fPIC` in its compile flags.

### `imported function does not match the expected type`

The extension and the DuckDB WASM runtime were built with emsdk versions that
generate different function signatures. A different core version does not get
this far; see the next entry.

### `The file was built specifically for DuckDB version ...`

DuckDB refuses the `LOAD` because the extension was built for another core
version. Levels 3 and 4 check the version before `LOAD` and stop with their own
error: Level 3 when its cached core in `build/wasm_shell/` is another version
(rerun with `--clean`), Level 4 at boot, naming the pin's file to bump.

### `function signature mismatch` from a Rust function

A Rust function reached through a path Emscripten's JS mishandles. The two
known causes need a runtime built without `WASM_BIGINT`, which the v2 runtimes
are not; see "Runtimes built without `WASM_BIGINT`" under "Emscripten exceptions" in
[embedded-tools.md](internals/embedded-tools.md#rust-crates-via-the-miint_rust_glue-umbrella)
for what to check.

### Extension loads but `miint_version()` fails

The extension loaded but initialization threw an error. Common cause: a
dependency extension (like `json`) isn't available. The test harness
pre-loads the json extension to avoid this.
