# Embedded Tools and Dependencies

How duckdb-miint embeds external libraries and tools. `CMakeLists.txt` is the source of truth; this file explains the *why* and the non-obvious gotchas.

Five embedding categories:
1. **Static libraries built from source** via `ExternalProject_Add`, linked into the extension
2. **Header-only libraries** included directly
3. **System/vcpkg libraries** via `find_package`
4. **Runtime binaries** invoked via `fork`/`exec` (not compiled in)
5. **Vendored source** compiled inline (a `.c`/`.h` pair added to the source lists)

## Feature Flags (CMake Options)

| Flag | Default | Auto-disabled when |
|---|---|---|
| `MIINT_ENABLE_CURL` | ON | Apple (vsearch/OpenSSL `MD5_*`/`SHA1_*` symbol clash), Emscripten (no network stack, and a loadable extension is a wasm *side module*, so undefined `curl_*` symbols become imports duckdb-wasm cannot supply — they trap when **called**, not at load). Also soft-disabled when neither the vcpkg CONFIG package nor the system `FindCURL` module locates libcurl, so the no-vcpkg tidy lane still configures. |
| `MIINT_ENABLE_HDF5` | ON | Emscripten (C++ static class members become unresolvable GOT.mem imports in WASM) |
| `MIINT_ENABLE_MAFFT` | ON | Windows (uses `mkdtemp` and other POSIX APIs; segfaults on MinGW) |
| `MIINT_ENABLE_ABPOA` | ON | Windows (POSIX APIs) |
| `MIINT_ENABLE_VSEARCH` | ON | Emscripten, Windows (autotools build not supported) |
| `MIINT_ENABLE_SORTMERNA` | ON | Emscripten (RocksDB vcpkg port not built for wasm32), Windows/MinGW (cmph assumes POSIX `<sys/time.h>`; MSVC-on-Windows would work if anyone wires it up) |
| `MIINT_ENABLE_GPL_BOUNDARY` | ON | Emscripten, Windows (subsystem uses POSIX shm + fork/exec) |
| `MIINT_ENABLE_UNIFRAC` | ON | Windows (libssu's inmem build assumes POSIX; first-class on Emscripten via the WASM target) |
| `MIINT_ENABLE_SYLPH`   | ON | Emscripten (WASM-incompatible), Windows/MinGW (POSIX-only sketch indexing) |

Corresponding preprocessor macros: `MIINT_HAS_CURL`, `MIINT_HAS_HDF5`, `MIINT_HAS_MAFFT`, `MIINT_HAS_ABPOA`, `MIINT_HAS_VSEARCH`, `MIINT_HAS_SORTMERNA`, `MIINT_HAS_GPL_BOUNDARY`, `MIINT_HAS_UNIFRAC`, `MIINT_HAS_SYLPH`. Also `MIINT_ASPERA_SUPPORTED=0` on Windows/WASM (POSIX-only runtime).

`MIINT_HAS_AVX2` and `MIINT_HAS_AVX512` are different in kind: not user-facing options but the results of a configure-time **compiler probe** (mirroring htslib's `hts_probe_cc.sh`), so an old toolchain builds fewer kernel variants instead of failing. Set only on x86-64 and never under Emscripten. See "Per-instruction-set compilation of the mmvec kernel" in §3.

Run-time / conditional: `MIINT_USE_JEMALLOC` is set when DuckDB's jemalloc is linked (not on musl/macOS/Windows). `MIINT_SIMD` selects the mmvec kernel at run time (`baseline` / `avx2` / `avx512` / `auto`); it is an environment variable rather than a DuckDB setting because the C++ unit tests need it and they never construct a database.

## 1. Statically Linked from Source (ExternalProject)

### HTSlib 1.22.1
- **Location:** `ext/htslib-1.22.1/` (release tarball, not a submodule)
- **Purpose:** SAM/BAM/CRAM parsing
- **Build:** autotools `configure` + `make`; produces `libhts.a`
- **Configure flags:** `--disable-libcurl --disable-bz2 --disable-lzma --without-libdeflate` (minimal feature set; zlib-only by default, zstd auto-detected via `PKG_CONFIG_PATH`)
- **Gotchas:**
  - **Windows/MinGW:** autotools configure is skipped — `config.h` and `config.mk` are generated manually (rtools42 CI lacks a full MSYS2 shell). `NONCONFIGURE_OBJS=` blanked to prevent `hfile_libcurl.o`.
  - **Emscripten:** `make lib-static` only (shared libs not supported on WASM); cross-compiled with `--host=wasm32-unknown-emscripten`.
  - **macOS cross-compile:** `--host=x86_64-apple-darwin` / `--host=aarch64-apple-darwin` selected from `OSX_BUILD_ARCH`.
  - **Explicit ZLIB linking:** `LIBS=<path-to-libz>` passed to configure because `-lz` doesn't always find vcpkg's zlib on macOS.

### minimap2 (v2.30)
- **Location:** `ext/minimap2/` (git submodule; version captured via `git describe` at configure time → `MINIMAP2_GIT_VERSION`)
- **Purpose:** Long-read / spliced alignment
- **Build:** Makefile; produces `libminimap2.a` with `HAVE_KALLOC`
- **SIMD / arch:**
  - Native ARM64 → `aarch64=1`
  - macOS cross → arch-specific `-arch` flag
  - Emscripten → `sse2only=1` with SIMDe via `-msimd128 -msse -msse2` (avoids SSE4.1 runtime CPU dispatch)
- **Jemalloc integration:** when `DUCKDB_EXTENSION_JEMALLOC_SHOULD_LINK` is set, minimap2 is built with `-DMIINT_USE_JEMALLOC -include src/include/mm_alloc.h` so its allocations route through DuckDB's jemalloc. On loadable-extension builds, `src/mm_alloc_stubs.c` forwards to system malloc because the host DuckDB loaded via Python uses `RTLD_LOCAL` and doesn't expose `duckdb_je_*` symbols.

### WFA2-lib v2.3.6
- **Location:** `ext/WFA2-lib/` (git submodule pinned to v2.3.6)
- **Purpose:** Wavefront Alignment Algorithm for pairwise alignment
- **Build:** Makefile (not CMake — WFA2's primary build is Makefile); `BUILD_TOOLS=0 BUILD_EXAMPLES=0 BUILD_WFA_PARALLEL=0`; produces `libwfa.a` (C core) + `libwfacpp.a` (C++ bindings)
- **Link order:** `wfa2cpp` **before** `wfa2` — C++ depends on C
- **Flag convention:** WFA2 uses `CC_FLAGS` (not `CFLAGS`/`CXXFLAGS`)
- **IUPAC matching:** `WFA2Aligner::align_cigar_semiglobal_iupac` uses the lambda-match callback API (`alignEndsFree(matchFunct, args, ...)`) for degenerate base matching in `extract_linked_amplicon`

### vsearch (v2.30.5-miint fork)
- **Location:** `ext/vsearch/` (fork)
- **Purpose:** search, clustering, UCHIME chimera detection, DUST masking, paired-end merge
- **Build:** autotools; produces `libvsearch.a` (PIC static archive)
- **Gotchas:**
  - Skipped entirely on Windows/Emscripten (autotools not supported)
  - Configure step `touch`es pre-generated autotools files (`aclocal.m4`, `configure`, `Makefile.in`, `config.h.in`, etc.) because fresh `git checkout` gives all files identical timestamps, which can trigger autotools regeneration rules — we don't want that.

### MAFFT PartTree
- **Location:** `ext/mafft/core/`
- **Purpose:** Multiple sequence alignment (PartTree algorithm)
- **Build:** Makefile; produces `libmafft_parttree.a` with `ENABLE_MULTITHREAD=-Denablemultithread`
- **Platform:** POSIX only (uses `mkdtemp`); auto-disabled on Windows

### abPOA (the-miint fork)
- **Location:** `ext/abpoa/` (git submodule at `embed-friction-fixes` on the `the-miint/abPOA` fork; version captured via `git describe` → `ABPOA_GIT_VERSION`)
- **Purpose:** Adaptive banded partial order alignment for MSA and consensus generation (`align_abpoa`, `consensus_abpoa`)
- **Gated by:** `MIINT_ENABLE_ABPOA` (defaults on; auto-disabled on Windows — POSIX APIs). `MIINT_HAS_ABPOA` compile define.
- **Build:** `add_subdirectory(ext/abpoa EXCLUDE_FROM_ALL)`. Unlike minimap2/MAFFT which use `ExternalProject_Add` with Makefiles, abPOA's fork CMakeLists.txt is fully subdirectory-friendly (no global flag pollution, gated install targets, proper multi-ISA SIMD dispatch including WASM SIMD128).
- **Linking:** `$<TARGET_FILE:abpoa>` + `add_dependencies` (not direct target linking) to avoid DuckDB export-set conflicts.
- **SIMD dispatch (x86_64):** CMake compiles `abpoa_align_simd.c` four times (SSE2, SSE4.1, AVX2, AVX512BW) as OBJECT libraries with per-ISA flags, plus a single CPUID dispatcher. The base library uses native `<immintrin.h>` at SSE2 level (no SIMDE on x86_64 dispatch mode).
- **ARM64:** single compile with SIMDE translating AVX2 intrinsics to NEON.
- **WASM:** SIMDE with `-msimd128` (128-bit SIMD, maps to SSE2 path).
- **Symbol namespacing:** the fork prefixes all externally-visible klib symbols (`km_init` → `abpoa_km_init`, `ks_resize` → `abpoa_ks_resize`, etc.) and minimap2-derived symbols (`mm_sketch` → `abpoa_mm_sketch`) to prevent collisions with minimap2 and htslib when statically linked into the same binary.
- **Nested submodule:** SIMDE is a submodule within `ext/abpoa/include/simde/`. Must be initialized (`git submodule update --init --recursive` in `ext/abpoa`).
- **Thread safety:** instance-based — each thread creates its own `abpoa_t`/`abpoa_para_t`. No global mutex needed. `max_threads_hint=0` in per-sample mode allows true parallel execution.
- **Known limitation:** abPOA's internal error handlers (`err_fatal_core` in `src/utils.c`) call `abort()` on OOM. Same class of issue as minimap2's embedded allocator — a very large or malformed input could kill the DuckDB server process.
- **macOS cross-compile:** `CMAKE_OSX_ARCHITECTURES` is set from `OSX_BUILD_ARCH` before `add_subdirectory` and unset after.

### SortMeRNA 4.4.0 (fork)
- **Location:** `ext/sortmerna/` (git submodule at `v4.4.0-miint` on the `the-miint/sortmerna` fork)
- **Purpose:** rRNA read filtering / alignment exposed as `align_sortmerna` (SAM schema) and `align_sortmerna_rrna` (identity/coverage/e-value schema)
- **Gated by:** `MIINT_ENABLE_SORTMERNA` (defaults on for Linux/macOS, off for Emscripten and MinGW — cmph's `cmph_time.h` assumes POSIX `<sys/time.h>`/`gettimeofday` which MinGW doesn't supply). `MIINT_HAS_SORTMERNA` compile define.
- **Build:** `ExternalProject_Add(sortmerna_build)` drives sortmerna's own CMake with `-DCONCURRENTQUEUE_HOME=<ext/concurrentqueue>`, `-DROCKSDB_DIST=<vcpkg_share>`, `-DROCKSDB_USE_STATIC_LIBS=ON`, `-DCMAKE_CXX_FLAGS=-Wno-register`, `-DWITH_TESTS=OFF`.
- **Bundle step:** sortmerna's internal OBJECT libraries (`cmph`, `alp`, `build_version`) are scoped to sortmerna's own CMake configuration. `ExternalProject_Add` isolates that scope from ours, so we cannot `target_link_libraries(... cmph)` against them directly. Additionally, when sortmerna links its STATIC `libsmr_api.a` against those OBJECT libraries, CMake does NOT physically embed the object files into the archive — they remain as loose `.o` files in sortmerna's build tree. `cmake/bundle_sortmerna.cmake` runs after the build, globs those object files, and `ar rs`-appends them into a copy named `libsortmerna_bundle.a`. That is the archive we link. Without this step the final link fails with undefined references to cmph/alp symbols.
- **Forcing a rebuild:** like `rype_build`, `sortmerna_build` is an ExternalProject and caches aggressively. Touching source files in `ext/sortmerna/` does not invalidate it. To force a rebuild: `touch build/release/extension/miint/sortmerna_build-prefix/src/sortmerna_build-stamp/sortmerna_build-configure`.
- **Submodule deps:**
  - `cmph` (minimal perfect hashing) — bundled transitively
  - `alp` (NCBI p-value library) — bundled transitively
  - `concurrentqueue` (lock-free SPMC queue) — header-only; vendored at `ext/concurrentqueue/concurrentqueue.h` rather than submoduled (upstream has no releases and submodule checkout was flaky)
- **RocksDB:** supplied by vcpkg (`find_package(RocksDB CONFIG REQUIRED)`). Pinned via `vcpkg.json` entry `{"name": "rocksdb", "platform": "!emscripten"}`.
- **Threading model:** sortmerna has a process-wide `g_run_mutex` serialising `smr_run_seqs_with_index`. `AlignSortMeRNA*TableFunction::GlobalState::MaxThreads()` returns `1` on the DuckDB side; sortmerna's internal thread pool does the real parallelism (controlled by `num_threads`).
- **Synthetic read IDs:** `SortMeRNAAligner::align` generates sequential CString IDs `"0".."N-1"` for each call and maps sortmerna's per-alignment output back to caller read IDs by position. Avoids upstream's per-batch ID-uniqueness contract and keeps sortmerna's in-memory CString table small.
- **Real-data oracle:** `data/sortmerna/real_oracle.blast.gz` is a gzipped native-CLI BLAST run used by `test/sql/align_sortmerna_realworld.test`. `data/sortmerna/real_oracle.submodule.sha` records the submodule HEAD that produced it; `run_tests.sh` refuses to run the regression test against a stale oracle. Regenerate with `MIINT_SORTMERNA_REAL_DATA=1 bash run_tests.sh`.
- **Known library/CLI divergence:** the streaming library path does not apply the CLI's internal minimum-score filter; our output is a documented strict superset of the CLI's. E-values also differ — library uses per-query Karlin-Altschul; CLI sums/corrects across the DB. Identity / coverage / score / CIGAR / edit_distance are bit-identical.

### unifrac-binaries + scikit-bio-binaries (UniFrac analytics)

Two coupled submodules implementing the UniFrac distance + Faith's PD + PERMANOVA + PCoA stack consumed by `unifrac_distances`, `unifrac_pcoa`, `unifrac_permanova`, and `unifrac_faith_pd`.

- **Locations:**
  - `ext/unifrac-binaries/` — libssu (the UniFrac and Faith PD engine). Version captured via `git describe` → `UNIFRAC_GIT_VERSION`.
  - `ext/scikit-bio-binaries/` — libskbb (PCoA via randomized FSVD and PERMANOVA pseudo-F). Version captured via `git describe` → `SKBB_GIT_VERSION`.
  - libskbb depends on libssu's `support_biom_t` / `support_bptree_t` types but the two compile independently. Both versions are reported in `miint_versions()` (rows `unifrac` and `scikit-bio-binaries`).
- **Build:**
  - `ExternalProject_Add(skbb_build)`: native → `make -C src inmem_static` (after `scripts/fetch_eigen.sh`); WASM → `make -C src wasm`. Produces `libskbb_inmem.a` or `libskbb_wasm.a`.
  - `ExternalProject_Add(unifrac_build)`: depends on `skbb_build`; receives `SKBB_DIR=${SKBB_SRC_DIR}` make-arg. Same target split as skbb; produces `libssu_inmem.a` or `libssu_wasm.a`.
  - **WASM `-fPIC` patch:** upstream's `wasm/emscripten_build.mk` (in both submodules) omits `-fPIC` from `WASM_CXXFLAGS`. DuckDB WASM extensions are side modules (`-sSIDE_MODULE=2`) so every TU must be position-independent — without the patch, `wasm-ld` fails with `R_WASM_MEMORY_ADDR_SLEB cannot be used against symbol ...; recompile with -fPIC`. Idempotent `sed` patches are applied as part of `skbb_build` / `unifrac_build`'s configure step; remove them once upstream lands the fix.
- **Native vs WASM divergence:**
  - **Native (inmem path):** uses cblas + LAPACKE for QR/SVD inside skbb's FSVD; OpenMP enabled via `INMEM_MPFLAG := -fopenmp` (see `ext/unifrac-binaries/src/inmem_build.mk`). `find_package(OpenMP REQUIRED COMPONENTS CXX)` is gated `MIINT_ENABLE_UNIFRAC AND NOT EMSCRIPTEN`; `OpenMP::OpenMP_CXX` is linked into extension/loadable_extension/tests.
  - **WASM:** uses Eigen for QR/SVD (`linalg_backend_eigen.cpp`); no OpenMP, no LAPACKE. Eigen is fetched on demand by `scripts/fetch_eigen.sh` at configure time. WASM-linked archives are added to `DUCKDB_EXTENSION_MIINT_LINKED_LIBS`.
- **Header inclusion contract:** libssu's `api.hpp` ships **without** include guards (only `task_parameters.hpp` and `status_enum.hpp` are properly guarded). All translation units include `src/include/unifrac_libssu.hpp` — a `#pragma once`-wrapped shim — never `api.hpp` directly.
- **`support_biom_t` orientation:** documented in `api.hpp:223` as "CSR-encoded input table" and consumed obs-major by `biom_inmem.cpp:144-156`. Layout: `indptr` length `n_obs+1`; `indices` stores **sample column indices** within each obs row; `data` holds count values. Misreading this as CSC corrupts every downstream call (the live `faith_pd_inmem` test in `test_UnifracSupportBiom.cpp` exists to keep this invariant honest).
- **Process-wide global skbb RNG:**
  - `one_off_matrix_inmem_fp32_v3` consumes libssu's global skbb RNG when `subsample_depth > 0` (via `ssu_set_random_seed(seed) → skbb_set_random_seed(new_seed_skbb)`).
  - `miint::unifrac::UnifracDistanceMatrix::Compute` serializes `ssu_set_random_seed + one_off_matrix_inmem_fp32_v3` behind a process-wide `std::mutex` so concurrent DuckDB connections invoking PCoA/PERMANOVA cannot interleave seed updates.
  - The Faith PD path uses `subsample_table_inmem_seeded` (takes a per-call seed; sidesteps the global RNG entirely) and `faith_pd_inmem` (no RNG); the mutex is not taken on that path. Confirmed by the `BridgeSubsample is reproducible under the same seed` C++ test running without coordination.
- **fp32 reproducibility quirk:** `skbb_pcoa_fsvd_fp32` on identical seeded inputs produces ~1e-7-magnitude bit-level differences across consecutive invocations even with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`. Almost certainly cblas/LAPACKE reduction ordering inside QR/SVD. PCoA SQL seed-reproducibility tests use a `< 1e-5` tolerance. A future fp64 path (not in v1) would give byte equality. Faith PD and PERMANOVA are deterministic (no FSVD in either, and skbb_permanova reconstructs its `RandomGeneratorArray` from the seed every call).
- **`skbb_pcoa_fsvd_fp32` output layout:** the API header comment says "Matrix of size (n_eighs × n_dims)", but the implementation in `principal_coordinate_analysis.cpp:574-578` writes `samples + row*n_eighs` with `row` iterating over **samples** — so the buffer is `(n_samples × n_dims)` sample-major. Index as `samples[sample_idx * n_dims + axis]`. The header comment is misleading; the implementation is canonical.
- **Subsampling drops samples:** when `subsample_depth > 0`, libssu drops samples whose total counts fall below the depth. `mat->n_samples` and `mat->sample_ids` reflect the **post-subsample** view; always trust them over the input feature-table's size. `UnifracDistanceMatrix::sample_ids()` exposes the post-subsample list as a `std::vector<std::string>` for downstream consumers.

### Rust crates via the `miint_rust_glue` umbrella

Two Rust crates — rype and sylph — are statically linked into the extension. They are bundled through a thin umbrella crate at `ext/miint-rust-glue/` that produces a single `libmiint_rust_glue.a`.

**Why an umbrella, not two independent staticlibs?** Each independent Rust staticlib embeds its own copy of the Rust standard library (`rust_eh_personality`, `std::panicking::EMPTY_PANIC`, ...). Linking two of them into the same binary trips duplicate-symbol errors on every supported linker. Previous workarounds (`-Wl,--allow-multiple-definition` on GNU ld, `-Wl,-ld_classic -Wl,-multiply_defined,suppress` on Apple) lived in `extension_config.cmake` until macOS 26 / Xcode 17's ld-prime stopped honoring `-multiply_defined,suppress` even via `-ld_classic`, breaking the duckdb shell, `libduckdb.dylib`, `plan_serializer`, and `unittest` link targets. The umbrella is the canonical fix per the Rust Reference "Linkage" chapter — a single staticlib emitted from one cargo invocation has one shared std and one set of symbols. Same pattern Mozilla has used in Firefox's `rul` super-crate since 2015.

- **Location:** `ext/miint-rust-glue/` (in-tree Cargo crate; `Cargo.toml` lists `rype` and `sylph` as path dependencies)
- **Build:** `ExternalProject_Add(miint_rust_glue_build)` drives `cargo build --release`; produces `libmiint_rust_glue.a`
- **Cargo features:**
  - `with-sylph` — on when `MIINT_ENABLE_SYLPH=ON` (default). Off on Emscripten and Windows/MinGW (sylph leans on POSIX APIs in its sketch indexer).
  - `rype-fastx` — enables rype's needletail-based FASTX reader. Off on Emscripten and Windows/MinGW because needletail's `cdylib` crate type triggers linker errors (standalone-linking failures on WASM; unresolved `___chkstk_ms` on MSVC).
  - `arrow-ffi` — always on for both crates; the extension consumes both via the Arrow C Data Interface for zero-copy FFI.
- **CMake exposure:** the umbrella is wired up as two `IMPORTED STATIC` targets — `rype` and `sylph` — both pointing at the same `libmiint_rust_glue.a`. Existing call sites like `target_link_libraries(... rype)` and `target_link_libraries(... sylph)` work unchanged; the linker dedupes the duplicate archive on the link line.
- **Cross-compilation gotchas:**
  - Cargo `--target` flag selected per platform; `rustup target add` called at configure time.
  - **Emscripten:** `--target wasm32-unknown-emscripten` + `RUSTFLAGS=-C relocation-model=pic` (default wasm32 is static; DuckDB WASM side modules need PIC).
  - **Windows/MinGW:** `--target x86_64-pc-windows-gnu` (NOT the MSVC default) so object files are compatible with MinGW `ld.exe`.
  - **macOS cross:** explicit `x86_64-apple-darwin` / `aarch64-apple-darwin` targets.
- **Forcing a rebuild:** ExternalProject caches `miint_rust_glue_build` aggressively. Touching source files in `ext/rype/` or `ext/sylph/` does not invalidate it. To force a rebuild: `touch build/release/extension/miint/miint_rust_glue_build-prefix/src/miint_rust_glue_build-stamp/miint_rust_glue_build-configure`.

#### rype subsystem

- **Location:** `ext/rype/` (git submodule on `the-miint/rype`; version captured via `git describe` → `RYPE_GIT_VERSION`)
- **Purpose:** Rust implementations of classify, extract, and log-ratio operations, exposed via Arrow C Data Interface
- **Compile definition:** `RYPE_ARROW`
- **Reported as:** `rype` row in `miint_versions()`

#### sylph subsystem

- **Location:** `ext/sylph/` (git submodule on `the-miint/sylph`, branch `v0.9.0-miint`; version captured via `git describe` → `SYLPH_GIT_VERSION`)
- **Purpose:** FracMinHash sketch-based relative-abundance profiling of microbial communities, exposed as the `sylph_profile` table function
- **Gated by:** `MIINT_ENABLE_SYLPH` (auto-off on Emscripten + Windows/MinGW per the POSIX-API requirement above). `MIINT_HAS_SYLPH` compile define when on; `SYLPH_GIT_VERSION` carries the configure-time `git describe` string.
- **Reported as:** `sylph` row in `miint_versions()` (only emitted when `MIINT_HAS_SYLPH` is defined)
- **C++ wrappers:**
  - `src/include/SylphDatabase.hpp` / `src/SylphDatabase.cpp` — `SylphDatabaseHandle`, a RAII wrapper around the C FFI `SylphDatabase*` from `ext/sylph/sylph.h`. Non-copyable, non-movable; owned by the `sylph_profile` GlobalState.
  - `src/include/sylph_profile.hpp` / `src/sylph_profile.cpp` — `SylphProfileTableFunction`, the DuckDB table-function binding. The `.syldb` is loaded once into GlobalState and shared read-only across worker threads; sylph treats the loaded database as immutable so no read-side mutex is required.

## 2. Header-Only

### kseq++
- **Location:** `ext/kseq++/`
- **Purpose:** Modern C++ FASTA/FASTQ parser
- **Integration:** header-only, included directly via `include_directories(ext)`
- **Requires C++20** (drives the project-wide `CMAKE_CXX_STANDARD 20` setting)

### LBFGS++
- **Location:** `ext/LBFGSpp/` (tracked file copies of upstream v0.4.0; per-file checksums and the omitted-files list in `ext/LBFGSpp/PROVENANCE.md`)
- **Purpose:** L-BFGS quasi-Newton optimizer; drives the MAP fit in `mmvec` (`src/mmvec.cpp`)
- **Integration:** header-only, no CMake change needed — upstream's two-level `include/LBFGS.h` + `include/LBFGSpp/*.h` is vendored **flat** into `ext/LBFGSpp/`, so `#include "LBFGSpp/LBFGS.h"` resolves through the existing `include_directories(ext)` and the headers' own relative includes resolve unpatched. Same flattening as kseq++.
- **Reported as:** `LBFGS++` row in `miint_versions()` (a literal — upstream ships no version macro)
- **Depends on Eigen** (`<Eigen/Core>`, `<Eigen/LU>` via `BFGSMat.h`), which the build already requires. Compiles clean under the globally-set `EIGEN_MPL2_ONLY`; nothing in the transitive include set trips that guard.
- **Only the unbounded solver is vendored.** `mmvec`'s parameters are unconstrained, so `LBFGSB.h` / `Cauchy.h` / `SubspaceMin.h` are omitted.
- **Gotchas for callers:** `LBFGSParam::ftol` and `wolfe` are line-search Armijo/Wolfe constants, *not* convergence tolerances (SciPy's convergence `ftol` maps to `past`/`delta`); and `minimize()` **throws** on line-search failure rather than returning a status, so callers that must survive non-convergence keep their own best-so-far snapshot.

## 3. System / vcpkg Libraries

| Library | How | Required | Notes |
|---|---|---|---|
| zlib | `find_package(ZLIB REQUIRED)` | yes | vcpkg or system |
| Threads | `find_package(Threads REQUIRED)` | yes | pthreads / Win32 |
| expat | `find_package(expat CONFIG REQUIRED)` | yes | XML parsing for mzML/mzXML |
| Eigen | `find_package(Eigen3 CONFIG QUIET)`; falls back to fetching the pinned, SHA-verified header-only 3.4.0 tarball | yes | Header-only, so nothing is linked. Backs `NewickTree.cpp`'s Mk rate matrices and the vendored LBFGS++ that drives `mmvec`. **Two globally-set defines, both load-bearing** — see below. |
| zstd | `find_package(zstd QUIET)` | optional | multiple target names tried (`zstd::libzstd_static`, `zstd::libzstd_shared`, `zstd::libzstd`, `zstd`); gates `HAVE_LIBZSTD` |
| HDF5 | `find_package(HDF5 REQUIRED COMPONENTS CXX)` | conditional | required when `MIINT_ENABLE_HDF5=ON`; links `hdf5-static`, `hdf5_hl-static`, `hdf5_cpp-static`, `hdf5_hl_cpp-static`. For `CLANG_TIDY` builds HDF5 is optional (tidy still runs without BIOM coverage). |
| RocksDB | `find_package(RocksDB CONFIG REQUIRED)` | conditional | required when `MIINT_ENABLE_SORTMERNA=ON`; vcpkg port gated via `"platform": "!emscripten"`. Linked against `RocksDB::rocksdb` (static archive). |
| libcurl | `find_package(CURL CONFIG QUIET)`, falling back to the system `FindCURL` module | conditional | required when `MIINT_ENABLE_CURL=ON`; powers the streaming HTTPS-PUT / FTP / FTPS transport in `ena_upload_reads`. Reported as the `libcurl` row in `miint_versions()`. |
| Catch2 | `FetchContent` @ v3.4.0 | yes (test) | used only for the `tests` C++ executable |

`VCPKG_TOOLCHAIN_PATH` must be set before CMake configure.

### The two Eigen defines

Both are `add_compile_definitions` in `CMakeLists.txt` — global, not per-source — and both are set for a reason that is easy to undo by accident.

- **`EIGEN_MPL2_ONLY`** confines the build to Eigen's pure-MPL2 subset, keeping miint's Modified-BSD licence clean. Nothing in the transitive include set of either Eigen user trips the guard.
- **`EIGEN_DONT_PARALLELIZE`** exists because the extension links OpenMP (libssu is built `-fopenmp`), so `_OPENMP` is defined for miint's own sources too and Eigen would otherwise enable its parallel GEMM path. `mmvec`'s objective contains matrix products reducing over a long dimension (d1 = 2720 at cystic-fibrosis scale), and thread-count-dependent blocking would make a **seeded fit non-reproducible**. Measured at CF scale, seed 0, 500 iterations: with the define removed the same fit returns a different model as soon as `OMP_NUM_THREADS > 1` (final loss 1709324046711880.2 against 1709329920514833.8). It costs 3–5% of the fit's wall clock. The fit *is* threaded, but by `miint::WorkerPool` (`src/include/miint_parallel.hpp`) around the products rather than inside them, which is exactly why its values do not move — so this define is more load-bearing now, not less.

Global rather than per-source because Eigen's macros configure templates instantiated in whichever TU includes the headers — a per-file define would leave two differently-configured instantiations of the same symbols for the linker to choose between.

### Per-instruction-set compilation of the mmvec kernel

Eigen picks its packet width at **compile** time from `__AVX__` / `__AVX512F__` and does not runtime-dispatch. `__attribute__((target))` does not help — under it those macros are undefined at preprocessing time and Eigen still emits SSE2. Reaching a wider instruction set therefore means compiling a source *again* with different flags and choosing between the results at load time, which is what `miint_add_isa_variants()` in `CMakeLists.txt` does. Its only caller today is `src/mmvec_kernel.cpp`; the mechanism (`src/include/miint_isa.hpp`, the `MIINT_SIMD` override, the `MIINT_HAS_AVX2` / `MIINT_HAS_AVX512` feature macros) is deliberately subsystem-agnostic so a future dense kernel need not re-derive it.

The precedent is htslib's htscodecs, which does the same thing and already ships through our distribution pipeline: the published `linux_amd64` artifact carries 78,583 AVX-512 and 25,511 AVX2 instructions (`rans_*_avx2` / `_avx512`, abPOA's aligner, OpenSSL's AES). Runtime dispatch needs **no CI or packaging change** — which matters because `DuckDBPlatform()` has no CPU-feature component and so cannot express a per-ISA download.

Three things here are load-bearing and easy to undo by accident:

- **The baseline variant takes no extra flags.** Every carved expected value in the repo is carved against it, and `run_tests.sh` pins `MIINT_SIMD=baseline` so CI asserts a fixed contract rather than whatever the runner's CPU offers.
- **The read path does not dispatch.** `ComputeLogits` calls `baseline::FillNonRefLogits` explicitly, so `mmvec_ranks` / `mmvec_predict` / `mmvec_score` stay bit-identical on every CPU for a given model. It costs nothing — that path runs once per call, not once per iteration.
- **`Workspace`'s Eigen-mapped buffers are SIMD-aligned** (`src/include/miint_aligned_vector.hpp`), and that is a *correctness* requirement, not tuning. Eigen peels leading scalars off an under-aligned buffer onto a different arithmetic path and decides how many from the runtime **address**, so without it the same seeded fit returns different models in one process as the allocator moves the scratch around. It reaches the matrix products, not just the exponential. The alternative — `EIGEN_MAX_ALIGN_BYTES=0`, which makes Eigen peel nothing — was measured at +75% (AVX2) and +113% (AVX-512), putting AVX-512 behind the SSE2 baseline; alignment costs 0.0025% of instructions retired.

A related trap in the fallback path: the pinned Eigen tarball's own `CMakeLists.txt` sets `-msse2` / `-mavx` / `-mavx512f` **globally** (lines 207–237). It is inert only because `FetchContent_Declare` points `SOURCE_SUBDIR` at a nonexistent directory so Eigen's CMake never runs. Removing that would inject arch flags into every target and silently move the baseline variant.

## 4. Runtime Binaries (not compiled in)

These are invoked via `fork`/`exec` at runtime; the extension links no code for them. Tests that exercise them are guarded in `run_tests.sh` by `command -v` detection.

### gpl-boundary (process-isolation host for GPL tools)
- **Gated by:** `MIINT_ENABLE_GPL_BOUNDARY` (controls whether `src/gpl_boundary/{process,session,shm,arrow_ipc}.cpp`, `src/phylogeny_fasttree.cpp`, `src/align_bowtie2.cpp`, `src/align_bowtie2_sharded.cpp`, `src/align_bowtie2_daemon_common.cpp`, and the vendored `third_party/nanoarrow_ipc/` object library are compiled). Auto-off on Emscripten and Windows.
- **Why a separate process?** FastTree and bowtie2 are GPL-licensed; miint is BSD. Statically linking either into the extension would cross the license boundary. Instead, they're statically linked into [`gpl-boundary`](https://github.com/the-miint/GPL-boundary), an independent GPL-licensed binary that miint launches as a child process and communicates with over a JSON-line control channel and POSIX shared memory.
- **Tools advertised by the daemon (v0.2.0):** `fasttree` (schema_version 2), `bowtie2-align` (schema_version 2), `bowtie2-build` (schema_version 1), `sortmerna`, `prodigal`. miint's per-tool wrappers fail loud at session boot if a required tool is missing or its schema_version doesn't match.
- **Wrapper:** `src/gpl_boundary/`
  - `process.{cpp,hpp}` — `FindGplBoundary` (PATH lookup), `ChildProcess` (fork/exec/wait with SIGTERM-then-SIGKILL graceful shutdown), `LineReader`/`WriteLine` (newline-framed JSON I/O over pipes).
  - `session.{cpp,hpp}` — `Session::Initialize` (handshake on `protocol_version=2`), `Session::Submit` (per-batch round trip with mandatory `shm_input_size` field), `Session::Shutdown`.
  - `shm.{cpp,hpp}` — `InputShmRegion` (created by miint, written by miint, **unlinked by miint**), `OutputShmRegion` (created by gpl-boundary, read by miint, **unlinked by miint** per gpl-boundary's README convention). Size is authoritative on both sides — neither side calls `fstat`.
  - `arrow_ipc.{cpp,hpp}` — `EncodeIpcStream` (DuckDB DataChunk → Arrow IPC stream bytes via vendored `nanoarrow_ipc`), `IpcStreamDecoder` (response bytes → `ArrowArrayWrapper`s).
- **Wire protocol invariants** (gpl-boundary v0.2.0+):
  - `protocol_version: 3` — the Init handshake hard-fails on mismatch so daemon-side bumps surface immediately at session boot rather than silently producing wrong data. The Init reply also returns the full tools registry (name + schema_version) so each wrapper can validate its expected schema upfront.
  - Every batch request includes `shm_input_size` (the exact byte count miint passed to `ftruncate`). Daemon does not `fstat` shm fds; explicit size is authoritative.
  - Output schema for `fasttree`: 8 columns (`node_index Int64 not-null`, `parent_index Int64 nullable`, `edge_id Int64 nullable`, `branch_length Float64 nullable`, `support Float64 nullable`, `n_children Int32 not-null`, `is_tip Boolean not-null`, `name Utf8 nullable`). Wire `n_children` is Int32; miint widens to Int64 and reorders to `is_tip, name, n_children` for SQL ergonomics. Schema-drift detection in `InitGlobal` checks each column name at every position — silent reorders fail loudly.
  - Output schema for `bowtie2-align`: 21 columns matching the `read_alignments` shape (`read_id`, `flags`, `reference`, `position`, `stop_position`, `mapq`, `cigar`, `mate_reference`, `mate_position`, `template_length`, then the optional SAM tags). `tag_sa` is always-NULL on this path (parity stub). Validated against `bt2_daemon::kOutputColumnNames` in `src/include/align_bowtie2_daemon_common.hpp`.
- **Lifecycle:**
  - One daemon spawned per table-function call (e.g. `phylogeny_fasttree(...)`, `align_bowtie2(...)`, `align_bowtie2_sharded(...)`); no cross-query caching yet — connection-scoped reuse via `ClientContext::registered_state` is a planned follow-up.
  - For `align_bowtie2`: `Bind` first submits the subjects to the daemon's `bowtie2-build` tool, writing index files into a per-call `mkdtemp("miint-bt2-XXXXXX")` directory under `$TMPDIR`. `Execute` streams the query batches to `bowtie2-align`. The `GlobalState` destructor calls `Session::Shutdown`, unlinks the `.bt2` files, and removes the temp directory.
  - Daemon shutdown via `~GlobalState` → `Session::Shutdown` → `~ChildProcess`: SIGTERM, then 30 × 10ms grace, then SIGKILL.
  - SIGPIPE is **blocked per-thread** via `pthread_sigmask` + drained with `sigtimedwait` — never `sigaction(SIGPIPE, SIG_IGN)` (would leak to other threads in the host process).
- **Vendored nanoarrow_ipc** (`third_party/nanoarrow_ipc/`):
  - DuckDB has no Arrow IPC byte-serialization at the C++ level (only the C Data Interface). gpl-boundary's wire format is FlatBuffers-framed Arrow IPC stream bytes, so we ship the official `nanoarrow_ipc` (~3k LOC, Apache 2.0) plus the `nanoarrow` C runtime and `flatcc` runtime it depends on.
  - **Pinned tag:** `apache-arrow-nanoarrow-0.8.0` (commit `a579fbf5...`).
  - **Symbol namespace:** `NANOARROW_NAMESPACE=miint` set in `third_party/nanoarrow_ipc/include/nanoarrow/nanoarrow_config.h`. Mangles every `Arrow*` and `ArrowIpc*` symbol to `miint_Arrow*` so we can't collide with DuckDB's bundled `duckdb_nanoarrow` (a C++ namespace).
  - **CMake integration:** declared as an OBJECT library `miint_nanoarrow_ipc`. The main `CMakeLists.txt` wires it via `target_sources(... PRIVATE $<TARGET_OBJECTS:miint_nanoarrow_ipc>)` rather than `target_link_libraries` — the latter triggers DuckDB's export-set machinery and demands the static lib be exported, which it has no business being.
- **Parity oracle:** `data/fasttree/{tiny,moderate}.golden.nwk` are bioconda FastTree binary outputs for committed input fixtures (`*.fa`); freshness is pinned via `*.golden.fixture.sha` (input-content hash) and `*.golden.fasttree.version` (binary version). `run_tests.sh` exports `MIINT_FASTTREE_{TINY,MODERATE}_PARITY_OK=1` only when the matching SHA matches. `test/sql/phylogeny_fasttree_parity.test` then verifies miint's daemon path produces topologically-identical trees (RF=0 via canonical bipartition multisets) with branch lengths matching within hybrid tolerance `max(1e-9, 1e-6 × max(|m|,|g|))`. Regenerate with `MIINT_FASTTREE_REGENERATE=1 bash run_tests.sh`.

### IBM Aspera ascp
- **Gated by:** `MIINT_ASPERA_SUPPORTED` (build-time; `0` on Windows/WASM)
- **Wrapper:** `src/aspera_utils.cpp`, `src/aspera_stream.cpp`
- Detected at runtime via `which ascp`; not required at build time
- SSH key auto-discovered at `~/.aspera/connect/etc/`, `$CONDA_PREFIX/etc/`, or downloaded from GitHub on demand
- Uses `stdio://` and `stdio-tar://` protocols to stream directly to stdout (confirmed working with ENA FASP servers)
- `read_ena_fastx` with `download_method='auto'` falls back to HTTP transparently when `ascp` isn't available

## 5. Vendored Source (compiled inline)

Small third-party sources dropped into the tree and compiled directly into the extension (added to both `EXTENSION_SOURCES` and `TEST_SOURCES`), rather than built as a separate `ExternalProject` or included header-only.

### microtar 0.1.0
- **Location:** `third_party/microtar/` (`microtar.c` + `microtar.h`, vendored verbatim from https://github.com/rxi/microtar; MIT, see `third_party/microtar/LICENSE` and `README.md`)
- **Purpose:** read `.tar` members. Used by `src/taxdump_archive.cpp` to extract `nodes.dmp`/`names.dmp`/`merged.dmp`/`delnodes.dmp` from NCBI's `taxdump.tar.gz` for `read_ncbi_taxdump`.
- **Build:** compiled as C directly into the extension; include path added via `include_directories(... third_party/microtar ...)`. No feature flag — always built (pure ANSI C, no platform deps).
- **Layering:** microtar owns only the `.tar` layer; the `.gz` layer is inflated separately by the existing zlib helper before microtar reads the resulting tar bytes from an in-memory buffer via its stream callbacks. Only the read side of the API is exercised.

## Platform Link-Deps Summary

Extra libraries pulled in per platform (beyond the main target_link_libraries set):

- **Apple:** `-framework CoreFoundation`
- **Windows/MinGW:** `regex` (HTSlib POSIX regex), `ws2_32 ntdll userenv bcrypt advapi32` (Rust x86_64-pc-windows-gnu stdlib deps — from `rustc_target` spec and `std #[link]` attributes)
- **Unix:** `dl m`

## WASM Loadable Extension Special Case

DuckDB's `build_loadable_extension` on Emscripten runs a POST_BUILD `emcc` command to create the `.wasm` side module. It reads extra libs to link from `DUCKDB_EXTENSION_MIINT_LINKED_LIBS`, **not** from `target_link_libraries` (which only records metadata for `STATIC IMPORTED` targets on Emscripten). Without this explicit list, symbols like `gzclose` become unresolved imports and fail at load time with `"bad export type: undefined"`. The list is assembled in `CMakeLists.txt` around line 707–727.

