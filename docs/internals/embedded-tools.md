# Embedded Tools and Dependencies

How duckdb-miint embeds external libraries and tools. `CMakeLists.txt` is the source of truth; this file explains the *why* and the non-obvious gotchas.

Four embedding categories:
1. **Static libraries built from source** via `ExternalProject_Add`, linked into the extension
2. **Header-only libraries** included directly
3. **System/vcpkg libraries** via `find_package`
4. **Runtime binaries** invoked via `fork`/`exec` (not compiled in)

## Feature Flags (CMake Options)

| Flag | Default | Auto-disabled when |
|---|---|---|
| `MIINT_ENABLE_HDF5` | ON | Emscripten (C++ static class members become unresolvable GOT.mem imports in WASM) |
| `MIINT_ENABLE_MAFFT` | ON | Windows (uses `mkdtemp` and other POSIX APIs; segfaults on MinGW) |
| `MIINT_ENABLE_VSEARCH` | ON | Emscripten, Windows (autotools build not supported) |
| `MIINT_ENABLE_SORTMERNA` | ON | Emscripten (RocksDB vcpkg port not built for wasm32), Windows/MinGW (cmph assumes POSIX `<sys/time.h>`; MSVC-on-Windows would work if anyone wires it up) |
| `MIINT_ENABLE_GPL_BOUNDARY` | ON | Emscripten, Windows (subsystem uses POSIX shm + fork/exec) |

Corresponding preprocessor macros: `MIINT_HAS_HDF5`, `MIINT_HAS_MAFFT`, `MIINT_HAS_VSEARCH`, `MIINT_HAS_SORTMERNA`, `MIINT_HAS_GPL_BOUNDARY`. Also `MIINT_ASPERA_SUPPORTED=0` on Windows/WASM (POSIX-only runtime).

Run-time / conditional: `MIINT_USE_JEMALLOC` is set when DuckDB's jemalloc is linked (not on musl/macOS/Windows).

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

### WFA2-lib v2.3.5
- **Location:** `ext/WFA2-lib/` (git submodule pinned to v2.3.5)
- **Purpose:** Wavefront Alignment Algorithm for pairwise alignment
- **Build:** Makefile (not CMake — WFA2's primary build is Makefile); `BUILD_TOOLS=0 BUILD_EXAMPLES=0 BUILD_WFA_PARALLEL=0`; produces `libwfa.a` (C core) + `libwfacpp.a` (C++ bindings)
- **Link order:** `wfa2cpp` **before** `wfa2` — C++ depends on C
- **Flag convention:** WFA2 uses `CC_FLAGS` (not `CFLAGS`/`CXXFLAGS`)
- **Known bug:** BiWFA alignment-scope score returns `INT32_MIN` for short sequences (see `ext/WFA2-lib/bialign_score_bug.c`)

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

### rype (Rust, Arrow FFI)
- **Location:** `ext/rype/` (git submodule; version captured via `git describe` → `RYPE_GIT_VERSION`)
- **Purpose:** Rust implementations (e.g., classify, extract, log-ratio) exposed via Arrow C Data Interface for zero-copy FFI
- **Build:** `cargo rustc --release --features arrow-ffi --lib --crate-type=staticlib`; produces `librype.a`
- **Compile definition:** `RYPE_ARROW`
- **Cross-compilation gotchas:**
  - Cargo `--target` flag selected per platform; `rustup target add` called at configure time
  - **Emscripten:** `--target wasm32-unknown-emscripten` + `RUSTFLAGS=-C relocation-model=pic` (default wasm32 is static; side modules need PIC)
  - **Windows/MinGW:** `--target x86_64-pc-windows-gnu` (NOT the MSVC default) so object files are compatible with MinGW `ld.exe`
  - **macOS cross:** explicit `x86_64-apple-darwin` / `aarch64-apple-darwin` targets
- **Feature gating:** on Emscripten and Windows, `--no-default-features --features arrow-ffi` disables the `fastx` (needletail) feature whose `cdylib` crate type triggers linker errors (standalone linking on WASM; unresolved `___chkstk_ms` on MSVC).
- **Forcing a rebuild:** ExternalProject caches `rype_build`. Touching source files doesn't invalidate — you must touch the configure stamp: `touch build/release/extension/miint/rype_build-prefix/src/rype_build-stamp/rype_build-configure`

## 2. Header-Only

### kseq++
- **Location:** `ext/kseq++/`
- **Purpose:** Modern C++ FASTA/FASTQ parser
- **Integration:** header-only, included directly via `include_directories(ext)`
- **Requires C++20** (drives the project-wide `CMAKE_CXX_STANDARD 20` setting)

## 3. System / vcpkg Libraries

| Library | How | Required | Notes |
|---|---|---|---|
| zlib | `find_package(ZLIB REQUIRED)` | yes | vcpkg or system |
| Threads | `find_package(Threads REQUIRED)` | yes | pthreads / Win32 |
| expat | `find_package(expat CONFIG REQUIRED)` | yes | XML parsing for mzML/mzXML |
| zstd | `find_package(zstd QUIET)` | optional | multiple target names tried (`zstd::libzstd_static`, `zstd::libzstd_shared`, `zstd::libzstd`, `zstd`); gates `HAVE_LIBZSTD` |
| HDF5 | `find_package(HDF5 REQUIRED COMPONENTS CXX)` | conditional | required when `MIINT_ENABLE_HDF5=ON`; links `hdf5-static`, `hdf5_hl-static`, `hdf5_cpp-static`, `hdf5_hl_cpp-static`. For `CLANG_TIDY` builds HDF5 is optional (tidy still runs without BIOM coverage). |
| RocksDB | `find_package(RocksDB CONFIG REQUIRED)` | conditional | required when `MIINT_ENABLE_SORTMERNA=ON`; vcpkg port gated via `"platform": "!emscripten"`. Linked against `RocksDB::rocksdb` (static archive). |
| Catch2 | `FetchContent` @ v3.4.0 | yes (test) | used only for the `tests` C++ executable |

`VCPKG_TOOLCHAIN_PATH` must be set before CMake configure.

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

## Platform Link-Deps Summary

Extra libraries pulled in per platform (beyond the main target_link_libraries set):

- **Apple:** `-framework CoreFoundation`
- **Windows/MinGW:** `regex` (HTSlib POSIX regex), `ws2_32 ntdll userenv bcrypt advapi32` (Rust x86_64-pc-windows-gnu stdlib deps — from `rustc_target` spec and `std #[link]` attributes)
- **Unix:** `dl m`

## WASM Loadable Extension Special Case

DuckDB's `build_loadable_extension` on Emscripten runs a POST_BUILD `emcc` command to create the `.wasm` side module. It reads extra libs to link from `DUCKDB_EXTENSION_MIINT_LINKED_LIBS`, **not** from `target_link_libraries` (which only records metadata for `STATIC IMPORTED` targets on Emscripten). Without this explicit list, symbols like `gzclose` become unresolved imports and fail at load time with `"bad export type: undefined"`. The list is assembled in `CMakeLists.txt` around line 707–727.

