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
| `MIINT_ENABLE_BOWTIE2` | ON | Emscripten, Windows (requires POSIX fork/exec) |
| `MIINT_ENABLE_MAFFT` | ON | Windows (uses `mkdtemp` and other POSIX APIs; segfaults on MinGW) |
| `MIINT_ENABLE_VSEARCH` | ON | Emscripten, Windows (autotools build not supported) |

Corresponding preprocessor macros: `MIINT_HAS_HDF5`, `MIINT_HAS_BOWTIE2`, `MIINT_HAS_MAFFT`, `MIINT_HAS_VSEARCH`. Also `MIINT_ASPERA_SUPPORTED=0` on Windows/WASM (POSIX-only runtime).

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
| Catch2 | `FetchContent` @ v3.4.0 | yes (test) | used only for the `tests` C++ executable |

`VCPKG_TOOLCHAIN_PATH` must be set before CMake configure.

## 4. Runtime Binaries (not compiled in)

These are invoked via `fork`/`exec` at runtime; the extension links no code for them. Tests that exercise them are guarded in `run_tests.sh` by `command -v` detection.

### bowtie2 / bowtie2-build
- **Gated by:** `MIINT_ENABLE_BOWTIE2` (controls whether the wrapper sources `Bowtie2Aligner.cpp`, `align_bowtie2.cpp`, `align_bowtie2_sharded.cpp` are compiled)
- **Wrapper:** `src/Bowtie2Aligner.cpp` — locates binaries on `PATH` via `fork`/`exec` of `which`; manages a per-process temp directory; calls `bowtie2-build` for indexing and `bowtie2` for alignment.
- **Platform:** POSIX only (auto-disabled on Windows/WASM)

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

