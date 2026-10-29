# Build system (CMake)

How miint's CMake is laid out, and what to edit when you add a source, a test or a dependency. Why each embedded library is built the way it is lives in [`embedded-tools.md`](embedded-tools.md).

## Layout

| File | Holds |
|---|---|
| `CMakeLists.txt` | The `MIINT_ENABLE_*` options, the **platform matrix** (the only place a platform turns a feature off), the includes below in order, the always-compiled sources, and the three targets. |
| `cmake/miint_helpers.cmake` | `miint_git_version()`, `miint_configure_target()`, `miint_check_unique()`, and the cross-compile flags (`MIINT_OSX_ARCH_CFLAGS`, `MIINT_AUTOTOOLS_HOST`, `EMSCRIPTEN_THREAD_CFLAGS`). |
| `cmake/system_deps.cmake` | Libraries found rather than built: ZLIB, Threads, expat, Eigen (and its two defines), libdeflate, zstd. |
| `cmake/openmp.cmake` | The `miint_openmp` target (UniFrac needs it; krepp uses it if present). |
| `cmake/isa_dispatch.cmake` | Per-instruction-set builds of a kernel (`miint_add_isa_variants`). |
| `cmake/deps/<name>.cmake` | One embedded dependency each: its build, and everything miint compiles or links because of it. The file name is the option's suffix: `MIINT_ENABLE_KREPP` → `krepp.cmake`. `htslib`, `minimap2`, `wfa2` and `rust_glue` are always built. Exceptions: sylph and st3 are cargo features of `rust_glue.cmake`, `skb_arrow.cmake` is included from `gpl_boundary.cmake` (it shares that option and transport), and libdeflate is in `system_deps.cmake`. |
| `cmake/scripts/` | Build-time scripts run with `cmake -P` (they see only their `-D` arguments), not `include()`d. |

## The registration lists

A deps file builds its library, then appends to these lists. `include()` shares the caller's scope, so plain `list(APPEND ...)` is all it takes.

| List | Holds | Applied to |
|---|---|---|
| `MIINT_EXTENSION_SOURCES` | `src/*.cpp` compiled only when the feature is on | static + loadable extension |
| `MIINT_TEST_SOURCES` | DuckDB-free `src/*.cpp` and `test/cpp/test_*.cpp` | `tests` |
| `MIINT_DEFINITIONS` | `MIINT_HAS_<NAME>`; `<NAME>_GIT_VERSION` via `miint_git_version()` | all three |
| `MIINT_INCLUDE_DIRS` | header directories | every target in the directory |
| `MIINT_LINK_LIBRARIES` | link items, **in link order** | all three |
| `MIINT_OBJECT_SOURCES` | `$<TARGET_OBJECTS:...>` | all three |
| `MIINT_TARGET_DEPENDENCIES` | targets that must build first (e.g. ones linked by path) | all three |
| `MIINT_WASM_LIBS` | archive paths for the wasm side-module link | loadable, under Emscripten |

`miint_configure_target()` applies them, so the static extension, the loadable extension and `tests` cannot drift apart. What only one target needs sits next to its call in `CMakeLists.txt`.

## Recipes

**Add a source that is always compiled.** Add it to `MIINT_EXTENSION_SOURCES` in `CMakeLists.txt`, keeping the list sorted. If it is DuckDB-free and has a C++ test, add it and the test to `MIINT_TEST_SOURCES` too.

**Add a source that belongs to an optional feature.** Append it in that feature's deps file, never in the core list: there it would compile with the feature off. `miint_check_unique()` stops configure if a source is in both.

**Add a C++ unit test.** Append `test/cpp/test_<Name>.cpp`, plus the `src/` files it exercises, to `MIINT_TEST_SOURCES`: in `CMakeLists.txt`, or in the feature's deps file if it needs that feature. A test that needs a running DuckDB goes in the `tests_duckdb` executable at the end of `CMakeLists.txt` instead.

**Add an embedded dependency.**
1. Write `cmake/deps/<name>.cmake`: build or import the library, then append to the lists above. Call `miint_git_version(<NAME> ext/<dir>)` if it is a submodule.
2. Add `option(MIINT_ENABLE_<NAME> ...)` and, if some platform cannot build it, a line in the platform matrix.
3. Add `<NAME>` to the feature `foreach` in `CMakeLists.txt`, at its link-order position. A dependency that is always built goes in `foreach(dep htslib minimap2 wfa2 rust_glue)` instead, with no option or matrix line.
4. Report the version and document the library: [`embedded-tools.md`](embedded-tools.md#reporting-a-new-dependencys-version).

## Rules

- **Include order is link order.** The linker resolves an archive's undefined symbols only from libraries after it, so a library goes after every library that uses it. zlib, expat, the compression codecs, Threads and the platform libraries are appended last, in `CMakeLists.txt`.
- **Link by path what is outside DuckDB's export set.** The static extension is `install(EXPORT)`ed with DuckDB, which rejects link dependencies on targets outside that set. Link such a library as `$<TARGET_FILE:x>` and add `x` to `MIINT_TARGET_DEPENDENCIES` (see `abpoa.cmake`, `krepp.cmake`); feed OBJECT libraries through `MIINT_OBJECT_SOURCES`.
- **Platform gates live only in the matrix.** A deps file may still turn its own feature off when a package is missing (`set(MIINT_ENABLE_<NAME> OFF)` then `return()`), as the tidy-check builds need.
- **Headers whose names collide get a narrower scope.** `MIINT_INCLUDE_DIRS` reaches every target in the directory, and its order decides which of two same-named headers an include finds when the header is not beside the including file. When a library ships headers named like another's, scope them to the files that need them with `set_source_files_properties` (see `krepp.cmake`).
- **`tests` sources must be DuckDB-free.** `tests` has no libduckdb and is built without `MIINT_STATIC_BUILD`, so a source in both lists is compiled twice with different `#ifdef`s. That is why the lists cannot share objects.
- **Format `cmake/` yourself.** `make format-check` covers the root `CMakeLists.txt` but nothing under `cmake/`; run `cmake-format -i` on what you change there, with cmake-format 0.6.13 (the version CI pins). CI's `cmake-format` job checks all of `cmake/`.
