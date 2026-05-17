# This file is included by DuckDB's build system. It specifies which extension to load

# MinGW's ld treats duplicate COMDAT symbols (C++ constexpr/inline from DuckDB headers)
# as errors when the extension's static lib and DuckDB's static lib are linked into the
# same binary. This is set here (not in the extension's CMakeLists.txt) because this file
# runs in DuckDB's CMake scope before libduckdb.dll/duckdb.exe targets are configured.
#
# Additionally: rype and sylph are both Rust staticlibs that statically embed the Rust
# standard library, so on platforms where ld pulls every archive member (or processes
# them with --whole-archive equivalent), each archive's copy of rust_eh_personality,
# ARGV_INIT_ARRAY, EMPTY_PANIC, ... causes "multiple definition" errors. This affects
# GNU ld (Linux) and MinGW ld; Apple's ld-prime handles duplicates from .a archives via
# Mach-O standard archive semantics (first-wins) without needing a flag, and in fact
# rejects --allow-multiple-definition as an unknown option. Emscripten's wasm-ld is
# similarly out of scope. Restrict the flag to GNU-ld targets only.
if(MINGW OR (NOT WIN32 AND NOT APPLE AND NOT EMSCRIPTEN))
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--allow-multiple-definition")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--allow-multiple-definition")
endif()

# Apple equivalent of --allow-multiple-definition: rype + sylph each carry
# their own copy of the Rust std (rust_eh_personality, EMPTY_PANIC, etc.),
# and the resulting cross-archive duplicates are errors by default.
#
# The classic linker option `-multiply_defined,suppress` says "pick the
# first definition" (both copies are bit-identical so this is correct).
# Apple's new ld-prime accepts the flag silently but no longer honors it —
# we have to fall back to the classic linker with `-ld_classic` to get the
# behavior we want. `-ld_classic` is deprecated and announced for removal
# in a future Xcode; when that happens, the fix is a Rust workspace glue
# crate that compiles rype + sylph as one compilation graph (sharing std),
# at which point this whole block goes away.
if(APPLE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-ld_classic -Wl,-multiply_defined,suppress")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-ld_classic -Wl,-multiply_defined,suppress")
endif()

# Extensions that miint depends on (must load before miint)
# Load order matters for LoadAllExtensions: macros reference functions from
# core_functions (FLOOR, map_from_entries) and json (read_json).
# Skip during tidy checks - these extensions have their own dependencies (e.g., CURL)
# that aren't available in the tidy CI environment, and tidy only analyzes our source code
if(NOT CLANG_TIDY)
    duckdb_extension_load(core_functions)
    duckdb_extension_load(json)
    # httpfs disabled: upstream API incompatibility between DuckDB v1.5-variegata
    # and httpfs main (GetHTTPUtil return type mismatch). Users can still
    # INSTALL httpfs; LOAD httpfs; at runtime for HTTP URL support.
    # Re-enable when httpfs is updated for the v1.5 API.
    # duckdb_extension_load(httpfs
    #     GIT_URL https://github.com/duckdb/duckdb-httpfs
    #     GIT_TAG main
    # )
endif()

# Extension from this repo
duckdb_extension_load(miint
    SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}
    LOAD_TESTS
)
