# This file is included by DuckDB's build system. It specifies which extension to load

# MinGW's ld treats duplicate COMDAT symbols (C++ constexpr/inline from DuckDB headers)
# as errors when the extension's static lib and DuckDB's static lib are linked into the
# same binary. This is set here (not in the extension's CMakeLists.txt) because this file
# runs in DuckDB's CMake scope before libduckdb.dll/duckdb.exe targets are configured.
if(MINGW)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--allow-multiple-definition")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--allow-multiple-definition")
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
