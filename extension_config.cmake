# This file is included by DuckDB's build system. It specifies which extension to load

# The rype + sylph duplicate-Rust-std issue used to require platform-specific
# linker workarounds here (-Wl,--allow-multiple-definition on GNU ld,
# -Wl,-ld_classic -Wl,-multiply_defined,suppress on Apple). Both crates are
# now bundled through the `miint_rust_glue` umbrella staticlib (see
# ext/miint-rust-glue/), which gives them a single shared std and removes
# the duplicate symbols at their source — no link flag needed.

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
