# This file is included by DuckDB's build system. It specifies which extension to load

# Extensions that miint depends on (must load before miint)
duckdb_extension_load(json)
duckdb_extension_load(httpfs
    GIT_URL https://github.com/duckdb/duckdb-httpfs
    GIT_TAG main
)

# Extension from this repo
duckdb_extension_load(miint
    SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}
    LOAD_TESTS
)
