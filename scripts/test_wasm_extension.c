// Minimal DuckDB WASM test harness for verifying extension loading.
// Compiled with emcc as MAIN_MODULE, loads the miint extension as a
// SIDE_MODULE, and runs basic queries to verify it works.
//
// Build: see scripts/build_wasm.sh --verify-load
// Run:   node test_wasm_harness.js

#include <stdio.h>
#include <string.h>
#include "duckdb.h"

static int run_query_check(duckdb_connection conn, const char *sql, const char *expected_prefix) {
    duckdb_result result;
    if (duckdb_query(conn, sql, &result) != DuckDBSuccess) {
        fprintf(stderr, "  FAIL: %s\n", duckdb_result_error(&result));
        duckdb_destroy_result(&result);
        return 1;
    }
    const char *val = duckdb_value_varchar(&result, 0, 0);
    if (expected_prefix && strncmp(val, expected_prefix, strlen(expected_prefix)) != 0) {
        fprintf(stderr, "  FAIL: expected prefix '%s', got '%s'\n", expected_prefix, val);
        duckdb_free((void *)val);
        duckdb_destroy_result(&result);
        return 1;
    }
    printf("  OK: %s\n", val);
    duckdb_free((void *)val);
    duckdb_destroy_result(&result);
    return 0;
}

int main(int argc, char **argv) {
    printf("=== miint WASM Extension Test ===\n\n");

    duckdb_database db;
    duckdb_connection conn;
    duckdb_config config;

    // Allow unsigned extensions (our local build isn't signed)
    duckdb_create_config(&config);
    duckdb_set_config(config, "allow_unsigned_extensions", "true");

    if (duckdb_open_ext(NULL, &db, config, NULL) != DuckDBSuccess) {
        fprintf(stderr, "Failed to open database\n");
        return 1;
    }
    duckdb_destroy_config(&config);
    duckdb_connect(db, &conn);

    // Check DuckDB version
    printf("1. DuckDB version\n");
    if (run_query_check(conn, "SELECT version()", "v")) return 1;

    // Load json dependency first (miint macros reference read_json)
    printf("\n2. Loading miint extension\n");
    duckdb_result result;
    // Load json dependency first (miint macros reference read_json)
    if (duckdb_query(conn, "LOAD '/extensions/json.duckdb_extension'", &result) != DuckDBSuccess) {
        fprintf(stderr, "  FAIL loading json: %s\n", duckdb_result_error(&result));
        duckdb_destroy_result(&result);
        return 1;
    }
    duckdb_destroy_result(&result);
    printf("  json extension loaded\n");
    if (duckdb_query(conn, "LOAD '/extensions/miint.duckdb_extension'", &result) != DuckDBSuccess) {
        fprintf(stderr, "  FAIL: %s\n", duckdb_result_error(&result));
        duckdb_destroy_result(&result);
        return 1;
    }
    duckdb_destroy_result(&result);
    printf("  OK: extension loaded\n");

    // Test miint_version()
    printf("\n3. miint_version()\n");
    if (run_query_check(conn, "SELECT miint_version()", NULL)) return 1;

    printf("\n=== All tests passed ===\n");

    duckdb_disconnect(&conn);
    duckdb_close(&db);
    return 0;
}
