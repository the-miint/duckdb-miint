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

    // Test UniFrac end-to-end: build a feature table + phylogeny inline (no
    // read_biom/read_newick — HDF5 and the filesystem are unavailable under
    // WASM) and compute Faith's phylogenetic diversity. Expected values are
    // hand-checked against test/sql/unifrac_faith_pd.test for the gg_otu tree
    // ((GG_OTU_1:0.1,GG_OTU_2:0.2)Int1:0.3,(GG_OTU_3:0.4,(GG_OTU_4:0.5,GG_OTU_5:0.6)Int2:0.7)Int3:0.8);
    //   Sample1 {OTU_2,OTU_4}      = 0.2+0.3+0.5+0.7+0.8 = 2.5
    //   Sample3 {OTU_1,OTU_3,OTU_4,OTU_5} = 0.1+0.3+0.4+0.5+0.6+0.7+0.8 = 3.4
    printf("\n4. UniFrac faith_pd (end-to-end libssu compute)\n");
    const char *setup[] = {
        "CREATE TABLE observations(sample_id VARCHAR, feature_id VARCHAR, value DOUBLE)",
        "INSERT INTO observations VALUES "
        "('Sample1','GG_OTU_2',1),('Sample1','GG_OTU_4',1),"
        "('Sample2','GG_OTU_2',1),('Sample2','GG_OTU_4',1),('Sample2','GG_OTU_5',1),"
        "('Sample3','GG_OTU_1',1),('Sample3','GG_OTU_3',1),('Sample3','GG_OTU_4',1),('Sample3','GG_OTU_5',1),"
        "('Sample4','GG_OTU_2',1),('Sample4','GG_OTU_3',1),"
        "('Sample5','GG_OTU_2',1),"
        "('Sample6','GG_OTU_2',1),('Sample6','GG_OTU_3',1),('Sample6','GG_OTU_4',1)",
        "CREATE TABLE tree(node_index BIGINT, name VARCHAR, branch_length DOUBLE, edge_id BIGINT, "
        "parent_index BIGINT, is_tip BOOLEAN)",
        "INSERT INTO tree VALUES "
        "(0,'GG_OTU_1',0.1,NULL,2,true),(1,'GG_OTU_2',0.2,NULL,2,true),(2,'Int1',0.3,NULL,8,false),"
        "(3,'GG_OTU_3',0.4,NULL,7,true),(4,'GG_OTU_4',0.5,NULL,6,true),(5,'GG_OTU_5',0.6,NULL,6,true),"
        "(6,'Int2',0.7,NULL,7,false),(7,'Int3',0.8,NULL,8,false),(8,'',NULL,NULL,NULL,false)",
    };
    for (size_t i = 0; i < sizeof(setup) / sizeof(setup[0]); i++) {
        if (duckdb_query(conn, setup[i], &result) != DuckDBSuccess) {
            fprintf(stderr, "  FAIL (setup): %s\n", duckdb_result_error(&result));
            duckdb_destroy_result(&result);
            return 1;
        }
        duckdb_destroy_result(&result);
    }
    // round() to 4dp so the check is robust to floating-point summation noise.
    if (run_query_check(conn,
                        "SELECT round(faith_pd, 4) FROM unifrac_faith_pd('observations', 'tree') "
                        "WHERE sample_id = 'Sample1'",
                        "2.5"))
        return 1;
    if (run_query_check(conn,
                        "SELECT round(faith_pd, 4) FROM unifrac_faith_pd('observations', 'tree') "
                        "WHERE sample_id = 'Sample3'",
                        "3.4"))
        return 1;

    printf("\n=== All tests passed ===\n");

    duckdb_disconnect(&conn);
    duckdb_close(&db);
    return 0;
}
