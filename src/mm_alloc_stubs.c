/*
 * Stubs for duckdb_je_* symbols used by minimap2 when built with
 * -DMIINT_USE_JEMALLOC. In the extension .so these resolve to DuckDB's
 * jemalloc; in the standalone test binary they forward to system malloc.
 *
 * WARNING: Because both allocation and deallocation use libc here, the test
 * binary cannot detect allocator mismatch bugs (e.g. malloc via jemalloc but
 * free via libc). Always validate the extension binary under AddressSanitizer
 * to catch such issues.
 */
#include <stdlib.h>

void *duckdb_je_malloc(size_t size) {
	return malloc(size);
}
void *duckdb_je_calloc(size_t num, size_t size) {
	return calloc(num, size);
}
void *duckdb_je_realloc(void *ptr, size_t size) {
	return realloc(ptr, size);
}
void duckdb_je_free(void *ptr) {
	free(ptr);
}
