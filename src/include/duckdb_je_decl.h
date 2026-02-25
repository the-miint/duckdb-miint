/*
 * Declarations of DuckDB's prefixed jemalloc symbols for use in C++ code.
 *
 * minimap2 is compiled with malloc/free redirected to duckdb_je_* via
 * mm_alloc.h (force-included with -include). C++ code that frees memory
 * allocated by minimap2 must call duckdb_je_free (not system free).
 *
 * In the extension .so, these resolve to DuckDB's statically-linked jemalloc.
 * In the test binary, they resolve to stubs in mm_alloc_stubs.c that forward
 * to system malloc.
 */
#ifndef DUCKDB_JE_DECL_H
#define DUCKDB_JE_DECL_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void *duckdb_je_malloc(size_t size);
extern void *duckdb_je_calloc(size_t num, size_t size);
extern void *duckdb_je_realloc(void *ptr, size_t size);
extern void duckdb_je_free(void *ptr);

#ifdef __cplusplus
}
#endif

#endif /* DUCKDB_JE_DECL_H */
