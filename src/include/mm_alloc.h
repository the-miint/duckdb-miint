/*
 * Redirects malloc/calloc/realloc/free to DuckDB's prefixed jemalloc symbols.
 *
 * This header is force-included into ALL minimap2 .c files via:
 *   CPPFLAGS="-DMIINT_USE_JEMALLOC -include mm_alloc.h"
 * in the ExternalProject_Add BUILD_COMMAND. It requires zero modifications
 * to minimap2 source files and survives submodule updates.
 *
 * The redirect covers all allocation paths in minimap2, including kalloc.c's
 * fallback to libc when km==NULL (e.g. krealloc with km=NULL calls realloc,
 * which this header redirects to duckdb_je_realloc at compile time).
 *
 * NOTE: posix_memalign/aligned_alloc are NOT redirected. minimap2 does not
 * use these directly (kalloc handles aligned allocations internally), but if
 * a future minimap2 version adds them, this header would need updating.
 */
#ifndef MM_ALLOC_H
#define MM_ALLOC_H

#ifdef __cplusplus
#error "mm_alloc.h must not be included in C++ code — use duckdb_je_decl.h instead"
#endif

#ifdef MIINT_USE_JEMALLOC

#include <stdlib.h> /* Process stdlib declarations before redefining */

extern void *duckdb_je_malloc(size_t size);
extern void *duckdb_je_calloc(size_t num, size_t size);
extern void *duckdb_je_realloc(void *ptr, size_t size);
extern void duckdb_je_free(void *ptr);

#define malloc(size)        duckdb_je_malloc(size)
#define calloc(count, size) duckdb_je_calloc(count, size)
#define realloc(ptr, size)  duckdb_je_realloc(ptr, size)
#define free(ptr)           duckdb_je_free(ptr)

#endif /* MIINT_USE_JEMALLOC */
#endif /* MM_ALLOC_H */
