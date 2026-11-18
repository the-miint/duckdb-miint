#pragma once

namespace miint {

// Keeps a statically linked LLVM libomp from shutting down while the process
// runs. Idempotent; LoadInternal calls it on every load.
//
// The macOS build links libomp.a (cmake/openmp.cmake). That runtime shuts down
// when the last thread that used OpenMP exits and initializes again on the next
// call, and the shutdown frees what an `omp critical` (unnamed, named, or the one
// behind a lock-based reduction) left in a variable in the caller's own binary:
// the variable keeps the stale pointer and the next critical section crashes
// (the-miint/duckdb-miint#482). DuckDB runs table functions on whichever thread
// executes the pipeline, and closing a database ends its workers, so the first
// thread to use OpenMP is often one that exits.
//
// One thread that has used OpenMP and never exits is enough for the runtime to
// stay up. This starts that thread the first time it is called and returns once
// the thread has registered with the runtime. Without a static libomp
// (MIINT_OPENMP_STATIC_LIBOMP unset: Linux's libgomp, a shared libomp, no
// OpenMP) it does nothing, since those runtimes do not shut down this way.
void KeepOpenMPRuntimeAlive();

} // namespace miint
