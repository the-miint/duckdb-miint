#pragma once

namespace miint {

// Keeps a statically linked libomp from shutting down while the process runs;
// why, in docs/internals/embedded-tools.md, "The OpenMP runtime". Idempotent;
// LoadInternal calls it on every load.
//
// The first call starts one thread that makes an OpenMP call and never exits,
// and returns once that thread has registered with the runtime. Without
// MIINT_OPENMP_STATIC_LIBOMP (Linux's libgomp, a shared libomp, no OpenMP) it
// does nothing.
void KeepOpenMPRuntimeAlive();

} // namespace miint
