#ifdef MIINT_HAS_HMMER
// The configuration HMMER and Easel are compiled with, pinned at compile time.
// cmake/deps/hmmer.cmake generates p7_config.h and esl_config.h; these checks
// fail the build if a setting that changes results, or that this build cannot
// support, ever turns on.
//
// Each header is checked right after it is included: both are autoconf
// templates whose unset settings remain as live "#undef NAME" lines, so
// including esl_config.h undefines, for instance, an HAVE_FLUSH_ZERO_MODE that
// p7_config.h had defined.
#include <catch2/catch_test_macros.hpp>

#include <string>

// Flush-to-zero and denormals-are-zero change floating-point results, so scores
// would no longer match the HMMER command-line programs, and they are set per
// thread, so they would leak into DuckDB's worker threads.
//
// HMMER_THREADS changes the layout of P7_HMMFILE, and the library's hot path
// needs no threads of its own: DuckDB runs it on its worker threads.

#include "p7_config.h"
#if defined(HAVE_FLUSH_ZERO_MODE)
#error "p7_config.h: HMMER must be built with FTZ off"
#endif
#if defined(HMMER_THREADS) || defined(HMMER_MPI)
#error "p7_config.h: HMMER must be built without its own threading or MPI"
#endif
#if defined(eslENABLE_SSE) == defined(eslENABLE_NEON)
#error "p7_config.h: exactly one of eslENABLE_SSE and eslENABLE_NEON must be defined"
#endif

#include "esl_config.h"
#if defined(HAVE_FLUSH_ZERO_MODE) || defined(HAVE_DENORMALS_ZERO_MODE)
#error "esl_config.h: Easel must be built with FTZ/DAZ off"
#endif
#if defined(HAVE_PTHREAD) || defined(HAVE_MPI)
#error "esl_config.h: Easel must be built without threads or MPI"
#endif
// HMMER never opens a file or pipe here: models are parsed from memory.
#if defined(HAVE_POPEN)
#error "esl_config.h: Easel must be built without popen"
#endif
#if defined(eslENABLE_SSE) == defined(eslENABLE_NEON)
#error "esl_config.h: exactly one of eslENABLE_SSE and eslENABLE_NEON must be defined"
#endif

// HMMER 3.4 has no scalar implementation: the SIMD implementation must match
// the target.
#if (defined(__aarch64__) || defined(__arm64__)) && !defined(eslENABLE_NEON)
#error "arm64 builds use HMMER's NEON implementation"
#endif
#if (defined(__x86_64__) || defined(__EMSCRIPTEN__)) && !defined(eslENABLE_SSE)
#error "x86-64 and wasm builds use HMMER's SSE implementation"
#endif

// The fork's limit on model length (upstream #350): the parsers reject a larger
// LENG before allocating anything.
static_assert(p7_MAXM == 100000, "p7_MAXM comes from the fork's p7_config.h.in");

TEST_CASE("HMMER config - versions of the vendored sources", "[hmmer]") {
	REQUIRE(std::string(HMMER_VERSION) == "3.4");
	REQUIRE(std::string(EASEL_VERSION) == "0.49");
}
#endif
