#include "unifrac_omp_scope.hpp"

#include <random>
#include <stdexcept>

#ifndef __EMSCRIPTEN__
#include <omp.h>
#endif

namespace miint::unifrac {

namespace {

// Seeds for callers that supplied none. Thread-local so concurrent calls draw
// without synchronization — the whole point is to keep them off the libraries'
// own process-global generators. Seeded once per thread from random_device, so
// an unseeded run stays non-reproducible exactly as it was when the library drew
// the seed itself. Masked to 31 bits because a negative value would send the
// call straight back to that global generator.
int DeriveSeed() {
	static thread_local std::mt19937 rng {std::random_device {}()};
	return static_cast<int>(rng() & 0x7fffffffu);
}

} // namespace

#ifdef __EMSCRIPTEN__
// WASM builds skip find_package(OpenMP); libssu_wasm.a / libskbb_wasm.a are
// compiled single-threaded with no OpenMP runtime to mutate.
OmpThreadPin::OmpThreadPin(int n_threads) : prev_threads_(1) {
	if (n_threads < 1) {
		throw std::invalid_argument("OmpThreadPin: n_threads must be >= 1");
	}
}

OmpThreadPin::~OmpThreadPin() = default;
#else
OmpThreadPin::OmpThreadPin(int n_threads) : prev_threads_(omp_get_max_threads()) {
	if (n_threads < 1) {
		throw std::invalid_argument("OmpThreadPin: n_threads must be >= 1");
	}
	omp_set_num_threads(n_threads);
}

OmpThreadPin::~OmpThreadPin() {
	omp_set_num_threads(prev_threads_);
}
#endif

ComputeCallScope::ComputeCallScope(int n_threads, int seed) : pin_(n_threads), seed_(seed >= 0 ? seed : DeriveSeed()) {
}

} // namespace miint::unifrac
