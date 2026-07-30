#include "unifrac_omp_scope.hpp"

#include <stdexcept>

#ifndef __EMSCRIPTEN__
#include <omp.h>
#endif

namespace miint::unifrac {

namespace {

std::mutex &GlobalLibssuMutex() {
	static std::mutex m;
	return m;
}

} // namespace

#ifdef __EMSCRIPTEN__
// WASM builds skip find_package(OpenMP); libssu_wasm.a / libskbb_wasm.a are
// compiled single-threaded with no OpenMP runtime to mutate. OmpThreadScope still
// takes the mutex, so the libssu/skbb global-state serialization contract holds.
OmpThreadPin::OmpThreadPin(int n_threads, OmpScopeHeld) : prev_threads_(1) {
	if (n_threads < 1) {
		throw std::invalid_argument("OmpThreadPin: n_threads must be >= 1");
	}
}

OmpThreadPin::~OmpThreadPin() = default;
#else
OmpThreadPin::OmpThreadPin(int n_threads, OmpScopeHeld) : prev_threads_(omp_get_max_threads()) {
	if (n_threads < 1) {
		throw std::invalid_argument("OmpThreadPin: n_threads must be >= 1");
	}
	omp_set_num_threads(n_threads);
}

OmpThreadPin::~OmpThreadPin() {
	omp_set_num_threads(prev_threads_);
}
#endif

// Lock first, then pin: the pin is what a concurrent caller would race on, so it
// is only touched under the lock. Members destruct in reverse, restoring the
// thread count before the lock is released.
OmpThreadScope::OmpThreadScope(int n_threads) : lock_(GlobalLibssuMutex()), pin_(n_threads, OmpScopeHeld()) {
}

OmpThreadScope::~OmpThreadScope() = default;

} // namespace miint::unifrac
