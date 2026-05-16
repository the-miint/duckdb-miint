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
// compiled single-threaded with no OpenMP runtime to mutate. Preserve the
// mutex so the libssu/skbb global-state serialization contract still holds.
OmpThreadScope::OmpThreadScope(int n_threads) : lock_(GlobalLibssuMutex()), prev_threads_(1) {
	if (n_threads < 1) {
		throw std::invalid_argument("OmpThreadScope: n_threads must be >= 1");
	}
}

OmpThreadScope::~OmpThreadScope() = default;
#else
OmpThreadScope::OmpThreadScope(int n_threads) : lock_(GlobalLibssuMutex()), prev_threads_(omp_get_max_threads()) {
	if (n_threads < 1) {
		throw std::invalid_argument("OmpThreadScope: n_threads must be >= 1");
	}
	omp_set_num_threads(n_threads);
}

OmpThreadScope::~OmpThreadScope() {
	omp_set_num_threads(prev_threads_);
}
#endif

} // namespace miint::unifrac
