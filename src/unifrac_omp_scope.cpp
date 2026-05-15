#include "unifrac_omp_scope.hpp"

#include <stdexcept>

#include <omp.h>

namespace miint::unifrac {

namespace {

std::mutex &GlobalLibssuMutex() {
	static std::mutex m;
	return m;
}

} // namespace

OmpThreadScope::OmpThreadScope(int n_threads) : lock_(GlobalLibssuMutex()), prev_threads_(omp_get_max_threads()) {
	if (n_threads < 1) {
		throw std::invalid_argument("OmpThreadScope: n_threads must be >= 1");
	}
	omp_set_num_threads(n_threads);
}

OmpThreadScope::~OmpThreadScope() {
	omp_set_num_threads(prev_threads_);
}

} // namespace miint::unifrac
