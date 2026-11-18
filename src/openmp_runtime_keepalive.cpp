#include "openmp_runtime_keepalive.hpp"

#ifdef MIINT_OPENMP_STATIC_LIBOMP
#include <chrono>
#include <future>
#include <mutex>
#include <thread>

#include <omp.h>
#endif

namespace miint {

void KeepOpenMPRuntimeAlive() {
#ifdef MIINT_OPENMP_STATIC_LIBOMP
	static std::once_flag started;
	std::call_once(started, [] {
		// The thread owns the promise and the caller holds a future on the same
		// reference-counted state, so neither touches memory the other can free.
		std::promise<void> registered;
		std::future<void> registered_future = registered.get_future();
		// Detached: it has to outlive every database, and a joinable std::thread
		// destroyed at exit would call std::terminate.
		std::thread([registered = std::move(registered)]() mutable {
			// Any OpenMP entry point registers the calling thread with the
			// runtime; this one has no side effects.
			omp_get_max_threads();
			registered.set_value();
			for (;;) {
				std::this_thread::sleep_for(std::chrono::hours(24));
			}
		}).detach();
		registered_future.wait();
	});
#endif
}

} // namespace miint
