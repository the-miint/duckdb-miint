#include "openmp_runtime_keepalive.hpp"

#ifdef MIINT_OPENMP_STATIC_LIBOMP
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>

// Declared rather than pulled in from <omp.h>: miint's own TUs are compiled
// without the OpenMP flag, so omp.h is not on their include path. The symbol
// comes from the libomp.a cmake/openmp.cmake links, and its signature is fixed by
// the OpenMP ABI.
extern "C" int omp_get_max_threads(void);
#endif

namespace miint {

void KeepOpenMPRuntimeAlive() {
#ifdef MIINT_OPENMP_STATIC_LIBOMP
	static std::once_flag started;
	std::call_once(started, [] {
		std::mutex mutex;
		std::condition_variable registered_cv;
		bool registered = false;
		// Detached: it has to outlive every database, and a joinable std::thread
		// destroyed at exit would call std::terminate.
		std::thread([&] {
			// Any OpenMP entry point registers the calling thread with the
			// runtime; this one has no side effects.
			omp_get_max_threads();
			{
				// Notified under the lock: the caller cannot return, and take
				// mutex, registered_cv and registered down with its frame, until
				// this scope has released it.
				std::lock_guard<std::mutex> lock(mutex);
				registered = true;
				registered_cv.notify_one();
			}
			// Nothing below may touch the caller's locals.
			for (;;) {
				std::this_thread::sleep_for(std::chrono::hours(24));
			}
		}).detach();
		std::unique_lock<std::mutex> lock(mutex);
		registered_cv.wait(lock, [&] { return registered; });
	});
#endif
}

} // namespace miint
