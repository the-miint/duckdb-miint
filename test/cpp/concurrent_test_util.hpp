#pragma once

#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace miint::test {

// Run `body` on `n_threads` threads that are all inside it at once, so a race
// between them has a chance to happen rather than being serialized by luck.
// Threads rendezvous with a bounded wait: a timeout degrades the test to "less
// overlap than intended", never to a hang.
//
// For any test whose point is that the threads overlap: a race it means to
// provoke, or work it means to see serialized.
inline void RunConcurrently(size_t n_threads, const std::function<void(size_t)> &body) {
	std::mutex mu;
	std::condition_variable cv;
	size_t arrived = 0;
	std::vector<std::thread> pool;
	pool.reserve(n_threads);
	for (size_t t = 0; t < n_threads; ++t) {
		pool.emplace_back([&, t]() {
			{
				std::unique_lock<std::mutex> lk(mu);
				++arrived;
				cv.notify_all();
				cv.wait_for(lk, std::chrono::seconds(5), [&] { return arrived == n_threads; });
			}
			body(t);
		});
	}
	for (auto &th : pool) {
		th.join();
	}
}

} // namespace miint::test
