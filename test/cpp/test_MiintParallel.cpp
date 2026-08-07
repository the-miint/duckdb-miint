#include <catch2/catch_test_macros.hpp>

#include "miint_parallel.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using miint::ParallelFor;
using miint::WorkerPool;

namespace {} // namespace

TEST_CASE("WorkerPool reports its budget and rejects a nonsensical one", "[parallel]") {
	// The budget INCLUDES the calling thread, so one is a legitimate request and
	// zero is not. Threads() can come back below what was asked for -- a build
	// without threads, or a pids cap, degrades to serial rather than failing -- so
	// the assertion is a range, not equality.
	REQUIRE(WorkerPool(1).Threads() == 1);
	REQUIRE(WorkerPool(4).Threads() >= 1);
	REQUIRE(WorkerPool(4).Threads() <= 4);
	REQUIRE_THROWS_AS(WorkerPool(0), std::invalid_argument);
	REQUIRE_THROWS_AS(WorkerPool(-1), std::invalid_argument);
}

TEST_CASE("ParallelFor visits every index exactly once", "[parallel]") {
	// Sizes deliberately straddle the thread counts -- fewer indices than threads,
	// exactly as many, and one more -- because that is where an off-by-one in the
	// chunk arithmetic hides.
	for (const int threads : {1, 2, 3, 4, 8}) {
		WorkerPool pool(threads);
		for (const int64_t n : {int64_t {0}, int64_t {1}, int64_t {7}, int64_t {8}, int64_t {9}, int64_t {1000}}) {
			std::vector<int> seen(static_cast<size_t>(n), 0);
			ParallelFor(&pool, 0, n, [&](int64_t begin, int64_t end) {
				for (int64_t i = begin; i < end; ++i) {
					seen[static_cast<size_t>(i)]++;
				}
			});
			REQUIRE(std::count(seen.begin(), seen.end(), 1) == static_cast<ptrdiff_t>(n));
		}
	}
}

TEST_CASE("ParallelFor's chunks are a contiguous ascending partition", "[parallel]") {
	// Stronger than "every index once", and the difference matters: the
	// decomposition mmvec relies on writes one slot per index and then sums those
	// slots in ASCENDING index order on a single thread. That reproduces the serial
	// sum only if the chunks tile the range with no gap and no overlap. A scheme
	// that interleaved indices round-robin would satisfy the previous test and
	// break this one.
	WorkerPool pool(4);
	// Collected under a lock because the chunks arrive from several threads, then
	// sorted -- the arrival ORDER is exactly the thing that is allowed to vary.
	std::mutex m;
	std::vector<std::pair<int64_t, int64_t>> ranges;
	ParallelFor(&pool, 5, 100, [&](int64_t begin, int64_t end) {
		std::lock_guard<std::mutex> lock(m);
		ranges.emplace_back(begin, end);
	});
	std::sort(ranges.begin(), ranges.end());
	REQUIRE_FALSE(ranges.empty());
	REQUIRE(ranges.front().first == 5);
	REQUIRE(ranges.back().second == 100);
	for (size_t i = 0; i < ranges.size(); ++i) {
		REQUIRE(ranges[i].first < ranges[i].second); // no empty chunk is handed out
		if (i > 0) {
			REQUIRE(ranges[i].first == ranges[i - 1].second);
		}
	}
}

TEST_CASE("a budget of one thread runs inline on the caller", "[parallel]") {
	// Not a micro-optimization. This is the Emscripten path and the path taken when
	// DuckDB gives the operator one thread, and it must create no OS thread at all:
	// a wasm build without -pthread throws from the std::thread constructor.
	WorkerPool pool(1);
	REQUIRE(pool.Threads() == 1);

	std::thread::id ran_on;
	ParallelFor(&pool, 0, 100, [&](int64_t, int64_t) { ran_on = std::this_thread::get_id(); });
	REQUIRE(ran_on == std::this_thread::get_id());
}

TEST_CASE("the calling thread is one of the workers", "[parallel]") {
	// A pool of T threads plus the caller is T+1 runnable threads on T cores. The
	// caller has to take a chunk itself, or every dispatch oversubscribes by one --
	// which is the shape of the 27x regression P4.5 measured.
	WorkerPool pool(4);
	std::mutex m;
	std::set<std::thread::id> ids;
	ParallelFor(&pool, 0, 4000, [&](int64_t, int64_t) {
		std::lock_guard<std::mutex> lock(m);
		ids.insert(std::this_thread::get_id());
	});
	REQUIRE(ids.count(std::this_thread::get_id()) == 1);
	REQUIRE(ids.size() <= 4);
}

TEST_CASE("an exception in the body reaches the caller, and the pool survives it", "[parallel]") {
	// A worker that lets an exception escape its thread function calls
	// std::terminate and takes the process with it, so this is a crash-or-error
	// test, not a tidiness one. Chunk 0 is the CALLER's, so throwing only for
	// begin > 0 guarantees the exception really did cross a thread boundary. On a
	// build with no threads at all (wasm without -pthread) there is no such
	// boundary, so the same intent is expressed by throwing from the only chunk
	// there is -- rather than letting the case silently stop being tested.
	WorkerPool pool(4);
	const bool threaded = pool.Threads() > 1;
	REQUIRE_THROWS_AS(ParallelFor(&pool, 0, 4000,
	                              [threaded](int64_t begin, int64_t) {
		                              if (!threaded || begin > 0) {
			                              throw std::runtime_error("boom");
		                              }
	                              }),
	                  std::runtime_error);

	// A pool that lost a worker or left its bookkeeping inconsistent on the way out
	// would hang here rather than fail.
	std::atomic<int64_t> total {0};
	ParallelFor(&pool, 0, 100, [&](int64_t begin, int64_t end) { total += end - begin; });
	REQUIRE(total.load() == 100);
}

TEST_CASE("the exception the caller sees is the lowest-numbered failing chunk", "[parallel]") {
	// Determinism has to extend to failures too. If the surviving exception were
	// whichever thread threw first, the same bad input would produce different
	// error messages from run to run and a user could not reproduce a report.
	//
	// Simply letting every chunk throw does NOT test this: chunk 0 runs on the
	// calling thread, which is already executing while the workers are still waking
	// from a condition variable, so chunk 0 wins on time as well as on number and a
	// first-thrower-wins implementation passes. Verified -- that mutation survived.
	// So chunk 0 is held back until another chunk has definitely thrown, which
	// forces a real choice between "earliest" and "lowest".
	WorkerPool pool(4);
	const bool threaded = pool.Threads() > 1;
	bool always_lowest = true;
	for (int rep = 0; rep < 20; ++rep) {
		std::atomic<int> others_thrown {0};
		try {
			ParallelFor(&pool, 0, 4000, [&](int64_t begin, int64_t) {
				if (begin > 0) {
					others_thrown.fetch_add(1);
					throw std::runtime_error(std::to_string(begin));
				}
				// Spin, not sleep: portable, and bounded by another chunk's
				// progress rather than by a duration. Guaranteed to end whenever
				// there is more than one chunk, since the workers own those.
				while (threaded && others_thrown.load() == 0) {
				}
				throw std::runtime_error("0");
			});
			always_lowest = false; // unreachable: every chunk throws
		} catch (const std::runtime_error &e) {
			always_lowest = always_lowest && std::string(e.what()) == "0";
		}
	}
	REQUIRE(always_lowest);
}

TEST_CASE("a nested ParallelFor runs inline instead of deadlocking", "[parallel]") {
	// Every worker is already inside the outer body, so a nested dispatch that
	// waited on them would wait forever. Running inline is a defined answer;
	// hanging is not. The assertion is that the work still all happens.
	WorkerPool pool(4);
	std::vector<int> seen(400, 0);
	ParallelFor(&pool, 0, 400, [&](int64_t begin, int64_t end) {
		ParallelFor(&pool, begin, end, [&](int64_t inner_begin, int64_t inner_end) {
			for (int64_t i = inner_begin; i < inner_end; ++i) {
				seen[static_cast<size_t>(i)]++;
			}
		});
	});
	REQUIRE(std::count(seen.begin(), seen.end(), 1) == 400);
}

TEST_CASE("scheduling cannot change a floating-point result", "[parallel]") {
	// The property the whole milestone rests on, stated as a test rather than as an
	// argument.
	constexpr int64_t kN = 10000;
	const auto value = [](int64_t i) {
		// 1.0 followed by values far below its ulp: adding them up first accumulates
		// a quantity that survives, adding them into 1.0 one at a time loses each
		// one. So the sum genuinely depends on the order.
		return i == 0 ? 1.0 : 1e-18 * static_cast<double>(i % 7 + 1);
	};

	std::vector<double> serial(static_cast<size_t>(kN));
	for (int64_t i = 0; i < kN; ++i) {
		serial[static_cast<size_t>(i)] = value(i);
	}
	double expected = 0.0;
	for (int64_t i = 0; i < kN; ++i) {
		expected += serial[static_cast<size_t>(i)];
	}

	// Guard the guard. If these values summed to the same double in any order, the
	// equality below would also pass for a decomposition that DID let a reduction
	// cross a thread, and the test would be worthless.
	double reversed = 0.0;
	for (int64_t i = kN - 1; i >= 0; --i) {
		reversed += serial[static_cast<size_t>(i)];
	}
	REQUIRE(reversed != expected);

	for (const int threads : {1, 2, 3, 4, 8}) {
		WorkerPool pool(threads);
		std::vector<double> slots(static_cast<size_t>(kN), 0.0);
		ParallelFor(&pool, 0, kN, [&](int64_t begin, int64_t end) {
			for (int64_t i = begin; i < end; ++i) {
				slots[static_cast<size_t>(i)] = value(i);
			}
		});
		double got = 0.0;
		for (int64_t i = 0; i < kN; ++i) {
			got += slots[static_cast<size_t>(i)];
		}
		REQUIRE(got == expected); // bit-identical, not within a tolerance
	}
}

TEST_CASE("an empty or inverted range does no work", "[parallel]") {
	WorkerPool pool(4);
	std::atomic<int> calls {0};
	ParallelFor(&pool, 10, 10, [&](int64_t, int64_t) { ++calls; });
	ParallelFor(&pool, 10, 5, [&](int64_t, int64_t) { ++calls; });
	REQUIRE(calls.load() == 0);
}

TEST_CASE("one pool serves many dispatches", "[parallel]") {
	// mmvec creates the pool once per fit and dispatches several times per objective
	// evaluation, a thousand evaluations deep, so the park/wake handshake has to
	// still be correct on the ten-thousandth call rather than only the first. A
	// generation counter that could be missed, or a worker that parked without
	// re-checking its predicate, shows up here as a hang or a dropped chunk.
	WorkerPool pool(4);
	bool all_ok = true;
	for (int rep = 0; rep < 2000; ++rep) {
		std::vector<int> seen(64, 0);
		ParallelFor(&pool, 0, 64, [&](int64_t begin, int64_t end) {
			for (int64_t i = begin; i < end; ++i) {
				seen[static_cast<size_t>(i)]++;
			}
		});
		all_ok = all_ok && std::count(seen.begin(), seen.end(), 1) == 64;
	}
	REQUIRE(all_ok);
}
