#pragma once
//
// A persistent, condition-variable-parked worker pool, and two ways to hand it
// work whose decomposition is deterministic.
//
// Written for mmvec's fit, whose requirement is unusual and shapes everything
// here: parallelism must change **no value**. The objective's hot loops are
// decomposed so that no floating-point reduction ever crosses a thread -- each
// thread writes its own slots, and the ordered sum happens afterwards on one
// thread -- so a fit is bit-identical at any thread count. That only holds if the
// chunks are a CONTIGUOUS ASCENDING PARTITION of the range, which is the one
// property of this file that callers depend on and tests pin down.
//
// Three choices, each measured in M6 P4.5 / M8 P0 rather than assumed:
//
//   * `std::thread` over OpenMP. `miint_openmp` exists only under
//     `if(MIINT_ENABLE_UNIFRAC AND NOT EMSCRIPTEN)`, so an OpenMP dependency would
//     make mmvec's parallelism vanish whenever UniFrac is off, and the ISA variant
//     objects never receive `-fopenmp` at all.
//   * **No spin-waiting.** Workers park on a condition variable; a spin-wait pool
//     measured 27x SLOWER once it oversubscribed a cpuset.
//   * **The calling thread takes a chunk**, or a pool of T plus the caller would
//     be T+1 runnable threads on T cores.
//
// A null `WorkerPool *` means "run serially", and it is the SAME state as a pool
// with a budget of one -- both run the body inline on the caller. Callers pass
// their pool pointer straight through without testing it.

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace miint {

//! A fixed set of worker threads, parked until there is work.
//!
//! Construct one per operation, not per call: the park/wake handshake is cheap
//! but thread creation is not, and mmvec dispatches thousands of times per fit.
//!
//! Not copyable, not movable, and not itself thread-safe -- one owner drives it.
class WorkerPool {
public:
	//! `n_threads` is the TOTAL budget INCLUDING the calling thread, so 1 means
	//! "run everything inline and create nothing". Throws std::invalid_argument
	//! below 1.
	//!
	//! Fewer threads than requested is not an error: a build without threads, or a
	//! host at its pids limit, degrades to serial. Ask `Threads()` what was
	//! actually obtained.
	explicit WorkerPool(int n_threads);
	~WorkerPool();

	WorkerPool(const WorkerPool &) = delete;
	WorkerPool &operator=(const WorkerPool &) = delete;

	//! The budget actually available, always >= 1 and never above the request.
	int Threads() const {
		return static_cast<int>(workers_.size()) + 1; // + the calling thread
	}

	//! How many times work has actually been handed to the worker threads.
	//!
	//! Counts only the threaded path -- an inline dispatch (budget of 1, a single
	//! chunk, a nested call) does not increment it. That is deliberate and is what
	//! makes it useful: it is the only way a test can tell "parallelised correctly"
	//! from "silently ran serially and therefore matched trivially", which
	//! bit-identity assertions by their nature cannot distinguish.
	uint64_t Dispatches() const {
		return dispatches_.load(std::memory_order_relaxed);
	}

	//! Run `chunk(c)` for every c in [0, n_chunks) and return once all have
	//! finished, rethrowing on the caller if any threw.
	//!
	//! The primitive both `ParallelFor` and `ParallelInvoke` are built on, and
	//! public because those two do not exhaust what it expresses: it is simply "run
	//! N independent things", with no assumption that they are slices of a range.
	void Dispatch(int n_chunks, const std::function<void(int)> &chunk);

private:
	//! Claim and run chunks until none are left. Called with `lock` HELD; releases
	//! it around each chunk and re-takes it. Shared by the caller and the workers so
	//! the next_chunk_/outstanding_ invariant is written once.
	void DrainChunks(std::unique_lock<std::mutex> &lock);

	//! Invoke one chunk, catching whatever it throws. Called with `m_` NOT held.
	void RunChunk(int c);

	void WorkerLoop();

	std::mutex m_;
	std::condition_variable work_ready_; //!< workers wait here for a new generation
	std::condition_variable all_done_;   //!< the caller waits here for outstanding_ == 0
	std::vector<std::thread> workers_;

	const std::function<void(int)> *chunk_ = nullptr;
	uint64_t generation_ = 0; //!< bumped per dispatch, so a worker cannot miss one
	int n_chunks_ = 0;
	int next_chunk_ = 0;
	int outstanding_ = 0;

	//! The surviving exception is the one from the LOWEST-numbered chunk, not the
	//! one that happened to be thrown first. Otherwise the same bad input would
	//! report differently from run to run and a user could not reproduce it.
	std::exception_ptr failure_;
	int failed_chunk_ = 0;

	bool stop_ = false;
	std::atomic<uint64_t> dispatches_ {0};
};

//! Split `[begin, end)` into `pool->Threads()` contiguous chunks and invoke
//! `fn(chunk_begin, chunk_end)` on each, returning once every chunk has finished.
//!
//! The chunks are a contiguous ascending partition and their boundaries depend
//! only on `begin`, `end` and the budget -- never on timing. An index therefore
//! lands in the same chunk on every run, which is what lets a caller write one
//! slot per index and reduce those slots afterwards in index order for a
//! bit-identical result.
//!
//! Runs inline on the calling thread, waking nothing, when `pool` is null, the
//! budget is 1, or the range holds fewer than two elements. A TEMPLATE rather than
//! a `std::function` parameter specifically so that path costs nothing: the
//! minibatch objective calls this with a null pool on every one of ~196000 Adam
//! updates, and type-erasing the body there was measured to heap-allocate a
//! closure per call for no parallelism at all.
//!
//! An exception from `fn` propagates to the caller once all chunks have finished --
//! never out of a worker's thread function, which would call std::terminate.
template <typename Fn>
void ParallelFor(WorkerPool *pool, int64_t begin, int64_t end, Fn &&fn) {
	const int64_t n = end - begin;
	if (n <= 0) {
		return;
	}
	if (pool == nullptr || pool->Threads() <= 1 || n == 1) {
		fn(begin, end);
		return;
	}
	// Never more chunks than indices, so every chunk handed out is non-empty.
	const int64_t chunks = std::min<int64_t>(pool->Threads(), n);
	const int64_t base = n / chunks;
	const int64_t remainder = n % chunks;
	pool->Dispatch(static_cast<int>(chunks), [&](int c) {
		// Spread the remainder over the first `remainder` chunks. Written as
		// products of the chunk index rather than `n * c / chunks` so it cannot
		// overflow for a large range, and so the boundaries are exact integers
		// rather than the result of two roundings.
		const int64_t lo = begin + c * base + std::min<int64_t>(c, remainder);
		const int64_t hi = lo + base + (c < remainder ? 1 : 0);
		fn(lo, hi);
	});
}

//! Run two independent pieces of work at the same time, each WHOLE.
//!
//! Not a special case of ParallelFor with two indices, even though that is how it
//! is implemented -- the distinction is the point. `ParallelFor` divides one
//! computation, which for an Eigen matrix product changes its result (M8 P0
//! measured Eigen re-blocking a sliced product and moving the fitted model). This
//! divides nothing: each callable computes exactly what it would have computed
//! alone, so the result is bit-identical by construction rather than by
//! measurement. That is what makes it usable on the two gradient GEMMs.
//!
//! Runs `a()` then `b()` inline when there is no second thread to put one on.
template <typename FnA, typename FnB>
void ParallelInvoke(WorkerPool *pool, FnA &&a, FnB &&b) {
	if (pool == nullptr || pool->Threads() <= 1) {
		a();
		b();
		return;
	}
	pool->Dispatch(2, [&](int task) {
		if (task == 0) {
			a();
		} else {
			b();
		}
	});
}

} // namespace miint
