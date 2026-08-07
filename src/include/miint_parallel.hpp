#pragma once
//
// A persistent, condition-variable-parked worker pool, and a `ParallelFor` over
// it whose decomposition is deterministic.
//
// Written for mmvec's fit, whose requirement is unusual and shapes everything
// here: parallelism must change **no value**. The objective's hot loops are
// decomposed so that no floating-point reduction ever crosses a thread -- each
// thread writes its own slots, and the ordered sum happens afterwards on one
// thread -- so a fit is bit-identical at any thread count. That only holds if the
// chunks are a CONTIGUOUS ASCENDING PARTITION of the range, which is the one
// property of this file that callers depend on and tests pin down.
//
// Three deliberate choices, each measured rather than assumed:
//
//   * `std::thread` rather than OpenMP. `miint_openmp` exists only under
//     `if(MIINT_ENABLE_UNIFRAC AND NOT EMSCRIPTEN)`, so an OpenMP dependency would
//     make mmvec's parallelism vanish whenever UniFrac is off, and the ISA variant
//     objects do not receive `-fopenmp` at all.
//   * **No spin-waiting.** Workers park on a condition variable. A spin-wait pool
//     measured 27x SLOWER in M6 P4.5 once it oversubscribed a cpuset.
//   * **The calling thread takes a chunk.** A pool of T plus the caller would be
//     T+1 runnable threads on T cores.

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
		return n_threads_;
	}

private:
	friend void ParallelFor(WorkerPool &pool, int64_t begin, int64_t end,
	                        const std::function<void(int64_t, int64_t)> &fn);

	//! Run `chunk(c)` for every c in [0, n_chunks), and return once all have
	//! finished. Rethrows on the caller if any chunk threw.
	void Dispatch(int n_chunks, const std::function<void(int)> &chunk);

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
	int n_threads_ = 1;
};

//! Split `[begin, end)` into `pool.Threads()` contiguous chunks and invoke
//! `fn(chunk_begin, chunk_end)` on each, returning once every chunk has finished.
//!
//! The chunks are a contiguous ascending partition and their boundaries depend
//! only on `begin`, `end` and the budget -- never on timing. An index therefore
//! lands in the same chunk on every run, which is what lets a caller write one
//! slot per index and reduce those slots afterwards in index order for a
//! bit-identical result.
//!
//! Runs inline on the calling thread, creating and waking nothing, when the budget
//! is 1, when the range is empty, or when called from inside another `ParallelFor`
//! body. That last case is not a nicety: every worker is already busy running the
//! outer body, so a nested dispatch that waited on them would wait forever.
//!
//! An exception from `fn` propagates to the caller once all chunks have finished
//! -- never out of a worker's thread function, which would call std::terminate.
void ParallelFor(WorkerPool &pool, int64_t begin, int64_t end, const std::function<void(int64_t, int64_t)> &fn);

} // namespace miint
