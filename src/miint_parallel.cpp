#include "miint_parallel.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <system_error>

namespace miint {

namespace {

//! Nonzero while this thread is running a ParallelFor body.
//!
//! A nested dispatch has to run inline, because by construction every worker is
//! already inside the outer body -- waiting for one to pick up the inner work
//! would wait forever. Running inline is a defined answer; deadlocking is not.
thread_local int g_in_parallel_body = 0;

struct BodyScope {
	BodyScope() {
		++g_in_parallel_body;
	}
	~BodyScope() {
		--g_in_parallel_body;
	}
	BodyScope(const BodyScope &) = delete;
	BodyScope &operator=(const BodyScope &) = delete;
};

} // namespace

WorkerPool::WorkerPool(int n_threads) {
	if (n_threads < 1) {
		throw std::invalid_argument("WorkerPool: n_threads must be >= 1 (got " + std::to_string(n_threads) + ")");
	}
#if defined(__EMSCRIPTEN__) && !defined(__EMSCRIPTEN_PTHREADS__)
	// A wasm build without -pthread has no threads to create at all, and there
	// std::thread's constructor aborts rather than throwing -- so the catch below
	// could not save us. Stay serial. WASM is a first-class CORRECTNESS target for
	// this code; it is deliberately not a performance one.
	(void)n_threads;
	return;
#else
	// One fewer than the budget: the calling thread is the last worker.
	workers_.reserve(static_cast<size_t>(n_threads - 1));
	try {
		for (int t = 1; t < n_threads; ++t) {
			workers_.emplace_back([this] { WorkerLoop(); });
		}
	} catch (const std::system_error &) {
		// The thread constructor failed, e.g. EAGAIN under a pids/ulimit cap. Do
		// NOT let it unwind past the joinable threads already in `workers_` -- a
		// joinable std::thread's destructor calls std::terminate() and would crash
		// the process. Degrade instead: whatever threads we did get are perfectly
		// usable, and Threads() reports the smaller budget. Same defence, and the
		// same reasoning, as community_distances.cpp.
		//
		// `workers_` is pre-reserved, so emplace_back itself cannot throw
		// bad_alloc; the thread constructor is the only thing here that throws.
	}
	n_threads_ = static_cast<int>(workers_.size()) + 1;
#endif
}

WorkerPool::~WorkerPool() {
	{
		std::lock_guard<std::mutex> lock(m_);
		stop_ = true;
	}
	work_ready_.notify_all();
	for (std::thread &worker : workers_) {
		worker.join();
	}
}

void WorkerPool::RunChunk(int c) {
	BodyScope scope;
	try {
		(*chunk_)(c);
	} catch (...) {
		std::lock_guard<std::mutex> lock(m_);
		if (!failure_ || c < failed_chunk_) {
			failure_ = std::current_exception();
			failed_chunk_ = c;
		}
	}
}

void WorkerPool::Dispatch(int n_chunks, const std::function<void(int)> &chunk) {
	if (n_chunks <= 0) {
		return;
	}
	if (workers_.empty() || n_chunks == 1 || g_in_parallel_body > 0) {
		// Nothing to hand out, or we are already inside a body. Exceptions
		// propagate straight to the caller here, which is the same observable
		// behaviour as the threaded path below.
		for (int c = 0; c < n_chunks; ++c) {
			BodyScope scope;
			chunk(c);
		}
		return;
	}

	{
		std::lock_guard<std::mutex> lock(m_);
		chunk_ = &chunk;
		n_chunks_ = n_chunks;
		next_chunk_ = 1; // chunk 0 is the caller's
		outstanding_ = n_chunks;
		failure_ = nullptr;
		failed_chunk_ = 0;
		++generation_;
	}
	work_ready_.notify_all();

	// The caller works too, then helps with anything the workers have not claimed,
	// so no thread sits idle while chunks remain.
	RunChunk(0);
	std::unique_lock<std::mutex> lock(m_);
	--outstanding_;
	while (next_chunk_ < n_chunks_) {
		const int c = next_chunk_++;
		lock.unlock();
		RunChunk(c);
		lock.lock();
		--outstanding_;
	}
	all_done_.wait(lock, [this] { return outstanding_ == 0; });

	if (failure_) {
		const std::exception_ptr failure = failure_;
		failure_ = nullptr;
		lock.unlock();
		std::rethrow_exception(failure);
	}
}

void WorkerPool::WorkerLoop() {
	uint64_t seen = 0;
	std::unique_lock<std::mutex> lock(m_);
	for (;;) {
		// Compared against a per-worker `seen` rather than a flag, so a worker that
		// was slow to wake cannot miss a dispatch and cannot re-run an old one.
		work_ready_.wait(lock, [this, &seen] { return stop_ || generation_ != seen; });
		if (stop_) {
			return;
		}
		seen = generation_;
		while (next_chunk_ < n_chunks_) {
			const int c = next_chunk_++;
			lock.unlock();
			RunChunk(c);
			lock.lock();
			if (--outstanding_ == 0) {
				all_done_.notify_one();
			}
		}
	}
}

void ParallelFor(WorkerPool &pool, int64_t begin, int64_t end, const std::function<void(int64_t, int64_t)> &fn) {
	const int64_t n = end - begin;
	if (n <= 0) {
		return;
	}
	// Never more chunks than indices, so every chunk handed out is non-empty.
	const int64_t chunks = std::min<int64_t>(pool.Threads(), n);
	const int64_t base = n / chunks;
	const int64_t remainder = n % chunks;

	pool.Dispatch(static_cast<int>(chunks), [&](int c) {
		// Spread the remainder over the first `remainder` chunks. Written as
		// products of the chunk index rather than `n * c / chunks` so it cannot
		// overflow for a large range, and so the boundaries are exact integers
		// rather than the result of two roundings.
		const int64_t index = c;
		const int64_t lo = begin + index * base + std::min<int64_t>(index, remainder);
		const int64_t hi = lo + base + (index < remainder ? 1 : 0);
		fn(lo, hi);
	});
}

} // namespace miint
