#include "miint_parallel.hpp"

#include <stdexcept>
#include <string>
#include <system_error>

namespace miint {

namespace {

//! Nonzero while this thread is running a dispatched body.
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
	// could not save us. Stay serial: Threads() reads workers_.size() + 1, which is
	// 1. WASM is a first-class CORRECTNESS target for this code; it is deliberately
	// not a performance one.
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

void WorkerPool::DrainChunks(std::unique_lock<std::mutex> &lock) {
	while (next_chunk_ < n_chunks_) {
		const int c = next_chunk_++;
		lock.unlock();
		RunChunk(c);
		lock.lock();
		if (--outstanding_ == 0) {
			// Inert when the CALLER runs this -- it is the only thread that ever
			// waits on all_done_, and it is not waiting yet. Kept in the shared
			// helper rather than split into two nearly-identical loops.
			all_done_.notify_one();
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
		BodyScope scope;
		for (int c = 0; c < n_chunks; ++c) {
			chunk(c);
		}
		return;
	}

	dispatches_.fetch_add(1, std::memory_order_relaxed);
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
	// Wake only as many workers as there is work for. `notify_all` for a two-chunk
	// dispatch wakes every worker so that all but one can take the mutex, find
	// nothing to do and re-park -- contending with the caller on its way to
	// all_done_.wait, and so slowing down the dispatch it is not helping.
	const auto needed = static_cast<size_t>(n_chunks - 1);
	if (needed >= workers_.size()) {
		work_ready_.notify_all();
	} else {
		for (size_t i = 0; i < needed; ++i) {
			work_ready_.notify_one();
		}
	}

	// The caller works too, then helps with anything the workers have not claimed,
	// so no thread sits idle while chunks remain.
	RunChunk(0);
	std::unique_lock<std::mutex> lock(m_);
	--outstanding_;
	DrainChunks(lock);
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
		DrainChunks(lock);
	}
}

} // namespace miint
