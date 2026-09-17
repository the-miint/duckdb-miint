#pragma once
#include "Minimap2Aligner.hpp"
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <string>

namespace miint {

// Owns a prebuilt .mmi (single- or multi-part) and walks N worker threads
// through its parts one at a time, keeping exactly one part resident.
//
// minimap2's own answer to "reference larger than RAM" is the multi-part index
// (`minimap2 -d out.mmi -I <batch_size> ref.fa`): every query is aligned against
// part 1, then part 2, and so on, so peak memory is one part rather than the
// whole index. This class is the thread-coordination half of that loop, shared
// by align_minimap2 (one cursor for the whole query) and align_minimap2_sharded
// (one cursor per shard). Minimap2IndexReader is the file half.
//
// Ownership of a part: the cursor's own `current_` reference plus one
// shared_ptr copy inside every Minimap2Aligner attached to it. A part is freed
// (mm_idx_destroy) when the LAST of those drops, so a transition has to drop
// both — Advance() detaches the calling thread's aligner first, and the leader
// resets `current_` before loading the next part. The only remaining overlap
// is a worker still inside align() on the outgoing part, which keeps its own
// copy alive until that one sub-batch finishes: the transient peak is
// (outgoing part still in flight + incoming part), and the steady-state peak
// tracks whichever part is biggest. Measured at 5.28GB peak RSS for an
// evenly-split 4-part ~9GB human genome index.
//
// Attachment is LAZY, on first real use via EnsureAttached, never eagerly at
// thread init: DuckDB may initialize more worker states than ever receive real
// work (observed with a handful of query rows and 12 threads — most threads got
// exactly one empty fetch and never ran again), and an eager attach would leave
// each such idle thread holding a live reference to whatever part was current
// at init time for the rest of the query, defeating the one-part-at-a-time
// bound. Measured against a 2-part ~9GB index with 12 threads and a tiny query
// table: peak RSS was 15.5GB — HIGHER than loading the same index single-part
// (11.75GB) — because 10-11 idle threads attached to part 1 and sat holding it
// while part 2 loaded fully.
//
// Protocol (caller side):
//   { std::lock_guard<std::mutex> g(cursor.Lock());
//     cursor.EnsureAttached(lstate.part, *lstate.aligner);
//     <claim this thread's next unit of work against the current part> }
//   ... align it ...
//   when the current part has no work left for this thread:
//     if (!cursor.Advance(lstate.part, *lstate.aligner, prepare, publish)) done;
//     else loop back (a newer part exists — re-attach and retry)
// Claiming work inside the same critical section as EnsureAttached is what
// makes the per-part reset in `publish` race-free: a thread can never take
// work belonging to part k+1 while its aligner is still attached to part k.
class Minimap2PartCursor {
public:
	// Per-thread record of which part of WHICH cursor that thread's aligner is
	// attached to. Lives in the caller's per-thread state.
	//
	// `owner` is what makes the token safe to carry between cursors: a sharded
	// worker moving from shard A to shard B arrives holding A's attachment, and
	// B's fresh cursor is also at generation 0, so comparing generations alone
	// would treat the thread as already attached and leave its aligner with no
	// index at all. Stamping the cursor here means a stale token is recognised
	// rather than relied on being cleared by the caller.
	struct Attachment {
		const Minimap2PartCursor *owner = nullptr;
		bool attached = false;
		uint64_t generation = 0;
	};

	// Opens the index and loads part 1, then peeks (4 bytes, see
	// Minimap2IndexReader::AtEof) for whether a part 2 exists. Throws
	// std::runtime_error on open/load failure or an index with no parts.
	//
	// `flush_freed_memory` is invoked on the calling thread right after every
	// point where this cursor may have just dropped the last reference to a part.
	// DuckDB does not return a busy worker's freed memory to the OS on its own —
	// see MakeFreedMemoryFlusher in align_common.hpp. May be empty.
	Minimap2PartCursor(const std::string &index_path, const Minimap2Config &config,
	                   std::function<void()> flush_freed_memory);

	Minimap2PartCursor(const Minimap2PartCursor &) = delete;
	Minimap2PartCursor &operator=(const Minimap2PartCursor &) = delete;

	// True if the file holds more than one part. Decided at construction.
	bool IsMultiPart() const {
		return is_multi_part_;
	}

	// Single-part fast path: hands part 1 to the caller and drops the reader, so
	// a caller that has its own single-index machinery (align_minimap2's
	// SharedMinimap2Index path) can use it unchanged. Only valid when
	// !IsMultiPart(); the cursor must not be used afterwards.
	std::shared_ptr<SharedMinimap2Index> ReleaseSinglePart();

	// Attaches `aligner` to the current part if needed, then runs `claim` with
	// the cursor still locked and returns whatever it returns.
	//
	// This is the ONLY way to claim a unit of work, and it is one call rather
	// than an exposed mutex because the two steps must not come apart: the work
	// source is reset every time a new part is published, so a claim made
	// outside this critical section could hand a worker still attached to part k
	// a unit belonging to part k+1 — aligned against the wrong part, and never
	// against the right one. `claim` runs under the lock, so it must not block
	// or run a DuckDB query; it is meant to be a counter bump or a pointer read.
	//
	// Attaching covers a fresh thread's lazy first attach and a thread returning
	// from Advance (which always detaches first). While a leader is
	// mid-transition no part is resident, so nothing is attached and `claim`
	// simply finds the outgoing part's source exhausted, sending the caller into
	// Advance() to wait for the leader.
	template <class Fn>
	auto WithCurrentPart(Attachment &att, Minimap2Aligner &aligner, Fn &&claim) -> decltype(claim()) {
		std::lock_guard<std::mutex> lock(lock_);
		EnsureAttached(att, aligner);
		return claim();
	}

	// Only valid from inside a WithCurrentPart `claim` (it reads state the lock
	// guards). True when no part follows the current one. Reads a flag refreshed
	// once per load, never the file, so a caller may ask per unit of work. Lets
	// a caller stop admitting new workers once the very last unit has been
	// claimed, without ever mistaking "last batch of part k" for "last batch of
	// the index". False while a leader is mid-transition (a next part exists —
	// it is being loaded), which is the conservative answer there.
	bool CurrentIsLastPart() const;

	// The calling thread has no work left against the part `att` is on. Detaches
	// its aligner, then either leads the transition to the next part or waits
	// for the thread already doing so.
	//
	// Exactly one thread performs the load per transition (the `advancing_` flag
	// elects a leader); every other thread blocks on the condition variable. Both
	// the index load and `prepare` run OUTSIDE the lock — every other thread is
	// parked or about to block on the lock, so no scheduler thread for this
	// pipeline is free to service a query issued while holding it — and inside
	// one try/catch, so a throw from either still notifies waiters instead of
	// stranding them.
	//
	//   prepare — leader only, outside the lock, after the next part loaded and
	//             before it is published. May run DuckDB queries (e.g. open a
	//             replay stream over the query snapshot). May throw: the query
	//             fails and every waiter is released.
	//   publish — leader only, UNDER the lock, immediately before the generation
	//             bump. Swaps the caller's per-part work source into place
	//             (install the stream opened in `prepare`, reset a batch offset
	//             to 0). Must not throw.
	//
	// Returns true once a part newer than `att.generation` exists to retry
	// against (whether this thread loaded it or another did), false once the
	// reader is exhausted with nothing left to try.
	bool Advance(Attachment &att, Minimap2Aligner &aligner, const std::function<void()> &prepare,
	             const std::function<void()> &publish);

	// Marks the cursor exhausted without touching the reader: a caller with zero
	// units of work (empty query relation) uses this so no later part is ever
	// loaded and decoded just to align nothing against it. The next Advance()
	// returns false.
	void MarkExhausted();

private:
	void FlushFreedMemory();
	void EnsureAttached(Attachment &att, Minimap2Aligner &aligner);

	std::unique_ptr<Minimap2IndexReader> reader_;
	std::shared_ptr<SharedMinimap2Index> current_;
	std::function<void()> flush_freed_memory_;
	bool is_multi_part_ = false;
	// Whether a part follows current_. Refreshed by the leader right after each
	// load, while it alone owns the reader, so CurrentIsLastPart() never has to
	// touch the file (it is asked once per claimed batch, not once per part).
	bool has_next_part_ = false;

	std::mutex lock_;
	std::condition_variable cv_;
	uint64_t generation_ = 0;
	bool advancing_ = false;       // true while one thread is loading the next part
	bool parts_exhausted_ = false; // reader returned no next part (or MarkExhausted)
};

} // namespace miint
