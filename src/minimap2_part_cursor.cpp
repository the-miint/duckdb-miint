#include "minimap2_part_cursor.hpp"
#include <cassert>
#include <exception>
#include <stdexcept>

namespace miint {

Minimap2PartCursor::Minimap2PartCursor(const std::string &index_path, const Minimap2Config &config,
                                       std::function<void()> flush_freed_memory)
    : reader_(std::make_unique<Minimap2IndexReader>(index_path, config)),
      flush_freed_memory_(std::move(flush_freed_memory)) {
	current_ = AdoptPart(reader_->ReadNextPart());
	if (!current_) {
		throw std::runtime_error("Index file '" + index_path + "' contains no parts");
	}
	// AtEof()'s probe can throw std::runtime_error too (fgetpos/fsetpos failure);
	// callers wrap everything from this constructor the same way.
	has_next_part_ = !reader_->AtEof();
	is_multi_part_ = has_next_part_;
}

std::shared_ptr<SharedMinimap2Index> Minimap2PartCursor::ReleaseSinglePart() {
	assert(!is_multi_part_);
	reader_.reset();
	return std::move(current_);
}

// The purge has to happen after the mm_idx_t is gone and on the thread that
// actually freed it, and neither is knowable from the outside: a part is read by
// whichever worker led the previous transition, and released by whichever one
// happens to drop the last reference — often a worker still inside align() when
// the leader has already moved on. Hanging the flush off the shared_ptr's
// deleter answers both questions for free, and replaces what would otherwise be
// a flush on every thread that might have been the last one (a process-wide
// purge per worker per transition, all but one of them purging nothing, most of
// them racing the leader's multi-GB read of the next part).
std::shared_ptr<SharedMinimap2Index> Minimap2PartCursor::AdoptPart(std::shared_ptr<SharedMinimap2Index> part) {
	if (!part || !flush_freed_memory_) {
		return part;
	}
	// Aliasing constructor: the new shared_ptr points at the same object but owns
	// `part` inside its deleter, so the index is destroyed (mm_idx_destroy) by the
	// reset below and the flush observes memory that is already free.
	auto flush = flush_freed_memory_;
	auto *raw = part.get();
	return std::shared_ptr<SharedMinimap2Index>(raw, [part, flush](SharedMinimap2Index *) mutable {
		part.reset();
		flush();
	});
}

void Minimap2PartCursor::EnsureAttached(Attachment &att, Minimap2Aligner &aligner) {
	if (att.owner == this && att.attached && att.generation == generation_) {
		return;
	}
	// current_ is null only while a leader is between resetting it and
	// installing the next part. A thread arriving in that window (a fresh
	// worker, or one whose attachment is stale from an earlier cursor) must not
	// attach a null index; recording the current generation without attaching
	// makes its next Advance() wait on this very transition rather than spin.
	att.owner = this;
	att.attached = false;
	att.generation = generation_;
	if (!current_) {
		return;
	}
	aligner.attach_shared_index(current_);
	att.attached = true;
}

bool Minimap2PartCursor::CurrentIsLastPart() const {
	return parts_exhausted_ || !has_next_part_;
}

void Minimap2PartCursor::MarkExhausted() {
	std::lock_guard<std::mutex> lock(lock_);
	parts_exhausted_ = true;
}

bool Minimap2PartCursor::Advance(Attachment &att, Minimap2Aligner &aligner, const std::function<void()> &prepare,
                                 const std::function<void()> &publish) {
	const uint64_t expected_generation = att.generation;
	// Detach FIRST, before either leading or waiting: a thread that waits while
	// still attached pins the outgoing part for the whole load, and an idle one
	// pins it for the rest of the query. See the class comment.
	aligner.detach_shared_index();
	att.attached = false;

	std::unique_lock<std::mutex> lock(lock_);

	while (true) {
		if (generation_ != expected_generation) {
			// Someone already advanced past the generation this thread was stuck on.
			return true;
		}
		if (parts_exhausted_) {
			return false;
		}
		if (advancing_) {
			cv_.wait(lock);
			continue;
		}

		// Become the leader for this transition. Clearing current_ under the lock
		// stops any thread arriving from here on from attaching to the outgoing
		// part; releasing it only after unlocking keeps the part's destructor —
		// and the process-wide purge its deleter runs — off the critical section
		// every other worker is about to block on.
		advancing_ = true;
		auto outgoing = std::move(current_);
		lock.unlock();
		outgoing.reset(); // free the just-finished part before loading the next

		std::shared_ptr<SharedMinimap2Index> next_index;
		bool next_has_successor = false;
		std::exception_ptr load_error;
		try {
			next_index = AdoptPart(reader_->ReadNextPart());
			if (next_index) {
				// Probe for the part after this one while this thread still has
				// the reader to itself; CurrentIsLastPart() then answers from the
				// flag without touching the file.
				next_has_successor = !reader_->AtEof();
				if (prepare) {
					prepare();
				}
			}
		} catch (...) {
			load_error = std::current_exception();
		}

		lock.lock();
		advancing_ = false;
		if (load_error) {
			// Stop every other thread too — the reader is now in an unknown
			// state and cannot be trusted for a retry.
			parts_exhausted_ = true;
			cv_.notify_all();
			lock.unlock();
			std::rethrow_exception(load_error);
		}
		if (!next_index) {
			parts_exhausted_ = true;
			cv_.notify_all();
			return false;
		}

		current_ = std::move(next_index);
		has_next_part_ = next_has_successor;
		if (publish) {
			publish();
		}
		generation_++;
		cv_.notify_all();
		return true;
	}
}

} // namespace miint
