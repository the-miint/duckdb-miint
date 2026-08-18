#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

namespace miint {

class IntervalCompressor {
public:
	std::vector<int64_t> starts;
	std::vector<int64_t> stops;

	IntervalCompressor();

	void Add(int64_t start, int64_t stop);
	void Compress();
	bool Empty() const;
	size_t Size() const;

	// How many times Compress() has run. Exposed so the amortization test can assert a
	// bound on total compression work deterministically instead of timing it.
	size_t CompressionCount() const {
		return compressions_;
	}

	// Add() compresses once the state reaches this many intervals. A FLOOR that grows --
	// see compress_floor_.
	static constexpr size_t COMPRESS_THRESHOLD = 1'000'000;

private:
	// Compression trigger. Must GROW after each compression, because Compress() merges
	// only OVERLAPPING intervals and reclaims nothing from a disjoint set. With a fixed
	// trigger, the first compression that fails to get below it makes every subsequent
	// Add() re-sort the entire state: measured ~10 ms PER ROW past 1e6, turning
	// O(n log n) into O(n^2 log n).
	//
	// Doubling off the post-compression size restores amortized O(n log n). Peak memory
	// becomes <= 2x the irreducible compressed state (or COMPRESS_THRESHOLD, whichever is
	// larger). No answer changes: Compress() computes the union of whatever is stored, and
	// replacing a subset by its union leaves the total union identical, so compressing more
	// or less often cannot change what a final Compress() produces.
	size_t compress_floor_ = COMPRESS_THRESHOLD;
	size_t compressions_ = 0;
};

} // namespace miint
