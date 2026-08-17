#pragma once

#include <cstdint>
#include <vector>

namespace miint {

struct CumulativePoint {
	int32_t rank;
	int64_t covered;
};

// Cumulative breadth of coverage in rank order: for each rank k, how many bases are
// covered by the UNION of every interval with rank <= k (issue #214).
//
// micov recompresses the entire accumulated interval set at every rank, which is
// O(n^2) in intervals and is the dominant cost of a micov run. This computes every
// rank in a single O(n log n) sweep using the min-covering-rank identity:
//
//   the union of ranks 0..k  ==  the bases whose LOWEST covering rank is <= k
//
// so histogramming bases by their minimum covering rank makes the curve that
// histogram's prefix sum. One sort, no per-rank recompression.
//
// Deliberately NOT built on IntervalCompressor: its Compress() re-sorts the whole
// accumulated set on every call, which is exactly the O(n^2) this exists to remove.
//
// Intervals are 1-based half-open, matching read_alignments and compress_intervals.
// Ranks are 0-based and must be contiguous -- see Curve().
class CumulativeCoverageAccumulator {
public:
	// Records an interval against a rank. An interval with start >= stop covers no
	// base and is dropped, but the rank is still registered. Note this DIVERGES from
	// IntervalCompressor, which silently swaps an inverted pair: swapping manufactures
	// coverage out of transposed columns, so the newer convention used by
	// region_presence and region_coverage -- drop it -- is followed here instead.
	void Add(int32_t rank, int64_t start, int64_t stop);

	// Registers a rank that has no coverage at all. Load-bearing: a sample in the
	// group with no coverage of the target contributes no interval rows, so without
	// this it would vanish from the curve, understating the group size and making
	// curves from groups with different detection rates incomparable.
	void AddEmptyRank(int32_t rank);

	// Merges another state into this one -- the DuckDB Combine path under a parallel
	// GROUP BY. Order-independent, like the sweep itself.
	void Absorb(const CumulativeCoverageAccumulator &other);

	bool Empty() const;

	// The curve, one point per rank, rank ascending.
	//
	// Throws miint::InvalidInputException unless the observed ranks are exactly
	// 0..R-1 with no gaps. A gap would silently shift the x-axis -- rank 1 reporting
	// what is really rank 2's coverage -- and ROW_NUMBER() OVER (...) - 1, the
	// intended way to produce ranks, always satisfies the constraint.
	std::vector<CumulativePoint> Curve() const;

	// Replaces each rank's own intervals with their union, and collapses a rank that
	// retains no interval to a single zero-coverage marker. Bounds the state without
	// changing any answer: merging a rank's intervals with each other cannot change any
	// base's MINIMUM covering rank, since that base is still covered by the same rank.
	//
	// Called automatically from Add() past a threshold and from Absorb(), mirroring
	// IntervalCompressor. It is public only so the unit tests can assert idempotence
	// and invariance directly.
	void Compact();

	// Number of stored observations. Exposed for the compaction tests -- Absorb() reads
	// the private vectors directly, being a member.
	size_t ObservationCount() const {
		return ranks_.size();
	}

	// How many times Compact() has run on this state. Exposed so the amortization test
	// can assert a BOUND on total compaction work deterministically, rather than timing
	// it: a fixed threshold degrades to one compaction per Add() on input that cannot be
	// merged, and a wall-clock assertion for that would be flaky.
	size_t CompactionCount() const {
		return compactions_;
	}

	// Compact() runs once the state reaches this many observations. It is a FLOOR that
	// grows: see the comment on compact_floor_.
	static constexpr size_t COMPACT_THRESHOLD = 1'000'000;

private:
	// Raw (rank, start, stop) observations, swept once at Curve() time. Compact()
	// shrinks them per rank; it deliberately does NOT merge across ranks, because the
	// curve needs to know which rank first covered each base.
	std::vector<int32_t> ranks_;
	std::vector<int64_t> starts_;
	std::vector<int64_t> stops_;

	// Compaction trigger. Must GROW after each compaction, because Compact() only merges
	// intervals WITHIN a rank and can therefore shrink the state by nothing at all -- and
	// the documented input shape is exactly the shape it cannot shrink, since positions
	// are normally compress_intervals output and already disjoint. With a fixed trigger,
	// the first compaction that fails to get below it makes every subsequent Add() re-sort
	// the entire state, turning O(n log n) into O(n^2 log n). Measured with a fixed floor at
	// COMPACT_THRESHOLD + 100 disjoint observations: 101 compactions of a ~1e6-observation
	// state, ~17.5 ms each. Extrapolating that per-row cost, the 1.2M rows a 100-sample
	// group of 12k intervals produces would take roughly an hour.
	//
	// Doubling the floor off the post-compaction size restores amortized O(n log n) -- the
	// state can at most double between compactions, so compactions are O(log n) and each
	// costs at most twice the irreducible size. Peak memory becomes <= 2x the irreducible
	// compacted state (or COMPACT_THRESHOLD, whichever is larger) rather than <=
	// COMPACT_THRESHOLD; that 2x is what buys the asymptotic fix.
	size_t compact_floor_ = COMPACT_THRESHOLD;
	size_t compactions_ = 0;
};

} // namespace miint
