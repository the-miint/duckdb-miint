#include "CumulativeCoverageAccumulator.hpp"
#include "alignment_functions_internal.hpp"

#include <algorithm>
#include <functional>
#include <queue>
#include <tuple>
#include <string>

namespace miint {

namespace {

struct SweepEvent {
	int64_t pos;
	int32_t rank;
	int32_t delta; // +1 opens an interval, -1 closes it
};

} // namespace

void CumulativeCoverageAccumulator::Add(int32_t rank, int64_t start, int64_t stop) {
	ranks_.push_back(rank);
	starts_.push_back(start);
	stops_.push_back(stop);

	if (ranks_.size() >= compact_floor_) {
		Compact();
	}
}

void CumulativeCoverageAccumulator::AddEmptyRank(int32_t rank) {
	// Stored as a degenerate interval: it registers the rank and covers nothing, which
	// is exactly the semantics wanted, with no second code path to keep in sync.
	Add(rank, 0, 0);
}

void CumulativeCoverageAccumulator::Absorb(const CumulativeCoverageAccumulator &other) {
	// Self-absorb would read the source range while insert() reallocates it -- undefined
	// behaviour. DuckDB never self-combines, so this is a latent case rather than a live
	// one, but Absorb() is public and documented as the Combine path.
	if (&other == this) {
		return;
	}

	ranks_.insert(ranks_.end(), other.ranks_.begin(), other.ranks_.end());
	starts_.insert(starts_.end(), other.starts_.begin(), other.starts_.end());
	stops_.insert(stops_.end(), other.stops_.begin(), other.stops_.end());

	// Same threshold Add() uses, NOT an unconditional compaction. compress_intervals
	// compacts on every Combine, but that is measurably the wrong policy here: it costs
	// O(m log m) per combine to reclaim nothing when the state is already small, so the
	// thread count made things WORSE. Measured on 500k rows / 50 ranks of disjoint
	// intervals (the shape Compact() cannot shrink), compacting unconditionally:
	// 0.019 s at threads=1 -> 0.082 s at threads=16, i.e. 4.3x slower for adding cores.
	// Guarded by the threshold, threads=16 is 0.024 s. The threshold alone bounds the state.
	if (ranks_.size() >= compact_floor_) {
		Compact();
	}
}

void CumulativeCoverageAccumulator::Compact() {
	if (ranks_.empty()) {
		return;
	}
	compactions_++;

	// (rank, start, stop) sorted by rank then start, so each rank's intervals are
	// contiguous and ascending and can be merged in one pass.
	std::vector<std::tuple<int32_t, int64_t, int64_t>> obs;
	obs.reserve(ranks_.size());
	for (size_t i = 0; i < ranks_.size(); i++) {
		obs.emplace_back(ranks_[i], starts_[i], stops_[i]);
	}
	std::sort(obs.begin(), obs.end());

	std::vector<int32_t> out_ranks;
	std::vector<int64_t> out_starts;
	std::vector<int64_t> out_stops;

	size_t i = 0;
	while (i < obs.size()) {
		const int32_t rank = std::get<0>(obs[i]);
		const size_t rank_begin = out_ranks.size();

		// Merge this rank's usable intervals. Degenerate ones contribute nothing and are
		// skipped here; the rank still survives via the marker below.
		bool open = false;
		int64_t cur_start = 0;
		int64_t cur_stop = 0;
		while (i < obs.size() && std::get<0>(obs[i]) == rank) {
			const int64_t s = std::get<1>(obs[i]);
			const int64_t e = std::get<2>(obs[i]);
			i++;
			if (s >= e) {
				continue;
			}
			if (!open) {
				open = true;
				cur_start = s;
				cur_stop = e;
			} else if (s <= cur_stop) {
				cur_stop = std::max(cur_stop, e);
			} else {
				out_ranks.push_back(rank);
				out_starts.push_back(cur_start);
				out_stops.push_back(cur_stop);
				cur_start = s;
				cur_stop = e;
			}
		}
		if (open) {
			out_ranks.push_back(rank);
			out_starts.push_back(cur_start);
			out_stops.push_back(cur_stop);
		}

		// A rank whose every interval was degenerate must NOT disappear -- it is a
		// zero-coverage sample, and losing it would shorten the curve and shift every
		// rank above it.
		if (out_ranks.size() == rank_begin) {
			out_ranks.push_back(rank);
			out_starts.push_back(0);
			out_stops.push_back(0);
		}
	}

	ranks_.swap(out_ranks);
	starts_.swap(out_starts);
	stops_.swap(out_stops);

	// Raise the floor above what we just failed to reclaim. Compact() merges only WITHIN a
	// rank, so on already-disjoint input -- the documented shape -- it reclaims nothing and
	// a fixed floor would re-trigger on the very next Add(). See compact_floor_'s comment.
	compact_floor_ = std::max(COMPACT_THRESHOLD, ranks_.size() * 2);
}

bool CumulativeCoverageAccumulator::Empty() const {
	return ranks_.empty();
}

std::vector<CumulativePoint> CumulativeCoverageAccumulator::Curve() const {
	if (ranks_.empty()) {
		return {};
	}

	// Validation is two linear passes plus a bounded presence array, NOT an O(n log n) sort
	// of a copy of every rank. The old form allocated a second int32 per observation and
	// sorted it on input the feeding hash join has already shuffled, purely to answer "are
	// these exactly 0..R-1"; that question does not need an ordering.
	//
	// Every branch stays reachable and the histogram allocation stays bounded WITHOUT a
	// separate "is max_rank plausible" guard -- an earlier version had one and it fired
	// first on every malformed input, leaving the gap and start-at-zero branches dead code
	// no test could reach. What replaces it is a pigeonhole argument: contiguous ranks
	// 0..R-1 need R distinct values, so a valid R is at most the number of observations.
	// Any rank above that count therefore GUARANTEES a gap at or below it, which is why a
	// rank of 1e9 is reported as a gap instead of sizing a billion buckets.
	int32_t min_rank = ranks_[0];
	int32_t max_rank = ranks_[0];
	for (const int32_t r : ranks_) {
		min_rank = std::min(min_rank, r);
		max_rank = std::max(max_rank, r);
	}

	if (min_rank != 0) {
		throw InvalidInputException(
		    "cumulative_coverage: ranks must start at 0 (lowest rank seen is " + std::to_string(min_rank) +
		    "). Produce ranks with ROW_NUMBER() OVER (...) - 1 -- note the - 1, ROW_NUMBER() itself is 1-based.");
	}

	// Every rank is >= 0 from here. Probing only the first min(max_rank, n) + 1 slots is
	// what bounds the allocation: past that the pigeonhole above already proves a gap lies
	// below, so there is nothing a bigger array could discover.
	const size_t probe = std::min(static_cast<size_t>(max_rank), ranks_.size());
	std::vector<char> present(probe + 1, 0);
	for (const int32_t r : ranks_) {
		if (static_cast<size_t>(r) <= probe) {
			present[static_cast<size_t>(r)] = 1;
		}
	}

	for (size_t r = 0; r <= probe; r++) {
		if (present[r] != 0) {
			continue;
		}
		// r is missing. r >= 1 here, because min_rank == 0 was just established and so
		// present[0] is set. Report the lowest rank actually observed above r, so the
		// message reads exactly as a sorted walk would have produced it.
		int32_t next_seen = max_rank;
		for (const int32_t seen : ranks_) {
			if (static_cast<size_t>(seen) > r) {
				next_seen = std::min(next_seen, seen);
			}
		}
		throw InvalidInputException("cumulative_coverage: ranks must be contiguous 0..N-1 (rank " + std::to_string(r) +
		                            " is missing; ranks jump from " + std::to_string(static_cast<int64_t>(r) - 1) +
		                            " to " + std::to_string(next_seen) +
		                            "). Produce ranks with ROW_NUMBER() OVER (...) - 1.");
	}

	// Contiguity is now established, so n_ranks <= the number of rows observed.
	const size_t n_ranks = static_cast<size_t>(max_rank) + 1;

	std::vector<SweepEvent> events;
	events.reserve(ranks_.size() * 2);
	for (size_t i = 0; i < ranks_.size(); i++) {
		if (starts_[i] >= stops_[i]) {
			// Covers no base. The rank is still counted, because n_ranks comes from the
			// sorted-ranks validation above, which reads every observation before this
			// filter runs.
			continue;
		}
		events.push_back({starts_[i], ranks_[i], 1});
		events.push_back({stops_[i], ranks_[i], -1});
	}

	std::sort(events.begin(), events.end(), [](const SweepEvent &a, const SweepEvent &b) { return a.pos < b.pos; });

	// hist[r] = number of bases whose minimum covering rank is exactly r.
	std::vector<int64_t> hist(n_ranks, 0);

	// The active set, as a count per rank plus a lazily-cleaned min-heap. A multiset
	// would be shorter but allocates a node per event; this is vector-backed, and each
	// rank is pushed at most once per interval and popped at most once overall, so the
	// sweep stays O(n log n) with a handful of allocations.
	std::vector<int32_t> active_count(n_ranks, 0);
	std::priority_queue<int32_t, std::vector<int32_t>, std::greater<int32_t>> min_active;
	int64_t n_active = 0;

	size_t i = 0;
	int64_t prev_pos = events.empty() ? 0 : events[0].pos;
	while (i < events.size()) {
		const int64_t pos = events[i].pos;

		// The segment [prev_pos, pos) is covered by the active set as it stood before
		// the events at `pos` are applied. Credit its whole width to the lowest rank
		// covering it -- that is the min-covering-rank identity, applied one
		// elementary segment at a time.
		if (n_active > 0 && pos > prev_pos) {
			while (!min_active.empty() && active_count[static_cast<size_t>(min_active.top())] == 0) {
				min_active.pop();
			}
			hist[static_cast<size_t>(min_active.top())] += pos - prev_pos;
		}

		while (i < events.size() && events[i].pos == pos) {
			const auto rank_idx = static_cast<size_t>(events[i].rank);
			if (events[i].delta > 0) {
				active_count[rank_idx]++;
				min_active.push(events[i].rank);
				n_active++;
			} else {
				active_count[rank_idx]--;
				n_active--;
			}
			i++;
		}

		prev_pos = pos;
	}

	std::vector<CumulativePoint> curve;
	curve.reserve(n_ranks);
	int64_t running = 0;
	for (size_t r = 0; r < n_ranks; r++) {
		running += hist[r];
		curve.push_back({static_cast<int32_t>(r), running});
	}
	return curve;
}

} // namespace miint
