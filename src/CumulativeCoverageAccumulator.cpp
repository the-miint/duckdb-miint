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

// Matches IntervalCompressor::COMPRESS_THRESHOLD. Past this many stored observations
// Add() folds each rank's intervals into their union, so a deep metagenomic group cannot
// grow the state without bound between combines.
constexpr size_t COMPACT_THRESHOLD = 1'000'000;

void CumulativeCoverageAccumulator::Add(int32_t rank, int64_t start, int64_t stop) {
	ranks_.push_back(rank);
	starts_.push_back(start);
	stops_.push_back(stop);

	if (ranks_.size() >= COMPACT_THRESHOLD) {
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

	// Compact after combining rather than before, mirroring compress_intervals: under a
	// parallel GROUP BY the combines are where the state would otherwise pile up.
	Compact();
}

void CumulativeCoverageAccumulator::Compact() {
	if (ranks_.empty()) {
		return;
	}

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
}

bool CumulativeCoverageAccumulator::Empty() const {
	return ranks_.empty();
}

std::vector<CumulativePoint> CumulativeCoverageAccumulator::Curve() const {
	if (ranks_.empty()) {
		return {};
	}

	// Validation walks the ranks in sorted order rather than marking a presence array.
	// Sorting first is what makes every branch below reachable AND bounds the histogram
	// allocation: a rank of 1e9 is rejected as a gap before anything is sized off it,
	// so there is no need for a separate "is max_rank plausible" guard. An earlier
	// version had one, and it fired first on every malformed input -- which left the
	// gap and start-at-zero branches dead code that no test could reach.
	std::vector<int32_t> sorted_ranks = ranks_;
	std::sort(sorted_ranks.begin(), sorted_ranks.end());

	if (sorted_ranks.front() != 0) {
		throw InvalidInputException(
		    "cumulative_coverage: ranks must start at 0 (lowest rank seen is " + std::to_string(sorted_ranks.front()) +
		    "). Produce ranks with ROW_NUMBER() OVER (...) - 1 -- note the - 1, ROW_NUMBER() itself is 1-based.");
	}

	int32_t prev_rank = sorted_ranks.front();
	for (size_t k = 1; k < sorted_ranks.size(); k++) {
		const int32_t r = sorted_ranks[k];
		if (r != prev_rank && r != prev_rank + 1) {
			throw InvalidInputException("cumulative_coverage: ranks must be contiguous 0..N-1 (rank " +
			                            std::to_string(prev_rank + 1) + " is missing; ranks jump from " +
			                            std::to_string(prev_rank) + " to " + std::to_string(r) +
			                            "). Produce ranks with ROW_NUMBER() OVER (...) - 1.");
		}
		prev_rank = r;
	}

	// Contiguity is now established, so n_ranks <= the number of rows observed.
	const size_t n_ranks = static_cast<size_t>(prev_rank) + 1;

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
