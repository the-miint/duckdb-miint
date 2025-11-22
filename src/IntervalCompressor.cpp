#include "IntervalCompressor.hpp"

namespace miint {

constexpr size_t COMPRESS_THRESHOLD = 1'000'000;

IntervalCompressor::IntervalCompressor() = default;

void IntervalCompressor::Add(int64_t start, int64_t stop) {
	// Ensure start <= stop; swap if inverted
	if (start > stop) {
		std::swap(start, stop);
	}

	starts.push_back(start);
	stops.push_back(stop);

	if (starts.size() >= COMPRESS_THRESHOLD) {
		Compress();
	}
}

void IntervalCompressor::Compress() {
	if (starts.empty()) {
		return;
	}

	std::vector<std::pair<int64_t, int64_t>> intervals;
	intervals.reserve(starts.size());
	for (size_t i = 0; i < starts.size(); i++) {
		intervals.emplace_back(starts[i], stops[i]);
	}

	std::sort(intervals.begin(), intervals.end());

	starts.clear();
	stops.clear();

	int64_t current_start = intervals[0].first;
	int64_t current_stop = intervals[0].second;

	for (size_t i = 1; i < intervals.size(); i++) {
		if (intervals[i].first <= current_stop) {
			current_stop = std::max(current_stop, intervals[i].second);
		} else {
			starts.push_back(current_start);
			stops.push_back(current_stop);
			current_start = intervals[i].first;
			current_stop = intervals[i].second;
		}
	}

	starts.push_back(current_start);
	stops.push_back(current_stop);
}

bool IntervalCompressor::Empty() const {
	return starts.empty();
}

size_t IntervalCompressor::Size() const {
	return starts.size();
}

} // namespace miint
