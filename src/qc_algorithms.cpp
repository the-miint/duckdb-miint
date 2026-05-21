// Portions of this file are derived from fastp
// (https://github.com/OpenGene/fastp), MIT-licensed,
// Copyright (c) 2016 OpenGene. See THIRD_PARTY_LICENSES.md.

#include "qc_algorithms.hpp"

#include <stdexcept>

namespace miint::qc {

namespace {

// Sum of the initial window [0, w). Caller guarantees len >= w.
inline std::uint32_t initial_window_sum(const std::uint8_t *q, std::size_t w) {
	std::uint32_t sum = 0;
	for (std::size_t i = 0; i < w; i++) {
		sum += q[i];
	}
	return sum;
}

// Sum of the final window [len-w, len). Caller guarantees len >= w.
inline std::uint32_t trailing_window_sum(const std::uint8_t *q, std::size_t len, std::size_t w) {
	std::uint32_t sum = 0;
	for (std::size_t i = len - w; i < len; i++) {
		sum += q[i];
	}
	return sum;
}

inline void require_window(std::size_t w) {
	if (w == 0) {
		throw std::invalid_argument("SlidingWindowTrimmer: window_size must be > 0");
	}
}

// Integer-safe mean-vs-threshold test: window_mean >= threshold  iff  sum >= threshold * w.
inline std::uint32_t threshold_sum(std::uint8_t mean_quality, std::size_t w) {
	return static_cast<std::uint32_t>(mean_quality) * static_cast<std::uint32_t>(w);
}

} // namespace

TrimResult SlidingWindowTrimmer::trim_3p(const std::uint8_t *q, std::size_t n, std::size_t w,
                                         std::uint8_t mean_quality) {
	require_window(w);
	if (n < w) {
		return {0, n};
	}

	std::uint32_t sum = trailing_window_sum(q, n, w);
	std::size_t end = n; // current window is [end-w, end)
	const std::uint32_t threshold = threshold_sum(mean_quality, w);

	while (sum < threshold) {
		if (end - w == 0) {
			return {0, 0};
		}
		// Slide back by 1: window becomes [end-w-1, end-1).
		sum -= q[end - 1];
		sum += q[end - w - 1];
		end--;
	}
	return {0, end};
}

TrimResult SlidingWindowTrimmer::trim_5p(const std::uint8_t *q, std::size_t n, std::size_t w,
                                         std::uint8_t mean_quality) {
	require_window(w);
	if (n < w) {
		return {0, n};
	}

	std::uint32_t sum = initial_window_sum(q, w);
	std::size_t start = 0; // current window is [start, start+w)
	const std::uint32_t threshold = threshold_sum(mean_quality, w);

	while (sum < threshold) {
		if (start + w == n) {
			return {n, n};
		}
		// Slide forward by 1: window becomes [start+1, start+w+1).
		sum -= q[start];
		sum += q[start + w];
		start++;
	}
	return {start, n};
}

TrimResult SlidingWindowTrimmer::trim_sliding(const std::uint8_t *q, std::size_t n, std::size_t w,
                                              std::uint8_t mean_quality) {
	require_window(w);
	if (n < w) {
		return {0, n};
	}

	std::uint32_t sum = initial_window_sum(q, w);
	std::size_t start = 0;
	const std::uint32_t threshold = threshold_sum(mean_quality, w);

	while (true) {
		if (sum < threshold) {
			return {0, start};
		}
		if (start + w == n) {
			return {0, n};
		}
		sum -= q[start];
		sum += q[start + w];
		start++;
	}
}

} // namespace miint::qc
