#pragma once

#include <cstddef>
#include <cstdint>

// Pure C++ implementations of the QC algorithms (adapter trimming, polyG/polyX
// trimming, sliding-window quality trimming, per-read filtering). Kept free of
// DuckDB types so they are unit-testable via Catch2 in isolation.
//
// Portions of these algorithms are derived from fastp
// (https://github.com/OpenGene/fastp), MIT-licensed, Copyright (c) 2016 OpenGene.
// See THIRD_PARTY_LICENSES.md for the full attribution.

namespace miint::qc {

// Half-open kept range [start, end) of the input. Use start as trimmed_5p,
// and (input_length - end) as trimmed_3p. If start == end, the entire read
// was trimmed.
struct TrimResult {
	std::size_t start;
	std::size_t end;
};

class SlidingWindowTrimmer {
public:
	// trim_3p (fastp's cut_tail): slide window from 3' toward 5'. While the
	// window mean is below `mean_quality`, recede; stop at the first window
	// that passes. Returns [0, end_of_last_passing_window).
	static TrimResult trim_3p(const std::uint8_t *quals, std::size_t len, std::size_t window_size,
	                          std::uint8_t mean_quality);

	// trim_5p (fastp's cut_front): slide window from 5' toward 3'. While the
	// window mean is below `mean_quality`, advance; stop at the first window
	// that passes. Returns [start_of_first_passing_window, len).
	static TrimResult trim_5p(const std::uint8_t *quals, std::size_t len, std::size_t window_size,
	                          std::uint8_t mean_quality);

	// trim_sliding (fastp's cut_right): slide window 5'→3'. At the first
	// window whose mean is below `mean_quality`, drop that window AND
	// everything to its right. Returns [0, first_bad_window_start).
	static TrimResult trim_sliding(const std::uint8_t *quals, std::size_t len, std::size_t window_size,
	                               std::uint8_t mean_quality);
};

} // namespace miint::qc
