// Portions of this file are derived from fastp
// (https://github.com/OpenGene/fastp), MIT-licensed,
// Copyright (c) 2016 OpenGene. See THIRD_PARTY_LICENSES.md.

#include "qc_algorithms.hpp"

#include <algorithm>
#include <cstdint>
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

// ---------------------------------------------------------------------------
// PolyXScanner
//
// Both scanners follow fastp's structure: walk from the 3' end, tracking
// either a single target base (polyG) or per-base counts (polyX), with a
// tolerance budget of one mismatch per 8 scanned bases plus a hard cap.
// Stops as soon as the tail no longer looks like a homopolymer.
// ---------------------------------------------------------------------------

namespace {

constexpr std::uint32_t POLY_ONE_PER_8 = 8;

// ACGT lookup (case-insensitive). Returns -1 for N or unknown bases.
inline int base_to_acgt_index(std::uint8_t c) {
	switch (c) {
	case 'A':
	case 'a':
		return 0;
	case 'C':
	case 'c':
		return 1;
	case 'G':
	case 'g':
		return 2;
	case 'T':
	case 't':
		return 3;
	default:
		return -1;
	}
}

// Compute mean Phred over [start, len). Caller guarantees len > start.
inline std::uint8_t mean_phred(const std::uint8_t *qual, std::size_t start, std::size_t len) {
	std::uint64_t sum = 0;
	for (std::size_t i = start; i < len; i++) {
		sum += qual[i];
	}
	return static_cast<std::uint8_t>(sum / (len - start));
}

} // namespace

TrimResult PolyXScanner::scan_polyg(const std::uint8_t *seq, const std::uint8_t *qual, std::size_t len,
                                    std::size_t min_len, std::uint32_t max_mismatch, std::uint8_t max_window_mean_q) {
	if (min_len == 0 || len < min_len) {
		return {0, len};
	}

	std::uint32_t mismatch = 0;
	std::size_t first_g_pos = len; // sentinel meaning "no G seen yet"
	for (std::size_t i = 0; i < len; i++) {
		const std::uint8_t c = seq[len - 1 - i];
		if (c == 'G' || c == 'g') {
			first_g_pos = len - 1 - i;
		} else {
			mismatch++;
		}
		const std::uint32_t allowed = std::min(max_mismatch, static_cast<std::uint32_t>((i + 1) / POLY_ONE_PER_8));
		if (mismatch > max_mismatch || (mismatch > allowed && i + 1 >= min_len)) {
			break;
		}
	}

	if (first_g_pos == len || len - first_g_pos < min_len) {
		return {0, len};
	}

	// Quality-aware gate. With max_window_mean_q==255 the check is a no-op
	// since any Phred 0..93 satisfies <= 255.
	if (mean_phred(qual, first_g_pos, len) > max_window_mean_q) {
		return {0, len};
	}
	return {0, first_g_pos};
}

TrimResult PolyXScanner::scan_polyx(const std::uint8_t *seq, std::size_t len, std::size_t min_len,
                                    std::uint32_t max_mismatch) {
	if (min_len == 0 || len < min_len) {
		return {0, len};
	}

	std::uint32_t counts[4] = {0, 0, 0, 0};
	std::size_t scanned = 0;
	for (std::size_t i = 0; i < len; i++) {
		const std::uint8_t c = seq[len - 1 - i];
		const int idx = base_to_acgt_index(c);
		if (idx >= 0) {
			counts[idx]++;
		}
		scanned = i + 1;

		const std::uint32_t allowed = std::min(max_mismatch, static_cast<std::uint32_t>(scanned / POLY_ONE_PER_8));
		bool any_base_survives = false;
		for (int b = 0; b < 4; b++) {
			if (scanned - counts[b] <= allowed) {
				any_base_survives = true;
				break;
			}
		}
		if (!any_base_survives && (i + 1 >= POLY_ONE_PER_8 || i + 1 >= min_len)) {
			break;
		}
	}

	// Dominant base — deterministic ACGT tie-break (first index with max count).
	int max_b = 0;
	for (int b = 1; b < 4; b++) {
		if (counts[b] > counts[max_b]) {
			max_b = b;
		}
	}
	if (counts[max_b] == 0) {
		return {0, len};
	}

	static constexpr char ACGT[] = "ACGT";
	const std::uint8_t target = static_cast<std::uint8_t>(ACGT[max_b]);

	// Leftmost occurrence of `target` within the scanned tail.
	std::size_t first_x_pos = len;
	const std::size_t scan_start = len - scanned;
	for (std::size_t i = scan_start; i < len; i++) {
		const std::uint8_t c = seq[i];
		// Case-insensitive: lowercase and uppercase target both qualify.
		if (c == target || c == static_cast<std::uint8_t>(target | 0x20)) {
			first_x_pos = i;
			break;
		}
	}

	if (first_x_pos == len || len - first_x_pos < min_len) {
		return {0, len};
	}
	return {0, first_x_pos};
}

// ---------------------------------------------------------------------------
// AdapterMatcher — 3-phase port of fastp's adapter trimming
// ---------------------------------------------------------------------------
//
// Phase 1: scan candidate positions left-to-right. At each position, compare
// adapter vs seq with a tolerance of cmplen/8 mismatches. Return first match.
//
// Phase 2: same scan but allow exactly one insertion in seq (seq has one
// extra base relative to adapter). Two-pointer walk; on first mismatch, skip
// a seq base and continue. Mismatches AFTER the indel still count.
//
// Phase 3: same as phase 2 but the indel is a deletion in seq (seq missing
// one base). On first mismatch, skip an adapter base.
//
// Phases are tried in order; first to find a match wins.

namespace {

constexpr std::uint32_t ADAPTER_ONE_PER_8 = 8;

// fastp's pre-start offset: -min(4, max(2, adapter_len/2)). For short
// adapters use a smaller offset; cap at -4 for typical 13bp adapters.
inline int pre_start_offset(std::size_t adapter_len) {
	int half = static_cast<int>(adapter_len / 2);
	if (half > 4) {
		half = 4;
	}
	if (half < 2) {
		half = 2;
	}
	return -half;
}

// Phase 1: scan and return first match by exact-Hamming-with-tolerance.
AdapterMatch phase1_hamming(const std::uint8_t *seq, std::size_t seq_len, const std::uint8_t *adapter,
                            std::size_t adapter_len, std::size_t min_match, int start_pos) {
	const int seq_len_i = static_cast<int>(seq_len);
	const int adapter_len_i = static_cast<int>(adapter_len);
	const int min_match_i = static_cast<int>(min_match);

	for (int pos = start_pos; pos + min_match_i <= seq_len_i; pos++) {
		const int adapter_off = pos < 0 ? -pos : 0;
		const int seq_off = pos < 0 ? 0 : pos;
		const int cmplen = std::min(seq_len_i - seq_off, adapter_len_i - adapter_off);
		if (cmplen < min_match_i) {
			continue;
		}

		const std::uint32_t allowed = static_cast<std::uint32_t>(cmplen / ADAPTER_ONE_PER_8);
		std::uint32_t mismatches = 0;
		bool ok = true;
		for (int i = 0; i < cmplen; i++) {
			if (adapter[adapter_off + i] != seq[seq_off + i]) {
				mismatches++;
				if (mismatches > allowed) {
					ok = false;
					break;
				}
			}
		}
		if (ok) {
			AdapterMatch m;
			m.matched = true;
			m.trim_start = static_cast<std::size_t>(seq_off);
			m.match_len = static_cast<std::size_t>(cmplen);
			m.mismatches = mismatches;
			m.indels = 0;
			return m;
		}
	}
	return {};
}

// Phase 2/3 shared body. `insertion_in_seq` selects direction:
//   true  → seq has one extra base (skip seq on first mismatch)
//   false → seq is missing one base (skip adapter on first mismatch)
AdapterMatch phase_indel(const std::uint8_t *seq, std::size_t seq_len, const std::uint8_t *adapter,
                         std::size_t adapter_len, std::size_t min_match, int start_pos, bool insertion_in_seq) {
	const int seq_len_i = static_cast<int>(seq_len);
	const int adapter_len_i = static_cast<int>(adapter_len);
	const int min_match_i = static_cast<int>(min_match);

	for (int pos = start_pos; pos + min_match_i <= seq_len_i; pos++) {
		const int adapter_off = pos < 0 ? -pos : 0;
		const int seq_off = pos < 0 ? 0 : pos;
		const int adapter_remain = adapter_len_i - adapter_off;
		const int seq_remain = seq_len_i - seq_off;
		if (adapter_remain < min_match_i) {
			continue;
		}

		const int seq_region_len = insertion_in_seq ? adapter_remain + 1 : adapter_remain - 1;
		if (seq_region_len < min_match_i) {
			continue;
		}
		if (seq_remain < seq_region_len) {
			continue;
		}

		const std::uint32_t allowed = static_cast<std::uint32_t>(adapter_remain / ADAPTER_ONE_PER_8);
		std::uint32_t mismatches = 0;
		bool indel_used = false;
		int ai = 0;
		int si = 0;
		bool ok = true;
		while (ai < adapter_remain && si < seq_region_len) {
			if (adapter[adapter_off + ai] == seq[seq_off + si]) {
				ai++;
				si++;
				continue;
			}
			if (!indel_used) {
				indel_used = true;
				if (insertion_in_seq) {
					si++;
				} else {
					ai++;
				}
				continue;
			}
			mismatches++;
			if (mismatches > allowed) {
				ok = false;
				break;
			}
			ai++;
			si++;
		}
		if (ok && indel_used) {
			AdapterMatch m;
			m.matched = true;
			m.trim_start = static_cast<std::size_t>(seq_off);
			m.match_len = static_cast<std::size_t>(seq_region_len);
			m.mismatches = mismatches;
			m.indels = 1;
			return m;
		}
	}
	return {};
}

} // namespace

AdapterMatch AdapterMatcher::find(const std::uint8_t *seq, std::size_t seq_len, const std::uint8_t *adapter,
                                  std::size_t adapter_len, std::size_t min_match, bool allow_pre_start) {
	if (min_match == 0 || adapter_len == 0 || seq_len == 0 || seq_len < min_match) {
		return {};
	}

	const int start_pos = allow_pre_start ? pre_start_offset(adapter_len) : 0;

	auto m = phase1_hamming(seq, seq_len, adapter, adapter_len, min_match, start_pos);
	if (m.matched) {
		return m;
	}
	m = phase_indel(seq, seq_len, adapter, adapter_len, min_match, start_pos, /*insertion_in_seq=*/true);
	if (m.matched) {
		return m;
	}
	return phase_indel(seq, seq_len, adapter, adapter_len, min_match, start_pos, /*insertion_in_seq=*/false);
}

} // namespace miint::qc
