// Portions of this file are derived from fastp
// (https://github.com/OpenGene/fastp), MIT-licensed,
// Copyright (c) 2016 OpenGene. See THIRD_PARTY_LICENSES.md.

#include "qc_algorithms.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

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

// fastp's adapter / poly tolerance: one mismatch budget per 8 scanned bases,
// used by both the polyG/polyX scanners and the adapter matcher's phase loops.
constexpr std::uint32_t MISMATCH_PER_8_BASES = 8;

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
		const std::uint32_t allowed =
		    std::min(max_mismatch, static_cast<std::uint32_t>((i + 1) / MISMATCH_PER_8_BASES));
		if (mismatch > max_mismatch || (mismatch > allowed && i + 1 >= min_len)) {
			break;
		}
	}

	if (first_g_pos == len || len - first_g_pos < min_len) {
		return {0, len};
	}

	// Quality-aware gate. Pass max_window_mean_q=93 (the max valid Phred) to
	// make the check a no-op since any real Phred 0..93 satisfies <= 93.
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

		const std::uint32_t allowed =
		    std::min(max_mismatch, static_cast<std::uint32_t>(scanned / MISMATCH_PER_8_BASES));
		bool any_base_survives = false;
		for (int b = 0; b < 4; b++) {
			if (scanned - counts[b] <= allowed) {
				any_base_survives = true;
				break;
			}
		}
		if (!any_base_survives && (i + 1 >= MISMATCH_PER_8_BASES || i + 1 >= min_len)) {
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

		const std::uint32_t allowed = static_cast<std::uint32_t>(cmplen / MISMATCH_PER_8_BASES);
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

// Phase 2/3 shared body. For each candidate starting position, try every
// possible indel position k in the adapter and pick the best (fewest
// mismatches) result. Exhaustive search across k is required for
// correctness: a greedy commit-on-first-mismatch approach produces false
// negatives when a substitution precedes the true indel (e.g. one leading
// sequencing error before the adapter's actual insertion site exhausts the
// mismatch budget if the indel slot is wasted on the substitution).
//
// `insertion_in_seq` selects direction:
//   true  → seq has one extra base at indel position k; skip seq[k]
//   false → seq is missing one base at adapter position k; skip adapter[k]
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
		// Defensive guard: adapter_remain == 0 would underflow seq_region_len in the deletion case.
		if (adapter_remain <= 0 || adapter_remain < min_match_i) {
			continue;
		}

		const int seq_region_len = insertion_in_seq ? adapter_remain + 1 : adapter_remain - 1;
		if (seq_region_len < min_match_i || seq_remain < seq_region_len) {
			continue;
		}

		const std::uint32_t allowed = static_cast<std::uint32_t>(adapter_remain / MISMATCH_PER_8_BASES);
		// max_k:
		//   insertion: k can be any "gap" in adapter from 0 to adapter_remain inclusive
		//              (inserted base at seq index k means adapter[j] maps to seq[j+1] for j>=k)
		//   deletion:  k can be any adapter index 0..adapter_remain-1 to be "missing" from seq
		const int max_k = insertion_in_seq ? adapter_remain : adapter_remain - 1;
		std::uint32_t best_mm = allowed + 1; // sentinel: "no valid match yet"

		for (int k = 0; k <= max_k; k++) {
			std::uint32_t mm = 0;
			bool ok = true;
			for (int j = 0; j < adapter_remain; j++) {
				if (!insertion_in_seq && j == k) {
					continue; // adapter[k] is the deleted base; not in seq
				}
				int si;
				if (insertion_in_seq) {
					si = j < k ? j : j + 1;
				} else {
					si = j < k ? j : j - 1;
				}
				if (adapter[adapter_off + j] != seq[seq_off + si]) {
					mm++;
					if (mm > allowed) {
						ok = false;
						break;
					}
				}
			}
			if (ok && mm < best_mm) {
				best_mm = mm;
			}
		}

		if (best_mm <= allowed) {
			AdapterMatch m;
			m.matched = true;
			m.trim_start = static_cast<std::size_t>(seq_off);
			m.match_len = static_cast<std::size_t>(seq_region_len);
			m.mismatches = best_mm;
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

std::size_t AdapterMatcher::find_leftmost(const std::uint8_t *seq, std::size_t seq_len,
                                          const std::vector<std::string> &candidates, std::size_t min_match,
                                          bool allow_pre_start) {
	std::size_t best = seq_len; // no match yet -> kept end is the whole read
	for (const auto &cand : candidates) {
		if (cand.size() < min_match) {
			continue;
		}
		const auto m = find(seq, seq_len, reinterpret_cast<const std::uint8_t *>(cand.data()), cand.size(), min_match,
		                    allow_pre_start);
		if (m.matched && m.trim_start < best) {
			best = m.trim_start;
			if (best == 0) {
				break; // 0 is the floor; no other candidate can trim further left
			}
		}
	}
	return best;
}

// ---------------------------------------------------------------------------
// OverlapAnalyzer — paired-end overlap analysis (fastp port)
// ---------------------------------------------------------------------------

namespace {

// fastp's complement: A<->T, C<->G; any other byte (incl. 'N', lowercase, IUPAC
// ambiguity codes) maps to 'N'. Overlap comparison is raw byte equality, so two
// 'N's match but 'N' vs a base does not — exactly fastp's behavior.
inline char overlap_complement(char b) {
	switch (b) {
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'C':
		return 'G';
	case 'G':
		return 'C';
	default:
		return 'N';
	}
}

// Mismatch count over [0, len), early-exiting once it exceeds `limit` (callers
// only test `> limit`, so the exact value past that point is irrelevant).
inline int count_mismatches_bounded(const char *a, const char *b, int len, int limit) {
	int mm = 0;
	for (int i = 0; i < len; i++) {
		if (a[i] != b[i] && ++mm > limit) {
			break;
		}
	}
	return mm;
}

inline int count_mismatches(const char *a, const char *b, int len) {
	int mm = 0;
	for (int i = 0; i < len; i++) {
		if (a[i] != b[i]) {
			mm++;
		}
	}
	return mm;
}

// fastp's complete_compare_require: acceptance is judged on at most this many
// leading bases of the overlap.
constexpr int kCompleteCompareRequire = 50;

// Accept the [a, a+len) vs [b, b+len) overlap iff its mismatch count within the
// first min(len, 50) bases does not exceed `diff_limit`, writing the mismatch
// count to `mm`. For an accepted overlap longer than 50 bases `mm` is updated to
// the FULL count (fastp judges acceptance on only the first 50 — see header).
inline bool accept_no_gap(const char *a, const char *b, int len, int diff_limit, int &mm) {
	const int prefix = std::min(len, kCompleteCompareRequire);
	mm = count_mismatches_bounded(a, b, prefix, diff_limit);
	if (mm > diff_limit) {
		return false;
	}
	if (len > kCompleteCompareRequire) {
		mm = count_mismatches(a, b, len);
	}
	return true;
}

} // namespace

OverlapResult OverlapAnalyzer::analyze(const char *seq1, std::size_t len1_sz, const char *seq2, std::size_t len2_sz,
                                       int diff_limit, int overlap_require, int diff_percent_limit) {
	OverlapResult ov;

	// Degenerate / out-of-range inputs have no meaningful overlap. The scalar
	// layer also validates user-supplied params, but guarding here keeps the
	// pure function safe: avoids signed overflow from the size_t->int casts and
	// the nonsense zero-length "overlap" a non-positive overlap_require would
	// otherwise accept in the reverse scan.
	if (overlap_require <= 0 || diff_limit < 0 || len1_sz > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
	    len2_sz > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
		return ov;
	}

	const int len1 = static_cast<int>(len1_sz);
	const int len2 = static_cast<int>(len2_sz);

	// fastp aligns R1 against the reverse complement of R2.
	std::string rc(static_cast<std::size_t>(len2), 'N');
	for (int i = 0; i < len2; i++) {
		rc[static_cast<std::size_t>(i)] = overlap_complement(seq2[len2 - 1 - i]);
	}
	const char *str1 = seq1;
	const char *str2 = rc.c_str();

	const double diff_percent = static_cast<double>(diff_percent_limit) / 100.0;

	// Forward: revcomp(R2) shifted right by `offset` along R1 (offset >= 0).
	for (int offset = 0; offset < len1 - overlap_require; offset++) {
		const int overlap_len = std::min(len1 - offset, len2);
		const int limit = std::min(diff_limit, static_cast<int>(overlap_len * diff_percent));
		int mm = 0;
		if (accept_no_gap(str1 + offset, str2, overlap_len, limit, mm)) {
			ov.overlapped = true;
			ov.offset = offset;
			ov.overlap_len = overlap_len;
			ov.diff = mm;
			return ov;
		}
	}

	// Reverse: revcomp(R2) shifted left (offset < 0) — the insert is shorter than
	// the read, so each mate reads through into adapter past the overlap.
	for (int offset = 0; offset > -(len2 - overlap_require); offset--) {
		const int overlap_len = std::min(len1, len2 + offset); // offset < 0: len2 - |offset|
		const int limit = std::min(diff_limit, static_cast<int>(overlap_len * diff_percent));
		int mm = 0;
		if (accept_no_gap(str1, str2 - offset, overlap_len, limit, mm)) { // str2 + |offset|
			ov.overlapped = true;
			ov.offset = offset;
			ov.overlap_len = overlap_len;
			ov.diff = mm;
			return ov;
		}
	}

	return ov;
}

// ---------------------------------------------------------------------------
// ReadFilter — single-pass per-read metric computation
// ---------------------------------------------------------------------------

FilterMetrics ReadFilter::measure(const std::uint8_t *seq, const std::uint8_t *qual, std::size_t len,
                                  std::uint8_t qualified_q) {
	FilterMetrics m {};
	m.length = static_cast<std::uint32_t>(len);
	for (std::size_t i = 0; i < len; i++) {
		const std::uint8_t s = seq[i];
		// Case-insensitive N detection: 'N' (0x4E) and 'n' (0x6E) differ only
		// in bit 5, so masking it off lets one comparison match both.
		if ((s & 0xDF) == 'N') {
			m.n_bases++;
		}
		const std::uint8_t q = qual[i];
		if (q < qualified_q) {
			m.low_qual_bases++;
		}
		m.qual_sum += q;
	}
	return m;
}

} // namespace miint::qc
