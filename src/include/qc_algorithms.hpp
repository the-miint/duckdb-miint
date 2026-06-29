#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

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

// Adapter match record. `trim_start` is the leftmost seq index removed
// (kept range is [0, trim_start)); `match_len` is the number of seq bases
// matched (which differs from adapter_len when an indel is used). For a
// pre-start match (adapter overlaps seq[0] from the left), `trim_start` is
// 0 and the entire read should be dropped.
struct AdapterMatch {
	bool matched = false;
	std::size_t trim_start = 0;
	std::size_t match_len = 0;
	std::uint32_t mismatches = 0;
	std::uint32_t indels = 0;
};

class AdapterMatcher {
public:
	// Hamming-tolerant adapter search ported from fastp's
	// AdapterTrimmer::trimBySequence: scan positions left-to-right, comparing
	// adapter vs seq with a tolerance of 1 mismatch per 8 compared bases, and
	// return the first match. No indel handling (fastp's by-sequence trim has
	// none) — a single inserted/deleted base in the adapter region is not
	// matched. The LEFTMOST match wins (fastp's behavior — the leftmost hit is
	// the most likely true adapter start; "best match" ranking would over-trim
	// on genomic chatter that happens to match late in the read). `indels` in
	// the result is therefore always 0; the field is kept for symmetry with the
	// other match metrics (only tests read it).
	//
	// `min_match` is the minimum length of compared region required for a
	// match. fastp auto-scales this 4..6 based on adapter list size; here we
	// expose it explicitly so the caller can pass an appropriate value.
	//
	// If `allow_pre_start`, the scan starts at a small negative offset so an
	// adapter that begins before seq[0] can still match. A pre-start match
	// returns `trim_start = 0`, indicating the entire read should be dropped.
	// Default off — fastp turns this on unconditionally for A-tailing, but
	// most non-Illumina protocols don't need it and it can over-trim.
	static AdapterMatch find(const std::uint8_t *seq, std::size_t seq_len, const std::uint8_t *adapter,
	                         std::size_t adapter_len, std::size_t min_match, bool allow_pre_start);

	// Leftmost trim point across a candidate adapter list: runs find() for each
	// candidate and returns the smallest matched trim_start, or `seq_len` if none
	// match (i.e. the kept end of the read). Candidates shorter than `min_match`
	// are skipped. Stops as soon as a candidate matches at position 0 — that is
	// the floor, so no remaining candidate can trim further left. The caller is
	// responsible for choosing `min_match` (e.g. fastp's list-size auto-scale)
	// BEFORE any deduplication, so that collapsing redundant candidates does not
	// change sensitivity.
	static std::size_t find_leftmost(const std::uint8_t *seq, std::size_t seq_len,
	                                 const std::vector<std::string> &candidates, std::size_t min_match,
	                                 bool allow_pre_start);
};

// Result of paired-end overlap analysis (ported from fastp's OverlapAnalysis).
// `offset` is the shift of reverse-complement(R2) relative to R1:
//   offset > 0  — R1 has bases 5' of the overlap; the insert is at least as long
//                 as the read, so neither mate reads into adapter.
//   offset < 0  — revcomp(R2) extends 5' of R1: the insert is SHORTER than the
//                 read length, so each mate reads through into adapter past
//                 `overlap_len`. This is the case fastp trims.
//   offset == 0 — the reads align with no shift.
// `overlap_len` is the length of the aligned region; `diff` the mismatch count
// in it (see OverlapAnalyzer::analyze for the >50bp reporting caveat).
struct OverlapResult {
	bool overlapped = false;
	int offset = 0;
	int overlap_len = 0;
	int diff = 0;
};

class OverlapAnalyzer {
public:
	// 3-parameter no-gap overlap analysis ported from fastp's
	// OverlapAnalysis::analyze. `seq2` is reverse-complemented internally using
	// fastp's complement (A<->T, C<->G; any other byte -> 'N'), so reads may
	// contain 'N'. Scans forward offsets [0, len1 - overlap_require) then reverse
	// offsets down to -(len2 - overlap_require) + 1, and returns the FIRST offset
	// whose mismatch count within the first min(overlap_len, 50) bases does not
	// exceed the per-offset limit min(diff_limit, overlap_len * diff_percent_limit / 100).
	//
	// Faithfully ported quirk: acceptance is decided on only the first 50 bases,
	// but for an accepted overlap longer than 50 the reported `diff` is the FULL
	// mismatch count (which may exceed the limit). fastp's one-gap path is not
	// ported (no-gap only).
	//
	// diff_limit / overlap_require / diff_percent_limit mirror fastp's
	// --overlap_diff_limit (default 5), --overlap_len_require (30),
	// --overlap_diff_percent_limit (20). Returns overlapped=false when no offset
	// qualifies (including when either read is shorter than overlap_require).
	static OverlapResult analyze(const char *seq1, std::size_t len1, const char *seq2, std::size_t len2, int diff_limit,
	                             int overlap_require, int diff_percent_limit);
};

// Per-read filter metrics — computed once in a single pass over seq+qual.
// All four threshold checks (low-qual %, avg quality, N count, length) are
// derived from these fields, so the user can audit failure modes via the
// scalar's STRUCT return without re-walking the bases.
struct FilterMetrics {
	std::uint32_t length;
	std::uint32_t n_bases;        // count of seq bases == 'N' (case-insensitive)
	std::uint32_t low_qual_bases; // count of qual bytes strictly < qualified_q
	std::uint64_t qual_sum;       // sum of Phred values (uint64 keeps headroom even for huge inputs)
};

class ReadFilter {
public:
	// Single pass over seq+qual computing all metrics. Threshold application
	// happens at the scalar layer so the metric struct can be reused for
	// auditing and for fastp-style precedence ordering.
	//
	// `seq` may be nullptr iff `len == 0`. Likewise for `qual`.
	static FilterMetrics measure(const std::uint8_t *seq, const std::uint8_t *qual, std::size_t len,
	                             std::uint8_t qualified_q);
};

class PolyXScanner {
public:
	// scan_polyg: identify a polyG run at the 3' end and return the trim point.
	// Ported from fastp's polyG logic with one improvement: a quality-aware
	// gate. After identifying a candidate trim region, the mean Phred of that
	// region must be <= max_window_mean_q or the trim is refused. Pass
	// max_window_mean_q = 93 (the maximum valid Phred score) to make the gate
	// a no-op — every real Phred value satisfies <= 93.
	//
	// `qual` must point to a buffer of `len` bytes parallel to `seq`.
	static TrimResult scan_polyg(const std::uint8_t *seq, const std::uint8_t *qual, std::size_t len,
	                             std::size_t min_len, std::uint32_t max_mismatch, std::uint8_t max_window_mean_q);

	// scan_polyx: identify the dominant base in the 3' tail and trim from its
	// leftmost position within the scanned region. Ports fastp's polyX scan
	// with one improvement: deterministic ACGT tie-break (earliest in ACGT
	// wins ties) so the result is reproducible across runs and platforms.
	// No quality gate — polyA/T/C runs are not platform artifacts the way
	// NextSeq polyG is.
	static TrimResult scan_polyx(const std::uint8_t *seq, std::size_t len, std::size_t min_len,
	                             std::uint32_t max_mismatch);
};

} // namespace miint::qc
