#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <qc_algorithms.hpp>

using miint::qc::PolyXScanner;
using miint::qc::SlidingWindowTrimmer;
using miint::qc::TrimResult;

// Helper: convert ASCII sequence string to uint8 buffer (no quality).
static std::vector<std::uint8_t> seq_bytes(const std::string &s) {
	return std::vector<std::uint8_t>(s.begin(), s.end());
}

// ---------------------------------------------------------------------------
// trim_3p — slide window from 3' toward 5'; while window mean < threshold,
// recede; stop at first passing window. Result spans [0, end).
// ---------------------------------------------------------------------------
TEST_CASE("SlidingWindowTrimmer::trim_3p", "[qc][sliding]") {
	SECTION("basic 3' tail trim (last 2 bases removed)") {
		// First 4 bases Q40, last 4 Q5. Window=4, threshold=Q20.
		// Window [4,8) sum=20 < 80 → slide back.
		// Window [3,7) sum=55 < 80 → slide back.
		// Window [2,6) sum=90 ≥ 80 → pass; end=6.
		std::vector<uint8_t> q = {40, 40, 40, 40, 5, 5, 5, 5};
		auto r = SlidingWindowTrimmer::trim_3p(q.data(), q.size(), 4, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 6);
	}

	SECTION("no trim when first 3' window already passes") {
		std::vector<uint8_t> q = {40, 40, 40, 40, 40, 40, 40, 40};
		auto r = SlidingWindowTrimmer::trim_3p(q.data(), q.size(), 4, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 8);
	}

	SECTION("entire read trimmed when no window passes") {
		std::vector<uint8_t> q = {5, 5, 5, 5};
		auto r = SlidingWindowTrimmer::trim_3p(q.data(), q.size(), 4, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 0);
	}

	SECTION("read shorter than window — no trim (fastp behavior)") {
		std::vector<uint8_t> q = {40, 40, 40};
		auto r = SlidingWindowTrimmer::trim_3p(q.data(), q.size(), 4, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 3);
	}

	SECTION("mean exactly equal to threshold passes") {
		// Sum 80 = threshold * window (20 * 4). Not strictly less than → pass.
		std::vector<uint8_t> q = {20, 20, 20, 20};
		auto r = SlidingWindowTrimmer::trim_3p(q.data(), q.size(), 4, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 4);
	}

	SECTION("empty input — no-op") {
		auto r = SlidingWindowTrimmer::trim_3p(nullptr, 0, 4, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 0);
	}

	SECTION("window of size 1, single low-Q base at 3' end") {
		std::vector<uint8_t> q = {30, 30, 10};
		auto r = SlidingWindowTrimmer::trim_3p(q.data(), q.size(), 1, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 2);
	}

	SECTION("window_size of 0 throws") {
		std::vector<uint8_t> q = {30, 30, 30, 30};
		CHECK_THROWS_AS(SlidingWindowTrimmer::trim_3p(q.data(), q.size(), 0, 20), std::invalid_argument);
	}
}

// ---------------------------------------------------------------------------
// trim_5p — slide window from 5' toward 3'; while window mean < threshold,
// advance; stop at first passing window. Result spans [start, n).
// ---------------------------------------------------------------------------
TEST_CASE("SlidingWindowTrimmer::trim_5p", "[qc][sliding]") {
	SECTION("basic 5' head trim (first 2 bases removed)") {
		// The first passing window is [2,6) with quals [5,5,40,40], mean 22.5.
		// Anchoring at the window's LEFT edge means index 2 is kept — even
		// though q[2]=5 is below the per-base threshold. This is fastp's
		// cut_front semantics: window-mean passes => keep the whole window.
		std::vector<uint8_t> q = {5, 5, 5, 5, 40, 40, 40, 40};
		auto r = SlidingWindowTrimmer::trim_5p(q.data(), q.size(), 4, 20);
		CHECK(r.start == 2);
		CHECK(r.end == 8);
	}

	SECTION("no trim when first 5' window already passes") {
		std::vector<uint8_t> q = {40, 40, 40, 40, 40, 40, 40, 40};
		auto r = SlidingWindowTrimmer::trim_5p(q.data(), q.size(), 4, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 8);
	}

	SECTION("entire read trimmed when no window passes") {
		std::vector<uint8_t> q = {5, 5, 5, 5};
		auto r = SlidingWindowTrimmer::trim_5p(q.data(), q.size(), 4, 20);
		CHECK(r.start == 4);
		CHECK(r.end == 4);
	}

	SECTION("read shorter than window — no trim") {
		std::vector<uint8_t> q = {40, 40, 40};
		auto r = SlidingWindowTrimmer::trim_5p(q.data(), q.size(), 4, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 3);
	}

	SECTION("empty input — no-op") {
		auto r = SlidingWindowTrimmer::trim_5p(nullptr, 0, 4, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 0);
	}

	SECTION("window_size of 0 throws") {
		std::vector<uint8_t> q = {30, 30, 30, 30};
		CHECK_THROWS_AS(SlidingWindowTrimmer::trim_5p(q.data(), q.size(), 0, 20), std::invalid_argument);
	}
}

// ---------------------------------------------------------------------------
// trim_sliding — slide window 5'→3'; AT FIRST failing window, drop window
// AND everything to the right. Result spans [0, first_bad_window_start).
// ---------------------------------------------------------------------------
TEST_CASE("SlidingWindowTrimmer::trim_sliding", "[qc][sliding]") {
	SECTION("drops everything from first failing window onward") {
		// All windows starting at 0,1,2 have sum ≥ 80; window at 3 has sum 55 < 80.
		// Result keeps positions [0, 3).
		std::vector<uint8_t> q = {40, 40, 40, 40, 5, 5, 5, 5};
		auto r = SlidingWindowTrimmer::trim_sliding(q.data(), q.size(), 4, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 3);
	}

	SECTION("no trim when every window passes") {
		std::vector<uint8_t> q = {40, 40, 40, 40, 40, 40, 40, 40};
		auto r = SlidingWindowTrimmer::trim_sliding(q.data(), q.size(), 4, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 8);
	}

	SECTION("entire read trimmed when first window fails") {
		std::vector<uint8_t> q = {5, 5, 5, 5};
		auto r = SlidingWindowTrimmer::trim_sliding(q.data(), q.size(), 4, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 0);
	}

	SECTION("low-quality middle region drops everything from that window") {
		// Window starts: [0,4)=90, [1,5)=55<80 → fail at start=1.
		std::vector<uint8_t> q = {40, 40, 5, 5, 5, 40, 40, 40};
		auto r = SlidingWindowTrimmer::trim_sliding(q.data(), q.size(), 4, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 1);
	}

	SECTION("read shorter than window — no trim") {
		std::vector<uint8_t> q = {5, 5, 5};
		auto r = SlidingWindowTrimmer::trim_sliding(q.data(), q.size(), 4, 20);
		CHECK(r.start == 0);
		CHECK(r.end == 3);
	}

	SECTION("window_size of 0 throws") {
		std::vector<uint8_t> q = {30, 30, 30, 30};
		CHECK_THROWS_AS(SlidingWindowTrimmer::trim_sliding(q.data(), q.size(), 0, 20), std::invalid_argument);
	}
}

// ---------------------------------------------------------------------------
// PolyXScanner::scan_polyg
// ---------------------------------------------------------------------------
// Sentinel quality value that disables the quality-aware check: the max
// valid Phred score is 93, so any real mean satisfies <= 93.
static constexpr std::uint8_t QUAL_AWARE_DISABLED = 93;

TEST_CASE("PolyXScanner::scan_polyg", "[qc][poly]") {
	SECTION("clean 10bp polyG tail — trimmed") {
		auto s = seq_bytes("AAGGGGGGGGGG"); // 2 A's + 10 G's = 12
		std::vector<std::uint8_t> q(12, 2); // Q2 throughout (matches NextSeq dark cycles)
		auto r = PolyXScanner::scan_polyg(s.data(), q.data(), s.size(), 10, 5, 5);
		CHECK(r.start == 0);
		CHECK(r.end == 2);
	}

	SECTION("polyG with one internal mismatch within tolerance — trimmed") {
		// Tail: G G G G A G G G G G G G   (12 bases; one A at position 7 from end)
		auto s = seq_bytes("AAGGGGGAGGGGGG"); // 2 A's + 6 G's + A + 5 G's + 1 G ... wait
		// Let me reconstruct: prefix "AA" (kept) + tail of 12 bases with one internal mismatch
		// Tail reading right-to-left: G,G,G,G,A,G,G,G,G,G,G,G => "GGGGGGGAGGGG"
		auto seq_str = std::string("AA") + std::string("GGGGGGGAGGGG"); // total 14
		s = seq_bytes(seq_str);
		std::vector<std::uint8_t> qual(s.size(), 2);
		auto r = PolyXScanner::scan_polyg(s.data(), qual.data(), s.size(), 10, 5, 5);
		CHECK(r.start == 0);
		CHECK(r.end == 2);
	}

	SECTION("polyG shorter than min_len — no trim") {
		auto s = seq_bytes("AAGGGGGGGGG"); // 2 A's + 9 G's = 11; min_len=10 means need >=10 trimmed
		std::vector<std::uint8_t> q(11, 2);
		auto r = PolyXScanner::scan_polyg(s.data(), q.data(), s.size(), 10, 5, QUAL_AWARE_DISABLED);
		CHECK(r.start == 0);
		CHECK(r.end == 11);
	}

	SECTION("no polyG present — no trim") {
		auto s = seq_bytes("ACGTACGTACGTACGT");
		std::vector<std::uint8_t> q(s.size(), 40);
		auto r = PolyXScanner::scan_polyg(s.data(), q.data(), s.size(), 10, 5, QUAL_AWARE_DISABLED);
		CHECK(r.start == 0);
		CHECK(r.end == 16);
	}

	SECTION("quality-aware refuses to trim high-quality G run") {
		// Tail is 10 G's BUT they are all Q40 (legitimate G-rich genome)
		auto s = seq_bytes("AAGGGGGGGGGG");
		std::vector<std::uint8_t> q(12, 40);
		// max_window_mean_q=5: mean of trim region (40) > 5 → refuse trim
		auto r = PolyXScanner::scan_polyg(s.data(), q.data(), s.size(), 10, 5, 5);
		CHECK(r.start == 0);
		CHECK(r.end == 12);
	}

	SECTION("quality-aware allows trim when quality is low") {
		auto s = seq_bytes("AAGGGGGGGGGG");
		std::vector<std::uint8_t> q(12, 2);
		auto r = PolyXScanner::scan_polyg(s.data(), q.data(), s.size(), 10, 5, 5);
		CHECK(r.start == 0);
		CHECK(r.end == 2);
	}

	SECTION("quality-aware disabled (sentinel 255) trims high-quality G run") {
		auto s = seq_bytes("AAGGGGGGGGGG");
		std::vector<std::uint8_t> q(12, 40);
		auto r = PolyXScanner::scan_polyg(s.data(), q.data(), s.size(), 10, 5, QUAL_AWARE_DISABLED);
		CHECK(r.start == 0);
		CHECK(r.end == 2);
	}

	SECTION("empty input — no-op") {
		auto r = PolyXScanner::scan_polyg(nullptr, nullptr, 0, 10, 5, QUAL_AWARE_DISABLED);
		CHECK(r.start == 0);
		CHECK(r.end == 0);
	}

	SECTION("read shorter than min_len — no trim") {
		auto s = seq_bytes("GGGGG");
		std::vector<std::uint8_t> q(s.size(), 2);
		auto r = PolyXScanner::scan_polyg(s.data(), q.data(), s.size(), 10, 5, QUAL_AWARE_DISABLED);
		CHECK(r.start == 0);
		CHECK(r.end == 5);
	}

	SECTION("custom min_len=4 trims short polyG") {
		auto s = seq_bytes("ACGTGGGG");
		std::vector<std::uint8_t> q(s.size(), 2);
		auto r = PolyXScanner::scan_polyg(s.data(), q.data(), s.size(), 4, 5, QUAL_AWARE_DISABLED);
		CHECK(r.start == 0);
		CHECK(r.end == 4);
	}
}

// ---------------------------------------------------------------------------
// PolyXScanner::scan_polyx
// ---------------------------------------------------------------------------
TEST_CASE("PolyXScanner::scan_polyx", "[qc][poly]") {
	SECTION("polyA tail — trimmed") {
		auto s = seq_bytes("CCAAAAAAAAAA"); // 2 C + 10 A
		auto r = PolyXScanner::scan_polyx(s.data(), s.size(), 10, 5);
		CHECK(r.start == 0);
		CHECK(r.end == 2);
	}

	SECTION("polyT tail — trimmed") {
		auto s = seq_bytes("CCTTTTTTTTTT");
		auto r = PolyXScanner::scan_polyx(s.data(), s.size(), 10, 5);
		CHECK(r.start == 0);
		CHECK(r.end == 2);
	}

	SECTION("polyC tail — trimmed") {
		auto s = seq_bytes("AACCCCCCCCCC");
		auto r = PolyXScanner::scan_polyx(s.data(), s.size(), 10, 5);
		CHECK(r.start == 0);
		CHECK(r.end == 2);
	}

	SECTION("polyG tail via polyX") {
		auto s = seq_bytes("AAGGGGGGGGGG");
		auto r = PolyXScanner::scan_polyx(s.data(), s.size(), 10, 5);
		CHECK(r.start == 0);
		CHECK(r.end == 2);
	}

	SECTION("no homopolymer — no trim") {
		auto s = seq_bytes("ACGTACGTACGTACGT");
		auto r = PolyXScanner::scan_polyx(s.data(), s.size(), 10, 5);
		CHECK(r.start == 0);
		CHECK(r.end == 16);
	}

	SECTION("ACGT tie-break: equal-count tail picks A (earliest in ACGT)") {
		// Tail "AACC": 2 A's and 2 C's tied. ACGT tie-break => A wins.
		// Build a 14-bp read so a polyx scan kicks in at min_len=4.
		auto s = seq_bytes("XXXXXXXXXXAACC"); // 10 X (other) + AACC tail (4)
		auto r = PolyXScanner::scan_polyx(s.data(), s.size(), 4, 5);
		// Dominant base = A; leftmost A within scanned tail is at index 10.
		CHECK(r.start == 0);
		CHECK(r.end == 10);
	}

	SECTION("empty input — no-op") {
		auto r = PolyXScanner::scan_polyx(nullptr, 0, 10, 5);
		CHECK(r.start == 0);
		CHECK(r.end == 0);
	}

	SECTION("read shorter than min_len — no trim") {
		auto s = seq_bytes("AAAA");
		auto r = PolyXScanner::scan_polyx(s.data(), s.size(), 10, 5);
		CHECK(r.start == 0);
		CHECK(r.end == 4);
	}

	SECTION("custom min_len=4 trims short polyA") {
		auto s = seq_bytes("CCGTAAAA");
		auto r = PolyXScanner::scan_polyx(s.data(), s.size(), 4, 5);
		CHECK(r.start == 0);
		CHECK(r.end == 4);
	}
}
