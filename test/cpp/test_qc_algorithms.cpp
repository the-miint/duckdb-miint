#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <qc_algorithms.hpp>

using miint::qc::AdapterMatch;
using miint::qc::AdapterMatcher;
using miint::qc::FilterMetrics;
using miint::qc::PolyXScanner;
using miint::qc::ReadFilter;
using miint::qc::SlidingWindowTrimmer;
using miint::qc::TrimResult;

// Helper: byte view over a string literal.
static const std::uint8_t *bp(const std::string &s) {
	return reinterpret_cast<const std::uint8_t *>(s.data());
}

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

	SECTION("quality-aware disabled (sentinel = max valid Phred 93) trims high-quality G run") {
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

// ---------------------------------------------------------------------------
// AdapterMatcher::find — 3-phase fastp port (exact Hamming, +1 insert, +1 delete)
// ---------------------------------------------------------------------------
TEST_CASE("AdapterMatcher::find phase 1 (exact Hamming)", "[qc][adapter]") {
	const std::string adapter = "AGATCGGAAGAGC"; // 13bp, TruSeq R1

	SECTION("exact match at 3' end") {
		const std::string seq = "ACGTACGT" + adapter; // 21bp; adapter starts at 8
		auto m = AdapterMatcher::find(bp(seq), seq.size(), bp(adapter), adapter.size(), 4, false);
		REQUIRE(m.matched);
		CHECK(m.trim_start == 8);
		CHECK(m.match_len == 13);
		CHECK(m.mismatches == 0);
		CHECK(m.indels == 0);
	}

	SECTION("exact match at 3' end with one allowed mismatch") {
		// Adapter at pos 8 with one mismatch in middle. cmplen=13, allowed=1.
		std::string seq = "ACGTACGTAGATCGGAAGAGC";
		seq[14] = 'T'; // was 'G' — one mismatch
		auto m = AdapterMatcher::find(bp(seq), seq.size(), bp(adapter), adapter.size(), 4, false);
		REQUIRE(m.matched);
		CHECK(m.trim_start == 8);
		CHECK(m.match_len == 13);
		CHECK(m.mismatches == 1);
	}

	SECTION("too many mismatches — no match") {
		// cmplen=13, allowed=1, but 3 mismatches.
		std::string seq = "ACGTACGTAGATCGGAAGAGC";
		seq[14] = 'T';
		seq[15] = 'T';
		seq[16] = 'T';
		auto m = AdapterMatcher::find(bp(seq), seq.size(), bp(adapter), adapter.size(), 4, false);
		CHECK_FALSE(m.matched);
	}

	SECTION("partial match at very end — at least min_match bases") {
		// Only the first 5 bases of the adapter present at the read end.
		const std::string seq = "ACGTACGTACGTAGATC";
		auto m = AdapterMatcher::find(bp(seq), seq.size(), bp(adapter), adapter.size(), 4, false);
		REQUIRE(m.matched);
		CHECK(m.trim_start == 12);
		CHECK(m.match_len == 5);
		CHECK(m.mismatches == 0);
	}

	SECTION("partial match below min_match — no match") {
		// Only 3 bases visible; min_match=4 — should NOT match.
		const std::string seq = "ACGTACGTACGTACAGA";
		auto m = AdapterMatcher::find(bp(seq), seq.size(), bp(adapter), adapter.size(), 4, false);
		CHECK_FALSE(m.matched);
	}

	SECTION("leftmost match wins when multiple possible positions") {
		// Plant the adapter twice; leftmost (most 5') wins — that's the
		// biologically correct adapter start.
		const std::string seq = "AAAAAAAAAGATCGGAAGAGCAAAGATCGGAAGAGC";
		auto m = AdapterMatcher::find(bp(seq), seq.size(), bp(adapter), adapter.size(), 4, false);
		REQUIRE(m.matched);
		CHECK(m.trim_start == 8); // earlier of the two matches
	}

	SECTION("no match found") {
		const std::string seq = "ACGTACGTACGTACGTACGT";
		auto m = AdapterMatcher::find(bp(seq), seq.size(), bp(adapter), adapter.size(), 4, false);
		CHECK_FALSE(m.matched);
	}
}

TEST_CASE("AdapterMatcher::find phase 2 (insertion in seq)", "[qc][adapter]") {
	const std::string adapter = "AGATCGGAAGAGC"; // 13bp

	SECTION("seq has one inserted base in the middle of the adapter") {
		// "AGATC" + 'X' + "GGAAGAGC" — adapter with extra X at position 5
		const std::string seq = "ACGTACGTAGATCXGGAAGAGC"; // 22bp; adapter region starts at 8
		auto m = AdapterMatcher::find(bp(seq), seq.size(), bp(adapter), adapter.size(), 4, false);
		REQUIRE(m.matched);
		CHECK(m.trim_start == 8);
		CHECK(m.indels == 1);
		CHECK(m.match_len == 14); // adapter_len + 1
	}

	SECTION("insertion + 1 mismatch within tolerance") {
		// Insert + one substitution; cmplen=14 (adapter+1), allowed=14/8=1.
		// Wait: allowed is based on adapter_len for indel phases. allowed=13/8=1.
		std::string seq = "ACGTACGTAGATCXGGAAGAGC";
		seq[16] = 'T'; // was 'G' — one mismatch after the insertion
		auto m = AdapterMatcher::find(bp(seq), seq.size(), bp(adapter), adapter.size(), 4, false);
		REQUIRE(m.matched);
		CHECK(m.trim_start == 8);
		CHECK(m.indels == 1);
		CHECK(m.mismatches == 1);
	}

	SECTION("substitution BEFORE true insertion site — exhaustive search required") {
		// Adapter (13): "AGATCGGAAGAGC"
		// Seq region (14): "TGATCXGGAAGAGC" — substitution T@0 + insertion X@5
		// A greedy commit-on-first-mismatch algorithm wastes the indel slot
		// on position 0 (the T-vs-A substitution) and then runs out of
		// mismatch budget. Exhaustive search across indel positions finds
		// k=5 with 1 mismatch (the T at position 0) and matches.
		const std::string seq = std::string("ACGT") + "TGATCXGGAAGAGC";
		auto m = AdapterMatcher::find(bp(seq), seq.size(), bp(adapter), adapter.size(), 4, false);
		REQUIRE(m.matched);
		CHECK(m.trim_start == 4);
		CHECK(m.indels == 1);
		CHECK(m.mismatches == 1);
	}
}

TEST_CASE("AdapterMatcher::find phase 3 (deletion in seq)", "[qc][adapter]") {
	const std::string adapter = "AGATCGGAAGAGC"; // 13bp

	SECTION("seq is missing one base from the adapter") {
		// Adapter minus the 'C' at position 5 → "AGATCGAAGAGC" (12bp)
		const std::string seq = "ACGTACGTAGATCGAAGAGC"; // 20bp; adapter region at 8
		auto m = AdapterMatcher::find(bp(seq), seq.size(), bp(adapter), adapter.size(), 4, false);
		REQUIRE(m.matched);
		CHECK(m.trim_start == 8);
		CHECK(m.indels == 1);
		CHECK(m.match_len == 12); // adapter_len - 1
	}
}

TEST_CASE("AdapterMatcher::find pre-start behavior", "[qc][adapter]") {
	const std::string adapter = "AGATCGGAAGAGC"; // 13bp

	SECTION("adapter starting exactly at seq[0] matches without pre-start") {
		// Adapter literally starts at the read's first base. This is a normal
		// phase-1 match at pos=0 (not a pre-start case), so allow_pre_start has
		// no effect here. trim_start=0 means the whole read is trimmed.
		const std::string seq = adapter + "ACGTACGT";
		auto m = AdapterMatcher::find(bp(seq), seq.size(), bp(adapter), adapter.size(), 4, false);
		REQUIRE(m.matched);
		CHECK(m.trim_start == 0);
	}

	SECTION("pre-start adapter overlap with allow_pre_start=true drops whole read") {
		// Adapter overlaps seq[0] by being partially in the seq: the first ~3 bases
		// of the adapter are missing (they conceptually live to the LEFT of seq[0]).
		// Read = adapter[3..] + tail
		const std::string read_visible = adapter.substr(3) + "ACGT"; // 10+4=14bp
		auto m = AdapterMatcher::find(bp(read_visible), read_visible.size(), bp(adapter), adapter.size(), 4, true);
		REQUIRE(m.matched);
		// Pre-start matches anchor at seq position 0 — whole read gets trimmed.
		CHECK(m.trim_start == 0);
	}

	SECTION("pre-start adapter overlap with allow_pre_start=false — no match") {
		// Same seq but without the flag. The visible adapter region (10bp) is
		// shorter than the full adapter, but pos=0 + cmplen=10 with adapter[0..10)
		// would need to match seq[0..10) which it doesn't (seq has adapter[3..]).
		// So no match found.
		const std::string read_visible = adapter.substr(3) + "ACGT";
		auto m = AdapterMatcher::find(bp(read_visible), read_visible.size(), bp(adapter), adapter.size(), 4, false);
		CHECK_FALSE(m.matched);
	}
}

// ---------------------------------------------------------------------------
// ReadFilter::measure — single-pass metric computation
// ---------------------------------------------------------------------------
TEST_CASE("ReadFilter::measure", "[qc][filter]") {
	SECTION("mixed quality buffer counts low-qual, N, and sums correctly") {
		// Seq:  A C G T N A C G T N (10 bases, 2 N's)
		// Qual: 5 5 5 5 5 40 40 40 40 40 (5 low-qual, 5 high-qual; qualified_q=15)
		// Sum = 25 + 200 = 225
		auto s = seq_bytes("ACGTNACGTN");
		std::vector<std::uint8_t> q = {5, 5, 5, 5, 5, 40, 40, 40, 40, 40};
		auto m = ReadFilter::measure(s.data(), q.data(), 10, 15);
		CHECK(m.length == 10);
		CHECK(m.n_bases == 2);
		CHECK(m.low_qual_bases == 5);
		CHECK(m.qual_sum == 225);
	}

	SECTION("all-N sequence sets n_bases == length") {
		auto s = seq_bytes("NNNNN");
		std::vector<std::uint8_t> q(5, 40);
		auto m = ReadFilter::measure(s.data(), q.data(), 5, 15);
		CHECK(m.n_bases == 5);
		CHECK(m.low_qual_bases == 0);
	}

	SECTION("all high quality: low_qual_bases == 0") {
		auto s = seq_bytes("ACGT");
		std::vector<std::uint8_t> q = {40, 40, 40, 40};
		auto m = ReadFilter::measure(s.data(), q.data(), 4, 15);
		CHECK(m.low_qual_bases == 0);
		CHECK(m.qual_sum == 160);
	}

	SECTION("quality exactly at threshold is NOT low-quality (>= passes)") {
		// qualified_q = 15 means qual < 15 is low-quality. Q15 is qualified.
		std::vector<std::uint8_t> q = {15, 15, 15};
		auto m = ReadFilter::measure(nullptr, q.data(), 0, 15);
		CHECK(m.length == 0); // empty seq path
		// Use a real seq for the actual threshold test:
		auto s = seq_bytes("ACG");
		auto m2 = ReadFilter::measure(s.data(), q.data(), 3, 15);
		CHECK(m2.low_qual_bases == 0);
		auto m3 = ReadFilter::measure(s.data(), q.data(), 3, 16);
		CHECK(m3.low_qual_bases == 3);
	}

	SECTION("case-insensitive N detection (n counts too)") {
		auto s = seq_bytes("AnNg");
		std::vector<std::uint8_t> q = {40, 40, 40, 40};
		auto m = ReadFilter::measure(s.data(), q.data(), 4, 15);
		CHECK(m.n_bases == 2);
	}

	SECTION("empty input — all zeros") {
		auto m = ReadFilter::measure(nullptr, nullptr, 0, 15);
		CHECK(m.length == 0);
		CHECK(m.n_bases == 0);
		CHECK(m.low_qual_bases == 0);
		CHECK(m.qual_sum == 0);
	}
}
