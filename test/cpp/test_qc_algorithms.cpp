#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <qc_algorithms.hpp>

using miint::qc::SlidingWindowTrimmer;
using miint::qc::TrimResult;

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
