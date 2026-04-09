#include <catch2/catch_test_macros.hpp>
#include <AlignmentSlicer.hpp>
#include <alignment_functions_internal.hpp>

using namespace miint;

// Helper to create a simple SliceInput (no sequence/quality)
static SliceInput make_input(const std::string &cigar, int64_t pos, int64_t stop_pos) {
	return {cigar, pos, stop_pos, "", {}};
}

// Helper to create SliceInput with sequence and quality
static SliceInput make_input_sq(const std::string &cigar, int64_t pos, int64_t stop_pos, const std::string &seq,
                                const std::vector<uint8_t> &qual) {
	return {cigar, pos, stop_pos, seq, qual};
}

// =====================================================================
// 1a. Overlap detection tests
// =====================================================================

TEST_CASE("AlignmentSlicer overlap - M operations", "[alignment_slicer][overlap]") {
	SECTION("M overlaps region") {
		// 10M at pos 5 covers [5,15). Region [8,12) overlaps.
		AlignmentSlicer slicer(8, 12, false);
		auto result = slicer.Slice(make_input("10M", 5, 15));
		REQUIRE(result.overlaps);
	}

	SECTION("M fully before region") {
		// 5M at pos 1 covers [1,6). Region [10,20) — no overlap.
		AlignmentSlicer slicer(10, 20, false);
		auto result = slicer.Slice(make_input("5M", 1, 6));
		REQUIRE_FALSE(result.overlaps);
	}

	SECTION("M fully after region") {
		// 5M at pos 25 covers [25,30). Region [10,20) — no overlap.
		AlignmentSlicer slicer(10, 20, false);
		auto result = slicer.Slice(make_input("5M", 25, 30));
		REQUIRE_FALSE(result.overlaps);
	}

	SECTION("M at left boundary") {
		// 5M at pos 10 covers [10,15). Region [10,20) — overlaps.
		AlignmentSlicer slicer(10, 20, false);
		auto result = slicer.Slice(make_input("5M", 10, 15));
		REQUIRE(result.overlaps);
	}

	SECTION("M just past right boundary") {
		// 5M at pos 20 covers [20,25). Region [10,20) — no overlap (stop is exclusive).
		AlignmentSlicer slicer(10, 20, false);
		auto result = slicer.Slice(make_input("5M", 20, 25));
		REQUIRE_FALSE(result.overlaps);
	}
}

TEST_CASE("AlignmentSlicer overlap - = and X operations", "[alignment_slicer][overlap]") {
	SECTION("= overlaps") {
		// 5=5M at pos 5 covers [5,15). Region [6,8) overlaps.
		AlignmentSlicer slicer(6, 8, false);
		auto result = slicer.Slice(make_input("5=5M", 5, 15));
		REQUIRE(result.overlaps);
	}

	SECTION("X overlaps") {
		// 3M5X at pos 1 covers [1,9). Region [5,8) overlaps via X at [4,9).
		AlignmentSlicer slicer(5, 8, false);
		auto result = slicer.Slice(make_input("3M5X", 1, 9));
		REQUIRE(result.overlaps);
	}

	SECTION("= fully before") {
		AlignmentSlicer slicer(10, 20, false);
		auto result = slicer.Slice(make_input("5=", 1, 6));
		REQUIRE_FALSE(result.overlaps);
	}

	SECTION("X fully after") {
		AlignmentSlicer slicer(10, 20, false);
		auto result = slicer.Slice(make_input("5X", 25, 30));
		REQUIRE_FALSE(result.overlaps);
	}
}

TEST_CASE("AlignmentSlicer overlap - I operations", "[alignment_slicer][overlap]") {
	SECTION("I at region boundary - M provides overlap") {
		// 5M2I5M at pos 1. 5M covers [1,6), I at ref_pos 6, 5M covers [6,11).
		// Region [6,12): the second 5M overlaps.
		AlignmentSlicer slicer(6, 12, false);
		auto result = slicer.Slice(make_input("5M2I5M", 1, 11));
		REQUIRE(result.overlaps);
	}

	SECTION("I before region - last M overlaps") {
		// 3M2I3M at pos 1. 3M covers [1,4), I at ref_pos 4, 3M covers [4,7).
		// Region [7,20): no M overlaps [7,20) since second 3M covers [4,7).
		// Actually second 3M covers ref [4,7), so stop is 7. Region starts at 7 (exclusive). No overlap!
		AlignmentSlicer slicer(7, 20, false);
		auto result = slicer.Slice(make_input("3M2I3M", 1, 7));
		REQUIRE_FALSE(result.overlaps);
	}

	SECTION("M at pos 5 is in region [5,6)") {
		// 1M2I at pos 5. M covers [5,6). Region [5,6) — M overlaps.
		AlignmentSlicer slicer(5, 6, false);
		auto result = slicer.Slice(make_input("1M2I", 5, 6));
		REQUIRE(result.overlaps);
	}
}

TEST_CASE("AlignmentSlicer overlap - D operations", "[alignment_slicer][overlap]") {
	SECTION("D-only overlap, exclude deletions") {
		// 3M5D3M at pos 1. 3M covers [1,4), 5D covers [4,9), 3M covers [9,12).
		// Region [5,9): only D overlaps. With include_deletions=false, no overlap.
		AlignmentSlicer slicer(5, 9, false);
		auto result = slicer.Slice(make_input("3M5D3M", 1, 12));
		REQUIRE_FALSE(result.overlaps);
	}

	SECTION("D-only overlap, include deletions") {
		AlignmentSlicer slicer(5, 9, true);
		auto result = slicer.Slice(make_input("3M5D3M", 1, 12));
		REQUIRE(result.overlaps);
	}

	SECTION("D with M overlap") {
		// 3M5D3M at pos 1. Region [2,5): M at [1,4) overlaps [2,5).
		AlignmentSlicer slicer(2, 5, false);
		auto result = slicer.Slice(make_input("3M5D3M", 1, 12));
		REQUIRE(result.overlaps);
	}
}

TEST_CASE("AlignmentSlicer overlap - N operations", "[alignment_slicer][overlap]") {
	SECTION("N-only overlap, exclude deletions") {
		// 3M10N3M at pos 1. 3M covers [1,4), 10N covers [4,14), 3M covers [14,17).
		// Region [5,12): only N overlaps. N never counts as overlap.
		AlignmentSlicer slicer(5, 12, false);
		auto result = slicer.Slice(make_input("3M10N3M", 1, 17));
		REQUIRE_FALSE(result.overlaps);
	}

	SECTION("N-only overlap, include deletions - still no overlap") {
		// N never counts as overlap regardless of include_deletions flag.
		AlignmentSlicer slicer(5, 12, true);
		auto result = slicer.Slice(make_input("3M10N3M", 1, 17));
		REQUIRE_FALSE(result.overlaps);
	}

	SECTION("N with M overlap") {
		// Region [1,5): first 3M at [1,4) overlaps.
		AlignmentSlicer slicer(1, 5, false);
		auto result = slicer.Slice(make_input("3M10N3M", 1, 17));
		REQUIRE(result.overlaps);
	}
}

TEST_CASE("AlignmentSlicer overlap - S operations", "[alignment_slicer][overlap]") {
	SECTION("S only overlap") {
		// 5S5M at pos 10. S doesn't consume ref. 5M covers [10,15).
		// Region [5,10): no M in region (M starts at 10, region stops at 10 exclusive).
		AlignmentSlicer slicer(5, 10, false);
		auto result = slicer.Slice(make_input("5S5M", 10, 15));
		REQUIRE_FALSE(result.overlaps);
	}

	SECTION("S before M in region") {
		// 3S5M at pos 4. S doesn't consume ref. 5M covers [4,9).
		// Region [4,9): M overlaps.
		AlignmentSlicer slicer(4, 9, false);
		auto result = slicer.Slice(make_input("3S5M", 4, 9));
		REQUIRE(result.overlaps);
	}
}

TEST_CASE("AlignmentSlicer overlap - H operations", "[alignment_slicer][overlap]") {
	SECTION("H only - no query overlap") {
		// 5H5M at pos 5. H doesn't consume ref. 5M covers [5,10).
		// Region [1,5): no M in region.
		AlignmentSlicer slicer(1, 5, false);
		auto result = slicer.Slice(make_input("5H5M", 5, 10));
		REQUIRE_FALSE(result.overlaps);
	}

	SECTION("H with M overlap") {
		// 5H5M at pos 5. 5M covers [5,10). Region [5,10) overlaps.
		AlignmentSlicer slicer(5, 10, false);
		auto result = slicer.Slice(make_input("5H5M", 5, 10));
		REQUIRE(result.overlaps);
	}
}

TEST_CASE("AlignmentSlicer overlap - P operation", "[alignment_slicer][overlap]") {
	// 5M1P5M at pos 1. P is padding, ignored. 5M covers [1,6), 5M covers [6,11).
	// Region [3,8) overlaps.
	AlignmentSlicer slicer(3, 8, false);
	auto result = slicer.Slice(make_input("5M1P5M", 1, 11));
	REQUIRE(result.overlaps);
}

TEST_CASE("AlignmentSlicer overlap - unmapped", "[alignment_slicer][overlap]") {
	AlignmentSlicer slicer(1, 100, false);
	auto result = slicer.Slice(make_input("*", 0, 0));
	REQUIRE_FALSE(result.overlaps);
}

TEST_CASE("AlignmentSlicer overlap - complex CIGAR with all ops", "[alignment_slicer][overlap]") {
	// 2H3S5M2I3M3D4M10N2M3S2H at pos 10
	// Walk: H(2) doesn't consume ref. S(3) doesn't consume ref.
	// 5M covers [10,15). 2I at ref 15. 3M covers [15,18). 3D covers [18,21).
	// 4M covers [21,25). 10N covers [25,35). 2M covers [35,37). S(3) doesn't consume ref. H(2) doesn't.
	// Region [12,25): 5M at [10,15) partially overlaps, 3M at [15,18) fully in, etc.
	AlignmentSlicer slicer(12, 25, false);
	auto result = slicer.Slice(make_input("2H3S5M2I3M3D4M10N2M3S2H", 10, 37));
	REQUIRE(result.overlaps);
}

// =====================================================================
// 1b. CIGAR trimming tests
// =====================================================================

TEST_CASE("AlignmentSlicer trim - M operations", "[alignment_slicer][trim]") {
	SECTION("No trim needed") {
		// 10M at pos 5, region [1,100). Fully within.
		AlignmentSlicer slicer(1, 100, false);
		auto result = slicer.Slice(make_input("10M", 5, 15));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "10M");
		REQUIRE(result.position == 5);
		REQUIRE(result.stop_position == 15);
	}

	SECTION("Trim left") {
		// 10M at pos 1 covers [1,11). Region [4,100). Trim 3 from left.
		AlignmentSlicer slicer(4, 100, false);
		auto result = slicer.Slice(make_input("10M", 1, 11));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "3H7M");
		REQUIRE(result.position == 4);
		REQUIRE(result.stop_position == 11);
	}

	SECTION("Trim right") {
		// 10M at pos 1 covers [1,11). Region [1,6). Trim 5 from right.
		AlignmentSlicer slicer(1, 6, false);
		auto result = slicer.Slice(make_input("10M", 1, 11));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "5M5H");
		REQUIRE(result.position == 1);
		REQUIRE(result.stop_position == 6);
	}

	SECTION("Trim both sides") {
		// 10M at pos 1 covers [1,11). Region [4,8). Trim 3 left, 3 right.
		AlignmentSlicer slicer(4, 8, false);
		auto result = slicer.Slice(make_input("10M", 1, 11));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "3H4M3H");
		REQUIRE(result.position == 4);
		REQUIRE(result.stop_position == 8);
	}

	SECTION("Trim removes entire first M op") {
		// 3M5M at pos 1 covers [1,9). Region [4,100). First 3M fully before.
		AlignmentSlicer slicer(4, 100, false);
		auto result = slicer.Slice(make_input("3M5M", 1, 9));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "3H5M");
		REQUIRE(result.position == 4);
		REQUIRE(result.stop_position == 9);
	}
}

TEST_CASE("AlignmentSlicer trim - = and X operations", "[alignment_slicer][trim]") {
	SECTION("= trim left") {
		// 10= at pos 1 covers [1,11). Region [4,100). Trim 3 from left.
		AlignmentSlicer slicer(4, 100, false);
		auto result = slicer.Slice(make_input("10=", 1, 11));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "3H7=");
		REQUIRE(result.position == 4);
		REQUIRE(result.stop_position == 11);
	}

	SECTION("X trim right") {
		// 10X at pos 1 covers [1,11). Region [1,6).
		AlignmentSlicer slicer(1, 6, false);
		auto result = slicer.Slice(make_input("10X", 1, 11));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "5X5H");
		REQUIRE(result.position == 1);
		REQUIRE(result.stop_position == 6);
	}

	SECTION("Mixed =/X trim both") {
		// 5=5X at pos 1 covers [1,11). Region [3,8).
		// 5= at [1,6): trim 2 from left → 3=. 5X at [6,11): trim 3 from right → 2X.
		// Wait: region [3,8). 5= covers [1,6). In region: [3,6) = 3=, before: 2.
		// 5X covers [6,11). In region: [6,8) = 2X, after: 3.
		AlignmentSlicer slicer(3, 8, false);
		auto result = slicer.Slice(make_input("5=5X", 1, 11));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "2H3=2X3H");
		REQUIRE(result.position == 3);
		REQUIRE(result.stop_position == 8);
	}
}

TEST_CASE("AlignmentSlicer trim - I operations", "[alignment_slicer][trim]") {
	SECTION("I within region kept") {
		// 5M3I5M at pos 1. 5M covers [1,6), I at ref 6, 5M covers [6,11).
		// Region [4,100). Trim 3 from left of first M.
		// Result: 3H2M3I5M, pos=4, stop=11
		AlignmentSlicer slicer(4, 100, false);
		auto result = slicer.Slice(make_input("5M3I5M", 1, 11));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "3H2M3I5M");
		REQUIRE(result.position == 4);
		REQUIRE(result.stop_position == 11);
	}

	SECTION("I at trim boundary") {
		// 5M2I5M at pos 1. 5M covers [1,6), I at ref 6, 5M covers [6,11).
		// Region [6,100). First 5M fully before. I at ref 6 — in region. Keep I and second 5M.
		// Left trim: 5 query bases from first M.
		// Result: 5H2I5M, pos=6, stop=11
		AlignmentSlicer slicer(6, 100, false);
		auto result = slicer.Slice(make_input("5M2I5M", 1, 11));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "5H2I5M");
		REQUIRE(result.position == 6);
		REQUIRE(result.stop_position == 11);
	}

	SECTION("I before region trimmed") {
		// 2M3I5M at pos 1. 2M covers [1,3), I at ref 3, 5M covers [3,8).
		// Region [3,8). 2M fully before → left_query_trim=2. I at ref 3 — in region [3,8)? Yes.
		// Keep I, keep 5M.
		// Result: 2H3I5M, pos=3, stop=8
		AlignmentSlicer slicer(3, 8, false);
		auto result = slicer.Slice(make_input("2M3I5M", 1, 8));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "2H3I5M");
		REQUIRE(result.position == 3);
		REQUIRE(result.stop_position == 8);
	}

	SECTION("I after region trimmed") {
		// 5M2I3M at pos 1. 5M covers [1,6), I at ref 6, 3M covers [6,9).
		// Region [1,6). 5M fully in. I at ref 6 — not in [1,6). Trimmed.
		// 3M fully after. right_query_trim = 2(I) + 3(M) = 5.
		// Result: 5M5H, pos=1, stop=6
		AlignmentSlicer slicer(1, 6, false);
		auto result = slicer.Slice(make_input("5M2I3M", 1, 9));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "5M5H");
		REQUIRE(result.position == 1);
		REQUIRE(result.stop_position == 6);
	}
}

TEST_CASE("AlignmentSlicer trim - D operations", "[alignment_slicer][trim]") {
	SECTION("D fully within") {
		// 5M3D5M at pos 1. 5M [1,6), 3D [6,9), 5M [9,14). Region [1,14).
		AlignmentSlicer slicer(1, 14, true);
		auto result = slicer.Slice(make_input("5M3D5M", 1, 14));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "5M3D5M");
		REQUIRE(result.position == 1);
		REQUIRE(result.stop_position == 14);
	}

	SECTION("D spans left boundary") {
		// 5M3D5M at pos 1. 5M [1,6), 3D [6,9), 5M [9,14).
		// Region [7,14). 5M before. 3D: before=1 (pos 6), within=2 (pos 7,8). 5M in region.
		// left_query_trim = 5 (from first M).
		// Result: 5H2D5M, pos=7, stop=14
		AlignmentSlicer slicer(7, 14, true);
		auto result = slicer.Slice(make_input("5M3D5M", 1, 14));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "5H2D5M");
		REQUIRE(result.position == 7);
		REQUIRE(result.stop_position == 14);
	}

	SECTION("D spans right boundary") {
		// 5M3D5M at pos 1. 5M [1,6), 3D [6,9), 5M [9,14).
		// Region [1,7). 5M in region. 3D: within=1 (pos 6), after=2 (pos 7,8). 5M after.
		// right_query_trim = 5 (from last M).
		// Result: 5M1D5H, pos=1, stop=7
		AlignmentSlicer slicer(1, 7, true);
		auto result = slicer.Slice(make_input("5M3D5M", 1, 14));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "5M1D5H");
		REQUIRE(result.position == 1);
		REQUIRE(result.stop_position == 7);
	}

	SECTION("D trimmed both sides") {
		// 3M10D3M at pos 1. 3M [1,4), 10D [4,14), 3M [14,17).
		// Region [5,12). 3M before → left_query_trim=3. 10D: before=1, within=7, after=2. 3M after →
		// right_query_trim=3. Result: 3H7D3H, pos=5, stop=12
		AlignmentSlicer slicer(5, 12, true);
		auto result = slicer.Slice(make_input("3M10D3M", 1, 17));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "3H7D3H");
		REQUIRE(result.position == 5);
		REQUIRE(result.stop_position == 12);
	}
}

TEST_CASE("AlignmentSlicer trim - N operations", "[alignment_slicer][trim]") {
	SECTION("N fully within") {
		// 5M100N5M at pos 1. 5M [1,6), 100N [6,106), 5M [106,111).
		// Region [1,111). All within.
		AlignmentSlicer slicer(1, 111, false);
		auto result = slicer.Slice(make_input("5M100N5M", 1, 111));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "5M100N5M");
		REQUIRE(result.position == 1);
		REQUIRE(result.stop_position == 111);
	}

	SECTION("N spans left boundary") {
		// 5M100N5M at pos 1. 5M [1,6), 100N [6,106), 5M [106,111).
		// Region [50,111). 5M before → left_query_trim=5. 100N: before=44 (6..49), within=56 (50..105).
		// 5M in region.
		// Result: 5H56N5M, pos=50, stop=111
		AlignmentSlicer slicer(50, 111, false);
		auto result = slicer.Slice(make_input("5M100N5M", 1, 111));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "5H56N5M");
		REQUIRE(result.position == 50);
		REQUIRE(result.stop_position == 111);
	}

	SECTION("N spans right boundary") {
		// 5M100N5M at pos 1. 5M [1,6), 100N [6,106), 5M [106,111).
		// Region [1,50). 5M in region. 100N: within=44 (6..49), after=56 (50..105). 5M after → right_query_trim=5.
		// Result: 5M44N5H, pos=1, stop=50
		AlignmentSlicer slicer(1, 50, false);
		auto result = slicer.Slice(make_input("5M100N5M", 1, 111));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "5M44N5H");
		REQUIRE(result.position == 1);
		REQUIRE(result.stop_position == 50);
	}

	SECTION("Region inside N only") {
		// 5M100N5M at pos 1. Region [20,80). 5M before, N spans, 5M after.
		// N-only overlap → not overlapping (N never counts).
		AlignmentSlicer slicer(20, 80, false);
		auto result = slicer.Slice(make_input("5M100N5M", 1, 111));
		REQUIRE_FALSE(result.overlaps);
	}
}

TEST_CASE("AlignmentSlicer trim - S operations", "[alignment_slicer][trim]") {
	SECTION("Leading S, M fully in region") {
		// 3S5M at pos 4. S doesn't consume ref. 5M covers [4,9).
		// Region [4,9). M fully in. S kept (it's at the start, M is in region).
		AlignmentSlicer slicer(4, 9, false);
		auto result = slicer.Slice(make_input("3S5M", 4, 9));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "3S5M");
		REQUIRE(result.position == 4);
		REQUIRE(result.stop_position == 9);
	}

	SECTION("Leading S, trim into M") {
		// 3S5M at pos 4. 5M covers [4,9). Region [6,9).
		// 2M trimmed from left → left_query_trim = 3(S) + 2(M) = 5. 3M remains.
		// Result: 5H3M, pos=6, stop=9
		AlignmentSlicer slicer(6, 9, false);
		auto result = slicer.Slice(make_input("3S5M", 4, 9));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "5H3M");
		REQUIRE(result.position == 6);
		REQUIRE(result.stop_position == 9);
	}

	SECTION("Trailing S, M fully in region") {
		// 5M3S at pos 1. 5M covers [1,6). S doesn't consume ref.
		// Region [1,6). M fully in. S kept (trailing, M was all in region).
		AlignmentSlicer slicer(1, 6, false);
		auto result = slicer.Slice(make_input("5M3S", 1, 6));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "5M3S");
		REQUIRE(result.position == 1);
		REQUIRE(result.stop_position == 6);
	}

	SECTION("Trailing S, trim from M") {
		// 5M3S at pos 1. 5M covers [1,6). Region [1,4).
		// 3M kept, 2M trimmed right. S after trimmed M → right_query_trim = 2(M) + 3(S) = 5.
		// Result: 3M5H, pos=1, stop=4
		AlignmentSlicer slicer(1, 4, false);
		auto result = slicer.Slice(make_input("5M3S", 1, 6));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "3M5H");
		REQUIRE(result.position == 1);
		REQUIRE(result.stop_position == 4);
	}
}

TEST_CASE("AlignmentSlicer trim - H operations (existing)", "[alignment_slicer][trim]") {
	SECTION("Existing H + additional trim") {
		// 3H10M2H at pos 5. 10M covers [5,15).
		// Region [8,12). Trim 3M left, 3M right. Plus existing 3H left, 2H right.
		// left H = 3(existing) + 3(trimmed) = 6. right H = 3(trimmed) + 2(existing) = 5.
		// Wait: 10M at [5,15). Region [8,12). Before: 3 (pos 5,6,7). After: 3 (pos 12,13,14).
		// Within: 4 (pos 8,9,10,11). left_query_trim=3. right_query_trim=3.
		// left H = 3 + 3 = 6. right H = 3 + 2 = 5.
		// Result: 6H4M5H, pos=8, stop=12
		AlignmentSlicer slicer(8, 12, false);
		auto result = slicer.Slice(make_input("3H10M2H", 5, 15));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "6H4M5H");
		REQUIRE(result.position == 8);
		REQUIRE(result.stop_position == 12);
	}

	SECTION("Existing H, no additional trim") {
		// 3H10M2H at pos 5. Region [1,100). Fully within.
		AlignmentSlicer slicer(1, 100, false);
		auto result = slicer.Slice(make_input("3H10M2H", 5, 15));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "3H10M2H");
		REQUIRE(result.position == 5);
		REQUIRE(result.stop_position == 15);
	}
}

TEST_CASE("AlignmentSlicer trim - P operation", "[alignment_slicer][trim]") {
	// 5M1P5M at pos 1. P is padding, doesn't consume ref or query.
	// 5M [1,6), P, 5M [6,11). Region [1,11). All within.
	AlignmentSlicer slicer(1, 11, false);
	auto result = slicer.Slice(make_input("5M1P5M", 1, 11));
	REQUIRE(result.overlaps);
	REQUIRE(result.cigar == "5M1P5M");
	REQUIRE(result.position == 1);
	REQUIRE(result.stop_position == 11);
}

TEST_CASE("AlignmentSlicer trim - complex multi-op", "[alignment_slicer][trim]") {
	SECTION("Multi-op left trim") {
		// 3M2I5M3D4M at pos 1.
		// 3M covers [1,4) (query 0-2). 2I at ref 4 (query 3-4). 5M covers [4,9) (query 5-9).
		// 3D covers [9,12). 4M covers [12,16) (query 10-13).
		// Region [5,100). 3M before (left_query_trim+=3). I at ref 4 — before region [5,100). left_query_trim+=2.
		// 5M at [4,9): split. Before: 1 (pos 4). Within: 4 (pos 5-8). left_query_trim+=1.
		// Total left_query_trim = 3+2+1 = 6. 4M within. 3D within. 4M within.
		// Result: 6H4M3D4M, pos=5, stop=16
		AlignmentSlicer slicer(5, 100, false);
		auto result = slicer.Slice(make_input("3M2I5M3D4M", 1, 16));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "6H4M3D4M");
		REQUIRE(result.position == 5);
		REQUIRE(result.stop_position == 16);
	}

	SECTION("Full complex trim with all op types") {
		// 2H3S5M2I3M3D4M10N2M3S2H at pos 10
		// Walk:
		// 2H: existing left H=2, no ref/query consumed
		// 3S: query 0-2, no ref consumed
		// 5M: ref [10,15), query 3-7
		// 2I: ref stays at 15, query 8-9
		// 3M: ref [15,18), query 10-12
		// 3D: ref [18,21), no query
		// 4M: ref [21,25), query 13-16
		// 10N: ref [25,35), no query
		// 2M: ref [35,37), query 17-18
		// 3S: query 19-21, no ref consumed
		// 2H: existing right H=2, no ref/query consumed
		//
		// Region [12,25).
		// 2H: left H existing = 2
		// 3S: before region (we haven't entered yet) → left_query_trim += 3
		// 5M [10,15): before=2 (10,11), within=3 (12,13,14). left_query_trim += 2. Keep 3M.
		// 2I at ref 15: in [12,25)? Yes. Keep 2I.
		// 3M [15,18): fully within. Keep 3M.
		// 3D [18,21): fully within. Keep 3D.
		// 4M [21,25): fully within [12,25)? [21,25) is within [12,25). Yes. Keep 4M.
		// 10N [25,35): starts at 25 which is >= 25 (region stop). After.
		// 2M [35,37): after → right_query_trim += 2
		// 3S: after → right_query_trim += 3
		// 2H: existing right H = 2
		//
		// left_query_trim = 3(S) + 2(M partial) = 5
		// right_query_trim = 2(M) + 3(S) = 5
		// left H = 2(existing) + 5(trimmed) = 7
		// right H = 5(trimmed) + 2(existing) = 7
		// CIGAR: 7H3M2I3M3D4M7H
		// pos = max(10, 12) = 12
		// stop = min(37, 25) = 25
		AlignmentSlicer slicer(12, 25, false);
		auto result = slicer.Slice(make_input("2H3S5M2I3M3D4M10N2M3S2H", 10, 37));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "7H3M2I3M3D4M7H");
		REQUIRE(result.position == 12);
		REQUIRE(result.stop_position == 25);
	}
}

// =====================================================================
// 1c. Sequence/quality trimming tests
// =====================================================================

TEST_CASE("AlignmentSlicer seq/qual trim - M operations", "[alignment_slicer][seqqual]") {
	SECTION("No trim") {
		AlignmentSlicer slicer(1, 6, false);
		auto result = slicer.Slice(make_input_sq("5M", 1, 6, "ACGTG", {30, 31, 32, 33, 34}));
		REQUIRE(result.overlaps);
		REQUIRE(result.sequence == "ACGTG");
		REQUIRE(result.quality == std::vector<uint8_t>({30, 31, 32, 33, 34}));
	}

	SECTION("Left trim 2") {
		AlignmentSlicer slicer(3, 6, false);
		auto result = slicer.Slice(make_input_sq("5M", 1, 6, "ACGTG", {30, 31, 32, 33, 34}));
		REQUIRE(result.overlaps);
		REQUIRE(result.sequence == "GTG");
		REQUIRE(result.quality == std::vector<uint8_t>({32, 33, 34}));
	}

	SECTION("Right trim 2") {
		AlignmentSlicer slicer(1, 4, false);
		auto result = slicer.Slice(make_input_sq("5M", 1, 6, "ACGTG", {30, 31, 32, 33, 34}));
		REQUIRE(result.overlaps);
		REQUIRE(result.sequence == "ACG");
		REQUIRE(result.quality == std::vector<uint8_t>({30, 31, 32}));
	}

	SECTION("Both trim") {
		AlignmentSlicer slicer(2, 4, false);
		auto result = slicer.Slice(make_input_sq("5M", 1, 6, "ACGTG", {30, 31, 32, 33, 34}));
		REQUIRE(result.overlaps);
		REQUIRE(result.sequence == "CG");
		REQUIRE(result.quality == std::vector<uint8_t>({31, 32}));
	}
}

TEST_CASE("AlignmentSlicer seq/qual trim - = and X operations", "[alignment_slicer][seqqual]") {
	SECTION("= op trim") {
		AlignmentSlicer slicer(3, 6, false);
		auto result = slicer.Slice(make_input_sq("5=", 1, 6, "ACGTG", {30, 31, 32, 33, 34}));
		REQUIRE(result.overlaps);
		REQUIRE(result.sequence == "GTG");
		REQUIRE(result.quality == std::vector<uint8_t>({32, 33, 34}));
	}

	SECTION("X op trim") {
		AlignmentSlicer slicer(3, 6, false);
		auto result = slicer.Slice(make_input_sq("5X", 1, 6, "ACGTG", {30, 31, 32, 33, 34}));
		REQUIRE(result.overlaps);
		REQUIRE(result.sequence == "GTG");
		REQUIRE(result.quality == std::vector<uint8_t>({32, 33, 34}));
	}
}

TEST_CASE("AlignmentSlicer seq/qual trim - I operations", "[alignment_slicer][seqqual]") {
	SECTION("I kept in region") {
		// 3M2I3M at pos 1. 3M [1,4) query 0-2. 2I at ref 4, query 3-4. 3M [4,7) query 5-7.
		// Region [2,7). Trim 1M left → left_query_trim=1.
		// Seq: ABCDEFGH (8 chars: 3M + 2I + 3M)
		AlignmentSlicer slicer(2, 7, false);
		auto result = slicer.Slice(make_input_sq("3M2I3M", 1, 7, "ABCDEFGH", {30, 31, 32, 33, 34, 35, 36, 37}));
		REQUIRE(result.overlaps);
		REQUIRE(result.sequence == "BCDEFGH");
		REQUIRE(result.quality == std::vector<uint8_t>({31, 32, 33, 34, 35, 36, 37}));
	}

	SECTION("I trimmed before region") {
		// 2M3I5M at pos 1. 2M [1,3) query 0-1. 3I at ref 3, query 2-4. 5M [3,8) query 5-9.
		// Region [3,8). 2M before → left_query_trim=2. I at ref 3 — in region.
		// Seq: ABCDEFGHIJ (10 chars: 2M + 3I + 5M). Trim 2 from left.
		AlignmentSlicer slicer(3, 8, false);
		auto result =
		    slicer.Slice(make_input_sq("2M3I5M", 1, 8, "ABCDEFGHIJ", {30, 31, 32, 33, 34, 35, 36, 37, 38, 39}));
		REQUIRE(result.overlaps);
		REQUIRE(result.sequence == "CDEFGHIJ");
		REQUIRE(result.quality == std::vector<uint8_t>({32, 33, 34, 35, 36, 37, 38, 39}));
	}
}

TEST_CASE("AlignmentSlicer seq/qual trim - D doesn't affect seq", "[alignment_slicer][seqqual]") {
	// 3M3D3M at pos 1. 3M [1,4), 3D [4,7), 3M [7,10). Seq: ABCDEF (6 chars, D has no query bases).
	// Region [1,10). All within.
	AlignmentSlicer slicer(1, 10, true);
	auto result = slicer.Slice(make_input_sq("3M3D3M", 1, 10, "ABCDEF", {30, 31, 32, 33, 34, 35}));
	REQUIRE(result.overlaps);
	REQUIRE(result.sequence == "ABCDEF");
	REQUIRE(result.quality == std::vector<uint8_t>({30, 31, 32, 33, 34, 35}));
}

TEST_CASE("AlignmentSlicer seq/qual trim - N doesn't affect seq", "[alignment_slicer][seqqual]") {
	// 3M5N3M at pos 1. 3M [1,4), 5N [4,9), 3M [9,12). Seq: ABCDEF (6 chars).
	// Region [1,12). All within.
	AlignmentSlicer slicer(1, 12, false);
	auto result = slicer.Slice(make_input_sq("3M5N3M", 1, 12, "ABCDEF", {30, 31, 32, 33, 34, 35}));
	REQUIRE(result.overlaps);
	REQUIRE(result.sequence == "ABCDEF");
	REQUIRE(result.quality == std::vector<uint8_t>({30, 31, 32, 33, 34, 35}));
}

TEST_CASE("AlignmentSlicer seq/qual trim - S operations", "[alignment_slicer][seqqual]") {
	SECTION("Leading S trimmed with partial M") {
		// 3S5M at pos 4. Seq: ABCDEFGH (8 chars: 3S + 5M).
		// Region [6,9). S is before, 2M trimmed. left_query_trim = 3(S) + 2(M) = 5.
		AlignmentSlicer slicer(6, 9, false);
		auto result = slicer.Slice(make_input_sq("3S5M", 4, 9, "ABCDEFGH", {30, 31, 32, 33, 34, 35, 36, 37}));
		REQUIRE(result.overlaps);
		REQUIRE(result.sequence == "FGH");
		REQUIRE(result.quality == std::vector<uint8_t>({35, 36, 37}));
	}
}

TEST_CASE("AlignmentSlicer seq/qual trim - H doesn't affect seq", "[alignment_slicer][seqqual]") {
	// 3H5M2H at pos 5. H not in seq. Seq: ABCDE (5 chars, only M).
	// Region [7,10). Trim 2M left.
	AlignmentSlicer slicer(7, 10, false);
	auto result = slicer.Slice(make_input_sq("3H5M2H", 5, 10, "ABCDE", {30, 31, 32, 33, 34}));
	REQUIRE(result.overlaps);
	REQUIRE(result.sequence == "CDE");
	REQUIRE(result.quality == std::vector<uint8_t>({32, 33, 34}));
}

TEST_CASE("AlignmentSlicer seq/qual trim - empty seq passthrough", "[alignment_slicer][seqqual]") {
	AlignmentSlicer slicer(3, 6, false);
	auto result = slicer.Slice(make_input_sq("5M", 1, 6, "", {}));
	REQUIRE(result.overlaps);
	REQUIRE(result.sequence.empty());
	REQUIRE(result.quality.empty());
}

TEST_CASE("AlignmentSlicer seq/qual trim - complex with all ops", "[alignment_slicer][seqqual]") {
	// 2H3S5M2I3M3D4M10N2M3S2H at pos 10
	// Query-consuming ops: S(3) M(5) I(2) M(3) M(4) M(2) S(3) = 22 query bases
	// H doesn't appear in sequence.
	// Region [12,25). From the trim test: left_query_trim=5 (3S + 2M), right_query_trim=5 (2M + 3S)
	// Seq has 22 chars. After trimming: 22 - 5 - 5 = 12 chars.
	std::string seq = "ABCDEFGHIJKLMNOPQRSTUV"; // 22 chars
	std::vector<uint8_t> qual;
	for (uint8_t i = 0; i < 22; i++) {
		qual.push_back(30 + i);
	}

	AlignmentSlicer slicer(12, 25, false);
	auto result = slicer.Slice(make_input_sq("2H3S5M2I3M3D4M10N2M3S2H", 10, 37, seq, qual));
	REQUIRE(result.overlaps);
	REQUIRE(result.sequence == "FGHIJKLMNOPQ"); // chars 5..16 (0-indexed)
	REQUIRE(result.quality.size() == 12);
	REQUIRE(result.quality[0] == 35);  // 30+5
	REQUIRE(result.quality[11] == 46); // 30+16
}

// =====================================================================
// 1d. Edge case tests
// =====================================================================

TEST_CASE("AlignmentSlicer edge cases", "[alignment_slicer][edge]") {
	SECTION("Region start > stop throws") {
		REQUIRE_THROWS_AS(AlignmentSlicer(10, 5, false), miint::InvalidInputException);
	}

	SECTION("Region start < 1 throws") {
		REQUIRE_THROWS_AS(AlignmentSlicer(0, 10, false), miint::InvalidInputException);
	}

	SECTION("Empty region (start == stop) - no reads overlap") {
		AlignmentSlicer slicer(5, 5, false);
		auto result = slicer.Slice(make_input("10M", 1, 11));
		REQUIRE_FALSE(result.overlaps);
	}

	SECTION("Position < 1 (unmapped) returns no overlap") {
		AlignmentSlicer slicer(1, 100, false);
		auto result = slicer.Slice(make_input("10M", 0, 10));
		REQUIRE_FALSE(result.overlaps);
	}

	SECTION("Negative position returns no overlap") {
		AlignmentSlicer slicer(1, 100, false);
		auto result = slicer.Slice(make_input("10M", -1, 9));
		REQUIRE_FALSE(result.overlaps);
	}
}

// =====================================================================
// Additional coverage: reviewer-identified gaps
// =====================================================================

TEST_CASE("AlignmentSlicer - pure insertion-only overlap", "[alignment_slicer][overlap]") {
	// Policy: an insertion at a ref position within the region counts as overlap.
	// Insertions don't consume reference but represent query bases inserted at
	// a specific reference coordinate. If that coordinate is in the region,
	// the insertion is evidence of the query interacting with the region.
	SECTION("I-only at region boundary - M before region") {
		// 5M2I at pos 1. 5M covers [1,6). I at ref 6. Region [6,10).
		// 5M fully before region. I at ref 6 is in [6,10). M is NOT in region.
		// Only the I is in the region. Does this overlap?
		// Yes — the insertion is at a reference coordinate within the region.
		AlignmentSlicer slicer(6, 10, false);
		auto result = slicer.Slice(make_input("5M2I", 1, 6));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "5H2I");
		REQUIRE(result.position == 6);
		REQUIRE(result.stop_position == 6);
	}

	SECTION("I-only NOT in region") {
		// 5M2I at pos 1. I at ref 6. Region [7,10). I not in region.
		AlignmentSlicer slicer(7, 10, false);
		auto result = slicer.Slice(make_input("5M2I", 1, 6));
		REQUIRE_FALSE(result.overlaps);
	}
}

TEST_CASE("AlignmentSlicer - sequence/quality length mismatch", "[alignment_slicer][edge]") {
	// When sequence length doesn't match CIGAR-implied query length, the trim
	// indices may exceed sequence bounds. In this case, sequence/quality are
	// returned empty (not crashed). This is a defensive behavior for corrupt input.
	SECTION("Sequence shorter than CIGAR implies") {
		// 10M implies 10 query bases, but only 3 provided
		AlignmentSlicer slicer(4, 8, false);
		auto result = slicer.Slice(make_input_sq("10M", 1, 11, "ABC", {30, 31, 32}));
		REQUIRE(result.overlaps);
		REQUIRE(result.cigar == "3H4M3H");
		// Sequence/quality should be empty (trim exceeds bounds)
		REQUIRE(result.sequence.empty());
		REQUIRE(result.quality.empty());
	}

	SECTION("Sequence matches CIGAR - normal case") {
		AlignmentSlicer slicer(4, 8, false);
		auto result = slicer.Slice(make_input_sq("10M", 1, 11, "ABCDEFGHIJ", {30, 31, 32, 33, 34, 35, 36, 37, 38, 39}));
		REQUIRE(result.overlaps);
		REQUIRE(result.sequence == "DEFG");
		REQUIRE(result.quality == std::vector<uint8_t>({33, 34, 35, 36}));
	}
}

TEST_CASE("AlignmentSlicer - leading S with I as first non-S op", "[alignment_slicer][overlap]") {
	// 3S2I5M at pos 5. S(3) pending. I(2) at ref 5. M(5) at [5,10).
	// Region [5,10). I at ref 5 is in region. S should be kept.
	AlignmentSlicer slicer(5, 10, false);
	auto result = slicer.Slice(make_input_sq("3S2I5M", 5, 10, "ABCDEFGHIJ", {30, 31, 32, 33, 34, 35, 36, 37, 38, 39}));
	REQUIRE(result.overlaps);
	REQUIRE(result.cigar == "3S2I5M");
	REQUIRE(result.position == 5);
	REQUIRE(result.stop_position == 10);
	REQUIRE(result.sequence == "ABCDEFGHIJ");
}
