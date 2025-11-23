#include "../../src/include/alignment_functions_internal.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("ParseCigar - Basic operations", "[alignment_functions]") {
	SECTION("Simple match") {
		auto stats = miint::ParseCigar("10M");
		REQUIRE((stats.matches == 10));
		REQUIRE((stats.alignment_columns == 10));
		REQUIRE((stats.gap_opens == 0));
	}

	SECTION("Match with explicit operations") {
		auto stats = miint::ParseCigar("5=2X3=");
		REQUIRE((stats.matches == 10));
		REQUIRE((stats.match_ops == 8));
		REQUIRE((stats.mismatch_ops == 2));
		REQUIRE((stats.alignment_columns == 10));
	}

	SECTION("Insertion") {
		auto stats = miint::ParseCigar("10M5I10M");
		REQUIRE((stats.matches == 20));
		REQUIRE((stats.insertions == 5));
		REQUIRE((stats.alignment_columns == 25));
		REQUIRE((stats.gap_opens == 1));
	}

	SECTION("Deletion") {
		auto stats = miint::ParseCigar("10M3D10M");
		REQUIRE((stats.matches == 20));
		REQUIRE((stats.deletions == 3));
		REQUIRE((stats.alignment_columns == 23));
		REQUIRE((stats.gap_opens == 1));
	}

	SECTION("Multiple gap opens") {
		auto stats = miint::ParseCigar("10M2I5M3D5M");
		REQUIRE((stats.matches == 20));
		REQUIRE((stats.insertions == 2));
		REQUIRE((stats.deletions == 3));
		REQUIRE((stats.gap_opens == 2));
		REQUIRE((stats.alignment_columns == 25));
	}

	SECTION("Consecutive insertions count as one gap open") {
		auto stats = miint::ParseCigar("10M2I3I5M");
		REQUIRE((stats.insertions == 5));
		REQUIRE((stats.gap_opens == 1));
	}

	SECTION("Consecutive deletions count as one gap open") {
		auto stats = miint::ParseCigar("10M2D3D5M");
		REQUIRE((stats.deletions == 5));
		REQUIRE((stats.gap_opens == 1));
	}
}

TEST_CASE("ParseCigar - Clipping and skipping", "[alignment_functions]") {
	SECTION("Soft clipping tracked") {
		auto stats = miint::ParseCigar("5S10M5S");
		REQUIRE((stats.matches == 10));
		REQUIRE((stats.alignment_columns == 10));
		REQUIRE((stats.soft_clips == 10));
		REQUIRE((stats.hard_clips == 0));
	}

	SECTION("Hard clipping tracked") {
		auto stats = miint::ParseCigar("5H10M5H");
		REQUIRE((stats.matches == 10));
		REQUIRE((stats.alignment_columns == 10));
		REQUIRE((stats.soft_clips == 0));
		REQUIRE((stats.hard_clips == 10));
	}

	SECTION("Reference skip (N) ignored") {
		auto stats = miint::ParseCigar("10M100N10M");
		REQUIRE((stats.matches == 20));
		REQUIRE((stats.alignment_columns == 20));
	}

	SECTION("Complex CIGAR with all operations") {
		auto stats = miint::ParseCigar("5H5S10M2I5M3D5M5S5H");
		REQUIRE((stats.matches == 20));
		REQUIRE((stats.insertions == 2));
		REQUIRE((stats.deletions == 3));
		REQUIRE((stats.gap_opens == 2));
		REQUIRE((stats.alignment_columns == 25));
		REQUIRE((stats.soft_clips == 10));
		REQUIRE((stats.hard_clips == 10));
	}

	SECTION("Only soft clipping") {
		auto stats = miint::ParseCigar("100S");
		REQUIRE((stats.matches == 0));
		REQUIRE((stats.alignment_columns == 0));
		REQUIRE((stats.soft_clips == 100));
		REQUIRE((stats.hard_clips == 0));
	}

	SECTION("Only hard clipping") {
		auto stats = miint::ParseCigar("100H");
		REQUIRE((stats.matches == 0));
		REQUIRE((stats.alignment_columns == 0));
		REQUIRE((stats.soft_clips == 0));
		REQUIRE((stats.hard_clips == 100));
	}

	SECTION("Mixed soft and hard clipping") {
		auto stats = miint::ParseCigar("10H20S10M20S10H");
		REQUIRE((stats.matches == 10));
		REQUIRE((stats.alignment_columns == 10));
		REQUIRE((stats.soft_clips == 40));
		REQUIRE((stats.hard_clips == 20));
	}

	SECTION("Complex CIGAR from Heng Li's blog post") {
		// The blog post states this CIGAR string has 43 matches and one mismatch,
		// however, it appears the string was produced without the use of -xeq
		// and thus does not differentiate =/X
		auto stats = miint::ParseCigar("18M3D2M2D2M1I22M");
		REQUIRE((stats.matches == 44));
		REQUIRE((stats.insertions == 1));
		REQUIRE((stats.deletions == 5));
		REQUIRE((stats.gap_opens == 3));
		REQUIRE((stats.alignment_columns == 50));
	}
}

TEST_CASE("ParseCigar - Edge cases", "[alignment_functions]") {
	SECTION("Empty string") {
		auto stats = miint::ParseCigar("");
		REQUIRE((stats.matches == 0));
		REQUIRE((stats.alignment_columns == 0));
	}

	SECTION("Unmapped (*)") {
		auto stats = miint::ParseCigar("*");
		REQUIRE((stats.matches == 0));
		REQUIRE((stats.alignment_columns == 0));
	}

	SECTION("Large numbers") {
		auto stats = miint::ParseCigar("150M");
		REQUIRE((stats.matches == 150));
		REQUIRE((stats.alignment_columns == 150));
	}
}

TEST_CASE("ParseCigar - Error handling", "[alignment_functions]") {
	SECTION("Invalid operation without length") {
		REQUIRE_THROWS_AS(miint::ParseCigar("M"), miint::InvalidInputException);
	}

	SECTION("Invalid operation character") {
		REQUIRE_THROWS_AS(miint::ParseCigar("10Z"), miint::InvalidInputException);
	}
}

TEST_CASE("ParseMd - Basic operations", "[alignment_functions]") {
	SECTION("All matches") {
		auto stats = miint::ParseMd("10");
		REQUIRE((stats.matches == 10));
		REQUIRE((stats.mismatches == 0));
	}

	SECTION("Single mismatch") {
		auto stats = miint::ParseMd("5A4");
		REQUIRE((stats.matches == 9));
		REQUIRE((stats.mismatches == 1));
	}

	SECTION("Multiple mismatches") {
		auto stats = miint::ParseMd("3A2T3");
		REQUIRE((stats.matches == 8));
		REQUIRE((stats.mismatches == 2));
	}

	SECTION("Deletion marker") {
		auto stats = miint::ParseMd("5^AC4");
		REQUIRE((stats.matches == 9));
		REQUIRE((stats.mismatches == 0));
	}

	SECTION("Multiple deletions") {
		auto stats = miint::ParseMd("3^A2^TG4");
		REQUIRE((stats.matches == 9));
		REQUIRE((stats.mismatches == 0));
	}

	SECTION("Mixed operations") {
		auto stats = miint::ParseMd("3A2^TG3C1");
		REQUIRE((stats.matches == 9));
		REQUIRE((stats.mismatches == 2));
	}
}

TEST_CASE("ParseMd - Edge cases", "[alignment_functions]") {
	SECTION("Empty string") {
		auto stats = miint::ParseMd("");
		REQUIRE((stats.matches == 0));
		REQUIRE((stats.mismatches == 0));
	}

	SECTION("Zero at start") {
		auto stats = miint::ParseMd("0A10");
		REQUIRE((stats.matches == 10));
		REQUIRE((stats.mismatches == 1));
	}

	SECTION("Zero at end") {
		auto stats = miint::ParseMd("10A0");
		REQUIRE((stats.matches == 10));
		REQUIRE((stats.mismatches == 1));
	}

	SECTION("Only mismatches") {
		auto stats = miint::ParseMd("0A0T0C0");
		REQUIRE((stats.matches == 0));
		REQUIRE((stats.mismatches == 3));
	}

	SECTION("Large numbers") {
		auto stats = miint::ParseMd("100A50");
		REQUIRE((stats.matches == 150));
		REQUIRE((stats.mismatches == 1));
	}
}

TEST_CASE("ParseMd - Complex patterns", "[alignment_functions]") {
	SECTION("Consecutive mismatches") {
		auto stats = miint::ParseMd("5AG3");
		REQUIRE((stats.matches == 8));
		REQUIRE((stats.mismatches == 2));
	}

	SECTION("Long deletion") {
		auto stats = miint::ParseMd("10^ACGTACGT10");
		REQUIRE((stats.matches == 20));
		REQUIRE((stats.mismatches == 0));
	}

	SECTION("All three types") {
		auto stats = miint::ParseMd("5A3^TG2C3");
		REQUIRE((stats.matches == 13));
		REQUIRE((stats.mismatches == 2));
	}

	SECTION("MD tag ending with deletion") {
		auto stats = miint::ParseMd("10^AC");
		REQUIRE((stats.matches == 10));
		REQUIRE((stats.mismatches == 0));
	}

	SECTION("MD tag with only deletions") {
		auto stats = miint::ParseMd("^AC^TG");
		REQUIRE((stats.matches == 0));
		REQUIRE((stats.mismatches == 0));
	}

	SECTION("MD tag with many consecutive mismatches") {
		auto stats = miint::ParseMd("5ACGT4");
		REQUIRE((stats.matches == 9));
		REQUIRE((stats.mismatches == 4));
	}

	SECTION("Very long deletion sequence") {
		auto stats = miint::ParseMd("10^ACGTACGTACGTACGT10");
		REQUIRE((stats.matches == 20));
		REQUIRE((stats.mismatches == 0));
	}
}

TEST_CASE("ComputeQueryLength - Basic operations", "[alignment_functions]") {
	SECTION("Simple match - no clipping") {
		auto stats = miint::ParseCigar("10M");
		REQUIRE((miint::ComputeQueryLength(stats, true) == 10));
		REQUIRE((miint::ComputeQueryLength(stats, false) == 10));
	}

	SECTION("Match with insertions") {
		auto stats = miint::ParseCigar("10M5I");
		REQUIRE((miint::ComputeQueryLength(stats, true) == 15));
		REQUIRE((miint::ComputeQueryLength(stats, false) == 15));
	}

	SECTION("Match with deletions - deletions don't consume query") {
		auto stats = miint::ParseCigar("10M5D");
		REQUIRE((miint::ComputeQueryLength(stats, true) == 10));
		REQUIRE((miint::ComputeQueryLength(stats, false) == 10));
	}

	SECTION("With soft clips only") {
		auto stats = miint::ParseCigar("5S10M5S");
		// Soft clips are always included in query length
		REQUIRE((miint::ComputeQueryLength(stats, true) == 20));
		REQUIRE((miint::ComputeQueryLength(stats, false) == 20));
	}

	SECTION("With hard clips only") {
		auto stats = miint::ParseCigar("5H10M5H");
		// Hard clips included only when parameter is true
		REQUIRE((miint::ComputeQueryLength(stats, true) == 20));
		REQUIRE((miint::ComputeQueryLength(stats, false) == 10));
	}

	SECTION("With both soft and hard clips") {
		auto stats = miint::ParseCigar("5H5S10M5S5H");
		// Soft clips always included, hard clips only when parameter is true
		REQUIRE((miint::ComputeQueryLength(stats, true) == 30));
		REQUIRE((miint::ComputeQueryLength(stats, false) == 20));
	}

	SECTION("Complex CIGAR with all operations") {
		auto stats = miint::ParseCigar("5H5S10M5I3D5M5S5H");
		// M=10+5=15, I=5, S=10, H=10
		// Total with H: 15+5+10+10=40
		// Total without H: 15+5+10=30
		REQUIRE((miint::ComputeQueryLength(stats, true) == 40));
		REQUIRE((miint::ComputeQueryLength(stats, false) == 30));
	}
}

TEST_CASE("ComputeQueryLength - Verify matches bam_cigar2qlen behavior", "[alignment_functions]") {
	// When include_hard_clips=false, should match bam_cigar2qlen
	// bam_cigar2qlen counts M, I, S, =, X operations (excludes H, D, N, P)

	SECTION("HTSlib example 1: 10M5I5S") {
		auto stats = miint::ParseCigar("10M5I5S");
		// bam_cigar2qlen would return 20 (10M + 5I + 5S)
		REQUIRE((miint::ComputeQueryLength(stats, false) == 20));
	}

	SECTION("HTSlib example 2: 5S10M3D5M5S") {
		auto stats = miint::ParseCigar("5S10M3D5M5S");
		// bam_cigar2qlen would return 25 (5S + 10M + 5M + 5S, excludes 3D)
		REQUIRE((miint::ComputeQueryLength(stats, false) == 25));
	}

	SECTION("HTSlib example 3: 5H10M5H") {
		auto stats = miint::ParseCigar("5H10M5H");
		// bam_cigar2qlen would return 10 (excludes hard clips)
		REQUIRE((miint::ComputeQueryLength(stats, false) == 10));
	}

	SECTION("HTSlib example 4: 5=2X3= (explicit match/mismatch)") {
		auto stats = miint::ParseCigar("5=2X3=");
		// bam_cigar2qlen would return 10 (= and X both consume query)
		REQUIRE((miint::ComputeQueryLength(stats, false) == 10));
	}

	SECTION("HTSlib example 5: 18M3D2M2D2M1I22M") {
		auto stats = miint::ParseCigar("18M3D2M2D2M1I22M");
		// bam_cigar2qlen: 18M + 2M + 2M + 1I + 22M = 45 (excludes deletions)
		REQUIRE((miint::ComputeQueryLength(stats, false) == 45));
	}
}

TEST_CASE("ComputeQueryLength - Edge cases", "[alignment_functions]") {
	SECTION("Empty CIGAR") {
		auto stats = miint::ParseCigar("");
		REQUIRE((miint::ComputeQueryLength(stats, true) == 0));
		REQUIRE((miint::ComputeQueryLength(stats, false) == 0));
	}

	SECTION("Unmapped read (*)") {
		auto stats = miint::ParseCigar("*");
		REQUIRE((miint::ComputeQueryLength(stats, true) == 0));
		REQUIRE((miint::ComputeQueryLength(stats, false) == 0));
	}

	SECTION("Only soft clipping") {
		auto stats = miint::ParseCigar("100S");
		REQUIRE((miint::ComputeQueryLength(stats, true) == 100));
		REQUIRE((miint::ComputeQueryLength(stats, false) == 100));
	}

	SECTION("Only hard clipping") {
		auto stats = miint::ParseCigar("100H");
		REQUIRE((miint::ComputeQueryLength(stats, true) == 100));
		REQUIRE((miint::ComputeQueryLength(stats, false) == 0));
	}

	SECTION("Reference skip (N) doesn't consume query") {
		auto stats = miint::ParseCigar("10M100N10M");
		REQUIRE((miint::ComputeQueryLength(stats, true) == 20));
		REQUIRE((miint::ComputeQueryLength(stats, false) == 20));
	}
}

TEST_CASE("ComputeQueryCoverage - Aligned type", "[alignment_functions]") {
	using Catch::Matchers::WithinRel;

	SECTION("Perfect alignment - no clipping") {
		auto stats = miint::ParseCigar("10M");
		// 10 aligned / 10 total = 1.0
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "aligned"), WithinRel(1.0, 0.001));
	}

	SECTION("Alignment with insertions") {
		auto stats = miint::ParseCigar("10M5I");
		// 10 aligned / 15 total = 0.6667
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "aligned"), WithinRel(0.6667, 0.001));
	}

	SECTION("Alignment with soft clips") {
		auto stats = miint::ParseCigar("5S10M5S");
		// 10 aligned / 20 total = 0.5
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "aligned"), WithinRel(0.5, 0.001));
	}

	SECTION("Alignment with hard clips") {
		auto stats = miint::ParseCigar("5H10M5H");
		// 10 aligned / 20 total = 0.5
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "aligned"), WithinRel(0.5, 0.001));
	}

	SECTION("Alignment with both soft and hard clips") {
		auto stats = miint::ParseCigar("5H5S10M5S5H");
		// 10 aligned / 30 total = 0.3333
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "aligned"), WithinRel(0.3333, 0.001));
	}

	SECTION("Complex alignment with insertions and clips") {
		auto stats = miint::ParseCigar("5H5S10M5I5S5H");
		// 10 aligned / 35 total = 0.2857
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "aligned"), WithinRel(0.2857, 0.001));
	}

	SECTION("Only soft clips - zero coverage") {
		auto stats = miint::ParseCigar("100S");
		// 0 aligned / 100 total = 0.0
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "aligned"), WithinRel(0.0, 0.001));
	}

	SECTION("Only hard clips - zero coverage") {
		auto stats = miint::ParseCigar("100H");
		// 0 aligned / 100 total = 0.0
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "aligned"), WithinRel(0.0, 0.001));
	}
}

TEST_CASE("ComputeQueryCoverage - Mapped type", "[alignment_functions]") {
	using Catch::Matchers::WithinRel;

	SECTION("Perfect alignment - no clipping") {
		auto stats = miint::ParseCigar("10M");
		// 10 mapped / 10 total = 1.0
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "mapped"), WithinRel(1.0, 0.001));
	}

	SECTION("Alignment with insertions") {
		auto stats = miint::ParseCigar("10M5I");
		// 15 mapped (M+I) / 15 total = 1.0
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "mapped"), WithinRel(1.0, 0.001));
	}

	SECTION("Alignment with soft clips") {
		auto stats = miint::ParseCigar("5S10M5S");
		// 10 mapped / 20 total = 0.5
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "mapped"), WithinRel(0.5, 0.001));
	}

	SECTION("Alignment with hard clips") {
		auto stats = miint::ParseCigar("5H10M5H");
		// 10 mapped / 20 total = 0.5
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "mapped"), WithinRel(0.5, 0.001));
	}

	SECTION("Alignment with insertions and soft clips") {
		auto stats = miint::ParseCigar("5S10M5I5S");
		// 15 mapped (M+I) / 25 total = 0.6
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "mapped"), WithinRel(0.6, 0.001));
	}

	SECTION("Complex alignment with all operations") {
		auto stats = miint::ParseCigar("5H5S10M5I5S5H");
		// 15 mapped (M+I) / 35 total = 0.4286
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "mapped"), WithinRel(0.4286, 0.001));
	}

	SECTION("Only soft clips - zero coverage") {
		auto stats = miint::ParseCigar("100S");
		// 0 mapped / 100 total = 0.0
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "mapped"), WithinRel(0.0, 0.001));
	}

	SECTION("Only hard clips - zero coverage") {
		auto stats = miint::ParseCigar("100H");
		// 0 mapped / 100 total = 0.0
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "mapped"), WithinRel(0.0, 0.001));
	}
}

TEST_CASE("ComputeQueryCoverage - Edge cases", "[alignment_functions]") {
	using Catch::Matchers::WithinRel;

	SECTION("Empty CIGAR - zero coverage") {
		auto stats = miint::ParseCigar("");
		// 0 / 0 = 0.0 (special case to avoid division by zero)
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "aligned"), WithinRel(0.0, 0.001));
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "mapped"), WithinRel(0.0, 0.001));
	}

	SECTION("Unmapped read (*) - zero coverage") {
		auto stats = miint::ParseCigar("*");
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "aligned"), WithinRel(0.0, 0.001));
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "mapped"), WithinRel(0.0, 0.001));
	}

	SECTION("Deletions don't affect coverage calculation") {
		auto stats = miint::ParseCigar("10M5D10M");
		// 20 aligned / 20 total = 1.0 (deletions don't consume query)
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "aligned"), WithinRel(1.0, 0.001));
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "mapped"), WithinRel(1.0, 0.001));
	}

	SECTION("Reference skip (N) doesn't affect coverage") {
		auto stats = miint::ParseCigar("10M100N10M");
		// 20 aligned / 20 total = 1.0
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "aligned"), WithinRel(1.0, 0.001));
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "mapped"), WithinRel(1.0, 0.001));
	}

	SECTION("Invalid type throws exception") {
		auto stats = miint::ParseCigar("10M");
		REQUIRE_THROWS_AS(miint::ComputeQueryCoverage(stats, "invalid"), miint::InvalidInputException);
	}
}

TEST_CASE("ParseCigar - Integer overflow protection", "[alignment_functions]") {
	SECTION("Huge CIGAR operation length should throw") {
		// This would overflow int64_t without bounds checking
		REQUIRE_THROWS_AS(miint::ParseCigar("999999999999999999999M"), miint::InvalidInputException);
	}

	SECTION("Maximum safe CIGAR operation length") {
		// Just below overflow threshold should work
		auto stats = miint::ParseCigar("922337203685477580M");
		REQUIRE((stats.matches == 922337203685477580));
	}

	SECTION("CIGAR with multiple operations parsing") {
		// Multiple operations are allowed - overflow check is per-operation
		// Note: Accumulated stats could overflow, but that's a usage issue
		auto stats = miint::ParseCigar("100M50I");
		REQUIRE((stats.matches == 100));
		REQUIRE((stats.insertions == 50));
	}
}

TEST_CASE("ParseMd - Integer overflow protection", "[alignment_functions]") {
	SECTION("Huge MD match length should throw") {
		// This would overflow int64_t without bounds checking
		REQUIRE_THROWS_AS(miint::ParseMd("999999999999999999999"), miint::InvalidInputException);
	}

	SECTION("Maximum safe MD match length") {
		// Just below overflow threshold should work
		auto stats = miint::ParseMd("922337203685477580");
		REQUIRE((stats.matches == 922337203685477580));
	}

	SECTION("MD tag with multiple large numbers") {
		// Each individual number is safe but total might overflow in usage
		auto stats = miint::ParseMd("500000000000000000A500000000000000000");
		REQUIRE((stats.matches == 1000000000000000000));
		REQUIRE((stats.mismatches == 1));
	}
}

TEST_CASE("ComputeQueryCoverage - Aligned vs Mapped comparison", "[alignment_functions]") {
	using Catch::Matchers::WithinRel;

	SECTION("Insertions increase mapped coverage but not aligned") {
		auto stats = miint::ParseCigar("10H10M10I10H");
		// aligned: 10 / 40 = 0.25
		// mapped: 20 / 40 = 0.5
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "aligned"), WithinRel(0.25, 0.001));
		REQUIRE_THAT(miint::ComputeQueryCoverage(stats, "mapped"), WithinRel(0.5, 0.001));
	}

	SECTION("Without insertions, aligned == mapped") {
		auto stats = miint::ParseCigar("5H5S10M5S5H");
		double aligned = miint::ComputeQueryCoverage(stats, "aligned");
		double mapped = miint::ComputeQueryCoverage(stats, "mapped");
		REQUIRE_THAT(aligned, WithinRel(mapped, 0.001));
	}
}
