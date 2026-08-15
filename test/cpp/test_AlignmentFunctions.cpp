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

TEST_CASE("ComputeCigarIdentity - Extended CIGAR ops", "[alignment_functions]") {
	using Catch::Matchers::WithinRel;

	SECTION("Perfect match (all =)") {
		auto stats = miint::ParseCigar("10=");
		auto result = miint::ComputeCigarIdentity(stats);
		REQUIRE(result.has_value());
		REQUIRE_THAT(*result, WithinRel(1.0, 0.001));
	}

	SECTION("Mixed =/X") {
		auto stats = miint::ParseCigar("5=5X");
		auto result = miint::ComputeCigarIdentity(stats);
		REQUIRE(result.has_value());
		REQUIRE_THAT(*result, WithinRel(0.5, 0.001));
	}

	SECTION("With insertions: gaps lower identity") {
		auto stats = miint::ParseCigar("5=2I3=");
		// match_ops=8, alignment_columns = 8 + 2 = 10 → 0.8
		auto result = miint::ComputeCigarIdentity(stats);
		REQUIRE(result.has_value());
		REQUIRE_THAT(*result, WithinRel(0.8, 0.001));
	}

	SECTION("With deletions: gaps lower identity") {
		auto stats = miint::ParseCigar("5=3D5=");
		// match_ops=10, alignment_columns = 10 + 3 = 13 → 10/13
		auto result = miint::ComputeCigarIdentity(stats);
		REQUIRE(result.has_value());
		REQUIRE_THAT(*result, WithinRel(10.0 / 13.0, 0.001));
	}

	SECTION("Soft clips excluded from identity") {
		auto stats = miint::ParseCigar("3S5=2X3=3S");
		// match_ops=8, mismatch_ops=2, alignment_columns=10 → 0.8
		auto result = miint::ComputeCigarIdentity(stats);
		REQUIRE(result.has_value());
		REQUIRE_THAT(*result, WithinRel(0.8, 0.001));
	}

	SECTION("Hard clips excluded from identity") {
		auto stats = miint::ParseCigar("5H5=5X5H");
		auto result = miint::ComputeCigarIdentity(stats);
		REQUIRE(result.has_value());
		REQUIRE_THAT(*result, WithinRel(0.5, 0.001));
	}

	SECTION("N op (RNA skip) ignored") {
		auto stats = miint::ParseCigar("5=100N5=");
		// match_ops=10, alignment_columns=10 (N not counted)
		auto result = miint::ComputeCigarIdentity(stats);
		REQUIRE(result.has_value());
		REQUIRE_THAT(*result, WithinRel(1.0, 0.001));
	}
}

TEST_CASE("ComputeCigarIdentity - Degenerate and ambiguous cases", "[alignment_functions]") {
	SECTION("Legacy M-only: no =/X, return no value") {
		auto stats = miint::ParseCigar("10M");
		REQUIRE(!miint::ComputeCigarIdentity(stats).has_value());
	}

	SECTION("Mixed M with =: inconsistent, return no value") {
		auto stats = miint::ParseCigar("5M5=");
		REQUIRE(!miint::ComputeCigarIdentity(stats).has_value());
	}

	SECTION("Mixed M with X: inconsistent, return no value") {
		auto stats = miint::ParseCigar("5M5X");
		REQUIRE(!miint::ComputeCigarIdentity(stats).has_value());
	}

	SECTION("Pure insertion: no =/X, return no value") {
		auto stats = miint::ParseCigar("10I");
		REQUIRE(!miint::ComputeCigarIdentity(stats).has_value());
	}

	SECTION("Pure deletion: no =/X, return no value") {
		auto stats = miint::ParseCigar("10D");
		REQUIRE(!miint::ComputeCigarIdentity(stats).has_value());
	}

	SECTION("Soft-clip only: no =/X, return no value") {
		auto stats = miint::ParseCigar("100S");
		REQUIRE(!miint::ComputeCigarIdentity(stats).has_value());
	}

	SECTION("Hard-clip only: no =/X, return no value") {
		auto stats = miint::ParseCigar("100H");
		REQUIRE(!miint::ComputeCigarIdentity(stats).has_value());
	}

	SECTION("Empty CIGAR: no =/X, return no value") {
		auto stats = miint::ParseCigar("");
		REQUIRE(!miint::ComputeCigarIdentity(stats).has_value());
	}

	SECTION("Unmapped (*): no =/X, return no value") {
		auto stats = miint::ParseCigar("*");
		REQUIRE(!miint::ComputeCigarIdentity(stats).has_value());
	}
}

// Helper: parse then compute, so the test bodies read as CIGAR -> intervals.
static std::vector<miint::QueryInterval> Intervals(const std::string &cigar, bool reverse, const std::string &type) {
	return miint::ComputeQueryIntervals(miint::ParseCigarOperations(cigar), reverse, type);
}

static bool Equals(const std::vector<miint::QueryInterval> &got, const std::vector<std::pair<int64_t, int64_t>> &want) {
	if (got.size() != want.size()) {
		return false;
	}
	for (size_t i = 0; i < got.size(); i++) {
		if (got[i].start != want[i].first || got[i].stop != want[i].second) {
			return false;
		}
	}
	return true;
}

TEST_CASE("ComputeQueryIntervals - Forward strand", "[alignment_functions]") {
	SECTION("Unclipped alignment covers the whole read") {
		REQUIRE(Equals(Intervals("6000=", false, "aligned"), {{0, 6000}}));
	}

	SECTION("Trailing soft clip is not covered") {
		REQUIRE(Equals(Intervals("3000=3000S", false, "aligned"), {{0, 3000}}));
	}

	SECTION("Leading hard clip offsets the covered block") {
		REQUIRE(Equals(Intervals("3000H3000=", false, "aligned"), {{3000, 6000}}));
	}

	SECTION("Trailing hard clip is not covered") {
		REQUIRE(Equals(Intervals("3000=3000H", false, "aligned"), {{0, 3000}}));
	}

	SECTION("Leading soft and trailing hard clip both offset and truncate") {
		REQUIRE(Equals(Intervals("10S30=5H", false, "aligned"), {{10, 40}}));
	}

	SECTION("Legacy M counts as aligned") {
		REQUIRE(Equals(Intervals("50M", false, "aligned"), {{0, 50}}));
	}

	SECTION("Mismatches are covered") {
		REQUIRE(Equals(Intervals("10=2X10=", false, "aligned"), {{0, 22}}));
	}
}

TEST_CASE("ComputeQueryIntervals - Operations that consume no query merge runs", "[alignment_functions]") {
	SECTION("Deletion leaves the query contiguous") {
		REQUIRE(Equals(Intervals("100=10D100=", false, "aligned"), {{0, 200}}));
	}

	SECTION("Reference skip leaves the query contiguous") {
		REQUIRE(Equals(Intervals("100=10N100=", false, "aligned"), {{0, 200}}));
	}

	SECTION("Padding leaves the query contiguous") {
		REQUIRE(Equals(Intervals("100=10P100=", false, "aligned"), {{0, 200}}));
	}
}

TEST_CASE("ComputeQueryIntervals - Insertions distinguish aligned from mapped", "[alignment_functions]") {
	SECTION("aligned: inserted bases are a gap in coverage") {
		REQUIRE(Equals(Intervals("100=10I100=", false, "aligned"), {{0, 100}, {110, 210}}));
	}

	SECTION("mapped: inserted bases are covered, so the block is contiguous") {
		REQUIRE(Equals(Intervals("100=10I100=", false, "mapped"), {{0, 210}}));
	}
}

TEST_CASE("ComputeQueryIntervals - Reverse strand mirrors onto the read axis", "[alignment_functions]") {
	// The CIGAR is written in reference orientation, so the leading clip is at the
	// read's 3' end. Without the mirror, intervals from a reverse fragment would be
	// placed on the wrong half of the read and could not be pooled with its mates.
	SECTION("Leading soft clip in reference orientation is trailing on the read") {
		REQUIRE(Equals(Intervals("3000S3000=", true, "aligned"), {{0, 3000}}));
	}

	SECTION("Trailing hard clip in reference orientation is leading on the read") {
		REQUIRE(Equals(Intervals("3000=3000H", true, "aligned"), {{3000, 6000}}));
	}

	SECTION("Multiple intervals stay ascending and non-overlapping after mirroring") {
		// Reference orientation: [0,100) and [110,160) of a 160 bp read.
		// Mirrored: [60,160) and [0,50).
		REQUIRE(Equals(Intervals("100=10I50=", true, "aligned"), {{0, 50}, {60, 160}}));
	}

	SECTION("An unclipped alignment is unchanged by mirroring") {
		REQUIRE(Equals(Intervals("6000=", true, "aligned"), {{0, 6000}}));
	}
}

TEST_CASE("ComputeQueryIntervals - Degenerate input", "[alignment_functions]") {
	SECTION("Unmapped CIGAR covers nothing") {
		REQUIRE(Intervals("*", false, "aligned").empty());
	}

	SECTION("Empty CIGAR covers nothing") {
		REQUIRE(Intervals("", false, "aligned").empty());
	}

	SECTION("Clip-only CIGAR covers nothing") {
		REQUIRE(Intervals("100S", false, "aligned").empty());
		REQUIRE(Intervals("100H", false, "aligned").empty());
	}

	SECTION("Pure insertion covers nothing under aligned but everything under mapped") {
		REQUIRE(Intervals("10I", false, "aligned").empty());
		REQUIRE(Equals(Intervals("10I", false, "mapped"), {{0, 10}}));
	}

	SECTION("Unrecognised type is rejected") {
		REQUIRE_THROWS_AS(Intervals("10=", false, "covered"), miint::InvalidInputException);
	}

	SECTION("Invalid CIGAR is rejected") {
		REQUIRE_THROWS_AS(Intervals("10Z", false, "aligned"), miint::InvalidInputException);
	}
}

TEST_CASE("ComputeQueryIntervals - Agrees with ComputeQueryCoverage on one fragment", "[alignment_functions]") {
	// The reason a caller can adopt the pooled coverage unconditionally: on a group of
	// one fragment the union of these intervals, over the full query length, is exactly
	// what cigar_query_coverage already reports. If this drifts, the two functions
	// disagree on non-wrapping reads and the whole design premise is broken.
	const std::vector<std::string> cigars = {"6000=",       "3000=3000S",  "3000H3000=", "10S30=5H",
	                                         "100=10I100=", "100=10D100=", "50M",        "10=2X10="};
	for (const auto &type : {std::string("aligned"), std::string("mapped")}) {
		for (const auto &cigar : cigars) {
			auto stats = miint::ParseCigar(cigar);
			auto query_length = miint::ComputeQueryLength(stats, /*include_hard_clips=*/true);
			int64_t covered = 0;
			for (const auto &iv : Intervals(cigar, false, type)) {
				covered += iv.stop - iv.start;
			}
			auto from_intervals = static_cast<double>(covered) / static_cast<double>(query_length);
			REQUIRE_THAT(from_intervals, Catch::Matchers::WithinAbs(miint::ComputeQueryCoverage(stats, type), 1e-12));
		}
	}
}

TEST_CASE("ComputeQueryIntervals - Origin-spanning fragments tile the read exactly once", "[alignment_functions]") {
	// The case from the issue: a 6 kb read across the origin of a 30 kb circular contig.
	// Neither fragment covers more than half, but together they must tile the read with
	// no overlap -- which is what makes the union, rather than the sum, the right
	// operator.
	auto primary = Intervals("3000=3000S", false, "aligned");
	auto supplementary = Intervals("3000H3000=", false, "aligned");
	REQUIRE(Equals(primary, {{0, 3000}}));
	REQUIRE(Equals(supplementary, {{3000, 6000}}));
	REQUIRE((primary[0].stop == supplementary[0].start));

	// Same read sequenced in the other orientation: both fragments carry FLAG 0x10 and
	// must still tile the read.
	auto rev_primary = Intervals("3000S3000=", true, "aligned");
	auto rev_supplementary = Intervals("3000=3000H", true, "aligned");
	REQUIRE(Equals(rev_primary, {{0, 3000}}));
	REQUIRE(Equals(rev_supplementary, {{3000, 6000}}));
}
