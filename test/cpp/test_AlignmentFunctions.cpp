#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <stdexcept>
#include <string>
// Include just the miint namespace parts we need for testing
// We need to define the exception type and the parsing functions

namespace miint {

// Copy the struct definitions and parsing functions from alignment_functions.cpp
struct CigarStats {
	int64_t matches = 0;
	int64_t match_ops = 0;
	int64_t mismatch_ops = 0;
	int64_t insertions = 0;
	int64_t deletions = 0;
	int64_t gap_opens = 0;
	int64_t alignment_columns = 0;
	int64_t soft_clips = 0;
	int64_t hard_clips = 0;
};

struct MdStats {
	int64_t matches = 0;
	int64_t mismatches = 0;
};

// Simple exception class for testing
class InvalidInputException : public std::runtime_error {
public:
	explicit InvalidInputException(const std::string &msg) : std::runtime_error(msg) {
	}
};

// Parse CIGAR string and extract statistics
static CigarStats ParseCigar(const std::string &cigar_str) {
	CigarStats stats;
	const char *cigar = cigar_str.data();
	size_t len = cigar_str.size();

	if (len == 0 || (len == 1 && cigar[0] == '*')) {
		return stats;
	}

	int64_t op_len = 0;
	char prev_op_type = '\0';

	for (size_t i = 0; i < len; i++) {
		char c = cigar[i];

		if (std::isdigit(c)) {
			op_len = op_len * 10 + (c - '0');
		} else {
			if (op_len == 0) {
				throw InvalidInputException("Invalid CIGAR string: operation without length");
			}

			switch (c) {
			case 'M':
				stats.matches += op_len;
				stats.alignment_columns += op_len;
				break;
			case '=':
				stats.matches += op_len;
				stats.match_ops += op_len;
				stats.alignment_columns += op_len;
				break;
			case 'X':
				stats.matches += op_len;
				stats.mismatch_ops += op_len;
				stats.alignment_columns += op_len;
				break;
			case 'I':
				stats.insertions += op_len;
				stats.alignment_columns += op_len;
				if (prev_op_type != 'I') {
					stats.gap_opens++;
				}
				break;
			case 'D':
				stats.deletions += op_len;
				stats.alignment_columns += op_len;
				if (prev_op_type != 'D') {
					stats.gap_opens++;
				}
				break;
			case 'N':
			case 'P':
				break;
			case 'S':
				stats.soft_clips += op_len;
				break;
			case 'H':
				stats.hard_clips += op_len;
				break;
			default:
				throw InvalidInputException("Invalid CIGAR operation: " + std::string(1, c));
			}

			prev_op_type = c;
			op_len = 0;
		}
	}

	return stats;
}

// Parse MD tag string and extract match/mismatch counts
static MdStats ParseMd(const std::string &md_str) {
	MdStats stats;
	const char *md = md_str.data();
	size_t len = md_str.size();

	if (len == 0) {
		return stats;
	}

	int64_t match_len = 0;

	for (size_t i = 0; i < len; i++) {
		char c = md[i];

		if (std::isdigit(c)) {
			match_len = match_len * 10 + (c - '0');
		} else if (c == '^') {
			if (match_len > 0) {
				stats.matches += match_len;
				match_len = 0;
			}
			i++;
			while (i < len && std::isalpha(md[i])) {
				i++;
			}
			i--;
		} else if (std::isalpha(c)) {
			if (match_len > 0) {
				stats.matches += match_len;
				match_len = 0;
			}
			stats.mismatches++;
		}
	}

	if (match_len > 0) {
		stats.matches += match_len;
	}

	return stats;
}

// Helper function to compute query length from CIGAR stats
// When include_hard_clips=false, result should match HTSlib's bam_cigar2qlen
static int64_t ComputeQueryLength(const CigarStats &stats, bool include_hard_clips) {
	int64_t length = stats.matches + stats.insertions + stats.soft_clips;
	if (include_hard_clips) {
		length += stats.hard_clips;
	}
	return length;
}

// Helper function to compute query coverage from CIGAR stats
static double ComputeQueryCoverage(const CigarStats &stats, const std::string &type) {
	// Total query length (always includes hard clips)
	int64_t query_length = stats.matches + stats.insertions + stats.soft_clips + stats.hard_clips;

	if (query_length == 0) {
		return 0.0; // Avoid division by zero, return 0% coverage
	}

	int64_t covered_bases = 0;
	if (type == "aligned") {
		// Only bases that align to reference (M, =, X)
		covered_bases = stats.matches;
	} else if (type == "mapped") {
		// Bases that are mapped (not clipped): M, I, =, X
		covered_bases = stats.matches + stats.insertions;
	} else {
		throw InvalidInputException("Invalid coverage type: " + type);
	}

	return static_cast<double>(covered_bases) / static_cast<double>(query_length);
}

} // namespace miint

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
