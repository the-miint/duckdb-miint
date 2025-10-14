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
			case 'S':
			case 'H':
			case 'P':
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
	SECTION("Soft clipping ignored") {
		auto stats = miint::ParseCigar("5S10M5S");
		REQUIRE((stats.matches == 10));
		REQUIRE((stats.alignment_columns == 10));
	}

	SECTION("Hard clipping ignored") {
		auto stats = miint::ParseCigar("5H10M5H");
		REQUIRE((stats.matches == 10));
		REQUIRE((stats.alignment_columns == 10));
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
