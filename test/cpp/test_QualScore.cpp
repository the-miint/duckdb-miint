#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <vector>
#include <cstdint>
#include <QualScore.hpp>

TEST_CASE("miint::find_low_quality_window", "[filter]") {
	SECTION("no low quality window returns none") {
		// 'M' => ASCII 77 → 77-33 = 44
		std::vector<uint8_t> qs(10, 44);
		auto obs = miint::find_low_quality_window(qs,
		                                          /*min_quality=*/20,
		                                          /*window_length=*/2);
		CHECK((obs == qs.size() + 1));
	}

	SECTION("bases exactly at min_quality are not low") {
		// all 44 ('M')
		std::vector<uint8_t> qs(10, 44);
		auto obs = miint::find_low_quality_window(qs,
		                                          /*min_quality=*/44,
		                                          /*window_length=*/2);
		CHECK((obs == qs.size() + 1));
	}

	SECTION("simple window of length 2 in MMM++MM is at index 3") {
		// 'M'→44, '+'→43-33=10
		std::vector<uint8_t> qs = {
		    44, 44, 44, // M, M, M
		    10, 10,     // +, +
		    44, 44      // M, M
		};
		auto obs = miint::find_low_quality_window(qs,
		                                          /*min_quality=*/15,
		                                          /*window_length=*/2);
		CHECK((obs == 3));
	}

	SECTION("window of length 3 in M++MM+++MM is at index 5") {
		// M, +, +, M, M, +, +, +, M, M
		// 44,10,10,44,44,10,10,10,44,44
		std::vector<uint8_t> qs = {44, 10, 10, 44, 44, 10, 10, 10, 44, 44};
		auto obs = miint::find_low_quality_window(qs,
		                                          /*min_quality=*/15,
		                                          /*window_length=*/3);
		CHECK((obs == 5));
	}

	SECTION("first two +'s in ++MMMM gives window at 0") {
		// +, +, M, M, M, M → 10,10,44,44,44,44
		std::vector<uint8_t> qs = {10, 10, 44, 44, 44, 44};
		auto obs = miint::find_low_quality_window(qs,
		                                          /*min_quality=*/11,
		                                          /*window_length=*/2);
		CHECK((obs == 0));
	}

	SECTION("first length-3 window in M++MMMM+++ is at index 7") {
		// M,+,+,M,M,M,M,+,+,+ → 44,10,10,44,44,44,44,10,10,10
		std::vector<uint8_t> qs = {44, 10, 10, 44, 44, 44, 44, 10, 10, 10};
		auto obs = miint::find_low_quality_window(qs,
		                                          /*min_quality=*/11,
		                                          /*window_length=*/3);
		CHECK((obs == 7));
	}

	SECTION("longer run still flags at its start in MMMMMMM+++") {
		// M×7 then +×3: 44×7,10×3
		std::vector<uint8_t> qs = {44, 44, 44, 44, 44, 44, 44, 10, 10, 10};
		auto obs = miint::find_low_quality_window(qs,
		                                          /*min_quality=*/11,
		                                          /*window_length=*/2);
		CHECK((obs == 7));
	}

	SECTION("when multiple windows exist, the first is returned") {
		// M,L,+,+,M,M,M,+,+,+
		// 'L'→76-33=43
		std::vector<uint8_t> qs1 = {44, 43, 10, 10, 44, 44, 44, 10, 10, 10};
		auto o1 = miint::find_low_quality_window(qs1,
		                                         /*min_quality=*/20,
		                                         /*window_length=*/2);
		CHECK((o1 == 2));

		// +,+,M,L,+,+,+,M,+,+,+,M,M,+,+
		// 10,10,44,43,10,10,10,44,10,10,10,44,44,10,10
		std::vector<uint8_t> qs2 = {10, 10, 44, 43, 10, 10, 10, 44, 10, 10, 10, 44, 44, 10, 10};
		auto o2 = miint::find_low_quality_window(qs2,
		                                         /*min_quality=*/20,
		                                         /*window_length=*/3);
		CHECK((o2 == 4));
	}

	SECTION("when all runs are too short, returns none") {
		// same as qs2 above but window_length=4
		std::vector<uint8_t> qs = {10, 10, 44, 43, 10, 10, 10, 44, 10, 10, 10, 44, 44, 10, 10};
		auto obs = miint::find_low_quality_window(qs,
		                                          /*min_quality=*/20,
		                                          /*window_length=*/4);
		CHECK((obs == qs.size() + 1));
	}
}
