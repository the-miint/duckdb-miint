// Unit tests for the pure codec used by id-column ingress/egress.
//
// These cover the BIGINT-side semantics that drive output emission when the
// caller's id column is BIGINT-typed. The DuckDB-aware wrappers
// (ExtractIdColumnAsStrings, EmitIdColumnFromStrings) live in
// id_column_utils.{hpp,cpp} and are covered by SQL integration tests because
// the unit-test binary does not link against the duckdb library.

#include "id_column_codec.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

using miint::FormatIdFromInt64;
using miint::ParseIdAsInt64;

TEST_CASE("ParseIdAsInt64 happy path", "[id_codec]") {
	SECTION("plain positive integer") {
		auto v = ParseIdAsInt64("123");
		REQUIRE(v.has_value());
		CHECK(*v == 123);
	}

	SECTION("zero") {
		auto v = ParseIdAsInt64("0");
		REQUIRE(v.has_value());
		CHECK(*v == 0);
	}

	SECTION("negative integer") {
		auto v = ParseIdAsInt64("-456");
		REQUIRE(v.has_value());
		CHECK(*v == -456);
	}

	SECTION("INT64_MAX boundary") {
		auto v = ParseIdAsInt64("9223372036854775807");
		REQUIRE(v.has_value());
		CHECK(*v == std::numeric_limits<int64_t>::max());
	}

	SECTION("INT64_MIN boundary") {
		auto v = ParseIdAsInt64("-9223372036854775808");
		REQUIRE(v.has_value());
		CHECK(*v == std::numeric_limits<int64_t>::min());
	}
}

TEST_CASE("ParseIdAsInt64 NULL sentinels", "[id_codec]") {
	SECTION("empty string returns nullopt (NULL)") {
		auto v = ParseIdAsInt64("");
		CHECK_FALSE(v.has_value());
	}

	SECTION("SAM unmapped sentinel '*' returns nullopt (NULL)") {
		auto v = ParseIdAsInt64("*");
		CHECK_FALSE(v.has_value());
	}
}

TEST_CASE("ParseIdAsInt64 invalid input throws", "[id_codec]") {
	SECTION("non-numeric string throws") {
		CHECK_THROWS_AS(ParseIdAsInt64("abc"), std::invalid_argument);
	}

	SECTION("float-like string throws (no implicit truncation)") {
		CHECK_THROWS_AS(ParseIdAsInt64("1.5"), std::invalid_argument);
	}

	SECTION("trailing garbage throws (no partial parse)") {
		CHECK_THROWS_AS(ParseIdAsInt64("123abc"), std::invalid_argument);
	}

	SECTION("leading whitespace throws (caller must clean input)") {
		CHECK_THROWS_AS(ParseIdAsInt64(" 123"), std::invalid_argument);
	}

	SECTION("overflow throws") {
		CHECK_THROWS_AS(ParseIdAsInt64("99999999999999999999"), std::out_of_range);
	}

	SECTION("SAM same-as-primary sentinel '=' throws (caller responsibility)") {
		// '=' must be resolved to the row's reference value by the caller
		// before invoking the codec. Failing loud here surfaces that contract.
		CHECK_THROWS_AS(ParseIdAsInt64("="), std::invalid_argument);
	}
}

TEST_CASE("FormatIdFromInt64", "[id_codec]") {
	SECTION("zero") {
		CHECK(FormatIdFromInt64(0) == "0");
	}

	SECTION("positive integer") {
		CHECK(FormatIdFromInt64(42) == "42");
	}

	SECTION("negative integer") {
		CHECK(FormatIdFromInt64(-1) == "-1");
	}

	SECTION("INT64_MAX") {
		CHECK(FormatIdFromInt64(std::numeric_limits<int64_t>::max()) == "9223372036854775807");
	}

	SECTION("INT64_MIN") {
		CHECK(FormatIdFromInt64(std::numeric_limits<int64_t>::min()) == "-9223372036854775808");
	}
}

TEST_CASE("ParseIdAsInt64 documented non-canonical inputs", "[id_codec]") {
	// These are accepted by from_chars but lose information through the
	// FormatIdFromInt64 round-trip. Asserting the behaviour explicitly so the
	// API contract is documented in code, not just folklore.
	SECTION("leading zeros: '01' parses to 1; round-trip canonicalises to '1'") {
		auto v = ParseIdAsInt64("01");
		REQUIRE(v.has_value());
		CHECK(*v == 1);
		CHECK(FormatIdFromInt64(*v) == "1");
	}

	SECTION("'-0' parses to 0; round-trip canonicalises to '0'") {
		auto v = ParseIdAsInt64("-0");
		REQUIRE(v.has_value());
		CHECK(*v == 0);
		CHECK(FormatIdFromInt64(*v) == "0");
	}

	SECTION("embedded NUL throws (no partial parse)") {
		// std::string with explicit length to preserve the embedded NUL.
		std::string s("12\x00"
		              "3",
		              4);
		REQUIRE(s.size() == 4);
		CHECK_THROWS_AS(ParseIdAsInt64(s), std::invalid_argument);
	}
}

TEST_CASE("ParseIdAsInt64 / FormatIdFromInt64 round-trip", "[id_codec]") {
	const int64_t values[] = {0,
	                          1,
	                          -1,
	                          42,
	                          -42,
	                          1'000'000'000'000LL,
	                          std::numeric_limits<int64_t>::max(),
	                          std::numeric_limits<int64_t>::min()};
	for (auto v : values) {
		auto round = ParseIdAsInt64(FormatIdFromInt64(v));
		REQUIRE(round.has_value());
		CHECK(*round == v);
	}
}
