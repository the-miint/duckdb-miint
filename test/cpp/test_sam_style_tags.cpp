#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "sam_style_tags.hpp"

#include <string>
#include <utility>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using miint::ParseSamStyleTags;
using miint::SamStyleTag;
using miint::TagSyntax;

namespace {

std::vector<std::pair<std::string, std::string>> Parse(const std::vector<std::string_view> &fields,
                                                       TagSyntax syntax = TagSyntax::Typed, size_t first = 0) {
	std::vector<SamStyleTag> tags;
	ParseSamStyleTags(fields, first, syntax, tags);
	std::vector<std::pair<std::string, std::string>> out;
	for (const auto &tag : tags) {
		out.emplace_back(tag.key, tag.value);
	}
	return out;
}

using Pairs = std::vector<std::pair<std::string, std::string>>;

} // namespace

TEST_CASE("Typed tags of every SAM/GFA type parse to key and value", "[sam_style_tags]") {
	// Values are the spec's own shapes; the type letter is validated and dropped.
	CHECK(Parse({"tp:A:P", "ST:A:+"}) == Pairs {{"tp", "P"}, {"ST", "+"}});
	CHECK(Parse({"NM:i:0", "KC:i:-42", "XX:i:+7"}) == Pairs {{"NM", "0"}, {"KC", "-42"}, {"XX", "+7"}});
	CHECK(Parse({"de:f:0", "fl:f:1.5e-3", "f2:f:.5", "f3:f:-2E+10"}) ==
	      Pairs {{"de", "0"}, {"fl", "1.5e-3"}, {"f2", ".5"}, {"f3", "-2E+10"}});
	CHECK(Parse({"cg:Z:2000M", "nm:Z:segment three", "em:Z:"}) ==
	      Pairs {{"cg", "2000M"}, {"nm", "segment three"}, {"em", ""}});
	CHECK(Parse({"js:J:{\"k\":[1,2]}"}) == Pairs {{"js", "{\"k\":[1,2]}"}});
	CHECK(Parse({"hx:H:1AE301", "h0:H:"}) == Pairs {{"hx", "1AE301"}, {"h0", ""}});
	CHECK(Parse({"ba:B:i,1,-2,3", "bf:B:f,0.5,1.25", "bc:B:c"}) ==
	      Pairs {{"ba", "i,1,-2,3"}, {"bf", "f,0.5,1.25"}, {"bc", "c"}});
}

TEST_CASE("Typed tags may have names longer than two characters", "[sam_style_tags]") {
	// The SAM and GFA specifications fix tag names at two characters, but myloasm's GFA
	// writes DP2:f:30.0 and MEDIAN_DP1:f:33.0. The type letter is still required and the
	// value still validated; reading DP2:f:30.0 as an untyped key with value "f:30.0"
	// would be wrong in either mode.
	const std::vector<std::string_view> myloasm = {"LN:i:196466", "DP2:f:30.0", "MEDIAN_DP1:f:33.0"};
	const Pairs expected = {{"LN", "196466"}, {"DP2", "30.0"}, {"MEDIAN_DP1", "33.0"}};
	CHECK(Parse(myloasm) == expected);
	CHECK(Parse(myloasm, TagSyntax::TypedOrUntyped) == expected);
	CHECK_THROWS_WITH(Parse({"MEDIAN_DP1:f:abc"}), ContainsSubstring("MEDIAN_DP1:f:abc"));
}

TEST_CASE("Tags keep field order and start at the requested field", "[sam_style_tags]") {
	CHECK(Parse({"mandatory", "also mandatory", "b2:i:2", "a1:i:1"}, TagSyntax::Typed, 2) ==
	      Pairs {{"b2", "2"}, {"a1", "1"}});
	CHECK(Parse({"only mandatory"}, TagSyntax::Typed, 1).empty());

	std::vector<SamStyleTag> tags {{"stale", "value"}};
	ParseSamStyleTags({"NM:i:1"}, 0, TagSyntax::Typed, tags);
	REQUIRE(tags.size() == 1);
	CHECK(tags[0].key == "NM");
}

TEST_CASE("A value that does not match its declared type is an error", "[sam_style_tags]") {
	// Accepting these would hand downstream casts a value the type letter promised was
	// valid; failing here names the field instead of a later, anonymous cast error.
	// Of strtod's non-numeric forms only nan and inf (below) are floats here: "nanx", a
	// truncated "infinit", strtod's nan(chars) and hex floats are not.
	for (const char *field :
	     {"NM:i:abc",   "NM:i:",      "NM:i:1.5", "de:f:nanx", "de:f:in",       "de:f:infinit", "de:f:nan(1)",
	      "de:f:0x1p3", "de:f:2.",    "de:f:",    "tp:A:PP",   "tp:A:",         "hx:H:1AE",     "hx:H:1ae3",
	      "ba:B:x,1",   "ba:B:i,1,a", "ba:B:i,",  "ba:B:",     "ba:B:f,1,nanx", "nm:Z:bad\x01"}) {
		INFO(field);
		CHECK_THROWS_WITH(Parse({field}), ContainsSubstring(field));
		// A typed-shaped field stays typed when untyped fields are also allowed.
		CHECK_THROWS_WITH(Parse({field}, TagSyntax::TypedOrUntyped), ContainsSubstring(field));
	}
}

TEST_CASE("A float may be nan or an infinity, as htslib reads one", "[sam_style_tags]") {
	// Each form htslib's strtod reads as nan or an infinity, in f values and B:f elements (why:
	// IsFloat). The value is kept as written.
	for (const char *value : {"nan", "NaN", "-nan", "inf", "-inf", "+Inf", "INF", "infinity", "-Infinity"}) {
		INFO(value);
		CHECK(Parse({std::string("de:f:") + value}) == Pairs {{"de", value}});
	}
	CHECK(Parse({"ba:B:f,1.5,nan,-inf"}) == Pairs {{"ba", "f,1.5,nan,-inf"}});
}

TEST_CASE("Untyped key:value fields are accepted only when asked for", "[sam_style_tags]") {
	// myloasm's PAF writes untyped tags; strict callers (GFA) must still reject them.
	const std::vector<std::string_view> myloasm = {"shared_mini:284", "diff_snp:0", "shared_snp:9"};
	CHECK_THROWS_WITH(Parse(myloasm), ContainsSubstring("shared_mini:284"));
	CHECK(Parse(myloasm, TagSyntax::TypedOrUntyped) ==
	      Pairs {{"shared_mini", "284"}, {"diff_snp", "0"}, {"shared_snp", "9"}});
	// Typed and untyped fields can sit on the same line; the value is everything after the
	// first ':' for an untyped field.
	CHECK(Parse({"NM:i:0", "note:a:b"}, TagSyntax::TypedOrUntyped) == Pairs {{"NM", "0"}, {"note", "a:b"}});
}

TEST_CASE("A field that is not a tag at all is an error in either mode", "[sam_style_tags]") {
	for (const char *field : {"abc", ":value", "1Z:Z:x", "a b:1", ""}) {
		INFO(field);
		CHECK_THROWS_WITH(Parse({field}), ContainsSubstring("is not a tag"));
		CHECK_THROWS_WITH(Parse({field}, TagSyntax::TypedOrUntyped), ContainsSubstring("is not a tag"));
	}
}

TEST_CASE("A key that appears twice is an error", "[sam_style_tags]") {
	// The tags become a MAP column, whose keys must be unique; keeping either value would
	// silently discard the other.
	CHECK_THROWS_WITH(Parse({"NM:i:0", "NM:i:1"}), ContainsSubstring("NM"));
	CHECK_THROWS_WITH(Parse({"xx:i:1", "xx:2"}, TagSyntax::TypedOrUntyped), ContainsSubstring("xx"));
}
