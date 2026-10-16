#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "gfa_parser.hpp"

#include <string>
#include <utility>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using miint::GfaParser;
using miint::GfaRecord;

namespace {

GfaRecord Parse(const std::string &line) {
	GfaParser parser;
	GfaRecord record;
	REQUIRE(parser.Parse(line, record));
	return record;
}

std::vector<std::pair<std::string, std::string>> Tags(const GfaRecord &record) {
	std::vector<std::pair<std::string, std::string>> out;
	for (const auto &tag : record.tags) {
		out.emplace_back(tag.key, tag.value);
	}
	return out;
}

using Fields = std::vector<std::string>;
using Pairs = std::vector<std::pair<std::string, std::string>>;

} // namespace

TEST_CASE("GfaParser splits each record type at its mandatory field count", "[gfa]") {
	// The specification fixes how many positional fields each record type has; everything
	// after them is tags. Getting a count wrong would either swallow a tag into the fields
	// or parse a positional field as a (malformed) tag.
	auto h = Parse("H\tVN:Z:1.1");
	CHECK(h.record_type == "H");
	CHECK(h.fields.empty());
	CHECK(Tags(h) == Pairs {{"VN", "1.1"}});

	auto s = Parse("S\ts1\tACGTACGT\tLN:i:8\tRC:i:120");
	CHECK(s.fields == Fields {"s1", "ACGTACGT"});
	CHECK(Tags(s) == Pairs {{"LN", "8"}, {"RC", "120"}});

	auto l = Parse("L\ts1\t+\ts2\t-\t4M\tID:Z:link1");
	CHECK(l.fields == Fields {"s1", "+", "s2", "-", "4M"});
	CHECK(Tags(l) == Pairs {{"ID", "link1"}});

	CHECK(Parse("C\ts1\t+\ts3\t-\t2\t6M").fields == Fields {"s1", "+", "s3", "-", "2", "6M"});
	CHECK(Parse("P\tp1\ts1+,s2-,s3+\t4M,5M").fields == Fields {"p1", "s1+,s2-,s3+", "4M,5M"});
	CHECK(Parse("W\tsample1\t1\tchr1\t0\t1016\t>s1<s2>s3").fields ==
	      Fields {"sample1", "1", "chr1", "0", "1016", ">s1<s2>s3"});

	// GFA 1.2's jump: like a link, with a distance (or *) in place of the overlap.
	auto j = Parse("J\ts1\t+\ts2\t-\t*\tSC:i:1");
	CHECK(j.fields == Fields {"s1", "+", "s2", "-", "*"});
	CHECK(Tags(j) == Pairs {{"SC", "1"}});
}

TEST_CASE("GfaParser reads myloasm's segment tags, whose names exceed two characters", "[gfa]") {
	auto s = Parse("S\tu491ctg\t*\tLN:i:196466\tDP:f:51.0\tDP2:f:30.0\tDP3:f:30.0\tMEDIAN_DP1:f:33.0");
	CHECK(s.fields == Fields {"u491ctg", "*"});
	CHECK(Tags(s) ==
	      Pairs {{"LN", "196466"}, {"DP", "51.0"}, {"DP2", "30.0"}, {"DP3", "30.0"}, {"MEDIAN_DP1", "33.0"}});
}

TEST_CASE("GfaParser gives a known record with no tags an empty tag list", "[gfa]") {
	auto l = Parse("L\ts2\t-\ts3\t+\t*");
	CHECK(l.has_tags);
	CHECK(l.tags.empty());
}

TEST_CASE("GfaParser keeps every field of an unknown record type positional", "[gfa]") {
	// A line in the format of myloasm's read placements: its last field looks tag-like but is
	// not a SAM-style tag. With no specification saying where positional fields end, none are
	// parsed.
	auto a = Parse("a\tu491ctg,0-2247\t517\t0\ts00517_linb_1009712_1023324_r+split0\t-\t2247\t"
	               "DP1:44,DP2:18,DP3:18,READ_LEN:13612");
	CHECK(a.record_type == "a");
	CHECK(a.fields == Fields {"u491ctg,0-2247", "517", "0", "s00517_linb_1009712_1023324_r+split0", "-", "2247",
	                          "DP1:44,DP2:18,DP3:18,READ_LEN:13612"});
	CHECK_FALSE(a.has_tags);
	CHECK(a.tags.empty());

	auto x = Parse("x\tcustom\trecord\twith:no:tags");
	CHECK(x.fields == Fields {"custom", "record", "with:no:tags"});
	CHECK_FALSE(x.has_tags);
}

TEST_CASE("GfaParser reports a comment line as carrying no record", "[gfa]") {
	GfaParser parser;
	GfaRecord record;
	CHECK_FALSE(parser.Parse("# Myloasm L tags: ol:i=overlap length (bp)", record));
}

TEST_CASE("GfaParser rejects malformed records with a message naming the defect", "[gfa]") {
	GfaParser parser;
	GfaRecord record;
	SECTION("an empty line") {
		CHECK_THROWS_WITH(parser.Parse("", record), ContainsSubstring("empty line"));
	}
	SECTION("an empty record type") {
		CHECK_THROWS_WITH(parser.Parse("\ts1\tACGT", record), ContainsSubstring("record type"));
	}
	SECTION("a known record type with too few fields") {
		CHECK_THROWS_WITH(parser.Parse("S\ts1", record), ContainsSubstring("S record needs 2"));
		CHECK_THROWS_WITH(parser.Parse("L\ts1\t+\ts2\t-", record), ContainsSubstring("L record needs 5"));
		CHECK_THROWS_WITH(parser.Parse("W\tsample1\t1\tchr1\t0\t1016", record), ContainsSubstring("W record needs 6"));
		CHECK_THROWS_WITH(parser.Parse("J\ts1\t+\ts2\t-", record), ContainsSubstring("a J record needs 5"));
	}
	SECTION("an orientation other than + or -") {
		// L, C and J join two oriented segments; any other orientation is not one GFA can
		// express. Each end is checked on its own.
		CHECK_THROWS_WITH(parser.Parse("L\ts1\t+\ts2\tx\t0M", record),
		                  ContainsSubstring("an L record's orientations must be + or -, found '+' and 'x'"));
		CHECK_THROWS_WITH(parser.Parse("L\ts1\t*\ts2\t-\t0M", record),
		                  ContainsSubstring("an L record's orientations must be + or -, found '*' and '-'"));
		CHECK_THROWS_WITH(parser.Parse("C\ts1\t+\ts3\t\t2\t6M", record),
		                  ContainsSubstring("a C record's orientations must be + or -, found '+' and ''"));
		CHECK_THROWS_WITH(parser.Parse("J\ts1\t++\ts2\t-\t*", record),
		                  ContainsSubstring("a J record's orientations must be + or -, found '++' and '-'"));
	}
	SECTION("a tag that is untyped, malformed, or repeated") {
		// GFA tags are always typed, so the lenient form PAF needs is not accepted here.
		CHECK_THROWS_WITH(parser.Parse("S\ts1\tACGT\tLN:4", record), ContainsSubstring("LN:4"));
		CHECK_THROWS_WITH(parser.Parse("S\ts1\tACGT\tLN:i:four", record), ContainsSubstring("LN:i:four"));
		CHECK_THROWS_WITH(parser.Parse("S\ts1\tACGT\tLN:i:4\tLN:i:4", record),
		                  ContainsSubstring("appears more than once"));
	}
	SECTION("a GFA 2 header") {
		// GFA 2 redefines S and adds record types; reading it as GFA 1 would mislabel fields.
		for (const char *version : {"2", "2.0", "2.1"}) {
			CAPTURE(version);
			CHECK_THROWS_WITH(parser.Parse(std::string("H\tVN:Z:") + version, record),
			                  ContainsSubstring("GFA 2 is not supported"));
		}
	}
	SECTION("a header naming any version other than GFA 1") {
		// Only GFA 1's layout is known; any other declared version is refused, not guessed at.
		// GFA 1 is "1" or "1.<minor>", so a version that merely starts with 1 is not it, and one
		// that merely starts with 2 is not GFA 2.
		for (const char *version : {"3.0", "abc", "10.0", "1abc", "", "20.0"}) {
			CAPTURE(version);
			CHECK_THROWS_WITH(parser.Parse(std::string("H\tVN:Z:") + version, record),
			                  ContainsSubstring("only GFA 1 is supported") &&
			                      ContainsSubstring(std::string("VN:Z:") + version + ")"));
		}
	}
	SECTION("a record type that is not text") {
		// Compressed bytes read as text, e.g. a gzipped graph without a .gz name. An unknown
		// record type is otherwise kept as written, so this would come back as rows.
		const char gzip_start[] = "\x1f\x8b\x08\x00\tx";
		CHECK_THROWS_WITH(parser.Parse(std::string(gzip_start, sizeof(gzip_start) - 1), record),
		                  ContainsSubstring("not text"));
		// Each end of printable ASCII on its own: a control byte, and DEL.
		CHECK_THROWS_WITH(parser.Parse("\x1f\tx", record), ContainsSubstring("not text"));
		CHECK_THROWS_WITH(parser.Parse("\x7f\tx", record), ContainsSubstring("not text"));
	}
}

TEST_CASE("GfaParser reads a GFA 1 header, or none", "[gfa]") {
	CHECK(Parse("H\tVN:Z:1.0").record_type == "H");
	CHECK(Parse("H\tVN:Z:1.2").record_type == "H");
	CHECK(Parse("H\tVN:Z:1").record_type == "H");
	CHECK(Parse("H\tTS:i:100").record_type == "H");
}

TEST_CASE("GfaParser keeps an unknown record type made of any printable text", "[gfa]") {
	// Space and '~' are the ends of printable ASCII; the not-text check must not refuse them.
	CHECK(Parse(" \tx").record_type == " ");
	CHECK(Parse("~\tx").record_type == "~");
}
