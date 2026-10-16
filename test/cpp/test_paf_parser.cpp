#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "paf_parser.hpp"

#include <string>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using miint::PafParser;
using miint::PafRecord;

namespace {

PafRecord Parse(const std::string &line) {
	PafParser parser;
	PafRecord record;
	parser.Parse(line, record);
	return record;
}

// kForward and kUnmapped are copied from data/paf/k12_minimap2.paf (minimap2 2.31), and
// kReverseOverhang is its ovhrev line with some tags dropped; data/paf/README.md says how the
// reads were built. kMyloasm is a made-up record in the format of myloasm 0.7.0's
// map_to_unitigs.paf.gz (untyped tags, a reverse-strand piece).
const std::string kForward = "fwd_1000_3000\t2000\t0\t2000\t+\tk12_100000_120000\t20000\t1000\t3000\t2000\t2000\t60\t"
                             "NM:i:0\tms:i:2000\tAS:i:2000\tnn:i:0\ttp:A:P\tcm:i:185\ts1:i:1973\ts2:i:0\tde:f:0\t"
                             "rl:i:0\tcg:Z:2000M";
const std::string kReverseOverhang = "ovhrev_17000_20000\t4500\t1500\t4500\t-\tk12_100000_120000\t20000\t17000\t"
                                     "20000\t3000\t3000\t60\tNM:i:0\ttp:A:P\tcg:Z:3000M";
const std::string kUnmapped = "unmapped\t3000\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\trl:i:0";
const std::string kMyloasm =
    "r953-s00933_ovh_1003000_1012000_4200000_r+split1\t8895\t6712\t8895\t-\tu517ctg\t172629\t0\t"
    "2183\t2159\t2159\t255\tshared_mini:186\tdiff_snp:6\tshared_snp:0";

// kForward with field `index` (0-based) replaced.
std::string WithField(size_t index, const std::string &value) {
	std::vector<std::string> fields;
	size_t start = 0;
	while (true) {
		size_t tab = kForward.find('\t', start);
		fields.push_back(kForward.substr(start, tab == std::string::npos ? std::string::npos : tab - start));
		if (tab == std::string::npos) {
			break;
		}
		start = tab + 1;
	}
	fields.at(index) = value;
	std::string out;
	for (size_t i = 0; i < fields.size(); i++) {
		out += (i ? "\t" : "") + fields[i];
	}
	return out;
}

} // namespace

TEST_CASE("PafParser converts the target interval to 1-based half-open", "[paf]") {
	// The read is reference bases [1000, 3000) in PAF's 0-based half-open terms. miint's
	// reference axis is 1-based half-open, so both ends move up by one and the length is
	// unchanged: [1001, 3001).
	auto r = Parse(kForward);
	CHECK(r.read_id == "fwd_1000_3000");
	CHECK(r.reference == "k12_100000_120000");
	CHECK(r.reference_length == 20000);
	CHECK(r.position == 1001);
	CHECK(r.stop_position == 3001);
	CHECK(r.stop_position - r.position == 2000);
	CHECK(r.strand == "+");
	CHECK(r.residue_matches == 2000);
	CHECK(r.block_length == 2000);
	CHECK(r.mapq == 60);
}

TEST_CASE("PafParser keeps the query interval 0-based on the read as sequenced", "[paf]") {
	// The reverse overhang read is foreign sequence [0, 1500) followed by the reverse
	// complement of the reference slice, so the aligned block is the read's last 3000
	// bases. PAF already reports forward-strand query coordinates; mirroring them would
	// give [0, 3000) and point at the foreign sequence.
	auto r = Parse(kReverseOverhang);
	CHECK(r.strand == "-");
	CHECK(r.query_length == 4500);
	CHECK(r.query_start == 1500);
	CHECK(r.query_stop == 4500);
	CHECK(r.position == 17001);
	CHECK(r.stop_position == 20001);
}

TEST_CASE("PafParser reports an unmapped record the way read_alignments does", "[paf]") {
	// minimap2 --paf-no-hit writes strand and target as '*' with zero coordinates.
	// read_alignments reports an unmapped SAM record as reference '*' at position 0, so a
	// query filtering unmapped rows works the same on both readers; position 1 would look
	// like an alignment at the reference's first base.
	auto r = Parse(kUnmapped);
	CHECK(r.read_id == "unmapped");
	CHECK(r.strand == "*");
	CHECK(r.reference == "*");
	CHECK(r.position == 0);
	CHECK(r.stop_position == 0);
	CHECK(r.query_length == 3000);
	CHECK(r.query_start == 0);
	CHECK(r.query_stop == 0);
}

TEST_CASE("PafParser keeps tags in order with their type letters dropped", "[paf]") {
	auto r = Parse(kForward);
	REQUIRE(r.tags.size() == 11);
	CHECK(r.tags.front().key == "NM");
	CHECK(r.tags.front().value == "0");
	CHECK(r.tags.back().key == "cg");
	CHECK(r.tags.back().value == "2000M");

	auto bare = Parse("bare\t2000\t0\t2000\t+\tk12\t20000\t1000\t3000\t2000\t2000\t60");
	CHECK(bare.tags.empty());
}

TEST_CASE("PafParser accepts myloasm's PAF, untyped tags and all, without interpreting it", "[paf]") {
	// myloasm writes untyped tags and puts '-' query coordinates on the reverse complement.
	// read_paf reports what the file says; fixing myloasm's coordinates is the myloasm
	// reader's job, so the values here are exactly as written (the target moved by one).
	auto r = Parse(kMyloasm);
	CHECK(r.read_id == "r953-s00933_ovh_1003000_1012000_4200000_r+split1");
	CHECK(r.query_start == 6712);
	CHECK(r.query_stop == 8895);
	CHECK(r.position == 1);
	CHECK(r.stop_position == 2184);
	CHECK(r.mapq == 255);
	REQUIRE(r.tags.size() == 3);
	CHECK(r.tags[0].key == "shared_mini");
	CHECK(r.tags[0].value == "186");
}

TEST_CASE("PafParser rejects malformed records with a message naming the defect", "[paf]") {
	SECTION("fewer than 12 fields") {
		CHECK_THROWS_WITH(Parse("a\t1\t0\t1\t+\tref\t10\t0\t1\t1\t1"), ContainsSubstring("12"));
	}
	SECTION("an empty line") {
		CHECK_THROWS_WITH(Parse(""), ContainsSubstring("empty line"));
	}
	SECTION("a number that is not a whole non-negative integer") {
		// A partial parse ("12x" as 12) or a clamped overflow would land on a plausible value.
		for (const char *bad : {"abc", "12x", "-1", "", "99999999999999999999", "+5"}) {
			INFO(bad);
			CHECK_THROWS_WITH(Parse(WithField(1, bad)), ContainsSubstring("query length"));
		}
		CHECK_THROWS_WITH(Parse(WithField(7, "x")), ContainsSubstring("target start"));
		CHECK_THROWS_WITH(Parse(WithField(10, "x")), ContainsSubstring("block length"));
	}
	SECTION("an interval outside its sequence or running backwards") {
		CHECK_THROWS_WITH(Parse(WithField(2, "2001")), ContainsSubstring("query"));   // start > stop
		CHECK_THROWS_WITH(Parse(WithField(3, "2001")), ContainsSubstring("query"));   // stop > length
		CHECK_THROWS_WITH(Parse(WithField(7, "3001")), ContainsSubstring("target"));  // start > stop
		CHECK_THROWS_WITH(Parse(WithField(8, "20001")), ContainsSubstring("target")); // stop > length
	}
	SECTION("a target end with no 1-based position") {
		// stop_position is target end + 1, which the largest integer does not have.
		CHECK_THROWS_WITH(Parse("r\t10\t0\t10\t+\tref\t9223372036854775807\t9223372036854775797\t"
		                        "9223372036854775807\t10\t10\t60"),
		                  ContainsSubstring("too large"));
	}
	SECTION("a strand that is not + or -") {
		CHECK_THROWS_WITH(Parse(WithField(4, "x")), ContainsSubstring("strand"));
		CHECK_THROWS_WITH(Parse(WithField(4, "")), ContainsSubstring("strand"));
	}
	SECTION("an unmapped record is '*' on both strand and target, or neither") {
		CHECK_THROWS_WITH(Parse(WithField(4, "*")), ContainsSubstring("unmapped"));
		CHECK_THROWS_WITH(Parse(WithField(5, "*")), ContainsSubstring("unmapped"));
		CHECK_THROWS_WITH(Parse("unmapped\t3000\t0\t0\t*\t*\t0\t5\t9\t0\t0\t0"), ContainsSubstring("unmapped"));
	}
	SECTION("a mapping quality above 255") {
		CHECK_THROWS_WITH(Parse(WithField(11, "256")), ContainsSubstring("mapping quality"));
	}
	SECTION("an empty read or target name") {
		CHECK_THROWS_WITH(Parse(WithField(0, "")), ContainsSubstring("query name"));
		CHECK_THROWS_WITH(Parse(WithField(5, "")), ContainsSubstring("target name"));
	}
	SECTION("a malformed tag") {
		CHECK_THROWS_WITH(Parse(kForward + "\tzz:i:abc"), ContainsSubstring("zz:i:abc"));
		CHECK_THROWS_WITH(Parse(kForward + "\tnotatag"), ContainsSubstring("notatag"));
	}
}
