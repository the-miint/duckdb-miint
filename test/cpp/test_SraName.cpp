#include <catch2/catch_test_macros.hpp>

#include "sra_name.hpp"
#include "sra_source.hpp"
#include "sra_test_helpers.hpp"

#include <atomic>
#include <string>
#include <thread>
#include <vector>

using namespace miint::sra;
using miint_test::OpenSraFixture;
using miint_test::ThrowsTypeMentioning;

// fasterq-dump prints each spot's name ($sn), and the oracles hold it byte for byte. Names are
// rebuilt from a format and the spot's coordinates; the format comes from the run's skey index.
// Expected keys below were read from the fixtures by a separate implementation of the index format
// (a script, not committed), whose names for every spot of 13 fixtures agreed with the oracles.

namespace {
SkeyIndex IndexOf(const std::string &run) {
	auto ar = OpenSraFixture(run + ".sra");
	return SkeyIndex(ar->ReadFile("tbl/SEQUENCE/idx/skey"));
}
} // namespace

TEST_CASE("SkeyIndex projects a row to its tile's name format", "[sra]") {
	// One key per row (a contiguous projection).
	SkeyIndex dense = IndexOf("SRR34459993");
	CHECK(dense.Project(1) == "NB501337:903:HM25KBGXJ:1:11110:$X:$Y");
	CHECK(dense.Project(39) == "NB501337:903:HM25KBGXJ:4:23612:$X:$Y");
	// Rows sharing a key are found by binary search over bit-packed run starts (a sparse
	// projection): 617 rows, 150 keys.
	SkeyIndex sparse = IndexOf("ERR16653849");
	CHECK(sparse.Project(1) == "A00877:1574:H7LGJDRX5:1:2101:$X:$Y");
	CHECK(sparse.Project(300) == "A00877:1574:H7LGJDRX5:1:2173:$X:$Y");
	CHECK(sparse.Project(617) == "A00877:1574:H7LGJDRX5:1:2278:$X:$Y");
	CHECK(sparse.Project(300) == "A00877:1574:H7LGJDRX5:1:2173:$X:$Y"); // after another key
	// Outside the index's rows there is no key.
	CHECK(sparse.Project(618).empty());
	CHECK(sparse.Project(0).empty());
}

TEST_CASE("SkeyIndex returns whole names stored as keys", "[sra]") {
	// Loaders that cannot tokenize a name store it whole, one key per spot.
	SkeyIndex idx = IndexOf("SRR24770212");
	CHECK(idx.Project(1) == "V300101884L1C001R0051159705");
	CHECK(idx.Project(1618) == "V300101879L4C004R0141343127");
}

TEST_CASE("SkeyIndex has no key for rows whose name is stored elsewhere", "[sra]") {
	// SRR32151205 keeps some formats in its NAME_FMT column instead: the index maps those rows to
	// node 0.
	SkeyIndex idx = IndexOf("SRR32151205");
	CHECK(idx.Project(1) == "M01639:326:000000000-LFHYR:1:1102:$X:$Y");
	CHECK(idx.Project(30).empty());
	CHECK(idx.Project(1592) == "M02485:125:000000000-LG2B2:1:1106:$X:$Y");
	CHECK(idx.Project(1593).empty());
}

TEST_CASE("SkeyIndex can be read from several threads at once", "[sra]") {
	// A run's readers share one index (it can be hundreds of MB), so projecting must not write to
	// it.
	const SkeyIndex sparse = IndexOf("ERR16653849");
	std::vector<std::string> expected;
	for (int64_t row = 1; row <= 617; row++) {
		expected.push_back(sparse.Project(row));
	}
	std::atomic<int> mismatches {0};
	std::vector<std::thread> threads;
	for (int t = 0; t < 4; t++) {
		threads.emplace_back([&, t] {
			for (int pass = 0; pass < 20; pass++) {
				for (int64_t i = 0; i < 617; i++) {
					int64_t row = (i * (t + 1) * 7919) % 617 + 1; // a different order on each thread
					mismatches += sparse.Project(row) != expected[size_t(row - 1)];
				}
			}
		});
	}
	for (auto &th : threads) {
		th.join();
	}
	CHECK(mismatches == 0);
}

TEST_CASE("SkeyIndex refuses index versions absent from public runs", "[sra]") {
	auto bytes = OpenSraFixture("SRR34459993.sra")->ReadFile("tbl/SEQUENCE/idx/skey");
	bytes[4] = 3; // version 3: the original trie id coding
	CHECK(ThrowsTypeMentioning<SraUnsupportedError>([&] { SkeyIndex idx(bytes); }, "version 3"));
}

TEST_CASE("FormatSpotName substitutes coordinates", "[sra]") {
	CHECK(FormatSpotName("M00763:72:1:1101:$X:$Y", 9430, 1976, 7, "") == "M00763:72:1:1101:9430:1976");
	CHECK(FormatSpotName("t$x-$y", 1, 2, 7, "") == "t1-2");
	// Coordinates print signed.
	CHECK(FormatSpotName("$X", -5, 0, 7, "") == "-5");
	// $X%d adds 24576 * d: coordinates split across a camera's sections.
	CHECK(FormatSpotName("$X%2:$Y", 10, 3, 7, "") == "49162:3");
	CHECK(FormatSpotName("$R", std::nullopt, std::nullopt, 42, "") == "42");
	// 454 names: base 36 of (X << 12) | Y.
	CHECK(FormatSpotName("ABCDE$Q", 1, 2, 7, "") == "ABCDEAADF4");
	// Five base-36 digits hold (X << 12) | Y only below 36^5: a larger X has no $Q name.
	CHECK(FormatSpotName("$Q", 14761, 0xFFF, 7, "") == "999HT");
	CHECK(ThrowsTypeMentioning<SraFormatError>([] { FormatSpotName("$Q", 14763, 0, 7, ""); }, "X is too large"));
}

TEST_CASE("FormatSpotName copies what is not a substitution", "[sra]") {
	// A '$' before any other character is copied, and so is a final '$'.
	CHECK(FormatSpotName("a$b", 1, 2, 7, "") == "a$b");
	CHECK(FormatSpotName("ab$", 1, 2, 7, "") == "ab$");
	CHECK(FormatSpotName("$", 1, 2, 7, "") == "$");
	CHECK(FormatSpotName("plain", std::nullopt, std::nullopt, 7, "") == "plain");
}

TEST_CASE("FormatSpotName falls back to the stored name, or fails the spot", "[sra]") {
	CHECK(FormatSpotName("", std::nullopt, std::nullopt, 7, "10346.Webster.2.H4.181_0") == "10346.Webster.2.H4.181_0");
	CHECK(ThrowsTypeMentioning<SraFormatError>([] { FormatSpotName("", std::nullopt, std::nullopt, 7, ""); },
	                                           "spot 7 has no name"));
	CHECK(ThrowsTypeMentioning<SraFormatError>([] { FormatSpotName("t:$X", std::nullopt, 2, 7, ""); }, "X"));
	CHECK(ThrowsTypeMentioning<SraFormatError>(
	    [] { FormatSpotName(std::string(1100, 'a'), std::nullopt, std::nullopt, 7, ""); }, "1023"));
}

TEST_CASE("FormatSpotNameNoCoord substitutes only the row id", "[sra]") {
	CHECK(FormatSpotNameNoCoord("$R", 274142, "") == "274142");
	CHECK(FormatSpotNameNoCoord("HWI-M00962:109:1:1101:14974:2234", 1, "") == "HWI-M00962:109:1:1101:14974:2234");
	CHECK(FormatSpotNameNoCoord("a$X", 1, "") == "a$X");
	CHECK(FormatSpotNameNoCoord("", 1, "kept") == "kept");
}

TEST_CASE("Spot names must be valid UTF-8", "[sra]") {
	// Names become VARCHAR values, which DuckDB requires to be UTF-8. fasterq-dump copies whatever
	// bytes a run stores; a name that is not text is reported, naming the spot, rather than
	// passed on as an invalid string.
	CHECK_NOTHROW(CheckSpotName("M00763:72:1:1101:9430:1976", 7));
	CHECK_NOTHROW(CheckSpotName("caf\xc3\xa9", 7));      // two-byte sequence
	CHECK_NOTHROW(CheckSpotName("\xf0\x9f\xa7\xac", 7)); // four-byte sequence
	CHECK(ThrowsTypeMentioning<SraFormatError>([] { CheckSpotName("a\xffz", 7); }, "spot 7"));
	CHECK(ThrowsTypeMentioning<SraFormatError>([] { CheckSpotName("\xc3", 7); }, "UTF-8"));             // truncated
	CHECK(ThrowsTypeMentioning<SraFormatError>([] { CheckSpotName("\xc0\xaf", 7); }, "UTF-8"));         // overlong
	CHECK(ThrowsTypeMentioning<SraFormatError>([] { CheckSpotName("\xed\xa0\x80", 7); }, "UTF-8"));     // surrogate
	CHECK(ThrowsTypeMentioning<SraFormatError>([] { CheckSpotName("\xf4\x90\x80\x80", 7); }, "UTF-8")); // > U+10FFFF
}
