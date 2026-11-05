#include <catch2/catch_test_macros.hpp>

#include "sra_run.hpp"
#include "sra_test_helpers.hpp"

#include <zlib.h>

#include <fstream>
#include <map>
#include <string>
#include <vector>

using namespace miint::sra;
using miint_test::OpenSraFixture;
using miint_test::ThrowsTypeMentioning;

// The reader must write exactly the reads fasterq-dump writes, in its order, to the same files.
// Routing rules follow tbl_join.c (sra-tools 3.4.1); the fixtures and oracles are described in
// data/sra/README.md. Spot names ($sn) are not compared here: they come in phase 6.

namespace {
Spot MakeSpot(std::string bases, std::vector<int32_t> start, std::vector<uint32_t> len, std::vector<uint8_t> type) {
	Spot s;
	s.row_id = 7;
	s.quality.assign(bases.size(), 30);
	s.bases = std::move(bases);
	s.read_start = std::move(start);
	s.read_len = std::move(len);
	s.read_type = std::move(type);
	s.has_read_type = !s.read_type.empty();
	return s;
}

constexpr uint8_t TECH = 0;
constexpr uint8_t BIO = 1;

std::vector<SraRead> Route(const Spot &s, SplitMode mode, bool tech = false) {
	std::vector<SraRead> out;
	RouteSpot(s, mode, tech, out);
	return out;
}

// One FASTQ record of an oracle file.
struct OracleRecord {
	std::string id; // "$ac.$si" or "$ac.$si/$ri"
	std::string length_field;
	std::string bases;
	std::string quality;
};

std::vector<OracleRecord> ReadOracle(const std::string &path) {
	gzFile f = gzopen(path.c_str(), "rb");
	REQUIRE(f != nullptr);
	std::string text;
	char buf[65536];
	int n;
	while ((n = gzread(f, buf, sizeof(buf))) > 0) {
		text.append(buf, size_t(n));
	}
	gzclose(f);
	std::vector<OracleRecord> out;
	size_t pos = 0;
	auto line = [&] {
		size_t nl = text.find('\n', pos);
		REQUIRE(nl != std::string::npos);
		std::string l = text.substr(pos, nl - pos);
		pos = nl + 1;
		return l;
	};
	while (pos < text.size()) {
		OracleRecord r;
		std::string header = line();
		r.id = header.substr(1, header.find(' ') - 1);
		r.length_field = header.substr(header.rfind(' ') + 1);
		r.bases = line();
		line();
		r.quality = line();
		out.push_back(std::move(r));
	}
	return out;
}

std::string QualityText(const std::vector<uint8_t> &q) {
	std::string s;
	for (uint8_t v : q) {
		s += char(v + 33);
	}
	return s;
}

// The accepted fixtures, and the oracle files of each mode as fasterq-dump named them.
const std::vector<std::string> &AcceptedRuns() {
	static const std::vector<std::string> RUNS = {
	    "SRR32151205", "SRR34459993", "ERR16653849", "ERR5009897",  "SRR4424282",  "SRR6016450",
	    "DRR501055",   "SRR4439672",  "SRR24770212", "ERR8523629",  "SRR30810369", "SRR30938008",
	    "ERR932865",   "ERR1160027",  "SRR14704581", "SRR36240868", "DRR880686",   "SRR37135232"};
	return RUNS;
}

// Every read the reader produces for a run, grouped by destination file, in order.
std::map<uint32_t, std::vector<SraRead>> ReadRun(const std::string &run, SplitMode mode, bool tech) {
	SpotReader reader(std::make_shared<const SraTable>(OpenSraFixture(run + ".sra")));
	std::map<uint32_t, std::vector<SraRead>> by_dest;
	Spot spot;
	std::vector<SraRead> reads;
	for (uint64_t i = 0; i < reader.RowCount(); i++) {
		reader.Read(reader.FirstRow() + int64_t(i), spot);
		reads.clear();
		RouteSpot(spot, mode, tech, reads);
		for (auto &r : reads) {
			by_dest[r.destination].push_back(std::move(r));
		}
	}
	return by_dest;
}

void CompareWithOracle(const std::string &run, const std::vector<SraRead> &ours, const std::string &oracle_path,
                       bool with_read_number) {
	INFO(oracle_path);
	auto expected = ReadOracle(oracle_path);
	REQUIRE(ours.size() == expected.size());
	for (size_t i = 0; i < ours.size(); i++) {
		const SraRead &r = ours[i];
		std::string id = run + "." + std::to_string(r.spot);
		if (with_read_number) {
			id += "/" + std::to_string(r.read_number);
		}
		INFO("record " << i << " " << expected[i].id);
		REQUIRE(id == expected[i].id);
		REQUIRE(r.bases == expected[i].bases);
		REQUIRE(QualityText(r.quality) == expected[i].quality);
		REQUIRE("length=" + std::to_string(r.bases.size()) == expected[i].length_field);
	}
}
} // namespace

TEST_CASE("split-3 sends mates of a two-read spot to files 1 and 2", "[sra]") {
	auto reads = Route(MakeSpot("AAAACCCC", {0, 4}, {4, 4}, {BIO, BIO}), SplitMode::Split3);
	REQUIRE(reads.size() == 2);
	CHECK(reads[0].destination == 1);
	CHECK(reads[0].bases == "AAAA");
	CHECK(reads[0].read_number == 1);
	CHECK(reads[1].destination == 2);
	CHECK(reads[1].bases == "CCCC");
}

TEST_CASE("split-3 numbers mates by the reads it writes, not by position", "[sra]") {
	// A technical read first: the biological reads still go to _1 and _2, keeping their
	// positions (2 and 3) as read numbers.
	auto reads = Route(MakeSpot("GGAAAACCCC", {0, 2, 6}, {2, 4, 4}, {TECH, BIO, BIO}), SplitMode::Split3);
	REQUIRE(reads.size() == 2);
	CHECK(reads[0].destination == 1);
	CHECK(reads[0].read_number == 2);
	CHECK(reads[1].destination == 2);
	CHECK(reads[1].read_number == 3);
}

TEST_CASE("split-3 sends a spot with one biological read to the unpaired file", "[sra]") {
	SECTION("the other read is empty") {
		auto reads = Route(MakeSpot("AAAA", {0, 4}, {4, 0}, {BIO, BIO}), SplitMode::Split3);
		REQUIRE(reads.size() == 1);
		CHECK(reads[0].destination == 0);
	}
	SECTION("the other read is technical") {
		auto reads = Route(MakeSpot("GGAAAA", {0, 2}, {2, 4}, {TECH, BIO}), SplitMode::Split3);
		REQUIRE(reads.size() == 1);
		CHECK(reads[0].destination == 0);
		CHECK(reads[0].bases == "AAAA");
	}
	SECTION("technical reads are skipped even when asked for") {
		auto reads = Route(MakeSpot("GGAAAA", {0, 2}, {2, 4}, {TECH, BIO}), SplitMode::Split3, true);
		CHECK(reads.size() == 1);
	}
}

TEST_CASE("split-3 gives a third biological read file 3", "[sra]") {
	auto reads = Route(MakeSpot("AACCGG", {0, 2, 4}, {2, 2, 2}, {BIO, BIO, BIO}), SplitMode::Split3);
	REQUIRE(reads.size() == 3);
	CHECK(reads[2].destination == 3);
}

TEST_CASE("split-3 writes nothing for a spot with no bases", "[sra]") {
	CHECK(Route(MakeSpot("", {0, 0}, {0, 0}, {BIO, BIO}), SplitMode::Split3).empty());
}

TEST_CASE("split-spot writes each read, technical ones only when asked", "[sra]") {
	Spot s = MakeSpot("GGAAAACCCC", {0, 2, 6}, {2, 4, 4}, {TECH, BIO, BIO});
	auto reads = Route(s, SplitMode::Spot);
	REQUIRE(reads.size() == 2);
	CHECK(reads[0].read_number == 2);
	CHECK(reads[0].destination == 0);
	reads = Route(s, SplitMode::Spot, true);
	REQUIRE(reads.size() == 3);
	CHECK(reads[0].bases == "GG");
	CHECK(reads[0].read_number == 1);
}

TEST_CASE("A one-read spot is written whole, as fasterq-dump does", "[sra]") {
	// fasterq-dump prints the spot's entire READ for a one-read spot, not the READ_START/READ_LEN
	// slice, and tests its whole length.
	Spot s = MakeSpot("AAAACC", {0}, {4}, {BIO});
	for (auto mode : {SplitMode::Split3, SplitMode::Spot}) {
		auto reads = Route(s, mode);
		REQUIRE(reads.size() == 1);
		CHECK(reads[0].bases == "AAAACC");
		CHECK(reads[0].read_number == 1);
	}
	CHECK(Route(MakeSpot("GGGG", {0}, {4}, {TECH}), SplitMode::Spot).empty());
	CHECK(Route(MakeSpot("GGGG", {0}, {4}, {TECH}), SplitMode::Spot, true).size() == 1);
}

TEST_CASE("Qualities above 93 are written as 93", "[sra]") {
	// fasterq-dump's table maps every q >= 93 to '~' (93 + 33).
	Spot s = MakeSpot("AC", {0}, {2}, {BIO});
	s.quality = {93, 120};
	auto reads = Route(s, SplitMode::Spot);
	CHECK(reads[0].quality == std::vector<uint8_t> {93, 93});
}

TEST_CASE("A spot whose READ_TYPE disagrees with READ_LEN stops the run", "[sra]") {
	Spot s = MakeSpot("AACC", {0, 2}, {2, 2}, {BIO});
	CHECK(ThrowsTypeMentioning<SraFormatError>([&] { Route(s, SplitMode::Split3); }, "READ_TYPE"));
}

TEST_CASE("Spots of every fixture match fasterq-dump's split-3 output", "[sra]") {
	for (const auto &run : AcceptedRuns()) {
		INFO(run);
		auto by_dest = ReadRun(run, SplitMode::Split3, false);
		std::string dir = "data/sra/oracle/" + run + "/split_3/";
		std::map<uint32_t, std::string> files;
		for (const auto &kv : by_dest) {
			files[kv.first] = dir + run + (kv.first == 0 ? "" : "_" + std::to_string(kv.first)) + ".fastq.gz";
		}
		// Every oracle file must be produced, and no other.
		for (const char *suffix : {"", "_1", "_2", "_3"}) {
			std::string path = dir + run + suffix + ".fastq.gz";
			bool exists = std::ifstream(path).good();
			bool produced = false;
			for (const auto &kv : files) {
				produced |= kv.second == path;
			}
			CHECK(exists == produced);
		}
		for (const auto &kv : by_dest) {
			CompareWithOracle(run, kv.second, files[kv.first], false);
		}
	}
}

TEST_CASE("Spots of every fixture match fasterq-dump's split-spot output", "[sra]") {
	for (bool tech : {false, true}) {
		for (const auto &run : AcceptedRuns()) {
			INFO(run << (tech ? " with technical reads" : ""));
			auto by_dest = ReadRun(run, SplitMode::Spot, tech);
			REQUIRE(by_dest.size() <= 1);
			std::vector<SraRead> reads = by_dest.empty() ? std::vector<SraRead>() : by_dest.begin()->second;
			// --include-technical wrote the same file as --split-spot for every fixture
			// (data/sra/README.md).
			CompareWithOracle(run, reads, "data/sra/oracle/" + run + "/split_spot/" + run + ".fastq.gz", true);
		}
	}
}
