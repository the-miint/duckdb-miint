#include <catch2/catch_test_macros.hpp>

#include "sra_run.hpp"
#include "sra_test_helpers.hpp"

#include <zlib.h>

#include <algorithm>
#include <atomic>
#include <cstring>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using namespace miint::sra;
using miint_test::OpenSraFixture;
using miint_test::SlurpBytes;
using miint_test::SraArchiveFromBytes;
using miint_test::ThrowsTypeMentioning;

// The reader must write exactly the reads fasterq-dump writes, in its order, to the same files,
// under the same spot names. Routing rules follow tbl_join.c (sra-tools 3.4.1); the fixtures and
// oracles are described in data/sra/README.md.

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
	bool has_name = false;
	std::string name; // $sn: between the id and " length=", and may contain spaces
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
		size_t first_space = header.find(' ');
		size_t length_at = header.rfind(' ');
		r.id = header.substr(1, first_space - 1);
		r.length_field = header.substr(length_at + 1);
		r.has_name = length_at > first_space;
		if (r.has_name) {
			r.name = header.substr(first_space + 1, length_at - first_space - 1);
		}
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
		REQUIRE(r.has_name == expected[i].has_name);
		REQUIRE(r.name == expected[i].name);
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
	REQUIRE(reads.size() == 1);
	CHECK(reads[0].quality == std::vector<uint8_t> {93, 93});
}

TEST_CASE("A spot whose READ_TYPE disagrees with READ_LEN stops the run", "[sra]") {
	Spot s = MakeSpot("AACC", {0, 2}, {2, 2}, {BIO});
	CHECK(ThrowsTypeMentioning<SraFormatError>([&] { Route(s, SplitMode::Split3); }, "READ_TYPE"));
}

TEST_CASE("split-3 rows pair a spot's mates and keep an unpaired read alone", "[sra]") {
	std::vector<SraRow> rows;
	auto reads = Route(MakeSpot("AAAACCCC", {0, 4}, {4, 4}, {BIO, BIO}), SplitMode::Split3);
	RowsForSpot(reads, SplitMode::Split3, rows);
	REQUIRE(rows.size() == 1);
	CHECK(rows[0].paired);
	CHECK(rows[0].read1.bases == "AAAA");
	CHECK(rows[0].read2.bases == "CCCC");

	rows.clear();
	reads = Route(MakeSpot("AAAA", {0, 4}, {4, 0}, {BIO, BIO}), SplitMode::Split3);
	RowsForSpot(reads, SplitMode::Split3, rows);
	REQUIRE(rows.size() == 1);
	CHECK_FALSE(rows[0].paired);
	CHECK(rows[0].read1.bases == "AAAA");
}

TEST_CASE("split-3 rows cannot hold a third biological read", "[sra]") {
	// fasterq-dump writes it to a file _3; a read_fastx row has two mates, so the reader refuses
	// rather than drop it, and names the mode that keeps it.
	std::vector<SraRow> rows;
	auto reads = Route(MakeSpot("AACCGG", {0, 2, 4}, {2, 2, 2}, {BIO, BIO, BIO}), SplitMode::Split3);
	CHECK(ThrowsTypeMentioning<SraUnsupportedError>([&] { RowsForSpot(reads, SplitMode::Split3, rows); },
	                                                "split => 'spot'"));
}

TEST_CASE("split-spot rows are one per read", "[sra]") {
	std::vector<SraRow> rows;
	auto reads = Route(MakeSpot("AACCGG", {0, 2, 4}, {2, 2, 2}, {BIO, BIO, BIO}), SplitMode::Spot);
	RowsForSpot(reads, SplitMode::Spot, rows);
	REQUIRE(rows.size() == 3);
	CHECK(rows[2].read_number == 3);
	CHECK_FALSE(rows[2].paired);
}

// A fixture with every occurrence of a name (in its table of contents and its schema text) renamed,
// which removes that object from the table without touching anything else.
namespace {
SraTable WithRenamed(const std::string &run, const std::string &from, const std::string &to) {
	auto bytes = SlurpBytes("data/sra/" + run + ".sra");
	auto at = bytes.begin();
	while ((at = std::search(at, bytes.end(), from.begin(), from.end())) != bytes.end()) {
		at = std::copy(to.begin(), to.end(), at);
	}
	return SraTable(SraArchiveFromBytes(std::move(bytes)));
}
} // namespace

namespace {
// An in-memory archive that counts the bytes read from it.
class CountingSource : public SraByteSource {
public:
	explicit CountingSource(std::vector<uint8_t> bytes) : bytes_(std::move(bytes)) {
	}
	uint64_t Size() const override {
		return bytes_.size();
	}
	void ReadAt(uint64_t pos, uint8_t *buf, size_t n) const override {
		std::memcpy(buf, bytes_.data() + pos, n);
		read_ += n;
	}
	uint64_t BytesRead() const {
		return read_;
	}

private:
	std::vector<uint8_t> bytes_;
	mutable std::atomic<uint64_t> read_ {0};
};
} // namespace

TEST_CASE("Readers of one run share its name index", "[sra]") {
	// Each decoding thread has a reader. A run that stores whole names as keys has an index of
	// tens of bytes per spot, so a reader made for another thread takes the first one's index
	// rather than reading its own copy.
	auto source = std::make_shared<CountingSource>(SlurpBytes("data/sra/SRR24770212.sra"));
	auto t = std::make_shared<const SraTable>(std::make_shared<const KarArchive>(source));
	const uint64_t skey_bytes = t->Archive().Find("tbl/SEQUENCE/idx/skey")->size;
	SpotReader first(t);
	uint64_t before = source->BytesRead();
	SpotReader alone(t);
	const uint64_t unshared = source->BytesRead() - before;
	before = source->BytesRead();
	SpotReader sharing(t, &first);
	const uint64_t shared = source->BytesRead() - before;
	CHECK(unshared - shared == skey_bytes);
	Spot a, b;
	for (int64_t row = first.FirstRow(); row < first.FirstRow() + int64_t(first.RowCount()); row++) {
		first.Read(row, a);
		sharing.Read(row, b);
		REQUIRE(a.name == b.name);
	}
}

TEST_CASE("A name format is only taken through the skey index", "[sra]") {
	// The schemas produce the format with idx:text:project<'skey'>(.NAME_FMT), which VDB cannot
	// build without the index: the name then falls back to the stored spot name, even when
	// NAME_FMT holds a format. (Every public run surveyed without the index had an empty NAME_FMT,
	// so either reading gives the same names there.)
	auto t = std::make_shared<const SraTable>(WithRenamed("SRR34459993", "skey", "skez"));
	REQUIRE(t->Static("SPOT_NAME") != nullptr);
	REQUIRE(t->Archive().Find("tbl/SEQUENCE/idx/skey") == nullptr);
	SpotReader reader(t);
	Spot spot;
	reader.Read(1, spot);
	CHECK(spot.has_name);
	CHECK(spot.name == t->Static("SPOT_NAME")->row);
}

TEST_CASE("An skey index without a NAME_FMT column is refused", "[sra]") {
	// The older skeyname schema, which projects with no column to fall back on and fails a row the
	// index has no key for; no surveyed public run has it.
	auto t = std::make_shared<const SraTable>(WithRenamed("SRR32151205", "NAME_FMT", "NAME_FMU"));
	REQUIRE_FALSE(t->HasPhysical("NAME_FMT"));
	CHECK(ThrowsTypeMentioning<SraUnsupportedError>([&] { SpotReader reader(t); }, "NAME_FMT"));
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
