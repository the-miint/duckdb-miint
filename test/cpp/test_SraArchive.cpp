#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "sra_archive.hpp"
#include "sra_test_helpers.hpp"

#include <cstring>
#include <memory>
#include <string>
#include <vector>

using namespace miint::sra;
using Catch::Matchers::ContainsSubstring;

// The fixtures and their provenance are described in data/sra/README.md. Expected values below
// were read from the files with an independent parser written while the format was worked out.
using miint_test::OpenSraFixture;
using miint_test::SlurpBytes;
using miint_test::SraArchiveFromBytes;
using miint_test::ThrowsTypeMentioning;

TEST_CASE("KarArchive reads the header and table of contents", "[sra]") {
	auto ar = OpenSraFixture("SRR4424282.sra");
	CHECK(ar->DataOffset() == 920);
	CHECK(ar->Entries().size() == 28);

	// A file's offset is relative to the data region in the TOC; Find reports it absolute, so a
	// wrong base would read some other file's bytes.
	const KarEntry *read_data = ar->Find("tbl/SEQUENCE/col/READ/data");
	REQUIRE(read_data != nullptr);
	CHECK(read_data->type == KarEntry::Type::File);
	CHECK(read_data->offset == 2612);
	CHECK(read_data->size == 482);
	CHECK(ar->ReadFile("tbl/SEQUENCE/col/READ/data").size() == 482);

	CHECK(ar->IsDirectory("tbl"));
	CHECK(ar->Children("tbl") == std::vector<std::string> {"SEQUENCE"});
	CHECK(ar->Find("lock")->type == KarEntry::Type::EmptyFile);
	CHECK(ar->Find("no/such/file") == nullptr);
}

TEST_CASE("KarArchive rejects bytes that are not an SRA archive", "[sra]") {
	std::vector<uint8_t> junk(64, 'x');
	CHECK(ThrowsTypeMentioning<SraFormatError>([&] { SraArchiveFromBytes(junk); }, "not an SRA archive"));
}

TEST_CASE("KarArchive reports a truncated table of contents as damage", "[sra]") {
	auto bytes = SlurpBytes("data/sra/SRR4424282.sra");
	bytes.resize(100); // the TOC runs to byte 920
	CHECK_THROWS_AS(SraArchiveFromBytes(bytes), SraFormatError);
}

TEST_CASE("KarArchive rejects archive variants absent from public runs", "[sra]") {
	auto bytes = SlurpBytes("data/sra/SRR4424282.sra");
	SECTION("big-endian") {
		// The byte-order tag as a big-endian writer stores it.
		bytes[8] = 0x05;
		bytes[9] = 0x03;
		bytes[10] = 0x19;
		bytes[11] = 0x88;
		CHECK(ThrowsTypeMentioning<SraUnsupportedError>([&] { SraArchiveFromBytes(bytes); }, "big-endian"));
	}
	SECTION("unknown byte-order tag") {
		bytes[8] = 0x00;
		CHECK_THROWS_AS(SraArchiveFromBytes(bytes), SraFormatError);
	}
	SECTION("archive version 2") {
		bytes[12] = 0x02;
		CHECK(ThrowsTypeMentioning<SraUnsupportedError>([&] { SraArchiveFromBytes(bytes); }, "version 2"));
	}
}

TEST_CASE("Metadata trees parse names, attributes, children and values", "[sra]") {
	auto ar = OpenSraFixture("SRR4424282.sra");
	MetadataNode md = ParseMetadata(ar->ReadFile("tbl/SEQUENCE/md/cur"));
	const MetadataNode *schema = md.Find("schema");
	REQUIRE(schema != nullptr);
	CHECK(schema->attrs.at("name") == "NCBI:align:tbl:seq#1.1");
	CHECK_THAT(schema->value, ContainsSubstring("NCBI:align:tbl:seq"));
	// SPOT_COUNT is stored in 8 bytes here; KDB writes counts in 1, 2, 4 or 8 bytes.
	REQUIRE(md.Find("STATS/TABLE/SPOT_COUNT") != nullptr);
	CHECK(md.Find("STATS/TABLE/SPOT_COUNT")->ValueAsU64() == 3);
	CHECK(md.Find("STATS/TABLE/NO_SUCH") == nullptr);
}

TEST_CASE("MetadataNode reads integers of every width KDB writes", "[sra]") {
	MetadataNode n;
	n.value = std::string("\x94", 1);
	CHECK(n.ValueAsU64() == 148);
	n.value = std::string("\x50\x2a", 2);
	CHECK(n.ValueAsU64() == 0x2a50);
	n.value = std::string("\x50\x2a\x01\x00", 4);
	CHECK(n.ValueAsU64() == 0x12a50);
	n.value = std::string("\x03\x00\x00\x00\x00\x00\x00\x00", 8);
	CHECK(n.ValueAsU64() == 3);
	n.value = std::string("\xff", 1);
	CHECK(n.ValueAsI64() == -1);
	// Three bytes is no integer width KDB uses: reading it as one would invent a value.
	n.value = std::string("\x01\x02\x03", 3);
	CHECK_THROWS_AS(n.ValueAsU64(), SraFormatError);
}

TEST_CASE("SraTable finds the SEQUENCE table of a database", "[sra]") {
	SraTable t(OpenSraFixture("SRR4424282.sra"));
	CHECK(t.InDatabase());
	CHECK(t.DatabaseTables() == std::vector<std::string> {"SEQUENCE"});
	CHECK(t.Prefix() == "tbl/SEQUENCE/");
	CHECK(t.TableType() == "NCBI:align:tbl:seq#1.1");
	CHECK(t.PhysicalColumns() == std::vector<std::string> {"QUALITY", "READ"});
	REQUIRE(t.DatabaseMetadata() != nullptr);
	CHECK(t.DatabaseMetadata()->Find("schema")->attrs.at("name") == "NCBI:align:db:alignment_sorted#1.3");

	// READ_LEN is the same for all three spots, so it is static: two reads of 250 bases.
	const StaticColumn *read_len = t.Static("READ_LEN");
	REQUIRE(read_len != nullptr);
	CHECK(read_len->type == "INSDC:coord:len");
	CHECK(read_len->row == std::string("\xfa\x00\x00\x00\xfa\x00\x00\x00", 8));
	CHECK(read_len->start_id == 1);
	CHECK(read_len->row_count == 3);
	CHECK_FALSE(t.HasPhysical("READ_LEN"));
	CHECK(t.Static("READ") == nullptr);
}

TEST_CASE("SraTable reads a flat table, which has no database around it", "[sra]") {
	// Flat tables are not only old runs (data/sra/README.md); missing them would refuse ~7% of
	// the runs surveyed.
	SraTable t(OpenSraFixture("ERR932865.sra"));
	CHECK_FALSE(t.InDatabase());
	CHECK(t.DatabaseTables().empty());
	CHECK(t.Prefix().empty());
	CHECK(t.DatabaseMetadata() == nullptr);
	CHECK(t.TableType() == "NCBI:SRA:Illumina:tbl:phred:v2#1.0.4");
	CHECK(t.HasPhysical("READ"));
	REQUIRE(t.Static("READ_LEN") != nullptr);
	CHECK(t.Static("READ_LEN")->row_count == 148);
}

TEST_CASE("SraTable lists every table of a database, not only SEQUENCE", "[sra]") {
	// Classification decides on the tables present, so they must all be reported.
	SraTable t(OpenSraFixture("SRR26779001.sra"));
	CHECK(t.DatabaseTables() == std::vector<std::string> {"PRIMARY_ALIGNMENT", "REFERENCE", "SEQUENCE"});
}

TEST_CASE("SraTable reads a run whose every column is static", "[sra]") {
	// A one-spot run stores even READ and its qualities in table metadata.
	SraTable t(OpenSraFixture("SRR30938008.sra"));
	CHECK(t.PhysicalColumns().empty());
	REQUIRE(t.Static("READ") != nullptr);
	REQUIRE(t.Static("ORIGINAL_QUALITY") != nullptr);
	CHECK(t.Static("READ")->row_count == 1);
}

TEST_CASE("An archive's table of contents is refused past a sane size", "[sra]") {
	// Public archives' tables of contents are a few KB; a header placing the data 64 MiB in is not
	// read whole into memory.
	const size_t size = (size_t(64) << 20) + 64;
	std::vector<uint8_t> bytes(size, 0);
	std::memcpy(bytes.data(), "NCBI.sra", 8);
	const uint32_t order = 0x05031988, version = 1;
	std::memcpy(bytes.data() + 8, &order, 4);
	std::memcpy(bytes.data() + 12, &version, 4);
	const uint64_t data_offset = size - 8;
	std::memcpy(bytes.data() + 16, &data_offset, 8);
	CHECK(ThrowsTypeMentioning<SraUnsupportedError>([&] { SraArchiveFromBytes(bytes); }, "table of contents"));
}
