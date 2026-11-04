#include <catch2/catch_test_macros.hpp>

#include "sra_run.hpp"
#include "sra_test_helpers.hpp"

#include <string>
#include <vector>

using namespace miint::sra;
using miint_test::OpenSraFixture;
using miint_test::ThrowsTypeMentioning;

// Which runs miint reads is decided by the tables present, the table type, the Lite markers and
// the column encodings: never by platform (Daniel, 2026-11-04). fasterq-dump reads PacBio and
// Nanopore runs loaded as ordinary sequence tables exactly as it reads Illumina runs, and nearly
// every current Nanopore run is one, so a platform rule would refuse runs that decode
// identically. data/sra/README.md describes the fixtures.

namespace {
RunFacts FactsOf(const std::string &file) {
	SraTable t(OpenSraFixture(file));
	return GatherRunFacts(t);
}

// A minimal accepted run: Illumina phred v2 database with READ and ORIGINAL_QUALITY.
RunFacts IlluminaFacts() {
	RunFacts f;
	f.in_database = true;
	f.database_tables = {"SEQUENCE"};
	f.table_type = "NCBI:SRA:Illumina:tbl:phred:v2#2";
	f.physical_columns = {"READ", "ORIGINAL_QUALITY", "READ_LEN", "READ_START"};
	f.encodings = {{"READ", ""},
	               {"ORIGINAL_QUALITY", "<INSDC:quality:phred>zip_encoding#1"},
	               {"READ_LEN", "<INSDC:coord:len>izip_encoding#1"},
	               {"READ_START", "<INSDC:coord:zero>izip_encoding#1"}};
	return f;
}
} // namespace

TEST_CASE("Every decodable fixture is accepted", "[sra]") {
	// The 15 runs in the README's accepted table and the 3 PacBio/Nanopore runs: flat tables,
	// every supported table type, PacBio HiFi and both Nanopore loads.
	const std::vector<std::string> accepted = {
	    "SRR32151205.sra", "SRR34459993.sra", "ERR16653849.sra", "ERR5009897.sra",  "SRR4424282.sra",
	    "SRR6016450.sra",  "DRR501055.sra",   "SRR4439672.sra",  "SRR24770212.sra", "ERR8523629.sra",
	    "SRR30810369.sra", "SRR30938008.sra", "ERR932865.sra",   "ERR1160027.sra",  "SRR14704581.sra",
	    "SRR36240868.sra", "DRR880686.sra",   "SRR37135232.sra"};
	for (const auto &file : accepted) {
		INFO(file);
		CHECK_NOTHROW(ClassifyRun(FactsOf(file)));
	}
}

TEST_CASE("Rejection fixtures are refused, each named for what it is", "[sra]") {
	CHECK(ThrowsTypeMentioning<SraUnsupportedError>([] { ClassifyRun(FactsOf("SRR26779001.sra")); }, "cSRA"));
	CHECK(ThrowsTypeMentioning<SraUnsupportedError>([] { ClassifyRun(FactsOf("SRR2513321.sra")); }, "CONSENSUS"));
	CHECK(ThrowsTypeMentioning<SraUnsupportedError>([] { ClassifyRun(FactsOf("SRR882123.sra")); },
	                                                "NCBI:SRA:_454_:tbl:v2#1.0.7"));
	CHECK(ThrowsTypeMentioning<SraUnsupportedError>([] { ClassifyRun(FactsOf("SRR4424282.sralite")); }, "SRA Lite"));
}

TEST_CASE("GatherRunFacts reports the Lite marker in table metadata", "[sra]") {
	// delite stamps SOFTWARE/delite into the table's metadata (checked on the committed .sralite)
	CHECK(FactsOf("SRR4424282.sralite").delite);
	CHECK_FALSE(FactsOf("SRR4424282.sra").delite);
}

TEST_CASE("A PacBio run counts as native only by its tables", "[sra]") {
	SECTION("PacBio HiFi as a sequence table is read") {
		CHECK_NOTHROW(ClassifyRun(IlluminaFacts()));
	}
	for (const char *table : {"CONSENSUS", "ZMW_METRICS", "PASSES"}) {
		SECTION(table) {
			RunFacts f = IlluminaFacts();
			f.database_tables.push_back(table);
			CHECK(ThrowsTypeMentioning<SraUnsupportedError>([&] { ClassifyRun(f); }, table));
		}
	}
}

TEST_CASE("Alignment data is cSRA however it shows", "[sra]") {
	SECTION("alignment tables") {
		RunFacts f = IlluminaFacts();
		f.database_tables = {"PRIMARY_ALIGNMENT", "REFERENCE", "SEQUENCE"};
		CHECK(ThrowsTypeMentioning<SraUnsupportedError>([&] { ClassifyRun(f); }, "cSRA"));
	}
	SECTION("reads stored against alignments (CMP_READ)") {
		RunFacts f = IlluminaFacts();
		f.physical_columns.insert("CMP_READ");
		CHECK(ThrowsTypeMentioning<SraUnsupportedError>([&] { ClassifyRun(f); }, "cSRA"));
	}
	SECTION("a physical PRIMARY_ALIGNMENT_ID") {
		RunFacts f = IlluminaFacts();
		f.physical_columns.insert("PRIMARY_ALIGNMENT_ID");
		CHECK(ThrowsTypeMentioning<SraUnsupportedError>([&] { ClassifyRun(f); }, "cSRA"));
	}
	SECTION("a static PRIMARY_ALIGNMENT_ID naming an alignment") {
		RunFacts f = IlluminaFacts();
		f.static_columns.insert("PRIMARY_ALIGNMENT_ID");
		f.static_primary_alignment_id = std::string("\x07\x00\x00\x00\x00\x00\x00\x00", 8);
		CHECK(ThrowsTypeMentioning<SraUnsupportedError>([&] { ClassifyRun(f); }, "cSRA"));
	}
	SECTION("a static PRIMARY_ALIGNMENT_ID of zeros means unaligned") {
		// latf-load writes this for runs with no alignments at all.
		RunFacts f = IlluminaFacts();
		f.static_columns.insert("PRIMARY_ALIGNMENT_ID");
		f.static_primary_alignment_id = std::string(16, '\0');
		CHECK_NOTHROW(ClassifyRun(f));
	}
}

TEST_CASE("A run without original qualities is SRA Lite", "[sra]") {
	SECTION("marked by delite") {
		RunFacts f = IlluminaFacts();
		f.delite = true;
		CHECK(ThrowsTypeMentioning<SraUnsupportedError>([&] { ClassifyRun(f); }, "SRA Lite"));
	}
	SECTION("unmarked but with no quality column") {
		// Its qualities would be synthesized (all 30), which fasterq-dump prints without warning.
		RunFacts f = IlluminaFacts();
		f.physical_columns.erase("ORIGINAL_QUALITY");
		f.encodings.erase("ORIGINAL_QUALITY");
		CHECK(ThrowsTypeMentioning<SraUnsupportedError>([&] { ClassifyRun(f); }, "SRA Lite"));
	}
	SECTION("a static quality column is still a quality column") {
		RunFacts f = IlluminaFacts();
		f.physical_columns.erase("ORIGINAL_QUALITY");
		f.encodings.erase("ORIGINAL_QUALITY");
		f.static_columns.insert("QUALITY");
		CHECK_NOTHROW(ClassifyRun(f));
	}
}

TEST_CASE("Table types are read only from the list that has been checked", "[sra]") {
	for (const char *type : {"NCBI:SRA:Illumina:tbl:phred:v2#1.0.4", "NCBI:SRA:GenericFastq:sequence#1",
	                         "NCBI:SRA:GenericFastq:sequence_nanopore#2", "NCBI:align:tbl:seq#2"}) {
		INFO(type);
		RunFacts f = IlluminaFacts();
		f.table_type = type;
		CHECK_NOTHROW(ClassifyRun(f));
	}
	for (const char *type : {"NCBI:SRA:Nanopore:sequence#1", "NCBI:SRA:GenericFastq:sequence_log_odds#2",
	                         "NCBI:SRA:IonTorrent:tbl:v2#1.0.3", "NCBI:sra:tbl:trace#1", ""}) {
		INFO(type);
		RunFacts f = IlluminaFacts();
		f.table_type = type;
		CHECK_THROWS_AS(ClassifyRun(f), SraUnsupportedError);
	}
}

TEST_CASE("An encoding outside the supported set is refused on a column that is read", "[sra]") {
	SECTION("seen in a public run: QUALITY with delta_average_zip_encoding (ERR263256)") {
		RunFacts f = IlluminaFacts();
		f.physical_columns.insert("QUALITY");
		f.encodings["QUALITY"] = "<INSDC:quality:phred>delta_average_zip_encoding#1";
		CHECK(ThrowsTypeMentioning<SraUnsupportedError>([&] { ClassifyRun(f); }, "delta_average_zip_encoding#1"));
	}
	SECTION("spacing and zip parameters are not part of the encoding's name") {
		RunFacts f = IlluminaFacts();
		f.physical_columns.insert("SPOT_NAME");
		f.encodings["SPOT_NAME"] = "< ascii > zip_encoding#1<0,1>";
		CHECK_NOTHROW(ClassifyRun(f));
	}
	SECTION("a column the reader never reads may use anything") {
		RunFacts f = IlluminaFacts();
		f.physical_columns.insert("SIGNAL");
		f.encodings["SIGNAL"] = "NCBI:SRA:_454_:encoding:SIGNAL#2";
		CHECK_NOTHROW(ClassifyRun(f));
	}
}

TEST_CASE("A database without a SEQUENCE table is refused", "[sra]") {
	RunFacts f = IlluminaFacts();
	f.database_tables = {"REFERENCE"};
	CHECK_THROWS_AS(ClassifyRun(f), SraUnsupportedError);
}
