// Tests for threading ENA's fastq_md5 field through to ENARunInfo.
//
// This is the blocker half of md5 verification: without fastq_md5
// present in BATCH_FIELDS (src/ena_resolver_cache.cpp) and correctly parsed /
// aligned here, read_ena_sequences would silently never have a digest to
// verify against, no matter how the streaming tap is wired up.
//
// NOTE on scope: this file intentionally does NOT unit test miint::StreamMd5
// or PerRunReader (src/include/stream_md5.hpp, src/per_run_reader.cpp)
// directly. Both use duckdb::MD5Context / duckdb::IOException / a real
// duckdb::FileSystem, none of which are linkable into this Catch2 binary
// without pulling in a large fraction of DuckDB's query engine (this target
// deliberately never links duckdb -- see the "tests" target's CMakeLists.txt
// comments). The sibling md5 pattern this mirrors, GzipMd5Stream in
// ena_upload_reads.cpp, has the same property and is likewise untested at
// this level; both are instead validated at the SQL level (see
// test/sql/ena_upload_reads_local.test for the write-side precedent and
// test/sql/read_ena_sequences.test for the read-side coverage added here).
#include "ena_run_info_extractor.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>

namespace {

miint::ENAColumnIndexMap StandardColsWithMd5() {
	std::vector<std::string> header = {"run_accession",    "sample_accession", "experiment_accession",
	                                   "fastq_ftp",        "fastq_aspera",     "fastq_bytes",
	                                   "fastq_md5",        "library_layout",   "submitted_ftp",
	                                   "submitted_aspera", "submitted_bytes",  "submitted_format"};
	return miint::ENAColumnIndexMap::FromHeader(header);
}

} // namespace

TEST_CASE("ENAColumnIndexMap maps fastq_md5", "[ena_md5]") {
	std::vector<std::string> header = {"run_accession", "fastq_ftp", "fastq_md5"};
	auto cols = miint::ENAColumnIndexMap::FromHeader(header);
	CHECK(cols.fastq_md5 == 2);
}

TEST_CASE("ENAColumnIndexMap leaves fastq_md5 at -1 when the column is absent", "[ena_md5]") {
	std::vector<std::string> header = {"run_accession", "fastq_ftp"};
	auto cols = miint::ENAColumnIndexMap::FromHeader(header);
	CHECK(cols.fastq_md5 == -1);
}

TEST_CASE("ENARunInfoExtractor threads a single-end fastq_md5", "[ena_md5]") {
	auto cols = StandardColsWithMd5();
	std::vector<std::string> row = {"ERR1074767",
	                                "SAMEA3271045",
	                                "ERX1106717",
	                                "ftp.sra.ebi.ac.uk/vol1/fastq/ERR107/007/ERR1074767/ERR1074767.fastq.gz",
	                                "",
	                                "12345",
	                                "d41d8cd98f00b204e9800998ecf8427e",
	                                "SINGLE",
	                                "",
	                                "",
	                                "",
	                                ""};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "auto");
	REQUIRE(runs.size() == 1);
	REQUIRE(runs[0].fastq_md5.size() == 1);
	CHECK(runs[0].fastq_md5[0] == "d41d8cd98f00b204e9800998ecf8427e");
}

TEST_CASE("ENARunInfoExtractor threads paired-end fastq_md5 aligned with fastq_urls", "[ena_md5]") {
	auto cols = StandardColsWithMd5();
	std::vector<std::string> row = {"ERR1000",
	                                "SAM1",
	                                "EXP1",
	                                "ftp.example.com/ERR1000_1.fastq.gz;ftp.example.com/ERR1000_2.fastq.gz",
	                                "",
	                                "100;200",
	                                "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa;bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
	                                "PAIRED",
	                                "",
	                                "",
	                                "",
	                                ""};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "auto");
	REQUIRE(runs.size() == 1);
	REQUIRE(runs[0].fastq_urls.size() == 2);
	REQUIRE(runs[0].fastq_md5.size() == 2);
	CHECK(runs[0].fastq_md5[0] == "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
	CHECK(runs[0].fastq_md5[1] == "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb");
}

TEST_CASE("ENARunInfoExtractor re-indexes fastq_md5 through the 3-file paired-end filter", "[ena_md5]") {
	auto cols = StandardColsWithMd5();
	// Same orphan-file shape as the existing fastq_bytes filter test
	// (test_ENAClient.cpp "filters 3-file paired-end to _1/_2"): an orphan
	// single file plus _1 and _2. fastq_md5 must survive the same fi-indexed
	// filtering as fastq_bytes, staying aligned with the filtered fastq_urls
	// (dropping the orphan's digest, keeping _1's and _2's).
	std::vector<std::string> row = {
	    "ERR2000",
	    "SAM2",
	    "EXP2",
	    "ftp.example.com/ERR2000.fastq.gz;ftp.example.com/ERR2000_1.fastq.gz;ftp.example.com/ERR2000_2.fastq.gz",
	    "",
	    "50;100;200",
	    "orphandigest00000000000000000000;aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa;bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
	    "PAIRED",
	    "",
	    "",
	    "",
	    ""};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "auto");
	REQUIRE(runs.size() == 1);
	REQUIRE(runs[0].fastq_urls.size() == 2);
	CHECK(runs[0].fastq_urls[0].find("_1.fast") != std::string::npos);
	CHECK(runs[0].fastq_urls[1].find("_2.fast") != std::string::npos);
	REQUIRE(runs[0].fastq_md5.size() == 2);
	CHECK(runs[0].fastq_md5[0] == "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
	CHECK(runs[0].fastq_md5[1] == "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb");
}

TEST_CASE("ENARunInfoExtractor leaves fastq_md5 empty when ENA omits the field", "[ena_md5]") {
	auto cols = StandardColsWithMd5();
	std::vector<std::string> row = {"ERR3000", "SAM3", "EXP3", "ftp.example.com/ERR3000.fastq.gz",
	                                "",        "999",
	                                "", // fastq_md5 empty — ENA didn't report one
	                                "SINGLE",  "",     "",     "",
	                                ""};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "auto");
	REQUIRE(runs.size() == 1);
	CHECK(runs[0].fastq_md5.empty());
}

TEST_CASE("ENARunInfoExtractor degrades to empty fastq_md5 on a count mismatch rather than misaligning", "[ena_md5]") {
	auto cols = StandardColsWithMd5();
	// Two files, but only one md5 token — a malformed/inconsistent ENA
	// response. Must not throw, and must not pair the single token with
	// either file (that would silently attribute the wrong digest to a
	// file). Degrades to "no md5 for this run".
	std::vector<std::string> row = {"ERR4000",
	                                "SAM4",
	                                "EXP4",
	                                "ftp.example.com/ERR4000_1.fastq.gz;ftp.example.com/ERR4000_2.fastq.gz",
	                                "",
	                                "100;200",
	                                "onlyonedigest0000000000000000000",
	                                "PAIRED",
	                                "",
	                                "",
	                                "",
	                                ""};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "auto");
	REQUIRE(runs.size() == 1);
	REQUIRE(runs[0].fastq_urls.size() == 2);
	CHECK(runs[0].fastq_md5.empty());
}
