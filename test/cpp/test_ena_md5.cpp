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
// The mismatch message text itself IS covered here, via the duckdb-free
// miint::BuildMd5MismatchMessage builder (read_ena_sequences_policy.hpp)
// that VerifyOrThrow delegates to.
#include "ena_run_info_extractor.hpp"
#include "read_ena_sequences_policy.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cctype>
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

TEST_CASE("ShouldSkipFailedRun skips a failed run only when siblings exist", "[ena_md5]") {
	// A run that fails to open (after its retry), or fails its md5 or ascp exit
	// check at completion, must not abort a multi-run scan and discard every
	// sibling already downloaded -- the
	// varchar[] / project-accession case that would otherwise make bulk fetches
	// unusable. A single-run scan keeps the hard throw: no sibling to protect,
	// and a downstream that resolves a scalar accession one run per call (the
	// Qiita ENA ingest job) relies on that error to fail the run.
	CHECK_FALSE(miint::ShouldSkipFailedRun(1)); // scalar / single run -> throw
	CHECK(miint::ShouldSkipFailedRun(2));       // varchar[] siblings -> skip the bad one
	CHECK(miint::ShouldSkipFailedRun(100));     // project expansion -> skip the bad one
}

TEST_CASE("EveryRunFailed holds only when every run of the scan was skipped", "[ena_md5]") {
	CHECK(miint::EveryRunFailed(1, 1));
	CHECK(miint::EveryRunFailed(3, 3));
	// A sibling's rows landed, so the scan returns them with a warning.
	CHECK_FALSE(miint::EveryRunFailed(2, 3));
	// A scan of no runs failed nothing.
	CHECK_FALSE(miint::EveryRunFailed(0, 0));
	CHECK_FALSE(miint::EveryRunFailed(0, 2));
}

TEST_CASE("BuildEveryRunFailedMessage names each run and the first failure", "[ena_md5]") {
	// The order is the order the runs failed in, which parallel workers decide.
	CHECK(miint::BuildEveryRunFailedMessage({"ERR2", "ERR1"}, "read_ena_sequences: ascp exited with code 1") ==
	      "read_ena_sequences: every run failed (2): ERR2 ERR1; first error: read_ena_sequences: ascp exited with "
	      "code 1");
}

TEST_CASE("BuildMd5MismatchMessage names the label and attributes each digest", "[ena_md5]") {
	auto message = miint::BuildMd5MismatchMessage("ERR1074767", "d41d8cd98f00b204e9800998ecf8427e",
	                                              "5eb63bbbe01eeed093cb22bb8f5acdc3");
	CHECK(message.find("md5 mismatch") != std::string::npos);
	CHECK(message.find("ERR1074767") != std::string::npos);
	auto expected_pos = message.find("ENA reported d41d8cd98f00b204e9800998ecf8427e");
	auto actual_pos = message.find("downloaded bytes hash to 5eb63bbbe01eeed093cb22bb8f5acdc3");
	CHECK(expected_pos != std::string::npos);
	CHECK(actual_pos != std::string::npos);
}

TEST_CASE("BuildMd5MismatchMessage never uses transport-failure wording", "[ena_md5]") {
	// This word list is miint's own published contract (docs/insdc_ena.md,
	// "Message contract"), not a quotation of any downstream classifier.
	static const std::vector<std::string> kExcludedWords = {
	    "connection", "timed out", "timeout", "network", "reset", "refused", "unreachable", "temporarily", "curl"};
	auto message = miint::BuildMd5MismatchMessage(
	    "ERR1074767 https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR107/007/ERR1074767/ERR1074767.fastq.gz",
	    "d41d8cd98f00b204e9800998ecf8427e", "5eb63bbbe01eeed093cb22bb8f5acdc3");
	std::string lowered = message;
	std::transform(lowered.begin(), lowered.end(), lowered.begin(), [](unsigned char c) { return std::tolower(c); });
	for (const auto &word : kExcludedWords) {
		CHECK(lowered.find(word) == std::string::npos);
	}
}

TEST_CASE("BuildMd5MismatchMessage does not transpose expected and actual", "[ena_md5]") {
	auto message = miint::BuildMd5MismatchMessage("ERR1074767", "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
	                                              "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb");
	auto ena_reported_pos = message.find("ENA reported");
	auto expected_digest_pos = message.find("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
	auto downloaded_pos = message.find("downloaded bytes hash to");
	auto actual_digest_pos = message.find("bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb");
	REQUIRE(ena_reported_pos != std::string::npos);
	REQUIRE(expected_digest_pos != std::string::npos);
	REQUIRE(downloaded_pos != std::string::npos);
	REQUIRE(actual_digest_pos != std::string::npos);
	CHECK(ena_reported_pos < expected_digest_pos);
	CHECK(expected_digest_pos < downloaded_pos);
	CHECK(downloaded_pos < actual_digest_pos);
}
