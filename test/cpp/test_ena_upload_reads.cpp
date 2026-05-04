// SPDX-License-Identifier: MIT
//
// Phase 6 RED/GREEN: pure-data unit tests for the helpers behind
// `ena_upload_reads`. Covers:
//   - Per-sample layout detection (single / paired / paired-interleaved /
//     mixed-error)
//   - Filename emission per layout
//   - Upload target URL parsing (aspera:// vs file://)
//   - `ascp --mode=send` argv construction
//
// All assertions are deterministic and require no network, no Aspera binary,
// and no DuckDB library — these helpers are explicitly DuckDB-free for that
// reason. The encoder tests in test_fastq_encoder.cpp cover the byte-exact
// FASTQ side; the wider integration (encode → gzip → MD5 → transport) is
// exercised by test/sql/ena_upload_reads_local.test (file://) and
// test/sql/ena_upload_reads.test (live Aspera, gated on credentials).

#include "ena_upload_helpers.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <string>
#include <vector>

using duckdb::AsperaSendOptions;
using duckdb::BuildAscpSendArgv;
using duckdb::FastqLayoutMode;
using duckdb::FastqLayoutModeName;
using duckdb::OutputFilenames;
using duckdb::ParseFastqLayoutMode;
using duckdb::ParseUploadTargetURL;
using duckdb::ResolveLayout;
using duckdb::UploadTargetURL;
using duckdb::UploadTransport;

// ---- Layout detection ----------------------------------------------------

TEST_CASE("ResolveLayout: AUTO collapses to SINGLE when no row has R2", "[ena_upload_reads]") {
	std::vector<bool> has_r2 = {false, false, false};
	REQUIRE(ResolveLayout("sampleA", FastqLayoutMode::AUTO, has_r2) == FastqLayoutMode::SINGLE);
}

TEST_CASE("ResolveLayout: AUTO collapses to PAIRED when every row has R2", "[ena_upload_reads]") {
	std::vector<bool> has_r2 = {true, true, true};
	REQUIRE(ResolveLayout("sampleA", FastqLayoutMode::AUTO, has_r2) == FastqLayoutMode::PAIRED);
}

TEST_CASE("ResolveLayout: mixed rows raise a clear error naming sample and row", "[ena_upload_reads]") {
	std::vector<bool> has_r2 = {true, true, false, true};
	REQUIRE_THROWS_WITH(ResolveLayout("sampleA", FastqLayoutMode::AUTO, has_r2),
	                    Catch::Matchers::ContainsSubstring("sampleA") && Catch::Matchers::ContainsSubstring("row 2"));

	std::vector<bool> has_r2_other = {false, false, true};
	REQUIRE_THROWS_WITH(ResolveLayout("anotherSample", FastqLayoutMode::AUTO, has_r2_other),
	                    Catch::Matchers::ContainsSubstring("anotherSample") &&
	                        Catch::Matchers::ContainsSubstring("row 2"));
}

TEST_CASE("ResolveLayout: explicit SINGLE rejects rows that supplied R2", "[ena_upload_reads]") {
	std::vector<bool> has_r2 = {false, true, false};
	REQUIRE_THROWS_WITH(ResolveLayout("sampleA", FastqLayoutMode::SINGLE, has_r2),
	                    Catch::Matchers::ContainsSubstring("layout=single") &&
	                        Catch::Matchers::ContainsSubstring("sampleA") &&
	                        Catch::Matchers::ContainsSubstring("row 1"));
}

TEST_CASE("ResolveLayout: SINGLE error reports first row with R2, even when row 0 has R2", "[ena_upload_reads]") {
	// Regression: an earlier implementation reused the AUTO "first mismatch
	// vs row 0" counter and reported the wrong row when the first row was
	// the offender.
	std::vector<bool> has_r2 = {true, false, true};
	REQUIRE_THROWS_WITH(ResolveLayout("sampleA", FastqLayoutMode::SINGLE, has_r2),
	                    Catch::Matchers::ContainsSubstring("layout=single") &&
	                        Catch::Matchers::ContainsSubstring("row 0"));
}

TEST_CASE("ResolveLayout: PAIRED error reports first row missing R2, even when row 0 has R2", "[ena_upload_reads]") {
	std::vector<bool> has_r2 = {true, true, false, true};
	REQUIRE_THROWS_WITH(ResolveLayout("sampleA", FastqLayoutMode::PAIRED, has_r2),
	                    Catch::Matchers::ContainsSubstring("layout=paired") &&
	                        Catch::Matchers::ContainsSubstring("row 2"));
}

TEST_CASE("ResolveLayout: explicit PAIRED rejects rows missing R2", "[ena_upload_reads]") {
	std::vector<bool> has_r2 = {true, true, false};
	REQUIRE_THROWS_WITH(ResolveLayout("sampleA", FastqLayoutMode::PAIRED, has_r2),
	                    Catch::Matchers::ContainsSubstring("layout=paired") &&
	                        Catch::Matchers::ContainsSubstring("sampleA"));
}

TEST_CASE("ResolveLayout: PAIRED_INTERLEAVED requires every row to have R2", "[ena_upload_reads]") {
	std::vector<bool> all_paired = {true, true};
	REQUIRE(ResolveLayout("sampleA", FastqLayoutMode::PAIRED_INTERLEAVED, all_paired) ==
	        FastqLayoutMode::PAIRED_INTERLEAVED);

	std::vector<bool> mixed = {true, false};
	REQUIRE_THROWS_WITH(ResolveLayout("sampleA", FastqLayoutMode::PAIRED_INTERLEAVED, mixed),
	                    Catch::Matchers::ContainsSubstring("paired_interleaved"));
}

TEST_CASE("ResolveLayout: empty rows is treated as a programming error", "[ena_upload_reads]") {
	std::vector<bool> empty;
	REQUIRE_THROWS(ResolveLayout("sampleA", FastqLayoutMode::AUTO, empty));
}

// ---- Filenames ------------------------------------------------------------

TEST_CASE("OutputFilenames: single layout produces one .fastq.gz", "[ena_upload_reads]") {
	auto names = OutputFilenames("sampleA", FastqLayoutMode::SINGLE);
	REQUIRE(names == std::vector<std::string> {"sampleA.fastq.gz"});
}

TEST_CASE("OutputFilenames: paired layout produces _1 and _2 .fastq.gz", "[ena_upload_reads]") {
	auto names = OutputFilenames("sampleA", FastqLayoutMode::PAIRED);
	REQUIRE(names == std::vector<std::string> {"sampleA_1.fastq.gz", "sampleA_2.fastq.gz"});
}

TEST_CASE("OutputFilenames: paired_interleaved produces one combined file", "[ena_upload_reads]") {
	auto names = OutputFilenames("sampleA", FastqLayoutMode::PAIRED_INTERLEAVED);
	REQUIRE(names == std::vector<std::string> {"sampleA.fastq.gz"});
}

TEST_CASE("OutputFilenames: AUTO is rejected — caller must resolve first", "[ena_upload_reads]") {
	REQUIRE_THROWS(OutputFilenames("sampleA", FastqLayoutMode::AUTO));
}

// ---- Layout name parsing --------------------------------------------------

TEST_CASE("ParseFastqLayoutMode: case-insensitive accepted values", "[ena_upload_reads]") {
	REQUIRE(ParseFastqLayoutMode("auto") == FastqLayoutMode::AUTO);
	REQUIRE(ParseFastqLayoutMode("AUTO") == FastqLayoutMode::AUTO);
	REQUIRE(ParseFastqLayoutMode("Single") == FastqLayoutMode::SINGLE);
	REQUIRE(ParseFastqLayoutMode("paired") == FastqLayoutMode::PAIRED);
	REQUIRE(ParseFastqLayoutMode("Paired_Interleaved") == FastqLayoutMode::PAIRED_INTERLEAVED);
}

TEST_CASE("ParseFastqLayoutMode: unknown value names the offender", "[ena_upload_reads]") {
	REQUIRE_THROWS_WITH(ParseFastqLayoutMode("triple"), Catch::Matchers::ContainsSubstring("triple"));
}

TEST_CASE("FastqLayoutModeName: round-trips with ParseFastqLayoutMode", "[ena_upload_reads]") {
	for (auto mode : {FastqLayoutMode::AUTO, FastqLayoutMode::SINGLE, FastqLayoutMode::PAIRED,
	                  FastqLayoutMode::PAIRED_INTERLEAVED}) {
		REQUIRE(ParseFastqLayoutMode(FastqLayoutModeName(mode)) == mode);
	}
}

// ---- Upload target URL parsing -------------------------------------------

TEST_CASE("ParseUploadTargetURL: aspera:// with default root", "[ena_upload_reads]") {
	auto target = ParseUploadTargetURL("aspera://webin2.ebi.ac.uk/");
	REQUIRE(target.transport == UploadTransport::ASPERA);
	REQUIRE(target.host == "webin2.ebi.ac.uk");
	REQUIRE(target.remote_dir == "/");
}

TEST_CASE("ParseUploadTargetURL: aspera:// with subdir; trailing slash normalised", "[ena_upload_reads]") {
	auto target = ParseUploadTargetURL("aspera://webin2.ebi.ac.uk/run42/");
	REQUIRE(target.transport == UploadTransport::ASPERA);
	REQUIRE(target.host == "webin2.ebi.ac.uk");
	REQUIRE(target.remote_dir == "/run42/");
}

TEST_CASE("ParseUploadTargetURL: aspera:// without path component defaults to /", "[ena_upload_reads]") {
	auto target = ParseUploadTargetURL("aspera://webin2.ebi.ac.uk");
	REQUIRE(target.transport == UploadTransport::ASPERA);
	REQUIRE(target.host == "webin2.ebi.ac.uk");
	REQUIRE(target.remote_dir == "/");
}

TEST_CASE("ParseUploadTargetURL: aspera:// missing slash before subdir gets one added", "[ena_upload_reads]") {
	auto target = ParseUploadTargetURL("aspera://webin2.ebi.ac.uk/run42");
	REQUIRE(target.remote_dir == "/run42/");
}

TEST_CASE("ParseUploadTargetURL: file:/// with absolute local path", "[ena_upload_reads]") {
	auto target = ParseUploadTargetURL("file:///tmp/ena-uploads/");
	REQUIRE(target.transport == UploadTransport::LOCAL_FILE);
	REQUIRE(target.host.empty());
	REQUIRE(target.remote_dir == "/tmp/ena-uploads/");
}

TEST_CASE("ParseUploadTargetURL: file:// with relative path is accepted", "[ena_upload_reads]") {
	// DuckDB's per-test __TEST_DIR__ macro substitutes a relative path; the
	// SQL test prefixes it with `file://` so this path-shape needs to work.
	auto target = ParseUploadTargetURL("file://duckdb_unittest_tempdir/X/uploads/");
	REQUIRE(target.transport == UploadTransport::LOCAL_FILE);
	REQUIRE(target.host.empty());
	REQUIRE(target.remote_dir == "duckdb_unittest_tempdir/X/uploads/");
}

TEST_CASE("ParseUploadTargetURL: ftp/ftps/http/https route to CURL transport", "[ena_upload_reads]") {
	for (const std::string scheme : {"ftp", "ftps", "http", "https"}) {
		const auto target = ParseUploadTargetURL(scheme + "://webin2.ebi.ac.uk/uploads/");
		REQUIRE(target.transport == UploadTransport::CURL);
		REQUIRE(target.scheme == scheme);
		REQUIRE(target.host == "webin2.ebi.ac.uk");
		REQUIRE(target.remote_dir == "/uploads/");
		REQUIRE(target.url_for_curl == scheme + "://webin2.ebi.ac.uk/uploads/");
	}
}

TEST_CASE("ParseUploadTargetURL: CURL scheme with no trailing slash gets normalised", "[ena_upload_reads]") {
	const auto target = ParseUploadTargetURL("https://example.org/run42");
	REQUIRE(target.transport == UploadTransport::CURL);
	REQUIRE(target.remote_dir == "/run42/");
	REQUIRE(target.url_for_curl == "https://example.org/run42/");
}

TEST_CASE("ParseUploadTargetURL: unknown scheme is rejected", "[ena_upload_reads]") {
	REQUIRE_THROWS_WITH(ParseUploadTargetURL("smb://server/share/"), Catch::Matchers::ContainsSubstring("smb"));
	REQUIRE_THROWS(ParseUploadTargetURL(""));
	REQUIRE_THROWS(ParseUploadTargetURL("aspera://"));
}

// ---- ascp send argv -------------------------------------------------------

namespace {

bool ContainsToken(const std::vector<std::string> &argv, const std::string &needle) {
	return std::any_of(argv.begin(), argv.end(), [&](const std::string &arg) { return arg == needle; });
}

// Returns the argv element immediately following `flag`, or "" if missing.
std::string ValueAfter(const std::vector<std::string> &argv, const std::string &flag) {
	for (size_t i = 0; i + 1 < argv.size(); i++) {
		if (argv[i] == flag) {
			return argv[i + 1];
		}
	}
	return "";
}

} // namespace

TEST_CASE("BuildAscpSendArgv: full flag set with no rate cap", "[ena_upload_reads]") {
	AsperaSendOptions opts;
	opts.ascp_path = "/usr/bin/ascp";
	opts.key_path = "/home/u/.aspera/connect/etc/asperaweb_id_dsa.openssh";
	opts.user = "Webin-12345";
	opts.host = "webin2.ebi.ac.uk";
	opts.port = 33001;
	opts.local_path = "/tmp/ena-upload-AbCdEf.fastq.gz";
	opts.remote_dir = "/";
	// max_rate left empty.

	auto argv = BuildAscpSendArgv(opts);

	// First argv element is always the binary path (POSIX convention).
	REQUIRE(argv.front() == opts.ascp_path);

	REQUIRE(ContainsToken(argv, "--mode=send"));
	REQUIRE(ContainsToken(argv, "--user=Webin-12345"));
	REQUIRE(ContainsToken(argv, "--host=webin2.ebi.ac.uk"));
	REQUIRE(ContainsToken(argv, "-Q"));
	REQUIRE(ContainsToken(argv, "-d"));
	REQUIRE(ValueAfter(argv, "-i") == opts.key_path);
	REQUIRE(ValueAfter(argv, "-P") == "33001");

	// Source then destination at the tail.
	REQUIRE(argv[argv.size() - 2] == opts.local_path);
	REQUIRE(argv.back() == "/");

	// No rate flag when not requested.
	REQUIRE_FALSE(ContainsToken(argv, "-l"));

	// The password must NEVER appear anywhere in argv — Aspera reads it from
	// ASPERA_SCP_PASS in the child env. Check a sentinel.
	for (const auto &a : argv) {
		REQUIRE(a.find("password") == std::string::npos);
	}
}

TEST_CASE("BuildAscpSendArgv: max_rate adds -l with the supplied value", "[ena_upload_reads]") {
	AsperaSendOptions opts;
	opts.ascp_path = "/usr/bin/ascp";
	opts.key_path = "/k";
	opts.user = "Webin-12345";
	opts.host = "webin2.ebi.ac.uk";
	opts.local_path = "/tmp/x.fastq.gz";
	opts.remote_dir = "/run42/";
	opts.max_rate = "100m";

	auto argv = BuildAscpSendArgv(opts);
	REQUIRE(ValueAfter(argv, "-l") == "100m");
	REQUIRE(argv.back() == "/run42/");
}

TEST_CASE("BuildAscpSendArgv: rejects empty mandatory fields", "[ena_upload_reads]") {
	AsperaSendOptions base;
	base.ascp_path = "/usr/bin/ascp";
	base.key_path = "/k";
	base.user = "Webin-12345";
	base.host = "webin2.ebi.ac.uk";
	base.local_path = "/tmp/x.fastq.gz";
	base.remote_dir = "/";

	auto missing_user = base;
	missing_user.user = "";
	REQUIRE_THROWS(BuildAscpSendArgv(missing_user));

	auto missing_host = base;
	missing_host.host = "";
	REQUIRE_THROWS(BuildAscpSendArgv(missing_host));

	auto missing_ascp = base;
	missing_ascp.ascp_path = "";
	REQUIRE_THROWS(BuildAscpSendArgv(missing_ascp));

	auto missing_key = base;
	missing_key.key_path = "";
	REQUIRE_THROWS(BuildAscpSendArgv(missing_key));

	auto missing_local = base;
	missing_local.local_path = "";
	REQUIRE_THROWS(BuildAscpSendArgv(missing_local));

	auto missing_remote = base;
	missing_remote.remote_dir = "";
	REQUIRE_THROWS(BuildAscpSendArgv(missing_remote));
}
