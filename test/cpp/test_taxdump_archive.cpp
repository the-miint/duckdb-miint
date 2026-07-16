#include <catch2/catch_test_macros.hpp>
#include "taxdump_archive.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

using namespace miint;

namespace {
// C++ tests are run from the repo root (see run_tests.sh), so data/ paths resolve.
std::string SlurpFile(const std::string &path) {
	std::ifstream in(path, std::ios::binary);
	REQUIRE_FALSE(in.fail()); // run the tests binary from the repo root
	std::ostringstream ss;
	ss << in.rdbuf();
	return ss.str();
}
} // namespace

// End-to-end: real gzip + real (Python-generated ustar) tar. A byte-exact match
// against the committed .dmp files proves both the gunzip and the microtar-driven
// in-memory extraction are correct.
TEST_CASE("TaxdumpArchive::ExtractTaxdump pulls byte-exact members from a .tar.gz", "[taxdump]") {
	auto files = TaxdumpArchive::ExtractTaxdump(SlurpFile("data/taxonomy/synthetic.tar.gz"));
	CHECK(files.nodes == SlurpFile("data/taxonomy/synthetic/nodes.dmp"));
	CHECK(files.names == SlurpFile("data/taxonomy/synthetic/names.dmp"));
	CHECK(files.merged == SlurpFile("data/taxonomy/synthetic/merged.dmp"));
	CHECK(files.delnodes == SlurpFile("data/taxonomy/synthetic/delnodes.dmp"));
}

// A half-downloaded taxdump must never be silently accepted as partial data
// (CLAUDE.md Rule 10 "Fail loud"). Truncating the valid fixture mid-stream means
// inflate can never reach the gzip trailer, so Gunzip must throw rather than
// return a truncated buffer that would later drop optional members undetected.
TEST_CASE("TaxdumpArchive::Gunzip rejects a truncated gzip stream", "[taxdump]") {
	std::string full = SlurpFile("data/taxonomy/synthetic.tar.gz");
	REQUIRE(full.size() > 20);
	std::string truncated = full.substr(0, full.size() / 2);
	CHECK_THROWS_AS(TaxdumpArchive::Gunzip(truncated), std::runtime_error);
}

TEST_CASE("TaxdumpArchive::ExtractTarMembers returns only present members", "[taxdump]") {
	std::string tar_bytes = TaxdumpArchive::Gunzip(SlurpFile("data/taxonomy/synthetic.tar.gz"));
	auto members = TaxdumpArchive::ExtractTarMembers(tar_bytes, {"nodes.dmp", "does-not-exist.dmp"});
	CHECK(members.count("nodes.dmp") == 1);
	CHECK(members.count("does-not-exist.dmp") == 0); // a requested-but-absent member is silently skipped
	CHECK(members["nodes.dmp"] == SlurpFile("data/taxonomy/synthetic/nodes.dmp"));
}
