#include <catch2/catch_test_macros.hpp>
#include "taxdump_archive.hpp"
#include "taxdump_test_helpers.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace miint;
using miint_test::ThrowsMentioning;

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
	auto files = TaxdumpArchive::ExtractTaxdump(SlurpFile("data/taxonomy/synthetic.tar.gz"), TaxdumpMemberSet::All());
	CHECK(files.nodes == SlurpFile("data/taxonomy/synthetic/nodes.dmp"));
	CHECK(files.names == SlurpFile("data/taxonomy/synthetic/names.dmp"));
	CHECK(files.merged == SlurpFile("data/taxonomy/synthetic/merged.dmp"));
	CHECK(files.delnodes == SlurpFile("data/taxonomy/synthetic/delnodes.dmp"));
}

// Each reader requests only what it parses, so extracting one member must not pay for
// the others: on the real dump, answering from delnodes.dmp (~7 MB) used to copy
// nodes.dmp and names.dmp (~500 MB) out of the archive as well.
TEST_CASE("TaxdumpArchive::ExtractTaxdump copies only the requested members", "[taxdump]") {
	TaxdumpMemberSet only_delnodes;
	only_delnodes.delnodes = true;
	auto files = TaxdumpArchive::ExtractTaxdump(SlurpFile("data/taxonomy/synthetic.tar.gz"), only_delnodes);
	CHECK(files.delnodes == SlurpFile("data/taxonomy/synthetic/delnodes.dmp"));
	CHECK(files.nodes.empty());
	CHECK(files.names.empty());
	CHECK(files.merged.empty());
}

// A requested member the archive does not carry is an error, not an empty string:
// "this release retires nothing" and "the extraction was incomplete" must not look the
// same to a consumer building a versioned release. noname.tar.gz carries only
// nodes.dmp and names.dmp, so it exercises the archive branch of that rule.
TEST_CASE("TaxdumpArchive::ExtractTaxdump rejects a missing requested member", "[taxdump]") {
	std::string archive = SlurpFile("data/taxonomy/noname.tar.gz");

	SECTION("requesting an absent member throws, naming it") {
		TaxdumpMemberSet wants_delnodes;
		wants_delnodes.delnodes = true;
		CHECK(ThrowsMentioning([&] { TaxdumpArchive::ExtractTaxdump(archive, wants_delnodes); }, {"delnodes.dmp"}));

		TaxdumpMemberSet wants_merged;
		wants_merged.merged = true;
		CHECK(ThrowsMentioning([&] { TaxdumpArchive::ExtractTaxdump(archive, wants_merged); }, {"merged.dmp"}));
	}

	SECTION("requesting only the members it does carry succeeds") {
		// The same archive is perfectly usable for the tree, which reads neither of the
		// missing members -- so the requirement is per-caller, not archive-wide.
		TaxdumpMemberSet wants_tree;
		wants_tree.nodes = true;
		wants_tree.names = true;
		auto files = TaxdumpArchive::ExtractTaxdump(archive, wants_tree);
		CHECK(files.nodes == SlurpFile("data/taxonomy/noname/nodes.dmp"));
		CHECK(files.names == SlurpFile("data/taxonomy/noname/names.dmp"));
	}
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
