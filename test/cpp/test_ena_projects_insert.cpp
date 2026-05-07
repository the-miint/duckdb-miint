// Unit tests for the pure-data layer of INSERT INTO ena.projects: assemble
// envelope, POST via injected functor, parse receipt, return row-shaped
// result tuples. Mirrors the test_ENAClient mock-fetcher pattern (string ->
// string functor) so we can drive end-to-end submit logic without TCP,
// threads, or DuckDB.

#include "ena_envelope_builder.hpp"
#include "ena_insert_test_helpers.hpp"
#include "ena_projects_insert.hpp"
#include "ena_receipt_parser.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <stdexcept>
#include <string>
#include <vector>

using namespace miint;
using miint_test::CapturedPost;

namespace {

std::string MakeReceipt(const std::vector<std::pair<std::string, std::string>> &alias_to_accession, bool success = true,
                        const std::string &error = "") {
	std::vector<miint_test::ReceiptObjectFixture> objects;
	objects.reserve(alias_to_accession.size());
	for (const auto &kv : alias_to_accession) {
		objects.push_back({"PROJECT", kv.first, kv.second, "study", "ERP" + kv.second.substr(5)});
	}
	return miint_test::MakeReceiptXML(objects, success, error);
}

} // namespace

TEST_CASE("ENA projects insert: single row builds correct envelope and parses receipt", "[ena_projects_insert]") {
	CapturedPost captured;
	auto post_fn = [&captured](const std::string &url, const std::string &body, const std::string &user,
	                           const std::string &password, const std::string &content_type) {
		captured = {url, body, user, password, content_type};
		return MakeReceipt({{"p1", "PRJEB42"}});
	};

	std::vector<ProjectSpec> projects = {{"p1", "Demo project", "An example", "METAGENOMIC", false}};
	ENAProjectInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	auto rows = SubmitProjectInsert(projects, opts, post_fn);

	REQUIRE(rows.size() == 1);
	REQUIRE(rows[0].alias == "p1");
	REQUIRE(rows[0].prjeb_accession == "PRJEB42");
	REQUIRE(rows[0].erp_accession == "ERP42");
	REQUIRE(rows[0].status == "PRIVATE");
	REQUIRE(rows[0].hold_until_date.empty());

	REQUIRE(captured.url == "http://mock.example/submit");
	REQUIRE(captured.user == "Webin-1");
	REQUIRE(captured.password == "pw");
	REQUIRE(captured.content_type.find("json") != std::string::npos);
	// Envelope contains the project alias and the ADD action
	REQUIRE(captured.body.find("\"alias\":\"p1\"") != std::string::npos);
	REQUIRE(captured.body.find("\"type\":\"ADD\"") != std::string::npos);
	REQUIRE(captured.body.find("\"sequencingProject\":{}") != std::string::npos);
}

TEST_CASE("ENA projects insert: multi-row envelope and order preservation", "[ena_projects_insert]") {
	CapturedPost captured;
	auto post_fn = [&captured](const std::string &url, const std::string &body, const std::string &user,
	                           const std::string &password, const std::string &content_type) {
		captured = {url, body, user, password, content_type};
		return MakeReceipt({{"a", "PRJEB10"}, {"b", "PRJEB11"}});
	};

	std::vector<ProjectSpec> projects = {{"a", "A", "", "METAGENOMIC", false}, {"b", "B", "", "WGS", false}};
	ENAProjectInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	auto rows = SubmitProjectInsert(projects, opts, post_fn);

	REQUIRE(rows.size() == 2);
	REQUIRE(rows[0].alias == "a");
	REQUIRE(rows[0].prjeb_accession == "PRJEB10");
	REQUIRE(rows[1].alias == "b");
	REQUIRE(rows[1].prjeb_accession == "PRJEB11");
}

TEST_CASE("ENA projects insert: receipt success=false throws with messages", "[ena_projects_insert]") {
	auto post_fn = [](const std::string &, const std::string &, const std::string &, const std::string &,
	                  const std::string &) {
		return MakeReceipt({}, false, "alias 'foo' already exists");
	};

	std::vector<ProjectSpec> projects = {{"foo", "Foo", "", "METAGENOMIC", false}};
	ENAProjectInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	REQUIRE_THROWS_WITH(SubmitProjectInsert(projects, opts, post_fn),
	                    Catch::Matchers::ContainsSubstring("alias 'foo' already exists"));
}

TEST_CASE("ENA projects insert: empty input is a no-op (no POST issued)", "[ena_projects_insert]") {
	bool called = false;
	auto post_fn = [&called](const std::string &, const std::string &, const std::string &, const std::string &,
	                         const std::string &) {
		called = true;
		return std::string();
	};

	std::vector<ProjectSpec> projects;
	ENAProjectInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	auto rows = SubmitProjectInsert(projects, opts, post_fn);
	REQUIRE(rows.empty());
	REQUIRE_FALSE(called);
}

TEST_CASE("ENA projects insert: receipt missing alias is reported clearly", "[ena_projects_insert]") {
	auto post_fn = [](const std::string &, const std::string &, const std::string &, const std::string &,
	                  const std::string &) {
		// receipt with no objects but success=true — pathological
		return MakeReceipt({}, true);
	};

	std::vector<ProjectSpec> projects = {{"orphan", "O", "", "METAGENOMIC", false}};
	ENAProjectInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	REQUIRE_THROWS_WITH(SubmitProjectInsert(projects, opts, post_fn), Catch::Matchers::ContainsSubstring("orphan"));
}

TEST_CASE("ENA projects insert: ADD failure on missing-alias clears any partially-populated rows",
          "[ena_projects_insert]") {
	// Two specs ('a' and 'b'); receipt only carries 'a'. The mid-loop bail
	// must drop the already-populated row for 'a' so callers reading
	// `outcome.rows` after `outcome.success == false` see a consistent empty
	// vector — not a half-filled prefix.
	auto post_fn = [](const std::string &, const std::string &, const std::string &, const std::string &,
	                  const std::string &) {
		return MakeReceipt({{"a", "PRJEB100"}}, true);
	};

	std::vector<ProjectSpec> projects = {{"a", "A", "", "METAGENOMIC", false}, {"b", "B", "", "WGS", false}};
	ENAProjectInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	auto outcome = SubmitProjectInsertOutcome(projects, opts, post_fn);
	REQUIRE_FALSE(outcome.success);
	REQUIRE(outcome.rows.empty());
	REQUIRE(outcome.error_messages.size() == 1);
	REQUIRE_THAT(outcome.error_messages[0], Catch::Matchers::ContainsSubstring("'b'"));
}

TEST_CASE("ENA projects insert: VALIDATE action emits VALIDATE in envelope and tolerates an empty per-object receipt",
          "[ena_projects_insert]") {
	// VALIDATE is a dry-run: ENA's wwwdev returns no per-object accessions on
	// validation receipts (just a <SUBMISSION/> stamp). The pure-data layer
	// must (a) emit "VALIDATE" in the envelope wire-format, (b) NOT error on
	// the absent per-spec receipt entries, (c) leave `outcome.rows` empty,
	// and (d) still surface era_accession so the submission_log row carries it.
	CapturedPost captured;
	auto post_fn = [&captured](const std::string &url, const std::string &body, const std::string &user,
	                           const std::string &password, const std::string &content_type) {
		captured = {url, body, user, password, content_type};
		return std::string("<?xml version=\"1.0\"?><RECEIPT receiptDate=\"2026-05-07T00:00:00Z\" "
		                   "submissionFile=\"mock\" success=\"true\">"
		                   "<SUBMISSION accession=\"ERA-VAL-1\" alias=\"mock\"/>"
		                   "<ACTIONS>VALIDATE</ACTIONS></RECEIPT>");
	};

	std::vector<ProjectSpec> projects = {{"vp1", "Validate me", "", "METAGENOMIC", false}};
	ENAProjectInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.action = ENAAction::VALIDATE;

	auto outcome = SubmitProjectInsertOutcome(projects, opts, post_fn);

	REQUIRE(outcome.success);
	REQUIRE(outcome.rows.empty());
	REQUIRE(outcome.era_accession == "ERA-VAL-1");
	REQUIRE(outcome.error_messages.empty());
	REQUIRE(captured.body.find("\"type\":\"VALIDATE\"") != std::string::npos);
	REQUIRE(captured.body.find("\"alias\":\"vp1\"") != std::string::npos);
}

TEST_CASE("ENA projects insert: VALIDATE failure path surfaces error_messages and success=false",
          "[ena_projects_insert]") {
	auto post_fn = [](const std::string &, const std::string &, const std::string &, const std::string &,
	                  const std::string &) {
		return std::string("<?xml version=\"1.0\"?><RECEIPT success=\"false\">"
		                   "<SUBMISSION accession=\"ERA-VAL-2\" alias=\"mock\"/>"
		                   "<ACTIONS>VALIDATE</ACTIONS>"
		                   "<MESSAGES><ERROR>missing required field 'taxon_id'</ERROR></MESSAGES>"
		                   "</RECEIPT>");
	};

	std::vector<ProjectSpec> projects = {{"badalias", "Bad", "", "METAGENOMIC", false}};
	ENAProjectInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.action = ENAAction::VALIDATE;

	auto outcome = SubmitProjectInsertOutcome(projects, opts, post_fn);
	REQUIRE_FALSE(outcome.success);
	REQUIRE(outcome.rows.empty());
	REQUIRE(outcome.error_messages.size() == 1);
	REQUIRE(outcome.error_messages[0] == "missing required field 'taxon_id'");
}

TEST_CASE("ENA projects insert: MODIFY round-trips the user-supplied accession back into the row",
          "[ena_projects_insert][modify]") {
	// MODIFY: user provides the existing PRJEB on the spec, the envelope
	// carries it on the project element, the receipt echoes it back.
	// `outcome.rows[0].prjeb_accession` must reflect the user's input (no
	// fabrication, no derivation from alias as the mock does for ADD).
	CapturedPost captured;
	auto post_fn = [&captured](const std::string &url, const std::string &body, const std::string &user,
	                           const std::string &password, const std::string &content_type) {
		captured = {url, body, user, password, content_type};
		return MakeReceipt({{"p1", "PRJEB55555"}});
	};

	std::vector<ProjectSpec> projects = {{"p1", "Updated title", "Updated abstract", "METAGENOMIC", false}};
	projects[0].accession = "PRJEB55555";
	ENAProjectInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.action = ENAAction::MODIFY;

	auto outcome = SubmitProjectInsertOutcome(projects, opts, post_fn);
	REQUIRE(outcome.success);
	REQUIRE(outcome.rows.size() == 1);
	REQUIRE(outcome.rows[0].alias == "p1");
	REQUIRE(outcome.rows[0].prjeb_accession == "PRJEB55555");
	REQUIRE(captured.body.find("\"type\":\"MODIFY\"") != std::string::npos);
	REQUIRE(captured.body.find("\"accession\":\"PRJEB55555\"") != std::string::npos);
}

TEST_CASE("ENA projects insert: MODIFY failure receipt surfaces the server error", "[ena_projects_insert][modify]") {
	auto post_fn = [](const std::string &, const std::string &, const std::string &, const std::string &,
	                  const std::string &) {
		return MakeReceipt({}, false, "PRJEB00000 not found in submission account");
	};
	std::vector<ProjectSpec> projects = {{"p1", "T", "D", "METAGENOMIC", false}};
	projects[0].accession = "PRJEB00000";
	ENAProjectInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.action = ENAAction::MODIFY;

	auto outcome = SubmitProjectInsertOutcome(projects, opts, post_fn);
	REQUIRE_FALSE(outcome.success);
	REQUIRE(outcome.rows.empty());
	REQUIRE(outcome.error_messages.size() == 1);
	REQUIRE_THAT(outcome.error_messages[0], Catch::Matchers::ContainsSubstring("PRJEB00000 not found"));
}

TEST_CASE("ENA projects insert: hold_until_date round-trips into the row", "[ena_projects_insert]") {
	auto post_fn = [](const std::string &, const std::string &body, const std::string &, const std::string &,
	                  const std::string &) {
		// Confirm the envelope carried the HOLD action
		REQUIRE(body.find("\"holdUntilDate\":\"2027-01-15\"") != std::string::npos);
		std::string receipt = "<RECEIPT success=\"true\">";
		receipt += "<PROJECT accession=\"PRJEB7\" alias=\"h\" status=\"PRIVATE\" "
		           "holdUntilDate=\"2027-01-15\">";
		receipt += "<EXT_ID accession=\"ERP7\" type=\"study\"/></PROJECT>";
		receipt += "</RECEIPT>";
		return receipt;
	};

	std::vector<ProjectSpec> projects = {{"h", "Held project", "", "METAGENOMIC", false}};
	ENAProjectInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.hold_until_date = "2027-01-15";

	auto rows = SubmitProjectInsert(projects, opts, post_fn);
	REQUIRE(rows.size() == 1);
	REQUIRE(rows[0].hold_until_date == "2027-01-15");
}
