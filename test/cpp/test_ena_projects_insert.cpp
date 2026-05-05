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
