// Phase 5 RED then GREEN: unit tests for the pure-data layer of
// INSERT INTO ena.samples — assemble envelope, POST via injected functor,
// parse receipt, return row-shaped result tuples (alias, ers_accession,
// samea_accession, status, hold_until_date).
//
// Mirrors test_ena_projects_insert.cpp's mock-fetcher pattern.

#include "ena_envelope_builder.hpp"
#include "ena_receipt_parser.hpp"
#include "ena_samples_insert.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <stdexcept>
#include <string>
#include <vector>

using namespace miint;

namespace {

struct CapturedPost {
	std::string url;
	std::string body;
	std::string user;
	std::string password;
	std::string content_type;
};

std::string MakeSampleReceipt(const std::vector<std::tuple<std::string, std::string, std::string>> &alias_ers_samea,
                              bool success = true, const std::string &error = "") {
	std::string out = "<?xml version=\"1.0\"?><RECEIPT receiptDate=\"2026-05-03T12:00:00Z\" "
	                  "submissionFile=\"mock\" success=\"";
	out += (success ? "true" : "false");
	out += "\">";
	for (const auto &t : alias_ers_samea) {
		out += "<SAMPLE accession=\"" + std::get<1>(t) + "\" alias=\"" + std::get<0>(t) + "\" status=\"PRIVATE\">";
		out += "<EXT_ID accession=\"" + std::get<2>(t) + "\" type=\"biosample\"/>";
		out += "</SAMPLE>";
	}
	out += "<SUBMISSION accession=\"ERA999\" alias=\"mock\"/><ACTIONS>ADD</ACTIONS>";
	if (!success && !error.empty()) {
		out += "<MESSAGES><ERROR>" + error + "</ERROR></MESSAGES>";
	}
	out += "</RECEIPT>";
	return out;
}

SampleSpec MinimalSample(const std::string &alias, int64_t taxon_id = 408170) {
	SampleSpec s;
	s.alias = alias;
	s.taxon_id = taxon_id;
	return s;
}

} // namespace

TEST_CASE("ENA samples insert: single row builds envelope and parses receipt", "[ena_samples_insert]") {
	CapturedPost captured;
	auto post_fn = [&captured](const std::string &url, const std::string &body, const std::string &user,
	                           const std::string &password, const std::string &content_type) {
		captured = {url, body, user, password, content_type};
		return MakeSampleReceipt({{"s1", "ERS42", "SAMEA42"}});
	};

	auto sample = MinimalSample("s1");
	sample.title = "Gut microbiota";
	sample.checklist = "ERC000015";
	sample.attributes = {{"collection date", "2026-04-01"},
	                     {"geographic location (country and/or sea)", "United States"},
	                     {"broad-scale environmental context", "human-associated habitat [ENVO:00009003]"}};
	std::vector<SampleSpec> samples = {sample};

	ENASampleInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	auto rows = SubmitSampleInsert(samples, opts, post_fn);

	REQUIRE(rows.size() == 1);
	REQUIRE(rows[0].alias == "s1");
	REQUIRE(rows[0].ers_accession == "ERS42");
	REQUIRE(rows[0].samea_accession == "SAMEA42");
	REQUIRE(rows[0].status == "PRIVATE");

	REQUIRE(captured.url == "http://mock.example/submit");
	REQUIRE(captured.user == "Webin-1");
	REQUIRE(captured.password == "pw");
	REQUIRE(captured.content_type.find("json") != std::string::npos);
	// Envelope contains the sample alias, the taxonId-as-string, and the
	// ENA-CHECKLIST attribute prepended to the user's attribute list.
	REQUIRE(captured.body.find("\"alias\":\"s1\"") != std::string::npos);
	REQUIRE(captured.body.find("\"taxonId\":\"408170\"") != std::string::npos);
	REQUIRE(captured.body.find("\"tag\":\"ENA-CHECKLIST\",\"value\":\"ERC000015\"") != std::string::npos);
	REQUIRE(captured.body.find("\"tag\":\"collection date\"") != std::string::npos);
}

TEST_CASE("ENA samples insert: multi-sample envelope and order preservation", "[ena_samples_insert]") {
	auto post_fn = [](const std::string &, const std::string &body, const std::string &, const std::string &,
	                  const std::string &) {
		REQUIRE(body.find("\"alias\":\"a\"") != std::string::npos);
		REQUIRE(body.find("\"alias\":\"b\"") != std::string::npos);
		return MakeSampleReceipt({{"a", "ERS10", "SAMEA10"}, {"b", "ERS11", "SAMEA11"}});
	};

	std::vector<SampleSpec> samples = {MinimalSample("a"), MinimalSample("b", 9606)};
	ENASampleInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	auto rows = SubmitSampleInsert(samples, opts, post_fn);

	REQUIRE(rows.size() == 2);
	REQUIRE(rows[0].alias == "a");
	REQUIRE(rows[0].ers_accession == "ERS10");
	REQUIRE(rows[0].samea_accession == "SAMEA10");
	REQUIRE(rows[1].alias == "b");
	REQUIRE(rows[1].ers_accession == "ERS11");
	REQUIRE(rows[1].samea_accession == "SAMEA11");
}

TEST_CASE("ENA samples insert: receipt success=false throws with messages", "[ena_samples_insert]") {
	auto post_fn = [](const std::string &, const std::string &, const std::string &, const std::string &,
	                  const std::string &) {
		return MakeSampleReceipt({}, false, "checklist ERC000015 missing mandatory attribute 'collection date'");
	};

	std::vector<SampleSpec> samples = {MinimalSample("foo")};
	ENASampleInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	REQUIRE_THROWS_WITH(SubmitSampleInsert(samples, opts, post_fn),
	                    Catch::Matchers::ContainsSubstring("checklist ERC000015"));
}

TEST_CASE("ENA samples insert: empty input is a no-op (no POST issued)", "[ena_samples_insert]") {
	bool called = false;
	auto post_fn = [&called](const std::string &, const std::string &, const std::string &, const std::string &,
	                         const std::string &) {
		called = true;
		return std::string();
	};

	std::vector<SampleSpec> samples;
	ENASampleInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	auto rows = SubmitSampleInsert(samples, opts, post_fn);
	REQUIRE(rows.empty());
	REQUIRE_FALSE(called);
}

TEST_CASE("ENA samples insert: receipt missing alias is reported clearly", "[ena_samples_insert]") {
	auto post_fn = [](const std::string &, const std::string &, const std::string &, const std::string &,
	                  const std::string &) {
		return MakeSampleReceipt({}, true);
	};

	std::vector<SampleSpec> samples = {MinimalSample("orphan")};
	ENASampleInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	REQUIRE_THROWS_WITH(SubmitSampleInsert(samples, opts, post_fn), Catch::Matchers::ContainsSubstring("orphan"));
}

TEST_CASE("ENA samples insert: hold_until_date round-trips via the submission action", "[ena_samples_insert]") {
	auto post_fn = [](const std::string &, const std::string &body, const std::string &, const std::string &,
	                  const std::string &) {
		REQUIRE(body.find("\"holdUntilDate\":\"2027-01-15\"") != std::string::npos);
		std::string receipt = "<RECEIPT success=\"true\">";
		receipt += "<SAMPLE accession=\"ERS7\" alias=\"h\" status=\"PRIVATE\" holdUntilDate=\"2027-01-15\">";
		receipt += "<EXT_ID accession=\"SAMEA7\" type=\"biosample\"/></SAMPLE>";
		receipt += "</RECEIPT>";
		return receipt;
	};

	std::vector<SampleSpec> samples = {MinimalSample("h")};
	ENASampleInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.hold_until_date = "2027-01-15";

	auto rows = SubmitSampleInsert(samples, opts, post_fn);
	REQUIRE(rows.size() == 1);
	REQUIRE(rows[0].hold_until_date == "2027-01-15");
}

TEST_CASE("ENA samples insert: SampleSpec exposes scientific_name and description in envelope",
          "[ena_samples_insert]") {
	auto post_fn = [](const std::string &, const std::string &body, const std::string &, const std::string &,
	                  const std::string &) {
		// SCIENTIFIC_NAME goes into organism, DESCRIPTION at the sample level.
		REQUIRE(body.find("\"scientificName\":\"human gut metagenome\"") != std::string::npos);
		REQUIRE(body.find("\"description\":\"clinical isolate from a healthy donor\"") != std::string::npos);
		return MakeSampleReceipt({{"sn", "ERS9", "SAMEA9"}});
	};

	SampleSpec sample = MinimalSample("sn");
	sample.scientific_name = "human gut metagenome";
	sample.description = "clinical isolate from a healthy donor";
	std::vector<SampleSpec> samples = {sample};

	ENASampleInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	auto rows = SubmitSampleInsert(samples, opts, post_fn);
	REQUIRE(rows.size() == 1);
	REQUIRE(rows[0].alias == "sn");
}
