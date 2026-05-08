// Unit tests for the pure-data layer of INSERT INTO ena.samples: assemble
// envelope, POST via injected functor, parse receipt, return row-shaped
// result tuples (alias, ers_accession, samea_accession, status,
// hold_until_date). Mirrors test_ena_projects_insert.cpp's mock-fetcher
// pattern.

#include "ena_envelope_builder.hpp"
#include "ena_insert_test_helpers.hpp"
#include "ena_receipt_parser.hpp"
#include "ena_samples_insert.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <stdexcept>
#include <string>
#include <vector>

using namespace miint;
using miint_test::CapturedPost;

namespace {

std::string MakeSampleReceipt(const std::vector<std::tuple<std::string, std::string, std::string>> &alias_ers_samea,
                              bool success = true, const std::string &error = "") {
	std::vector<miint_test::ReceiptObjectFixture> objects;
	objects.reserve(alias_ers_samea.size());
	for (const auto &t : alias_ers_samea) {
		objects.push_back({"SAMPLE", std::get<0>(t), std::get<1>(t), "biosample", std::get<2>(t)});
	}
	return miint_test::MakeReceiptXML(objects, success, error);
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
	// Sample envelopes flipped to XML in L4b-fix to bypass the V2 JSON
	// dispatcher's <DESCRIPTION>-before-<SAMPLE_NAME> ordering bug.
	REQUIRE(captured.content_type.find("xml") != std::string::npos);
	// Envelope contains the sample alias, the TAXON_ID, and the ENA-CHECKLIST
	// attribute prepended to the user's attribute list.
	REQUIRE(captured.body.find(R"X(<SAMPLE alias="s1">)X") != std::string::npos);
	REQUIRE(captured.body.find("<TAXON_ID>408170</TAXON_ID>") != std::string::npos);
	REQUIRE(captured.body.find("<TAG>ENA-CHECKLIST</TAG><VALUE>ERC000015</VALUE>") != std::string::npos);
	REQUIRE(captured.body.find("<TAG>collection date</TAG>") != std::string::npos);
}

TEST_CASE("ENA samples insert: multi-sample envelope and order preservation", "[ena_samples_insert]") {
	auto post_fn = [](const std::string &, const std::string &body, const std::string &, const std::string &,
	                  const std::string &) {
		REQUIRE(body.find(R"X(<SAMPLE alias="a">)X") != std::string::npos);
		REQUIRE(body.find(R"X(<SAMPLE alias="b">)X") != std::string::npos);
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

TEST_CASE("ENA samples insert: MODIFY round-trips the user-supplied accession back into the row",
          "[ena_samples_insert][modify]") {
	CapturedPost captured;
	auto post_fn = [&captured](const std::string &url, const std::string &body, const std::string &user,
	                           const std::string &password, const std::string &content_type) {
		captured = {url, body, user, password, content_type};
		// Real wwwdev MODIFY echoes the user-supplied ERS verbatim; the EXT_ID
		// SAMEA accession remains stable across MODIFY (BioSample is the
		// permanent identifier; ENA's ERS is a re-versioning of the same
		// biosample). Mirror that here for fidelity.
		return MakeSampleReceipt({{"s1", "ERS9999100", "SAMEA9999100"}});
	};

	auto sample = MinimalSample("s1");
	sample.accession = "ERS9999100";
	sample.title = "Updated sample title";
	sample.checklist = "ERC000015";
	sample.attributes = {{"collection date", "2026-05-07"}};
	std::vector<SampleSpec> samples = {sample};

	ENASampleInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.action = ENAAction::MODIFY;

	auto outcome = SubmitSampleInsertOutcome(samples, opts, post_fn);
	REQUIRE(outcome.success);
	REQUIRE(outcome.rows.size() == 1);
	REQUIRE(outcome.rows[0].alias == "s1");
	REQUIRE(outcome.rows[0].ers_accession == "ERS9999100");
	REQUIRE(outcome.rows[0].samea_accession == "SAMEA9999100");
	REQUIRE(captured.body.find("<MODIFY/>") != std::string::npos);
	REQUIRE(captured.body.find(R"X(<SAMPLE alias="s1" accession="ERS9999100">)X") != std::string::npos);
}

TEST_CASE("ENA samples insert: MODIFY failure receipt surfaces the server error", "[ena_samples_insert][modify]") {
	auto post_fn = [](const std::string &, const std::string &, const std::string &, const std::string &,
	                  const std::string &) {
		return MakeSampleReceipt({}, false, "ERS00000 not found in submission account");
	};
	auto sample = MinimalSample("s1");
	sample.accession = "ERS00000";
	std::vector<SampleSpec> samples = {sample};

	ENASampleInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.action = ENAAction::MODIFY;

	auto outcome = SubmitSampleInsertOutcome(samples, opts, post_fn);
	REQUIRE_FALSE(outcome.success);
	REQUIRE(outcome.rows.empty());
	REQUIRE(outcome.error_messages.size() == 1);
	REQUIRE_THAT(outcome.error_messages[0], Catch::Matchers::ContainsSubstring("ERS00000 not found"));
}

TEST_CASE("ENA samples insert: hold_until_date round-trips via the submission action", "[ena_samples_insert]") {
	auto post_fn = [](const std::string &, const std::string &body, const std::string &, const std::string &,
	                  const std::string &) {
		REQUIRE(body.find(R"X(<HOLD HoldUntilDate="2027-01-15"/>)X") != std::string::npos);
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
		// SCIENTIFIC_NAME nests inside SAMPLE_NAME; DESCRIPTION is a sibling
		// AFTER SAMPLE_NAME (SRA.sample.xsd ordering — the bug L4b-fix exists
		// to enforce). Pin the relative position so a future reorder of
		// AppendXmlSample regresses loudly.
		REQUIRE(body.find("<SAMPLE_NAME><TAXON_ID>408170</TAXON_ID>"
		                  "<SCIENTIFIC_NAME>human gut metagenome</SCIENTIFIC_NAME></SAMPLE_NAME>") !=
		        std::string::npos);
		REQUIRE(body.find("<DESCRIPTION>clinical isolate from a healthy donor</DESCRIPTION>") != std::string::npos);
		const auto sn_pos = body.find("<SAMPLE_NAME>");
		const auto desc_pos = body.find("<DESCRIPTION>");
		REQUIRE(sn_pos != std::string::npos);
		REQUIRE(desc_pos != std::string::npos);
		REQUIRE(sn_pos < desc_pos);
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
