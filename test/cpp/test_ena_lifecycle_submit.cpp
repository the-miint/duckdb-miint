// SPDX-License-Identifier: MIT
//
// Pure-data tests for SubmitLifecycle: the shared CANCEL / RELEASE /
// targeted-HOLD core that L1's ena_cancel/ena_release/ena_hold table
// functions and PlanDelete dispatch through. Mirrors the
// test_ena_projects_insert.cpp pattern: inject a fake ENAPostFn,
// inspect the captured request, return a canned receipt, assert the
// outcome.

#include "ena_envelope_builder.hpp"
#include "ena_insert_test_helpers.hpp"
#include "ena_lifecycle_submit.hpp"

#include <catch2/catch_all.hpp>

#include <cstring>

using miint_test::CapturedPost;
using miint_test::StubPost;

namespace {

// Build a targeted-action success receipt that mirrors what wwwdev actually
// returns (verified live 2026-05-07 against PRJEB112641): RECEIPT
// @success="false" because no new accessions are assigned, with an <INFO>
// child confirming the action took effect. SubmitLifecycle classifies this
// as success because there are no <ERROR> elements. No <SUBMISSION> child
// — lifecycle ops on existing objects don't carry a submission accession.
std::string LifecycleSuccessReceipt(const char *action) {
	const char *status_word = std::strcmp(action, "CANCEL") == 0    ? "cancelled"
	                          : std::strcmp(action, "RELEASE") == 0 ? "public"
	                          : std::strcmp(action, "HOLD") == 0    ? "private"
	                                                                : "updated";
	std::string out = R"(<?xml version="1.0" encoding="UTF-8"?>
<RECEIPT receiptDate="2026-05-06T12:00:00.000Z" submissionFile="mock.xml" success="false">
    <MESSAGES><INFO>accession is set to )";
	out += status_word;
	out += " status.</INFO></MESSAGES>\n    <ACTIONS>";
	out += action;
	out += "</ACTIONS></RECEIPT>";
	return out;
}

} // namespace

TEST_CASE("SubmitLifecycle: CANCEL success emits target= and parses ERA accession", "[ena_lifecycle_submit][cancel]") {
	miint::LifecycleSubmitOptions opts;
	opts.endpoint_url = "https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.target.accession = "ERS123456";

	CapturedPost captured;
	auto outcome =
	    miint::SubmitLifecycle(miint::ENAAction::CANCEL, opts, StubPost(captured, LifecycleSuccessReceipt("CANCEL")));

	CHECK(outcome.success == true);
	CHECK(outcome.action == miint::ENAAction::CANCEL);
	CHECK(outcome.target == "ERS123456");
	// Lifecycle receipts on existing objects don't carry a submission
	// accession — the receipt has no <SUBMISSION> child.
	CHECK(outcome.era_accession.empty());
	CHECK(outcome.error_messages.empty());
	// Wire-form sanity: XML payload includes the right <CANCEL target=…/>
	CHECK(captured.content_type == "application/xml");
	CHECK(captured.url == opts.endpoint_url);
	CHECK(captured.user == "Webin-1");
	CHECK(captured.body.find(R"(<CANCEL target="ERS123456"/>)") != std::string::npos);
	// Defensive: no PROJECT_SET/SAMPLE_SET/etc. in the body
	CHECK(captured.body.find("PROJECT_SET") == std::string::npos);
	CHECK(captured.body.find("SAMPLE_SET") == std::string::npos);
}

TEST_CASE("SubmitLifecycle: CANCEL of nonexistent target reports failure with errors",
          "[ena_lifecycle_submit][cancel]") {
	miint::LifecycleSubmitOptions opts;
	opts.endpoint_url = "https://test/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.target.accession = "ERS-doesnotexist";

	const std::string response = R"(<?xml version="1.0" encoding="UTF-8"?>
<RECEIPT receiptDate="2026-05-06T12:00:00.000Z" submissionFile="mock.xml" success="false">
    <MESSAGES><ERROR>Object 'ERS-doesnotexist' not found in submission account</ERROR></MESSAGES>
    <ACTIONS>CANCEL</ACTIONS>
</RECEIPT>)";

	CapturedPost captured;
	auto outcome = miint::SubmitLifecycle(miint::ENAAction::CANCEL, opts, StubPost(captured, response));

	CHECK(outcome.success == false);
	REQUIRE(outcome.error_messages.size() == 1);
	CHECK(outcome.error_messages[0].find("not found") != std::string::npos);
}

TEST_CASE("SubmitLifecycle: RELEASE emits target=", "[ena_lifecycle_submit][release]") {
	miint::LifecycleSubmitOptions opts;
	opts.endpoint_url = "https://test/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.target.accession = "ERS123456";

	CapturedPost captured;
	auto outcome =
	    miint::SubmitLifecycle(miint::ENAAction::RELEASE, opts, StubPost(captured, LifecycleSuccessReceipt("RELEASE")));

	CHECK(outcome.success == true);
	CHECK(outcome.action == miint::ENAAction::RELEASE);
	CHECK(outcome.target == "ERS123456");
	CHECK(captured.body.find(R"(<RELEASE target="ERS123456"/>)") != std::string::npos);
}

TEST_CASE("SubmitLifecycle: HOLD emits target= and HoldUntilDate", "[ena_lifecycle_submit][hold]") {
	miint::LifecycleSubmitOptions opts;
	opts.endpoint_url = "https://test/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.target.accession = "ERS123456";
	opts.hold_until_date = "2027-12-31";

	const std::string response = R"(<?xml version="1.0" encoding="UTF-8"?>
<RECEIPT receiptDate="2026-05-06T12:00:00.000Z" submissionFile="mock.xml" success="true">
    <SAMPLE accession="ERS123456" status="PRIVATE" holdUntilDate="2027-12-31Z"/>
    <SUBMISSION accession="ERA77777"/>
    <ACTIONS>HOLD</ACTIONS>
</RECEIPT>)";

	CapturedPost captured;
	auto outcome = miint::SubmitLifecycle(miint::ENAAction::HOLD, opts, StubPost(captured, response));

	CHECK(outcome.success == true);
	CHECK(outcome.action == miint::ENAAction::HOLD);
	CHECK(outcome.target == "ERS123456");
	CHECK(outcome.hold_until_date == "2027-12-31");
	CHECK(outcome.era_accession == "ERA77777");
	CHECK(captured.body.find(R"(target="ERS123456")") != std::string::npos);
	CHECK(captured.body.find(R"(HoldUntilDate="2027-12-31")") != std::string::npos);
}

TEST_CASE("SubmitLifecycle: refname is used when accession is empty", "[ena_lifecycle_submit][cancel]") {
	miint::LifecycleSubmitOptions opts;
	opts.endpoint_url = "https://test/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.target.refname = "my-sample-alias";

	CapturedPost captured;
	auto outcome =
	    miint::SubmitLifecycle(miint::ENAAction::CANCEL, opts, StubPost(captured, LifecycleSuccessReceipt("CANCEL")));

	CHECK(outcome.success == true);
	CHECK(outcome.target == "my-sample-alias");
	CHECK(captured.body.find(R"(<CANCEL target="my-sample-alias"/>)") != std::string::npos);
}

TEST_CASE("SubmitLifecycle: wwwdev returns success=\"false\" with <INFO> on a successful CANCEL",
          "[ena_lifecycle_submit][cancel]") {
	// Real wwwdev shape verified live on 2026-05-07 against PRJEB112641:
	// the RECEIPT @success attribute is always "false" for cross-submission
	// lifecycle (no new accessions assigned). The actual outcome is in
	// <INFO> children. SubmitLifecycle treats the action as successful iff
	// the receipt has no <ERROR> elements.
	miint::LifecycleSubmitOptions opts;
	opts.endpoint_url = "https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.target.accession = "PRJEB112641";

	const std::string response = R"(<?xml version="1.0" encoding="UTF-8"?>
<RECEIPT receiptDate="2026-05-07T16:42:55.855+01:00" success="false">
    <MESSAGES>
        <INFO>STUDY accession "ERP193165" is set to cancelled status.</INFO>
        <INFO>PROJECT accession "PRJEB112641" is set to cancelled status.</INFO>
    </MESSAGES>
    <ACTIONS>CANCEL</ACTIONS>
</RECEIPT>)";

	CapturedPost captured;
	auto outcome = miint::SubmitLifecycle(miint::ENAAction::CANCEL, opts, StubPost(captured, response));

	CHECK(outcome.success == true);
	CHECK(outcome.error_messages.empty());
	CHECK(outcome.target == "PRJEB112641");
	// Lifecycle receipts on existing objects carry no <SUBMISSION> child;
	// pin this so a future change that started forwarding era_accession
	// from somewhere else would surface here.
	CHECK(outcome.era_accession.empty());
}

TEST_CASE("SubmitLifecycle: malformed receipt yields outcome.success=false with parse error",
          "[ena_lifecycle_submit]") {
	miint::LifecycleSubmitOptions opts;
	opts.endpoint_url = "https://test/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.target.accession = "ERS123456";

	CapturedPost captured;
	auto outcome = miint::SubmitLifecycle(miint::ENAAction::CANCEL, opts, StubPost(captured, "<unclosed-tag"));

	CHECK(outcome.success == false);
	REQUIRE(outcome.error_messages.size() == 1);
	CHECK(outcome.error_messages[0].find("parse") != std::string::npos);
	// Raw receipt is preserved even when parsing fails — needed for audit log.
	CHECK(outcome.raw_receipt == "<unclosed-tag");
}

TEST_CASE("SubmitLifecycle: missing target throws before POST", "[ena_lifecycle_submit][cancel]") {
	miint::LifecycleSubmitOptions opts;
	opts.endpoint_url = "https://test/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	// no target

	bool post_called = false;
	auto post_fn = [&post_called](const std::string &, const std::string &, const std::string &, const std::string &,
	                              const std::string &) -> std::string {
		post_called = true;
		return "";
	};
	CHECK_THROWS_WITH(miint::SubmitLifecycle(miint::ENAAction::CANCEL, opts, post_fn),
	                  Catch::Matchers::ContainsSubstring("CANCEL"));
	CHECK_FALSE(post_called);
}

TEST_CASE("SubmitLifecycle: ADD/MODIFY/VALIDATE rejected (wrong API for those actions)", "[ena_lifecycle_submit]") {
	miint::LifecycleSubmitOptions opts;
	opts.endpoint_url = "https://test/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.target.accession = "ERS123456";

	CapturedPost captured;
	for (auto bad_action : {miint::ENAAction::ADD, miint::ENAAction::MODIFY, miint::ENAAction::VALIDATE}) {
		CHECK_THROWS_WITH(miint::SubmitLifecycle(bad_action, opts, StubPost(captured, "")),
		                  Catch::Matchers::ContainsSubstring("SubmitLifecycle"));
	}
}

TEST_CASE("SubmitLifecycle: duration_ms is sub-second under stubbed POST", "[ena_lifecycle_submit]") {
	miint::LifecycleSubmitOptions opts;
	opts.endpoint_url = "https://test/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.target.accession = "ERS123456";

	CapturedPost captured;
	auto outcome =
	    miint::SubmitLifecycle(miint::ENAAction::CANCEL, opts, StubPost(captured, LifecycleSuccessReceipt("CANCEL")));
	// Catches an accidental sleep-unit-of-measurement bug — a stubbed POST
	// should never take a full second. (>= 0 alone is vacuous on a signed
	// type with a monotonic clock.)
	CHECK(outcome.duration_ms < 1000);
}

TEST_CASE("SubmitLifecycle: hold_until_date on a non-HOLD action is rejected", "[ena_lifecycle_submit]") {
	// hold_until_date is only meaningful with ADD or HOLD. Setting it on
	// CANCEL/RELEASE is silently ignored by the wire format, which would
	// otherwise let the outcome.hold_until_date echo back into submission_log
	// and mislead audit consumers into thinking the action set a date.
	miint::LifecycleSubmitOptions opts;
	opts.endpoint_url = "https://test/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.target.accession = "ERS123456";
	opts.hold_until_date = "2027-12-31"; // meaningless for CANCEL

	bool post_called = false;
	auto post_fn = [&post_called](const std::string &, const std::string &, const std::string &, const std::string &,
	                              const std::string &) -> std::string {
		post_called = true;
		return "";
	};
	CHECK_THROWS_WITH(miint::SubmitLifecycle(miint::ENAAction::CANCEL, opts, post_fn),
	                  Catch::Matchers::ContainsSubstring("hold_until_date"));
	CHECK_FALSE(post_called);

	// Same for RELEASE.
	bool post_called_release = false;
	auto post_fn2 = [&post_called_release](const std::string &, const std::string &, const std::string &,
	                                       const std::string &, const std::string &) -> std::string {
		post_called_release = true;
		return "";
	};
	CHECK_THROWS_WITH(miint::SubmitLifecycle(miint::ENAAction::RELEASE, opts, post_fn2),
	                  Catch::Matchers::ContainsSubstring("hold_until_date"));
	CHECK_FALSE(post_called_release);
}
