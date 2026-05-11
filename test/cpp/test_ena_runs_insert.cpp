// Unit tests for the pure-data layer of INSERT INTO ena.runs: assemble
// envelope, POST via injected functor, parse receipt, return row-shaped
// result tuples (alias, err_accession, status). Mirrors the
// projects/samples/experiments insert tests (mock-fetcher pattern).

#include "ena_envelope_builder.hpp"
#include "ena_insert_test_helpers.hpp"
#include "ena_receipt_parser.hpp"
#include "ena_runs_insert.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <stdexcept>
#include <string>
#include <vector>

using namespace miint;
using miint_test::CapturedPost;

namespace {

std::string MakeRunReceipt(const std::vector<std::pair<std::string, std::string>> &alias_to_err, bool success = true,
                           const std::string &error = "") {
	std::vector<miint_test::ReceiptObjectFixture> objects;
	objects.reserve(alias_to_err.size());
	for (const auto &kv : alias_to_err) {
		objects.push_back({"RUN", kv.first, kv.second, "", ""});
	}
	return miint_test::MakeReceiptXML(objects, success, error);
}

RunSpec PairedRun(const std::string &alias, const std::string &experiment_alias) {
	RunSpec r;
	r.alias = alias;
	r.experiment_ref.refname = experiment_alias;
	r.files.push_back({alias + "_1.fastq.gz", "fastq", "9b8932f85caa54e687eba62fca3edce2"});
	r.files.push_back({alias + "_2.fastq.gz", "fastq", "183d6a24e0c3704e993bebe75bbbd989"});
	return r;
}

} // namespace

TEST_CASE("ENA runs insert: paired-fastq run builds envelope and parses receipt", "[ena_runs_insert]") {
	CapturedPost captured;
	auto post_fn = [&captured](const std::string &url, const std::string &body, const std::string &user,
	                           const std::string &password, const std::string &content_type) {
		captured = {url, body, user, password, content_type};
		return MakeRunReceipt({{"r1", "ERR42"}});
	};

	std::vector<RunSpec> runs = {PairedRun("r1", "e1")};
	ENARunInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	auto rows = SubmitRunInsert(runs, opts, post_fn);

	REQUIRE(rows.size() == 1);
	REQUIRE(rows[0].alias == "r1");
	REQUIRE(rows[0].err_accession == "ERR42");
	REQUIRE(rows[0].status == "PRIVATE");

	REQUIRE(captured.url == "http://mock.example/submit");
	REQUIRE(captured.user == "Webin-1");
	REQUIRE(captured.password == "pw");
	// V2 dispatch quirk: SRA-side objects (experiment / run) require XML.
	REQUIRE(captured.content_type.find("xml") != std::string::npos);
	REQUIRE(captured.body.find("<RUN alias=\"r1\">") != std::string::npos);
	REQUIRE(captured.body.find("<EXPERIMENT_REF refname=\"e1\"/>") != std::string::npos);
	REQUIRE(captured.body.find("checksum=\"9b8932f85caa54e687eba62fca3edce2\"") != std::string::npos);
	REQUIRE(captured.body.find("checksum=\"183d6a24e0c3704e993bebe75bbbd989\"") != std::string::npos);
}

TEST_CASE("ENA runs insert: multi-run preserves order", "[ena_runs_insert]") {
	auto post_fn = [](const std::string &, const std::string &body, const std::string &, const std::string &,
	                  const std::string &) {
		REQUIRE(body.find("<RUN alias=\"a\">") != std::string::npos);
		REQUIRE(body.find("<RUN alias=\"b\">") != std::string::npos);
		return MakeRunReceipt({{"a", "ERR10"}, {"b", "ERR11"}});
	};

	std::vector<RunSpec> runs = {PairedRun("a", "e1"), PairedRun("b", "e1")};
	ENARunInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	auto rows = SubmitRunInsert(runs, opts, post_fn);
	REQUIRE(rows.size() == 2);
	REQUIRE(rows[0].alias == "a");
	REQUIRE(rows[0].err_accession == "ERR10");
	REQUIRE(rows[1].alias == "b");
	REQUIRE(rows[1].err_accession == "ERR11");
}

TEST_CASE("ENA runs insert: receipt success=false throws with messages", "[ena_runs_insert]") {
	auto post_fn = [](const std::string &, const std::string &, const std::string &, const std::string &,
	                  const std::string &) {
		return MakeRunReceipt({}, false, "experiment reference 'e-missing' not found");
	};

	std::vector<RunSpec> runs = {PairedRun("r1", "e-missing")};
	ENARunInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	REQUIRE_THROWS_WITH(SubmitRunInsert(runs, opts, post_fn),
	                    Catch::Matchers::ContainsSubstring("experiment reference 'e-missing' not found"));
}

TEST_CASE("ENA runs insert: empty input is a no-op (no POST issued)", "[ena_runs_insert]") {
	bool called = false;
	auto post_fn = [&called](const std::string &, const std::string &, const std::string &, const std::string &,
	                         const std::string &) {
		called = true;
		return std::string();
	};

	std::vector<RunSpec> runs;
	ENARunInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	auto rows = SubmitRunInsert(runs, opts, post_fn);
	REQUIRE(rows.empty());
	REQUIRE_FALSE(called);
}

TEST_CASE("ENA runs insert: receipt missing alias is reported clearly", "[ena_runs_insert]") {
	auto post_fn = [](const std::string &, const std::string &, const std::string &, const std::string &,
	                  const std::string &) {
		return MakeRunReceipt({}, true);
	};

	std::vector<RunSpec> runs = {PairedRun("orphan", "e1")};
	ENARunInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	REQUIRE_THROWS_WITH(SubmitRunInsert(runs, opts, post_fn), Catch::Matchers::ContainsSubstring("orphan"));
}

TEST_CASE("ENA runs insert: MODIFY round-trips the user-supplied accession back into the row",
          "[ena_runs_insert][modify]") {
	CapturedPost captured;
	auto post_fn = [&captured](const std::string &url, const std::string &body, const std::string &user,
	                           const std::string &password, const std::string &content_type) {
		captured = {url, body, user, password, content_type};
		// Real wwwdev MODIFY echoes the user-supplied ERR verbatim. Mirror.
		return MakeRunReceipt({{"r1", "ERR9999100"}});
	};

	auto run = PairedRun("r1", "e1");
	run.accession = "ERR9999100";
	run.title = "Updated run title";
	std::vector<RunSpec> runs = {run};

	ENARunInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.action = ENAAction::MODIFY;

	auto outcome = SubmitRunInsertOutcome(runs, opts, post_fn);
	REQUIRE(outcome.success);
	REQUIRE(outcome.rows.size() == 1);
	REQUIRE(outcome.rows[0].alias == "r1");
	REQUIRE(outcome.rows[0].err_accession == "ERR9999100");
	REQUIRE(captured.body.find("<MODIFY/>") != std::string::npos);
	REQUIRE(captured.body.find(R"X(<RUN alias="r1" accession="ERR9999100">)X") != std::string::npos);
	// Pin TITLE emission so a future AppendXmlRun regression that drops the
	// optional element doesn't slip past the envelope-only test layer.
	REQUIRE(captured.body.find("<TITLE>Updated run title</TITLE>") != std::string::npos);
}

TEST_CASE("ENA runs insert: MODIFY with accession-form experiment_ref emits accession on EXPERIMENT_REF",
          "[ena_runs_insert][modify]") {
	// Production users will typically pull `experiment_ref` from
	// `ena.submission_log` after the experiment ADD round-trip; that's an
	// ERX accession, not a refname. The refname-form path is covered by
	// the basic happy path; this case pins the accession-form wire.
	CapturedPost captured;
	auto post_fn = [&captured](const std::string &url, const std::string &body, const std::string &user,
	                           const std::string &password, const std::string &content_type) {
		captured = {url, body, user, password, content_type};
		return MakeRunReceipt({{"r1", "ERR9999100"}});
	};

	auto run = PairedRun("r1", "e1");
	run.experiment_ref.accession = "ERX42";
	run.experiment_ref.refname.clear();
	run.accession = "ERR9999100";
	std::vector<RunSpec> runs = {run};

	ENARunInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.action = ENAAction::MODIFY;

	auto outcome = SubmitRunInsertOutcome(runs, opts, post_fn);
	REQUIRE(outcome.success);
	REQUIRE(captured.body.find(R"X(<EXPERIMENT_REF accession="ERX42"/>)X") != std::string::npos);
	// Refname must NOT appear when accession is set (RefDescriptor's
	// accession-wins precedence).
	REQUIRE(captured.body.find("refname=") == std::string::npos);
}

TEST_CASE("ENA runs insert: MODIFY failure receipt surfaces the server error", "[ena_runs_insert][modify]") {
	auto post_fn = [](const std::string &, const std::string &, const std::string &, const std::string &,
	                  const std::string &) {
		return MakeRunReceipt({}, false, "ERR00000 not found in submission account");
	};
	auto run = PairedRun("r1", "e1");
	run.accession = "ERR00000";
	std::vector<RunSpec> runs = {run};

	ENARunInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";
	opts.action = ENAAction::MODIFY;

	auto outcome = SubmitRunInsertOutcome(runs, opts, post_fn);
	REQUIRE_FALSE(outcome.success);
	REQUIRE(outcome.rows.empty());
	REQUIRE(outcome.error_messages.size() == 1);
	REQUIRE_THAT(outcome.error_messages[0], Catch::Matchers::ContainsSubstring("ERR00000 not found"));
}
