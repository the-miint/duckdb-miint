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
