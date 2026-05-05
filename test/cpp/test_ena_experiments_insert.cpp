// Unit tests for the pure-data layer of INSERT INTO ena.experiments: assemble
// envelope, POST via injected functor, parse receipt, return row-shaped
// result tuples (alias, erx_accession, status). Mirrors the projects/samples
// insert tests (mock-fetcher pattern).

#include "ena_envelope_builder.hpp"
#include "ena_experiments_insert.hpp"
#include "ena_insert_test_helpers.hpp"
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

std::string MakeExperimentReceipt(const std::vector<std::pair<std::string, std::string>> &alias_to_erx,
                                  bool success = true, const std::string &error = "") {
	std::vector<miint_test::ReceiptObjectFixture> objects;
	objects.reserve(alias_to_erx.size());
	for (const auto &kv : alias_to_erx) {
		objects.push_back({"EXPERIMENT", kv.first, kv.second, "", ""});
	}
	return miint_test::MakeReceiptXML(objects, success, error);
}

ExperimentSpec MinimalExperiment(const std::string &alias) {
	ExperimentSpec e;
	e.alias = alias;
	e.study_ref.refname = "p1";
	e.sample_ref.refname = "s1";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	// SRA `typeIlluminaModel` enum requires the `Illumina` prefix
	// (verified live 2026-05-04 — bare "NovaSeq 6000" is rejected).
	e.instrument_model = "Illumina NovaSeq 6000";
	return e;
}

} // namespace

TEST_CASE("ENA experiments insert: single row builds envelope and parses receipt", "[ena_experiments_insert]") {
	CapturedPost captured;
	auto post_fn = [&captured](const std::string &url, const std::string &body, const std::string &user,
	                           const std::string &password, const std::string &content_type) {
		captured = {url, body, user, password, content_type};
		return MakeExperimentReceipt({{"e1", "ERX42"}});
	};

	std::vector<ExperimentSpec> exps = {MinimalExperiment("e1")};
	ENAExperimentInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	auto rows = SubmitExperimentInsert(exps, opts, post_fn);

	REQUIRE(rows.size() == 1);
	REQUIRE(rows[0].alias == "e1");
	REQUIRE(rows[0].erx_accession == "ERX42");
	REQUIRE(rows[0].status == "PRIVATE");

	REQUIRE(captured.url == "http://mock.example/submit");
	REQUIRE(captured.user == "Webin-1");
	REQUIRE(captured.password == "pw");
	// V2 dispatch quirk: SRA-side objects (experiment / run) require XML.
	REQUIRE(captured.content_type.find("xml") != std::string::npos);
	REQUIRE(captured.body.find("<EXPERIMENT alias=\"e1\">") != std::string::npos);
	REQUIRE(captured.body.find("<LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>") != std::string::npos);
	REQUIRE(captured.body.find("<INSTRUMENT_MODEL>Illumina NovaSeq 6000</INSTRUMENT_MODEL>") != std::string::npos);
}

TEST_CASE("ENA experiments insert: multi-row preserves order", "[ena_experiments_insert]") {
	auto post_fn = [](const std::string &, const std::string &body, const std::string &, const std::string &,
	                  const std::string &) {
		REQUIRE(body.find("<EXPERIMENT alias=\"a\">") != std::string::npos);
		REQUIRE(body.find("<EXPERIMENT alias=\"b\">") != std::string::npos);
		return MakeExperimentReceipt({{"a", "ERX10"}, {"b", "ERX11"}});
	};

	std::vector<ExperimentSpec> exps = {MinimalExperiment("a"), MinimalExperiment("b")};
	ENAExperimentInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	auto rows = SubmitExperimentInsert(exps, opts, post_fn);
	REQUIRE(rows.size() == 2);
	REQUIRE(rows[0].alias == "a");
	REQUIRE(rows[0].erx_accession == "ERX10");
	REQUIRE(rows[1].alias == "b");
	REQUIRE(rows[1].erx_accession == "ERX11");
}

TEST_CASE("ENA experiments insert: receipt success=false throws with messages", "[ena_experiments_insert]") {
	auto post_fn = [](const std::string &, const std::string &, const std::string &, const std::string &,
	                  const std::string &) {
		return MakeExperimentReceipt({}, false, "study reference 'p-missing' not found");
	};

	std::vector<ExperimentSpec> exps = {MinimalExperiment("e1")};
	ENAExperimentInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	REQUIRE_THROWS_WITH(SubmitExperimentInsert(exps, opts, post_fn),
	                    Catch::Matchers::ContainsSubstring("study reference 'p-missing' not found"));
}

TEST_CASE("ENA experiments insert: empty input is a no-op (no POST issued)", "[ena_experiments_insert]") {
	bool called = false;
	auto post_fn = [&called](const std::string &, const std::string &, const std::string &, const std::string &,
	                         const std::string &) {
		called = true;
		return std::string();
	};

	std::vector<ExperimentSpec> exps;
	ENAExperimentInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	auto rows = SubmitExperimentInsert(exps, opts, post_fn);
	REQUIRE(rows.empty());
	REQUIRE_FALSE(called);
}

TEST_CASE("ENA experiments insert: receipt missing alias is reported clearly", "[ena_experiments_insert]") {
	auto post_fn = [](const std::string &, const std::string &, const std::string &, const std::string &,
	                  const std::string &) {
		return MakeExperimentReceipt({}, true);
	};

	std::vector<ExperimentSpec> exps = {MinimalExperiment("orphan")};
	ENAExperimentInsertOptions opts;
	opts.endpoint_url = "http://mock.example/submit";
	opts.user = "Webin-1";
	opts.password = "pw";

	REQUIRE_THROWS_WITH(SubmitExperimentInsert(exps, opts, post_fn), Catch::Matchers::ContainsSubstring("orphan"));
}
