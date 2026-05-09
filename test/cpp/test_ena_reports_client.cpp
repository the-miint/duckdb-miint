// SPDX-License-Identifier: MIT
//
// Pure-data tests for the Webin Reports API client used by L5 to translate
// user-supplied aliases into server-assigned accessions before lifecycle
// envelope build. URL construction + JSON parse + miss/hit/transport-error
// behaviour all driven by a mock URLFetcher (the same testability pattern
// `ena_alias_check` uses).
//
// Wire shape pinned by `localdocs/ena-live-reports-probe.sh` (live-verified
// 2026-05-08 against wwwdev):
//   GET <base>/<kind>/{id}?format=json
//     hit  → HTTP 200 + [{"report":{"id":"PRJEB…","alias":"…", …}}]
//     miss → HTTP 200 + []
// The {id} segment accepts BOTH alias and primary accession; the server
// resolves either to the same record. Lifecycle ops always target the
// primary accession (decision #16), so for each AliasObjectKind we hit the
// path that returns the primary accession in `report.id`:
//   STUDY      → /projects (id = PRJEB…)
//   SAMPLE     → /samples  (id = ERS…)
//   EXPERIMENT → /experiments (id = ERX…)
//   RUN        → /runs (id = ERR…)
// (NOT /studies — that returns id = ERP, the secondary accession.)

#include "ena_reports_client.hpp"

#include <catch2/catch_all.hpp>

#include <stdexcept>
#include <string>
#include <vector>

using miint::AliasObjectKind;
using miint::BuildReportsLookupURL;
using miint::DEFAULT_REPORTS_BASE_PROD;
using miint::DEFAULT_REPORTS_BASE_TEST;
using miint::LookupAccessionByAlias;
using miint::ResolveReportsBaseForEndpoint;

namespace {

// Records every URL it receives and answers from a fixture map. URLs that
// don't match a fixture get the canonical "alias not found" response shape
// (HTTP 200, empty array). Mirrors the `ena_alias_check` MockFetcher.
struct MockFetcher {
	std::vector<std::string> urls;
	std::vector<std::pair<std::string, std::string>> responses;

	std::string operator()(const std::string &url) {
		urls.push_back(url);
		for (const auto &p : responses) {
			if (url.find(p.first) != std::string::npos) {
				return p.second;
			}
		}
		// Default: live wwwdev returns HTTP 200 with `[]` for any {id} the
		// account doesn't own. The fetcher contract is "throw on transport
		// failure"; HTTP 200 + empty array is a normal response, so we
		// echo it.
		return "[]";
	}
};

const std::string REPORTS_BASE = "https://wwwdev.ebi.ac.uk/ena/submit/report";

// Small helper to build a hit-response from a (alias, primary-accession) pair.
// The live wwwdev shape is `[{"report":{"id":"…","alias":"…", …}}]`; we keep
// a few extra fields to mirror the real payload so tests reading additional
// fields get realistic input.
std::string HitResponse(const std::string &kind_id, const std::string &alias) {
	return "[{\"report\":{\"id\":\"" + kind_id + "\",\"alias\":\"" + alias +
	       "\",\"firstCreated\":\"2026-05-08T03:27:35\",\"releaseStatus\":\"PRIVATE\"," +
	       "\"submissionAccountId\":\"Webin-34\"},\"links\":[]}]";
}

} // namespace

// -----------------------------------------------------------------------------
// URL construction
// -----------------------------------------------------------------------------

TEST_CASE("BuildReportsLookupURL maps each kind to the primary-accession path", "[ena_reports_client]") {
	// L5 design: STUDY routes to /projects (id = PRJEB), NOT /studies (id = ERP).
	// /studies would resolve aliases too, but its `id` field is the secondary
	// accession; lifecycle targets the primary. Keeping these in lock-step is
	// the load-bearing invariant of L5.
	struct Case {
		AliasObjectKind kind;
		const char *path_segment;
	};
	const std::vector<Case> cases = {
	    {AliasObjectKind::STUDY, "/projects/"},
	    {AliasObjectKind::SAMPLE, "/samples/"},
	    {AliasObjectKind::EXPERIMENT, "/experiments/"},
	    {AliasObjectKind::RUN, "/runs/"},
	};
	for (const auto &c : cases) {
		const auto url = BuildReportsLookupURL(REPORTS_BASE, c.kind, "alpha");
		INFO(url);
		CHECK(url.find(REPORTS_BASE) == 0);
		CHECK(url.find(c.path_segment) != std::string::npos);
		// Always ask for JSON — the parser only consumes the JSON shape.
		CHECK(url.find("format=json") != std::string::npos);
		// The {id} segment is always alpha here.
		CHECK(url.find("/alpha?") != std::string::npos);
	}
}

TEST_CASE("BuildReportsLookupURL strips a single trailing slash from the base", "[ena_reports_client]") {
	const auto a = BuildReportsLookupURL(REPORTS_BASE, AliasObjectKind::SAMPLE, "alpha");
	const auto b = BuildReportsLookupURL(REPORTS_BASE + "/", AliasObjectKind::SAMPLE, "alpha");
	CHECK(a == b);
}

TEST_CASE("BuildReportsLookupURL percent-encodes the {id} path segment", "[ena_reports_client]") {
	// User aliases may contain spaces / ampersands / equals if a careless
	// caller uses them; each survives the round-trip via percent encoding so
	// the server sees the literal characters.
	const auto url = BuildReportsLookupURL(REPORTS_BASE, AliasObjectKind::SAMPLE, "with space&other=val");
	INFO(url);
	CHECK(url.find("with%20space%26other%3Dval") != std::string::npos);
}

TEST_CASE("BuildReportsLookupURL rejects empty {id}", "[ena_reports_client]") {
	REQUIRE_THROWS_AS(BuildReportsLookupURL(REPORTS_BASE, AliasObjectKind::SAMPLE, ""), std::invalid_argument);
}

TEST_CASE("BuildReportsLookupURL rejects {id} with newline", "[ena_reports_client]") {
	// HTTP requests would split on a newline (header injection). Reject up
	// front rather than relying on the URL encoder to sanitise — the encoder
	// would produce %0A but a defensive caller is clearer.
	REQUIRE_THROWS_AS(BuildReportsLookupURL(REPORTS_BASE, AliasObjectKind::SAMPLE, "line1\nline2"),
	                  std::invalid_argument);
	REQUIRE_THROWS_AS(BuildReportsLookupURL(REPORTS_BASE, AliasObjectKind::SAMPLE, "line1\rline2"),
	                  std::invalid_argument);
}

// -----------------------------------------------------------------------------
// Lookup driver — single-call, hit / miss / parse
// -----------------------------------------------------------------------------

TEST_CASE("LookupAccessionByAlias: hit returns report.id from the JSON array", "[ena_reports_client]") {
	MockFetcher fetch;
	fetch.responses.push_back({"/projects/my-study", HitResponse("PRJEB99999", "my-study")});
	const auto acc = LookupAccessionByAlias(REPORTS_BASE, AliasObjectKind::STUDY, "my-study", std::ref(fetch));
	CHECK(acc == "PRJEB99999");
	REQUIRE(fetch.urls.size() == 1);
	CHECK(fetch.urls[0].find("/projects/my-study") != std::string::npos);
}

TEST_CASE("LookupAccessionByAlias: hit on each kind", "[ena_reports_client]") {
	struct Case {
		AliasObjectKind kind;
		const char *path;
		const char *primary;
	};
	const std::vector<Case> cases = {
	    {AliasObjectKind::STUDY, "/projects/", "PRJEB12345"},
	    {AliasObjectKind::SAMPLE, "/samples/", "ERS12345"},
	    {AliasObjectKind::EXPERIMENT, "/experiments/", "ERX12345"},
	    {AliasObjectKind::RUN, "/runs/", "ERR12345"},
	};
	for (const auto &c : cases) {
		MockFetcher fetch;
		fetch.responses.push_back({c.path, HitResponse(c.primary, "alpha")});
		const auto acc = LookupAccessionByAlias(REPORTS_BASE, c.kind, "alpha", std::ref(fetch));
		INFO("kind path " << c.path);
		CHECK(acc == c.primary);
	}
}

TEST_CASE("LookupAccessionByAlias: miss returns empty string (HTTP 200, body=[])", "[ena_reports_client]") {
	// Live wwwdev returns HTTP 200 + `[]` for an alias the account doesn't
	// own. We treat empty array as "not found" and return empty string;
	// callers (the lifecycle table fns) surface this as a friendly
	// "alias not found" diagnostic rather than throwing.
	MockFetcher fetch; // default response is "[]"
	const auto acc = LookupAccessionByAlias(REPORTS_BASE, AliasObjectKind::STUDY, "not-mine", std::ref(fetch));
	CHECK(acc.empty());
	CHECK(fetch.urls.size() == 1);
}

TEST_CASE("LookupAccessionByAlias: tolerates whitespace + trailing newline in JSON", "[ena_reports_client]") {
	// A pretty-printed body should still parse — Spring REST sometimes adds
	// a trailing newline; our parser must not be position-sensitive.
	MockFetcher fetch;
	fetch.responses.push_back(
	    {"/projects/alpha", "  [  { \"report\" :  { \"id\"  : \"PRJEB1\" , \"alias\":\"alpha\" } } ]  \n"});
	const auto acc = LookupAccessionByAlias(REPORTS_BASE, AliasObjectKind::STUDY, "alpha", std::ref(fetch));
	CHECK(acc == "PRJEB1");
}

TEST_CASE("LookupAccessionByAlias: skips non-string sibling values when scanning for 'id'", "[ena_reports_client]") {
	// Live wwwdev only carries string-typed fields at the outer object level
	// today, but the inner `report` may carry numeric / null / boolean /
	// nested array / nested object siblings if ENA extends the schema.
	// SkipJsonValue must handle all four JSON value classes (string is
	// covered by the existing tests). Earlier impl had a "consumed_first_token"
	// gate that silently mis-skipped bare scalars at depth 0 — pin all
	// classes here so a regression that mishandles any class surfaces
	// directly. Test order: number, null, bool, nested array, nested
	// object, all preceding the `id` we want.
	MockFetcher fetch;
	fetch.responses.push_back({"/samples/alpha", "[{\"report\":{\"taxId\":408170,"
	                                             "\"commonName\":null,"
	                                             "\"isReleased\":false,"
	                                             "\"externalIds\":[\"SAMEA1\",\"SAMEA2\"],"
	                                             "\"submitter\":{\"id\":\"X\",\"role\":\"PI\"},"
	                                             "\"id\":\"ERS9\","
	                                             "\"alias\":\"alpha\"}}]"});
	const auto acc = LookupAccessionByAlias(REPORTS_BASE, AliasObjectKind::SAMPLE, "alpha", std::ref(fetch));
	CHECK(acc == "ERS9");
}

TEST_CASE("LookupAccessionByAlias: throws on transport error (fetcher exception)", "[ena_reports_client]") {
	struct ThrowingFetcher {
		std::string operator()(const std::string &) {
			throw std::runtime_error("network down");
		}
	};
	ThrowingFetcher fetch;
	REQUIRE_THROWS_AS(LookupAccessionByAlias(REPORTS_BASE, AliasObjectKind::STUDY, "alpha", std::ref(fetch)),
	                  std::runtime_error);
}

TEST_CASE("LookupAccessionByAlias: malformed JSON surfaces as runtime_error", "[ena_reports_client]") {
	MockFetcher fetch;
	fetch.responses.push_back({"/projects/alpha", "not even close to JSON"});
	REQUIRE_THROWS_AS(LookupAccessionByAlias(REPORTS_BASE, AliasObjectKind::STUDY, "alpha", std::ref(fetch)),
	                  std::runtime_error);
}

TEST_CASE("LookupAccessionByAlias: array element missing report.id is a parse error", "[ena_reports_client]") {
	// The server contract is "id is always present on a hit". If the response
	// has an entry but no id, that's a server-side bug we should surface
	// loudly rather than silently treat as miss.
	MockFetcher fetch;
	fetch.responses.push_back({"/projects/alpha", "[{\"report\":{\"alias\":\"alpha\"}}]"});
	REQUIRE_THROWS_AS(LookupAccessionByAlias(REPORTS_BASE, AliasObjectKind::STUDY, "alpha", std::ref(fetch)),
	                  std::runtime_error);
}

TEST_CASE("LookupAccessionByAlias: empty alias is a programming error", "[ena_reports_client]") {
	// Same defensive contract as BuildReportsLookupURL — surfaces at the
	// pure-data layer so the lifecycle table fn can rely on bind-time
	// whitespace guards.
	MockFetcher fetch;
	REQUIRE_THROWS_AS(LookupAccessionByAlias(REPORTS_BASE, AliasObjectKind::STUDY, "", std::ref(fetch)),
	                  std::invalid_argument);
}

// -----------------------------------------------------------------------------
// ResolveReportsBaseForEndpoint — env-var override and endpoint-aware default
// (POSIX setenv/unsetenv: gate on non-WASM since Emscripten's libc doesn't
// expose them; matches the project-wide WASM-first-class build constraint).
// -----------------------------------------------------------------------------

#ifndef __EMSCRIPTEN__
TEST_CASE("ResolveReportsBaseForEndpoint defaults to wwwdev for non-production endpoints", "[ena_reports_client]") {
	::unsetenv("MIINT_ENA_REPORTS_URL_BASE");
	// Mirrors ResolveDefaultENAEndpointURL: only "production" routes to
	// the prod base; everything else (including "test", empty, or any
	// unrecognised label) defaults to wwwdev.
	CHECK(ResolveReportsBaseForEndpoint("production") == DEFAULT_REPORTS_BASE_PROD);
	CHECK(ResolveReportsBaseForEndpoint("test") == DEFAULT_REPORTS_BASE_TEST);
	CHECK(ResolveReportsBaseForEndpoint("") == DEFAULT_REPORTS_BASE_TEST);
	CHECK(ResolveReportsBaseForEndpoint("anything-else") == DEFAULT_REPORTS_BASE_TEST);
}

TEST_CASE("ResolveReportsBaseForEndpoint honours MIINT_ENA_REPORTS_URL_BASE override", "[ena_reports_client]") {
	// Env var wins regardless of endpoint label — the override is the test
	// fixture mechanism (mock server / private wwwdev mirror) and must
	// short-circuit the per-endpoint default.
	::setenv("MIINT_ENA_REPORTS_URL_BASE", "https://example.test/report", /*overwrite=*/1);
	CHECK(ResolveReportsBaseForEndpoint("production") == "https://example.test/report");
	CHECK(ResolveReportsBaseForEndpoint("test") == "https://example.test/report");
	// Trailing slash stripped — same convention as the portal base.
	::setenv("MIINT_ENA_REPORTS_URL_BASE", "https://example.test/report/", /*overwrite=*/1);
	CHECK(ResolveReportsBaseForEndpoint("test") == "https://example.test/report");
	::unsetenv("MIINT_ENA_REPORTS_URL_BASE");
}
#endif
