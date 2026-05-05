// SPDX-License-Identifier: MIT
//
// Pre-INSERT alias collision check. Pure-data tests of the
// portal-API URL builder and the chunked-fetch driver, with a mock fetcher
// substituting for the live https://www.ebi.ac.uk/ena/portal/api/search
// endpoint.

#include "ena_alias_check.hpp"

#include <catch2/catch_all.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

using miint::AliasObjectKind;
using miint::AliasObjectKindFromTableName;
using miint::AliasObjectKindName;
using miint::BuildAliasCollisionURL;
using miint::CheckAliasCollisions;

namespace {

// Records every URL it receives and answers from a fixture map. URLs that
// don't match a fixture get a header-only TSV (i.e. no collision).
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
		return "alias\n";
	}

	int CountContaining(const std::string &needle) const {
		int n = 0;
		for (const auto &u : urls) {
			if (u.find(needle) != std::string::npos) {
				n++;
			}
		}
		return n;
	}
};

const std::string PORTAL_BASE = "https://www.ebi.ac.uk/ena/portal/api";

} // namespace

TEST_CASE("AliasObjectKindName returns the canonical lowercase token", "[ena_alias_check]") {
	CHECK(std::string(AliasObjectKindName(AliasObjectKind::STUDY)) == "study");
	CHECK(std::string(AliasObjectKindName(AliasObjectKind::SAMPLE)) == "sample");
	CHECK(std::string(AliasObjectKindName(AliasObjectKind::EXPERIMENT)) == "experiment");
	CHECK(std::string(AliasObjectKindName(AliasObjectKind::RUN)) == "run");
}

TEST_CASE("AliasObjectKindFromTableName covers virtual-catalog names", "[ena_alias_check]") {
	CHECK(AliasObjectKindFromTableName("projects") == AliasObjectKind::STUDY);
	CHECK(AliasObjectKindFromTableName("samples") == AliasObjectKind::SAMPLE);
	CHECK(AliasObjectKindFromTableName("experiments") == AliasObjectKind::EXPERIMENT);
	CHECK(AliasObjectKindFromTableName("runs") == AliasObjectKind::RUN);
	REQUIRE_THROWS_AS(AliasObjectKindFromTableName("unknown_table"), std::invalid_argument);
	REQUIRE_THROWS_AS(AliasObjectKindFromTableName(""), std::invalid_argument);
}

TEST_CASE("BuildAliasCollisionURL builds one-alias query for STUDY", "[ena_alias_check]") {
	const auto url = BuildAliasCollisionURL(PORTAL_BASE, AliasObjectKind::STUDY, {"my-study"});
	// result=read_study; field=study_alias; query=study_alias IN ("my-study"); format=tsv
	CHECK(url.find("https://www.ebi.ac.uk/ena/portal/api/search?") == 0);
	CHECK(url.find("result=read_study") != std::string::npos);
	CHECK(url.find("fields=study_alias") != std::string::npos);
	CHECK(url.find("study_alias%20IN%20%28%22my-study%22%29") != std::string::npos);
	CHECK(url.find("format=tsv") != std::string::npos);
}

TEST_CASE("BuildAliasCollisionURL builds correctly for each kind", "[ena_alias_check]") {
	struct Case {
		AliasObjectKind kind;
		const char *result;
		const char *field;
	};
	const std::vector<Case> cases = {
	    {AliasObjectKind::STUDY, "result=read_study", "fields=study_alias"},
	    {AliasObjectKind::SAMPLE, "result=read_sample", "fields=sample_alias"},
	    {AliasObjectKind::EXPERIMENT, "result=read_experiment", "fields=experiment_alias"},
	    {AliasObjectKind::RUN, "result=read_run", "fields=run_alias"},
	};
	for (const auto &c : cases) {
		const auto url = BuildAliasCollisionURL(PORTAL_BASE, c.kind, {"alpha"});
		INFO(url);
		CHECK(url.find(c.result) != std::string::npos);
		CHECK(url.find(c.field) != std::string::npos);
	}
}

TEST_CASE("BuildAliasCollisionURL chains multiple aliases via IN list", "[ena_alias_check]") {
	const auto url = BuildAliasCollisionURL(PORTAL_BASE, AliasObjectKind::SAMPLE, {"a1", "a2", "a3"});
	// IN ("a1","a2","a3") -> sample_alias%20IN%20%28%22a1%22%2C%22a2%22%2C%22a3%22%29
	CHECK(url.find("sample_alias%20IN%20%28%22a1%22%2C%22a2%22%2C%22a3%22%29") != std::string::npos);
}

TEST_CASE("BuildAliasCollisionURL rejects empty alias list", "[ena_alias_check]") {
	REQUIRE_THROWS_AS(BuildAliasCollisionURL(PORTAL_BASE, AliasObjectKind::SAMPLE, {}), std::invalid_argument);
}

TEST_CASE("BuildAliasCollisionURL rejects aliases with embedded quote/newline", "[ena_alias_check]") {
	REQUIRE_THROWS_AS(BuildAliasCollisionURL(PORTAL_BASE, AliasObjectKind::SAMPLE, {"bad\"alias"}),
	                  std::invalid_argument);
	REQUIRE_THROWS_AS(BuildAliasCollisionURL(PORTAL_BASE, AliasObjectKind::SAMPLE, {"line1\nline2"}),
	                  std::invalid_argument);
	REQUIRE_THROWS_AS(BuildAliasCollisionURL(PORTAL_BASE, AliasObjectKind::SAMPLE, {""}), std::invalid_argument);
}

TEST_CASE("CheckAliasCollisions: empty input → no fetch", "[ena_alias_check]") {
	MockFetcher fetch;
	auto hits = CheckAliasCollisions(PORTAL_BASE, AliasObjectKind::SAMPLE, {}, std::ref(fetch));
	CHECK(hits.empty());
	CHECK(fetch.urls.empty());
}

TEST_CASE("CheckAliasCollisions: header-only response → no collisions", "[ena_alias_check]") {
	MockFetcher fetch; // default response is "alias\n" (header only)
	auto hits = CheckAliasCollisions(PORTAL_BASE, AliasObjectKind::SAMPLE, {"fresh-1", "fresh-2"}, std::ref(fetch));
	CHECK(hits.empty());
	CHECK(fetch.urls.size() == 1);
}

TEST_CASE("CheckAliasCollisions: TSV with column header + collision rows returns hits", "[ena_alias_check]") {
	MockFetcher fetch;
	// Portal API returns a TSV with header row "<kind>_alias" and one alias per line.
	fetch.responses.push_back({"sample_alias%20IN", "sample_alias\nexisting-1\nexisting-3\n"});
	auto hits = CheckAliasCollisions(PORTAL_BASE, AliasObjectKind::SAMPLE, {"existing-1", "fresh-2", "existing-3"},
	                                 std::ref(fetch));
	REQUIRE(hits.size() == 2);
	// Returned in input order, not response order.
	CHECK(hits[0] == "existing-1");
	CHECK(hits[1] == "existing-3");
}

TEST_CASE("CheckAliasCollisions: ignores aliases the server returned that we didn't ask about", "[ena_alias_check]") {
	// Defensive: a malicious / stale portal response naming an alias outside our
	// query must not poison the result.
	MockFetcher fetch;
	fetch.responses.push_back({"sample_alias%20IN", "sample_alias\nexisting-1\nrogue-record\n"});
	auto hits = CheckAliasCollisions(PORTAL_BASE, AliasObjectKind::SAMPLE, {"existing-1", "fresh-2"}, std::ref(fetch));
	REQUIRE(hits.size() == 1);
	CHECK(hits[0] == "existing-1");
}

TEST_CASE("CheckAliasCollisions: dedups when an alias appears twice in the input", "[ena_alias_check]") {
	// Same alias twice: caller may do this if upstream produced duplicates. We
	// fetch once and report one collision (preserving first occurrence's
	// position).
	MockFetcher fetch;
	fetch.responses.push_back({"sample_alias%20IN", "sample_alias\ndupe\n"});
	auto hits = CheckAliasCollisions(PORTAL_BASE, AliasObjectKind::SAMPLE, {"dupe", "dupe"}, std::ref(fetch));
	REQUIRE(hits.size() == 1);
	CHECK(hits[0] == "dupe");
}

TEST_CASE("CheckAliasCollisions: chunks queries when alias list exceeds the per-URL limit", "[ena_alias_check]") {
	// 120 aliases at chunk-size 50 → ceil(120/50) = 3 fetches.
	std::vector<std::string> aliases;
	for (int i = 0; i < 120; i++) {
		aliases.push_back("alias-" + std::to_string(i));
	}
	MockFetcher fetch;
	// Mark alias-7, alias-77, alias-119 as colliding across different chunks.
	// Each fixture is anchored with a trailing `%2C` (URL-encoded comma) or
	// `%29` (URL-encoded close paren) so `%22alias-7%22` doesn't accidentally
	// match `%22alias-77%22` or `%22alias-78%22` if the chunk boundaries ever
	// shift.
	fetch.responses.push_back({"%22alias-7%22%2C", "sample_alias\nalias-7\n"});
	fetch.responses.push_back({"%22alias-77%22%2C", "sample_alias\nalias-77\n"});
	fetch.responses.push_back({"%22alias-119%22%29", "sample_alias\nalias-119\n"});

	auto hits = CheckAliasCollisions(PORTAL_BASE, AliasObjectKind::SAMPLE, aliases, std::ref(fetch));
	CHECK(fetch.urls.size() == 3);
	REQUIRE(hits.size() == 3);
	// Returned in input order.
	CHECK(hits[0] == "alias-7");
	CHECK(hits[1] == "alias-77");
	CHECK(hits[2] == "alias-119");
}

TEST_CASE("BuildAliasCollisionURL percent-encodes special characters in alias values", "[ena_alias_check]") {
	// `%`, `&`, `=` would corrupt the URL query parameters once the portal
	// API decodes them; they must be %HH-encoded inside the IN-list.
	const auto url = BuildAliasCollisionURL(PORTAL_BASE, AliasObjectKind::SAMPLE, {"weird%alias&v=1"});
	CHECK(url.find("%22weird%25alias%26v%3D1%22") != std::string::npos);
	// Common safe characters survive verbatim.
	const auto url2 = BuildAliasCollisionURL(PORTAL_BASE, AliasObjectKind::SAMPLE, {"safe-alias_1.0~tilde"});
	CHECK(url2.find("%22safe-alias_1.0~tilde%22") != std::string::npos);
}

TEST_CASE("CheckAliasCollisions: response matching uses the original alias, not the encoded form",
          "[ena_alias_check]") {
	// The portal API returns the alias in the response body in its original
	// (non-percent-encoded) form. Round-trip a `%`-bearing alias to confirm
	// we don't double-decode or accidentally compare encoded against decoded.
	MockFetcher fetch;
	fetch.responses.push_back({"%22weird%25alias%22", "sample_alias\nweird%alias\n"});
	auto hits =
	    CheckAliasCollisions(PORTAL_BASE, AliasObjectKind::SAMPLE, {"weird%alias", "fresh-other"}, std::ref(fetch));
	REQUIRE(hits.size() == 1);
	CHECK(hits[0] == "weird%alias");
}

TEST_CASE("CheckAliasCollisions: tolerates trailing CRLF + blank lines", "[ena_alias_check]") {
	MockFetcher fetch;
	fetch.responses.push_back({"sample_alias%20IN", "sample_alias\r\nexisting-1\r\n\r\n"});
	auto hits = CheckAliasCollisions(PORTAL_BASE, AliasObjectKind::SAMPLE, {"existing-1"}, std::ref(fetch));
	REQUIRE(hits.size() == 1);
	CHECK(hits[0] == "existing-1");
}
