// SPDX-License-Identifier: MIT
//
// Phase 3 RED -> GREEN tests for the ENA Webin V2 JSON submission envelope
// builder. Fixtures derive from the canonical ENA tutorial format documented
// in localdocs/ena-research-webin-v2-deep.md §1 and the V2 OpenAPI behaviour
// (Accept-driven response shape; request body is opaque string).

#include "ena_envelope_builder.hpp"

#include <catch2/catch_all.hpp>

namespace {

// Convenience: compact-JSON byte-equal compare.
void CheckEqual(const std::string &actual, const std::string &expected) {
	CHECK(actual == expected);
}

} // namespace

TEST_CASE("ENA envelope: VALIDATE action only, no objects", "[ena_envelope]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::VALIDATE;
	auto json = miint::BuildEnvelopeJSON(env);
	CheckEqual(json, R"X({"submission":{"actions":[{"type":"VALIDATE"}]}})X");
}

TEST_CASE("ENA envelope: ADD with HOLD until a future date", "[ena_envelope]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::ADD;
	env.hold_until_date = "2027-01-01";
	auto json = miint::BuildEnvelopeJSON(env);
	CheckEqual(json, R"X({"submission":{"actions":[{"type":"ADD"},{"type":"HOLD","holdUntilDate":"2027-01-01"}]}})X");
}

TEST_CASE("ENA envelope: single submission_project (PRJEB will be issued)", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ProjectSpec p;
	p.alias = "gut-cohort-2026";
	p.title = "Adult gut microbiome cohort";
	p.description = "Phase 1 collection";
	p.project_type = "METAGENOMIC";
	env.projects.push_back(p);

	auto json = miint::BuildEnvelopeJSON(env);
	CheckEqual(json, R"X({"submission":{"actions":[{"type":"ADD"}]},)X"
	                 R"X("projects":[{"alias":"gut-cohort-2026","title":"Adult gut microbiome cohort",)X"
	                 R"X("description":"Phase 1 collection","sequencingProject":{}}]})X");
}

TEST_CASE("ENA envelope: project with no description omits the field", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ProjectSpec p;
	p.alias = "x";
	p.title = "y";
	p.project_type = "METAGENOMIC";
	env.projects.push_back(p);

	auto json = miint::BuildEnvelopeJSON(env);
	// Compact form — note no "description" key
	CHECK(json.find("\"description\"") == std::string::npos);
	CHECK(json.find("\"sequencingProject\":{}") != std::string::npos);
}

TEST_CASE("ENA envelope: umbrella project uses umbrellaProject marker", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ProjectSpec p;
	p.alias = "consortium-2026";
	p.title = "Consortium umbrella";
	p.is_umbrella = true;
	env.projects.push_back(p);

	auto json = miint::BuildEnvelopeJSON(env);
	CHECK(json.find("\"umbrellaProject\":{}") != std::string::npos);
	CHECK(json.find("sequencingProject") == std::string::npos);
}

TEST_CASE("ENA envelope: single sample with checklist and attributes", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::SampleSpec s;
	s.alias = "gut-001";
	s.title = "Adult subject 001 stool";
	s.taxon_id = 408170;
	s.checklist = "ERC000015";
	// Stable, deterministic order — builder emits attributes in caller's order.
	s.attributes.emplace_back("collection date", "2024-06-15");
	s.attributes.emplace_back("geographic location (country and/or sea)", "United States");
	env.samples.push_back(s);

	auto json = miint::BuildEnvelopeJSON(env);
	// Sanity: contains the right top-level keys
	CHECK(json.find(R"X("submission":{"actions":[{"type":"ADD"}]})X") != std::string::npos);
	CHECK(json.find(R"X("samples":)X") != std::string::npos);
	// Sample shape: alias, title, organism.taxonId (string per V2 convention)
	CHECK(json.find(R"X("alias":"gut-001")X") != std::string::npos);
	CHECK(json.find(R"X("title":"Adult subject 001 stool")X") != std::string::npos);
	CHECK(json.find(R"X("organism":{"taxonId":"408170"})X") != std::string::npos);
	// Checklist becomes an attribute with tag="ENA-CHECKLIST"
	CHECK(json.find(R"X({"tag":"ENA-CHECKLIST","value":"ERC000015"})X") != std::string::npos);
	// User attributes preserve order
	auto col_pos = json.find(R"X({"tag":"collection date","value":"2024-06-15"})X");
	auto geo_pos = json.find(R"X({"tag":"geographic location (country and/or sea)","value":"United States"})X");
	REQUIRE(col_pos != std::string::npos);
	REQUIRE(geo_pos != std::string::npos);
	CHECK(col_pos < geo_pos);
}

TEST_CASE("ENA envelope: multi-sample emits one element per sample", "[ena_envelope]") {
	miint::SubmissionSpec env;
	for (int i = 1; i <= 3; ++i) {
		miint::SampleSpec s;
		s.alias = "sample-" + std::to_string(i);
		s.taxon_id = 408170;
		s.checklist = "ERC000015";
		env.samples.push_back(s);
	}
	auto json = miint::BuildEnvelopeJSON(env);
	// Three sample objects, each with a unique alias
	CHECK(json.find(R"X("alias":"sample-1")X") != std::string::npos);
	CHECK(json.find(R"X("alias":"sample-2")X") != std::string::npos);
	CHECK(json.find(R"X("alias":"sample-3")X") != std::string::npos);
}

TEST_CASE("ENA envelope: project + sample mixed bundle", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ProjectSpec p;
	p.alias = "p1";
	p.title = "P1";
	p.project_type = "METAGENOMIC";
	env.projects.push_back(p);
	miint::SampleSpec s;
	s.alias = "s1";
	s.taxon_id = 408170;
	env.samples.push_back(s);

	auto json = miint::BuildEnvelopeJSON(env);
	// Both arrays present, in spec'd order: submission, projects, samples
	auto sub_pos = json.find(R"X("submission":)X");
	auto proj_pos = json.find(R"X("projects":)X");
	auto samp_pos = json.find(R"X("samples":)X");
	REQUIRE(sub_pos != std::string::npos);
	REQUIRE(proj_pos != std::string::npos);
	REQUIRE(samp_pos != std::string::npos);
	CHECK(sub_pos < proj_pos);
	CHECK(proj_pos < samp_pos);
}

TEST_CASE("ENA envelope: JSON string escaping (quote, backslash, control chars)", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ProjectSpec p;
	p.alias = "x";
	// Deliberately nasty: contains " \ \n \t and a high-bit byte (passes through as UTF-8 byte).
	p.title = std::string("Quoted \"thing\" \\ with\nnewline\tand tab");
	p.project_type = "METAGENOMIC";
	env.projects.push_back(p);

	auto json = miint::BuildEnvelopeJSON(env);
	// Each forbidden character must be escaped per RFC 8259 §7.
	CHECK(json.find(R"X("title":"Quoted \"thing\" \\ with\nnewline\tand tab")X") != std::string::npos);
}

TEST_CASE("ENA envelope: empty project alias is rejected at build time", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ProjectSpec p;
	p.title = "x";
	p.project_type = "METAGENOMIC";
	env.projects.push_back(p);
	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("alias"));
}

TEST_CASE("ENA envelope: sample with taxon_id <= 0 rejected", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::SampleSpec s;
	s.alias = "s1";
	s.taxon_id = 0;
	env.samples.push_back(s);
	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("taxon"));
}

TEST_CASE("ENA envelope: control characters escape as \\u00XX", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ProjectSpec p;
	p.alias = "x";
	// 0x01 (SOH) and 0x1F (US) — must escape as  /  per RFC 8259 §7.
	p.title = std::string("a\x01"
	                      "b\x1f"
	                      "c");
	p.project_type = "METAGENOMIC";
	env.projects.push_back(p);

	auto json = miint::BuildEnvelopeJSON(env);
	// Look for the literal six-char escape sequence in the produced JSON.
	CHECK(json.find("\\u0001") != std::string::npos);
	CHECK(json.find("\\u001f") != std::string::npos);
}

TEST_CASE("ENA envelope: action=HOLD without hold_until_date is rejected", "[ena_envelope]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::HOLD;
	// hold_until_date deliberately empty
	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("hold_until_date"));
}

TEST_CASE("ENA envelope: action=HOLD AND hold_until_date is rejected (avoid double-HOLD)", "[ena_envelope]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::HOLD;
	env.hold_until_date = "2027-01-01";
	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("automatically"));
}
