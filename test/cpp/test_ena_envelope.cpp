// SPDX-License-Identifier: MIT
//
// Tests for the ENA Webin V2 JSON + XML submission envelope builders.
// Fixtures derive from the canonical ENA tutorial format and the V2 OpenAPI
// behaviour (Accept-driven response shape; request body is opaque string).

#include "ena_envelope_builder.hpp"

#include <catch2/catch_all.hpp>

#include <algorithm>

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

TEST_CASE("ENA envelope: project with no description falls back to the title", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ProjectSpec p;
	p.alias = "x";
	p.title = "y";
	p.project_type = "METAGENOMIC";
	env.projects.push_back(p);

	auto json = miint::BuildEnvelopeJSON(env);
	// `description` is always emitted; when the user didn't provide one we
	// reuse `title` so the XML <DESCRIPTION> element has non-empty content
	// (wwwdev intermittently rejects PROJECT documents missing a populated
	// DESCRIPTION even though the XSD says minOccurs=0).
	CHECK(json.find("\"description\":\"y\"") != std::string::npos);
	CHECK(json.find("\"sequencingProject\":{}") != std::string::npos);
}

TEST_CASE("ENA envelope: MODIFY project emits the existing accession on the project element",
          "[ena_envelope][modify]") {
	// MODIFY in Webin V2 is "re-submit the full updated XML with the
	// existing accession" (`<PROJECT alias="..." accession="PRJEB...">`).
	// The accession identifies which already-registered object to update;
	// without it the server has no way to disambiguate against the alias's
	// possible reuse.
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::MODIFY;
	miint::ProjectSpec p;
	p.alias = "gut-cohort-2026";
	p.accession = "PRJEB123456";
	p.title = "Adult gut microbiome cohort (revised)";
	p.description = "Phase 1 collection — updated abstract";
	p.project_type = "METAGENOMIC";
	env.projects.push_back(p);

	auto json = miint::BuildEnvelopeJSON(env);
	CHECK(json.find("\"type\":\"MODIFY\"") != std::string::npos);
	// Accession must appear right next to the alias (same JSON object) so
	// the per-element binding is unambiguous to any reasonable parser.
	CHECK(json.find("\"alias\":\"gut-cohort-2026\",\"accession\":\"PRJEB123456\"") != std::string::npos);
	CHECK(json.find("\"title\":\"Adult gut microbiome cohort (revised)\"") != std::string::npos);
}

TEST_CASE("ENA envelope: ADD with accession on a project is rejected at build time", "[ena_envelope][modify]") {
	// On ADD the server assigns the accession; setting one on the spec is
	// almost certainly a programmer error (mistakenly reusing a MODIFY-shaped
	// spec on the ADD path). The server-side rejection would be confusing
	// ("alias exists with accession X"); throw here so the diagnostic names
	// the offending alias and the cause directly.
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::ADD;
	miint::ProjectSpec p;
	p.alias = "fresh-project";
	p.accession = "PRJEB999"; // bogus pre-fill
	p.title = "t";
	p.project_type = "METAGENOMIC";
	env.projects.push_back(p);

	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env),
	                  Catch::Matchers::ContainsSubstring("ADD must not set an accession") &&
	                      Catch::Matchers::ContainsSubstring("fresh-project"));
}

TEST_CASE("ENA envelope: MODIFY without accession on a project is rejected at build time", "[ena_envelope][modify]") {
	// MODIFY needs the accession to identify the object — without it Webin V2
	// can't disambiguate against a re-used alias. Catch missing accession at
	// build time so the user gets a useful "you forgot the accession" message
	// instead of a server-side "alias not found" half a round-trip later.
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::MODIFY;
	miint::ProjectSpec p;
	p.alias = "needs-accession";
	p.title = "t";
	p.project_type = "METAGENOMIC";
	env.projects.push_back(p);

	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("MODIFY requires accession") &&
	                                                     Catch::Matchers::ContainsSubstring("needs-accession"));
}

TEST_CASE("ENA envelope XML: MODIFY sample emits the existing accession on the sample element",
          "[ena_envelope][modify]") {
	// Production sample envelopes flipped to XML in L4b-fix to bypass the V2
	// JSON dispatcher's <DESCRIPTION>-before-<SAMPLE_NAME> ordering bug.
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::MODIFY;
	miint::SampleSpec s;
	s.alias = "s-2026";
	s.accession = "ERS9999001";
	s.taxon_id = 408170;
	s.title = "Updated sample title";
	s.checklist = "ERC000015";
	s.attributes = {{"collection date", "2026-05-07"}};
	env.samples.push_back(s);

	auto xml = miint::BuildEnvelopeXML(env);
	CHECK(xml.find("<MODIFY/>") != std::string::npos);
	CHECK(xml.find(R"X(<SAMPLE alias="s-2026" accession="ERS9999001">)X") != std::string::npos);
	CHECK(xml.find("<TAXON_ID>408170</TAXON_ID>") != std::string::npos);
}

TEST_CASE("ENA envelope XML: ADD with accession on a sample is rejected at build time", "[ena_envelope][modify]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::ADD;
	miint::SampleSpec s;
	s.alias = "fresh-sample";
	s.accession = "ERS999"; // bogus pre-fill
	s.taxon_id = 408170;
	env.samples.push_back(s);

	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env),
	                  Catch::Matchers::ContainsSubstring("ADD must not set an accession") &&
	                      Catch::Matchers::ContainsSubstring("fresh-sample"));
}

TEST_CASE("ENA envelope XML: MODIFY without accession on a sample is rejected at build time",
          "[ena_envelope][modify]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::MODIFY;
	miint::SampleSpec s;
	s.alias = "needs-accession-sample";
	s.taxon_id = 408170;
	env.samples.push_back(s);

	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("MODIFY requires accession") &&
	                                                    Catch::Matchers::ContainsSubstring("needs-accession-sample"));
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

TEST_CASE("ENA envelope: sample attribute units emit alongside value", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::SampleSpec s;
	s.alias = "gut-002";
	s.taxon_id = 408170;
	s.checklist = "ERC000015";
	s.attributes.emplace_back("collection date", "2026-04-01");
	s.attributes.emplace_back("geographic location (latitude)", "32.7157");
	s.attributes.emplace_back("geographic location (longitude)", "-117.1611");
	// Sparse units: lat/lon get DD; collection date has no units entry.
	s.attribute_units.emplace_back("geographic location (latitude)", "DD");
	s.attribute_units.emplace_back("geographic location (longitude)", "DD");
	env.samples.push_back(s);

	auto json = miint::BuildEnvelopeJSON(env);
	// Lat/lon attribute objects carry the `unit` field. JSON key is singular
	// (`unit`) even though the SQL column / spec field name is plural
	// (`attribute_units`) — V2 validator rejects `units` plural with HTTP 500.
	CHECK(json.find(R"X({"tag":"geographic location (latitude)","value":"32.7157","unit":"DD"})X") !=
	      std::string::npos);
	CHECK(json.find(R"X({"tag":"geographic location (longitude)","value":"-117.1611","unit":"DD"})X") !=
	      std::string::npos);
	// Date attribute does NOT carry units (no entry in attribute_units, no
	// silent default).
	CHECK(json.find(R"X({"tag":"collection date","value":"2026-04-01"})X") != std::string::npos);
	// Empty units string for an entry omits the field entirely (caller can use
	// this to suppress units even if the tag appears in the units map).
	miint::SampleSpec s2 = s;
	s2.alias = "gut-003";
	s2.attribute_units.clear();
	s2.attribute_units.emplace_back("geographic location (latitude)", "");
	miint::SubmissionSpec env2;
	env2.samples.push_back(s2);
	auto json2 = miint::BuildEnvelopeJSON(env2);
	CHECK(json2.find(R"X("unit":)X") == std::string::npos);
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

// =====================================================================
// Experiments + runs
// =====================================================================
//
// Wire shape mirrors the SRA.experiment.xsd / SRA.run.xsd nesting documented
// in localdocs/ena-research-webin-v2-deep.md §4.4 / §4.5, expressed with the
// camelCase JSON keys used elsewhere in the V2 envelope. Cross-references
// use refname (the parent's alias) — V2 also accepts accession, exposed via
// `accession` on the ref struct.

TEST_CASE("ENA envelope: minimal experiment with paired layout and ILLUMINA platform", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ExperimentSpec e;
	e.alias = "exp-001";
	e.study_ref.refname = "p1";
	e.sample_ref.refname = "s1";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "Illumina HiSeq 4000";
	env.experiments.push_back(e);

	auto json = miint::BuildEnvelopeJSON(env);
	CheckEqual(json, R"X({"submission":{"actions":[{"type":"ADD"}]},)X"
	                 R"X("experiments":[{"alias":"exp-001",)X"
	                 R"X("studyRef":{"refname":"p1"},)X"
	                 R"X("design":{"sampleDescriptor":{"refname":"s1"},)X"
	                 R"X("libraryDescriptor":{"libraryStrategy":"WGS","librarySource":"METAGENOMIC",)X"
	                 R"X("librarySelection":"RANDOM","libraryLayout":{"paired":{}}}},)X"
	                 R"X("platform":{"ILLUMINA":{"instrumentModel":"Illumina HiSeq 4000"}}}]})X");
}

TEST_CASE("ENA envelope: experiment with optional fields (title, design description, library_name, single layout)",
          "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ExperimentSpec e;
	e.alias = "exp-002";
	e.title = "Stool 16S amplicon";
	e.study_ref.refname = "p1";
	e.sample_ref.refname = "s1";
	e.design_description = "16S V4 amplicon, 250 bp paired-end";
	e.library_name = "lib-002";
	e.library_strategy = "AMPLICON";
	e.library_source = "METAGENOMIC";
	e.library_selection = "PCR";
	e.library_layout = miint::ENALibraryLayout::SINGLE;
	e.platform = "OXFORD_NANOPORE";
	e.instrument_model = "MinION";
	env.experiments.push_back(e);

	auto json = miint::BuildEnvelopeJSON(env);
	CHECK(json.find(R"X("alias":"exp-002")X") != std::string::npos);
	CHECK(json.find(R"X("title":"Stool 16S amplicon")X") != std::string::npos);
	CHECK(json.find(R"X("designDescription":"16S V4 amplicon, 250 bp paired-end")X") != std::string::npos);
	CHECK(json.find(R"X("libraryName":"lib-002")X") != std::string::npos);
	CHECK(json.find(R"X("libraryLayout":{"single":{}})X") != std::string::npos);
	CHECK(json.find(R"X("OXFORD_NANOPORE":{"instrumentModel":"MinION"})X") != std::string::npos);
}

TEST_CASE("ENA envelope: experiment study_ref accepts accession or refname", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ExperimentSpec e;
	e.alias = "exp-acc";
	e.study_ref.accession = "PRJEB42";
	e.sample_ref.accession = "SAMEA42";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "NovaSeq 6000";
	env.experiments.push_back(e);

	auto json = miint::BuildEnvelopeJSON(env);
	// `accession` wins when set; refname not emitted.
	CHECK(json.find(R"X("studyRef":{"accession":"PRJEB42"})X") != std::string::npos);
	CHECK(json.find(R"X("sampleDescriptor":{"accession":"SAMEA42"})X") != std::string::npos);
}

TEST_CASE("ENA envelope: experiment empty alias rejected", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ExperimentSpec e;
	e.study_ref.refname = "p1";
	e.sample_ref.refname = "s1";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "NovaSeq 6000";
	env.experiments.push_back(e);
	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("alias"));
}

TEST_CASE("ENA envelope: experiment requires study_ref", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ExperimentSpec e;
	e.alias = "exp-no-study";
	e.sample_ref.refname = "s1";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "NovaSeq 6000";
	env.experiments.push_back(e);
	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("study_ref"));
}

TEST_CASE("ENA envelope: experiment requires sample_ref", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ExperimentSpec e;
	e.alias = "exp-no-sample";
	e.study_ref.refname = "p1";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "NovaSeq 6000";
	env.experiments.push_back(e);
	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("sample_ref"));
}

TEST_CASE("ENA envelope: experiment unknown library_strategy rejected", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ExperimentSpec e;
	e.alias = "exp-bad";
	e.study_ref.refname = "p1";
	e.sample_ref.refname = "s1";
	e.library_strategy = "DEFINITELY_NOT_A_REAL_STRATEGY";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "NovaSeq 6000";
	env.experiments.push_back(e);
	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("library_strategy"));
}

TEST_CASE("ENA envelope: experiment unknown platform rejected", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ExperimentSpec e;
	e.alias = "exp-bad-platform";
	e.study_ref.refname = "p1";
	e.sample_ref.refname = "s1";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "MADE_UP_PLATFORM";
	e.instrument_model = "NovaSeq 6000";
	env.experiments.push_back(e);
	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("platform"));
}

TEST_CASE("ENA envelope: minimal run with two paired-fastq files", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::RunSpec r;
	r.alias = "run-001";
	r.experiment_ref.refname = "exp-001";
	r.files.push_back({"sample-A_1.fastq.gz", "fastq", "9b8932f85caa54e687eba62fca3edce2"});
	r.files.push_back({"sample-A_2.fastq.gz", "fastq", "183d6a24e0c3704e993bebe75bbbd989"});
	env.runs.push_back(r);

	auto json = miint::BuildEnvelopeJSON(env);
	CheckEqual(json, R"X({"submission":{"actions":[{"type":"ADD"}]},)X"
	                 R"X("runs":[{"alias":"run-001","experimentRef":{"refname":"exp-001"},)X"
	                 R"X("files":[)X"
	                 R"X({"filename":"sample-A_1.fastq.gz","filetype":"fastq","checksumMethod":"MD5",)X"
	                 R"X("checksum":"9b8932f85caa54e687eba62fca3edce2"},)X"
	                 R"X({"filename":"sample-A_2.fastq.gz","filetype":"fastq","checksumMethod":"MD5",)X"
	                 R"X("checksum":"183d6a24e0c3704e993bebe75bbbd989"}]}]})X");
}

TEST_CASE("ENA envelope: single-end run emits one file", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::RunSpec r;
	r.alias = "run-single";
	r.experiment_ref.refname = "exp-single";
	r.files.push_back({"sample-B.fastq.gz", "fastq", "abcdef0123456789abcdef0123456789"});
	env.runs.push_back(r);
	auto json = miint::BuildEnvelopeJSON(env);
	CHECK(json.find(R"X("files":[{"filename":"sample-B.fastq.gz")X") != std::string::npos);
	CHECK(std::count(json.begin(), json.end(), '{') > 0);
	// Only one FILE entry — search for "filename" once.
	auto first = json.find("\"filename\":");
	REQUIRE(first != std::string::npos);
	auto second = json.find("\"filename\":", first + 1);
	CHECK(second == std::string::npos);
}

TEST_CASE("ENA envelope: run requires alias", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::RunSpec r;
	r.experiment_ref.refname = "exp-001";
	r.files.push_back({"a.fastq.gz", "fastq", "deadbeef"});
	env.runs.push_back(r);
	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("alias"));
}

TEST_CASE("ENA envelope: run requires experiment_ref", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::RunSpec r;
	r.alias = "run-no-exp";
	r.files.push_back({"a.fastq.gz", "fastq", "deadbeef"});
	env.runs.push_back(r);
	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("experiment_ref"));
}

TEST_CASE("ENA envelope: run requires at least one file", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::RunSpec r;
	r.alias = "run-no-files";
	r.experiment_ref.refname = "exp-001";
	env.runs.push_back(r);
	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("file"));
}

TEST_CASE("ENA envelope: run rejects empty checksum", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::RunSpec r;
	r.alias = "run-bad-checksum";
	r.experiment_ref.refname = "exp-001";
	r.files.push_back({"a.fastq.gz", "fastq", ""});
	env.runs.push_back(r);
	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("checksum"));
}

TEST_CASE("ENA envelope: empty filetype rejected (no silent default)", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::RunSpec r;
	r.alias = "run-empty-filetype";
	r.experiment_ref.refname = "exp-001";
	r.files.push_back({"a.fastq.gz", "", "abcd"});
	env.runs.push_back(r);
	CHECK_THROWS_WITH(miint::BuildEnvelopeJSON(env), Catch::Matchers::ContainsSubstring("filetype"));
}

// =====================================================================
// XML envelope (experiments + runs) — V2 server requires XML for SRA-side
// objects; verified live 2026-05-04. See BuildEnvelopeXML in
// src/ena_envelope_builder.cpp.
// =====================================================================

// =====================================================================
// ResolveENARefDescriptor + per-kind wrappers — accession-vs-refname
// disambiguation shared by the INSERT and MODIFY paths.
// =====================================================================

TEST_CASE("ENA RefDescriptor: accession prefix + digits routes to accession", "[ena_envelope]") {
	auto study = miint::ResolveENAStudyRef("PRJEB123456");
	CHECK(study.accession == "PRJEB123456");
	CHECK(study.refname.empty());

	auto sample = miint::ResolveENASampleRef("SAMEA42");
	CHECK(sample.accession == "SAMEA42");
	CHECK(sample.refname.empty());
}

TEST_CASE("ENA RefDescriptor: exact prefix with no digits routes to refname", "[ena_envelope]") {
	// A bare prefix ("PRJEB" with nothing after) is not a valid accession — at
	// least one digit must follow. The `<= prefix.size()` guard in
	// ResolveENARefDescriptor pins this; without it a user alias literally
	// equal to a prefix would be silently misclassified.
	auto study = miint::ResolveENAStudyRef("PRJEB");
	CHECK(study.accession.empty());
	CHECK(study.refname == "PRJEB");

	auto sample = miint::ResolveENASampleRef("ERS");
	CHECK(sample.accession.empty());
	CHECK(sample.refname == "ERS");
}

TEST_CASE("ENA RefDescriptor: prefix followed by non-digits routes to refname", "[ena_envelope]") {
	// A user alias like "ERPmycoolstudy" or "PRJEB-2026-cohort" must NOT be
	// silently classified as an accession the server can't find. The
	// "digits only after prefix" check is what protects against this.
	auto a = miint::ResolveENAStudyRef("ERPmycoolstudy");
	CHECK(a.accession.empty());
	CHECK(a.refname == "ERPmycoolstudy");

	auto b = miint::ResolveENAStudyRef("PRJEB-2026-cohort");
	CHECK(b.accession.empty());
	CHECK(b.refname == "PRJEB-2026-cohort");

	auto c = miint::ResolveENASampleRef("SAMEA_with_underscore");
	CHECK(c.accession.empty());
	CHECK(c.refname == "SAMEA_with_underscore");
}

TEST_CASE("ENA RefDescriptor: non-matching prefix routes to refname", "[ena_envelope]") {
	auto a = miint::ResolveENAStudyRef("my-study-alias");
	CHECK(a.accession.empty());
	CHECK(a.refname == "my-study-alias");

	auto b = miint::ResolveENASampleRef("alpha");
	CHECK(b.accession.empty());
	CHECK(b.refname == "alpha");

	// Empty string also routes to refname (caller's responsibility to reject
	// empties before calling — the helper itself doesn't validate non-emptiness).
	auto c = miint::ResolveENAStudyRef("");
	CHECK(c.accession.empty());
	CHECK(c.refname.empty());
}

TEST_CASE("ENA RefDescriptor: study and sample prefix lists are disjoint", "[ena_envelope]") {
	auto a = miint::ResolveENAStudyRef("SAMEA42");
	CHECK(a.accession.empty());
	CHECK(a.refname == "SAMEA42");

	auto b = miint::ResolveENASampleRef("PRJEB123456");
	CHECK(b.accession.empty());
	CHECK(b.refname == "PRJEB123456");
}

TEST_CASE("ENA RefDescriptor: every canonical study prefix is recognised", "[ena_envelope]") {
	for (const char *p : {"PRJEB", "PRJNA", "PRJDB", "ERP"}) {
		auto ref = miint::ResolveENAStudyRef(std::string(p) + "1");
		CHECK(ref.accession == std::string(p) + "1");
		CHECK(ref.refname.empty());
	}
}

TEST_CASE("ENA RefDescriptor: every canonical sample prefix is recognised", "[ena_envelope]") {
	for (const char *p : {"ERS", "SAMEA", "SAMN", "SAMD"}) {
		auto ref = miint::ResolveENASampleRef(std::string(p) + "1");
		CHECK(ref.accession == std::string(p) + "1");
		CHECK(ref.refname.empty());
	}
}

TEST_CASE("ENA RefDescriptor: generic primitive accepts arbitrary prefix lists", "[ena_envelope]") {
	// Named wrappers cover the L4c-shipped kinds; the generic primitive stays
	// exposed so future kinds (analyses, runs) can supply their own prefix
	// lists without creating a wrapper that's later abandoned.
	auto a = miint::ResolveENARefDescriptor("ERX42", {"ERX"});
	CHECK(a.accession == "ERX42");

	auto b = miint::ResolveENARefDescriptor("not-an-accession", {"ERX"});
	CHECK(b.accession.empty());
	CHECK(b.refname == "not-an-accession");
}

TEST_CASE("ENA envelope XML: minimal experiment with paired layout and ILLUMINA platform", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ExperimentSpec e;
	e.alias = "exp-001";
	e.study_ref.accession = "PRJEB42";
	e.sample_ref.accession = "SAMEA42";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "Illumina NovaSeq 6000";
	env.experiments.push_back(e);

	auto xml = miint::BuildEnvelopeXML(env);
	CheckEqual(
	    xml, R"X(<?xml version="1.0" encoding="UTF-8"?>)X"
	         R"X(<WEBIN><SUBMISSION><ACTIONS><ACTION><ADD/></ACTION></ACTIONS></SUBMISSION>)X"
	         R"X(<EXPERIMENT_SET>)X"
	         R"X(<EXPERIMENT alias="exp-001">)X"
	         R"X(<STUDY_REF accession="PRJEB42"/>)X"
	         R"X(<DESIGN><DESIGN_DESCRIPTION></DESIGN_DESCRIPTION>)X"
	         R"X(<SAMPLE_DESCRIPTOR accession="SAMEA42"/>)X"
	         R"X(<LIBRARY_DESCRIPTOR><LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>)X"
	         R"X(<LIBRARY_SOURCE>METAGENOMIC</LIBRARY_SOURCE>)X"
	         R"X(<LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION>)X"
	         R"X(<LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT>)X"
	         R"X(</LIBRARY_DESCRIPTOR></DESIGN>)X"
	         R"X(<PLATFORM><ILLUMINA><INSTRUMENT_MODEL>Illumina NovaSeq 6000</INSTRUMENT_MODEL></ILLUMINA></PLATFORM>)X"
	         R"X(</EXPERIMENT></EXPERIMENT_SET></WEBIN>)X");
}

TEST_CASE("ENA envelope XML: experiment with refname-style cross-references and optional fields", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ExperimentSpec e;
	e.alias = "exp-002";
	e.title = "Stool 16S amplicon";
	e.study_ref.refname = "p1";
	e.sample_ref.refname = "s1";
	e.design_description = "16S V4 amplicon, 250 bp paired-end";
	e.library_name = "lib-002";
	e.library_strategy = "AMPLICON";
	e.library_source = "METAGENOMIC";
	e.library_selection = "PCR";
	e.library_layout = miint::ENALibraryLayout::SINGLE;
	e.platform = "OXFORD_NANOPORE";
	e.instrument_model = "MinION";
	env.experiments.push_back(e);

	auto xml = miint::BuildEnvelopeXML(env);
	CHECK(xml.find(R"X(<EXPERIMENT alias="exp-002">)X") != std::string::npos);
	CHECK(xml.find(R"X(<TITLE>Stool 16S amplicon</TITLE>)X") != std::string::npos);
	CHECK(xml.find(R"X(<STUDY_REF refname="p1"/>)X") != std::string::npos);
	CHECK(xml.find(R"X(<SAMPLE_DESCRIPTOR refname="s1"/>)X") != std::string::npos);
	CHECK(xml.find(R"X(<DESIGN_DESCRIPTION>16S V4 amplicon, 250 bp paired-end</DESIGN_DESCRIPTION>)X") !=
	      std::string::npos);
	CHECK(xml.find(R"X(<LIBRARY_NAME>lib-002</LIBRARY_NAME>)X") != std::string::npos);
	CHECK(xml.find(R"X(<LIBRARY_LAYOUT><SINGLE/></LIBRARY_LAYOUT>)X") != std::string::npos);
	CHECK(xml.find(R"X(<OXFORD_NANOPORE><INSTRUMENT_MODEL>MinION</INSTRUMENT_MODEL></OXFORD_NANOPORE>)X") !=
	      std::string::npos);
}

TEST_CASE("ENA envelope XML: minimal run with paired files", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::RunSpec r;
	r.alias = "run-001";
	r.experiment_ref.refname = "exp-001";
	r.files.push_back({"sample-A_1.fastq.gz", "fastq", "9b8932f85caa54e687eba62fca3edce2"});
	r.files.push_back({"sample-A_2.fastq.gz", "fastq", "183d6a24e0c3704e993bebe75bbbd989"});
	env.runs.push_back(r);

	auto xml = miint::BuildEnvelopeXML(env);
	CheckEqual(xml, R"X(<?xml version="1.0" encoding="UTF-8"?>)X"
	                R"X(<WEBIN><SUBMISSION><ACTIONS><ACTION><ADD/></ACTION></ACTIONS></SUBMISSION>)X"
	                R"X(<RUN_SET><RUN alias="run-001">)X"
	                R"X(<EXPERIMENT_REF refname="exp-001"/>)X"
	                R"X(<DATA_BLOCK><FILES>)X"
	                R"X(<FILE filename="sample-A_1.fastq.gz" filetype="fastq" checksum_method="MD5" )X"
	                R"X(checksum="9b8932f85caa54e687eba62fca3edce2"/>)X"
	                R"X(<FILE filename="sample-A_2.fastq.gz" filetype="fastq" checksum_method="MD5" )X"
	                R"X(checksum="183d6a24e0c3704e993bebe75bbbd989"/>)X"
	                R"X(</FILES></DATA_BLOCK></RUN></RUN_SET></WEBIN>)X");
}

TEST_CASE("ENA envelope XML: combined experiments + runs in one envelope", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ExperimentSpec e;
	e.alias = "e1";
	e.study_ref.refname = "p1";
	e.sample_ref.refname = "s1";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "Illumina NovaSeq 6000";
	env.experiments.push_back(e);
	miint::RunSpec r;
	r.alias = "r1";
	r.experiment_ref.refname = "e1";
	r.files.push_back({"a.fastq.gz", "fastq", "deadbeef"});
	env.runs.push_back(r);

	auto xml = miint::BuildEnvelopeXML(env);
	auto submission_pos = xml.find("<SUBMISSION>");
	auto exp_set_pos = xml.find("<EXPERIMENT_SET>");
	auto run_set_pos = xml.find("<RUN_SET>");
	REQUIRE(submission_pos != std::string::npos);
	REQUIRE(exp_set_pos != std::string::npos);
	REQUIRE(run_set_pos != std::string::npos);
	// Spec order: SUBMISSION, EXPERIMENT_SET, RUN_SET.
	CHECK(submission_pos < exp_set_pos);
	CHECK(exp_set_pos < run_set_pos);
}

TEST_CASE("ENA envelope XML: HOLD action emits HoldUntilDate sibling action", "[ena_envelope]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::ADD;
	env.hold_until_date = "2027-01-15";
	miint::ExperimentSpec e;
	e.alias = "e-hold";
	e.study_ref.refname = "p1";
	e.sample_ref.refname = "s1";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "Illumina NovaSeq 6000";
	env.experiments.push_back(e);
	auto xml = miint::BuildEnvelopeXML(env);
	CHECK(xml.find(R"X(<ACTION><HOLD HoldUntilDate="2027-01-15"/></ACTION>)X") != std::string::npos);
}

TEST_CASE("ENA envelope XML: XML escaping of attribute and text values", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::RunSpec r;
	r.alias = "run<&\"'>";
	r.experiment_ref.refname = "exp&1";
	r.files.push_back({"a&b.fastq.gz", "fastq", "deadbeef"});
	env.runs.push_back(r);
	auto xml = miint::BuildEnvelopeXML(env);
	// Alias attribute value is escaped (used inside double-quoted attr).
	CHECK(xml.find(R"X(<RUN alias="run&lt;&amp;&quot;&apos;&gt;">)X") != std::string::npos);
	CHECK(xml.find(R"X(<EXPERIMENT_REF refname="exp&amp;1"/>)X") != std::string::npos);
	CHECK(xml.find(R"X(filename="a&amp;b.fastq.gz")X") != std::string::npos);
}

TEST_CASE("ENA envelope XML: empty experiments + runs yields just the SUBMISSION block", "[ena_envelope]") {
	miint::SubmissionSpec env;
	auto xml = miint::BuildEnvelopeXML(env);
	CheckEqual(xml, R"X(<?xml version="1.0" encoding="UTF-8"?>)X"
	                R"X(<WEBIN><SUBMISSION><ACTIONS><ACTION><ADD/></ACTION></ACTIONS></SUBMISSION></WEBIN>)X");
}

TEST_CASE("ENA envelope XML: MODIFY experiment emits the existing accession on the experiment element",
          "[ena_envelope][modify]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::MODIFY;
	miint::ExperimentSpec e;
	e.alias = "exp-2026";
	e.accession = "ERX9999001";
	e.title = "Updated experiment title";
	e.study_ref.accession = "PRJEB42";
	e.sample_ref.accession = "SAMEA42";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "Illumina NovaSeq 6000";
	env.experiments.push_back(e);

	auto xml = miint::BuildEnvelopeXML(env);
	CHECK(xml.find("<MODIFY/>") != std::string::npos);
	CHECK(xml.find(R"X(<EXPERIMENT alias="exp-2026" accession="ERX9999001">)X") != std::string::npos);
	CHECK(xml.find("<TITLE>Updated experiment title</TITLE>") != std::string::npos);
}

TEST_CASE("ENA envelope XML: ADD with accession on an experiment is rejected at build time", "[ena_envelope][modify]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::ADD;
	miint::ExperimentSpec e;
	e.alias = "fresh-experiment";
	e.accession = "ERX999"; // bogus pre-fill
	e.study_ref.refname = "p1";
	e.sample_ref.refname = "s1";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "Illumina NovaSeq 6000";
	env.experiments.push_back(e);

	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env),
	                  Catch::Matchers::ContainsSubstring("ADD must not set an accession") &&
	                      Catch::Matchers::ContainsSubstring("fresh-experiment"));
}

TEST_CASE("ENA envelope XML: MODIFY without accession on an experiment is rejected at build time",
          "[ena_envelope][modify]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::MODIFY;
	miint::ExperimentSpec e;
	e.alias = "needs-accession-experiment";
	// e.accession deliberately empty
	e.study_ref.refname = "p1";
	e.sample_ref.refname = "s1";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "Illumina NovaSeq 6000";
	env.experiments.push_back(e);

	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env),
	                  Catch::Matchers::ContainsSubstring("MODIFY requires accession") &&
	                      Catch::Matchers::ContainsSubstring("needs-accession-experiment"));
}

TEST_CASE("ENA envelope XML: experiment empty alias rejected", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ExperimentSpec e;
	e.study_ref.refname = "p1";
	e.sample_ref.refname = "s1";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "Illumina NovaSeq 6000";
	env.experiments.push_back(e);
	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("alias"));
}

TEST_CASE("ENA envelope XML: run with no files rejected", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::RunSpec r;
	r.alias = "run-no-files";
	r.experiment_ref.refname = "e1";
	env.runs.push_back(r);
	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("file"));
}

// =====================================================================
// XML envelope: samples — flipped from JSON to XML in L4b-fix to bypass the
// V2 JSON dispatcher's element-ordering bug where <DESCRIPTION> is emitted
// before <SAMPLE_NAME> regardless of JSON-key order, violating
// SRA.sample.xsd. Verified live on wwwdev 2026-05-07.
// =====================================================================

TEST_CASE("ENA envelope XML: minimal sample with checklist and attributes", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::SampleSpec s;
	s.alias = "gut-001";
	s.title = "Adult subject 001 stool";
	s.taxon_id = 408170;
	s.checklist = "ERC000015";
	s.attributes.emplace_back("collection date", "2024-06-15");
	s.attributes.emplace_back("geographic location (country and/or sea)", "United States");
	env.samples.push_back(s);

	auto xml = miint::BuildEnvelopeXML(env);
	CheckEqual(xml, R"X(<?xml version="1.0" encoding="UTF-8"?>)X"
	                R"X(<WEBIN><SUBMISSION><ACTIONS><ACTION><ADD/></ACTION></ACTIONS></SUBMISSION>)X"
	                R"X(<SAMPLE_SET>)X"
	                R"X(<SAMPLE alias="gut-001">)X"
	                R"X(<TITLE>Adult subject 001 stool</TITLE>)X"
	                R"X(<SAMPLE_NAME><TAXON_ID>408170</TAXON_ID></SAMPLE_NAME>)X"
	                R"X(<SAMPLE_ATTRIBUTES>)X"
	                R"X(<SAMPLE_ATTRIBUTE><TAG>ENA-CHECKLIST</TAG><VALUE>ERC000015</VALUE></SAMPLE_ATTRIBUTE>)X"
	                R"X(<SAMPLE_ATTRIBUTE><TAG>collection date</TAG><VALUE>2024-06-15</VALUE></SAMPLE_ATTRIBUTE>)X"
	                R"X(<SAMPLE_ATTRIBUTE><TAG>geographic location (country and/or sea)</TAG>)X"
	                R"X(<VALUE>United States</VALUE></SAMPLE_ATTRIBUTE>)X"
	                R"X(</SAMPLE_ATTRIBUTES>)X"
	                R"X(</SAMPLE></SAMPLE_SET></WEBIN>)X");
}

TEST_CASE("ENA envelope XML: minimal MODIFY shape with no optional fields set", "[ena_envelope][modify]") {
	// The smallest legal MODIFY a caller could send: alias + accession +
	// taxon_id, no title, no description, no scientific_name, no checklist,
	// no attributes. Confirms the optional-field gating in AppendXmlSample
	// produces a clean envelope rather than empty <TITLE/>, <DESCRIPTION/>,
	// or <SAMPLE_ATTRIBUTES/> elements (each of which would be rejected by
	// SRA.sample.xsd's content-model on those children).
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::MODIFY;
	miint::SampleSpec s;
	s.alias = "minimal-mod";
	s.accession = "ERS9999777";
	s.taxon_id = 408170;
	env.samples.push_back(s);

	auto xml = miint::BuildEnvelopeXML(env);
	CheckEqual(xml, R"X(<?xml version="1.0" encoding="UTF-8"?>)X"
	                R"X(<WEBIN><SUBMISSION><ACTIONS><ACTION><MODIFY/></ACTION></ACTIONS></SUBMISSION>)X"
	                R"X(<SAMPLE_SET>)X"
	                R"X(<SAMPLE alias="minimal-mod" accession="ERS9999777">)X"
	                R"X(<SAMPLE_NAME><TAXON_ID>408170</TAXON_ID></SAMPLE_NAME>)X"
	                R"X(</SAMPLE></SAMPLE_SET></WEBIN>)X");
}

TEST_CASE("ENA envelope XML: SAMPLE_NAME precedes DESCRIPTION (SRA.sample.xsd ordering)", "[ena_envelope]") {
	// The bug that drove this phase: V2's JSON dispatcher emitted <DESCRIPTION>
	// BEFORE <SAMPLE_NAME>, violating SRA.sample.xsd. Going through XML gives
	// us direct control over element order — pin that here so a future
	// reordering of AppendXmlSample doesn't silently regress the XSD ordering.
	miint::SubmissionSpec env;
	miint::SampleSpec s;
	s.alias = "with-description";
	s.taxon_id = 408170;
	s.scientific_name = "human gut metagenome";
	s.description = "clinical isolate from a healthy donor";
	env.samples.push_back(s);

	auto xml = miint::BuildEnvelopeXML(env);
	const auto sample_name_pos = xml.find("<SAMPLE_NAME>");
	const auto description_pos = xml.find("<DESCRIPTION>");
	REQUIRE(sample_name_pos != std::string::npos);
	REQUIRE(description_pos != std::string::npos);
	CHECK(sample_name_pos < description_pos);
	CHECK(xml.find("<SCIENTIFIC_NAME>human gut metagenome</SCIENTIFIC_NAME>") != std::string::npos);
	CHECK(xml.find("<DESCRIPTION>clinical isolate from a healthy donor</DESCRIPTION>") != std::string::npos);
}

TEST_CASE("ENA envelope XML: sample attribute units emit a <UNITS> child sibling", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::SampleSpec s;
	s.alias = "gut-002";
	s.taxon_id = 408170;
	s.checklist = "ERC000015";
	s.attributes.emplace_back("collection date", "2026-04-01");
	s.attributes.emplace_back("geographic location (latitude)", "32.7157");
	s.attributes.emplace_back("geographic location (longitude)", "-117.1611");
	s.attribute_units.emplace_back("geographic location (latitude)", "DD");
	s.attribute_units.emplace_back("geographic location (longitude)", "DD");
	env.samples.push_back(s);

	auto xml = miint::BuildEnvelopeXML(env);
	CHECK(xml.find(R"X(<SAMPLE_ATTRIBUTE><TAG>geographic location (latitude)</TAG>)X"
	               R"X(<VALUE>32.7157</VALUE><UNITS>DD</UNITS></SAMPLE_ATTRIBUTE>)X") != std::string::npos);
	CHECK(xml.find(R"X(<SAMPLE_ATTRIBUTE><TAG>geographic location (longitude)</TAG>)X"
	               R"X(<VALUE>-117.1611</VALUE><UNITS>DD</UNITS></SAMPLE_ATTRIBUTE>)X") != std::string::npos);
	CHECK(xml.find(R"X(<SAMPLE_ATTRIBUTE><TAG>collection date</TAG>)X"
	               R"X(<VALUE>2026-04-01</VALUE></SAMPLE_ATTRIBUTE>)X") != std::string::npos);

	// Empty units string for an entry suppresses the <UNITS> sibling entirely
	// (caller can use this to opt out without removing the tag from
	// attribute_units).
	miint::SampleSpec s2 = s;
	s2.alias = "gut-003";
	s2.attribute_units.clear();
	s2.attribute_units.emplace_back("geographic location (latitude)", "");
	miint::SubmissionSpec env2;
	env2.samples.push_back(s2);
	auto xml2 = miint::BuildEnvelopeXML(env2);
	CHECK(xml2.find("<UNITS>") == std::string::npos);
}

TEST_CASE("ENA envelope XML: multi-sample emits one element per sample", "[ena_envelope]") {
	miint::SubmissionSpec env;
	for (int i = 1; i <= 3; ++i) {
		miint::SampleSpec s;
		s.alias = "sample-" + std::to_string(i);
		s.taxon_id = 408170;
		s.checklist = "ERC000015";
		env.samples.push_back(s);
	}
	auto xml = miint::BuildEnvelopeXML(env);
	CHECK(xml.find(R"X(<SAMPLE alias="sample-1">)X") != std::string::npos);
	CHECK(xml.find(R"X(<SAMPLE alias="sample-2">)X") != std::string::npos);
	CHECK(xml.find(R"X(<SAMPLE alias="sample-3">)X") != std::string::npos);
}

TEST_CASE("ENA envelope XML: SAMPLE_SET precedes EXPERIMENT_SET and RUN_SET", "[ena_envelope]") {
	// Spec order per SRA.submission.xsd: SAMPLE_SET → EXPERIMENT_SET → RUN_SET.
	miint::SubmissionSpec env;
	miint::SampleSpec s;
	s.alias = "s1";
	s.taxon_id = 408170;
	env.samples.push_back(s);
	miint::ExperimentSpec e;
	e.alias = "e1";
	e.study_ref.refname = "p1";
	e.sample_ref.refname = "s1";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "NovaSeq 6000";
	env.experiments.push_back(e);
	miint::RunSpec r;
	r.alias = "r1";
	r.experiment_ref.refname = "e1";
	r.files.push_back({"r1_1.fastq.gz", "fastq", "deadbeef"});
	env.runs.push_back(r);

	auto xml = miint::BuildEnvelopeXML(env);
	const auto submission_pos = xml.find("<SUBMISSION>");
	const auto sample_pos = xml.find("<SAMPLE_SET>");
	const auto exp_pos = xml.find("<EXPERIMENT_SET>");
	const auto run_pos = xml.find("<RUN_SET>");
	REQUIRE(submission_pos != std::string::npos);
	REQUIRE(sample_pos != std::string::npos);
	REQUIRE(exp_pos != std::string::npos);
	REQUIRE(run_pos != std::string::npos);
	CHECK(submission_pos < sample_pos);
	CHECK(sample_pos < exp_pos);
	CHECK(exp_pos < run_pos);
}

TEST_CASE("ENA envelope XML: sample with taxon_id <= 0 rejected", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::SampleSpec s;
	s.alias = "s1";
	s.taxon_id = 0;
	env.samples.push_back(s);
	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("taxon"));
}

TEST_CASE("ENA envelope XML: sample empty alias rejected", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::SampleSpec s;
	s.taxon_id = 408170;
	env.samples.push_back(s);
	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("alias"));
}

TEST_CASE("ENA envelope XML: sample XML escaping of alias, attribute tag/value, and text", "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::SampleSpec s;
	s.alias = R"X(s<&"'>1)X";
	s.taxon_id = 408170;
	s.title = R"X(<unsafe & "title")X";
	s.attributes.emplace_back(R"X(tag<&)X", R"X(value>")X");
	env.samples.push_back(s);

	auto xml = miint::BuildEnvelopeXML(env);
	CHECK(xml.find(R"X(<SAMPLE alias="s&lt;&amp;&quot;&apos;&gt;1">)X") != std::string::npos);
	CHECK(xml.find(R"X(<TITLE>&lt;unsafe &amp; &quot;title&quot;</TITLE>)X") != std::string::npos);
	CHECK(xml.find(R"X(<TAG>tag&lt;&amp;</TAG><VALUE>value&gt;&quot;</VALUE>)X") != std::string::npos);
}

TEST_CASE("ENA envelope: full graph (project + sample + experiment + run) emits arrays in spec order",
          "[ena_envelope]") {
	miint::SubmissionSpec env;
	miint::ProjectSpec p;
	p.alias = "p1";
	p.title = "P1";
	env.projects.push_back(p);
	miint::SampleSpec s;
	s.alias = "s1";
	s.taxon_id = 408170;
	env.samples.push_back(s);
	miint::ExperimentSpec e;
	e.alias = "e1";
	e.study_ref.refname = "p1";
	e.sample_ref.refname = "s1";
	e.library_strategy = "WGS";
	e.library_source = "METAGENOMIC";
	e.library_selection = "RANDOM";
	e.library_layout = miint::ENALibraryLayout::PAIRED;
	e.platform = "ILLUMINA";
	e.instrument_model = "NovaSeq 6000";
	env.experiments.push_back(e);
	miint::RunSpec r;
	r.alias = "r1";
	r.experiment_ref.refname = "e1";
	r.files.push_back({"r1_1.fastq.gz", "fastq", "deadbeef"});
	r.files.push_back({"r1_2.fastq.gz", "fastq", "cafebabe"});
	env.runs.push_back(r);

	auto json = miint::BuildEnvelopeJSON(env);
	const auto sub_pos = json.find("\"submission\":");
	const auto proj_pos = json.find("\"projects\":");
	const auto samp_pos = json.find("\"samples\":");
	const auto exp_pos = json.find("\"experiments\":");
	const auto run_pos = json.find("\"runs\":");
	REQUIRE(sub_pos != std::string::npos);
	REQUIRE(proj_pos != std::string::npos);
	REQUIRE(samp_pos != std::string::npos);
	REQUIRE(exp_pos != std::string::npos);
	REQUIRE(run_pos != std::string::npos);
	CHECK(sub_pos < proj_pos);
	CHECK(proj_pos < samp_pos);
	CHECK(samp_pos < exp_pos);
	CHECK(exp_pos < run_pos);
}

// =====================================================================
// Lifecycle (targeted) actions: CANCEL / HOLD / RELEASE on existing accessions
// =====================================================================
//
// These actions reference an already-registered object via `target=` on the
// action element. Body sets (PROJECT_SET / SAMPLE_SET / EXPERIMENT_SET /
// RUN_SET) are not emitted — the action itself is the entire payload.
// MODIFY also uses `action=MODIFY` but with a body identifying the object
// (covered separately under the existing per-object envelope tests).

TEST_CASE("ENA envelope: CANCEL targets an accession (XML)", "[ena_envelope][lifecycle]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::CANCEL;
	env.target_accession = "ERS123456";

	auto xml = miint::BuildEnvelopeXML(env);
	CheckEqual(xml, R"X(<?xml version="1.0" encoding="UTF-8"?>)X"
	                R"X(<WEBIN><SUBMISSION>)X"
	                R"X(<ACTIONS><ACTION><CANCEL target="ERS123456"/></ACTION></ACTIONS>)X"
	                R"X(</SUBMISSION></WEBIN>)X");
}

TEST_CASE("ENA envelope: RELEASE targets an accession (XML)", "[ena_envelope][lifecycle]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::RELEASE;
	env.target_accession = "ERS123456";

	auto xml = miint::BuildEnvelopeXML(env);
	CheckEqual(xml, R"X(<?xml version="1.0" encoding="UTF-8"?>)X"
	                R"X(<WEBIN><SUBMISSION>)X"
	                R"X(<ACTIONS><ACTION><RELEASE target="ERS123456"/></ACTION></ACTIONS>)X"
	                R"X(</SUBMISSION></WEBIN>)X");
}

TEST_CASE("ENA envelope: HOLD targets an accession with HoldUntilDate (XML)", "[ena_envelope][lifecycle]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::HOLD;
	env.target_accession = "ERS123456";
	env.hold_until_date = "2027-01-01";

	auto xml = miint::BuildEnvelopeXML(env);
	CheckEqual(xml, R"X(<?xml version="1.0" encoding="UTF-8"?>)X"
	                R"X(<WEBIN><SUBMISSION>)X"
	                R"X(<ACTIONS><ACTION><HOLD target="ERS123456" HoldUntilDate="2027-01-01"/></ACTION></ACTIONS>)X"
	                R"X(</SUBMISSION></WEBIN>)X");
}

TEST_CASE("ENA envelope: target_refname on a lifecycle action is rejected", "[ena_envelope][lifecycle]") {
	// Webin V2 only resolves refname/alias on `<*_REF>` children inside an
	// ADD body where the alias is also defined locally; for cross-submission
	// lifecycle ops on already-registered objects it requires the server-
	// assigned accession. Verified live 2026-05-07.
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::CANCEL;
	env.target_refname = "my-sample-alias";

	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env),
	                  Catch::Matchers::ContainsSubstring("by refname/alias is not supported"));
}

TEST_CASE("ENA envelope: CANCEL without target is rejected", "[ena_envelope][lifecycle]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::CANCEL;
	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("CANCEL"));
}

TEST_CASE("ENA envelope: RELEASE without target is rejected", "[ena_envelope][lifecycle]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::RELEASE;
	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("RELEASE"));
}

TEST_CASE("ENA envelope: HOLD with target requires hold_until_date", "[ena_envelope][lifecycle]") {
	// action=HOLD with a target is the post-hoc embargo path; the date is
	// required (server uses it as the new hold-until-date).
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::HOLD;
	env.target_accession = "ERS123456";
	// hold_until_date deliberately empty
	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("hold_until_date"));
}

TEST_CASE("ENA envelope: ADD with target_accession is rejected (ADD doesn't target)", "[ena_envelope][lifecycle]") {
	// target_accession is exclusive to lifecycle actions. Setting it on an
	// ADD/MODIFY/VALIDATE catches a category-of-mistake at envelope-build time.
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::ADD;
	env.target_accession = "ERS123456";
	miint::ProjectSpec p;
	p.alias = "p1";
	p.title = "t";
	env.projects.push_back(p);
	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("target"));
}

TEST_CASE("ENA envelope: CANCEL with target_refname set is rejected even when accession is also set",
          "[ena_envelope][lifecycle]") {
	// Pre-tighten this verified an "accession wins" precedence for callers
	// that supplied both. After 2026-05-07, target_refname is rejected
	// outright on lifecycle actions because it has no usable wire-form
	// (Webin V2 cannot resolve refname cross-submission). Carrying both
	// is now a programming mistake — fail loudly so the caller drops the
	// refname rather than relying on silent precedence.
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::CANCEL;
	env.target_accession = "ERS123456";
	env.target_refname = "my-alias";

	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env),
	                  Catch::Matchers::ContainsSubstring("by refname/alias is not supported"));
}

TEST_CASE("ENA envelope: CANCEL with body content is rejected (no silent drop)", "[ena_envelope][lifecycle]") {
	// A targeted CANCEL is the entire payload; if a caller accidentally also
	// fills in projects/samples (e.g. reusing a SubmissionSpec across calls),
	// silently dropping the body would look like a successful submit-then-cancel.
	// Reject loudly instead so the caller fixes their spec construction.
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::CANCEL;
	env.target_accession = "ERS123456";
	miint::ProjectSpec p;
	p.alias = "stray";
	p.title = "should-not-appear";
	env.projects.push_back(p);
	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("body content"));
}

TEST_CASE("ENA envelope: HOLD-with-target with body content is rejected", "[ena_envelope][lifecycle]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::HOLD;
	env.target_accession = "ERS123456";
	env.hold_until_date = "2027-01-01";
	miint::SampleSpec s;
	s.alias = "stray";
	s.taxon_id = 408170;
	env.samples.push_back(s);
	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("body content"));
}

TEST_CASE("ENA envelope: MODIFY with target_accession is rejected", "[ena_envelope][lifecycle]") {
	// MODIFY identifies its object via the body, not via target=. A target
	// here would be a programming mistake.
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::MODIFY;
	env.target_accession = "ERS123456";
	miint::SampleSpec s;
	s.alias = "s1";
	s.taxon_id = 408170;
	env.samples.push_back(s);
	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("target"));
}

TEST_CASE("ENA envelope: VALIDATE with target_refname is rejected", "[ena_envelope][lifecycle]") {
	// Two distinct ValidateActions checks could throw here:
	//   (a) L1d refname-rejection: "...VALIDATE by refname/alias is not
	//       supported by Webin V2 for cross-submission lifecycle ops..."
	//   (b) Older "ADD/MODIFY/VALIDATE with target": "...VALIDATE action does
	//       not take a target accession or refname"
	// (a) currently fires first because `HasNonWhitespace(target_refname)` is
	// checked before the action-vs-target compatibility check. Match on the
	// substring "refname" — present in both messages — so this test stays
	// green if the rejection order ever swaps back.
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::VALIDATE;
	env.target_refname = "my-alias";
	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("refname"));
}

TEST_CASE("ENA envelope: whitespace-only target_accession is rejected", "[ena_envelope][lifecycle]") {
	// "   " has non-zero length but is meaningless to the server. Catch it
	// at envelope-build time rather than letting the round-trip discover it.
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::CANCEL;
	env.target_accession = "   ";
	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("target_accession"));
}

TEST_CASE("ENA envelope: whitespace-only target_refname is rejected", "[ena_envelope][lifecycle]") {
	miint::SubmissionSpec env;
	env.action = miint::ENAAction::RELEASE;
	env.target_refname = "\t\n";
	CHECK_THROWS_WITH(miint::BuildEnvelopeXML(env), Catch::Matchers::ContainsSubstring("target_refname"));
}
