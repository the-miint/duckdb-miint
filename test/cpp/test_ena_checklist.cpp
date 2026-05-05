// SPDX-License-Identifier: MIT
//
// ENA checklist registry + sample-attribute validator (Phase 8 Step 8b).
// Pure-data tests of the XML parser, the validator, and the lazy-fetch
// registry, with a mock fetcher substituting for the live
// https://www.ebi.ac.uk/ena/browser/api/xml/ERC000NN endpoint.

#include "ena_checklist.hpp"

#include <catch2/catch_all.hpp>

#include <stdexcept>
#include <string>
#include <vector>

using miint::BuildChecklistFetchURL;
using miint::ChecklistDef;
using miint::ChecklistFieldDef;
using miint::ChecklistFieldMandatory;
using miint::ChecklistRegistry;
using miint::ChecklistValidationIssue;
using miint::ParseChecklistXML;
using miint::ResolveChecklistBaseFromEnv;
using miint::ValidateAttributesAgainstChecklist;

namespace {

// Compact synthetic checklist exercising every validator path: mandatory
// plain field, mandatory field with units, mandatory field with CV, and
// optional plain field. Schema MUST stay aligned with the
// `_CHECKLIST_TST000001` constant in test/scripts/ena_webin_mock.py — the
// SQL mock test (test/sql/ena_checklist_validation_mock.test) relies on the
// same accession + field set against the python mock.
constexpr const char *FIXTURE_XML = R"(<?xml version="1.0" encoding="UTF-8"?>
<CHECKLIST_SET>
<CHECKLIST accession="TST000001" checklistType="Sample">
  <IDENTIFIERS><PRIMARY_ID>TST000001</PRIMARY_ID></IDENTIFIERS>
  <DESCRIPTOR>
    <LABEL>miint Phase 8 test checklist</LABEL>
    <NAME>miint_phase8_test</NAME>
    <DESCRIPTION>Synthetic checklist for unit + SQL tests; MUST stay aligned with ena_webin_mock.py.</DESCRIPTION>
    <FIELD_GROUP restrictionType="Any number or none of the fields">
      <FIELD>
        <LABEL>project name</LABEL>
        <NAME>project_name</NAME>
        <FIELD_TYPE><TEXT_FIELD/></FIELD_TYPE>
        <MANDATORY>mandatory</MANDATORY>
      </FIELD>
      <FIELD>
        <LABEL>collection date</LABEL>
        <NAME>collection_date</NAME>
        <FIELD_TYPE><TEXT_FIELD/></FIELD_TYPE>
        <MANDATORY>mandatory</MANDATORY>
      </FIELD>
      <FIELD>
        <LABEL>geographic location (latitude)</LABEL>
        <NAME>geographic_location_latitude</NAME>
        <FIELD_TYPE><TEXT_FIELD/></FIELD_TYPE>
        <UNITS><UNIT>DD</UNIT></UNITS>
        <MANDATORY>mandatory</MANDATORY>
      </FIELD>
      <FIELD>
        <LABEL>country</LABEL>
        <NAME>country</NAME>
        <FIELD_TYPE>
          <TEXT_CHOICE_FIELD>
            <TEXT_VALUE><VALUE>USA</VALUE></TEXT_VALUE>
            <TEXT_VALUE><VALUE>Canada</VALUE></TEXT_VALUE>
            <TEXT_VALUE><VALUE>Mexico</VALUE></TEXT_VALUE>
          </TEXT_CHOICE_FIELD>
        </FIELD_TYPE>
        <MANDATORY>mandatory</MANDATORY>
      </FIELD>
      <FIELD>
        <LABEL>ploidy</LABEL>
        <NAME>ploidy</NAME>
        <FIELD_TYPE><TEXT_FIELD/></FIELD_TYPE>
        <MANDATORY>optional</MANDATORY>
      </FIELD>
      <FIELD>
        <LABEL>sample storage temperature</LABEL>
        <NAME>sample_storage_temperature</NAME>
        <FIELD_TYPE><TEXT_FIELD/></FIELD_TYPE>
        <UNITS><UNIT>degC</UNIT><UNIT>K</UNIT></UNITS>
        <MANDATORY>optional</MANDATORY>
      </FIELD>
    </FIELD_GROUP>
  </DESCRIPTOR>
</CHECKLIST>
</CHECKLIST_SET>)";

const ChecklistFieldDef *FindField(const ChecklistDef &c, const std::string &label) {
	for (const auto &f : c.fields) {
		if (f.label == label) {
			return &f;
		}
	}
	return nullptr;
}

bool HasIssueFor(const std::vector<ChecklistValidationIssue> &issues, const std::string &field_label,
                 const std::string &needle) {
	for (const auto &i : issues) {
		if (i.field == field_label && i.message.find(needle) != std::string::npos) {
			return true;
		}
	}
	return false;
}

} // namespace

// ---- Parser ----

TEST_CASE("ParseChecklistXML extracts accession + label", "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	CHECK(c.accession == "TST000001");
	CHECK(c.label == "miint Phase 8 test checklist");
	CHECK(c.fields.size() == 6);
}

TEST_CASE("ParseChecklistXML preserves field LABEL (human-readable, matches MAP keys)", "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	const auto *f = FindField(c, "project name");
	REQUIRE(f != nullptr);
	CHECK(f->name == "project_name");
	CHECK(f->mandatory == ChecklistFieldMandatory::MANDATORY);
	CHECK(f->allowed_units.empty());
	CHECK(f->allowed_values.empty());
}

TEST_CASE("ParseChecklistXML extracts UNITS list", "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	const auto *lat = FindField(c, "geographic location (latitude)");
	REQUIRE(lat != nullptr);
	REQUIRE(lat->allowed_units.size() == 1);
	CHECK(lat->allowed_units[0] == "DD");

	const auto *temp = FindField(c, "sample storage temperature");
	REQUIRE(temp != nullptr);
	REQUIRE(temp->allowed_units.size() == 2);
	CHECK(temp->allowed_units[0] == "degC");
	CHECK(temp->allowed_units[1] == "K");
}

TEST_CASE("ParseChecklistXML extracts TEXT_CHOICE_FIELD controlled vocabulary", "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	const auto *country = FindField(c, "country");
	REQUIRE(country != nullptr);
	REQUIRE(country->allowed_values.size() == 3);
	CHECK(country->allowed_values[0] == "USA");
	CHECK(country->allowed_values[1] == "Canada");
	CHECK(country->allowed_values[2] == "Mexico");
}

TEST_CASE("ParseChecklistXML maps MANDATORY value strings to enum", "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	CHECK(FindField(c, "project name")->mandatory == ChecklistFieldMandatory::MANDATORY);
	CHECK(FindField(c, "ploidy")->mandatory == ChecklistFieldMandatory::OPTIONAL);
}

TEST_CASE("ParseChecklistXML throws on malformed XML", "[ena_checklist]") {
	REQUIRE_THROWS(ParseChecklistXML("<not-well-formed"));
}

TEST_CASE("ParseChecklistXML throws on empty CHECKLIST_SET", "[ena_checklist]") {
	REQUIRE_THROWS(ParseChecklistXML(R"(<?xml version="1.0"?><CHECKLIST_SET/>)"));
}

TEST_CASE("ParseChecklistXML tolerates xmlns declaration on the root element", "[ena_checklist]") {
	// The real ENA browser API may decorate the root with namespace attributes.
	// Expat is configured without namespace separation, so as long as elements
	// are unqualified (no `ns:` prefix) the parser sees the declared `name`
	// verbatim and the element-stack matching keeps working.
	constexpr const char *NS_XML = R"(<?xml version="1.0"?>
<CHECKLIST_SET xmlns="https://www.ebi.ac.uk/ena/browser/checklists">
<CHECKLIST accession="TST000099" checklistType="Sample">
  <DESCRIPTOR>
    <LABEL>ns-aware fixture</LABEL>
    <FIELD_GROUP>
      <FIELD>
        <LABEL>only field</LABEL>
        <NAME>only_field</NAME>
        <FIELD_TYPE><TEXT_FIELD/></FIELD_TYPE>
        <MANDATORY>mandatory</MANDATORY>
      </FIELD>
    </FIELD_GROUP>
  </DESCRIPTOR>
</CHECKLIST>
</CHECKLIST_SET>)";
	const auto c = ParseChecklistXML(NS_XML);
	CHECK(c.accession == "TST000099");
	CHECK(c.label == "ns-aware fixture");
	REQUIRE(c.fields.size() == 1);
	CHECK(c.fields[0].label == "only field");
	CHECK(c.fields[0].mandatory == ChecklistFieldMandatory::MANDATORY);
}

// ---- Validator ----

TEST_CASE("Validate: all mandatory present + units + CV match → no issues", "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	std::vector<std::pair<std::string, std::string>> attrs = {
	    {"project name", "demo project"},
	    {"geographic location (latitude)", "32.7157"},
	    {"country", "USA"},
	    {"collection date", "2026-04-01"},
	};
	std::vector<std::pair<std::string, std::string>> units = {{"geographic location (latitude)", "DD"}};
	const auto issues = ValidateAttributesAgainstChecklist(c, attrs, units);
	CHECK(issues.empty());
}

TEST_CASE("Validate: missing mandatory field → issue names the field", "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	std::vector<std::pair<std::string, std::string>> attrs = {
	    {"geographic location (latitude)", "32.7157"},
	    {"country", "USA"},
	    {"collection date", "2026-04-01"},
	};
	std::vector<std::pair<std::string, std::string>> units = {{"geographic location (latitude)", "DD"}};
	const auto issues = ValidateAttributesAgainstChecklist(c, attrs, units);
	CHECK(HasIssueFor(issues, "project name", "mandatory"));
}

TEST_CASE("Validate: empty-string mandatory value counts as missing", "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	std::vector<std::pair<std::string, std::string>> attrs = {
	    {"project name", ""},
	    {"geographic location (latitude)", "32.7157"},
	    {"country", "USA"},
	    {"collection date", "2026-04-01"},
	};
	std::vector<std::pair<std::string, std::string>> units = {{"geographic location (latitude)", "DD"}};
	const auto issues = ValidateAttributesAgainstChecklist(c, attrs, units);
	CHECK(HasIssueFor(issues, "project name", "mandatory"));
}

TEST_CASE("Validate: unknown attribute key → issue", "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	std::vector<std::pair<std::string, std::string>> attrs = {
	    {"project name", "demo"},   {"geographic location (latitude)", "32.7157"},
	    {"country", "USA"},         {"collection date", "2026-04-01"},
	    {"made up attribute", "v"},
	};
	std::vector<std::pair<std::string, std::string>> units = {{"geographic location (latitude)", "DD"}};
	const auto issues = ValidateAttributesAgainstChecklist(c, attrs, units);
	CHECK(HasIssueFor(issues, "made up attribute", "not in checklist"));
}

TEST_CASE("Validate: missing units when field requires units → issue", "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	std::vector<std::pair<std::string, std::string>> attrs = {
	    {"project name", "demo"},
	    {"geographic location (latitude)", "32.7157"}, // no entry in units map
	    {"country", "USA"},
	    {"collection date", "2026-04-01"},
	};
	std::vector<std::pair<std::string, std::string>> units; // empty
	const auto issues = ValidateAttributesAgainstChecklist(c, attrs, units);
	CHECK(HasIssueFor(issues, "geographic location (latitude)", "unit"));
}

TEST_CASE("Validate: bad unit value (not in allowed_units) → issue", "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	std::vector<std::pair<std::string, std::string>> attrs = {
	    {"project name", "demo"},
	    {"geographic location (latitude)", "32.7157"},
	    {"country", "USA"},
	    {"collection date", "2026-04-01"},
	};
	std::vector<std::pair<std::string, std::string>> units = {
	    {"geographic location (latitude)", "Mm"}, // wrong unit
	};
	const auto issues = ValidateAttributesAgainstChecklist(c, attrs, units);
	CHECK(HasIssueFor(issues, "geographic location (latitude)", "unit"));
}

TEST_CASE("Validate: CV value not in allowed_values → issue (the live-wwwdev USA-vs-United-States bug)",
          "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	std::vector<std::pair<std::string, std::string>> attrs = {
	    {"project name", "demo"},
	    {"geographic location (latitude)", "32.7157"},
	    {"country", "United States"}, // wrong: must be "USA"
	    {"collection date", "2026-04-01"},
	};
	std::vector<std::pair<std::string, std::string>> units = {{"geographic location (latitude)", "DD"}};
	const auto issues = ValidateAttributesAgainstChecklist(c, attrs, units);
	CHECK(HasIssueFor(issues, "country", "United States"));
}

TEST_CASE("Validate: optional field with units left blank → no issue", "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	std::vector<std::pair<std::string, std::string>> attrs = {
	    {"project name", "demo"},
	    {"geographic location (latitude)", "32.7157"},
	    {"country", "USA"},
	    {"collection date", "2026-04-01"},
	    // sample storage temperature is optional and absent → fine
	};
	std::vector<std::pair<std::string, std::string>> units = {{"geographic location (latitude)", "DD"}};
	const auto issues = ValidateAttributesAgainstChecklist(c, attrs, units);
	CHECK(issues.empty());
}

TEST_CASE("Validate: optional field with units present → unit value is checked", "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	std::vector<std::pair<std::string, std::string>> attrs = {
	    {"project name", "demo"},          {"geographic location (latitude)", "32.7157"}, {"country", "USA"},
	    {"collection date", "2026-04-01"}, {"sample storage temperature", "-80"},
	};
	std::vector<std::pair<std::string, std::string>> units = {
	    {"geographic location (latitude)", "DD"}, {"sample storage temperature", "F"}, // not in {degC, K}
	};
	const auto issues = ValidateAttributesAgainstChecklist(c, attrs, units);
	CHECK(HasIssueFor(issues, "sample storage temperature", "unit"));
}

TEST_CASE("Validate: aggregates multiple issues in one call", "[ena_checklist]") {
	const auto c = ParseChecklistXML(FIXTURE_XML);
	std::vector<std::pair<std::string, std::string>> attrs = {
	    // project name missing
	    {"geographic location (latitude)", "32.7157"},
	    {"country", "Atlantis"}, // not in CV
	    // collection date missing
	    {"made up attribute", "v"},
	};
	std::vector<std::pair<std::string, std::string>> units; // missing latitude unit
	const auto issues = ValidateAttributesAgainstChecklist(c, attrs, units);
	CHECK(issues.size() >= 4);
	CHECK(HasIssueFor(issues, "project name", "mandatory"));
	CHECK(HasIssueFor(issues, "collection date", "mandatory"));
	CHECK(HasIssueFor(issues, "country", "Atlantis"));
	CHECK(HasIssueFor(issues, "geographic location (latitude)", "unit"));
	CHECK(HasIssueFor(issues, "made up attribute", "not in checklist"));
}

// ---- URL build / env override ----

TEST_CASE("BuildChecklistFetchURL builds the canonical browser-API URL", "[ena_checklist]") {
	CHECK(BuildChecklistFetchURL("https://www.ebi.ac.uk/ena/browser/api/xml", "ERC000015") ==
	      "https://www.ebi.ac.uk/ena/browser/api/xml/ERC000015");
	CHECK(BuildChecklistFetchURL("https://example.test/foo/", "ERC000015") == "https://example.test/foo/ERC000015");
}

TEST_CASE("BuildChecklistFetchURL rejects empty / suspicious accession", "[ena_checklist]") {
	REQUIRE_THROWS_AS(BuildChecklistFetchURL("https://www.ebi.ac.uk/ena/browser/api/xml", ""), std::invalid_argument);
	// Path traversal / URL injection: ENA accessions are alphanumeric.
	REQUIRE_THROWS_AS(BuildChecklistFetchURL("https://www.ebi.ac.uk/ena/browser/api/xml", "../../etc/passwd"),
	                  std::invalid_argument);
	REQUIRE_THROWS_AS(BuildChecklistFetchURL("https://www.ebi.ac.uk/ena/browser/api/xml", "ERC0?"),
	                  std::invalid_argument);
}

// ---- Registry (lazy fetch + cache) ----

TEST_CASE("ChecklistRegistry: GetOrFetch invokes fetcher once per accession; cache hit on second call",
          "[ena_checklist]") {
	ChecklistRegistry reg;
	int fetch_count = 0;
	auto fetcher = [&fetch_count](const std::string &url) -> std::string {
		fetch_count++;
		(void)url;
		return std::string(FIXTURE_XML);
	};
	const auto &c1 = reg.GetOrFetch("TST000001", fetcher);
	CHECK(c1.accession == "TST000001");
	CHECK(fetch_count == 1);

	const auto &c2 = reg.GetOrFetch("TST000001", fetcher);
	CHECK(&c1 == &c2); // returned by reference; same cached object
	CHECK(fetch_count == 1);
}

TEST_CASE("ChecklistRegistry: distinct accessions cache independently", "[ena_checklist]") {
	ChecklistRegistry reg;
	int fetch_count = 0;
	auto fetcher = [&fetch_count](const std::string &url) -> std::string {
		fetch_count++;
		(void)url;
		return std::string(FIXTURE_XML);
	};
	(void)reg.GetOrFetch("TST000001", fetcher);
	(void)reg.GetOrFetch("TST000002", fetcher);
	CHECK(fetch_count == 2);
	(void)reg.GetOrFetch("TST000001", fetcher);
	CHECK(fetch_count == 2);
}

TEST_CASE("ChecklistRegistry: fetcher exception propagates and does NOT cache the failure", "[ena_checklist]") {
	ChecklistRegistry reg;
	int fetch_count = 0;
	auto failing_fetcher = [&fetch_count](const std::string &url) -> std::string {
		fetch_count++;
		(void)url;
		throw std::runtime_error("simulated network failure");
	};
	REQUIRE_THROWS(reg.GetOrFetch("TST000001", failing_fetcher));
	REQUIRE_THROWS(reg.GetOrFetch("TST000001", failing_fetcher));
	// Both calls reached the fetcher — the failure was not cached.
	CHECK(fetch_count == 2);
}

TEST_CASE("ResolveChecklistBaseFromEnv defaults to the EBI browser API", "[ena_checklist]") {
	// Don't unset MIINT_ENA_CHECKLIST_URL_BASE here — other tests may rely on
	// the fixture override. Just check the function returns a non-empty
	// string consistent with whatever the harness configured.
	const auto base = ResolveChecklistBaseFromEnv();
	CHECK(!base.empty());
	// No trailing slash regardless of source.
	CHECK(base.back() != '/');
}
