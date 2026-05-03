// SPDX-License-Identifier: MIT
//
// Phase 3 RED -> GREEN tests for ENA Webin V2 receipt parsing.
// Fixtures derive from the official ENA programmatic-submission tutorial
// captured in localdocs/ena-research-webin-v2-deep.md §1 and SRA.receipt.xsd.

#include "ena_receipt_parser.hpp"

#include <catch2/catch_all.hpp>

TEST_CASE("ENA receipt: success XML with project, sample, experiment, run, submission",
          "[ena_receipt]") {
	const char *xml = R"(<?xml version="1.0" encoding="UTF-8"?>
<RECEIPT receiptDate="2022-07-27T09:54:37.869+01:00"
         submissionFile="submission-EMBL-EBI_1658912077869.xml"
         success="true">
    <EXPERIMENT accession="ERX9535365" alias="illumina-hiSeq" status="PRIVATE"/>
    <RUN        accession="ERR9994219" alias="paired-data"    status="PRIVATE"/>
    <SAMPLE accession="ERS12520704" alias="gut-microbiota" status="PRIVATE"
            holdUntilDate="2024-07-12+01:00">
        <EXT_ID accession="SAMEA110422334" type="biosample"/>
    </SAMPLE>
    <PROJECT accession="PRJEB55033" alias="comparative-analysis"
             status="PRIVATE" holdUntilDate="2024-07-12+01:00">
        <EXT_ID accession="ERP139895" type="study"/>
    </PROJECT>
    <SUBMISSION accession="ERA16500666"
                alias="SUBMISSION-27-07-2022-09:54:36:278"/>
    <MESSAGES><INFO>All objects in this submission are set to private status (HOLD).</INFO></MESSAGES>
    <ACTIONS>ADD</ACTIONS>
    <ACTIONS>HOLD</ACTIONS>
</RECEIPT>)";

	auto receipt = miint::ParseReceiptXML(xml);

	CHECK(receipt.success == true);
	CHECK(receipt.receipt_date == "2022-07-27T09:54:37.869+01:00");
	CHECK(receipt.submission_file == "submission-EMBL-EBI_1658912077869.xml");
	CHECK(receipt.submission_accession == "ERA16500666");

	REQUIRE(receipt.actions.size() == 2);
	CHECK(receipt.actions[0] == "ADD");
	CHECK(receipt.actions[1] == "HOLD");

	REQUIRE(receipt.info_messages.size() == 1);
	CHECK(receipt.info_messages[0].find("private status (HOLD)") != std::string::npos);
	CHECK(receipt.errors.empty());

	// Find each object by type+alias (parse order is element order)
	REQUIRE(receipt.objects.size() == 4); // EXPERIMENT, RUN, SAMPLE, PROJECT (SUBMISSION goes to submission_accession)

	auto by_type_alias = [&](const std::string &type, const std::string &alias) -> const miint::ENAObjectReceipt * {
		for (const auto &o : receipt.objects) {
			if (o.object_type == type && o.alias == alias) {
				return &o;
			}
		}
		return nullptr;
	};

	auto *exp = by_type_alias("EXPERIMENT", "illumina-hiSeq");
	REQUIRE(exp != nullptr);
	CHECK(exp->accession == "ERX9535365");
	CHECK(exp->status == "PRIVATE");
	CHECK(exp->ext_ids.empty());

	auto *run = by_type_alias("RUN", "paired-data");
	REQUIRE(run != nullptr);
	CHECK(run->accession == "ERR9994219");

	auto *samp = by_type_alias("SAMPLE", "gut-microbiota");
	REQUIRE(samp != nullptr);
	CHECK(samp->accession == "ERS12520704");
	CHECK(samp->hold_until_date == "2024-07-12+01:00");
	REQUIRE(samp->ext_ids.size() == 1);
	CHECK(samp->ext_ids[0].accession == "SAMEA110422334");
	CHECK(samp->ext_ids[0].type == "biosample");

	auto *proj = by_type_alias("PROJECT", "comparative-analysis");
	REQUIRE(proj != nullptr);
	CHECK(proj->accession == "PRJEB55033");
	CHECK(proj->hold_until_date == "2024-07-12+01:00");
	REQUIRE(proj->ext_ids.size() == 1);
	CHECK(proj->ext_ids[0].accession == "ERP139895");
	CHECK(proj->ext_ids[0].type == "study");
}

TEST_CASE("ENA receipt: failure XML with success=false and ERROR messages",
          "[ena_receipt]") {
	const char *xml = R"(<?xml version="1.0" encoding="UTF-8"?>
<RECEIPT receiptDate="2022-01-01T17:05:01.114+01:00"
         submissionFile="failed_submission.xml"
         success="false">
    <SAMPLE alias="gut-001"/>
    <MESSAGES>
        <ERROR>Sample alias 'gut-001' is missing required attribute 'collection date'</ERROR>
        <ERROR>Sample alias 'gut-001' has invalid taxon_id</ERROR>
    </MESSAGES>
    <ACTIONS>ADD</ACTIONS>
</RECEIPT>)";

	auto receipt = miint::ParseReceiptXML(xml);

	CHECK(receipt.success == false);
	CHECK(receipt.submission_accession.empty());
	REQUIRE(receipt.errors.size() == 2);
	CHECK(receipt.errors[0].find("collection date") != std::string::npos);
	CHECK(receipt.errors[1].find("taxon_id") != std::string::npos);
	CHECK(receipt.info_messages.empty());

	// SAMPLE element is present but with no accession (failed objects can appear
	// in the receipt without server-assigned ids).
	REQUIRE(receipt.objects.size() == 1);
	CHECK(receipt.objects[0].object_type == "SAMPLE");
	CHECK(receipt.objects[0].alias == "gut-001");
	CHECK(receipt.objects[0].accession.empty());
}

TEST_CASE("ENA receipt: VALIDATE action receipt has no accessions", "[ena_receipt]") {
	const char *xml = R"(<?xml version="1.0" encoding="UTF-8"?>
<RECEIPT receiptDate="2022-01-01T00:00:00.000+01:00"
         submissionFile="validate.xml"
         success="true">
    <SAMPLE alias="gut-001"/>
    <ACTIONS>VALIDATE</ACTIONS>
</RECEIPT>)";

	auto receipt = miint::ParseReceiptXML(xml);

	CHECK(receipt.success == true);
	REQUIRE(receipt.actions.size() == 1);
	CHECK(receipt.actions[0] == "VALIDATE");
	REQUIRE(receipt.objects.size() == 1);
	CHECK(receipt.objects[0].alias == "gut-001");
	CHECK(receipt.objects[0].accession.empty());
}

TEST_CASE("ENA receipt: malformed XML throws", "[ena_receipt]") {
	CHECK_THROWS_WITH(miint::ParseReceiptXML("<unclosed-tag"),
	                  Catch::Matchers::ContainsSubstring("XML"));
}

TEST_CASE("ENA receipt: missing 'success' attribute defaults to false", "[ena_receipt]") {
	// Defensive: per XSD, success is required, but if a malformed server
	// response omits it we should default conservatively to failure.
	const char *xml = R"(<RECEIPT receiptDate="2022-01-01T00:00:00.000Z" submissionFile="x.xml">
        <MESSAGES><ERROR>broken receipt</ERROR></MESSAGES>
    </RECEIPT>)";
	auto receipt = miint::ParseReceiptXML(xml);
	CHECK(receipt.success == false);
}

TEST_CASE("ENA receipt: extracts SAMEA from a single sample", "[ena_receipt]") {
	// Most common consumer path: I just want SAMEA out of the receipt.
	const char *xml = R"(<RECEIPT receiptDate="2022-01-01T00:00:00.000Z"
                        submissionFile="x.xml" success="true">
        <SAMPLE accession="ERS27605861" alias="stomach_microbiota" status="PRIVATE"
                holdUntilDate="2023-01-01Z">
            <EXT_ID accession="SAMEA130793922" type="biosample"/>
        </SAMPLE>
        <SUBMISSION accession="ERA12956757" alias="SUBMISSION-X"/>
        <ACTIONS>ADD</ACTIONS>
    </RECEIPT>)";

	auto receipt = miint::ParseReceiptXML(xml);
	REQUIRE(receipt.success);
	REQUIRE(receipt.objects.size() == 1);
	const auto &samp = receipt.objects[0];
	CHECK(samp.accession == "ERS27605861");
	REQUIRE(samp.ext_ids.size() == 1);
	CHECK(samp.ext_ids[0].accession == "SAMEA130793922");
	CHECK(samp.ext_ids[0].type == "biosample");
}
