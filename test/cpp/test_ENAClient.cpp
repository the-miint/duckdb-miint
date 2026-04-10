#include <catch2/catch_test_macros.hpp>
#include "ena_parser.hpp"

// ---- Step 1.1: Accession type detection ----

TEST_CASE("ENAClient accession type detection", "[ena]") {
	SECTION("Study accessions") {
		CHECK(miint::ENAParser::DetectAccessionType("PRJNA555783") == miint::ENAAccessionType::STUDY);
		CHECK(miint::ENAParser::DetectAccessionType("PRJEB5291") == miint::ENAAccessionType::STUDY);
		CHECK(miint::ENAParser::DetectAccessionType("PRJDB1234") == miint::ENAAccessionType::STUDY);
		CHECK(miint::ENAParser::DetectAccessionType("ERP005989") == miint::ENAAccessionType::STUDY);
		CHECK(miint::ENAParser::DetectAccessionType("SRP228341") == miint::ENAAccessionType::STUDY);
		CHECK(miint::ENAParser::DetectAccessionType("DRP000001") == miint::ENAAccessionType::STUDY);
	}

	SECTION("Sample accessions") {
		CHECK(miint::ENAParser::DetectAccessionType("SAMN12329471") == miint::ENAAccessionType::SAMPLE);
		CHECK(miint::ENAParser::DetectAccessionType("SAMEA3271045") == miint::ENAAccessionType::SAMPLE);
		CHECK(miint::ENAParser::DetectAccessionType("SAMD00000001") == miint::ENAAccessionType::SAMPLE);
	}

	SECTION("Run accessions") {
		CHECK(miint::ENAParser::DetectAccessionType("ERR1074767") == miint::ENAAccessionType::RUN);
		CHECK(miint::ENAParser::DetectAccessionType("SRR19396253") == miint::ENAAccessionType::RUN);
		CHECK(miint::ENAParser::DetectAccessionType("DRR000001") == miint::ENAAccessionType::RUN);
	}

	SECTION("Experiment accessions") {
		CHECK(miint::ENAParser::DetectAccessionType("ERX1106717") == miint::ENAAccessionType::EXPERIMENT);
		CHECK(miint::ENAParser::DetectAccessionType("SRX7123456") == miint::ENAAccessionType::EXPERIMENT);
		CHECK(miint::ENAParser::DetectAccessionType("DRX000001") == miint::ENAAccessionType::EXPERIMENT);
	}

	SECTION("Unknown accessions") {
		CHECK(miint::ENAParser::DetectAccessionType("") == miint::ENAAccessionType::UNKNOWN);
		CHECK(miint::ENAParser::DetectAccessionType("garbage") == miint::ENAAccessionType::UNKNOWN);
		CHECK(miint::ENAParser::DetectAccessionType("NC_001416.1") == miint::ENAAccessionType::UNKNOWN);
		CHECK(miint::ENAParser::DetectAccessionType("GCF_000005845.2") == miint::ENAAccessionType::UNKNOWN);
	}
}

// ---- Step 1.2: URL construction ----

TEST_CASE("ENAClient URL construction", "[ena]") {
	SECTION("Study accession uses study_accession query param") {
		auto url = miint::ENAParser::BuildSearchURL("PRJEB5291", "read_run", "run_accession,fastq_ftp");
		CHECK(url == "https://www.ebi.ac.uk/ena/portal/api/search?"
		             "result=read_run"
		             "&query=study_accession%3D%22PRJEB5291%22"
		             "&fields=run_accession,fastq_ftp"
		             "&limit=0&format=tsv");
	}

	SECTION("Run accession uses run_accession query param") {
		auto url = miint::ENAParser::BuildSearchURL("ERR1074767", "read_run", "run_accession,fastq_ftp");
		CHECK(url.find("run_accession%3D%22ERR1074767%22") != std::string::npos);
	}

	SECTION("Sample accession uses sample_accession query param") {
		auto url = miint::ENAParser::BuildSearchURL("SAMEA3271045", "read_run", "run_accession");
		CHECK(url.find("sample_accession%3D%22SAMEA3271045%22") != std::string::npos);
	}

	SECTION("Experiment accession uses experiment_accession query param") {
		auto url = miint::ENAParser::BuildSearchURL("ERX1106717", "read_run", "run_accession");
		CHECK(url.find("experiment_accession%3D%22ERX1106717%22") != std::string::npos);
	}

	SECTION("SRP study accession") {
		auto url = miint::ENAParser::BuildSearchURL("SRP228341", "read_run", "run_accession");
		CHECK(url.find("study_accession%3D%22SRP228341%22") != std::string::npos);
	}

	SECTION("Different result types are embedded") {
		auto url = miint::ENAParser::BuildSearchURL("ERR1074767", "sample", "sample_accession");
		CHECK(url.find("result=sample") != std::string::npos);
	}
}

// ---- Step 1.3: TSV parsing ----

TEST_CASE("ENAClient TSV parsing", "[ena]") {
	SECTION("Basic TSV with headers and one row") {
		std::string tsv = "run_accession\tsample_accession\tfastq_ftp\n"
		                  "ERR1074767\tSAMEA3271045\tftp.sra.ebi.ac.uk/vol1/fastq/ERR107/007/ERR1074767/"
		                  "ERR1074767.fastq.gz\n";

		auto result = miint::ENAParser::ParseTSV(tsv);

		REQUIRE(result.column_names.size() == 3);
		CHECK(result.column_names[0] == "run_accession");
		CHECK(result.column_names[1] == "sample_accession");
		CHECK(result.column_names[2] == "fastq_ftp");

		REQUIRE(result.rows.size() == 1);
		CHECK(result.rows[0][0] == "ERR1074767");
		CHECK(result.rows[0][1] == "SAMEA3271045");
	}

	SECTION("Empty TSV (headers only, no data rows)") {
		std::string tsv = "run_accession\tsample_accession\n";
		auto result = miint::ENAParser::ParseTSV(tsv);
		CHECK(result.column_names.size() == 2);
		CHECK(result.rows.empty());
	}

	SECTION("Multiple rows") {
		std::string tsv = "col1\tcol2\n"
		                  "a\tb\n"
		                  "c\td\n";
		auto result = miint::ENAParser::ParseTSV(tsv);
		REQUIRE(result.rows.size() == 2);
		CHECK(result.rows[0][0] == "a");
		CHECK(result.rows[0][1] == "b");
		CHECK(result.rows[1][0] == "c");
		CHECK(result.rows[1][1] == "d");
	}

	SECTION("Empty fields between tabs") {
		std::string tsv = "col1\tcol2\tcol3\n"
		                  "a\t\tc\n";
		auto result = miint::ENAParser::ParseTSV(tsv);
		REQUIRE(result.rows.size() == 1);
		CHECK(result.rows[0][0] == "a");
		CHECK(result.rows[0][1] == "");
		CHECK(result.rows[0][2] == "c");
	}

	SECTION("No trailing newline") {
		std::string tsv = "col1\n"
		                  "val1";
		auto result = miint::ENAParser::ParseTSV(tsv);
		REQUIRE(result.rows.size() == 1);
		CHECK(result.rows[0][0] == "val1");
	}

	SECTION("Empty body") {
		std::string tsv = "";
		auto result = miint::ENAParser::ParseTSV(tsv);
		CHECK(result.column_names.empty());
		CHECK(result.rows.empty());
	}
}

// ---- Step 1.4: XML sample attribute parsing ----

TEST_CASE("ENAClient XML sample attribute parsing", "[ena]") {
	SECTION("Basic sample attributes") {
		std::string xml =
		    R"(<?xml version="1.0" encoding="UTF-8"?>
<ROOT>
  <SAMPLE accession="SAMEA3271045" alias="test">
    <SAMPLE_ATTRIBUTES>
      <SAMPLE_ATTRIBUTE>
        <TAG>collection date</TAG>
        <VALUE>2013</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>geographic location (country and/or sea)</TAG>
        <VALUE>United Kingdom</VALUE>
      </SAMPLE_ATTRIBUTE>
    </SAMPLE_ATTRIBUTES>
  </SAMPLE>
</ROOT>)";

		auto attrs = miint::ENAParser::ParseSampleAttributesXML(xml);

		REQUIRE(attrs.size() == 2);
		CHECK(attrs[0].sample_accession == "SAMEA3271045");
		CHECK(attrs[0].tag == "collection date");
		CHECK(attrs[0].value == "2013");
		CHECK(attrs[1].sample_accession == "SAMEA3271045");
		CHECK(attrs[1].tag == "geographic location (country and/or sea)");
		CHECK(attrs[1].value == "United Kingdom");
	}

	SECTION("Multiple samples in one XML") {
		std::string xml =
		    R"(<?xml version="1.0" encoding="UTF-8"?>
<ROOT>
  <SAMPLE accession="SAMEA001">
    <SAMPLE_ATTRIBUTES>
      <SAMPLE_ATTRIBUTE><TAG>depth</TAG><VALUE>10m</VALUE></SAMPLE_ATTRIBUTE>
    </SAMPLE_ATTRIBUTES>
  </SAMPLE>
  <SAMPLE accession="SAMEA002">
    <SAMPLE_ATTRIBUTES>
      <SAMPLE_ATTRIBUTE><TAG>depth</TAG><VALUE>20m</VALUE></SAMPLE_ATTRIBUTE>
    </SAMPLE_ATTRIBUTES>
  </SAMPLE>
</ROOT>)";

		auto attrs = miint::ENAParser::ParseSampleAttributesXML(xml);
		REQUIRE(attrs.size() == 2);
		CHECK(attrs[0].sample_accession == "SAMEA001");
		CHECK(attrs[0].value == "10m");
		CHECK(attrs[1].sample_accession == "SAMEA002");
		CHECK(attrs[1].value == "20m");
	}

	SECTION("Empty attributes section") {
		std::string xml =
		    R"(<ROOT><SAMPLE accession="SAMEA001"><SAMPLE_ATTRIBUTES></SAMPLE_ATTRIBUTES></SAMPLE></ROOT>)";
		auto attrs = miint::ENAParser::ParseSampleAttributesXML(xml);
		CHECK(attrs.empty());
	}

	SECTION("Attribute with empty value") {
		std::string xml =
		    R"(<ROOT><SAMPLE accession="SAMEA001"><SAMPLE_ATTRIBUTES>
        <SAMPLE_ATTRIBUTE><TAG>host</TAG><VALUE></VALUE></SAMPLE_ATTRIBUTE>
    </SAMPLE_ATTRIBUTES></SAMPLE></ROOT>)";
		auto attrs = miint::ENAParser::ParseSampleAttributesXML(xml);
		REQUIRE(attrs.size() == 1);
		CHECK(attrs[0].tag == "host");
		CHECK(attrs[0].value == "");
	}

	SECTION("Attribute with no VALUE element") {
		std::string xml =
		    R"(<ROOT><SAMPLE accession="SAMEA001"><SAMPLE_ATTRIBUTES>
        <SAMPLE_ATTRIBUTE><TAG>orphan_tag</TAG></SAMPLE_ATTRIBUTE>
    </SAMPLE_ATTRIBUTES></SAMPLE></ROOT>)";
		auto attrs = miint::ENAParser::ParseSampleAttributesXML(xml);
		REQUIRE(attrs.size() == 1);
		CHECK(attrs[0].tag == "orphan_tag");
		CHECK(attrs[0].value == "");
	}
}

// ---- Step 1.5: FTP-to-HTTPS URL conversion ----

TEST_CASE("ENAClient FTP to HTTPS conversion", "[ena]") {
	SECTION("Single-end FTP path") {
		auto urls =
		    miint::ENAParser::FTPtoHTTPS("ftp.sra.ebi.ac.uk/vol1/fastq/ERR107/007/ERR1074767/ERR1074767.fastq.gz");
		REQUIRE(urls.size() == 1);
		CHECK(urls[0] == "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR107/007/ERR1074767/ERR1074767.fastq.gz");
	}

	SECTION("Paired-end FTP paths (semicolon-separated)") {
		auto urls =
		    miint::ENAParser::FTPtoHTTPS("ftp.sra.ebi.ac.uk/vol1/fastq/ERR107/007/ERR1074767/ERR1074767_1.fastq.gz;"
		                                 "ftp.sra.ebi.ac.uk/vol1/fastq/ERR107/007/ERR1074767/ERR1074767_2.fastq.gz");
		REQUIRE(urls.size() == 2);
		CHECK(urls[0] == "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR107/007/ERR1074767/ERR1074767_1.fastq.gz");
		CHECK(urls[1] == "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR107/007/ERR1074767/ERR1074767_2.fastq.gz");
	}

	SECTION("Empty FTP field") {
		auto urls = miint::ENAParser::FTPtoHTTPS("");
		CHECK(urls.empty());
	}
}

// ---- Step 1.6: Default fields ----

TEST_CASE("ENAClient default fields", "[ena]") {
	SECTION("Default read_run fields include key columns") {
		auto fields = miint::ENAParser::DefaultFields("read_run");
		CHECK(fields.find("run_accession") != std::string::npos);
		CHECK(fields.find("sample_accession") != std::string::npos);
		CHECK(fields.find("experiment_accession") != std::string::npos);
		CHECK(fields.find("study_accession") != std::string::npos);
		CHECK(fields.find("fastq_ftp") != std::string::npos);
		CHECK(fields.find("fastq_bytes") != std::string::npos);
		CHECK(fields.find("fastq_md5") != std::string::npos);
		CHECK(fields.find("library_layout") != std::string::npos);
		CHECK(fields.find("library_strategy") != std::string::npos);
		CHECK(fields.find("instrument_model") != std::string::npos);
		CHECK(fields.find("read_count") != std::string::npos);
		CHECK(fields.find("base_count") != std::string::npos);
	}

	SECTION("Default sample fields include key columns") {
		auto fields = miint::ENAParser::DefaultFields("sample");
		CHECK(fields.find("sample_accession") != std::string::npos);
		CHECK(fields.find("scientific_name") != std::string::npos);
		CHECK(fields.find("tax_id") != std::string::npos);
	}

	SECTION("Default study fields include key columns") {
		auto fields = miint::ENAParser::DefaultFields("study");
		CHECK(fields.find("study_accession") != std::string::npos);
		CHECK(fields.find("study_title") != std::string::npos);
	}

	SECTION("Unknown result type returns empty string") {
		auto fields = miint::ENAParser::DefaultFields("nonsense");
		CHECK(fields.empty());
	}
}

// ---- Step 1.7: XML URL construction ----

TEST_CASE("ENAClient XML URL construction", "[ena]") {
	SECTION("Single accession") {
		auto url = miint::ENAParser::BuildXMLURL({"SAMEA001"});
		CHECK(url == "https://www.ebi.ac.uk/ena/browser/api/xml/SAMEA001");
	}

	SECTION("Multiple accessions batched") {
		auto url = miint::ENAParser::BuildXMLURL({"SAMEA001", "SAMEA002", "SAMEA003"});
		CHECK(url == "https://www.ebi.ac.uk/ena/browser/api/xml/SAMEA001,SAMEA002,SAMEA003");
	}
}

// ---- Step 1.8: Input validation ----

TEST_CASE("ENAParser input validation", "[ena]") {
	SECTION("Valid accessions pass validation") {
		CHECK_NOTHROW(miint::ENAParser::ValidateAccession("ERR1074767"));
		CHECK_NOTHROW(miint::ENAParser::ValidateAccession("PRJNA555783"));
		CHECK_NOTHROW(miint::ENAParser::ValidateAccession("GCF_000005845.2"));
	}

	SECTION("Accessions with special characters are rejected") {
		CHECK_THROWS(miint::ENAParser::ValidateAccession("ERR107&limit=1"));
		CHECK_THROWS(miint::ENAParser::ValidateAccession("ERR107 extra"));
		CHECK_THROWS(miint::ENAParser::ValidateAccession("ERR107%22"));
		CHECK_THROWS(miint::ENAParser::ValidateAccession("ERR107=bad"));
	}

	SECTION("Valid fields pass validation") {
		CHECK_NOTHROW(miint::ENAParser::ValidateFields("run_accession,fastq_ftp"));
		CHECK_NOTHROW(miint::ENAParser::ValidateFields("sample_accession"));
	}

	SECTION("Fields with special characters are rejected") {
		CHECK_THROWS(miint::ENAParser::ValidateFields("run_accession&limit=1"));
		CHECK_THROWS(miint::ENAParser::ValidateFields("field with spaces"));
	}

	SECTION("FTPtoHTTPS strips ftp:// prefix") {
		auto urls = miint::ENAParser::FTPtoHTTPS("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR107/file.gz");
		REQUIRE(urls.size() == 1);
		CHECK(urls[0] == "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR107/file.gz");
	}
}
