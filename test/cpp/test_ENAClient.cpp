#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "ena_client.hpp"
#include "ena_parser.hpp"
#include "ena_resolver_cache.hpp"
#include "ena_run_info_extractor.hpp"

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

// ---- FilterSubmittedByFormat ----

TEST_CASE("FilterSubmittedByFormat extracts matching files", "[ena]") {
	SECTION("Two SFF files among mixed formats") {
		auto result = miint::ENAParser::FilterSubmittedByFormat(
		    "ftp.example.com/a.sff;ftp.example.com/b.fastq.gz;ftp.example.com/c.sff", "", "SFF;FASTQ;SFF",
		    "100;200;300", "SFF");
		REQUIRE(result.urls.size() == 2);
		CHECK(result.urls[0] == "https://ftp.example.com/a.sff");
		CHECK(result.urls[1] == "https://ftp.example.com/c.sff");
		CHECK(result.total_bytes == 400);
		CHECK(result.aspera_raw.empty());
	}

	SECTION("No matches returns empty") {
		auto result =
		    miint::ENAParser::FilterSubmittedByFormat("ftp.example.com/a.fastq.gz", "", "FASTQ", "100", "SFF");
		CHECK(result.urls.empty());
		CHECK(result.total_bytes == 0);
	}

	SECTION("Empty inputs returns empty") {
		auto result = miint::ENAParser::FilterSubmittedByFormat("", "", "", "", "SFF");
		CHECK(result.urls.empty());
		CHECK(result.total_bytes == 0);
	}

	SECTION("Mismatched field counts uses minimum") {
		// More formats than FTP entries — only process indices where all fields exist
		auto result =
		    miint::ENAParser::FilterSubmittedByFormat("ftp.example.com/a.sff", "", "SFF;SFF", "100;200", "SFF");
		REQUIRE(result.urls.size() == 1);
		CHECK(result.total_bytes == 100);
	}

	SECTION("Aspera paths extracted alongside URLs") {
		auto result = miint::ENAParser::FilterSubmittedByFormat(
		    "ftp.example.com/a.sff;ftp.example.com/b.fastq.gz",
		    "fasp.example.com:/vol1/a.sff;fasp.example.com:/vol1/b.fastq.gz", "SFF;FASTQ", "100;200", "SFF");
		REQUIRE(result.urls.size() == 1);
		CHECK(result.urls[0] == "https://ftp.example.com/a.sff");
		REQUIRE(result.aspera_raw.size() == 1);
		CHECK(result.aspera_raw[0] == "fasp.example.com:/vol1/a.sff");
		CHECK(result.total_bytes == 100);
	}

	SECTION("Case-sensitive format matching") {
		auto result = miint::ENAParser::FilterSubmittedByFormat("ftp.example.com/a.sff", "", "sff", "100", "SFF");
		CHECK(result.urls.empty());
	}

	SECTION("Bytes field with invalid number treated as 0") {
		auto result =
		    miint::ENAParser::FilterSubmittedByFormat("ftp.example.com/a.sff", "", "SFF", "notanumber", "SFF");
		REQUIRE(result.urls.size() == 1);
		CHECK(result.total_bytes == 0);
	}
}

// ---- BuildSearchURLBatch ----

TEST_CASE("BuildSearchURLBatch with single run accession uses IN clause", "[ena]") {
	auto url = miint::ENAParser::BuildSearchURLBatch({"ERR1074767"}, miint::ENAAccessionType::RUN, "read_run",
	                                                 "run_accession,fastq_ftp");
	// Single-element: query=run_accession IN ("ERR1074767")
	CHECK(url.find("run_accession%20IN%20%28%22ERR1074767%22%29") != std::string::npos);
	CHECK(url.find("result=read_run") != std::string::npos);
	CHECK(url.find("fields=run_accession,fastq_ftp") != std::string::npos);
}

TEST_CASE("BuildSearchURLBatch with multiple run accessions", "[ena]") {
	auto url = miint::ENAParser::BuildSearchURLBatch({"ERR1074767", "ERR1074768"}, miint::ENAAccessionType::RUN,
	                                                 "read_run", "run_accession,fastq_ftp");
	// query=run_accession IN ("ERR1074767","ERR1074768")
	CHECK(url.find("run_accession%20IN%20%28%22ERR1074767%22%2C%22ERR1074768%22%29") != std::string::npos);
}

TEST_CASE("BuildSearchURLBatch with study accessions uses study_accession column", "[ena]") {
	auto url = miint::ENAParser::BuildSearchURLBatch({"PRJEB5291", "PRJNA555783"}, miint::ENAAccessionType::STUDY,
	                                                 "read_run", "run_accession,study_accession");
	CHECK(url.find("study_accession%20IN%20%28%22PRJEB5291%22%2C%22PRJNA555783%22%29") != std::string::npos);
}

TEST_CASE("BuildSearchURLBatch URL-encodes parens, quotes, commas", "[ena]") {
	auto url = miint::ENAParser::BuildSearchURLBatch({"ERR1", "ERR2", "ERR3"}, miint::ENAAccessionType::RUN, "read_run",
	                                                 "run_accession");
	// %20 (space), %28 '(', %22 '"', %2C ',', %29 ')'
	CHECK(url.find("%20IN%20%28") != std::string::npos); // IN (
	CHECK(url.find("%22ERR1%22%2C%22ERR2%22%2C%22ERR3%22") != std::string::npos);
	CHECK(url.find("%29&fields=") != std::string::npos); // ) before fields
}

TEST_CASE("BuildSearchURLBatch throws on empty accession vector", "[ena]") {
	CHECK_THROWS(miint::ENAParser::BuildSearchURLBatch({}, miint::ENAAccessionType::RUN, "read_run", "run_accession"));
}

TEST_CASE("BuildSearchURLBatch throws on UNKNOWN accession type", "[ena]") {
	CHECK_THROWS(miint::ENAParser::BuildSearchURLBatch({"ERR1074767"}, miint::ENAAccessionType::UNKNOWN, "read_run",
	                                                   "run_accession"));
}

TEST_CASE("BuildSearchURLBatch validates every accession", "[ena]") {
	CHECK_THROWS(miint::ENAParser::BuildSearchURLBatch({"ERR1", "bad accession"}, miint::ENAAccessionType::RUN,
	                                                   "read_run", "run_accession"));
	CHECK_THROWS(miint::ENAParser::BuildSearchURLBatch({"ERR1", "ERR2&limit=1"}, miint::ENAAccessionType::RUN,
	                                                   "read_run", "run_accession"));
}

TEST_CASE("BuildSearchURLBatch starts with portal base", "[ena]") {
	auto url = miint::ENAParser::BuildSearchURLBatch({"ERR1074767"}, miint::ENAAccessionType::RUN, "read_run",
	                                                 "run_accession");
	CHECK(url.rfind("https://www.ebi.ac.uk/ena/portal/api/search?", 0) == 0);
	CHECK(url.find("&limit=0&format=tsv") != std::string::npos);
}

// ---- ENAColumnIndexMap ----

TEST_CASE("ENAColumnIndexMap maps column names to indices", "[ena]") {
	SECTION("All known columns are mapped") {
		std::vector<std::string> header = {"run_accession",    "sample_accession", "experiment_accession",
		                                   "study_accession",  "fastq_ftp",        "fastq_aspera",
		                                   "fastq_bytes",      "library_layout",   "submitted_ftp",
		                                   "submitted_aspera", "submitted_bytes",  "submitted_format"};
		auto cols = miint::ENAColumnIndexMap::FromHeader(header);
		CHECK(cols.run_accession == 0);
		CHECK(cols.sample_accession == 1);
		CHECK(cols.experiment_accession == 2);
		CHECK(cols.study_accession == 3);
		CHECK(cols.fastq_ftp == 4);
		CHECK(cols.fastq_aspera == 5);
		CHECK(cols.fastq_bytes == 6);
		CHECK(cols.library_layout == 7);
		CHECK(cols.submitted_ftp == 8);
		CHECK(cols.submitted_aspera == 9);
		CHECK(cols.submitted_bytes == 10);
		CHECK(cols.submitted_format == 11);
	}

	SECTION("Missing columns remain -1") {
		std::vector<std::string> header = {"run_accession", "fastq_ftp"};
		auto cols = miint::ENAColumnIndexMap::FromHeader(header);
		CHECK(cols.run_accession == 0);
		CHECK(cols.fastq_ftp == 1);
		CHECK(cols.sample_accession == -1);
		CHECK(cols.library_layout == -1);
	}

	SECTION("Unknown column names are ignored") {
		std::vector<std::string> header = {"run_accession", "irrelevant", "fastq_ftp"};
		auto cols = miint::ENAColumnIndexMap::FromHeader(header);
		CHECK(cols.run_accession == 0);
		CHECK(cols.fastq_ftp == 2);
	}

	SECTION("Get returns empty string for missing column") {
		std::vector<std::string> row = {"ERR1", "SAM1"};
		CHECK(miint::ENAColumnIndexMap::Get(row, -1) == "");
		CHECK(miint::ENAColumnIndexMap::Get(row, 99) == "");
		CHECK(miint::ENAColumnIndexMap::Get(row, 0) == "ERR1");
	}
}

// ---- ENARunInfoExtractor ----

static miint::ENAColumnIndexMap StandardCols() {
	std::vector<std::string> header = {"run_accession",    "sample_accession", "experiment_accession", "fastq_ftp",
	                                   "fastq_aspera",     "fastq_bytes",      "library_layout",       "submitted_ftp",
	                                   "submitted_aspera", "submitted_bytes",  "submitted_format"};
	return miint::ENAColumnIndexMap::FromHeader(header);
}

TEST_CASE("ENARunInfoExtractor extracts single-end FASTQ", "[ena]") {
	auto cols = StandardCols();
	std::vector<std::string> row = {"ERR1074767",
	                                "SAMEA3271045",
	                                "ERX1106717",
	                                "ftp.sra.ebi.ac.uk/vol1/fastq/ERR107/007/ERR1074767/ERR1074767.fastq.gz",
	                                "fasp.sra.ebi.ac.uk:/vol1/fastq/ERR107/007/ERR1074767/ERR1074767.fastq.gz",
	                                "12345",
	                                "SINGLE",
	                                "",
	                                "",
	                                "",
	                                ""};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "auto");
	REQUIRE(runs.size() == 1);
	CHECK(runs[0].run_accession == "ERR1074767");
	CHECK(runs[0].sample_accession == "SAMEA3271045");
	CHECK(runs[0].experiment_accession == "ERX1106717");
	CHECK(runs[0].is_paired == false);
	REQUIRE(runs[0].fastq_urls.size() == 1);
	CHECK(runs[0].fastq_urls[0] == "https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR107/007/ERR1074767/ERR1074767.fastq.gz");
	CHECK(runs[0].total_bytes == 12345);
	CHECK(runs[0].format == miint::ENASequenceFormat::FASTX);
	CHECK(runs[0].sff_url.empty());
	REQUIRE(runs[0].aspera_paths.size() == 1);
	CHECK(runs[0].aspera_paths[0].host == "fasp.sra.ebi.ac.uk");
}

TEST_CASE("ENARunInfoExtractor extracts paired-end FASTQ with 2 files", "[ena]") {
	auto cols = StandardCols();
	std::vector<std::string> row = {"ERR1000",
	                                "SAM1",
	                                "EXP1",
	                                "ftp.example.com/ERR1000_1.fastq.gz;ftp.example.com/ERR1000_2.fastq.gz",
	                                "fasp.example.com:/a_1.fastq.gz;fasp.example.com:/a_2.fastq.gz",
	                                "100;200",
	                                "PAIRED",
	                                "",
	                                "",
	                                "",
	                                ""};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "auto");
	REQUIRE(runs.size() == 1);
	CHECK(runs[0].is_paired == true);
	REQUIRE(runs[0].fastq_urls.size() == 2);
	CHECK(runs[0].total_bytes == 300);
	REQUIRE(runs[0].aspera_paths.size() == 2);
}

TEST_CASE("ENARunInfoExtractor filters 3-file paired-end to _1/_2", "[ena]") {
	auto cols = StandardCols();
	// ENA sometimes returns 3 files for PAIRED layout: an orphan single plus _1 and _2.
	// The extractor must drop the orphan.
	std::vector<std::string> row = {
	    "ERR2000",
	    "SAM2",
	    "EXP2",
	    "ftp.example.com/ERR2000.fastq.gz;ftp.example.com/ERR2000_1.fastq.gz;ftp.example.com/ERR2000_2.fastq.gz",
	    "fasp.example.com:/ERR2000.fastq.gz;fasp.example.com:/ERR2000_1.fastq.gz;fasp.example.com:/ERR2000_2.fastq.gz",
	    "50;100;200",
	    "PAIRED",
	    "",
	    "",
	    "",
	    ""};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "auto");
	REQUIRE(runs.size() == 1);
	CHECK(runs[0].is_paired == true);
	REQUIRE(runs[0].fastq_urls.size() == 2);
	CHECK(runs[0].fastq_urls[0].find("_1.fast") != std::string::npos);
	CHECK(runs[0].fastq_urls[1].find("_2.fast") != std::string::npos);
	// Total bytes should reflect only the _1/_2 files (100+200), not the orphan
	CHECK(runs[0].total_bytes == 300);
	REQUIRE(runs[0].aspera_paths.size() == 2);
	CHECK(runs[0].aspera_paths[0].remote_path.find("_1.fast") != std::string::npos);
	CHECK(runs[0].aspera_paths[1].remote_path.find("_2.fast") != std::string::npos);
}

TEST_CASE("ENARunInfoExtractor picks SFF when prefer_format='sff' and SFF available", "[ena]") {
	auto cols = StandardCols();
	std::vector<std::string> row = {"ERR3000",
	                                "SAM3",
	                                "EXP3",
	                                "ftp.example.com/ERR3000.fastq.gz", // FASTQ also present
	                                "",
	                                "500",
	                                "SINGLE",
	                                "ftp.example.com/a.sff",
	                                "",
	                                "1000",
	                                "SFF"};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "sff");
	REQUIRE(runs.size() == 1);
	CHECK(runs[0].format == miint::ENASequenceFormat::SFF);
	CHECK(runs[0].sff_url == "https://ftp.example.com/a.sff");
	CHECK(runs[0].total_bytes == 1000);
	CHECK(runs[0].fastq_urls.empty());
}

TEST_CASE("ENARunInfoExtractor prefer_format='sff' falls back to FASTQ when no SFF", "[ena]") {
	auto cols = StandardCols();
	std::vector<std::string> row = {
	    "ERR3001", "SAM3", "EXP3", "ftp.example.com/ERR3001.fastq.gz", "", "500", "SINGLE", "", "", "", ""};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "sff");
	REQUIRE(runs.size() == 1);
	CHECK(runs[0].format == miint::ENASequenceFormat::FASTX);
	REQUIRE(runs[0].fastq_urls.size() == 1);
}

TEST_CASE("ENARunInfoExtractor prefer_format='fastq' skips SFF-only row", "[ena]") {
	auto cols = StandardCols();
	std::vector<std::string> row = {"ERR3002", "SAM3", "EXP3", "", "", "", "SINGLE", "ftp.example.com/a.sff",
	                                "",        "1000", "SFF"};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "fastq");
	CHECK(runs.empty());
}

TEST_CASE("ENARunInfoExtractor prefer_format='auto' picks FASTQ when both present", "[ena]") {
	auto cols = StandardCols();
	std::vector<std::string> row = {"ERR3003", "SAM3", "EXP3",   "ftp.example.com/ERR3003.fastq.gz",
	                                "",        "500",  "SINGLE", "ftp.example.com/a.sff",
	                                "",        "1000", "SFF"};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "auto");
	REQUIRE(runs.size() == 1);
	CHECK(runs[0].format == miint::ENASequenceFormat::FASTX);
}

TEST_CASE("ENARunInfoExtractor flattens multiple SFF URLs into one RunInfo each", "[ena]") {
	auto cols = StandardCols();
	std::vector<std::string> row = {"ERR4000", "SAM4",    "EXP4",   "",
	                                "",        "",        "SINGLE", "ftp.example.com/a.sff;ftp.example.com/b.sff",
	                                "",        "400;600", "SFF;SFF"};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "auto");
	REQUIRE(runs.size() == 2);
	CHECK(runs[0].format == miint::ENASequenceFormat::SFF);
	CHECK(runs[0].sff_url == "https://ftp.example.com/a.sff");
	CHECK(runs[1].format == miint::ENASequenceFormat::SFF);
	CHECK(runs[1].sff_url == "https://ftp.example.com/b.sff");
	// Total SFF bytes (1000) split equally across files
	CHECK(runs[0].total_bytes == 500);
	CHECK(runs[1].total_bytes == 500);
	// SFF uses same metadata for all flattened entries
	CHECK(runs[0].run_accession == "ERR4000");
	CHECK(runs[1].run_accession == "ERR4000");
}

TEST_CASE("ENARunInfoExtractor returns empty for row with no usable format", "[ena]") {
	auto cols = StandardCols();
	std::vector<std::string> row = {"ERR5000", "SAM5", "EXP5", "", "", "", "SINGLE", "", "", "", ""};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "auto");
	CHECK(runs.empty());
}

TEST_CASE("ENARunInfoExtractor bytes field with invalid number treated as 0", "[ena]") {
	auto cols = StandardCols();
	std::vector<std::string> row = {
	    "ERR6000", "SAM6", "EXP6", "ftp.example.com/ERR6000.fastq.gz", "", "notanumber", "SINGLE", "", "", "", ""};

	auto runs = miint::ENARunInfoExtractor::FromTSVRow(row, cols, "auto");
	REQUIRE(runs.size() == 1);
	CHECK(runs[0].total_bytes == 0);
}

// ---- ResolveRunsBatch ----
//
// MockFetcher records every URL it receives and returns TSV bodies keyed by
// substring match. Lets tests assert which URLs were fetched and how many times.
// State lives in a std::shared_ptr so that std::function's copy-on-wrap still
// mutates the same state the test harness inspects.
struct MockFetcher {
	struct State {
		std::vector<std::string> urls;
		std::vector<std::pair<std::string, std::string>> responses;
	};
	std::shared_ptr<State> state {std::make_shared<State>()};

	std::string operator()(const std::string &url) const {
		state->urls.push_back(url);
		for (const auto &p : state->responses) {
			if (url.find(p.first) != std::string::npos) {
				return p.second;
			}
		}
		return "run_accession\tstudy_accession\tfastq_ftp\n";
	}

	void AddResponse(const std::string &needle, const std::string &body) {
		state->responses.push_back({needle, body});
	}

	const std::vector<std::string> &urls() const {
		return state->urls;
	}

	int CountContaining(const std::string &needle) const {
		int n = 0;
		for (const auto &u : state->urls) {
			if (u.find(needle) != std::string::npos) {
				n++;
			}
		}
		return n;
	}
};

// Build a minimal single-row TSV response matching BATCH_FIELDS.
// Columns (12): run_accession,sample_accession,experiment_accession,study_accession,
//               fastq_ftp,fastq_aspera,fastq_bytes,library_layout,
//               submitted_ftp,submitted_aspera,submitted_bytes,submitted_format
static std::string RunRowTSV(const std::string &run_accession, const std::string &study_accession = "") {
	std::string header = "run_accession\tsample_accession\texperiment_accession\tstudy_accession\t"
	                     "fastq_ftp\tfastq_aspera\tfastq_bytes\tlibrary_layout\t"
	                     "submitted_ftp\tsubmitted_aspera\tsubmitted_bytes\tsubmitted_format\n";
	std::string ftp = "ftp.sra.ebi.ac.uk/vol1/fastq/" + run_accession + "/" + run_accession + ".fastq.gz";
	return header + run_accession + "\tSAM1\tEXP1\t" + study_accession + "\t" + ftp + "\t\t100\tSINGLE\t\t\t\t\n";
}

TEST_CASE("ResolveRunsBatch groups by accession type", "[ena]") {
	miint::ENAResolverCache cache(16);
	MockFetcher fetch;
	fetch.AddResponse("run_accession%20IN", RunRowTSV("ERR1") + "ERR2\tSAM1\tEXP1\t\t" +
	                                            "ftp.example.com/ERR2.fastq.gz\t\t200\tSINGLE\t\t\t\t\n");
	// For study query, return a row that has study_accession so grouping works
	fetch.AddResponse("study_accession%20IN",
	                  "run_accession\tsample_accession\texperiment_accession\tstudy_accession\tfastq_ftp\t"
	                  "fastq_aspera\tfastq_bytes\tlibrary_layout\tsubmitted_ftp\tsubmitted_aspera\t"
	                  "submitted_bytes\tsubmitted_format\n"
	                  "ERR9\tSAM9\tEXP9\tPRJNA3\tftp.example.com/ERR9.fastq.gz\t\t300\tSINGLE\t\t\t\t\n");

	auto result = miint::ResolveRunsBatch(fetch, cache, {"ERR1", "ERR2", "PRJNA3"}, "auto");

	// Exactly two fetches: one run batch, one study batch.
	CHECK(fetch.urls().size() == 2);
	CHECK(fetch.CountContaining("run_accession%20IN") == 1);
	CHECK(fetch.CountContaining("study_accession%20IN") == 1);

	REQUIRE(result.count("ERR1") == 1);
	REQUIRE(result.count("ERR2") == 1);
	REQUIRE(result.count("PRJNA3") == 1);
	CHECK(result["ERR1"].size() == 1);
	CHECK(result["ERR2"].size() == 1);
	CHECK(result["PRJNA3"].size() == 1);
	CHECK(result["PRJNA3"][0].run_accession == "ERR9");
}

TEST_CASE("ResolveRunsBatch respects max_batch_size", "[ena]") {
	miint::ENAResolverCache cache(256);
	MockFetcher fetch;

	std::vector<std::string> inputs;
	for (int i = 0; i < 75; i++) {
		inputs.push_back("ERR" + std::to_string(1000 + i));
	}

	miint::ENAFetcher fetcher = [&fetch](const std::string &url) {
		return fetch(url);
	};
	miint::ResolveRunsBatch(fetcher, cache, inputs, "auto", 50);

	// 75 accessions / 50 per batch = 2 fetches
	CHECK(fetch.urls().size() == 2);
	CHECK(fetch.CountContaining("run_accession%20IN") == 2);
}

TEST_CASE("ResolveRunsBatch populates cache for all input accessions", "[ena]") {
	miint::ENAResolverCache cache(16);
	MockFetcher fetch;
	fetch.AddResponse("run_accession%20IN", RunRowTSV("ERR1"));

	miint::ENAFetcher fetcher = [&fetch](const std::string &url) {
		return fetch(url);
	};
	miint::ResolveRunsBatch(fetcher, cache, {"ERR1", "ERR2"}, "auto");

	// Both should be cached (ERR1 with data, ERR2 as negative cache entry)
	std::vector<miint::ENARunInfo> tmp;
	CHECK(cache.Get({"ERR1", "auto"}, tmp));
	CHECK(tmp.size() == 1);
	CHECK(cache.Get({"ERR2", "auto"}, tmp));
	CHECK(tmp.empty());
}

TEST_CASE("ResolveRunsBatch serves from cache without calling fetcher", "[ena]") {
	miint::ENAResolverCache cache(16);
	miint::ENARunInfo pre;
	pre.run_accession = "ERR1";
	cache.Put({"ERR1", "auto"}, {pre});

	MockFetcher fetch;
	auto result = miint::ResolveRunsBatch(fetch, cache, {"ERR1"}, "auto");

	CHECK(fetch.urls().empty()); // no HTTP fetch
	REQUIRE(result.count("ERR1") == 1);
	CHECK(result["ERR1"].size() == 1);
}

TEST_CASE("ResolveRunsBatch mixes cached + uncached accessions", "[ena]") {
	miint::ENAResolverCache cache(16);
	miint::ENARunInfo pre;
	pre.run_accession = "ERR1";
	cache.Put({"ERR1", "auto"}, {pre});

	MockFetcher fetch;
	fetch.AddResponse("run_accession%20IN", RunRowTSV("ERR2"));

	auto result = miint::ResolveRunsBatch(fetch, cache, {"ERR1", "ERR2"}, "auto");

	// Only one fetch for ERR2
	CHECK(fetch.urls().size() == 1);
	// Fetch URL should contain ERR2 but not ERR1
	CHECK(fetch.urls()[0].find("ERR2") != std::string::npos);
	CHECK(fetch.urls()[0].find("ERR1") == std::string::npos);
	CHECK(result["ERR1"].size() == 1);
	CHECK(result["ERR2"].size() == 1);
}

TEST_CASE("ResolveRunsBatch negative-caches empty results", "[ena]") {
	miint::ENAResolverCache cache(16);
	MockFetcher fetch;
	// Fetcher returns header-only TSV → no rows → empty

	miint::ResolveRunsBatch(fetch, cache, {"ERR_MISSING"}, "auto");
	CHECK(fetch.urls().size() == 1);

	// Second call should hit cache, not fetch
	miint::ResolveRunsBatch(fetch, cache, {"ERR_MISSING"}, "auto");
	CHECK(fetch.urls().size() == 1);
}

TEST_CASE("ResolveRunsBatch deduplicates input accessions", "[ena]") {
	miint::ENAResolverCache cache(16);
	MockFetcher fetch;
	fetch.AddResponse("run_accession%20IN", RunRowTSV("ERR1"));

	auto result = miint::ResolveRunsBatch(fetch, cache, {"ERR1", "ERR1", "ERR1"}, "auto");

	CHECK(fetch.urls().size() == 1);
	// URL should include ERR1 only once in the IN clause
	size_t n = 0;
	const auto &url = fetch.urls()[0];
	auto pos = url.find("%22ERR1%22");
	while (pos != std::string::npos) {
		n++;
		pos = url.find("%22ERR1%22", pos + 1);
	}
	CHECK(n == 1);
	REQUIRE(result.count("ERR1") == 1);
}

TEST_CASE("ResolveRunsBatch drops rows whose grouping key matches no input", "[ena]") {
	// ENA sometimes returns extra rows (e.g., alternative run accessions under the same
	// sample or study). FetchBatch must NOT cache those under their grouping keys, because
	// doing so would pollute the cache with entries the caller never asked for and could
	// create false-positive hits in future queries.
	miint::ENAResolverCache cache(16);
	MockFetcher fetch;

	// Batch query for ERR1 — server returns an extra row for ERR_UNEXPECTED.
	std::string tsv = "run_accession\tsample_accession\texperiment_accession\tstudy_accession\t"
	                  "fastq_ftp\tfastq_aspera\tfastq_bytes\tlibrary_layout\t"
	                  "submitted_ftp\tsubmitted_aspera\tsubmitted_bytes\tsubmitted_format\n"
	                  "ERR1\tSAM1\tEXP1\t\tftp.example.com/ERR1.fastq.gz\t\t100\tSINGLE\t\t\t\t\n"
	                  "ERR_UNEXPECTED\tSAM2\tEXP2\t\tftp.example.com/ERR_UNEXPECTED.fastq.gz\t\t200\tSINGLE\t\t\t\t\n";
	fetch.AddResponse("run_accession%20IN", tsv);

	auto result = miint::ResolveRunsBatch(fetch, cache, {"ERR1"}, "auto");

	// Only ERR1 should be in the result
	CHECK(result.count("ERR1") == 1);
	CHECK(result["ERR1"].size() == 1);
	CHECK(result.count("ERR_UNEXPECTED") == 0);

	// Cache must not hold an entry for ERR_UNEXPECTED
	std::vector<miint::ENARunInfo> tmp;
	CHECK_FALSE(cache.Get({"ERR_UNEXPECTED", "auto"}, tmp));
	CHECK(cache.Get({"ERR1", "auto"}, tmp));
}

TEST_CASE("ResolveRunsBatch throws when batch response is missing grouping column", "[ena]") {
	// If ENA returns a TSV that doesn't include the column we grouped by (e.g., an API
	// change or server error), we must throw rather than silently report no results.
	miint::ENAResolverCache cache(16);
	MockFetcher fetch;
	// Response with columns but none that match run_accession (the grouping column for RUN type)
	fetch.AddResponse("run_accession%20IN", "sample_accession\tfastq_ftp\nSAM1\tftp.example.com/a.fastq.gz\n");

	CHECK_THROWS(miint::ResolveRunsBatch(fetch, cache, {"ERR1"}, "auto"));
}

TEST_CASE("ResolveRunsBatch: UNKNOWN-type falls back to per-accession fetch", "[ena]") {
	miint::ENAResolverCache cache(16);
	MockFetcher fetch;
	// Unknown-type accessions go through BuildSearchURL (single-accession path), which
	// uses %3D%22 (equals) rather than %20IN%20%28 (compound).

	miint::ENAFetcher fetcher = [&fetch](const std::string &url) {
		return fetch(url);
	};
	miint::ResolveRunsBatch(fetcher, cache, {"unknown-acc1", "unknown-acc2"}, "auto");

	CHECK(fetch.urls().size() == 2); // one fetch per accession
	for (const auto &u : fetch.urls()) {
		CHECK(u.find("%3D%22") != std::string::npos);      // equality query
		CHECK(u.find("%20IN%20%28") == std::string::npos); // NOT a batched query
	}
}

// ---- HTTP Basic auth header construction (Phase 2) ----
//
// Round-trip POST testing (request body, response handling, retry on 5xx)
// requires a real HTTP endpoint and lives in the Phase 4 mock-server tier.
// Here we cover only the pure-function piece: RFC 7617 Base64 encoding.

TEST_CASE("ENAClient BuildBasicAuthHeader: RFC 7617 canonical example", "[ena_basic_auth]") {
	// RFC 7617 §2: "Aladdin:open sesame" -> "QWxhZGRpbjpvcGVuIHNlc2FtZQ=="
	auto header = miint::BuildBasicAuthHeader("Aladdin", "open sesame");
	CHECK(header == "Basic QWxhZGRpbjpvcGVuIHNlc2FtZQ==");
}

TEST_CASE("ENAClient BuildBasicAuthHeader: simple no-padding case", "[ena_basic_auth]") {
	// "user:pass" is 9 bytes = 12 base64 chars exactly, no padding
	auto header = miint::BuildBasicAuthHeader("user", "pass");
	CHECK(header == "Basic dXNlcjpwYXNz");
}

TEST_CASE("ENAClient BuildBasicAuthHeader: padding cases", "[ena_basic_auth]") {
	// 4-byte input -> 2 padding chars: "ab:c" -> "YWI6Yw=="
	CHECK(miint::BuildBasicAuthHeader("ab", "c") == "Basic YWI6Yw==");
	// 5-byte input -> 1 padding char: "abc:d" -> "YWJjOmQ="
	CHECK(miint::BuildBasicAuthHeader("abc", "d") == "Basic YWJjOmQ=");
	// 6-byte input -> 0 padding chars: "abc:de" -> "YWJjOmRl"
	CHECK(miint::BuildBasicAuthHeader("abc", "de") == "Basic YWJjOmRl");
}

TEST_CASE("ENAClient BuildBasicAuthHeader: realistic Webin id", "[ena_basic_auth]") {
	// Sanity: prefix is 'Basic ', payload is exactly the right length.
	auto header = miint::BuildBasicAuthHeader("Webin-12345", "shh");
	CHECK_THAT(header, Catch::Matchers::StartsWith("Basic "));
	// "Webin-12345:shh" is 15 bytes -> ceil(15/3)*4 = 20 base64 chars, no padding
	CHECK(header.size() == 6 + 20);
	CHECK(header.find('=') == std::string::npos);
}

TEST_CASE("ENAClient BuildBasicAuthHeader: bytes outside printable ASCII", "[ena_basic_auth]") {
	// Webin passwords can contain any ASCII; the Base64 encoder MUST treat
	// inputs as raw bytes, not text. Use \x80 (high bit set) to verify
	// no UTF-8 / signed-char surprise.
	auto header = miint::BuildBasicAuthHeader("u", std::string("\x80\x81\x82", 3));
	// "u:\x80\x81\x82" = bytes 0x75 0x3A 0x80 0x81 0x82 (5 bytes)
	// Expected base64: "dTqAgYI=" (1 padding char)
	CHECK(header == "Basic dTqAgYI=");
}

TEST_CASE("ENAClient BuildBasicAuthHeader: empty user or password rejected", "[ena_basic_auth]") {
	CHECK_THROWS_WITH(miint::BuildBasicAuthHeader("", "pw"), Catch::Matchers::ContainsSubstring("user"));
	CHECK_THROWS_WITH(miint::BuildBasicAuthHeader("user", ""), Catch::Matchers::ContainsSubstring("password"));
}

TEST_CASE("ENAClient BuildBasicAuthHeader: colon in user rejected (RFC 7617)", "[ena_basic_auth]") {
	// User containing ':' would split incorrectly server-side and silently
	// authenticate as the wrong identity. RFC 7617 §2.1 forbids it.
	CHECK_THROWS_WITH(miint::BuildBasicAuthHeader("us:er", "pw"), Catch::Matchers::ContainsSubstring("':'"));
	// Colons in passwords are fine.
	CHECK_NOTHROW(miint::BuildBasicAuthHeader("user", "p:ass"));
}

TEST_CASE("ENAClient Base64Encode: empty input returns empty", "[ena_basic_auth]") {
	CHECK(miint::Base64Encode("") == "");
}

TEST_CASE("ENAClient IsRetryableStatus: 4xx not retried, 5xx and 429 retried", "[ena_retry]") {
	CHECK(miint::ENAClient::IsRetryableStatus(429));
	CHECK(miint::ENAClient::IsRetryableStatus(500));
	CHECK(miint::ENAClient::IsRetryableStatus(502));
	CHECK(miint::ENAClient::IsRetryableStatus(503));
	CHECK(miint::ENAClient::IsRetryableStatus(504));

	CHECK_FALSE(miint::ENAClient::IsRetryableStatus(200));
	CHECK_FALSE(miint::ENAClient::IsRetryableStatus(201));
	CHECK_FALSE(miint::ENAClient::IsRetryableStatus(400));
	CHECK_FALSE(miint::ENAClient::IsRetryableStatus(401));
	CHECK_FALSE(miint::ENAClient::IsRetryableStatus(403));
	CHECK_FALSE(miint::ENAClient::IsRetryableStatus(404));
	CHECK_FALSE(miint::ENAClient::IsRetryableStatus(415));
	CHECK_FALSE(miint::ENAClient::IsRetryableStatus(505));
}
