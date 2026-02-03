#include <catch2/catch_test_macros.hpp>
#include "ncbi_parser.hpp"
#include "zip_utils.hpp"
#include "miniz.hpp"
#include <cstring>

using namespace miint;

TEST_CASE("NCBIParser accession type detection", "[ncbi]") {
	SECTION("Assembly accessions (GCF_/GCA_)") {
		CHECK(NCBIParser::DetectAccessionType("GCF_000001405.40") == AccessionType::ASSEMBLY);
		CHECK(NCBIParser::DetectAccessionType("GCA_000001635.9") == AccessionType::ASSEMBLY);
		CHECK(NCBIParser::DetectAccessionType("GCF_009858895.2") == AccessionType::ASSEMBLY);

		CHECK(NCBIParser::IsAssemblyAccession("GCF_000001405.40") == true);
		CHECK(NCBIParser::IsAssemblyAccession("GCA_000001635.9") == true);
	}

	SECTION("Sequence accessions (NC_/NM_/NP_/etc)") {
		CHECK(NCBIParser::DetectAccessionType("NC_000001.11") == AccessionType::SEQUENCE);
		CHECK(NCBIParser::DetectAccessionType("NM_001101.5") == AccessionType::SEQUENCE);
		CHECK(NCBIParser::DetectAccessionType("NP_000509.1") == AccessionType::SEQUENCE);
		CHECK(NCBIParser::DetectAccessionType("AC_000001.1") == AccessionType::SEQUENCE);
		CHECK(NCBIParser::DetectAccessionType("NG_000001.1") == AccessionType::SEQUENCE);
		CHECK(NCBIParser::DetectAccessionType("NT_000001.1") == AccessionType::SEQUENCE);
		CHECK(NCBIParser::DetectAccessionType("NW_000001.1") == AccessionType::SEQUENCE);
		CHECK(NCBIParser::DetectAccessionType("NZ_CP001234.1") == AccessionType::SEQUENCE);

		CHECK(NCBIParser::IsAssemblyAccession("NC_000001.11") == false);
		CHECK(NCBIParser::IsAssemblyAccession("NM_001101.5") == false);
	}

	SECTION("Unknown accessions") {
		CHECK(NCBIParser::DetectAccessionType("") == AccessionType::UNKNOWN);
		CHECK(NCBIParser::DetectAccessionType("invalid") == AccessionType::UNKNOWN);
		CHECK(NCBIParser::DetectAccessionType("ABC123") == AccessionType::UNKNOWN);
	}
}

TEST_CASE("NCBIParser FASTA parsing", "[ncbi]") {
	SECTION("Single sequence") {
		std::string fasta = ">NC_001416.1 Escherichia phage Lambda, complete genome\n"
		                    "GGGCGGCGACCTCGCGGGTTTTCGCTATTTTAT\n"
		                    "GAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTT\n";

		auto batch = NCBIParser::ParseFasta(fasta);

		REQUIRE(batch.size() == 1);
		CHECK(batch.read_ids[0] == "NC_001416.1");
		CHECK(batch.comments[0] == "Escherichia phage Lambda, complete genome");
		CHECK(batch.sequences1[0] == "GGGCGGCGACCTCGCGGGTTTTCGCTATTTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTT");
		CHECK(batch.is_paired == false);
		CHECK(batch.quals1[0].as_string().empty()); // No quality for FASTA
	}

	SECTION("Multiple sequences") {
		std::string fasta = ">seq1 first sequence\n"
		                    "ATCGATCG\n"
		                    ">seq2 second sequence\n"
		                    "GCTAGCTA\n"
		                    ">seq3\n"
		                    "NNNNNNNN\n";

		auto batch = NCBIParser::ParseFasta(fasta);

		REQUIRE(batch.size() == 3);

		CHECK(batch.read_ids[0] == "seq1");
		CHECK(batch.comments[0] == "first sequence");
		CHECK(batch.sequences1[0] == "ATCGATCG");

		CHECK(batch.read_ids[1] == "seq2");
		CHECK(batch.comments[1] == "second sequence");
		CHECK(batch.sequences1[1] == "GCTAGCTA");

		CHECK(batch.read_ids[2] == "seq3");
		CHECK(batch.comments[2] == ""); // No comment
		CHECK(batch.sequences1[2] == "NNNNNNNN");
	}

	SECTION("Empty FASTA") {
		std::string fasta = "";
		auto batch = NCBIParser::ParseFasta(fasta);
		CHECK(batch.empty());
	}

	SECTION("FASTA with only whitespace") {
		std::string fasta = "  \n\n  \n";
		auto batch = NCBIParser::ParseFasta(fasta);
		CHECK(batch.empty());
	}

	SECTION("Multiline sequence") {
		std::string fasta = ">long_seq description here\n"
		                    "AAAAAAAAAA\n"
		                    "CCCCCCCCCC\n"
		                    "GGGGGGGGGG\n"
		                    "TTTTTTTTTT\n";

		auto batch = NCBIParser::ParseFasta(fasta);

		REQUIRE(batch.size() == 1);
		CHECK(batch.sequences1[0] == "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT");
	}
}

TEST_CASE("NCBIParser GenBank XML parsing", "[ncbi]") {
	SECTION("Basic GenBank XML") {
		// Simplified GenBank XML structure
		std::string xml = R"(<?xml version="1.0" encoding="UTF-8"?>
<GBSet>
  <GBSeq>
    <GBSeq_locus>NC_001416</GBSeq_locus>
    <GBSeq_length>48502</GBSeq_length>
    <GBSeq_moltype>DNA</GBSeq_moltype>
    <GBSeq_update-date>23-JAN-2024</GBSeq_update-date>
    <GBSeq_definition>Escherichia phage Lambda, complete genome</GBSeq_definition>
    <GBSeq_primary-accession>NC_001416</GBSeq_primary-accession>
    <GBSeq_accession-version>NC_001416.1</GBSeq_accession-version>
    <GBSeq_organism>Enterobacteria phage lambda</GBSeq_organism>
    <GBSeq_taxonomy>Viruses; Duplodnaviria; Heunggongvirae; Uroviricota</GBSeq_taxonomy>
    <GBSeq_source>Enterobacteria phage lambda</GBSeq_source>
    <GBSeq_feature-table>
      <GBFeature>
        <GBFeature_key>source</GBFeature_key>
        <GBFeature_location>1..48502</GBFeature_location>
        <GBFeature_quals>
          <GBQualifier>
            <GBQualifier_name>organism</GBQualifier_name>
            <GBQualifier_value>Enterobacteria phage lambda</GBQualifier_value>
          </GBQualifier>
          <GBQualifier>
            <GBQualifier_name>db_xref</GBQualifier_name>
            <GBQualifier_value>taxon:10710</GBQualifier_value>
          </GBQualifier>
        </GBFeature_quals>
      </GBFeature>
    </GBSeq_feature-table>
  </GBSeq>
</GBSet>)";

		auto metadata = NCBIParser::ParseGenBankXML(xml);

		CHECK(metadata.accession == "NC_001416.1");
		CHECK(metadata.version == 1);
		CHECK(metadata.description == "Escherichia phage Lambda, complete genome");
		CHECK(metadata.organism == "Enterobacteria phage lambda");
		CHECK(metadata.taxonomy_id == 10710);
		CHECK(metadata.length == 48502);
		CHECK(metadata.molecule_type == "DNA");
	}

	SECTION("Empty XML") {
		std::string xml = "";
		auto metadata = NCBIParser::ParseGenBankXML(xml);
		CHECK(metadata.accession.empty());
	}
}

TEST_CASE("NCBIParser Feature Table parsing", "[ncbi]") {
	SECTION("Basic feature table") {
		std::string ft = ">Feature ref|NC_001416.1|\n"
		                 "191\t736\tgene\n"
		                 "\t\t\tgene\tnu1\n"
		                 "\t\t\tlocus_tag\tlambdap01\n"
		                 "191\t736\tCDS\n"
		                 "\t\t\tproduct\tterminase small subunit\n"
		                 "\t\t\tprotein_id\tref|NP_040580.1|\n";

		auto batch = NCBIParser::ParseFeatureTable(ft);

		REQUIRE(batch.size() == 2);

		// First feature: gene
		CHECK(batch.features[0].seqid == "NC_001416.1");
		CHECK(batch.features[0].source == "RefSeq");
		CHECK(batch.features[0].type == "gene");
		CHECK(batch.features[0].position == 191);
		CHECK(batch.features[0].stop_position == 736);
		CHECK(batch.features[0].strand == "+");
		CHECK(batch.features[0].phase == -1); // Not CDS
		REQUIRE(batch.features[0].attrs.size() == 2);
		CHECK(batch.features[0].attrs[0].first == "gene");
		CHECK(batch.features[0].attrs[0].second == "nu1");

		// Second feature: CDS
		CHECK(batch.features[1].type == "CDS");
		CHECK(batch.features[1].phase == 0); // CDS has phase
		REQUIRE(batch.features[1].attrs.size() == 2);
		CHECK(batch.features[1].attrs[0].first == "product");
		CHECK(batch.features[1].attrs[0].second == "terminase small subunit");
	}

	SECTION("Complement strand (reversed positions)") {
		std::string ft = ">Feature NC_001416.1\n"
		                 "5000\t4000\tgene\n"
		                 "\t\t\tgene\tcomplement_gene\n";

		auto batch = NCBIParser::ParseFeatureTable(ft);

		REQUIRE(batch.size() == 1);
		CHECK(batch.features[0].position == 4000); // Swapped to smaller
		CHECK(batch.features[0].stop_position == 5000);
		CHECK(batch.features[0].strand == "-");
	}

	SECTION("Empty feature table") {
		std::string ft = "";
		auto batch = NCBIParser::ParseFeatureTable(ft);
		CHECK(batch.empty());
	}

	SECTION("Feature table with only header") {
		std::string ft = ">Feature ref|NC_001416.1|\n";
		auto batch = NCBIParser::ParseFeatureTable(ft);
		CHECK(batch.empty());
	}

	SECTION("Partial feature indicators") {
		std::string ft = ">Feature NC_001416.1\n"
		                 "<1\t100\tgene\n"
		                 "\t\t\tgene\tpartial5\n"
		                 "200\t>300\tgene\n"
		                 "\t\t\tgene\tpartial3\n";

		auto batch = NCBIParser::ParseFeatureTable(ft);

		REQUIRE(batch.size() == 2);
		CHECK(batch.features[0].position == 1);
		CHECK(batch.features[0].stop_position == 100);
		CHECK(batch.features[1].position == 200);
		CHECK(batch.features[1].stop_position == 300);
	}
}

// Helper function to create a ZIP archive in memory with a single file
static std::string CreateZipInMemory(const std::string &filename, const std::string &content) {
	duckdb_miniz::mz_zip_archive zip;
	memset(&zip, 0, sizeof(zip));

	// Initialize writer for heap allocation
	if (!duckdb_miniz::mz_zip_writer_init_heap(&zip, 0, 0)) {
		throw std::runtime_error("Failed to initialize ZIP writer");
	}

	// Add the file to the archive (use level 6 compression)
	if (!duckdb_miniz::mz_zip_writer_add_mem(&zip, filename.c_str(), content.data(), content.size(),
	                                         duckdb_miniz::MZ_DEFAULT_COMPRESSION)) {
		duckdb_miniz::mz_zip_writer_end(&zip);
		throw std::runtime_error("Failed to add file to ZIP");
	}

	// Finalize and get the buffer
	void *buf = nullptr;
	size_t size = 0;
	if (!duckdb_miniz::mz_zip_writer_finalize_heap_archive(&zip, &buf, &size)) {
		duckdb_miniz::mz_zip_writer_end(&zip);
		throw std::runtime_error("Failed to finalize ZIP archive");
	}

	// Create string from buffer
	std::string result(static_cast<char *>(buf), size);

	// Cleanup
	duckdb_miniz::mz_zip_writer_end(&zip);
	free(buf);

	return result;
}

TEST_CASE("ExtractFromZip", "[ncbi]") {
	SECTION("Single file extraction - no corruption") {
		// Create a test FASTA with multiple sequences (simulating assembly)
		std::string original_fasta = ">seq1 first sequence\n"
		                             "ATCGATCGATCG\n"
		                             ">seq2 second sequence\n"
		                             "GCTAGCTAGCTA\n"
		                             ">seq3 last sequence\n"
		                             "NNNNNNNNNNNN\n";

		// Create ZIP in memory
		std::string zip_data = CreateZipInMemory("test_genomic.fna", original_fasta);

		// Extract using miint::ExtractFromZip
		std::string extracted = ExtractFromZip(zip_data, "_genomic.fna");

		// Verify exact match - no extra null bytes
		CHECK(extracted.size() == original_fasta.size());
		CHECK(extracted == original_fasta);

		// Explicitly check for null bytes at the end (the bug we're fixing)
		CHECK(extracted.find('\0') == std::string::npos);
	}

	SECTION("Extraction preserves exact content size") {
		// Test with a specific size that might cause alignment issues
		std::string content = "ABCDEFGHIJ"; // 10 bytes
		std::string zip_data = CreateZipInMemory("test.txt", content);

		std::string extracted = ExtractFromZip(zip_data, ".txt");

		CHECK(extracted.size() == 10);
		CHECK(extracted == content);
	}

	SECTION("Large content extraction") {
		// Create a larger content to test (1MB)
		std::string large_content;
		large_content.reserve(1024 * 1024);
		for (int i = 0; i < 1024 * 1024; i++) {
			large_content += "ACGT"[i % 4];
		}

		std::string zip_data = CreateZipInMemory("large_genomic.fna", large_content);
		std::string extracted = ExtractFromZip(zip_data, "_genomic.fna");

		CHECK(extracted.size() == large_content.size());
		CHECK(extracted == large_content);
		CHECK(extracted.find('\0') == std::string::npos);
	}

	SECTION("Pattern matching") {
		// Verify pattern matching works
		std::string content = "test content";
		std::string zip_data = CreateZipInMemory("ncbi_dataset/data/GCF_000001/GCF_000001_genomic.fna", content);

		std::string extracted = ExtractFromZip(zip_data, "_genomic.fna");

		CHECK(extracted == content);
	}

	SECTION("No matching file throws") {
		std::string content = "test content";
		std::string zip_data = CreateZipInMemory("other.txt", content);

		// Should throw because no file matches "_genomic.fna"
		CHECK_THROWS_AS(ExtractFromZip(zip_data, "_genomic.fna"), std::runtime_error);
	}

	SECTION("Empty file extraction succeeds") {
		// An empty file in the ZIP should extract successfully (return empty string)
		std::string empty_content = "";
		std::string zip_data = CreateZipInMemory("empty_genomic.fna", empty_content);

		std::string extracted = ExtractFromZip(zip_data, "_genomic.fna");

		CHECK(extracted.empty());
		CHECK(extracted.size() == 0);
	}

	SECTION("Error message includes file count") {
		std::string content = "test content";
		std::string zip_data = CreateZipInMemory("other.txt", content);

		try {
			ExtractFromZip(zip_data, "_genomic.fna");
			FAIL("Expected exception was not thrown");
		} catch (const std::runtime_error &e) {
			std::string msg = e.what();
			// Error should mention the substring and file count
			CHECK(msg.find("_genomic.fna") != std::string::npos);
			CHECK(msg.find("1 files checked") != std::string::npos);
		}
	}
}
