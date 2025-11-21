#include "htslib-1.22.1/htslib/hts_log.h"
#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <string>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include "SAMReader.hpp"

// Test fixture for RAII-based temp file management
class TempFileFixture {
public:
	~TempFileFixture() {
		for (const auto &path : temp_files_) {
			std::filesystem::remove(path);
		}
	}

	void write_temp_sam(const std::string &path, const std::vector<std::string> &records) {
		std::ofstream out(path);
		if (!out) {
			throw std::runtime_error("Failed to create temp file: " + path);
		}
		for (const auto &r : records) {
			out << r;
		}
		out.close();
		temp_files_.push_back(path);
	}

private:
	std::vector<std::string> temp_files_;
};

// Helper functions for generating SAM content
std::string header_reference_line(const std::string &id, const std::string &length) {
	return "@SQ\tSN:" + id + "\tLN:" + length + "\n";
}

std::string unpaired_record(const std::string &id, const std::string &ref, const std::string &flags,
                            const std::string &cigar, const std::string &fields = "") {
	return id + "\t" + flags + "\t" + ref + "\t" + "1\t0\t" + cigar + "\t*\t0\t0\tATGC\tIIII\n";
}

TEST_CASE("SAMReader with header unpaired no flags", "[SAMReader]") {
	TempFileFixture fixture;
	auto path = "test_header.sam";
	fixture.write_temp_sam(path, {header_reference_line("g1", "1000"), unpaired_record("r1", "g1", "0", "4M"),
	                              unpaired_record("r2", "g1", "256", "4M")});

	miint::SAMReader reader(path);
	auto batch = reader.read(5);

	REQUIRE((batch.size() == 2));
	REQUIRE((batch.read_ids[0] == "r1"));
	REQUIRE((batch.flags[0] == 0));
	REQUIRE((batch.references[0] == "g1"));
	REQUIRE((batch.positions[0] == 1));
	REQUIRE((batch.mapqs[0] == 0));
	REQUIRE((batch.cigars[0] == "4M"));
	REQUIRE((batch.mate_references[0] == "*"));
	REQUIRE((batch.mate_positions[0] == 0));
	REQUIRE((batch.template_lengths[0] == 0));
	REQUIRE((batch.tag_as_values[0] == -1));
	REQUIRE((batch.tag_xs_values[0] == -1));
	REQUIRE((batch.tag_ys_values[0] == -1));
	REQUIRE((batch.tag_xn_values[0] == -1));
	REQUIRE((batch.tag_xm_values[0] == -1));
	REQUIRE((batch.tag_xo_values[0] == -1));
	REQUIRE((batch.tag_xg_values[0] == -1));
	REQUIRE((batch.tag_nm_values[0] == -1));
	REQUIRE((batch.tag_yt_values[0] == ""));
	REQUIRE((batch.tag_md_values[0] == ""));
	REQUIRE((batch.tag_sa_values[0] == ""));

	REQUIRE((batch.read_ids[1] == "r2"));
	REQUIRE((batch.flags[1] == 256));
	REQUIRE((batch.references[1] == "g1"));
	REQUIRE((batch.positions[1] == 1));
	REQUIRE((batch.mapqs[1] == 0));
	REQUIRE((batch.cigars[1] == "4M"));
	REQUIRE((batch.mate_references[1] == "*"));
	REQUIRE((batch.mate_positions[1] == 0));
	REQUIRE((batch.template_lengths[1] == 0));
	REQUIRE((batch.tag_as_values[1] == -1));
	REQUIRE((batch.tag_xs_values[1] == -1));
	REQUIRE((batch.tag_ys_values[1] == -1));
	REQUIRE((batch.tag_xn_values[1] == -1));
	REQUIRE((batch.tag_xm_values[1] == -1));
	REQUIRE((batch.tag_xo_values[1] == -1));
	REQUIRE((batch.tag_xg_values[1] == -1));
	REQUIRE((batch.tag_nm_values[1] == -1));
	REQUIRE((batch.tag_yt_values[1] == ""));
	REQUIRE((batch.tag_md_values[1] == ""));
	REQUIRE((batch.tag_sa_values[1] == ""));
}

TEST_CASE("Headerless constructor with single reference", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_headerless_single.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "ref1", "0", "4M")});

	std::unordered_map<std::string, uint64_t> refs = {{"ref1", 1000}};
	miint::SAMReader reader(path, refs);
	auto batch = reader.read(1);

	REQUIRE((batch.size() == 1));
	REQUIRE((batch.read_ids[0] == "r1"));
	REQUIRE((batch.references[0] == "ref1"));
}

TEST_CASE("Headerless constructor with multiple references", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_headerless_multi.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "chr1", "0", "4M"), unpaired_record("r2", "chr2", "0", "4M"),
	                              unpaired_record("r3", "chr3", "0", "4M")});

	std::unordered_map<std::string, uint64_t> refs = {{"chr1", 1000}, {"chr2", 2000}, {"chr3", 3000}};
	miint::SAMReader reader(path, refs);
	auto batch = reader.read(5);

	REQUIRE((batch.size() == 3));
	REQUIRE((batch.references[0] == "chr1"));
	REQUIRE((batch.references[1] == "chr2"));
	REQUIRE((batch.references[2] == "chr3"));
}

TEST_CASE("Headerless constructor reads records correctly", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_headerless_integrity.sam";
	fixture.write_temp_sam(
	    path, {unpaired_record("read1", "genome1", "0", "4M"), unpaired_record("read2", "genome1", "256", "4M")});

	std::unordered_map<std::string, uint64_t> refs = {{"genome1", 5000}};
	miint::SAMReader reader(path, refs);
	auto batch = reader.read(5);

	REQUIRE((batch.size() == 2));
	REQUIRE((batch.read_ids[0] == "read1"));
	REQUIRE((batch.flags[0] == 0));
	REQUIRE((batch.references[0] == "genome1"));
	REQUIRE((batch.cigars[0] == "4M"));

	REQUIRE((batch.read_ids[1] == "read2"));
	REQUIRE((batch.flags[1] == 256));
	REQUIRE((batch.references[1] == "genome1"));
}

TEST_CASE("Headerless constructor with large reference length", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_headerless_large.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "longref", "0", "4M")});

	std::unordered_map<std::string, uint64_t> refs = {{"longref", 18446744073709551615ULL}};
	miint::SAMReader reader(path, refs);
	auto batch = reader.read(1);

	REQUIRE((batch.size() == 1));
	REQUIRE((batch.references[0] == "longref"));
}

// Note: We cannot reliably detect if a SAM file has a header without consuming the file position,
// which breaks stdin support. Users must use the correct constructor for their file type.
// The header-based constructor will fail appropriately if the file lacks a header.
//
// TEST_CASE("Headerless constructor throws on file with header", "[SAMReader][headerless]") {
// 	TempFileFixture fixture;
// 	auto path = "test_has_header.sam";
// 	fixture.write_temp_sam(path, {header_reference_line("ref1", "1000"), unpaired_record("r1", "ref1", "0", "4M")});
//
// 	std::unordered_map<std::string, uint64_t> refs = {{"ref1", 1000}};
// 	REQUIRE_THROWS_WITH(miint::SAMReader(path, refs), "File contains header but headerless mode was requested");
// }

TEST_CASE("Headerless constructor throws on empty reference map", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_empty_refs.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "ref1", "0", "4M")});

	std::unordered_map<std::string, uint64_t> refs;
	REQUIRE_THROWS_WITH(miint::SAMReader(path, refs), "Reference map cannot be empty");
}

TEST_CASE("Headerless constructor throws on empty reference name", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_empty_name.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "ref1", "0", "4M")});

	std::unordered_map<std::string, uint64_t> refs = {{"", 1000}};
	REQUIRE_THROWS_WITH(miint::SAMReader(path, refs), "Reference name cannot be empty");
}

TEST_CASE("Headerless constructor throws on zero length reference", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_zero_length.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "ref1", "0", "4M")});

	std::unordered_map<std::string, uint64_t> refs = {{"ref1", 0}};
	REQUIRE_THROWS_WITH(miint::SAMReader(path, refs), "Reference length cannot be zero");
}

TEST_CASE("Headerless constructor throws on nonexistent file", "[SAMReader][headerless]") {
	std::unordered_map<std::string, uint64_t> refs = {{"ref1", 1000}};
	hts_set_log_level(HTS_LOG_OFF);
	REQUIRE_THROWS_AS(miint::SAMReader("nonexistent_file.sam", refs), std::runtime_error);
}

TEST_CASE("Header constructor throws on headerless file", "[SAMReader][header]") {
	TempFileFixture fixture;
	auto path = "test_no_header.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "ref1", "0", "4M")});

	REQUIRE_THROWS_WITH(miint::SAMReader(path), "SAM file missing required header");
}

TEST_CASE("Header constructor throws on nonexistent file", "[SAMReader][header]") {
	hts_set_log_level(HTS_LOG_OFF);
	REQUIRE_THROWS_AS(miint::SAMReader("nonexistent_file.sam"), std::runtime_error);
}

TEST_CASE("Headerless constructor with unknown reference in data", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_unknown_ref.sam";
	hts_set_log_level(HTS_LOG_OFF);
	fixture.write_temp_sam(path, {unpaired_record("r1", "unknown_ref", "0", "4M")});

	std::unordered_map<std::string, uint64_t> refs = {{"ref1", 1000}};
	miint::SAMReader reader(path, refs);

	// LIMITATION: htslib automatically sets the unmapped flag (0x4) when it can't find
	// a reference in the header. This means we cannot distinguish between genuinely
	// unmapped reads and reads with unknown references. The record will appear as
	// unmapped with reference="*".
	// To catch unknown references, users must ensure their reference_lengths table
	// includes all references present in the data files.
	auto batch = reader.read(1);
	REQUIRE((batch.size() == 1));
	REQUIRE((batch.references[0] == "*")); // Treated as unmapped by htslib
	REQUIRE((batch.flags[0] & 0x4) != 0);  // Unmapped flag is set by htslib
}

TEST_CASE("Headerless constructor with empty file", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_empty.sam";
	fixture.write_temp_sam(path, {});

	std::unordered_map<std::string, uint64_t> refs = {{"ref1", 1000}};
	miint::SAMReader reader(path, refs);
	auto batch = reader.read(1);
	REQUIRE((batch.size() == 0));
}

TEST_CASE("Headerless constructor with whitespace only file", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_whitespace.sam";
	fixture.write_temp_sam(path, {"\n\n  \n\t\n"});

	std::unordered_map<std::string, uint64_t> refs = {{"ref1", 1000}};
	miint::SAMReader reader(path, refs);
	auto batch = reader.read(1);
	REQUIRE((batch.size() == 0));
}

TEST_CASE("Headerless constructor multiple sequential reads", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_sequential.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "ref1", "0", "4M"), unpaired_record("r2", "ref1", "0", "4M"),
	                              unpaired_record("r3", "ref1", "0", "4M"), unpaired_record("r4", "ref1", "0", "4M")});

	std::unordered_map<std::string, uint64_t> refs = {{"ref1", 1000}};
	miint::SAMReader reader(path, refs);

	auto batch1 = reader.read(2);
	REQUIRE((batch1.size() == 2));
	REQUIRE((batch1.read_ids[0] == "r1"));
	REQUIRE((batch1.read_ids[1] == "r2"));

	auto batch2 = reader.read(2);
	REQUIRE((batch2.size() == 2));
	REQUIRE((batch2.read_ids[0] == "r3"));
	REQUIRE((batch2.read_ids[1] == "r4"));

	auto batch3 = reader.read(2);
	REQUIRE((batch3.size() == 0));
}

TEST_CASE("Header constructor multiple sequential reads", "[SAMReader][header]") {
	TempFileFixture fixture;
	auto path = "test_header_sequential.sam";
	fixture.write_temp_sam(path, {header_reference_line("ref1", "1000"), unpaired_record("r1", "ref1", "0", "4M"),
	                              unpaired_record("r2", "ref1", "0", "4M"), unpaired_record("r3", "ref1", "0", "4M")});

	miint::SAMReader reader(path);

	auto batch1 = reader.read(2);
	REQUIRE((batch1.size() == 2));
	REQUIRE((batch1.read_ids[0] == "r1"));
	REQUIRE((batch1.read_ids[1] == "r2"));

	auto batch2 = reader.read(2);
	REQUIRE((batch2.size() == 1));
	REQUIRE((batch2.read_ids[0] == "r3"));
}

TEST_CASE("Headerless constructor throws on reference name with tab", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_tab_in_name.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "ref1", "0", "4M")});

	std::unordered_map<std::string, uint64_t> refs = {{"ref\t1", 1000}};
	REQUIRE_THROWS_WITH(miint::SAMReader(path, refs), "Reference name contains invalid characters");
}

TEST_CASE("Headerless constructor throws on reference name with newline", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_newline_in_name.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "ref1", "0", "4M")});

	std::unordered_map<std::string, uint64_t> refs = {{"ref\n1", 1000}};
	REQUIRE_THROWS_WITH(miint::SAMReader(path, refs), "Reference name contains invalid characters");
}

TEST_CASE("Headerless constructor with reference map superset", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_superset.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "ref1", "0", "4M")});

	std::unordered_map<std::string, uint64_t> refs = {{"ref1", 1000}, {"ref2", 2000}, {"ref3", 3000}};
	miint::SAMReader reader(path, refs);
	auto batch = reader.read(1);

	REQUIRE((batch.size() == 1));
	REQUIRE((batch.references[0] == "ref1"));
}

TEST_CASE("Headerless constructor with existing test data", "[SAMReader][headerless]") {
	std::unordered_map<std::string, uint64_t> refs = {{"G1234", 20}, {"G000144735", 90}};
	miint::SAMReader reader("data/sam/foo_no_header.sam", refs);
	auto batch = reader.read(10);

	REQUIRE((batch.size() == 4));
	REQUIRE((batch.read_ids[0] == "foo-1"));
	REQUIRE((batch.references[0] == "G1234"));
	REQUIRE((batch.read_ids[2] == "foo-3"));
	REQUIRE((batch.references[2] == "G000144735"));
}

TEST_CASE("Header constructor with existing test data", "[SAMReader][header]") {
	miint::SAMReader reader("data/sam/foo_has_header.sam");
	auto batch = reader.read(10);

	REQUIRE((batch.size() == 4));
	REQUIRE((batch.read_ids[0] == "foo-1"));
	REQUIRE((batch.references[0] == "G1234"));
	REQUIRE((batch.read_ids[2] == "foo-3"));
	REQUIRE((batch.references[2] == "G000144735"));
}

TEST_CASE("Headerless constructor throws on reference name starting with asterisk", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_asterisk_start.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "ref1", "0", "4M")});

	std::unordered_map<std::string, uint64_t> refs = {{"*ref1", 1000}};
	REQUIRE_THROWS_WITH(miint::SAMReader(path, refs), "Reference name cannot start with '*' or '='");
}

TEST_CASE("Headerless constructor throws on reference name starting with equals", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_equals_start.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "ref1", "0", "4M")});

	std::unordered_map<std::string, uint64_t> refs = {{"=ref1", 1000}};
	REQUIRE_THROWS_WITH(miint::SAMReader(path, refs), "Reference name cannot start with '*' or '='");
}

TEST_CASE("Headerless constructor throws on reference name exceeding max length", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_long_name.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "ref1", "0", "4M")});

	std::string long_name(1025, 'A');
	std::unordered_map<std::string, uint64_t> refs = {{long_name, 1000}};
	REQUIRE_THROWS_WITH(miint::SAMReader(path, refs), "Reference name exceeds maximum length of 1024 characters");
}

TEST_CASE("Headerless constructor throws on reference name with position-like pattern", "[SAMReader][headerless]") {
	TempFileFixture fixture;
	auto path = "test_position_pattern.sam";
	fixture.write_temp_sam(path, {unpaired_record("r1", "ref1", "0", "4M")});

	SECTION("Pattern with single position") {
		std::unordered_map<std::string, uint64_t> refs = {{"HLA-A:01", 1000}};
		REQUIRE_THROWS_WITH(miint::SAMReader(path, refs),
		                    "Reference name ends with position-like pattern (:<digits> or :<digits>-<digits>)");
	}

	SECTION("Pattern with range") {
		std::unordered_map<std::string, uint64_t> refs = {{"chr1:100-200", 1000}};
		REQUIRE_THROWS_WITH(miint::SAMReader(path, refs),
		                    "Reference name ends with position-like pattern (:<digits> or :<digits>-<digits>)");
	}

	SECTION("Pattern with multiple colons") {
		std::unordered_map<std::string, uint64_t> refs = {{"HLA-A*01:01:01:01", 1000}};
		REQUIRE_THROWS_WITH(miint::SAMReader(path, refs),
		                    "Reference name ends with position-like pattern (:<digits> or :<digits>-<digits>)");
	}

	SECTION("Valid name with colon not at end") {
		std::unordered_map<std::string, uint64_t> refs = {{"prefix:123suffix", 1000}};
		miint::SAMReader reader(path, refs);
		REQUIRE_NOTHROW(reader.read(1));
	}
}
