#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <string>
#include <fstream>
#include <filesystem>
#include "QualScore.hpp"
#include "SequenceRecord.hpp"
#include "SequenceReader.hpp"

// Utilities for test data
std::string simple_read(const std::string &id, const std::string &seq, const std::string &qual,
                        const std::string &comment = "") {
	return "@" + id + (comment.empty() ? "" : " " + comment) + "\n" + seq + "\n+\n" + qual + "\n";
}

void write_temp_file(const std::string &path, const std::vector<std::string> &records) {
	std::ofstream out(path);
	for (const auto &r : records) {
		out << r;
	}
}

TEST_CASE("SequenceReader single-end", "[SequenceReader]") {
	auto path = "test_R1.fastq";
	write_temp_file(path, {simple_read("r1", "ACGT", "IIII"), simple_read("r2", "TGCA", "HHHH")});

	miint::SequenceReader reader(path);
	auto batch = reader.read(5);

	REQUIRE((batch.size() == 2));
	REQUIRE((batch[0].read_id == "r1"));
	REQUIRE((batch[1].read_id == "r2"));
	std::filesystem::remove(path);
}

TEST_CASE("SequenceReader paired-end valid", "[SequenceReader]") {
	auto r1 = "test_p1.fastq";
	auto r2 = "test_p2.fastq";

	write_temp_file(r1, {simple_read("x/1", "ACGT", "IIII"), simple_read("y", "TGCA", "HHHH")});
	write_temp_file(r2, {simple_read("x/2", "AAAA", "DDDD"), simple_read("y", "CCCC", "EEEE")});

	miint::SequenceReader reader(r1, r2);
	auto batch = reader.read(2);
	REQUIRE((batch.size() == 2));
	REQUIRE((batch[0].read_id == "x"));
	REQUIRE((batch[0].read2.value() == "AAAA"));
	REQUIRE((batch[1].qual2->as_string() == "EEEE"));

	std::filesystem::remove(r1);
	std::filesystem::remove(r2);
}

TEST_CASE("SequenceReader mismatched IDs", "[SequenceReader]") {
	auto r1 = "bad_p1.fastq";
	auto r2 = "bad_p2.fastq";

	write_temp_file(r1, {simple_read("r1/1", "ACGT", "IIII")});
	write_temp_file(r2, {simple_read("r2/2", "ACGT", "IIII")});

	REQUIRE_THROWS_WITH(miint::SequenceReader(r1, r2).read(1), "Mismatched read IDs: r1/1 vs r2/2");

	std::filesystem::remove(r1);
	std::filesystem::remove(r2);
}

TEST_CASE("SequenceReader partial batch", "[SequenceReader]") {
	auto path = "few_reads.fastq";
	write_temp_file(path, {simple_read("only_one", "ACGT", "IIII")});

	miint::SequenceReader reader(path);
	auto batch = reader.read(5);
	REQUIRE((batch.size() == 1));
	std::filesystem::remove(path);
}

TEST_CASE("QualScore from string and back", "[QualScore]") {
	std::string qstr = "IJKLMNOP"; // ASCII 73+
	miint::QualScore qs(qstr);
	REQUIRE((qs.as_string() == qstr));

	auto vec = qs.as_vec();
	REQUIRE((vec.size() == qstr.size()));
	for (size_t i = 0; i < vec.size(); ++i) {
		REQUIRE((static_cast<char>(vec[i] + 33) == qstr[i]));
	}

	miint::QualScore from_vec(vec);
	REQUIRE((from_vec.as_string() == qstr));
}

TEST_CASE("SequenceRecord single-end", "[SequenceRecord]") {
	miint::SequenceRecord rec("read1", "abc", "ATGC", "IIII");
	REQUIRE((rec.read_id == "read1"));
	REQUIRE((rec.comment == "abc"));
	REQUIRE((rec.read1 == "ATGC"));
	REQUIRE_FALSE(rec.read2.has_value());
	REQUIRE((rec.qual1.as_string() == "IIII"));
	REQUIRE_FALSE(rec.qual2.has_value());
}

TEST_CASE("SequenceRecord paired-end", "[SequenceRecord]") {
	miint::SequenceRecord rec("read1", "xyz", "ACGT", "HHHH", "TGCA", "BBBB");
	REQUIRE((rec.read2.value() == "TGCA"));
	REQUIRE((rec.qual2->as_string() == "BBBB"));
}
