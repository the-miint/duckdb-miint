#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <string>
#include <fstream>
#include <filesystem>
#include "QualScore.hpp"
#include "SequenceRecord.hpp"
#include "SequenceReader.hpp"

// Test fixture for RAII-based temp file management
class TempFileFixture {
public:
	~TempFileFixture() {
		for (const auto &path : temp_files_) {
			std::filesystem::remove(path);
		}
	}

	void write_temp_fastq(const std::string &path, const std::vector<std::string> &records) {
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

	std::string simple_read(const std::string &id, const std::string &seq, const std::string &qual,
	                        const std::string &comment = "") {
		return "@" + id + (comment.empty() ? "" : " " + comment) + "\n" + seq + "\n+\n" + qual + "\n";
	}

	std::string simple_fasta(const std::string &id, const std::string &seq, const std::string &comment = "") {
		return ">" + id + (comment.empty() ? "" : " " + comment) + "\n" + seq + "\n";
	}

private:
	std::vector<std::string> temp_files_;
};

TEST_CASE("SequenceReader single-end", "[SequenceReader]") {
	TempFileFixture fixture;
	auto path = "test_R1.fastq";
	fixture.write_temp_fastq(path,
	                         {fixture.simple_read("r1", "ACGT", "IIII"), fixture.simple_read("r2", "TGCA", "HHHH")});

	miint::SequenceReader reader(path);
	auto batch = reader.read(5);

	REQUIRE((batch.size() == 2));
	REQUIRE((batch[0].read_id == "r1"));
	REQUIRE((batch[1].read_id == "r2"));
}

TEST_CASE("SequenceReader paired-end valid", "[SequenceReader]") {
	TempFileFixture fixture;
	auto r1 = "test_p1.fastq";
	auto r2 = "test_p2.fastq";

	fixture.write_temp_fastq(r1,
	                         {fixture.simple_read("x/1", "ACGT", "IIII"), fixture.simple_read("y", "TGCA", "HHHH")});
	fixture.write_temp_fastq(r2,
	                         {fixture.simple_read("x/2", "AAAA", "DDDD"), fixture.simple_read("y", "CCCC", "EEEE")});

	miint::SequenceReader reader(r1, r2);
	auto batch = reader.read(2);
	REQUIRE((batch.size() == 2));
	REQUIRE((batch[0].read_id == "x"));
	REQUIRE((batch[0].read2.value() == "AAAA"));
	REQUIRE((batch[1].qual2->as_string() == "EEEE"));
}

TEST_CASE("SequenceReader mismatched IDs", "[SequenceReader]") {
	TempFileFixture fixture;
	auto r1 = "bad_p1.fastq";
	auto r2 = "bad_p2.fastq";

	fixture.write_temp_fastq(r1, {fixture.simple_read("r1/1", "ACGT", "IIII")});
	fixture.write_temp_fastq(r2, {fixture.simple_read("r2/2", "ACGT", "IIII")});

	REQUIRE_THROWS_WITH(miint::SequenceReader(r1, r2).read(1), "Mismatched read IDs: r1/1 vs r2/2");
}

TEST_CASE("SequenceReader partial batch", "[SequenceReader]") {
	TempFileFixture fixture;
	auto path = "few_reads.fastq";
	fixture.write_temp_fastq(path, {fixture.simple_read("only_one", "ACGT", "IIII")});

	miint::SequenceReader reader(path);
	auto batch = reader.read(5);
	REQUIRE((batch.size() == 1));
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

TEST_CASE("SequenceReader FASTA format single-end", "[SequenceReader][FASTA]") {
	TempFileFixture fixture;
	auto path = "test_fasta.fa";
	fixture.write_temp_fastq(path,
	                         {fixture.simple_fasta("seq1", "ATGC", "comment1"), fixture.simple_fasta("seq2", "GGCC")});

	miint::SequenceReader reader(path);
	auto batch = reader.read(5);

	REQUIRE((batch.size() == 2));
	REQUIRE((batch[0].read_id == "seq1"));
	REQUIRE((batch[0].comment == "comment1"));
	REQUIRE((batch[0].read1 == "ATGC"));
	REQUIRE((batch[0].qual1.as_string().empty()));
	REQUIRE((batch[1].read_id == "seq2"));
	REQUIRE((batch[1].comment == ""));
	REQUIRE((batch[1].qual1.as_string().empty()));
}

TEST_CASE("SequenceReader FASTA format paired-end", "[SequenceReader][FASTA]") {
	TempFileFixture fixture;
	auto r1 = "test_r1.fa";
	auto r2 = "test_r2.fa";

	fixture.write_temp_fastq(r1, {fixture.simple_fasta("read1/1", "ATGC"), fixture.simple_fasta("read2", "GGCC")});
	fixture.write_temp_fastq(r2, {fixture.simple_fasta("read1/2", "TGCA"), fixture.simple_fasta("read2", "CCGG")});

	miint::SequenceReader reader(r1, r2);
	auto batch = reader.read(5);

	REQUIRE((batch.size() == 2));
	REQUIRE((batch[0].read_id == "read1"));
	REQUIRE((batch[0].read1 == "ATGC"));
	REQUIRE((batch[0].read2.value() == "TGCA"));
	REQUIRE((batch[0].qual1.as_string().empty()));
	REQUIRE((batch[0].qual2->as_string().empty()));
}

TEST_CASE("SequenceReader mixed FASTA/FASTQ paired-end throws", "[SequenceReader][FASTA][FASTQ][error]") {
	TempFileFixture fixture;
	auto r1 = "test_mixed_r1.fa";
	auto r2 = "test_mixed_r2.fq";

	fixture.write_temp_fastq(r1, {fixture.simple_fasta("read1/1", "ATGC"), fixture.simple_fasta("read2", "GGCC")});
	fixture.write_temp_fastq(
	    r2, {fixture.simple_read("read1/2", "TGCA", "IIII"), fixture.simple_read("read2", "CCGG", "HHHH")});

	REQUIRE_THROWS_WITH(miint::SequenceReader(r1, r2),
	                    Catch::Matchers::ContainsSubstring("Cannot mix FASTA and FASTQ"));
}

TEST_CASE("SequenceReader empty file throws", "[SequenceReader][error]") {
	TempFileFixture fixture;
	auto path = "empty.fq";
	fixture.write_temp_fastq(path, {});

	REQUIRE_THROWS_WITH(miint::SequenceReader(path), Catch::Matchers::ContainsSubstring("Empty file"));
}

TEST_CASE("SequenceReader gzipped file", "[SequenceReader][compression]") {
	miint::SequenceReader reader("data/fastq/foo.r1.fastq.gz");
	auto batch = reader.read(5);

	REQUIRE((batch.size() == 2));
	REQUIRE((batch[0].read_id == "foo1"));
	REQUIRE((batch[0].read1 == "ATGC"));
}

TEST_CASE("SequenceReader multiple sequential exhaustive reads", "[SequenceReader]") {
	TempFileFixture fixture;
	auto path = "multi_batch.fq";
	fixture.write_temp_fastq(path,
	                         {fixture.simple_read("r1", "AAAA", "IIII"), fixture.simple_read("r2", "TTTT", "HHHH"),
	                          fixture.simple_read("r3", "GGGG", "GGGG"), fixture.simple_read("r4", "CCCC", "FFFF")});

	miint::SequenceReader reader(path);

	auto batch1 = reader.read(2);
	REQUIRE((batch1.size() == 2));
	REQUIRE((batch1[0].read_id == "r1"));
	REQUIRE((batch1[1].read_id == "r2"));

	auto batch2 = reader.read(2);
	REQUIRE((batch2.size() == 2));
	REQUIRE((batch2[0].read_id == "r3"));
	REQUIRE((batch2[1].read_id == "r4"));

	auto batch3 = reader.read(2);
	REQUIRE((batch3.size() == 0));
}

TEST_CASE("SequenceReader large batch size", "[SequenceReader]") {
	TempFileFixture fixture;
	auto path = "large_batch.fq";
	std::vector<std::string> records;
	records.reserve(10000);
	for (int i = 0; i < 10000; i++) {
		records.push_back(fixture.simple_read("read" + std::to_string(i), "ACGT", "IIII"));
	}
	fixture.write_temp_fastq(path, records);

	miint::SequenceReader reader(path);
	auto batch = reader.read(10000);

	REQUIRE((batch.size() == 10000));
	REQUIRE((batch[0].read_id == "read0"));
	REQUIRE((batch[9999].read_id == "read9999"));
}

TEST_CASE("SequenceReader partial read at EOF", "[SequenceReader]") {
	TempFileFixture fixture;
	auto path = "partial_eof.fq";
	fixture.write_temp_fastq(path,
	                         {fixture.simple_read("r1", "ACGT", "IIII"), fixture.simple_read("r2", "TGCA", "HHHH"),
	                          fixture.simple_read("r3", "GGCC", "GGGG")});

	miint::SequenceReader reader(path);
	auto batch = reader.read(10);

	REQUIRE((batch.size() == 3));
}

TEST_CASE("QualScore edge cases", "[QualScore]") {
	std::string low_qual(10, '!');
	miint::QualScore qs_low(low_qual);
	auto vec_low = qs_low.as_vec();
	for (auto q : vec_low) {
		REQUIRE((q == 0));
	}

	std::string high_qual(10, '~');
	miint::QualScore qs_high(high_qual);
	auto vec_high = qs_high.as_vec();
	for (auto q : vec_high) {
		REQUIRE((q == 93));
	}
}

TEST_CASE("SequenceReader paired-end mismatched count with context", "[SequenceReader][error]") {
	TempFileFixture fixture;
	auto r1 = "mismatch_r1.fq";
	auto r2 = "mismatch_r2.fq";

	fixture.write_temp_fastq(r1,
	                         {fixture.simple_read("r1", "ACGT", "IIII"), fixture.simple_read("r2", "TGCA", "HHHH")});
	fixture.write_temp_fastq(r2, {fixture.simple_read("r1", "AAAA", "DDDD")});

	miint::SequenceReader reader(r1, r2);
	REQUIRE_THROWS_WITH(reader.read(5), Catch::Matchers::ContainsSubstring("missing mate"));
}
