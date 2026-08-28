#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <string>
#include <fstream>
#include <filesystem>
#include "QualScore.hpp"
#include "SequenceRecord.hpp"
#include "SequenceReader.hpp"
#include <vector>

// getrlimit/geteuid and permission bits that actually deny the owner are POSIX-only; the
// open-failure test below is compiled out elsewhere. Same guard as test_NewickParser.cpp.
#if !defined(_WIN32) && !defined(__EMSCRIPTEN__)
#include <cerrno>
#include <fcntl.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

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
	REQUIRE((batch.read_ids[0] == "r1"));
	REQUIRE((batch.read_ids[1] == "r2"));
}

TEST_CASE("SequenceReader single-end / bug", "[SequenceReader]") {
	TempFileFixture fixture;
	auto path = "test_R1.fastq";
	fixture.write_temp_fastq(
	    path, {fixture.simple_read("r1/foo/bar", "ACGT", "IIII"), fixture.simple_read("r2/123/456", "TGCA", "HHHH")});

	miint::SequenceReader reader(path);
	auto batch = reader.read(5);

	REQUIRE((batch.size() == 2));
	REQUIRE((batch.read_ids[0] == "r1/foo/bar"));
	REQUIRE((batch.read_ids[1] == "r2/123/456"));
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
	REQUIRE((batch.read_ids[0] == "x"));
	REQUIRE((batch.sequences2[0] == "AAAA"));
	REQUIRE((batch.quals2[1].as_string() == "EEEE"));
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
	REQUIRE((batch.read_ids[0] == "seq1"));
	REQUIRE((batch.comments[0] == "comment1"));
	REQUIRE((batch.sequences1[0] == "ATGC"));
	REQUIRE((batch.quals1[0].as_string().empty()));
	REQUIRE((batch.read_ids[1] == "seq2"));
	REQUIRE((batch.comments[1] == ""));
	REQUIRE((batch.quals1[1].as_string().empty()));
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
	REQUIRE((batch.read_ids[0] == "read1"));
	REQUIRE((batch.sequences1[0] == "ATGC"));
	REQUIRE((batch.sequences2[0] == "TGCA"));
	REQUIRE((batch.quals1[0].as_string().empty()));
	REQUIRE((batch.quals2[0].as_string().empty()));
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

TEST_CASE("SequenceReader mixed FASTA/FASTQ paired-end allowed when flag set",
          "[SequenceReader][FASTA][FASTQ][paired-end]") {
	TempFileFixture fixture;
	auto r1 = "test_allowed_r1.fa";
	auto r2 = "test_allowed_r2.fq";

	fixture.write_temp_fastq(r1, {fixture.simple_fasta("read1/1", "ATGC"), fixture.simple_fasta("read2", "GGCC")});
	fixture.write_temp_fastq(
	    r2, {fixture.simple_read("read1/2", "TGCA", "IIII"), fixture.simple_read("read2", "CCGG", "HHHH")});

	miint::SequenceReader reader(r1, r2, true);
	auto batch = reader.read(5);

	REQUIRE((batch.size() == 2));
	REQUIRE((batch.read_ids[0] == "read1"));
	REQUIRE((batch.sequences1[0] == "ATGC"));
	REQUIRE((batch.sequences2[0] == "TGCA"));
	REQUIRE((batch.quals1[0].as_string().empty()));   // FASTA: no quality
	REQUIRE((batch.quals2[0].as_string() == "IIII")); // FASTQ: has quality
}

TEST_CASE("SequenceReader mixed FASTA/FASTQ with longer sequences", "[SequenceReader][FASTA][FASTQ][paired-end]") {
	TempFileFixture fixture;
	auto r1 = "test_long_r1.fa";
	auto r2 = "test_long_r2.fq";

	fixture.write_temp_fastq(
	    r1, {fixture.simple_fasta("read1/1", "ATGCATGC"), fixture.simple_fasta("read2/1", "GGCCGGCC")});
	fixture.write_temp_fastq(r2, {fixture.simple_read("read1/2", "TGCATGCA", "IIIIIIII"),
	                              fixture.simple_read("read2/2", "CCGGCCGG", "HHHHHHHH")});

	miint::SequenceReader reader(r1, r2, true);
	auto batch = reader.read(5);

	REQUIRE((batch.size() == 2));
	REQUIRE((batch.read_ids[0] == "read1"));
	REQUIRE((batch.read_ids[1] == "read2"));
	REQUIRE((batch.sequences1[0] == "ATGCATGC"));
	REQUIRE((batch.sequences2[0] == "TGCATGCA"));
	REQUIRE((batch.quals1[0].as_string().empty()));
	REQUIRE((batch.quals2[0].as_string() == "IIIIIIII"));
}

TEST_CASE("SequenceReader empty file throws", "[SequenceReader][error]") {
	TempFileFixture fixture;
	auto path = "empty.fq";
	fixture.write_temp_fastq(path, {});

	REQUIRE_THROWS_WITH(miint::SequenceReader(path), Catch::Matchers::ContainsSubstring("Empty file"));
}

#if !defined(_WIN32) && !defined(__EMSCRIPTEN__)
// kseq++'s SeqStreamIn does not check gzopen()'s return value. A NULL handle reads as a
// stream with zero records, so before CheckedSeqStreamIn every failed open was reported as
// "Empty file: <path>" -- naming a file that is present, non-empty and perfectly valid.
// That is what made descriptor exhaustion over a large VARCHAR[] scan (issue #260)
// undiagnosable: the error blamed the data. These two conditions must stay distinguishable,
// because the message is the only signal a user gets about which one actually happened.
namespace {
// Returns the message SequenceReader fails with, or "" if it does not throw.
std::string open_failure_message(const std::string &path) {
	try {
		miint::SequenceReader reader(path);
	} catch (const std::exception &e) {
		return e.what();
	}
	return "";
}

// Lowers RLIMIT_NOFILE for its lifetime and always puts the soft limit back, including on
// the stack unwind of a failed REQUIRE.
class FdLimit {
public:
	explicit FdLimit(rlim_t soft) {
		if (getrlimit(RLIMIT_NOFILE, &saved_) != 0) {
			throw std::runtime_error("getrlimit failed");
		}
		struct rlimit lowered = saved_;
		lowered.rlim_cur = soft;
		if (setrlimit(RLIMIT_NOFILE, &lowered) != 0) {
			throw std::runtime_error("setrlimit failed");
		}
	}
	~FdLimit() {
		setrlimit(RLIMIT_NOFILE, &saved_);
	}
	FdLimit(const FdLimit &) = delete;
	FdLimit &operator=(const FdLimit &) = delete;

private:
	struct rlimit saved_;
};

// Opens `path` until the process runs out of descriptors, recording the errno that stopped
// it, and closes them all again on destruction.
class HeldFds {
public:
	explicit HeldFds(const char *path) {
		for (int fd = open(path, O_RDONLY); fd >= 0; fd = open(path, O_RDONLY)) {
			fds_.push_back(fd);
		}
		exhausted_errno = errno;
	}
	~HeldFds() {
		for (int fd : fds_) {
			close(fd);
		}
	}
	HeldFds(const HeldFds &) = delete;
	HeldFds &operator=(const HeldFds &) = delete;

	int exhausted_errno = 0;

private:
	std::vector<int> fds_;
};
} // namespace

TEST_CASE("SequenceReader reports a failed open as a failed open, not an empty file", "[SequenceReader][error]") {
	TempFileFixture fixture;

	SECTION("unreadable file") {
		auto path = "unreadable.fq";
		fixture.write_temp_fastq(path, {fixture.simple_read("read1", "ATGC", "IIII")});
		// Root ignores permission bits, so this case cannot be exercised as root.
		if (geteuid() == 0) {
			SUCCEED("skipped: running as root, permission bits are not enforced");
			return;
		}
		std::filesystem::permissions(path, std::filesystem::perms::none);

		const std::string msg = open_failure_message(path);
		std::filesystem::permissions(path, std::filesystem::perms::owner_all); // let the fixture clean up

		INFO("error was: " << msg);
		REQUIRE(msg.find("Failed to open file") != std::string::npos);
		REQUIRE(msg.find(path) != std::string::npos);
		REQUIRE(msg.find("Empty file") == std::string::npos);
	}

	SECTION("descriptor exhaustion") {
		auto path = "fd_exhausted.fq";
		fixture.write_temp_fastq(path, {fixture.simple_read("read1", "ATGC", "IIII")});

		// Hold every remaining descriptor so the next gzopen() has to fail with EMFILE -- the
		// exact condition a large-VARCHAR[] scan hits once readers are not released. Both the
		// lowered limit and the held descriptors unwind via RAII: a REQUIRE that fails (or one
		// a later edit adds) inside this scope would otherwise leak a 96-descriptor limit into
		// every test that runs after it in this binary, turning one failure into a cascade.
		FdLimit limit(96); // low enough to fill quickly, high enough for Catch2 itself
		HeldFds held(path);

		REQUIRE(held.exhausted_errno == EMFILE);
		const std::string msg = open_failure_message(path);
		INFO("error was: " << msg);
		REQUIRE(msg.find("Failed to open file") != std::string::npos);
		REQUIRE(msg.find("Empty file") == std::string::npos);
	}

	SECTION("a genuinely empty file still reports emptiness") {
		auto path = "really_empty.fq";
		fixture.write_temp_fastq(path, {});

		const std::string msg = open_failure_message(path);
		INFO("error was: " << msg);
		REQUIRE(msg.find("Empty file") != std::string::npos);
		REQUIRE(msg.find("Failed to open file") == std::string::npos);
	}
}
#endif // !_WIN32 && !__EMSCRIPTEN__

TEST_CASE("SequenceReader gzipped file", "[SequenceReader][compression]") {
	miint::SequenceReader reader("data/fastq/foo.r1.fastq.gz");
	auto batch = reader.read(5);

	REQUIRE((batch.size() == 2));
	REQUIRE((batch.read_ids[0] == "foo1"));
	REQUIRE((batch.sequences1[0] == "ATGC"));
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
	REQUIRE((batch1.read_ids[0] == "r1"));
	REQUIRE((batch1.read_ids[1] == "r2"));

	auto batch2 = reader.read(2);
	REQUIRE((batch2.size() == 2));
	REQUIRE((batch2.read_ids[0] == "r3"));
	REQUIRE((batch2.read_ids[1] == "r4"));

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
	REQUIRE((batch.read_ids[0] == "read0"));
	REQUIRE((batch.read_ids[9999] == "read9999"));
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

TEST_CASE("SequenceReader byte budget SE: generous budget returns all records", "[SequenceReader][byte_budget]") {
	TempFileFixture fixture;
	auto path = "budget_generous.fq";
	std::vector<std::string> records;
	for (int i = 0; i < 10; i++) {
		records.push_back(fixture.simple_read("r" + std::to_string(i), "ACGTACGT", "IIIIIIII"));
	}
	fixture.write_temp_fastq(path, records);

	miint::SequenceReader reader(path);
	auto batch = reader.read(100, 1024 * 1024); // 1 MiB -- way more than needed

	REQUIRE((batch.size() == 10));
}

TEST_CASE("SequenceReader byte budget SE: tight budget caps records per chunk", "[SequenceReader][byte_budget]") {
	TempFileFixture fixture;
	auto path = "budget_tight.fq";
	std::vector<std::string> records;
	for (int i = 0; i < 10; i++) {
		// Each record: name=~3 chars, comment=0, seq=8, qual=8 -> ~19 bytes/record
		records.push_back(fixture.simple_read("r" + std::to_string(i), "ACGTACGT", "IIIIIIII"));
	}
	fixture.write_temp_fastq(path, records);

	miint::SequenceReader reader(path);

	// Budget ~25 bytes -> first record (~19 B) fits, adding a second (~38 B total) exceeds.
	// Expect 1 record per call.
	size_t total = 0;
	for (int call = 0; call < 20; call++) {
		auto batch = reader.read(100, 25);
		if (batch.empty()) {
			break;
		}
		REQUIRE((batch.size() == 1));
		total += batch.size();
	}
	REQUIRE((total == 10));
}

TEST_CASE("SequenceReader byte budget SE: starvation guard returns at least one record",
          "[SequenceReader][byte_budget]") {
	TempFileFixture fixture;
	auto path = "budget_starve.fq";
	fixture.write_temp_fastq(path, {fixture.simple_read("big", std::string(1000, 'A'), std::string(1000, 'I'))});

	miint::SequenceReader reader(path);
	auto batch = reader.read(100, 1); // budget absurdly small; single record massively over

	// Starvation guard: always return at least one record when stream has data.
	REQUIRE((batch.size() == 1));
	REQUIRE((batch.read_ids[0] == "big"));
}

TEST_CASE("SequenceReader byte budget SE: budget preserves full content across calls",
          "[SequenceReader][byte_budget]") {
	TempFileFixture fixture;
	auto path = "budget_preserve.fq";
	std::vector<std::string> expected_ids;
	std::vector<std::string> records;
	for (int i = 0; i < 7; i++) {
		std::string id = "rec" + std::to_string(i);
		expected_ids.push_back(id);
		records.push_back(fixture.simple_read(id, "ACGT", "IIII"));
	}
	fixture.write_temp_fastq(path, records);

	miint::SequenceReader reader(path);

	std::vector<std::string> seen;
	for (int call = 0; call < 20; call++) {
		auto batch = reader.read(100, 20); // tight budget; ~1 record per call
		if (batch.empty()) {
			break;
		}
		for (auto &id : batch.read_ids) {
			seen.push_back(id);
		}
	}

	REQUIRE((seen == expected_ids));
}

TEST_CASE("SequenceReader byte budget PE: budget counts both mates", "[SequenceReader][byte_budget]") {
	TempFileFixture fixture;
	auto r1 = "pe_budget_r1.fq";
	auto r2 = "pe_budget_r2.fq";

	std::vector<std::string> recs1, recs2;
	for (int i = 0; i < 5; i++) {
		std::string id = "pair" + std::to_string(i);
		recs1.push_back(fixture.simple_read(id, "ACGT", "IIII"));
		recs2.push_back(fixture.simple_read(id, "TGCA", "HHHH"));
	}
	fixture.write_temp_fastq(r1, recs1);
	fixture.write_temp_fastq(r2, recs2);

	miint::SequenceReader reader(r1, r2);

	// Each mate is ~5 (name) + 4 (seq) + 4 (qual) = 13 bytes; pair ~26 bytes.
	// Budget 40 bytes -> pair 1 fits, adding pair 2 (>52 bytes) exceeds -> 1 pair per call.
	size_t total = 0;
	for (int call = 0; call < 20; call++) {
		auto batch = reader.read(100, 40);
		if (batch.empty()) {
			break;
		}
		REQUIRE((batch.size() == 1));
		REQUIRE((batch.is_paired));
		total += batch.size();
	}
	REQUIRE((total == 5));
}

TEST_CASE("SequenceReader byte budget SE: large records do not over-prefetch", "[SequenceReader][byte_budget]") {
	// Each record alone exceeds the budget, so the byte budget caps the batch at one record.
	// The poll that feeds that batch must NOT eagerly materialize a fixed chunk of large
	// records: with a 500 B budget and ~4 KB records, a single poll should pull one record,
	// not POLL_CHUNK_SIZE of them. This is the memory bound the byte budget is meant to give.
	TempFileFixture fixture;
	auto path = "budget_large_prefetch.fq";
	std::vector<std::string> records;
	for (int i = 0; i < 12; i++) {
		records.push_back(fixture.simple_read("r" + std::to_string(i), std::string(2000, 'A'), std::string(2000, 'I')));
	}
	fixture.write_temp_fastq(path, records);

	miint::SequenceReader reader(path);

	size_t total = 0;
	for (int call = 0; call < 50; call++) {
		auto batch = reader.read(100, 500); // budget << one record
		if (batch.empty()) {
			break;
		}
		REQUIRE((batch.size() == 1));
		total += batch.size();
	}
	REQUIRE((total == 12));
	// The dynamic poll must size to the budget: at most one ~4 KB record per poll.
	REQUIRE((reader.MaxPollCount() == 1));
}

TEST_CASE("SequenceReader byte budget PE: large pairs do not over-prefetch", "[SequenceReader][byte_budget]") {
	TempFileFixture fixture;
	auto r1 = "budget_large_pe_r1.fq";
	auto r2 = "budget_large_pe_r2.fq";
	std::vector<std::string> recs1, recs2;
	for (int i = 0; i < 8; i++) {
		std::string id = "pair" + std::to_string(i);
		recs1.push_back(fixture.simple_read(id, std::string(2000, 'A'), std::string(2000, 'I')));
		recs2.push_back(fixture.simple_read(id, std::string(2000, 'C'), std::string(2000, 'H')));
	}
	fixture.write_temp_fastq(r1, recs1);
	fixture.write_temp_fastq(r2, recs2);

	miint::SequenceReader reader(r1, r2);

	size_t total = 0;
	for (int call = 0; call < 50; call++) {
		auto batch = reader.read(100, 500); // budget << one pair (~8 KB)
		if (batch.empty()) {
			break;
		}
		REQUIRE((batch.size() == 1));
		total += batch.size();
	}
	REQUIRE((total == 8));
	REQUIRE((reader.MaxPollCount() == 1));
}

TEST_CASE("SequenceReader byte budget SE: small records still poll in chunks", "[SequenceReader][byte_budget]") {
	// Guard against the dynamic poll collapsing to 1 for the common short-read case: with a
	// generous budget and tiny records, the poll must still amortize (reach POLL_CHUNK_SIZE).
	TempFileFixture fixture;
	auto path = "budget_small_chunked.fq";
	std::vector<std::string> records;
	for (int i = 0; i < 200; i++) {
		records.push_back(fixture.simple_read("r" + std::to_string(i), "ACGT", "IIII"));
	}
	fixture.write_temp_fastq(path, records);

	miint::SequenceReader reader(path);
	auto batch = reader.read(100, 1024 * 1024); // generous budget, fills to n
	REQUIRE((batch.size() == 100));
	// Must still batch the poll, not degrade to one-record-at-a-time.
	REQUIRE((reader.MaxPollCount() > 1));
}

TEST_CASE("SequenceReader byte budget default (SIZE_MAX) unchanged behavior", "[SequenceReader][byte_budget]") {
	// Regression: callers that don't pass a budget see the same rows-only behavior.
	TempFileFixture fixture;
	auto path = "budget_default.fq";
	std::vector<std::string> records;
	for (int i = 0; i < 50; i++) {
		records.push_back(fixture.simple_read("r" + std::to_string(i), "ACGT", "IIII"));
	}
	fixture.write_temp_fastq(path, records);

	miint::SequenceReader reader(path);
	auto batch = reader.read(100); // no budget argument -> default SIZE_MAX
	REQUIRE((batch.size() == 50));
}
