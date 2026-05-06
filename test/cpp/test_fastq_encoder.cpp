// SPDX-License-Identifier: MIT
//
// Phase 6 RED/GREEN: byte-exact FASTQ encoder coverage.
//
// The encoder is the only piece of the upload pipeline that's deterministic
// without any HTTP/FTP harness, so we pin its output here. Quality offset
// validation, comment-on/comment-off, paired-style interleaving (single
// stream with two records), and the Phred+33 happy path all live here.

#include "fastq_encoder.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cstdint>
#include <string>
#include <vector>

using duckdb::FastqEncoder;

namespace {

struct StringSink {
	std::string buf;
	void operator()(const char *data, std::size_t size) {
		buf.append(data, size);
	}
};

// Phred+33 inverse: pass raw quality values 0..93. Encoder adds 33 → '!'..'~'.
std::vector<std::uint8_t> RawQuality(const std::string &phred33_chars) {
	std::vector<std::uint8_t> out;
	out.reserve(phred33_chars.size());
	for (char c : phred33_chars) {
		out.push_back(static_cast<std::uint8_t>(c) - 33);
	}
	return out;
}

} // namespace

TEST_CASE("FastqEncoder: single Phred+33 record byte-exact", "[fastq_encoder]") {
	FastqEncoder encoder;
	StringSink sink;
	const std::string id = "read1";
	const std::string seq = "ACGTACGTAC";
	const auto qual = RawQuality("IIIIIIIIII"); // 'I' = Q40

	encoder.Encode(std::ref(sink), id.data(), id.size(), nullptr, 0, seq.data(), seq.size(), qual.data(), qual.size());

	REQUIRE(sink.buf == "@read1\nACGTACGTAC\n+\nIIIIIIIIII\n");
}

TEST_CASE("FastqEncoder: comment is space-prefixed on the header line", "[fastq_encoder]") {
	FastqEncoder encoder;
	StringSink sink;
	const std::string id = "read1";
	const std::string comment = "1:N:0:ATCG";
	const std::string seq = "ACGT";
	const auto qual = RawQuality("ABCD");

	encoder.Encode(std::ref(sink), id.data(), id.size(), comment.data(), comment.size(), seq.data(), seq.size(),
	               qual.data(), qual.size());

	REQUIRE(sink.buf == "@read1 1:N:0:ATCG\nACGT\n+\nABCD\n");
}

TEST_CASE("FastqEncoder: empty comment emits no separator", "[fastq_encoder]") {
	FastqEncoder encoder;
	StringSink sink_a, sink_b;
	const std::string id = "r";
	const std::string seq = "AC";
	const auto qual = RawQuality("!!");

	encoder.Encode(std::ref(sink_a), id.data(), id.size(), "", 0, seq.data(), seq.size(), qual.data(), qual.size());
	encoder.Encode(std::ref(sink_b), id.data(), id.size(), nullptr, 0, seq.data(), seq.size(), qual.data(),
	               qual.size());

	REQUIRE(sink_a.buf == sink_b.buf);
	REQUIRE(sink_a.buf == "@r\nAC\n+\n!!\n");
}

TEST_CASE("FastqEncoder: three reads paired-end interleaved into one stream", "[fastq_encoder]") {
	FastqEncoder encoder;
	StringSink sink;

	struct Record {
		std::string id;
		std::string seq;
		std::string qual_chars;
	};
	const std::vector<std::pair<Record, Record>> reads = {
	    {{"r1", "AAAA", "IIII"}, {"r1", "TTTT", "IIII"}},
	    {{"r2", "CCCC", "????"}, {"r2", "GGGG", "????"}},
	    {{"r3", "ACGT", "AAAA"}, {"r3", "TGCA", "BBBB"}},
	};

	for (const auto &p : reads) {
		const auto q1 = RawQuality(p.first.qual_chars);
		const auto q2 = RawQuality(p.second.qual_chars);
		encoder.Encode(std::ref(sink), p.first.id.data(), p.first.id.size(), nullptr, 0, p.first.seq.data(),
		               p.first.seq.size(), q1.data(), q1.size());
		encoder.Encode(std::ref(sink), p.second.id.data(), p.second.id.size(), nullptr, 0, p.second.seq.data(),
		               p.second.seq.size(), q2.data(), q2.size());
	}

	const std::string expected = "@r1\nAAAA\n+\nIIII\n"
	                             "@r1\nTTTT\n+\nIIII\n"
	                             "@r2\nCCCC\n+\n????\n"
	                             "@r2\nGGGG\n+\n????\n"
	                             "@r3\nACGT\n+\nAAAA\n"
	                             "@r3\nTGCA\n+\nBBBB\n";
	REQUIRE(sink.buf == expected);
}

TEST_CASE("FastqEncoder: rejects out-of-range quality scores", "[fastq_encoder]") {
	FastqEncoder encoder;
	StringSink sink;
	const std::string id = "r";
	const std::string seq = "A";
	// 94 + 33 = 127, exceeds the printable-ASCII upper bound (126).
	const std::vector<std::uint8_t> qual = {94};

	REQUIRE_THROWS_WITH(encoder.Encode(std::ref(sink), id.data(), id.size(), nullptr, 0, seq.data(), seq.size(),
	                                   qual.data(), qual.size()),
	                    Catch::Matchers::ContainsSubstring("Quality score overflow"));
}

TEST_CASE("FastqEncoder: Phred+64 offset routes through the same validator", "[fastq_encoder]") {
	FastqEncoder encoder(64);
	StringSink sink;
	const std::string id = "r";
	const std::string seq = "A";
	// Raw 0 + 64 = 64 = '@'; valid.
	const std::vector<std::uint8_t> qual = {0};
	encoder.Encode(std::ref(sink), id.data(), id.size(), nullptr, 0, seq.data(), seq.size(), qual.data(), qual.size());
	REQUIRE(sink.buf == "@r\nA\n+\n@\n");

	// Raw 63 + 64 = 127, must overflow.
	StringSink sink_bad;
	const std::vector<std::uint8_t> qual_bad = {63};
	REQUIRE_THROWS_WITH(encoder.Encode(std::ref(sink_bad), id.data(), id.size(), nullptr, 0, seq.data(), seq.size(),
	                                   qual_bad.data(), qual_bad.size()),
	                    Catch::Matchers::ContainsSubstring("Quality score overflow"));
}

TEST_CASE("FastqEncoder: encoder is reusable across many calls without state leak", "[fastq_encoder]") {
	FastqEncoder encoder;
	StringSink sink;
	const std::string id = "r";
	const auto qual_short = RawQuality("I");
	const auto qual_long = RawQuality("IIIIIIIIIIIIIIII");

	// Encode a long record then a short one: the short call must NOT see
	// trailing bytes from the long record's quality buffer.
	encoder.Encode(std::ref(sink), id.data(), id.size(), nullptr, 0, "AAAAAAAAAAAAAAAA", 16, qual_long.data(),
	               qual_long.size());
	encoder.Encode(std::ref(sink), id.data(), id.size(), nullptr, 0, "C", 1, qual_short.data(), qual_short.size());

	REQUIRE(sink.buf == "@r\nAAAAAAAAAAAAAAAA\n+\nIIIIIIIIIIIIIIII\n"
	                    "@r\nC\n+\nI\n");
}
