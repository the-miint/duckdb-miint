#ifdef MIINT_HAS_ABPOA
#include "AbpoaAligner.hpp"
#include <catch2/catch_test_macros.hpp>
#include <algorithm>

TEST_CASE("AbpoaAligner - Construction", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	(void)aligner;
}

TEST_CASE("AbpoaAligner - Two identical sequences", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	std::vector<std::string> names = {"s1", "s2"};
	std::vector<std::string> seqs = {"ACGTACGT", "ACGTACGT"};

	auto result = aligner.align(names, seqs);

	REQUIRE(result.names.size() == 2);
	REQUIRE(result.aligned_sequences.size() == 2);
	REQUIRE(result.aligned_length == 8);
	CHECK(result.aligned_sequences[0] == "ACGTACGT");
	CHECK(result.aligned_sequences[1] == "ACGTACGT");
	CHECK(result.original_lengths[0] == 8);
	CHECK(result.original_lengths[1] == 8);
}

TEST_CASE("AbpoaAligner - Three similar sequences (tiny.fa)", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	std::vector<std::string> names = {"s1", "s2", "s3"};
	std::vector<std::string> seqs = {"ACGTACGT", "ACGAACGT", "ACGTACGA"};

	auto result = aligner.align(names, seqs);

	REQUIRE(result.names.size() == 3);
	REQUIRE(result.aligned_length == 8);
	CHECK(result.aligned_sequences[0] == "ACGTACGT");
	CHECK(result.aligned_sequences[1] == "ACGAACGT");
	CHECK(result.aligned_sequences[2] == "ACGTACGA");
}

TEST_CASE("AbpoaAligner - Error: fewer than 2 sequences", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	CHECK_THROWS_AS(aligner.align({"s1"}, {"ACGTACGT"}), std::invalid_argument);
	CHECK_THROWS_AS(aligner.align({}, {}), std::invalid_argument);
}

TEST_CASE("AbpoaAligner - Error: empty sequence", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	CHECK_THROWS_AS(aligner.align({"s1", "s2"}, {"ACGT", ""}), std::invalid_argument);
}

TEST_CASE("AbpoaAligner - Error: mismatched names/sequences", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	CHECK_THROWS_AS(aligner.align({"s1"}, {"ACGT", "TGCA"}), std::invalid_argument);
}

TEST_CASE("AbpoaAligner - Consensus: identical sequences", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	std::vector<std::string> names = {"s1", "s2"};
	std::vector<std::string> seqs = {"ACGTACGT", "ACGTACGT"};

	auto result = aligner.consensus(names, seqs);

	REQUIRE(result.entries.size() == 1);
	CHECK(result.entries[0].consensus_id == 0);
	CHECK(result.entries[0].sequence == "ACGTACGT");
	CHECK(result.entries[0].length == 8);
	CHECK(result.entries[0].num_reads == 2);
}

TEST_CASE("AbpoaAligner - Consensus: tiny.fa ground truth", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	std::vector<std::string> names = {"s1", "s2", "s3"};
	std::vector<std::string> seqs = {"ACGTACGT", "ACGAACGT", "ACGTACGA"};

	auto result = aligner.consensus(names, seqs);

	REQUIRE(result.entries.size() == 1);
	CHECK(result.entries[0].sequence == "ACGTACGT");
	CHECK(result.entries[0].length == 8);
	CHECK(result.entries[0].num_reads == 3);
}

TEST_CASE("AbpoaAligner - Consensus: valid nucleotides only", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	std::vector<std::string> names = {"s1", "s2", "s3"};
	std::vector<std::string> seqs = {"ACGTACGTACGT", "ACGAACGAACGA", "TGCATGCATGCA"};

	auto result = aligner.consensus(names, seqs);

	REQUIRE(!result.entries.empty());
	for (auto &entry : result.entries) {
		for (char c : entry.sequence) {
			bool valid = (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
			CHECK(valid);
		}
	}
}

TEST_CASE("AbpoaAligner - Parameter forwarding", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	std::vector<std::string> names = {"s1", "s2"};
	std::vector<std::string> seqs = {"ACGTACGT", "ACGAACGT"};

	miint::AbpoaAlignParams params;
	params.match = 1;
	params.mismatch = 2;
	auto result = aligner.align(names, seqs, params);

	REQUIRE(result.names.size() == 2);
	REQUIRE(result.aligned_length > 0);
}

TEST_CASE("AbpoaAligner - Reuse: call align twice", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	std::vector<std::string> names = {"s1", "s2"};
	std::vector<std::string> seqs = {"ACGTACGT", "ACGAACGT"};

	auto result1 = aligner.align(names, seqs);
	auto result2 = aligner.align(names, seqs);

	CHECK(result1.aligned_length == result2.aligned_length);
	CHECK(result1.aligned_sequences == result2.aligned_sequences);
}

TEST_CASE("AbpoaAligner - Multi-consensus", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	std::vector<std::string> names = {"a", "b", "c", "d", "e", "f"};
	std::vector<std::string> seqs = {"ACGTACGTACGT", "ACGTTTTTACGT", "ACGAACGAACGA",
	                                 "ACGATTTTACGA", "ACGTACGTACGT", "ACGTTTTTACGT"};

	miint::AbpoaAlignParams params;
	params.max_num_cons = 2;
	auto result = aligner.consensus(names, seqs, params);

	REQUIRE(result.entries.size() >= 1);
	REQUIRE(result.entries.size() <= 2);
	for (auto &entry : result.entries) {
		CHECK(entry.length > 0);
		CHECK(entry.num_reads > 0);
	}
}

TEST_CASE("AbpoaAligner - Progressive mode", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	std::vector<std::string> names = {"s1", "s2", "s3"};
	std::vector<std::string> seqs = {"ACGTACGT", "ACGAACGT", "ACGTACGA"};

	miint::AbpoaAlignParams params;
	params.progressive = true;
	auto result = aligner.align(names, seqs, params);

	REQUIRE(result.names.size() == 3);
	REQUIRE(result.aligned_length == 8);
}

TEST_CASE("AbpoaAligner - Majority voting algorithm", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	std::vector<std::string> names = {"s1", "s2", "s3"};
	std::vector<std::string> seqs = {"ACGTACGT", "ACGAACGT", "ACGTACGA"};

	miint::AbpoaAlignParams params;
	params.algorithm = "majority_voting";
	auto result = aligner.consensus(names, seqs, params);

	REQUIRE(result.entries.size() == 1);
	CHECK(result.entries[0].length > 0);
}

TEST_CASE("AbpoaAligner - Invalid align mode", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	std::vector<std::string> names = {"s1", "s2"};
	std::vector<std::string> seqs = {"ACGTACGT", "ACGAACGT"};

	miint::AbpoaAlignParams params;
	params.align_mode = "invalid";
	CHECK_THROWS_AS(aligner.align(names, seqs, params), std::invalid_argument);
}

TEST_CASE("AbpoaAligner - Ungapped content preserved", "[AbpoaAligner]") {
	miint::AbpoaAligner aligner;
	std::vector<std::string> names = {"s1", "s2", "s3"};
	std::vector<std::string> seqs = {"ACGTACGT", "ACGAACGT", "ACGTACGA"};

	auto result = aligner.align(names, seqs);

	for (size_t i = 0; i < seqs.size(); i++) {
		std::string ungapped;
		for (char c : result.aligned_sequences[i]) {
			if (c != '-') {
				ungapped += c;
			}
		}
		CHECK(ungapped == seqs[i]);
	}
}

#endif // MIINT_HAS_ABPOA
