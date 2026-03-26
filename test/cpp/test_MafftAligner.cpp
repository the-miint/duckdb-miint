#include "../../src/include/MafftAligner.hpp"
#include <catch2/catch_test_macros.hpp>
#include <fstream>
#include <sstream>
#include <tuple>

using miint::MafftAligner;
using miint::MafftAlignResult;

// Helper: parse a FASTA file into names, comments, sequences
static std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>>
parse_fasta(const std::string &path) {
	std::vector<std::string> names, comments, seqs;
	std::ifstream f(path);
	REQUIRE(f.is_open());
	std::string line, current_seq;
	std::string current_name, current_comment;
	bool in_seq = false;
	while (std::getline(f, line)) {
		if (line.empty()) {
			continue;
		}
		if (line[0] == '>') {
			if (in_seq) {
				names.push_back(current_name);
				comments.push_back(current_comment);
				seqs.push_back(current_seq);
				current_seq.clear();
			}
			// Parse header: >name comment
			std::string header = line.substr(1);
			auto space = header.find(' ');
			if (space != std::string::npos) {
				current_name = header.substr(0, space);
				current_comment = header.substr(space + 1);
			} else {
				current_name = header;
				current_comment = "";
			}
			in_seq = true;
		} else {
			current_seq += line;
		}
	}
	if (in_seq) {
		names.push_back(current_name);
		comments.push_back(current_comment);
		seqs.push_back(current_seq);
	}
	return {names, comments, seqs};
}

// Cycle 2.1: Construction
TEST_CASE("MafftAligner - Construction", "[MafftAligner]") {
	SECTION("Default construction succeeds") {
		REQUIRE_NOTHROW(MafftAligner());
	}
}

// Cycle 2.2: Two trivial identical sequences
TEST_CASE("MafftAligner - Two identical sequences", "[MafftAligner]") {
	MafftAligner aligner;
	auto result = aligner.align({"seq1", "seq2"}, {"", ""}, {"ACGTACGT", "ACGTACGT"});

	REQUIRE(result.sequences.size() == 2);
	REQUIRE(result.aligned_length > 0);
	REQUIRE(result.sequences[0] == result.sequences[1]);
	REQUIRE(result.original_lengths[0] == 8);
	REQUIRE(result.original_lengths[1] == 8);
	REQUIRE(result.sequences[0].size() == static_cast<size_t>(result.aligned_length));
}

// Cycle 2.3: Real sequences from test.fa
TEST_CASE("MafftAligner - Real sequences", "[MafftAligner]") {
	MafftAligner aligner;

	auto [names, comments, seqs] = parse_fasta("data/fastq/test.fa");
	REQUIRE(names.size() >= 2);

	auto result = aligner.align(names, comments, seqs);

	REQUIRE(result.sequences.size() == names.size());
	// All aligned sequences have the same length
	for (const auto &s : result.sequences) {
		REQUIRE(s.size() == static_cast<size_t>(result.aligned_length));
	}
	// Names preserved
	REQUIRE(result.names[0] == names[0]);
	// Original lengths preserved
	for (size_t i = 0; i < seqs.size(); i++) {
		REQUIRE(result.original_lengths[i] == static_cast<int>(seqs[i].size()));
	}
	// Non-gap chars match original (case-insensitive since MAFFT may change case,
	// but our wrapper restores it)
	for (size_t i = 0; i < result.sequences.size(); i++) {
		std::string ungapped;
		for (char c : result.sequences[i]) {
			if (c != '-') {
				ungapped += c;
			}
		}
		REQUIRE(ungapped == seqs[i]);
	}
}

// Cycle 2.4: Error handling
TEST_CASE("MafftAligner - Error handling", "[MafftAligner]") {
	MafftAligner aligner;

	SECTION("Fewer than 2 sequences throws") {
		REQUIRE_THROWS_AS(aligner.align({"seq1"}, {""}, {"ACGTACGT"}), std::invalid_argument);
	}
	SECTION("Empty sequences vector throws") {
		REQUIRE_THROWS_AS(aligner.align({}, {}, {}), std::invalid_argument);
	}
	SECTION("Mismatched names/sequences sizes throws") {
		REQUIRE_THROWS_AS(aligner.align({"seq1"}, {""}, {"ACGT", "TGCA"}), std::invalid_argument);
	}
	SECTION("Empty sequence in input throws") {
		REQUIRE_THROWS_AS(aligner.align({"seq1", "seq2"}, {"", ""}, {"ACGTACGT", ""}), std::invalid_argument);
	}
	SECTION("Sequence shorter than 6 chars throws") {
		REQUIRE_THROWS_AS(aligner.align({"seq1", "seq2"}, {"", ""}, {"ACGT", "ACGTACGT"}), std::invalid_argument);
	}
	SECTION("Sequence of exactly 5 chars throws") {
		REQUIRE_THROWS_AS(aligner.align({"seq1", "seq2"}, {"", ""}, {"ACGTA", "ACGTACGT"}), std::invalid_argument);
	}
}

// Cycle 2.5: Case preservation
TEST_CASE("MafftAligner - Case preservation", "[MafftAligner]") {
	MafftAligner aligner;

	auto result = aligner.align({"s1", "s2"}, {"", ""}, {"AcGtAcGt", "aCgTaCgT"});

	REQUIRE(result.sequences.size() == 2);

	// Non-gap characters must match original case exactly
	for (size_t i = 0; i < result.sequences.size(); i++) {
		std::string ungapped;
		for (char c : result.sequences[i]) {
			if (c != '-') {
				ungapped += c;
			}
		}
		std::vector<std::string> originals = {"AcGtAcGt", "aCgTaCgT"};
		REQUIRE(ungapped == originals[i]);
	}
}

// Cycle 2.6: Reuse (multiple calls)
TEST_CASE("MafftAligner - Reuse", "[MafftAligner]") {
	MafftAligner aligner;

	auto r1 = aligner.align({"a", "b"}, {"", ""}, {"AAAAAAAA", "TTTTTTTT"});
	auto r2 = aligner.align({"c", "d", "e"}, {"", "", ""}, {"ACGTACGT", "ACGTACGT", "TGCATGCA"});

	REQUIRE(r1.sequences.size() == 2);
	REQUIRE(r2.sequences.size() == 3);
	// Second call should not be corrupted by first
	REQUIRE(r2.names[0] == "c");
	REQUIRE(r2.names[1] == "d");
	REQUIRE(r2.names[2] == "e");
}

// Cycle 2.7: Reference correctness (larger dataset)
TEST_CASE("MafftAligner - Sample dataset", "[MafftAligner]") {
	MafftAligner aligner;

	auto [names, comments, seqs] = parse_fasta("data/mafft/sample.fa");
	REQUIRE(names.size() == 36);

	auto result = aligner.align(names, comments, seqs);

	REQUIRE(result.sequences.size() == 36);
	REQUIRE(result.aligned_length > 0);
	// All aligned sequences have the same length
	for (const auto &s : result.sequences) {
		REQUIRE(s.size() == static_cast<size_t>(result.aligned_length));
	}
}
