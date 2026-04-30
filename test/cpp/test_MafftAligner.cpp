/*
 * MafftAligner unit tests.
 *
 * Reference alignments were generated using the splittbfast binary compiled
 * from the same MAFFT source (scratch/mafft/core/, library-api branch) with:
 *
 *   scratch/mafft/core/splittbfast -D -f -1.53 -Q 100 -h 0 -p 50 -s -1 -i <input>
 *
 * For protein (sample.fa):
 *   scratch/mafft/core/splittbfast -P -f -1.53 -Q 100 -h 0 -p 50 -s -1 -i <input>
 *
 * These flags match the argv constructed by MafftAligner::align() and
 * correspond to the PartTree algorithm invoked by:
 *   mafft --quiet --preservecase --parttree <input>
 *
 * The binary output is lowercased (DNA) or uppercased (protein) and
 * reordered by similarity. MafftAligner restores original case and
 * preserves input order. Tests below account for both differences.
 *
 * Correctness was cross-verified against native mafft v7.526 (bioconda,
 * conda env "mafft"). The sorted set of ungapped sequences is identical
 * between the binary, the library, and native mafft. Aligned sequences
 * (with gaps) are identical between binary and library when sorted.
 * Native mafft v7.526 produces a different alignment width (752 vs 732
 * for sample.fa) because it is a different release with potentially
 * different default parameters; however the ungapped content matches.
 */

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

// Cycle 2.2: Two trivial identical sequences — exact expected output
// Reference: splittbfast -D ... produces "acgtacgt" for both (lowercased).
// MafftAligner restores original case → "ACGTACGT".
// No gaps needed (identical sequences, same length).
TEST_CASE("MafftAligner - Two identical sequences", "[MafftAligner]") {
	MafftAligner aligner;
	auto result = aligner.align({"seq1", "seq2"}, {"", ""}, {"ACGTACGT", "ACGTACGT"});

	REQUIRE(result.sequences.size() == 2);
	REQUIRE(result.aligned_length == 8);
	// Exact expected alignment: no gaps, original case
	REQUIRE(result.sequences[0] == "ACGTACGT");
	REQUIRE(result.sequences[1] == "ACGTACGT");
	REQUIRE(result.original_lengths[0] == 8);
	REQUIRE(result.original_lengths[1] == 8);
}

// Cycle 2.3: test.fa — exact expected output
// Input: seq1=ATGCATGCATGC (12bp), seq2=GGCCGGCCGGCC (12bp)
// Reference: splittbfast -D ... produces lowercased, no gaps (same length).
// MafftAligner restores case → original sequences unchanged.
TEST_CASE("MafftAligner - Real sequences", "[MafftAligner]") {
	MafftAligner aligner;

	auto [names, comments, seqs] = parse_fasta("data/fastq/test.fa");
	REQUIRE(names.size() == 2);

	auto result = aligner.align(names, comments, seqs);

	REQUIRE(result.sequences.size() == 2);
	REQUIRE(result.aligned_length == 12);
	// Exact expected alignment: no gaps (same-length sequences)
	REQUIRE(result.sequences[0] == "ATGCATGCATGC");
	REQUIRE(result.sequences[1] == "GGCCGGCCGGCC");
	// Names preserved in input order
	REQUIRE(result.names[0] == "seq1");
	REQUIRE(result.names[1] == "seq2");
	REQUIRE(result.original_lengths[0] == 12);
	REQUIRE(result.original_lengths[1] == 12);
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
		REQUIRE_THROWS_AS(aligner.align({"seq1"}, {""}, {"ACGTACGT", "TGCATGCA"}), std::invalid_argument);
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

// Cycle 2.5: Case preservation — exact expected output
// Input: s1=AcGtAcGt, s2=aCgTaCgT (mixed case)
// Reference: splittbfast -D ... lowercases to "acgtacgt" for both.
// MafftAligner restores original case character-by-character (fillorichar algorithm).
TEST_CASE("MafftAligner - Case preservation", "[MafftAligner]") {
	MafftAligner aligner;

	auto result = aligner.align({"s1", "s2"}, {"", ""}, {"AcGtAcGt", "aCgTaCgT"});

	REQUIRE(result.sequences.size() == 2);
	REQUIRE(result.aligned_length == 8);
	// Exact expected output: original case restored, no gaps
	REQUIRE(result.sequences[0] == "AcGtAcGt");
	REQUIRE(result.sequences[1] == "aCgTaCgT");
}

// Cycle 2.6: Reuse (multiple calls) — exact expected output
// Verifies global state cleanup between calls (the static buffer reset
// sentinels in localcommonsextet_p, pairalign, splitseq_mq).
TEST_CASE("MafftAligner - Reuse", "[MafftAligner]") {
	MafftAligner aligner;

	auto r1 = aligner.align({"a", "b"}, {"", ""}, {"AAAAAAAA", "TTTTTTTT"});
	REQUIRE(r1.sequences.size() == 2);
	REQUIRE(r1.aligned_length == 8);
	REQUIRE(r1.sequences[0] == "AAAAAAAA");
	REQUIRE(r1.sequences[1] == "TTTTTTTT");

	auto r2 = aligner.align({"c", "d", "e"}, {"", "", ""}, {"ACGTACGT", "ACGTACGT", "TGCATGCA"});
	REQUIRE(r2.sequences.size() == 3);
	// Second call produces valid results not corrupted by first
	REQUIRE(r2.names[0] == "c");
	REQUIRE(r2.names[1] == "d");
	REQUIRE(r2.names[2] == "e");
	// All sequences must have the same aligned length
	for (const auto &s : r2.sequences) {
		REQUIRE(s.size() == static_cast<size_t>(r2.aligned_length));
	}
}

// Cycle 2.7: Reference correctness — 36 protein opsin sequences
// Input: data/mafft/sample.fa (36 opsin sequences from MAFFT test suite)
// MafftAligner now drives MAFFT through MAFFT_STRATEGY_AUTO (FFT-NS-2 for
// inputs of this size, not PartTree). The exact aligned width is therefore
// algorithm-dependent and may differ from the historical PartTree value of
// 732. The invariants we still assert are: (a) 36 output rows; (b) all rows
// share a common width; (c) ungapped content is byte-identical to the input
// after case restoration.
TEST_CASE("MafftAligner - Sample dataset", "[MafftAligner]") {
	MafftAligner aligner;

	auto [names, comments, seqs] = parse_fasta("data/mafft/sample.fa");
	REQUIRE(names.size() == 36);

	auto result = aligner.align(names, comments, seqs);

	REQUIRE(result.sequences.size() == 36);
	REQUIRE(result.aligned_length > 0);
	// All aligned sequences share the reported aligned_length
	for (const auto &s : result.sequences) {
		REQUIRE(s.size() == static_cast<size_t>(result.aligned_length));
	}
	// Ungapped content must match original input (case-preserved)
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

// Cycle 2.8: tiny.fa — exact expected output
// Input: s1=ACGTACGT, s2=ACGAACGT, s3=ACGTACGA (3 DNA seqs, 8bp each)
// Reference: splittbfast -D ... produces same sequences lowercased, no gaps.
// MafftAligner restores case. Verified against binary output.
TEST_CASE("MafftAligner - Tiny exact alignment", "[MafftAligner]") {
	MafftAligner aligner;

	auto result =
	    aligner.align({"s1", "s2", "s3"}, {"first sequence", "", "third one"}, {"ACGTACGT", "ACGAACGT", "ACGTACGA"});

	REQUIRE(result.sequences.size() == 3);
	REQUIRE(result.aligned_length == 8);
	// Exact expected alignment: no gaps needed (all same length, high similarity)
	REQUIRE(result.sequences[0] == "ACGTACGT");
	REQUIRE(result.sequences[1] == "ACGAACGT");
	REQUIRE(result.sequences[2] == "ACGTACGA");
	// Names and comments preserved in input order
	REQUIRE(result.names[0] == "s1");
	REQUIRE(result.comments[0] == "first sequence");
	REQUIRE(result.names[1] == "s2");
	REQUIRE(result.comments[1] == "");
	REQUIRE(result.names[2] == "s3");
	REQUIRE(result.comments[2] == "third one");
}
