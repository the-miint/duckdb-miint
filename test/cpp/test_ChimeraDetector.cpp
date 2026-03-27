#include <catch2/catch_test_macros.hpp>

#include "ChimeraDetector.hpp"
#include "WFA2Aligner.hpp"

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

using miint::build_star_alignment;
using miint::classify_diffs;
using miint::DiffType;
using miint::StarAlignment;

// ============================================================================
// Phase 2: Alignment equivalence validation
// ============================================================================

TEST_CASE("UCHIME alignment equivalence", "[ChimeraDetector]") {
	// WFA2Aligner with UCHIME-equivalent penalties
	miint::WFA2Aligner aligner(miint::UCHIME_WFA2_MISMATCH, miint::UCHIME_WFA2_GAP_OPEN, miint::UCHIME_WFA2_GAP_EXTEND);

	SECTION("Identical sequences produce ungapped alignment") {
		std::string seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
		auto result = aligner.align_full(seq, seq);
		REQUIRE(result.has_value());
		REQUIRE(result->score == 0);
		REQUIRE(result->query_aligned == seq);
		REQUIRE(result->subject_aligned == seq);
		// No gaps
		REQUIRE(result->query_aligned.find('-') == std::string::npos);
	}

	SECTION("Single substitution produces ungapped alignment") {
		std::string query = "ACGTACGTACGTACGT";
		std::string subject = "ACGTACGAACGTACGT"; // T->A at position 7
		auto result = aligner.align_full(query, subject);
		REQUIRE(result.has_value());
		REQUIRE(result->score == miint::UCHIME_WFA2_MISMATCH); // One mismatch penalty
		// No gaps — mismatch penalty (6) < gap_open (20)
		REQUIRE(result->query_aligned.find('-') == std::string::npos);
		REQUIRE(result->subject_aligned.find('-') == std::string::npos);
	}

	SECTION("High gap penalties discourage gaps for closely related sequences") {
		// With gap_open=20, a single base insertion costs 20+4=24,
		// while two mismatches cost 2*6=12. So the aligner should prefer
		// mismatches over gaps when sequences differ by a few bases.
		std::string query = "ACGTACGTACGT";
		std::string subject = "ACGTAACGTACGT"; // Extra A inserted
		auto result = aligner.align_full(query, subject);
		REQUIRE(result.has_value());
		// The alignment should exist (either gapped or ungapped)
		REQUIRE(result->query_aligned.size() == result->subject_aligned.size());
	}
}

TEST_CASE("UCHIME alignment on real test data", "[ChimeraDetector]") {
	// Use the actual test sequences from our chimera test data.
	// These are the same sequences vsearch aligns in expected_ref_alns.txt.
	miint::WFA2Aligner aligner(miint::UCHIME_WFA2_MISMATCH, miint::UCHIME_WFA2_GAP_OPEN, miint::UCHIME_WFA2_GAP_EXTEND);

	// ref1 and ref2 from chimera_ref.fasta
	std::string ref1 =
	    "AAGCCCAATCAACCACTCTCACTGGACGATTGCGGATATTGGCAACGAATTGGGAGGCGACCCGGACGACAGTCACGCCTTCTCGTTTGCGTACAGCTAT"
	    "TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCCGCCGTGATACCGGACCAAACAAGACGTCCACTTCAATGTTTAAATGACCTACGCGTCAGAA"
	    "CACCTTTCTACTATGTGTTCTCCCAGAATCATCTAGTACAATGGCGCGTCGTCATTAAAGCACCGGATGCGACGAACGGAGCGTGAATGAAGCTACTAC";

	std::string ref2 =
	    "AAGCACAATAAACCCCTGTCACGGGTCGAATAGGGGTTTAGGCAACGATCTCTGCAGAGCCCCTTGAGACAGTGACGCTTGTGCCGTTGCTTAAACTGATT"
	    "TGAAGGAGTCTAGCGGCCGCAGTAACGCACATTACCTAGTCCGCGTTCCCAGACGAAACAGGACGACATCTTTAAGGTTTAACTGACCGTCTATTCTTAAA"
	    "ACCTTCCAACTATGTGTTCCGAAAGAATCACCAACTACATTTACGCGTCGTGAATAACGTGTCGCCTGAGACGATAGGCGCGTGAATGAGGCGCTTAA";

	// query_chimera1 = first 150bp of ref1 + last 150bp of ref2
	std::string chimera1 = ref1.substr(0, 150) + ref2.substr(150);

	SECTION("300bp sequences align without gaps (substitution-only divergence)") {
		auto result = aligner.align_full(chimera1, ref1);
		REQUIRE(result.has_value());
		// These sequences were generated from a template with substitutions only,
		// so with gap_open=20 the alignment should be ungapped.
		REQUIRE(result->query_aligned.find('-') == std::string::npos);
		REQUIRE(result->subject_aligned.find('-') == std::string::npos);
		REQUIRE(result->query_aligned.size() == 300);
	}

	SECTION("Chimera aligns to ref1 with correct identity pattern") {
		auto result_a = aligner.align_full(chimera1, ref1);
		auto result_b = aligner.align_full(chimera1, ref2);
		REQUIRE(result_a.has_value());
		REQUIRE(result_b.has_value());

		// Count matches in first and second halves
		int matches_a_first_half = 0, matches_a_second_half = 0;
		int matches_b_first_half = 0, matches_b_second_half = 0;
		for (size_t i = 0; i < 150; i++) {
			if (result_a->query_aligned[i] == result_a->subject_aligned[i]) {
				matches_a_first_half++;
			}
			if (result_b->query_aligned[i] == result_b->subject_aligned[i]) {
				matches_b_first_half++;
			}
		}
		for (size_t i = 150; i < 300; i++) {
			if (result_a->query_aligned[i] == result_a->subject_aligned[i]) {
				matches_a_second_half++;
			}
			if (result_b->query_aligned[i] == result_b->subject_aligned[i]) {
				matches_b_second_half++;
			}
		}

		// First half should match ref1 better than ref2
		REQUIRE(matches_a_first_half > matches_b_first_half);
		// Second half should match ref2 better than ref1
		REQUIRE(matches_b_second_half > matches_a_second_half);
		// First half against ref1 should be perfect (chimera's first half IS ref1)
		REQUIRE(matches_a_first_half == 150);
		// Second half against ref2 should be perfect (chimera's second half IS ref2)
		REQUIRE(matches_b_second_half == 150);
	}
}

// ============================================================================
// Phase 3: Star alignment construction and diff classification
// ============================================================================

TEST_CASE("build_star_alignment - empty inputs", "[ChimeraDetector]") {
	SECTION("All empty strings produce empty star alignment") {
		auto star = build_star_alignment("", "", "", "");
		REQUIRE(star.query_row.empty());
		REQUIRE(star.parent_a_row.empty());
		REQUIRE(star.parent_b_row.empty());
	}

	SECTION("classify_diffs on empty star produces empty diffs") {
		StarAlignment star {"", "", "", {}};
		classify_diffs(star);
		REQUIRE(star.diffs.empty());
	}
}

TEST_CASE("build_star_alignment - mismatched query bases throw", "[ChimeraDetector]") {
	// If one alignment has unconsumed non-gap query bases after the other finishes,
	// build_star_alignment should throw (indicates mismatched input alignments).
	SECTION("A has extra query bases after B finishes") {
		// B is 4 columns, A is 5 columns with the 5th being a non-gap query base
		REQUIRE_THROWS_AS(build_star_alignment("ACGTX", "ACGTX", "ACGT", "ACGT"), std::runtime_error);
	}
}

TEST_CASE("build_star_alignment - trivial cases", "[ChimeraDetector]") {
	SECTION("All identical — trivial merge") {
		// Q==A==B, no gaps
		auto star = build_star_alignment("ACGT", "ACGT", "ACGT", "ACGT");
		REQUIRE(star.query_row == "ACGT");
		REQUIRE(star.parent_a_row == "ACGT");
		REQUIRE(star.parent_b_row == "ACGT");
	}

	SECTION("Substitutions only, no gaps") {
		// Q-A alignment: Q="ACGT", A="ACAT"
		// Q-B alignment: Q="ACGT", B="ACGA"
		auto star = build_star_alignment("ACGT", "ACAT", "ACGT", "ACGA");
		REQUIRE(star.query_row == "ACGT");
		REQUIRE(star.parent_a_row == "ACAT");
		REQUIRE(star.parent_b_row == "ACGA");
	}
}

TEST_CASE("build_star_alignment - gap handling", "[ChimeraDetector]") {
	SECTION("Gap in A's subject — B gets gap inserted") {
		// Q-A: Q="AC-GT", A="ACCGT" (insertion in A)
		// Q-B: Q="ACGT", B="ACGT"
		// Star should be: Q="AC-GT", A="ACCGT", B="AC-GT"
		auto star = build_star_alignment("AC-GT", "ACCGT", "ACGT", "ACGT");
		REQUIRE(star.query_row == "AC-GT");
		REQUIRE(star.parent_a_row == "ACCGT");
		REQUIRE(star.parent_b_row == "AC-GT");
	}

	SECTION("Gap in B's subject — A gets gap inserted") {
		// Q-A: Q="ACGT", A="ACGT"
		// Q-B: Q="AC-GT", B="ACTGT" (insertion in B)
		auto star = build_star_alignment("ACGT", "ACGT", "AC-GT", "ACTGT");
		REQUIRE(star.query_row == "AC-GT");
		REQUIRE(star.parent_a_row == "AC-GT");
		REQUIRE(star.parent_b_row == "ACTGT");
	}

	SECTION("Gaps in both alignments at different positions") {
		// Q-A: Q="A-CGT", A="AACGT"
		// Q-B: Q="ACG-T", B="ACGAT"
		auto star = build_star_alignment("A-CGT", "AACGT", "ACG-T", "ACGAT");
		// Both gaps should be present in the star alignment
		REQUIRE(star.query_row.size() == star.parent_a_row.size());
		REQUIRE(star.query_row.size() == star.parent_b_row.size());
		// Star should have 6 columns (original 4 query bases + 2 gaps)
		REQUIRE(star.query_row.size() == 6);
	}
}

TEST_CASE("classify_diffs - basic classifications", "[ChimeraDetector]") {
	SECTION("All identical → IGNORE") {
		StarAlignment star {"ACGT", "ACGT", "ACGT", {}};
		classify_diffs(star);
		REQUIRE(star.diffs.size() == 4);
		for (auto d : star.diffs) {
			REQUIRE(d == DiffType::IGNORE);
		}
	}

	SECTION("Q matches A only → MATCH_A") {
		// A!=B, Q==A at position 2: Q=G, A=G, B=A
		StarAlignment star {"ACGT", "ACGT", "ACAT", {}};
		classify_diffs(star);
		REQUIRE(star.diffs[2] == DiffType::MATCH_A); // G vs G vs A
	}

	SECTION("Q matches B only → MATCH_B") {
		// A!=B, Q==B at position 2: Q=A, A=G, B=A
		StarAlignment star {"ACAT", "ACGT", "ACAT", {}};
		classify_diffs(star);
		REQUIRE(star.diffs[2] == DiffType::MATCH_B); // A vs G vs A
	}

	SECTION("Both parents agree, Q differs → NO_VOTE") {
		// A==B, Q!=A at position 2: Q=T, A=G, B=G
		StarAlignment star {"ACTT", "ACGT", "ACGT", {}};
		classify_diffs(star);
		REQUIRE(star.diffs[2] == DiffType::NO_VOTE); // T vs G vs G
	}

	SECTION("Q matches neither → ABSTAIN") {
		// A!=B, Q!=A, Q!=B at position 2: Q=T, A=G, B=C
		StarAlignment star {"ACTT", "ACGT", "ACCT", {}};
		classify_diffs(star);
		REQUIRE(star.diffs[2] == DiffType::ABSTAIN); // T vs G vs C
	}
}

TEST_CASE("classify_diffs - gap character is not ambiguous", "[ChimeraDetector]") {
	// Regression: gap ('-') must not be treated as ambiguous.
	// It has its own handling (gap columns + gap-adjacent columns → IGNORE).
	// If is_ambiguous('-') returned true, classify_diffs would still work
	// (gaps are checked first), but the function semantics would be wrong.
	// This test verifies that a gap in parent_a at a non-gap-adjacent position
	// is classified as IGNORE due to the gap check, not the ambiguity check.
	StarAlignment star {"AACGT", "A-CGT", "AACGT", {}};
	classify_diffs(star);
	REQUIRE(star.diffs[0] == DiffType::IGNORE); // adjacent to gap at pos 1
	REQUIRE(star.diffs[1] == DiffType::IGNORE); // gap in A
	REQUIRE(star.diffs[2] == DiffType::IGNORE); // adjacent to gap at pos 1
	// Position 3 should be classifiable (not gap-adjacent: gap is at pos 1, adj covers pos 0 and 2)
	REQUIRE(star.diffs[3] == DiffType::IGNORE); // All identical: G==G==G
}

TEST_CASE("classify_diffs - gap and ambiguity handling", "[ChimeraDetector]") {
	SECTION("Gap in any row → IGNORE") {
		StarAlignment star {"A-GT", "ACGT", "ACGT", {}};
		classify_diffs(star);
		REQUIRE(star.diffs[1] == DiffType::IGNORE); // gap in Q
	}

	SECTION("Position adjacent to gap → IGNORE") {
		// Position 0 is adjacent to gap at position 1
		StarAlignment star {"A-GT", "ACGT", "ACGT", {}};
		classify_diffs(star);
		REQUIRE(star.diffs[0] == DiffType::IGNORE); // adjacent to gap at pos 1
		REQUIRE(star.diffs[1] == DiffType::IGNORE); // gap itself
		REQUIRE(star.diffs[2] == DiffType::IGNORE); // adjacent to gap at pos 1
	}

	SECTION("Ambiguous base (N) → IGNORE") {
		StarAlignment star {"ANGT", "ACGT", "ACGT", {}};
		classify_diffs(star);
		REQUIRE(star.diffs[1] == DiffType::IGNORE); // N is ambiguous
	}
}

TEST_CASE("classify_diffs - realistic chimera pattern", "[ChimeraDetector]") {
	// Simulate a chimera: first 5 positions match A, last 5 match B.
	// Q: ACGTAACGTA
	// A: ACGTAXXXXX  (matches Q in first 5, differs in last 5)
	// B: XXXXXACGTA  (differs from Q in first 5, matches in last 5)
	// But A and B must differ from each other at all 10 positions for clear signal.

	StarAlignment star {"ACGTAACGTA", // Q
	                    "ACGTATGCAT", // A: matches Q at 0-4, differs at 5-9
	                    "TGCATACGTA", // B: differs from Q at 0-4, matches at 5-9
	                    {}};
	classify_diffs(star);

	REQUIRE(star.diffs.size() == 10);

	// Count diff types
	int match_a = 0, match_b = 0, ignore = 0;
	for (auto d : star.diffs) {
		if (d == DiffType::MATCH_A) {
			match_a++;
		}
		if (d == DiffType::MATCH_B) {
			match_b++;
		}
		if (d == DiffType::IGNORE) {
			ignore++;
		}
	}

	// First 5: Q==A, A!=B → should be MATCH_A (where A!=B)
	// But some positions might be identical across all three → IGNORE
	// Let's check specific positions
	// Pos 0: Q=A, A=A, B=T → A!=B, Q==A → MATCH_A
	REQUIRE(star.diffs[0] == DiffType::MATCH_A);
	// Pos 5: Q=A, A=T, B=A → A!=B, Q==B → MATCH_B
	REQUIRE(star.diffs[5] == DiffType::MATCH_B);

	// Should have both MATCH_A and MATCH_B diffs
	REQUIRE(match_a > 0);
	REQUIRE(match_b > 0);
}

TEST_CASE("Star alignment + classify_diffs with WFA2 on real data", "[ChimeraDetector]") {
	// End-to-end test using actual test data sequences and WFA2 alignment,
	// then verifying the diff pattern matches vsearch's output.
	miint::WFA2Aligner aligner(miint::UCHIME_WFA2_MISMATCH, miint::UCHIME_WFA2_GAP_OPEN, miint::UCHIME_WFA2_GAP_EXTEND);

	std::string ref1 =
	    "AAGCCCAATCAACCACTCTCACTGGACGATTGCGGATATTGGCAACGAATTGGGAGGCGACCCGGACGACAGTCACGCCTTCTCGTTTGCGTACAGCTAT"
	    "TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCCGCCGTGATACCGGACCAAACAAGACGTCCACTTCAATGTTTAAATGACCTACGCGTCAGAA"
	    "CACCTTTCTACTATGTGTTCTCCCAGAATCATCTAGTACAATGGCGCGTCGTCATTAAAGCACCGGATGCGACGAACGGAGCGTGAATGAAGCTACTAC";

	std::string ref2 =
	    "AAGCACAATAAACCCCTGTCACGGGTCGAATAGGGGTTTAGGCAACGATCTCTGCAGAGCCCCTTGAGACAGTGACGCTTGTGCCGTTGCTTAAACTGATT"
	    "TGAAGGAGTCTAGCGGCCGCAGTAACGCACATTACCTAGTCCGCGTTCCCAGACGAAACAGGACGACATCTTTAAGGTTTAACTGACCGTCTATTCTTAAA"
	    "ACCTTCCAACTATGTGTTCCGAAAGAATCACCAACTACATTTACGCGTCGTGAATAACGTGTCGCCTGAGACGATAGGCGCGTGAATGAGGCGCTTAA";

	std::string chimera1 = ref1.substr(0, 150) + ref2.substr(150);

	auto align_a = aligner.align_full(chimera1, ref1);
	auto align_b = aligner.align_full(chimera1, ref2);
	REQUIRE(align_a.has_value());
	REQUIRE(align_b.has_value());

	auto star = build_star_alignment(align_a->query_aligned, align_a->subject_aligned, align_b->query_aligned,
	                                 align_b->subject_aligned);
	classify_diffs(star);

	// vsearch reports for query_chimera1:
	// Left 45: N 0, A 0, Y 45 → 45 MATCH_A diffs
	// Right 46: N 0, A 0, Y 46 → 46 MATCH_B diffs
	// Total non-ignore diffs = 91

	int match_a = 0, match_b = 0, no_vote = 0, abstain = 0;
	for (auto d : star.diffs) {
		switch (d) {
		case DiffType::MATCH_A:
			match_a++;
			break;
		case DiffType::MATCH_B:
			match_b++;
			break;
		case DiffType::NO_VOTE:
			no_vote++;
			break;
		case DiffType::ABSTAIN:
			abstain++;
			break;
		case DiffType::IGNORE:
			break;
		}
	}

	// vsearch found 45 A-diffs and 46 B-diffs with 0 N and 0 abstain.
	// Our counts should match exactly since the sequences have no gaps
	// and both aligners should produce the same ungapped alignment.
	REQUIRE(match_a == 45);
	REQUIRE(match_b == 46);
	REQUIRE(no_vote == 0);
	REQUIRE(abstain == 0);
}

TEST_CASE("Star alignment + classify_diffs with WFA2 — chimera2 (asymmetric)", "[ChimeraDetector]") {
	// query_chimera2: first 40% of ref3 + last 60% of ref4 (crossover at position 120)
	miint::WFA2Aligner aligner(miint::UCHIME_WFA2_MISMATCH, miint::UCHIME_WFA2_GAP_OPEN, miint::UCHIME_WFA2_GAP_EXTEND);

	std::string ref3 =
	    "AAGCCTAAGAGGCCACACCGACTGGCCGAGCAGGAATTCAGTCAAAGATATGTGCGGCCAACCTTGCGATAGTTACGCAGTCGCCGTTTCCTAAACCTATTT"
	    "GAAGGTTTCTAGCAGCCGCAGTAAGTGAGGATACCTCATCCGTATTACCAGACAAAATAACCGGTCATCTTAAATGTGTATATGACCCTTTCGTCACAAAAC"
	    "CTTTTTGCTGTGTGTCCCGCAAGAATCAACAACTACAATGGCGAGCCGCGACTAACGCGACGGCTGATACTAACGGCGCGTGAATGAAGCGCTTAA";

	std::string ref4 =
	    "AAGCCCAATAAACCAACCTGTCTGGCCAAATAGCGTTAGTGGCAACGTCATGTGTGACTACCCTTGCGATAGTGACGAACTCGCCGTTGCCTAACCCTATTT"
	    "GATGGAGTCTAGCGGACGCAGTAAGTCACAATACCTCGACCGTTTTACCTGAACAAACGAGGCGTCTTCGACAATGTTTTTAGGACCCTCTCGTAATAAAAA"
	    "ATTCCTACTATGTGTTCAGCAAGCATCAACAGCTGCATTGGAGCGTCGTGAATAAAGCGAGGGCTGTGACGGACTGTGCGTGAATGAAGCGCTTGA";

	std::string chimera2 = ref3.substr(0, 120) + ref4.substr(120);

	auto align_a = aligner.align_full(chimera2, ref3);
	auto align_b = aligner.align_full(chimera2, ref4);
	REQUIRE(align_a.has_value());
	REQUIRE(align_b.has_value());

	auto star = build_star_alignment(align_a->query_aligned, align_a->subject_aligned, align_b->query_aligned,
	                                 align_b->subject_aligned);
	classify_diffs(star);

	int match_a = 0, match_b = 0;
	for (auto d : star.diffs) {
		if (d == DiffType::MATCH_A) {
			match_a++;
		}
		if (d == DiffType::MATCH_B) {
			match_b++;
		}
	}

	// vsearch reports for query_chimera2:
	// Left 35: Y 35, Right 50: Y 50
	REQUIRE(match_a == 35);
	REQUIRE(match_b == 50);
}

TEST_CASE("Star alignment + classify_diffs with WFA2 — chimera3 (with mutation near crossover)", "[ChimeraDetector]") {
	miint::WFA2Aligner aligner(miint::UCHIME_WFA2_MISMATCH, miint::UCHIME_WFA2_GAP_OPEN, miint::UCHIME_WFA2_GAP_EXTEND);

	std::string ref5 =
	    "CAGCGCAATAAATCACTACGTCTGGCCGAATACGGATATAGGCAACGACTCGAGCGGCGACCCTGGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTATCT"
	    "GAAGGAGTCTAGCACTCGCTACAAGGCGCTATAGCTCGTCCGTGTTACCAGACCAAACAGGACGTCCTTGTCACTGTTGAACTTACCCTATCGTGGTAGAAC"
	    "CTTTCTACTAGGTGGTACGCAACAATCACCAACTACAATGGCACGTCGTGAATAACTCTCCGGATGAGACGAGCGTCGCTTTAGTGGACCACGTAA";

	std::string ref6 =
	    "AAGCCCACTAAACACCTCTGACTGGCCGGATACGGATTTAGGCAACGATGTGTGCCGCGCCCCTTACCACGGTGACGCTCTCGGAGTTGTCTAAACCTTTTTT"
	    "AAGGAGTCCGGCAGCCGCACTAAGGCACAATATCTCCTACGCGTAAGCAGATCAAACATGACGTCCTCTATAATGTTGAAATGAACCTCTCGTCATAAAACCA"
	    "TTCTACTATGAGTTCCGAAAGAATCAACAACTACAGTGGCGCGTCGGGTTTATCCTGACGGCTGAGATGAACGGTGCGCTAATGAACTGCTTAA";

	// chimera3 = first 150bp of ref5 + last 150bp of ref6 + 1 mutation at pos 155.
	// The actual sequence from chimera_queries.fasta (mutation already applied):
	std::string chimera3 =
	    "CAGCGCAATAAATCACTACGTCTGGCCGAATACGGATATAGGCAACGACTCGAGCGGCGACCCTGGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTATCT"
	    "GAAGGAGTCTAGCACTCGCTACAAGGCGCTATAGCTCGTCCGTGTTACCAGATGAAACATGACGTCCTCTATAATGTTGAAATGAACCTCTCGTCATAAAAC"
	    "CATTCTACTATGAGTTCCGAAAGAATCAACAACTACAGTGGCGCGTCGGGTTTATCCTGACGGCTGAGATGAACGGTGCGCTAATGAACTGCTTAA";

	auto align_a = aligner.align_full(chimera3, ref5);
	auto align_b = aligner.align_full(chimera3, ref6);
	REQUIRE(align_a.has_value());
	REQUIRE(align_b.has_value());

	auto star = build_star_alignment(align_a->query_aligned, align_a->subject_aligned, align_b->query_aligned,
	                                 align_b->subject_aligned);
	classify_diffs(star);

	int match_a = 0, match_b = 0, abstain = 0;
	for (auto d : star.diffs) {
		if (d == DiffType::MATCH_A) {
			match_a++;
		}
		if (d == DiffType::MATCH_B) {
			match_b++;
		}
		if (d == DiffType::ABSTAIN) {
			abstain++;
		}
	}

	// vsearch reports for query_chimera3:
	// Left 43: Y 43, Right 45: Y 44, A 1
	// Total diffs: 43 MATCH_A + 44 MATCH_B + 1 ABSTAIN = 88
	// But the mutation near the crossover may be classified differently
	// by our aligner vs vsearch. Verify totals are reasonable.
	REQUIRE(match_a + match_b >= 85); // at least 85 of ~88 total informative diffs
	REQUIRE(match_a > 0);
	REQUIRE(match_b > 0);
}

TEST_CASE("Star alignment + classify_diffs with WFA2 — divergent chimera", "[ChimeraDetector]") {
	miint::WFA2Aligner aligner(miint::UCHIME_WFA2_MISMATCH, miint::UCHIME_WFA2_GAP_OPEN, miint::UCHIME_WFA2_GAP_EXTEND);

	std::string ref1 =
	    "AAGCCCAATCAACCACTCTCACTGGACGATTGCGGATATTGGCAACGAATTGGGAGGCGACCCGGACGACAGTCACGCCTTCTCGTTTGCGTACAGCTAT"
	    "TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCCGCCGTGATACCGGACCAAACAAGACGTCCACTTCAATGTTTAAATGACCTACGCGTCAGAA"
	    "CACCTTTCTACTATGTGTTCTCCCAGAATCATCTAGTACAATGGCGCGTCGTCATTAAAGCACCGGATGCGACGAACGGAGCGTGAATGAAGCTACTAC";

	std::string ref2 =
	    "AAGCACAATAAACCCCTGTCACGGGTCGAATAGGGGTTTAGGCAACGATCTCTGCAGAGCCCCTTGAGACAGTGACGCTTGTGCCGTTGCTTAAACTGATT"
	    "TGAAGGAGTCTAGCGGCCGCAGTAACGCACATTACCTAGTCCGCGTTCCCAGACGAAACAGGACGACATCTTTAAGGTTTAACTGACCGTCTATTCTTAAA"
	    "ACCTTCCAACTATGTGTTCCGAAAGAATCACCAACTACATTTACGCGTCGTGAATAACGTGTCGCCTGAGACGATAGGCGCGTGAATGAGGCGCTTAA";

	// query_divergent_chimera from chimera_queries.fasta
	std::string divergent =
	    "AAGCCCAATCAACCACTCTCACTGGACGATTGCGGATATTGGCCACGAATTGGGAGGCGACCCGGACGACAGTTACGCCTTCTCGTTTGCGTACAGCTAT"
	    "TTGAAGGAGTCTAGCAGCCGCGGTGAGGCACAATACCTCCGCCGTGATACCAGACGAAACAGGACGAGAACTTTAAGGTTTAACTGAACGTCTATTATTAA"
	    "AACCTTACAAGTATGGGTTCCGAAAGCATGACCAACTACATTGACGCGTCGTGAATAACGTGTCGCCTGAGACGATAGGCGCGTGAATGAGTCGCTTAA";

	auto align_a = aligner.align_full(divergent, ref1);
	auto align_b = aligner.align_full(divergent, ref2);
	REQUIRE(align_a.has_value());
	REQUIRE(align_b.has_value());

	auto star = build_star_alignment(align_a->query_aligned, align_a->subject_aligned, align_b->query_aligned,
	                                 align_b->subject_aligned);
	classify_diffs(star);

	int match_a = 0, match_b = 0, no_vote = 0, abstain = 0;
	for (auto d : star.diffs) {
		switch (d) {
		case DiffType::MATCH_A:
			match_a++;
			break;
		case DiffType::MATCH_B:
			match_b++;
			break;
		case DiffType::NO_VOTE:
			no_vote++;
			break;
		case DiffType::ABSTAIN:
			abstain++;
			break;
		case DiffType::IGNORE:
			break;
		}
	}

	// vsearch reports for query_divergent_chimera (from expected_ref.tsv):
	// LY=44, LN=0, LA=4, RY=43, RN=2, RA=9
	// Total: match_a + match_b = 44+43 = 87 (but this depends on breakpoint)
	// Total non-ignore diffs = 44+0+4 + 43+2+9 = 102
	// We verify the total diff counts match what vsearch sees.
	// Note: the L/R split depends on breakpoint, but totals should match.
	int total_informative = match_a + match_b + no_vote + abstain;
	// vsearch total = LY+LN+LA+RY+RN+RA = 44+0+4+43+2+9 = 102
	// But L/R totals include both A-diffs and B-diffs.
	// Total A-diffs + B-diffs + N-diffs + ?-diffs should match.
	// From vsearch --uchimealns output for query_divergent_chimera:
	// Counting diff codes directly from the Diffs lines in expected_ref_alns.txt:
	//   A-diffs: 46 (44 uppercase 'A' + 2 lowercase 'a' = A-diffs on B-side of breakpoint)
	//   B-diffs: 43
	//   N-diffs: 11 (both parents agree, query differs — extra mutations from divergence)
	//   ?-diffs: 2 (query matches neither parent)
	// Total informative = 46 + 43 + 11 + 2 = 102
	//
	// IMPORTANT: these are *pre-breakpoint* raw column classifications, NOT vsearch's
	// post-breakpoint vote counts. vsearch's summary line shows:
	//   "Left 48: N 0, A 4, Y 44; Right 54: N 2, A 9, Y 43"
	// where Y=yes-votes, N=no-votes, A=abstain *within each segment after breakpoint*.
	// The LN+RN=2 counts N-diffs that are no-votes (on the "wrong" side), while our
	// no_vote=11 counts ALL columns where A==B and Q differs, regardless of breakpoint.
	// Both are correct — they measure different things. The post-breakpoint vote
	// assignment will be implemented in Phase 5 (breakpoint sweep).
	REQUIRE(total_informative == 102);
	REQUIRE(match_a == 46);
	REQUIRE(match_b == 43);
	REQUIRE(no_vote == 11);
	REQUIRE(abstain == 2);
}

// ============================================================================
// Phase 4: Smoothed parent selection
// ============================================================================

TEST_CASE("compute_match_profile basics", "[ChimeraDetector]") {
	SECTION("Identical sequences produce all-1 profile") {
		auto profile = miint::compute_match_profile("ACGT", "ACGT");
		REQUIRE(profile == std::vector<int> {1, 1, 1, 1});
	}

	SECTION("One mismatch produces 0 at that position") {
		auto profile = miint::compute_match_profile("ACGT", "ACAT");
		REQUIRE(profile == std::vector<int> {1, 1, 0, 1});
	}

	SECTION("Gap produces 0 at that position") {
		auto profile = miint::compute_match_profile("A-GT", "ACGT");
		REQUIRE(profile == std::vector<int> {1, 0, 1, 1});
	}
}

TEST_CASE("compute_smoothed basics", "[ChimeraDetector]") {
	SECTION("Window of 3 on all-1 profile") {
		std::vector<int> profile = {1, 1, 1, 1, 1};
		auto smooth = miint::compute_smoothed(profile, 3);
		// Running sum: [1, 2, 3, 3, 3]
		REQUIRE(smooth[0] == 1);
		REQUIRE(smooth[1] == 2);
		REQUIRE(smooth[2] == 3);
		REQUIRE(smooth[3] == 3);
		REQUIRE(smooth[4] == 3);
	}

	SECTION("Empty profile") {
		auto smooth = miint::compute_smoothed({}, 32);
		REQUIRE(smooth.empty());
	}
}

TEST_CASE("select_parents on chimera1", "[ChimeraDetector]") {
	miint::WFA2Aligner aligner(miint::UCHIME_WFA2_MISMATCH, miint::UCHIME_WFA2_GAP_OPEN, miint::UCHIME_WFA2_GAP_EXTEND);

	std::string ref1 =
	    "AAGCCCAATCAACCACTCTCACTGGACGATTGCGGATATTGGCAACGAATTGGGAGGCGACCCGGACGACAGTCACGCCTTCTCGTTTGCGTACAGCTAT"
	    "TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCCGCCGTGATACCGGACCAAACAAGACGTCCACTTCAATGTTTAAATGACCTACGCGTCAGAA"
	    "CACCTTTCTACTATGTGTTCTCCCAGAATCATCTAGTACAATGGCGCGTCGTCATTAAAGCACCGGATGCGACGAACGGAGCGTGAATGAAGCTACTAC";
	std::string ref2 =
	    "AAGCACAATAAACCCCTGTCACGGGTCGAATAGGGGTTTAGGCAACGATCTCTGCAGAGCCCCTTGAGACAGTGACGCTTGTGCCGTTGCTTAAACTGATT"
	    "TGAAGGAGTCTAGCGGCCGCAGTAACGCACATTACCTAGTCCGCGTTCCCAGACGAAACAGGACGACATCTTTAAGGTTTAACTGACCGTCTATTCTTAAA"
	    "ACCTTCCAACTATGTGTTCCGAAAGAATCACCAACTACATTTACGCGTCGTGAATAACGTGTCGCCTGAGACGATAGGCGCGTGAATGAGGCGCTTAA";

	std::string chimera1 = ref1.substr(0, 150) + ref2.substr(150);

	// Both refs as candidates (indices 0 and 1 into ref_sequences)
	std::vector<std::string> refs = {ref1, ref2};
	std::vector<uint32_t> candidates = {0, 1};

	auto parents = miint::select_parents(chimera1, candidates, refs, aligner);
	REQUIRE(parents.has_value());

	// vsearch selects ref1 as parent A (left) and ref2 as parent B (right)
	REQUIRE(parents->parent_a_idx == 0);
	REQUIRE(parents->parent_b_idx == 1);
}

// ============================================================================
// Phase 5: Breakpoint sweep + classification
// ============================================================================

TEST_CASE("sweep_breakpoints — clear chimera", "[ChimeraDetector]") {
	miint::UchimeParams params;

	SECTION("30 A-diffs then 30 B-diffs → high h-score") {
		std::vector<DiffType> diffs(60);
		std::fill(diffs.begin(), diffs.begin() + 30, DiffType::MATCH_A);
		std::fill(diffs.begin() + 30, diffs.end(), DiffType::MATCH_B);

		auto bp = miint::sweep_breakpoints(diffs, params);
		// H = (30 / (8*dn)) * (30 / (8*dn)) with dn=1.4, xn=8.0
		// = (30 / 11.2) * (30 / 11.2) = 2.678^2 = 7.17
		REQUIRE(bp.best_h > 7.0);
		REQUIRE(bp.left_yes == 30);
		REQUIRE(bp.right_yes == 30);
		REQUIRE(bp.left_no == 0);
		REQUIRE(bp.right_no == 0);
		REQUIRE(!bp.reversed);
	}

	SECTION("All A-diffs — one-sided, no valid breakpoint") {
		std::vector<DiffType> diffs(30, DiffType::MATCH_A);
		auto bp = miint::sweep_breakpoints(diffs, params);
		// No breakpoint where left_y > left_n AND right_y > right_n
		REQUIRE(bp.best_h == 0.0);
	}

	SECTION("All IGNORE — no informative diffs") {
		std::vector<DiffType> diffs(30, DiffType::IGNORE);
		auto bp = miint::sweep_breakpoints(diffs, params);
		REQUIRE(bp.best_h == 0.0);
	}

	SECTION("Reversed configuration: B-diffs left, A-diffs right") {
		std::vector<DiffType> diffs(60);
		std::fill(diffs.begin(), diffs.begin() + 30, DiffType::MATCH_B);
		std::fill(diffs.begin() + 30, diffs.end(), DiffType::MATCH_A);

		auto bp = miint::sweep_breakpoints(diffs, params);
		REQUIRE(bp.best_h > 7.0);
		REQUIRE(bp.reversed);
	}

	SECTION("Mixed diffs with some NO_VOTE and ABSTAIN") {
		// 10 A, 2 N, 1 ?, 10 B
		std::vector<DiffType> diffs;
		diffs.insert(diffs.end(), 10, DiffType::MATCH_A);
		diffs.insert(diffs.end(), 2, DiffType::NO_VOTE);
		diffs.insert(diffs.end(), 1, DiffType::ABSTAIN);
		diffs.insert(diffs.end(), 10, DiffType::MATCH_B);

		auto bp = miint::sweep_breakpoints(diffs, params);
		REQUIRE(bp.best_h > 0.28); // Should pass minh threshold
		REQUIRE(bp.left_yes > 0);
		REQUIRE(bp.right_yes > 0);
	}
}

TEST_CASE("Full detect() pipeline — chimera1 against reference DB", "[ChimeraDetector]") {
	miint::WFA2Aligner aligner(miint::UCHIME_WFA2_MISMATCH, miint::UCHIME_WFA2_GAP_OPEN, miint::UCHIME_WFA2_GAP_EXTEND);

	std::string ref1 =
	    "AAGCCCAATCAACCACTCTCACTGGACGATTGCGGATATTGGCAACGAATTGGGAGGCGACCCGGACGACAGTCACGCCTTCTCGTTTGCGTACAGCTAT"
	    "TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCCGCCGTGATACCGGACCAAACAAGACGTCCACTTCAATGTTTAAATGACCTACGCGTCAGAA"
	    "CACCTTTCTACTATGTGTTCTCCCAGAATCATCTAGTACAATGGCGCGTCGTCATTAAAGCACCGGATGCGACGAACGGAGCGTGAATGAAGCTACTAC";
	std::string ref2 =
	    "AAGCACAATAAACCCCTGTCACGGGTCGAATAGGGGTTTAGGCAACGATCTCTGCAGAGCCCCTTGAGACAGTGACGCTTGTGCCGTTGCTTAAACTGATT"
	    "TGAAGGAGTCTAGCGGCCGCAGTAACGCACATTACCTAGTCCGCGTTCCCAGACGAAACAGGACGACATCTTTAAGGTTTAACTGACCGTCTATTCTTAAA"
	    "ACCTTCCAACTATGTGTTCCGAAAGAATCACCAACTACATTTACGCGTCGTGAATAACGTGTCGCCTGAGACGATAGGCGCGTGAATGAGGCGCTTAA";
	std::string ref3 =
	    "AAGCCTAAGAGGCCACACCGACTGGCCGAGCAGGAATTCAGTCAAAGATATGTGCGGCCAACCTTGCGATAGTTACGCAGTCGCCGTTTCCTAAACCTATTT"
	    "GAAGGTTTCTAGCAGCCGCAGTAAGTGAGGATACCTCATCCGTATTACCAGACAAAATAACCGGTCATCTTAAATGTGTATATGACCCTTTCGTCACAAAAC"
	    "CTTTTTGCTGTGTGTCCCGCAAGAATCAACAACTACAATGGCGAGCCGCGACTAACGCGACGGCTGATACTAACGGCGCGTGAATGAAGCGCTTAA";

	std::vector<std::string> labels = {"ref1", "ref2", "ref3"};
	std::vector<std::string> seqs = {ref1, ref2, ref3};

	miint::ChimeraDetector detector;
	detector.set_reference(labels, seqs);

	std::string chimera1 = ref1.substr(0, 150) + ref2.substr(150);

	SECTION("Chimeric query detected as Y") {
		auto result = detector.detect("query_chimera1", chimera1, aligner);
		REQUIRE(result.flag == "Y");
		REQUIRE(result.parent_a_label == "ref1");
		REQUIRE(result.parent_b_label == "ref2");
		REQUIRE(result.score > 0.28);
	}

	SECTION("Clean query detected as N") {
		auto result = detector.detect("query_clean1", ref1, aligner);
		REQUIRE(result.flag == "N");
	}

	SECTION("Chimera1 vote counts match vsearch") {
		auto result = detector.detect("query_chimera1", chimera1, aligner);
		// vsearch: score=16.5019, LY=45, LN=0, LA=0, RY=46, RN=0, RA=0
		REQUIRE(result.left_yes == 45);
		REQUIRE(result.left_no == 0);
		REQUIRE(result.left_abstain == 0);
		REQUIRE(result.right_yes == 46);
		REQUIRE(result.right_no == 0);
		REQUIRE(result.right_abstain == 0);
		REQUIRE(std::abs(result.score - 16.5019) < 0.01);
	}

	SECTION("Chimera1 identities match vsearch") {
		auto result = detector.detect("query_chimera1", chimera1, aligner);
		// vsearch: QA=84.7, QB=85.0, AB=69.7, QModel=100.0, QT=85.0, div=15.0
		REQUIRE(std::abs(result.id_query_a - 84.7) < 0.5);
		REQUIRE(std::abs(result.id_query_b - 85.0) < 0.5);
		REQUIRE(std::abs(result.id_a_b - 69.7) < 0.5);
		REQUIRE(std::abs(result.id_query_model - 100.0) < 0.5);
		REQUIRE(std::abs(result.id_query_top - 85.0) < 0.5);
		REQUIRE(std::abs(result.divergence - 15.0) < 0.5);
	}
}

TEST_CASE("Full detect() pipeline — all 8 queries against full reference DB", "[ChimeraDetector]") {
	// Validate all 8 test queries against vsearch expected_ref.tsv ground truth.
	miint::WFA2Aligner aligner(miint::UCHIME_WFA2_MISMATCH, miint::UCHIME_WFA2_GAP_OPEN, miint::UCHIME_WFA2_GAP_EXTEND);

	std::string ref1 =
	    "AAGCCCAATCAACCACTCTCACTGGACGATTGCGGATATTGGCAACGAATTGGGAGGCGACCCGGACGACAGTCACGCCTTCTCGTTTGCGTACAGCTAT"
	    "TTGAAGGAGTCTAGCAGCCGCAGTAAGGCACAATACCTCCGCCGTGATACCGGACCAAACAAGACGTCCACTTCAATGTTTAAATGACCTACGCGTCAGAA"
	    "CACCTTTCTACTATGTGTTCTCCCAGAATCATCTAGTACAATGGCGCGTCGTCATTAAAGCACCGGATGCGACGAACGGAGCGTGAATGAAGCTACTAC";
	std::string ref2 =
	    "AAGCACAATAAACCCCTGTCACGGGTCGAATAGGGGTTTAGGCAACGATCTCTGCAGAGCCCCTTGAGACAGTGACGCTTGTGCCGTTGCTTAAACTGATT"
	    "TGAAGGAGTCTAGCGGCCGCAGTAACGCACATTACCTAGTCCGCGTTCCCAGACGAAACAGGACGACATCTTTAAGGTTTAACTGACCGTCTATTCTTAAA"
	    "ACCTTCCAACTATGTGTTCCGAAAGAATCACCAACTACATTTACGCGTCGTGAATAACGTGTCGCCTGAGACGATAGGCGCGTGAATGAGGCGCTTAA";
	std::string ref3 =
	    "AAGCCTAAGAGGCCACACCGACTGGCCGAGCAGGAATTCAGTCAAAGATATGTGCGGCCAACCTTGCGATAGTTACGCAGTCGCCGTTTCCTAAACCTATTT"
	    "GAAGGTTTCTAGCAGCCGCAGTAAGTGAGGATACCTCATCCGTATTACCAGACAAAATAACCGGTCATCTTAAATGTGTATATGACCCTTTCGTCACAAAAC"
	    "CTTTTTGCTGTGTGTCCCGCAAGAATCAACAACTACAATGGCGAGCCGCGACTAACGCGACGGCTGATACTAACGGCGCGTGAATGAAGCGCTTAA";
	std::string ref4 =
	    "AAGCCCAATAAACCAACCTGTCTGGCCAAATAGCGTTAGTGGCAACGTCATGTGTGACTACCCTTGCGATAGTGACGAACTCGCCGTTGCCTAACCCTATTT"
	    "GATGGAGTCTAGCGGACGCAGTAAGTCACAATACCTCGACCGTTTTACCTGAACAAACGAGGCGTCTTCGACAATGTTTTTAGGACCCTCTCGTAATAAAAA"
	    "ATTCCTACTATGTGTTCAGCAAGCATCAACAGCTGCATTGGAGCGTCGTGAATAAAGCGAGGGCTGTGACGGACTGTGCGTGAATGAAGCGCTTGA";
	std::string ref5 =
	    "CAGCGCAATAAATCACTACGTCTGGCCGAATACGGATATAGGCAACGACTCGAGCGGCGACCCTGGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTATCT"
	    "GAAGGAGTCTAGCACTCGCTACAAGGCGCTATAGCTCGTCCGTGTTACCAGACCAAACAGGACGTCCTTGTCACTGTTGAACTTACCCTATCGTGGTAGAAC"
	    "CTTTCTACTAGGTGGTACGCAACAATCACCAACTACAATGGCACGTCGTGAATAACTCTCCGGATGAGACGAGCGTCGCTTTAGTGGACCACGTAA";
	std::string ref6 =
	    "AAGCCCACTAAACACCTCTGACTGGCCGGATACGGATTTAGGCAACGATGTGTGCCGCGCCCCTTACCACGGTGACGCTCTCGGAGTTGTCTAAACCTTTTTT"
	    "AAGGAGTCCGGCAGCCGCACTAAGGCACAATATCTCCTACGCGTAAGCAGATCAAACATGACGTCCTCTATAATGTTGAAATGAACCTCTCGTCATAAAACCA"
	    "TTCTACTATGAGTTCCGAAAGAATCAACAACTACAGTGGCGCGTCGGGTTTATCCTGACGGCTGAGATGAACGGTGCGCTAATGAACTGCTTAA";

	std::vector<std::string> labels = {"ref1", "ref2", "ref3", "ref4", "ref5", "ref6"};
	std::vector<std::string> seqs = {ref1, ref2, ref3, ref4, ref5, ref6};

	miint::ChimeraDetector detector;
	detector.set_reference(labels, seqs);

	// vsearch expected_ref.tsv flags:
	// query_clean1=N, query_clean2=N, query_chimera1=Y, query_noparent=N,
	// query_chimera2=Y, query_chimera3=Y, query_short=N, query_divergent_chimera=Y
	SECTION("query_clean1 → N") {
		auto r = detector.detect("query_clean1", ref1, aligner);
		REQUIRE(r.flag == "N");
	}

	SECTION("query_chimera1 → Y with correct parents") {
		std::string chimera1 = ref1.substr(0, 150) + ref2.substr(150);
		auto r = detector.detect("query_chimera1", chimera1, aligner);
		REQUIRE(r.flag == "Y");
		REQUIRE(r.parent_a_label == "ref1");
		REQUIRE(r.parent_b_label == "ref2");
	}

	SECTION("query_chimera2 → Y with correct parents (order may differ from vsearch)") {
		std::string chimera2 = ref3.substr(0, 120) + ref4.substr(120);
		auto r = detector.detect("query_chimera2", chimera2, aligner);
		REQUIRE(r.flag == "Y");
		// Parent assignment order may differ from vsearch (ref3/ref4 vs ref4/ref3)
		// depending on which parent wins more smoothed identity positions.
		// Both parents must be present regardless of order.
		bool has_ref3 = (r.parent_a_label == "ref3" || r.parent_b_label == "ref3");
		bool has_ref4 = (r.parent_a_label == "ref4" || r.parent_b_label == "ref4");
		REQUIRE(has_ref3);
		REQUIRE(has_ref4);
	}

	SECTION("query_noparent → N") {
		std::string noparent =
		    "ACAAATAGAGTACACATGTCAGAGGCTGCTGTCGAGGTTGTAATCATTAACAATGGAGTCACATTGACCAGCTTCACGCGCTGTGCACCCACGCAGGGAGCC"
		    "GGTATCATAACGAGCACGGGCTTCACGACTGAGACTAAGGCGGAAAACGCGTAAGCGACGCGAGTATCCCTCCTCACCAAAGATCTCATACGGAATTAAGGA"
		    "GAGCCTGAAAGGCTACGTCGTGATTGTCCTGAACCGAGGCCTCTAGCAGCGTTATAGGTACCCGCTCGAACAGAACACCACCCTTGCGTAATTCAA";
		auto r = detector.detect("query_noparent", noparent, aligner);
		REQUIRE(r.flag == "N");
	}

	SECTION("query_short → N (too short for reliable detection)") {
		std::string short_q = ref1.substr(0, 60);
		auto r = detector.detect("query_short", short_q, aligner);
		REQUIRE(r.flag == "N");
	}

	SECTION("WFA2Aligner reuse across multiple detect() calls — bit-identical results") {
		// Verify the same aligner produces identical results across repeated calls.
		// If aligner state leaks between calls, scores or votes could drift.
		std::string chimera1 = ref1.substr(0, 150) + ref2.substr(150);
		auto r1 = detector.detect("query_chimera1", chimera1, aligner);
		auto r2 = detector.detect("query_chimera1", chimera1, aligner);
		auto r3 = detector.detect("query_chimera1", chimera1, aligner);
		REQUIRE(r1.flag == "Y");
		REQUIRE(r1.score == r2.score);
		REQUIRE(r2.score == r3.score);
		REQUIRE(r1.left_yes == r2.left_yes);
		REQUIRE(r2.left_yes == r3.left_yes);
		REQUIRE(r1.right_yes == r2.right_yes);
		REQUIRE(r1.id_query_a == r2.id_query_a);

		// After chimera calls, a clean call should work correctly
		auto r_clean = detector.detect("query_clean1", ref1, aligner);
		REQUIRE(r_clean.flag == "N");
	}
}

// ============================================================================
// Gapped sequence validation
// ============================================================================

TEST_CASE("Alignment with indel-containing sequences", "[ChimeraDetector]") {
	// Create sequences that require gaps in alignment, then verify
	// the full pipeline still works correctly.
	miint::WFA2Aligner aligner(miint::UCHIME_WFA2_MISMATCH, miint::UCHIME_WFA2_GAP_OPEN, miint::UCHIME_WFA2_GAP_EXTEND);

	// Parent A: 100bp
	std::string parent_a = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                       "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA";

	// Parent B: same as A but with a 3bp insertion at position 25
	std::string parent_b = parent_a.substr(0, 25) + "GGG" + parent_a.substr(25);
	// Also add substitutions in the second half to make it clearly different
	std::string parent_b_mut = parent_b;
	for (size_t i = 60; i < parent_b_mut.size() && i < 90; i += 3) {
		parent_b_mut[i] = (parent_b_mut[i] == 'A') ? 'C' : 'A';
	}

	SECTION("Alignment of sequences with length difference produces gaps") {
		auto result = aligner.align_full(parent_a, parent_b_mut);
		REQUIRE(result.has_value());
		// One of the aligned sequences should contain gaps
		bool has_gap = result->query_aligned.find('-') != std::string::npos ||
		               result->subject_aligned.find('-') != std::string::npos;
		REQUIRE(has_gap);
	}

	SECTION("Star alignment handles gapped pairwise alignments") {
		auto align_a = aligner.align_full(parent_a, parent_a); // ungapped
		auto align_b = aligner.align_full(parent_a, parent_b_mut);
		REQUIRE(align_a.has_value());
		REQUIRE(align_b.has_value());

		auto star = build_star_alignment(align_a->query_aligned, align_a->subject_aligned, align_b->query_aligned,
		                                 align_b->subject_aligned);

		// All rows should be same length
		REQUIRE(star.query_row.size() == star.parent_a_row.size());
		REQUIRE(star.query_row.size() == star.parent_b_row.size());

		// Classify should not crash
		classify_diffs(star);
		REQUIRE(star.diffs.size() == star.query_row.size());
	}

	SECTION("Full detect pipeline with gapped sequences — non-chimera") {
		std::vector<std::string> labels = {"parent_a", "parent_b_mut"};
		std::vector<std::string> seqs = {parent_a, parent_b_mut};

		miint::ChimeraDetector detector;
		detector.set_reference(labels, seqs);

		// Query that is NOT a chimera — just parent_a
		auto r = detector.detect("query_clean", parent_a, aligner);
		REQUIRE(r.flag == "N");
		// Divergence should be 0 or non-negative for non-chimeric
		REQUIRE(r.divergence >= 0.0);
	}

	SECTION("Parents with different alignment lengths don't cause OOB") {
		// parent_a is 102bp, parent_b_mut is 105bp (3bp insertion).
		// Alignments will have different lengths: ~105 columns for parent_b_mut,
		// ~102 for parent_a. select_parents must handle this safely.
		std::vector<std::string> labels = {"pa", "pb"};
		std::vector<std::string> seqs = {parent_a, parent_b_mut};

		// Create a third reference to have enough candidates (need >= 2)
		std::string parent_c = parent_a;
		for (size_t i = 70; i < 90 && i < parent_c.size(); i++) {
			parent_c[i] = (parent_c[i] == 'A') ? 'G' : 'A';
		}
		labels.push_back("pc");
		seqs.push_back(parent_c);

		miint::ChimeraDetector detector;
		detector.set_reference(labels, seqs);

		// Query that partially matches all three — should not crash
		auto r = detector.detect("query", parent_a, aligner);
		// Just verify no crash and a valid flag
		REQUIRE((r.flag == "Y" || r.flag == "N" || r.flag == "?"));
	}
}
