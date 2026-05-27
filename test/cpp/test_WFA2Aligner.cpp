#include "../../src/include/WFA2Aligner.hpp"
#include "../../src/include/sequence_utils.hpp"

#include <catch2/catch_test_macros.hpp>
#include <string>

using miint::WFA2Aligner;
using miint::WFA2CigarResult;
using miint::WFA2FullResult;

TEST_CASE("WFA2Aligner - Construction", "[WFA2Aligner]") {
	SECTION("Default construction succeeds") {
		REQUIRE_NOTHROW(WFA2Aligner());
	}
	SECTION("Custom penalties construction succeeds") {
		REQUIRE_NOTHROW(WFA2Aligner(2, 4, 1));
	}
}

TEST_CASE("WFA2Aligner - Validation", "[WFA2Aligner]") {
	SECTION("Rejects mismatch <= 0") {
		REQUIRE_THROWS_AS(WFA2Aligner(0, 6, 2), std::invalid_argument);
		REQUIRE_THROWS_AS(WFA2Aligner(-1, 6, 2), std::invalid_argument);
	}
	SECTION("Rejects gap_extend <= 0") {
		REQUIRE_THROWS_AS(WFA2Aligner(4, 6, 0), std::invalid_argument);
		REQUIRE_THROWS_AS(WFA2Aligner(4, 6, -1), std::invalid_argument);
	}
	SECTION("Rejects negative gap_open") {
		REQUIRE_THROWS_AS(WFA2Aligner(4, -1, 2), std::invalid_argument);
	}
	SECTION("Allows gap_open = 0") {
		REQUIRE_NOTHROW(WFA2Aligner(4, 0, 2));
	}
}

TEST_CASE("WFA2Aligner - align_score", "[WFA2Aligner]") {
	WFA2Aligner aligner;

	SECTION("Identical sequences -> 0") {
		auto result = aligner.align_score("ACGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 0);
	}
	SECTION("Single mismatch -> mismatch penalty (4)") {
		auto result = aligner.align_score("ACGT", "ACAT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 4);
	}
	SECTION("Single insertion -> gap_open + gap_extend (6+2=8)") {
		// query has extra base compared to subject
		auto result = aligner.align_score("ACGGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 8);
	}
	SECTION("Single deletion -> gap_open + gap_extend (6+2=8)") {
		auto result = aligner.align_score("ACGT", "ACGGT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 8);
	}
	SECTION("Custom penalties: mismatch=2") {
		WFA2Aligner custom(2, 6, 2);
		auto result = custom.align_score("ACGT", "ACAT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 2);
	}
}

TEST_CASE("WFA2Aligner - align_cigar", "[WFA2Aligner]") {
	WFA2Aligner aligner;

	SECTION("Identical sequences") {
		auto result = aligner.align_cigar("ACGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->score == 0);
		REQUIRE(result->cigar == "4=");
	}
	SECTION("Single mismatch") {
		auto result = aligner.align_cigar("ACGT", "ACAT");
		REQUIRE(result.has_value());
		REQUIRE(result->score == 4);
		// ACGT vs ACAT: positions 0,1 match, position 2 mismatches (G vs A), position 3 matches
		REQUIRE(result->cigar == "2=1X1=");
	}
	SECTION("CIGAR contains I for insertion") {
		auto result = aligner.align_cigar("ACGGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->cigar.find('I') != std::string::npos);
	}
	SECTION("CIGAR contains D for deletion") {
		auto result = aligner.align_cigar("ACGT", "ACGGT");
		REQUIRE(result.has_value());
		REQUIRE(result->cigar.find('D') != std::string::npos);
	}
}

TEST_CASE("WFA2Aligner - align_full", "[WFA2Aligner]") {
	WFA2Aligner aligner;

	SECTION("Identical sequences") {
		auto result = aligner.align_full("ACGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->score == 0);
		REQUIRE(result->cigar == "4=");
		REQUIRE(result->query_aligned == "ACGT");
		REQUIRE(result->subject_aligned == "ACGT");
	}
	SECTION("Single mismatch - no gaps in aligned seqs") {
		auto result = aligner.align_full("ACGT", "ACAT");
		REQUIRE(result.has_value());
		REQUIRE(result->query_aligned == "ACGT");
		REQUIRE(result->subject_aligned == "ACAT");
	}
	SECTION("Insertion - gap in subject_aligned") {
		auto result = aligner.align_full("ACGGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->subject_aligned.find('-') != std::string::npos);
		REQUIRE(result->query_aligned.size() == result->subject_aligned.size());
	}
	SECTION("Deletion - gap in query_aligned") {
		auto result = aligner.align_full("ACGT", "ACGGT");
		REQUIRE(result.has_value());
		REQUIRE(result->query_aligned.find('-') != std::string::npos);
		REQUIRE(result->query_aligned.size() == result->subject_aligned.size());
	}
}

TEST_CASE("WFA2Aligner - Reuse", "[WFA2Aligner]") {
	WFA2Aligner aligner;

	auto r1 = aligner.align_score("ACGT", "ACGT");
	auto r2 = aligner.align_score("ACGT", "ACAT");
	auto r3 = aligner.align_score("ACGT", "ACGT");

	REQUIRE(r1.has_value());
	REQUIRE(r2.has_value());
	REQUIRE(r3.has_value());
	REQUIRE(r1.value() == 0);
	REQUIRE(r2.value() == 4);
	REQUIRE(r3.value() == 0);
}

TEST_CASE("WFA2Aligner - Empty sequences", "[WFA2Aligner]") {
	WFA2Aligner aligner;

	SECTION("Both empty -> score 0") {
		auto result = aligner.align_score("", "");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 0);
	}
	SECTION("Both empty -> empty cigar") {
		auto result = aligner.align_cigar("", "");
		REQUIRE(result.has_value());
		REQUIRE(result->score == 0);
		REQUIRE(result->cigar.empty());
	}
	SECTION("Both empty -> full with empty aligned seqs") {
		auto result = aligner.align_full("", "");
		REQUIRE(result.has_value());
		REQUIRE(result->score == 0);
		REQUIRE(result->query_aligned.empty());
		REQUIRE(result->subject_aligned.empty());
	}
	SECTION("Query non-empty, subject empty -> all-insertion alignment") {
		// Score = gap_open(6) + 4 * gap_extend(2) = 14
		auto score = aligner.align_score("ACGT", "");
		REQUIRE(score.has_value());
		REQUIRE(score.value() == 14);

		auto cigar = aligner.align_cigar("ACGT", "");
		REQUIRE(cigar.has_value());
		REQUIRE(cigar->cigar == "4I");

		auto full = aligner.align_full("ACGT", "");
		REQUIRE(full.has_value());
		REQUIRE(full->query_aligned == "ACGT");
		REQUIRE(full->subject_aligned == "----");
	}
	SECTION("Query empty, subject non-empty -> all-deletion alignment") {
		auto score = aligner.align_score("", "ACGT");
		REQUIRE(score.has_value());
		REQUIRE(score.value() == 14);

		auto cigar = aligner.align_cigar("", "ACGT");
		REQUIRE(cigar.has_value());
		REQUIRE(cigar->cigar == "4D");

		auto full = aligner.align_full("", "ACGT");
		REQUIRE(full.has_value());
		REQUIRE(full->query_aligned == "----");
		REQUIRE(full->subject_aligned == "ACGT");
	}
}

TEST_CASE("WFA2Aligner - reconstruct_aligned edge cases via align_full", "[WFA2Aligner]") {
	WFA2Aligner aligner;

	SECTION("Multi-op CIGAR: insertion + deletion + matches") {
		// ACGGT vs ACGT: expect insertion in query
		auto result = aligner.align_full("ACGGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->query_aligned.size() == result->subject_aligned.size());
		// Subject should have a gap character
		REQUIRE(result->subject_aligned.find('-') != std::string::npos);
	}
	SECTION("Single character match") {
		auto result = aligner.align_full("A", "A");
		REQUIRE(result.has_value());
		REQUIRE(result->query_aligned == "A");
		REQUIRE(result->subject_aligned == "A");
	}
	SECTION("Single character mismatch") {
		auto result = aligner.align_full("A", "C");
		REQUIRE(result.has_value());
		REQUIRE(result->query_aligned == "A");
		REQUIRE(result->subject_aligned == "C");
	}
	SECTION("Aligned lengths always equal for various indel patterns") {
		// Multiple insertions
		auto r1 = aligner.align_full("AACCGGTT", "ACGT");
		REQUIRE(r1.has_value());
		REQUIRE(r1->query_aligned.size() == r1->subject_aligned.size());

		// Multiple deletions
		auto r2 = aligner.align_full("ACGT", "AACCGGTT");
		REQUIRE(r2.has_value());
		REQUIRE(r2->query_aligned.size() == r2->subject_aligned.size());
	}
}

TEST_CASE("WFA2Aligner - Realistic lengths", "[WFA2Aligner]") {
	WFA2Aligner aligner;

	SECTION("150bp pair with ~5% divergence") {
		std::string query = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";
		// Introduce ~5% mismatches (7-8 positions)
		std::string subject = query;
		subject[10] = 'T';
		subject[30] = 'T';
		subject[50] = 'T';
		subject[70] = 'T';
		subject[90] = 'T';
		subject[110] = 'T';
		subject[130] = 'T';

		auto result = aligner.align_score(query, subject);
		REQUIRE(result.has_value());
		REQUIRE(result.value() > 0);
	}
	SECTION("1000bp pair") {
		std::string query(1000, 'A');
		std::string subject(1000, 'A');
		subject[500] = 'C';

		auto result = aligner.align_score(query, subject);
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 4); // single mismatch with default penalty
	}
}

TEST_CASE("WFA2Aligner - Ambiguity codes (N)", "[WFA2Aligner]") {
	WFA2Aligner aligner;

	SECTION("N treated as mismatch against other bases") {
		auto result = aligner.align_score("ACGT", "ACNT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() > 0);
	}
	SECTION("N vs N treated as match") {
		auto result = aligner.align_score("ANGT", "ANGT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 0);
	}
}

TEST_CASE("WFA2Aligner - Edge cases", "[WFA2Aligner]") {
	WFA2Aligner aligner;

	SECTION("Single character sequences - identical") {
		auto result = aligner.align_score("A", "A");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 0);
	}
	SECTION("Single character sequences - different") {
		auto result = aligner.align_score("A", "C");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 4);
	}
	SECTION("All same character") {
		auto result = aligner.align_score("AAAA", "AAAA");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 0);
	}
}

// ---------------------------------------------------------------------------
// IUPAC bitmask utilities
// ---------------------------------------------------------------------------

TEST_CASE("IupacMatch - concrete bases", "[WFA2Aligner][iupac]") {
	REQUIRE(miint::IupacMatch('A', 'A'));
	REQUIRE(miint::IupacMatch('C', 'C'));
	REQUIRE(miint::IupacMatch('G', 'G'));
	REQUIRE(miint::IupacMatch('T', 'T'));
	REQUIRE_FALSE(miint::IupacMatch('A', 'C'));
	REQUIRE_FALSE(miint::IupacMatch('A', 'G'));
	REQUIRE_FALSE(miint::IupacMatch('A', 'T'));
}

TEST_CASE("IupacMatch - N matches everything", "[WFA2Aligner][iupac]") {
	REQUIRE(miint::IupacMatch('N', 'A'));
	REQUIRE(miint::IupacMatch('N', 'C'));
	REQUIRE(miint::IupacMatch('N', 'G'));
	REQUIRE(miint::IupacMatch('N', 'T'));
	REQUIRE(miint::IupacMatch('A', 'N'));
	REQUIRE(miint::IupacMatch('N', 'N'));
}

TEST_CASE("IupacMatch - two-base codes", "[WFA2Aligner][iupac]") {
	// R = A|G
	REQUIRE(miint::IupacMatch('R', 'A'));
	REQUIRE(miint::IupacMatch('R', 'G'));
	REQUIRE_FALSE(miint::IupacMatch('R', 'C'));
	REQUIRE_FALSE(miint::IupacMatch('R', 'T'));
	// Y = C|T
	REQUIRE(miint::IupacMatch('Y', 'C'));
	REQUIRE(miint::IupacMatch('Y', 'T'));
	REQUIRE_FALSE(miint::IupacMatch('Y', 'A'));
	REQUIRE_FALSE(miint::IupacMatch('Y', 'G'));
	// Symmetry
	REQUIRE(miint::IupacMatch('A', 'R'));
	REQUIRE(miint::IupacMatch('G', 'R'));
	// Two degenerate codes that overlap
	REQUIRE(miint::IupacMatch('R', 'M'));       // R=A|G, M=A|C → share A
	REQUIRE_FALSE(miint::IupacMatch('R', 'Y')); // R=A|G, Y=C|T → no overlap
}

TEST_CASE("IupacMatch - three-base codes", "[WFA2Aligner][iupac]") {
	// B = C|G|T (not A)
	REQUIRE_FALSE(miint::IupacMatch('B', 'A'));
	REQUIRE(miint::IupacMatch('B', 'C'));
	REQUIRE(miint::IupacMatch('B', 'G'));
	REQUIRE(miint::IupacMatch('B', 'T'));
}

TEST_CASE("IupacMatch - case insensitive", "[WFA2Aligner][iupac]") {
	REQUIRE(miint::IupacMatch('a', 'A'));
	REQUIRE(miint::IupacMatch('n', 'T'));
	REQUIRE(miint::IupacMatch('r', 'a'));
}

TEST_CASE("IupacMatch - invalid returns false", "[WFA2Aligner][iupac]") {
	REQUIRE_FALSE(miint::IupacMatch('Z', 'A'));
	REQUIRE_FALSE(miint::IupacMatch('A', '!'));
}

// ---------------------------------------------------------------------------
// IUPAC-aware semi-global alignment
// ---------------------------------------------------------------------------

TEST_CASE("align_cigar_semiglobal_iupac - ACGT-only matches literal", "[WFA2Aligner][iupac]") {
	WFA2Aligner aligner;
	// anchor = "ACGT", window = "XXACGTXX" — should find anchor in window
	auto literal = aligner.align_full_semiglobal("ACGT", "XXACGTXX");
	auto iupac = aligner.align_cigar_semiglobal_iupac("ACGT", "XXACGTXX");
	REQUIRE(literal.has_value());
	REQUIRE(iupac.has_value());
	REQUIRE(literal->cigar == iupac->cigar);
}

TEST_CASE("align_cigar_semiglobal_iupac - N in query matches any base", "[WFA2Aligner][iupac]") {
	WFA2Aligner aligner;
	// anchor = "ACNT" (N at pos 2), window = "ACGT"
	// N should match G at zero cost → pure match CIGAR
	auto result = aligner.align_cigar_semiglobal_iupac("ACNT", "ACGT");
	REQUIRE(result.has_value());
	REQUIRE(result->cigar == "4=");
}

TEST_CASE("align_cigar_semiglobal_iupac - R in query matches A or G", "[WFA2Aligner][iupac]") {
	WFA2Aligner aligner;
	// R = A|G → R vs A should match
	auto match = aligner.align_cigar_semiglobal_iupac("ACRT", "ACAT");
	REQUIRE(match.has_value());
	REQUIRE(match->cigar == "4=");

	// R vs C should mismatch
	auto mismatch = aligner.align_cigar_semiglobal_iupac("ACRT", "ACCT");
	REQUIRE(mismatch.has_value());
	REQUIRE(mismatch->cigar == "2=1X1=");
}

TEST_CASE("align_cigar_semiglobal_iupac - degenerate in subject too", "[WFA2Aligner][iupac]") {
	WFA2Aligner aligner;
	// Both sides degenerate: R vs M → R=A|G, M=A|C → share A → match
	auto result = aligner.align_cigar_semiglobal_iupac("R", "M");
	REQUIRE(result.has_value());
	REQUIRE(result->cigar == "1=");
}
