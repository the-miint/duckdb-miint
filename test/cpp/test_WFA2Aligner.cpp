#include "../../src/include/WFA2Aligner.hpp"

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
