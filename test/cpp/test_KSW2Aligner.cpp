#include "../../src/include/KSW2Aligner.hpp"

#include <catch2/catch_test_macros.hpp>
#include <string>

using miint::KSW2Aligner;
using miint::KSW2CigarResult;
using miint::KSW2FullResult;

// Scoring semantics for KSW2 (different from WFA2!):
//   identical 4-bp ACGT with default match=2  ->  score 8 (4 * match)
//   single mismatch (ACGT vs ACAT)            ->  score 8 - 2 - 4 = 2  (3 matches, minus one mismatch penalty)
//   single indel (gap_open=6 + 1*gap_extend=2) ->  matches recovery wins or loses depending on length
//   empty inputs                              ->  score 0
//   all-insertion of length 4                 ->  score -(gap_open + 4*gap_extend) = -14

TEST_CASE("KSW2Aligner - Construction", "[KSW2Aligner]") {
	SECTION("Default construction succeeds") {
		REQUIRE_NOTHROW(KSW2Aligner());
	}
	SECTION("Custom penalties (basic 4-arg) succeeds") {
		REQUIRE_NOTHROW(KSW2Aligner(2, 4, 6, 2));
	}
	SECTION("Advanced (6-arg, w + zdrop) succeeds") {
		REQUIRE_NOTHROW(KSW2Aligner(2, 4, 6, 2, -1, -1));
	}
	SECTION("Advanced with explicit band and zdrop") {
		REQUIRE_NOTHROW(KSW2Aligner(2, 4, 6, 2, 100, 400));
	}
}

TEST_CASE("KSW2Aligner - Validation", "[KSW2Aligner]") {
	SECTION("Rejects match <= 0") {
		REQUIRE_THROWS_AS(KSW2Aligner(0, 4, 6, 2), std::invalid_argument);
		REQUIRE_THROWS_AS(KSW2Aligner(-1, 4, 6, 2), std::invalid_argument);
	}
	SECTION("Rejects mismatch <= 0") {
		REQUIRE_THROWS_AS(KSW2Aligner(2, 0, 6, 2), std::invalid_argument);
		REQUIRE_THROWS_AS(KSW2Aligner(2, -1, 6, 2), std::invalid_argument);
	}
	SECTION("Rejects gap_extend <= 0") {
		REQUIRE_THROWS_AS(KSW2Aligner(2, 4, 6, 0), std::invalid_argument);
		REQUIRE_THROWS_AS(KSW2Aligner(2, 4, 6, -1), std::invalid_argument);
	}
	SECTION("Rejects negative gap_open") {
		REQUIRE_THROWS_AS(KSW2Aligner(2, 4, -1, 2), std::invalid_argument);
	}
	SECTION("Allows gap_open = 0") {
		REQUIRE_NOTHROW(KSW2Aligner(2, 4, 0, 2));
	}
}

TEST_CASE("KSW2Aligner - align_extz_score", "[KSW2Aligner]") {
	KSW2Aligner aligner;

	SECTION("Identical sequences -> 4 * match = 8") {
		auto result = aligner.align_extz_score("ACGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 8);
	}
	SECTION("Single mismatch -> 3*match - mismatch = 6 - 4 = 2") {
		auto result = aligner.align_extz_score("ACGT", "ACAT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 2);
	}
	SECTION("Single insertion -> 4*match - (gap_open + gap_extend) = 8 - 8 = 0") {
		// query ACGGT (5bp), subject ACGT (4bp): 4 matches, 1-base insertion
		auto result = aligner.align_extz_score("ACGGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 0);
	}
	SECTION("Single deletion -> 4*match - (gap_open + gap_extend) = 8 - 8 = 0") {
		auto result = aligner.align_extz_score("ACGT", "ACGGT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 0);
	}
	SECTION("Custom penalties: match=1") {
		KSW2Aligner custom(1, 4, 6, 2);
		auto result = custom.align_extz_score("ACGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 4); // 4 * match(1)
	}
}

TEST_CASE("KSW2Aligner - align_extz_cigar", "[KSW2Aligner]") {
	KSW2Aligner aligner;

	SECTION("Identical sequences -> 4M") {
		auto result = aligner.align_extz_cigar("ACGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->score == 8);
		REQUIRE(result->cigar == "4M");
	}
	SECTION("Single mismatch -> still 4M (KSW2 lumps match+mismatch)") {
		auto result = aligner.align_extz_cigar("ACGT", "ACAT");
		REQUIRE(result.has_value());
		REQUIRE(result->score == 2);
		REQUIRE(result->cigar == "4M");
	}
	SECTION("CIGAR contains I for insertion") {
		auto result = aligner.align_extz_cigar("ACGGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->cigar.find('I') != std::string::npos);
	}
	SECTION("CIGAR contains D for deletion") {
		auto result = aligner.align_extz_cigar("ACGT", "ACGGT");
		REQUIRE(result.has_value());
		REQUIRE(result->cigar.find('D') != std::string::npos);
	}
}

TEST_CASE("KSW2Aligner - align_extz_full", "[KSW2Aligner]") {
	KSW2Aligner aligner;

	SECTION("Identical sequences") {
		auto result = aligner.align_extz_full("ACGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->score == 8);
		REQUIRE(result->cigar == "4M");
		REQUIRE(result->query_aligned == "ACGT");
		REQUIRE(result->subject_aligned == "ACGT");
	}
	SECTION("Single mismatch - no gaps in aligned seqs") {
		auto result = aligner.align_extz_full("ACGT", "ACAT");
		REQUIRE(result.has_value());
		REQUIRE(result->query_aligned == "ACGT");
		REQUIRE(result->subject_aligned == "ACAT");
	}
	SECTION("Insertion - gap in subject_aligned") {
		auto result = aligner.align_extz_full("ACGGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->subject_aligned.find('-') != std::string::npos);
		REQUIRE(result->query_aligned.size() == result->subject_aligned.size());
	}
	SECTION("Deletion - gap in query_aligned") {
		auto result = aligner.align_extz_full("ACGT", "ACGGT");
		REQUIRE(result.has_value());
		REQUIRE(result->query_aligned.find('-') != std::string::npos);
		REQUIRE(result->query_aligned.size() == result->subject_aligned.size());
	}
}

TEST_CASE("KSW2Aligner - Reuse", "[KSW2Aligner]") {
	KSW2Aligner aligner;

	auto r1 = aligner.align_extz_score("ACGT", "ACGT");
	auto r2 = aligner.align_extz_score("ACGT", "ACAT");
	auto r3 = aligner.align_extz_score("ACGT", "ACGT");

	REQUIRE(r1.has_value());
	REQUIRE(r2.has_value());
	REQUIRE(r3.has_value());
	REQUIRE(r1.value() == 8);
	REQUIRE(r2.value() == 2);
	REQUIRE(r3.value() == 8);
}

TEST_CASE("KSW2Aligner - Empty sequences", "[KSW2Aligner]") {
	KSW2Aligner aligner;

	SECTION("Both empty -> score 0") {
		auto result = aligner.align_extz_score("", "");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 0);
	}
	SECTION("Both empty -> empty cigar") {
		auto result = aligner.align_extz_cigar("", "");
		REQUIRE(result.has_value());
		REQUIRE(result->score == 0);
		REQUIRE(result->cigar.empty());
	}
	SECTION("Both empty -> full with empty aligned seqs") {
		auto result = aligner.align_extz_full("", "");
		REQUIRE(result.has_value());
		REQUIRE(result->score == 0);
		REQUIRE(result->query_aligned.empty());
		REQUIRE(result->subject_aligned.empty());
	}
	SECTION("Query non-empty, subject empty -> all-insertion: -(gap_open + 4*gap_extend) = -14") {
		auto score = aligner.align_extz_score("ACGT", "");
		REQUIRE(score.has_value());
		REQUIRE(score.value() == -14);

		auto cigar = aligner.align_extz_cigar("ACGT", "");
		REQUIRE(cigar.has_value());
		REQUIRE(cigar->cigar == "4I");

		auto full = aligner.align_extz_full("ACGT", "");
		REQUIRE(full.has_value());
		REQUIRE(full->query_aligned == "ACGT");
		REQUIRE(full->subject_aligned == "----");
	}
	SECTION("Query empty, subject non-empty -> all-deletion: -14") {
		auto score = aligner.align_extz_score("", "ACGT");
		REQUIRE(score.has_value());
		REQUIRE(score.value() == -14);

		auto cigar = aligner.align_extz_cigar("", "ACGT");
		REQUIRE(cigar.has_value());
		REQUIRE(cigar->cigar == "4D");

		auto full = aligner.align_extz_full("", "ACGT");
		REQUIRE(full.has_value());
		REQUIRE(full->query_aligned == "----");
		REQUIRE(full->subject_aligned == "ACGT");
	}
}

TEST_CASE("KSW2Aligner - reconstruct edge cases via align_extz_full", "[KSW2Aligner]") {
	KSW2Aligner aligner;

	SECTION("Multi-op CIGAR: insertion + matches") {
		auto result = aligner.align_extz_full("ACGGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->query_aligned.size() == result->subject_aligned.size());
		REQUIRE(result->subject_aligned.find('-') != std::string::npos);
	}
	SECTION("Single character match") {
		auto result = aligner.align_extz_full("A", "A");
		REQUIRE(result.has_value());
		REQUIRE(result->query_aligned == "A");
		REQUIRE(result->subject_aligned == "A");
	}
	SECTION("Single character mismatch") {
		auto result = aligner.align_extz_full("A", "C");
		REQUIRE(result.has_value());
		REQUIRE(result->query_aligned == "A");
		REQUIRE(result->subject_aligned == "C");
	}
}

TEST_CASE("KSW2Aligner - Realistic lengths", "[KSW2Aligner]") {
	KSW2Aligner aligner;

	SECTION("150bp identical pair -> 150*match = 300") {
		std::string s(150, 'A');
		auto result = aligner.align_extz_score(s, s);
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 300);
	}
	SECTION("1000bp pair with single mismatch") {
		std::string query(1000, 'A');
		std::string subject(1000, 'A');
		subject[500] = 'C';
		auto result = aligner.align_extz_score(query, subject);
		REQUIRE(result.has_value());
		// 999 matches * 2 - 1 mismatch * 4 = 1998 - 4 = 1994
		REQUIRE(result.value() == 1994);
	}
}

TEST_CASE("KSW2Aligner - Ambiguity codes (N)", "[KSW2Aligner]") {
	KSW2Aligner aligner;

	SECTION("Sequences with N do not crash") {
		auto result = aligner.align_extz_score("ACGT", "ACNT");
		REQUIRE(result.has_value());
		// N is encoded as 4 in seq_nt4_table; with default match/mismatch matrix (no GENERIC_SC flag),
		// position-3 N vs T is treated under KSW2's fast path with wildcard-as-match. So we don't assert
		// exact value here; just that no crash and result is valid.
	}
}

TEST_CASE("KSW2Aligner - Edge cases", "[KSW2Aligner]") {
	KSW2Aligner aligner;

	SECTION("Single character sequences - identical -> match = 2") {
		auto result = aligner.align_extz_score("A", "A");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 2);
	}
	SECTION("Single character sequences - different -> -mismatch = -4") {
		auto result = aligner.align_extz_score("A", "C");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == -4);
	}
	SECTION("All same character -> 4*match") {
		auto result = aligner.align_extz_score("AAAA", "AAAA");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 8);
	}
}

// ---------------------------------------------------------------------------
// Dual-affine (extd) mode: ksw_extd2_sse.
// Param surface adds gap_open2 + gap_extend2 (the "second" gap-penalty pair).
// For any gap of length L, KSW2 chooses the cheaper of:
//   first affine:  gap_open  + L * gap_extend
//   second affine: gap_open2 + L * gap_extend2
// Typical use: first affine is cheap-open/expensive-extend (short indels);
//              second affine is expensive-open/cheap-extend (long indels).
// ---------------------------------------------------------------------------

TEST_CASE("KSW2Aligner - Dual-affine construction", "[KSW2Aligner]") {
	SECTION("8-arg dual-affine ctor succeeds") {
		REQUIRE_NOTHROW(KSW2Aligner(2, 4, 6, 2, 24, 1, -1, -1));
	}
	SECTION("Dual-affine with explicit bandwidth and zdrop") {
		REQUIRE_NOTHROW(KSW2Aligner(2, 4, 6, 2, 24, 1, 100, 400));
	}
}

TEST_CASE("KSW2Aligner - Dual-affine validation", "[KSW2Aligner]") {
	SECTION("Rejects negative gap_open2") {
		REQUIRE_THROWS_AS(KSW2Aligner(2, 4, 6, 2, -1, 1, -1, -1), std::invalid_argument);
	}
	SECTION("Rejects gap_extend2 <= 0") {
		REQUIRE_THROWS_AS(KSW2Aligner(2, 4, 6, 2, 24, 0, -1, -1), std::invalid_argument);
		REQUIRE_THROWS_AS(KSW2Aligner(2, 4, 6, 2, 24, -1, -1, -1), std::invalid_argument);
	}
	SECTION("Allows gap_open2 = 0") {
		REQUIRE_NOTHROW(KSW2Aligner(2, 4, 6, 2, 0, 1, -1, -1));
	}
	SECTION("Common-param validation still applies") {
		REQUIRE_THROWS_AS(KSW2Aligner(0, 4, 6, 2, 24, 1, -1, -1), std::invalid_argument); // match
		REQUIRE_THROWS_AS(KSW2Aligner(2, 0, 6, 2, 24, 1, -1, -1), std::invalid_argument); // mismatch
		REQUIRE_THROWS_AS(KSW2Aligner(2, 4, 6, 0, 24, 1, -1, -1), std::invalid_argument); // gap_extend
	}
}

TEST_CASE("KSW2Aligner - align_extd_score: identical and short gaps match extz", "[KSW2Aligner]") {
	// With short gaps, the first affine pair is cheaper than the second pair (typical
	// minimap2-style defaults), so dual-affine should score identically to extz.
	KSW2Aligner aligner(2, 4, 6, 2, 24, 1, -1, -1);

	SECTION("Identical sequences -> 4*match = 8") {
		auto result = aligner.align_extd_score("ACGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 8);
	}
	SECTION("Single mismatch -> 3*match - mismatch = 2") {
		auto result = aligner.align_extd_score("ACGT", "ACAT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 2);
	}
	SECTION("Single 1bp insertion -> first affine wins (cheaper for L=1)") {
		// Gap of 1: first = 6 + 1*2 = 8; second = 24 + 1*1 = 25. First wins.
		// 4 matches * 2 - 8 = 0.
		auto result = aligner.align_extd_score("ACGGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 0);
	}
}

TEST_CASE("KSW2Aligner - align_extd_score: long gap uses second affine", "[KSW2Aligner]") {
	// Sequences with a 30-base insertion in the query.
	// Common params: gap_open=6, gap_extend=2 (first); gap_open2=24, gap_extend2=1 (second).
	// First-affine cost  for 30bp gap: 6  + 30*2 = 66
	// Second-affine cost for 30bp gap: 24 + 30*1 = 54  (cheaper -> chosen)
	// 20 matches * 2 = 40; score under dual-affine = 40 - 54 = -14.
	// Compare with pure extz (no dual): would be 40 - 66 = -26.
	std::string subject(20, 'A');                       // 20 A's
	std::string query = subject + std::string(30, 'A'); // 50 A's = subject + 30 insertion

	KSW2Aligner dual_affine(2, 4, 6, 2, 24, 1, -1, -1);
	KSW2Aligner extz_only(2, 4, 6, 2, -1, -1);

	auto dual_score = dual_affine.align_extd_score(query, subject);
	auto extz_score = extz_only.align_extz_score(query, subject);

	REQUIRE(dual_score.has_value());
	REQUIRE(extz_score.has_value());
	REQUIRE(dual_score.value() == -14);
	REQUIRE(extz_score.value() == -26);
	REQUIRE(dual_score.value() > extz_score.value()); // dual is better (less negative)
}

TEST_CASE("KSW2Aligner - align_extd_cigar", "[KSW2Aligner]") {
	KSW2Aligner aligner(2, 4, 6, 2, 24, 1, -1, -1);

	SECTION("Identical sequences -> 4M") {
		auto result = aligner.align_extd_cigar("ACGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->score == 8);
		REQUIRE(result->cigar == "4M");
	}
	SECTION("Single mismatch -> still 4M (KSW2 lumps M)") {
		auto result = aligner.align_extd_cigar("ACGT", "ACAT");
		REQUIRE(result.has_value());
		REQUIRE(result->score == 2);
		REQUIRE(result->cigar == "4M");
	}
	SECTION("CIGAR contains I for insertion") {
		auto result = aligner.align_extd_cigar("ACGGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->cigar.find('I') != std::string::npos);
	}
}

TEST_CASE("KSW2Aligner - align_extd_full", "[KSW2Aligner]") {
	KSW2Aligner aligner(2, 4, 6, 2, 24, 1, -1, -1);

	SECTION("Identical sequences -> aligned seqs equal originals") {
		auto result = aligner.align_extd_full("ACGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->query_aligned == "ACGT");
		REQUIRE(result->subject_aligned == "ACGT");
	}
	SECTION("Insertion -> gap in subject_aligned") {
		auto result = aligner.align_extd_full("ACGGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->subject_aligned.find('-') != std::string::npos);
		REQUIRE(result->query_aligned.size() == result->subject_aligned.size());
	}
}

TEST_CASE("KSW2Aligner - align_extd empty inputs", "[KSW2Aligner]") {
	KSW2Aligner aligner(2, 4, 6, 2, 24, 1, -1, -1);

	SECTION("Both empty -> score 0, empty cigar") {
		auto s = aligner.align_extd_score("", "");
		REQUIRE(s.has_value());
		REQUIRE(s.value() == 0);
		auto c = aligner.align_extd_cigar("", "");
		REQUIRE(c.has_value());
		REQUIRE(c->cigar.empty());
	}
	SECTION("Query non-empty, subject empty -> all-insertion with cheaper affine for length 4") {
		// L=4: first = 6 + 4*2 = 14; second = 24 + 4*1 = 28. First wins. Score = -14.
		auto s = aligner.align_extd_score("ACGT", "");
		REQUIRE(s.has_value());
		REQUIRE(s.value() == -14);
		auto c = aligner.align_extd_cigar("ACGT", "");
		REQUIRE(c.has_value());
		REQUIRE(c->cigar == "4I");
	}
}

TEST_CASE("KSW2Aligner - extd and extz coexist on one instance", "[KSW2Aligner]") {
	// Same instance, both modes callable.
	KSW2Aligner aligner(2, 4, 6, 2, 24, 1, -1, -1);

	auto z = aligner.align_extz_score("ACGT", "ACGT");
	auto d = aligner.align_extd_score("ACGT", "ACGT");
	REQUIRE(z.has_value());
	REQUIRE(d.has_value());
	REQUIRE(z.value() == 8);
	REQUIRE(d.value() == 8);
}

// ---------------------------------------------------------------------------
// Splice-aware (exts) mode: ksw_exts2_sse.
// Param surface adds gap_open2 (intron-open penalty) + noncan (non-canonical splice
// penalty). No second extension penalty — introns extend at gap_extend rate.
// No `w` bandwidth — ksw_exts2_sse omits it. Junction array is passed as NULL (no
// guidance).
// ---------------------------------------------------------------------------

TEST_CASE("KSW2Aligner - Splice construction", "[KSW2Aligner]") {
	SECTION("7-arg splice ctor succeeds") {
		REQUIRE_NOTHROW(KSW2Aligner(2, 4, 6, 2, 24, 9, -1));
	}
	SECTION("Splice with explicit zdrop") {
		REQUIRE_NOTHROW(KSW2Aligner(2, 4, 6, 2, 24, 9, 400));
	}
}

TEST_CASE("KSW2Aligner - Splice validation", "[KSW2Aligner]") {
	SECTION("Rejects negative gap_open2") {
		REQUIRE_THROWS_AS(KSW2Aligner(2, 4, 6, 2, -1, 9, -1), std::invalid_argument);
	}
	SECTION("Rejects negative noncan") {
		REQUIRE_THROWS_AS(KSW2Aligner(2, 4, 6, 2, 24, -1, -1), std::invalid_argument);
	}
	SECTION("Allows noncan = 0") {
		REQUIRE_NOTHROW(KSW2Aligner(2, 4, 6, 2, 24, 0, -1));
	}
	SECTION("Common-param validation still applies") {
		REQUIRE_THROWS_AS(KSW2Aligner(0, 4, 6, 2, 24, 9, -1), std::invalid_argument);  // match
		REQUIRE_THROWS_AS(KSW2Aligner(2, 0, 6, 2, 24, 9, -1), std::invalid_argument);  // mismatch
		REQUIRE_THROWS_AS(KSW2Aligner(2, 4, 6, 0, 24, 9, -1), std::invalid_argument);  // gap_extend
		REQUIRE_THROWS_AS(KSW2Aligner(2, 4, -1, 2, 24, 9, -1), std::invalid_argument); // gap_open
	}
}

TEST_CASE("KSW2Aligner - align_exts_score: identical and short gaps behave like extz", "[KSW2Aligner]") {
	// For sequences without splice-pattern boundaries (GT-AG dinucleotides), exts gives the
	// same score as extz with the same gap_open / gap_extend.
	KSW2Aligner aligner(2, 4, 6, 2, 24, 9, -1);

	SECTION("Identical sequences -> 4*match = 8") {
		auto result = aligner.align_exts_score("ACGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 8);
	}
	SECTION("Single mismatch -> 3*match - mismatch = 2") {
		auto result = aligner.align_exts_score("ACGT", "ACAT");
		REQUIRE(result.has_value());
		REQUIRE(result.value() == 2);
	}
}

TEST_CASE("KSW2Aligner - align_exts_cigar", "[KSW2Aligner]") {
	KSW2Aligner aligner(2, 4, 6, 2, 24, 9, -1);

	SECTION("Identical -> 4M, score 8") {
		auto result = aligner.align_exts_cigar("ACGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->score == 8);
		REQUIRE(result->cigar == "4M");
	}
	SECTION("Mismatch -> 4M, score 2") {
		auto result = aligner.align_exts_cigar("ACGT", "ACAT");
		REQUIRE(result.has_value());
		REQUIRE(result->score == 2);
		REQUIRE(result->cigar == "4M");
	}
}

TEST_CASE("KSW2Aligner - align_exts_full reconstructs aligned strings", "[KSW2Aligner]") {
	KSW2Aligner aligner(2, 4, 6, 2, 24, 9, -1);

	SECTION("Identical") {
		auto result = aligner.align_exts_full("ACGT", "ACGT");
		REQUIRE(result.has_value());
		REQUIRE(result->query_aligned == "ACGT");
		REQUIRE(result->subject_aligned == "ACGT");
	}
}

TEST_CASE("KSW2Aligner - align_exts empty inputs", "[KSW2Aligner]") {
	KSW2Aligner aligner(2, 4, 6, 2, 24, 9, -1);

	SECTION("Both empty -> score 0") {
		auto s = aligner.align_exts_score("", "");
		REQUIRE(s.has_value());
		REQUIRE(s.value() == 0);
	}
	SECTION("Query empty -> all-deletion using first affine") {
		// Splice mode uses first affine for non-splice gaps in the empty-input case.
		auto s = aligner.align_exts_score("", "ACGT");
		REQUIRE(s.has_value());
		REQUIRE(s.value() == -14);
		auto c = aligner.align_exts_cigar("", "ACGT");
		REQUIRE(c.has_value());
		REQUIRE(c->cigar == "4D");
	}
}

TEST_CASE("KSW2Aligner - exts coexists with extz on same splice-ctor instance", "[KSW2Aligner]") {
	// Splice ctor stores gap_open2 + noncan. align_extz_* reads only gap_open/gap_extend (which
	// the splice ctor sets correctly) and ignores noncan, so it is safe to call on a splice
	// instance. align_extd_* is NOT safe here (gap_extend2 is a sentinel) and is intentionally
	// not exercised — see the contract comment in KSW2Aligner.hpp.
	KSW2Aligner aligner(2, 4, 6, 2, 24, 9, -1);

	auto z = aligner.align_extz_score("ACGT", "ACGT");
	auto s = aligner.align_exts_score("ACGT", "ACGT");
	REQUIRE(z.has_value());
	REQUIRE(s.has_value());
	REQUIRE(z.value() == 8);
	REQUIRE(s.value() == 8);
}
