#include <catch2/catch_test_macros.hpp>

#include "PileupWalker.hpp"
#include "alignment_functions_internal.hpp"

#include <cstdint>
#include <string>
#include <vector>

using namespace miint;

// Test fixture: 10-bp ref, common alignment vehicles.
static const std::string REF_10 = "ACGTACGTAC";
static const std::string REF_ID = "chr1";
static const std::string READ_ID = "r1";

static std::vector<PileupRow> walk(const std::string &cigar, std::int64_t start, const std::string &seq,
                                   const std::vector<std::uint8_t> &qual, bool qual_null = false) {
	std::vector<PileupRow> rows;
	const std::uint8_t *q = qual.empty() ? nullptr : qual.data();
	PileupWalker::Walk(cigar, REF_ID, REF_10, start, READ_ID, seq, q, qual.size(), qual_null, rows);
	return rows;
}

TEST_CASE("PileupWalker emits one row per ref pos on all-match", "[pileup][cigar]") {
	auto rows = walk("10=", 1, "ACGTACGTAC", {40, 40, 40, 40, 40, 40, 40, 40, 40, 40});
	REQUIRE(rows.size() == 10);
	for (std::size_t i = 0; i < rows.size(); ++i) {
		REQUIRE(rows[i].ref_pos == static_cast<std::int64_t>(i + 1));
		REQUIRE(rows[i].ref_base == REF_10[i]);
		REQUIRE(rows[i].query_base == REF_10[i]);
		REQUIRE(rows[i].query_qual == 40);
		REQUIRE_FALSE(rows[i].query_is_null);
		REQUIRE_FALSE(rows[i].qual_is_null);
		REQUIRE(rows[i].insert_pos == 0);
		REQUIRE_FALSE(rows[i].ref_base_is_null);
	}
}

TEST_CASE("PileupWalker emits X rows with mismatched query_base", "[pileup][cigar]") {
	auto rows = walk("3=1X6=", 1, "ACGGACGTAC", {40, 40, 40, 30, 40, 40, 40, 40, 40, 40});
	REQUIRE(rows.size() == 10);
	REQUIRE(rows[3].ref_base == 'T');
	REQUIRE(rows[3].query_base == 'G');
	REQUIRE(rows[3].query_qual == 30);
}

TEST_CASE("PileupWalker drops soft-clipped bases", "[pileup][cigar]") {
	// 3S7= — first 3 bases of query are clipped; ref_pos starts at the
	// alignment position regardless of clipping.
	auto rows = walk("3S7=", 1, "NNNACGTACG", {10, 10, 10, 40, 40, 40, 40, 40, 40, 40});
	REQUIRE(rows.size() == 7);
	REQUIRE(rows[0].ref_pos == 1);
	REQUIRE(rows[0].query_base == 'A');
	REQUIRE(rows[0].query_qual == 40);
}

TEST_CASE("PileupWalker ignores hard-clip ops entirely", "[pileup][cigar]") {
	// 3H7= — H consumes no query (SAM seq is pre-trimmed).
	auto rows = walk("3H7=", 1, "ACGTACG", {40, 40, 40, 40, 40, 40, 40});
	REQUIRE(rows.size() == 7);
	REQUIRE(rows[0].ref_pos == 1);
	REQUIRE(rows[0].query_base == 'A');
}

TEST_CASE("PileupWalker emits NULL query_base on D positions", "[pileup][cigar]") {
	// 3=2D5= — 2 reference positions skipped in query.
	auto rows = walk("3=2D5=", 1, "ACGCGTAC", {40, 40, 40, 40, 40, 40, 40, 40});
	REQUIRE(rows.size() == 10);
	REQUIRE(rows[3].query_is_null);
	REQUIRE(rows[4].query_is_null);
	REQUIRE(rows[3].qual_is_null);
	REQUIRE(rows[4].qual_is_null);
	REQUIRE(rows[3].ref_base == 'T');
	REQUIRE(rows[4].ref_base == 'A');
	for (const auto &r : rows) {
		REQUIRE(r.insert_pos == 0);
		REQUIRE_FALSE(r.ref_base_is_null);
	}
}

TEST_CASE("PileupWalker treats N like D (skipped reference span)", "[pileup][cigar]") {
	// 3=2N5= — same emit behavior as 3=2D5= in v1.
	auto rows = walk("3=2N5=", 1, "ACGCGTAC", {40, 40, 40, 40, 40, 40, 40, 40});
	REQUIRE(rows.size() == 10);
	REQUIRE(rows[3].query_is_null);
	REQUIRE(rows[4].query_is_null);
}

TEST_CASE("PileupWalker emits insertion rows for I ops", "[pileup][cigar]") {
	std::vector<PileupRow> rows;
	const std::string ref8 = "ACGTACGT";
	const std::string query = "ACGNNTACGT"; // 3 match + 2 insert + 5 match
	const std::vector<std::uint8_t> qual = {40, 40, 40, 5, 5, 40, 40, 40, 40, 40};
	PileupWalker::Walk("3=2I5=", REF_ID, ref8, 1, READ_ID, query, qual.data(), qual.size(), false, rows);
	REQUIRE(rows.size() == 10);

	// First 3 rows: ref-aligned match
	for (std::size_t i = 0; i < 3; ++i) {
		REQUIRE(rows[i].insert_pos == 0);
		REQUIRE_FALSE(rows[i].ref_base_is_null);
	}

	// Rows 3-4: insertion rows
	REQUIRE(rows[3].ref_pos == 3); // preceding ref position (SAM convention)
	REQUIRE(rows[3].insert_pos == 1);
	REQUIRE(rows[3].query_base == 'N');
	REQUIRE(rows[3].query_qual == 5);
	REQUIRE(rows[3].ref_base_is_null);
	REQUIRE_FALSE(rows[3].query_is_null);

	REQUIRE(rows[4].ref_pos == 3);
	REQUIRE(rows[4].insert_pos == 2);
	REQUIRE(rows[4].query_base == 'N');
	REQUIRE(rows[4].query_qual == 5);
	REQUIRE(rows[4].ref_base_is_null);

	// Rows 5-9: ref-aligned match (post-insertion)
	for (std::size_t i = 5; i < 10; ++i) {
		REQUIRE(rows[i].insert_pos == 0);
		REQUIRE(rows[i].query_qual == 40);
		REQUIRE_FALSE(rows[i].ref_base_is_null);
	}
}

TEST_CASE("PileupWalker emits insertion rows for leading insertion", "[pileup][cigar]") {
	std::vector<PileupRow> rows;
	const std::string ref8 = "ACGTACGT";
	const std::string query = "NNACGTACGT"; // 2 insert + 8 match
	const std::vector<std::uint8_t> qual = {5, 5, 40, 40, 40, 40, 40, 40, 40, 40};
	PileupWalker::Walk("2I8=", REF_ID, ref8, 1, READ_ID, query, qual.data(), qual.size(), false, rows);
	REQUIRE(rows.size() == 10);

	// First 2 rows: insertion at ref_pos 0 (before align_start=1)
	REQUIRE(rows[0].ref_pos == 0);
	REQUIRE(rows[0].insert_pos == 1);
	REQUIRE(rows[0].query_base == 'N');
	REQUIRE(rows[0].ref_base_is_null);
	REQUIRE(rows[1].ref_pos == 0);
	REQUIRE(rows[1].insert_pos == 2);

	// Remaining 8 rows: ref-aligned
	for (std::size_t i = 2; i < 10; ++i) {
		REQUIRE(rows[i].insert_pos == 0);
		REQUIRE(rows[i].ref_pos == static_cast<std::int64_t>(i - 1));
	}
}

TEST_CASE("PileupWalker emits insertion rows for trailing insertion", "[pileup][cigar]") {
	std::vector<PileupRow> rows;
	const std::string ref8 = "ACGTACGT";
	const std::string query = "ACGTACGTNN"; // 8 match + 2 insert
	const std::vector<std::uint8_t> qual = {40, 40, 40, 40, 40, 40, 40, 40, 5, 5};
	PileupWalker::Walk("8=2I", REF_ID, ref8, 1, READ_ID, query, qual.data(), qual.size(), false, rows);
	REQUIRE(rows.size() == 10);

	// First 8 rows: ref-aligned
	for (std::size_t i = 0; i < 8; ++i) {
		REQUIRE(rows[i].insert_pos == 0);
	}

	// Last 2 rows: insertion at ref_pos 8 (last consumed ref position)
	REQUIRE(rows[8].ref_pos == 8);
	REQUIRE(rows[8].insert_pos == 1);
	REQUIRE(rows[8].query_base == 'N');
	REQUIRE(rows[8].query_qual == 5);
	REQUIRE(rows[8].ref_base_is_null);
	REQUIRE(rows[9].ref_pos == 8);
	REQUIRE(rows[9].insert_pos == 2);
}

TEST_CASE("PileupWalker handles 1-based ref position correctly mid-reference", "[pileup][cigar]") {
	auto rows = walk("5=", 3, "GTACG", {40, 40, 40, 40, 40});
	REQUIRE(rows.size() == 5);
	REQUIRE(rows[0].ref_pos == 3);
	REQUIRE(rows[0].ref_base == 'G');
	REQUIRE(rows[4].ref_pos == 7);
	REQUIRE(rows[4].ref_base == 'G');
}

TEST_CASE("PileupWalker propagates qual_is_null to every emitted row", "[pileup][cigar]") {
	std::vector<PileupRow> rows;
	PileupWalker::Walk("3=", REF_ID, REF_10, 1, READ_ID, "ACG", /*qual_data=*/nullptr, /*qual_len=*/0,
	                   /*qual_is_null=*/true, rows);
	REQUIRE(rows.size() == 3);
	for (const auto &r : rows) {
		REQUIRE(r.qual_is_null);
		REQUIRE_FALSE(r.query_is_null);
		REQUIRE(r.query_qual == 0);
	}
}

TEST_CASE("PileupWalker emits zero rows on empty / unmapped CIGAR", "[pileup][cigar]") {
	auto rows = walk("*", 0, "", {});
	REQUIRE(rows.empty());
	auto rows2 = walk("", 0, "", {});
	REQUIRE(rows2.empty());
}

TEST_CASE("PileupWalker rejects alignment past reference end", "[pileup][cigar]") {
	// 5= starting at pos 8 needs positions 8..12 but ref is 10 bp.
	REQUIRE_THROWS_AS(walk("5=", 8, "ACGTA", {40, 40, 40, 40, 40}), miint::InvalidInputException);
}

TEST_CASE("PileupWalker rejects CIGAR/seq length mismatch", "[pileup][cigar]") {
	// 5= demands 5 query bases; seq has 3.
	REQUIRE_THROWS_AS(walk("5=", 1, "ACG", {40, 40, 40}), miint::InvalidInputException);
}

TEST_CASE("PileupWalker rejects qual length mismatch when qual present", "[pileup][cigar]") {
	REQUIRE_THROWS_AS(walk("5=", 1, "ACGTA", {40, 40, 40}), miint::InvalidInputException);
}

TEST_CASE("PileupWalker rejects unsupported CIGAR op", "[pileup][cigar]") {
	// CIGAR with '?' is rejected at the parse step; the walker validation
	// catches non-standard ops via ParseCigarOperations throwing.
	REQUIRE_THROWS_AS(walk("3?", 1, "ACG", {40, 40, 40}), miint::InvalidInputException);
}
