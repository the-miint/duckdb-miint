#include <catch2/catch_test_macros.hpp>

#include "MsaConsensus.hpp"
#include "alignment_functions_internal.hpp"

#include <cstdint>
#include <string>
#include <vector>

using namespace miint;

// ---------------------------------------------------------------------------
// VoteColumn: per-column 5-state log-likelihood reduction
// ---------------------------------------------------------------------------

TEST_CASE("VoteColumn: single high-Q observation picks the observed base", "[consensus][vote]") {
	std::vector<Observation> obs = {{'A', 40}};
	auto v = VoteColumn(obs);
	REQUIRE(v.base == 'A');
	REQUIRE(v.qual > 0);
}

TEST_CASE("VoteColumn: all-same high-Q observations cap qual at 60", "[consensus][vote]") {
	std::vector<Observation> obs = {{'A', 40}, {'A', 40}, {'A', 40}, {'A', 40}, {'A', 40}};
	auto v = VoteColumn(obs);
	REQUIRE(v.base == 'A');
	REQUIRE(v.qual == 60);
}

TEST_CASE("VoteColumn: higher-Q base beats lower-Q disagreement", "[consensus][vote]") {
	// A@40 vs C@20 → A wins
	std::vector<Observation> obs = {{'A', 40}, {'C', 20}};
	auto v = VoteColumn(obs);
	REQUIRE(v.base == 'A');
}

TEST_CASE("VoteColumn: majority gap suppresses the column", "[consensus][vote]") {
	// 4 gaps, 1 A@40 → gap wins
	std::vector<Observation> obs = {{'-', 0}, {'-', 0}, {'-', 0}, {'-', 0}, {'A', 40}};
	auto v = VoteColumn(obs);
	REQUIRE(v.base == '-');
}

TEST_CASE("VoteColumn: minority gap with strong base evidence keeps the base", "[consensus][vote]") {
	// 1 gap, 4 A@40 → A wins clearly
	std::vector<Observation> obs = {{'-', 0}, {'A', 40}, {'A', 40}, {'A', 40}, {'A', 40}};
	auto v = VoteColumn(obs);
	REQUIRE(v.base == 'A');
}

TEST_CASE("VoteColumn: equal-evidence tie breaks alphabetically", "[consensus][vote]") {
	// 2 A@30, 2 C@30 → exactly tied log-likelihood; A wins by alphabet
	std::vector<Observation> obs = {{'A', 30}, {'A', 30}, {'C', 30}, {'C', 30}};
	auto v = VoteColumn(obs);
	REQUIRE(v.base == 'A');
}

TEST_CASE("VoteColumn: posterior Q reflects mixed-evidence uncertainty", "[consensus][vote]") {
	// A@40 vs C@20: A wins but with reduced confidence vs all-agree case
	std::vector<Observation> all_agree = {{'A', 40}, {'A', 40}};
	std::vector<Observation> conflict = {{'A', 40}, {'C', 20}};
	auto v_all = VoteColumn(all_agree);
	auto v_mix = VoteColumn(conflict);
	REQUIRE(v_all.base == 'A');
	REQUIRE(v_mix.base == 'A');
	REQUIRE(v_mix.qual < v_all.qual);
}

TEST_CASE("VoteColumn: empty observation list throws", "[consensus][vote]") {
	std::vector<Observation> obs;
	REQUIRE_THROWS_AS(VoteColumn(obs), miint::InvalidInputException);
}

// ---------------------------------------------------------------------------
// HpCorrect: median-per-read homopolymer-length correction
// ---------------------------------------------------------------------------

TEST_CASE("HpCorrect: HP run with consistent length is preserved", "[consensus][hp]") {
	// Naive consensus = "TAAAAAAT" (6 A's between T's). All reads have 6 A's
	// ungapped. HP correction should keep length 6.
	std::string cons = "TAAAAAAT";
	std::vector<std::uint8_t> q(cons.size(), 60);
	std::vector<std::size_t> msa_col_of_cons_pos = {0, 1, 2, 3, 4, 5, 6, 7};
	std::vector<std::string> aligned = {"TAAAAAAT", "TAAAAAAT", "TAAAAAAT", "TAAAAAAT", "TAAAAAAT"};
	auto out = HpCorrect(cons, q, msa_col_of_cons_pos, aligned);
	REQUIRE(out.first == "TAAAAAAT");
	REQUIRE(out.second.size() == 8);
}

TEST_CASE("HpCorrect: column-vote overcounts HP length, median per-read wins", "[consensus][hp]") {
	// Reads: 2× "TAAAA-T" (4 A's) + 3× shifted-gap patterns (4 A's each).
	// Per-column vote → "TAAAAAT" (5 A's). HP correction → "TAAAAT" (4 A's).
	std::string cons = "TAAAAAT";
	std::vector<std::uint8_t> q(cons.size(), 60);
	std::vector<std::size_t> msa_col_of_cons_pos = {0, 1, 2, 3, 4, 5, 6};
	std::vector<std::string> aligned = {"TAAAA-T", "TAAAA-T", "TAA-AAT", "TA-AAAT", "T-AAAAT"};
	auto out = HpCorrect(cons, q, msa_col_of_cons_pos, aligned);
	REQUIRE(out.first == "TAAAAT");
	REQUIRE(out.second.size() == out.first.size());
}

TEST_CASE("HpCorrect: single-base 'runs' are not corrected", "[consensus][hp]") {
	// "ACGT" — every char unique; HP correction is a no-op.
	std::string cons = "ACGT";
	std::vector<std::uint8_t> q = {60, 60, 60, 60};
	std::vector<std::size_t> msa_col_of_cons_pos = {0, 1, 2, 3};
	std::vector<std::string> aligned = {"ACGT", "ACGT"};
	auto out = HpCorrect(cons, q, msa_col_of_cons_pos, aligned);
	REQUIRE(out.first == "ACGT");
}

TEST_CASE("HpCorrect: HP correction preserves the qual scalar across the run", "[consensus][hp]") {
	// Single 4-A consensus run; all reads have 4 A's → run length stays 4.
	// Qual at every emitted run position should equal the input qual at the
	// run's first position (we don't synthesise new posterior Qs in the HP pass).
	std::string cons = "AAAA";
	std::vector<std::uint8_t> q = {55, 55, 55, 55};
	std::vector<std::size_t> msa_col_of_cons_pos = {0, 1, 2, 3};
	std::vector<std::string> aligned = {"AAAA", "AAAA"};
	auto out = HpCorrect(cons, q, msa_col_of_cons_pos, aligned);
	REQUIRE(out.first == "AAAA");
	REQUIRE(out.second == std::vector<std::uint8_t>({55, 55, 55, 55}));
}

TEST_CASE("HpCorrect: HP shorter than column-vote when most reads agree it's shorter", "[consensus][hp]") {
	// Consensus "AAAA" (column vote), but 4 of 5 reads actually have 3 A's
	// ungapped. Median = 3 → correct down to "AAA".
	std::string cons = "AAAA";
	std::vector<std::uint8_t> q = {50, 50, 50, 50};
	std::vector<std::size_t> msa_col_of_cons_pos = {0, 1, 2, 3};
	std::vector<std::string> aligned = {"AAA-", "AAA-", "AAA-", "AAA-", "AAAA"};
	auto out = HpCorrect(cons, q, msa_col_of_cons_pos, aligned);
	REQUIRE(out.first == "AAA");
}

// ---------------------------------------------------------------------------
// VoteColumn: boundary behavior at the kGapErrProb tuning point. Pinning
// these so the constant doesn't silently drift on future tweaks.
// ---------------------------------------------------------------------------

TEST_CASE("VoteColumn: 3 base + 2 gap keeps the base", "[consensus][vote][gap-boundary]") {
	// At kGapErrProb=0.05 with 3 strong base obs vs 2 gap obs, the base
	// should win cleanly — gap evidence is real but minority and per-gap
	// evidence is weaker than per-Q40 base evidence.
	std::vector<Observation> obs = {{'A', 40}, {'A', 40}, {'A', 40}, {'-', 0}, {'-', 0}};
	auto v = VoteColumn(obs);
	REQUIRE(v.base == 'A');
}

TEST_CASE("VoteColumn: 2 base + 3 gap still keeps the base at Q=40", "[consensus][vote][gap-boundary]") {
	// 2 strong base obs vs 3 gaps: with Q=40 (very confident) base evidence
	// outweighs majority-gap evidence. If kGapErrProb is raised much past
	// ~0.10 this test will start to fail — guards the boundary.
	std::vector<Observation> obs = {{'A', 40}, {'A', 40}, {'-', 0}, {'-', 0}, {'-', 0}};
	auto v = VoteColumn(obs);
	REQUIRE(v.base == 'A');
}

TEST_CASE("VoteColumn: 1 base + 4 gap suppresses the column", "[consensus][vote][gap-boundary]") {
	// 1 strong base obs vs 4 gaps: gap wins. Already covered by an earlier
	// test but kept here in the gap-boundary group for cohesion.
	std::vector<Observation> obs = {{'A', 40}, {'-', 0}, {'-', 0}, {'-', 0}, {'-', 0}};
	auto v = VoteColumn(obs);
	REQUIRE(v.base == '-');
}

// ---------------------------------------------------------------------------
// BuildConsensus: partition-invariance test (Combine semantics).
//
// DuckDB parallelises GROUP BY aggregates by accumulating partial states
// per thread, then merging via Combine. We need the result to be invariant
// under that partitioning. The aggregate's Combine just concatenates the
// per-row vectors; if BuildConsensus is order-insensitive, the merged-state
// path yields the same consensus as the single-batch path.
// ---------------------------------------------------------------------------

TEST_CASE("BuildConsensus is invariant under row partitioning (Combine semantics)", "[consensus][build][combine]") {
	std::vector<std::string> all_seqs = {"ACGT", "ACGT", "ACAT", "ACGT", "ACGT"};
	std::vector<std::vector<std::uint8_t>> all_quals(5, {40, 40, 40, 40});

	const auto full = BuildConsensus(all_seqs, all_quals);

	// Build two partials, then concat them (Combine's semantics: target =
	// target ++ source). Round-trip through BuildConsensus and compare.
	std::vector<std::string> partA(all_seqs.begin(), all_seqs.begin() + 2);
	std::vector<std::vector<std::uint8_t>> qualA(all_quals.begin(), all_quals.begin() + 2);
	std::vector<std::string> partB(all_seqs.begin() + 2, all_seqs.end());
	std::vector<std::vector<std::uint8_t>> qualB(all_quals.begin() + 2, all_quals.end());

	partA.insert(partA.end(), partB.begin(), partB.end());
	qualA.insert(qualA.end(), qualB.begin(), qualB.end());

	const auto combined = BuildConsensus(partA, qualA);
	REQUIRE(full.first == combined.first);
	REQUIRE(full.second == combined.second);
}

TEST_CASE("BuildConsensus is invariant under partitioning with HP correction in play",
          "[consensus][build][combine][hp]") {
	// HP-correction depends on the per-row aligned_seqs vector. If Combine
	// concatenates in the wrong order or drops rows, the HP median changes
	// and we'd see "TAAAAAT" instead of "TAAAAT".
	std::vector<std::string> all_seqs = {"TAAAA-T", "TAAAA-T", "TAA-AAT", "TA-AAAT", "T-AAAAT"};
	std::vector<std::vector<std::uint8_t>> all_quals(5, {40, 40, 40, 40, 40, 40});

	const auto full = BuildConsensus(all_seqs, all_quals);

	std::vector<std::string> partA(all_seqs.begin(), all_seqs.begin() + 3);
	std::vector<std::vector<std::uint8_t>> qualA(all_quals.begin(), all_quals.begin() + 3);
	std::vector<std::string> partB(all_seqs.begin() + 3, all_seqs.end());
	std::vector<std::vector<std::uint8_t>> qualB(all_quals.begin() + 3, all_quals.end());
	partA.insert(partA.end(), partB.begin(), partB.end());
	qualA.insert(qualA.end(), qualB.begin(), qualB.end());

	const auto combined = BuildConsensus(partA, qualA);
	REQUIRE(full.first == combined.first);
	REQUIRE(full.second == combined.second);
	REQUIRE(full.first == "TAAAAT");
}

TEST_CASE("BuildConsensus: n=1 bypass yields gap-stripped seq + unchanged qual", "[consensus][build]") {
	std::vector<std::string> aligned = {"A-CGT"};
	std::vector<std::vector<std::uint8_t>> quals = {{40, 40, 40, 40}};
	auto out = BuildConsensus(aligned, quals);
	REQUIRE(out.first == "ACGT");
	REQUIRE(out.second == std::vector<std::uint8_t>({40, 40, 40, 40}));
}

TEST_CASE("BuildConsensus: inconsistent widths throw", "[consensus][build]") {
	std::vector<std::string> aligned = {"ACGT", "ACGTA"};
	std::vector<std::vector<std::uint8_t>> quals = {{40, 40, 40, 40}, {40, 40, 40, 40, 40}};
	REQUIRE_THROWS_AS(BuildConsensus(aligned, quals), miint::InvalidInputException);
}
