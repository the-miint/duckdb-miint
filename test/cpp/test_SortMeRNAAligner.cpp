#include <catch2/catch_test_macros.hpp>

#include "SortMeRNAAligner.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>

// Test-local constants and helpers. All golden values below are regenerated
// from the native sortmerna CLI via tools/sortmerna_oracle.sh — do NOT
// edit by hand; rerun the script if test inputs or the embedded sortmerna
// version change. Identity/coverage match native exactly at CLI's 1-decimal
// display precision; e_value is the only field we expect to differ (documented
// per-query Karlin-Altschul divergence in smr_api.h).
namespace {

// Path relative to the repo root (where tests are executed from). Same
// approach used by test_Minimap2Aligner.cpp.
const std::string kTestRef = "ext/sortmerna/data/test_ref.fasta";

// BLAST reports identity/coverage at 1-decimal precision. The library emits
// full precision from the same computation, so rounding the library value
// to one decimal must yield exactly the CLI value.
int64_t round_1dp(double v) {
	return static_cast<int64_t>(std::llround(v * 10.0));
}

// AB271211 16S rRNA read from ext/sortmerna/data/test_read.fasta (1487 nt).
constexpr const char *kAB271211 = "TCCAACGCGTTGGGAGCTCTCCCATATGGTCGACCTGCAGGCGGCCGCACTAGTGATTAGAGTTTGATCCTGGCTCAG"
                                  "GATGAACGCTGGCGGCGTGCCTAACACATGCAAGTCGAACGGGAATCTTCGGATTCTAGTGGCGGACGGGTGAGTAAC"
                                  "GCGTAAGAATCTAACTTCAGGACGGGGACAACAGTGGGAAACGACTGCTAATACCCGATGTGCCGCGAGGTGAAACCT"
                                  "AATTGGCCTGAAGAGGAGCTTGCGTCTGATTAGCTAGTTGGTGGGGTAAGAGCCTACCAAGGCGACGATCAGTAGCTG"
                                  "GTCTGAGAGGATGAGCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATTTTCC"
                                  "GCAATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGGGAGGAAGGTCTTTGGATTGTAAACCTCTTTTCTCAAGGA"
                                  "AGAAGTTCTGACGGTACTTGAGGAATCAGCCTCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGGGGAGGCAAGC"
                                  "GTTATCCGGAATTATTGGGCGTAAAGCGTCCGCAGGTGGTCAGCCAAGTCTGCCGTCAAATCAGGTTGCTTAACGACC"
                                  "TAAAGGCGGTGGAAACTGGCAGACTAGAGAGCAGTAGGGGTAGCAGGAATTCCCAGTGTAGCGGTGAAATGCGTAGAG"
                                  "ATTGGGAAGAACATCGGTGGCGAAAGCGTGCTACTGGGCTGTATCTGACACTCAGGGACGAAAGCTAGGGGAGCGAAA"
                                  "GGGATTAGATACCCCTGTAGTCCTAGCCGTAAACGATGGATACTAGGCGTGGCTTGTATCGACCCGAGCCGTGCCGAA"
                                  "GCTAACGCGTTAAGTATCCCGCCTGGGGAGTACGCACGCAAGTGTGAAACTCAAAGGAATTGACGGGGGCCCGCACAA"
                                  "GCGGTGGAGTATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCAAGACTTGACATGTCGCGAACCCTGGTGAAA"
                                  "GCTGGGGGTGCCTTCGGGAGCGCGAACACAGGTGGTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAG"
                                  "TCCCGCAACGAGCGCAACCCTCGTTCTTAGTTGCCAGCATTAAGTTGGGGACTCTAAGGAGACTGCCGGTGACAAACC"
                                  "GGAGGAAGGTGGGGATGACGTCAAGTCAGCATGCCCCTTACGTCTTGGGCGACACACGTACTACAATGGTCGGGACAA"
                                  "AGGGCAGCGAACTTGCGAGAGCCAGCGAATCCCAGCAAACCCGGCCTCAGTTCAGATTGCAGGCTGCAACTCGCCTGC"
                                  "ATGAAGGAGGAATCGCTAGTAATCGCCGGTCAGCATACGGCGGTGAATTCGTTCCCGGGCCTTGTACACACCGCCCGT"
                                  "CACACCATGGAAGCTGGTCACGCCCGAAGTCATTACCTCAACCGCAAGGAGGGGGATGCCTAAGGCAGGGCTAGTGAC"
                                  "TGGGG";

constexpr const char *kAB271211_CIGAR = "57S57M2I12M2D4M2I29M1D11M2I3M2D11M1I7M1D13M5D4M3D9M2D3M7D1260M";

// Reverse-complement of kAB271211, generated with `tr ACGT TGCA | rev`.
constexpr const char *kAB271211_RC = "CCCCAGTCACTAGCCCTGCCTTAGGCATCCCCCTCCTTGCGGTTGAGGTAATGACTTCGGGCGTGACCAGCTTCCATG"
                                     "GTGTGACGGGCGGTGTGTACAAGGCCCGGGAACGAATTCACCGCCGTATGCTGACCGGCGATTACTAGCGATTCCTCC"
                                     "TTCATGCAGGCGAGTTGCAGCCTGCAATCTGAACTGAGGCCGGGTTTGCTGGGATTCGCTGGCTCTCGCAAGTTCGCT"
                                     "GCCCTTTGTCCCGACCATTGTAGTACGTGTGTCGCCCAAGACGTAAGGGGCATGCTGACTTGACGTCATCCCCACCTT"
                                     "CCTCCGGTTTGTCACCGGCAGTCTCCTTAGAGTCCCCAACTTAATGCTGGCAACTAAGAACGAGGGTTGCGCTCGTTG"
                                     "CGGGACTTAACCCAACATCTCACGACACGAGCTGACGACAGCCATGCACCACCTGTGTTCGCGCTCCCGAAGGCACCC"
                                     "CCAGCTTTCACCAGGGTTCGCGACATGTCAAGTCTTGGTAAGGTTCTTCGCGTTGCATCGAATTAAACCACATACTCC"
                                     "ACCGCTTGTGCGGGCCCCCGTCAATTCCTTTGAGTTTCACACTTGCGTGCGTACTCCCCAGGCGGGATACTTAACGCG"
                                     "TTAGCTTCGGCACGGCTCGGGTCGATACAAGCCACGCCTAGTATCCATCGTTTACGGCTAGGACTACAGGGGTATCTA"
                                     "ATCCCTTTCGCTCCCCTAGCTTTCGTCCCTGAGTGTCAGATACAGCCCAGTAGCACGCTTTCGCCACCGATGTTCTTC"
                                     "CCAATCTCTACGCATTTCACCGCTACACTGGGAATTCCTGCTACCCCTACTGCTCTCTAGTCTGCCAGTTTCCACCGC"
                                     "CTTTAGGTCGTTAAGCAACCTGATTTGACGGCAGACTTGGCTGACCACCTGCGGACGCTTTACGCCCAATAATTCCGG"
                                     "ATAACGCTTGCCTCCCCCGTATTACCGCGGCTGCTGGCACGGAGTTAGCCGAGGCTGATTCCTCAAGTACCGTCAGAA"
                                     "CTTCTTCCTTGAGAAAAGAGGTTTACAATCCAAAGACCTTCCTCCCTCACGCGGCGTTGCTCCGTCAGGCTTTCGCCC"
                                     "ATTGCGGAAAATTCCCCACTGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAGTGTGGCTGCTCATCCTCT"
                                     "CAGACCAGCTACTGATCGTCGCCTTGGTAGGCTCTTACCCCACCAACTAGCTAATCAGACGCAAGCTCCTCTTCAGGC"
                                     "CAATTAGGTTTCACCTCGCGGCACATCGGGTATTAGCAGTCGTTTCCCACTGTTGTCCCCGTCCTGAAGTTAGATTCT"
                                     "TACGCGTTACTCACCCGTCCGCCACTAGAATCCGAAGATTCCCGTTCGACTTGCATGTGTTAGGCACGCCGCCAGCGT"
                                     "TCATCCTGAGCCAGGATCAAACTCTAATCACTAGTGCGGCCGCCTGCAGGTCGACCATATGGGAGAGCTCCCAACGCG"
                                     "TTGGA";

} // namespace

TEST_CASE("SortMeRNAAligner constructs with a valid reference", "[sortmerna][wrapper]") {
	REQUIRE(std::filesystem::exists(kTestRef));
	miint::SortMeRNAConfig cfg;
	cfg.num_threads = 1;
	miint::SortMeRNAAligner aligner(cfg, {kTestRef});
	SUCCEED();
}

TEST_CASE("SortMeRNAAligner throws on missing reference", "[sortmerna][wrapper]") {
	miint::SortMeRNAConfig cfg;
	cfg.num_threads = 1;
	REQUIRE_THROWS_AS(miint::SortMeRNAAligner(cfg, {"/nonexistent/reference.fasta"}), std::runtime_error);
}

TEST_CASE("SortMeRNAAligner throws on empty ref_paths", "[sortmerna][wrapper]") {
	miint::SortMeRNAConfig cfg;
	cfg.num_threads = 1;
	REQUIRE_THROWS_AS(miint::SortMeRNAAligner(cfg, {}), std::invalid_argument);
}

TEST_CASE("SortMeRNAAligner reports library version", "[sortmerna][wrapper]") {
	REQUIRE(miint::SortMeRNAAligner::version() == "4.4.0");
}

TEST_CASE("SortMeRNAAligner aligns single-end forward-strand read", "[sortmerna][wrapper]") {
	miint::SortMeRNAConfig cfg;
	cfg.num_threads = 1;
	miint::SortMeRNAAligner aligner(cfg, {kTestRef});

	miint::SortMeRNAQueryBatch queries;
	queries.read_ids.emplace_back("AB271211");
	queries.sequences.emplace_back(kAB271211);

	miint::SortMeRNAResultBatch out;
	aligner.align(queries, out);

	REQUIRE(out.size() == 1);
	REQUIRE(out.read_ids[0] == "AB271211");
	REQUIRE(out.aligned[0] == 1);
	REQUIRE(out.strands[0] == 1); // forward (SAM FLAG=0 → no 0x10 bit)
	REQUIRE(out.ref_names[0] == "Unc49508");
	REQUIRE(out.ref_starts[0] == 1);
	REQUIRE(out.ref_ends[0] == 1446);
	REQUIRE(out.cigars[0] == kAB271211_CIGAR);
	REQUIRE(out.scores[0] == 2430);       // SAM AS tag
	REQUIRE(out.edit_distances[0] == 94); // SAM NM tag
	REQUIRE(round_1dp(out.identities[0]) == 935);
	REQUIRE(round_1dp(out.coverages[0]) == 962);
	REQUIRE(out.segment_indices[0] == 0);
}

TEST_CASE("SortMeRNAAligner paired-end emits two rows per input with segment_idx 0/1", "[sortmerna][wrapper]") {
	miint::SortMeRNAConfig cfg;
	cfg.num_threads = 1;
	cfg.paired = true;
	miint::SortMeRNAAligner aligner(cfg, {kTestRef});

	miint::SortMeRNAQueryBatch queries;
	queries.read_ids.emplace_back("AB271211");
	queries.sequences.emplace_back(kAB271211);
	queries.sequences2.emplace_back(kAB271211_RC);

	miint::SortMeRNAResultBatch out;
	aligner.align(queries, out);

	REQUIRE(out.size() == 2);

	// Forward segment
	REQUIRE(out.read_ids[0] == "AB271211");
	REQUIRE(out.segment_indices[0] == 0);
	REQUIRE(out.aligned[0] == 1);
	REQUIRE(out.strands[0] == 1);
	REQUIRE(out.ref_names[0] == "Unc49508");
	REQUIRE(out.ref_starts[0] == 1);
	REQUIRE(out.ref_ends[0] == 1446);
	REQUIRE(out.cigars[0] == kAB271211_CIGAR);
	REQUIRE(out.scores[0] == 2430);
	REQUIRE(out.edit_distances[0] == 94);
	REQUIRE(round_1dp(out.identities[0]) == 935);
	REQUIRE(round_1dp(out.coverages[0]) == 962);

	// Reverse segment maps to same input row
	REQUIRE(out.read_ids[1] == "AB271211");
	REQUIRE(out.segment_indices[1] == 1);
	REQUIRE(out.aligned[1] == 1);
	REQUIRE(out.strands[1] == 0);
	REQUIRE(out.ref_names[1] == "Unc49508");
	REQUIRE(out.ref_starts[1] == 1);
	REQUIRE(out.ref_ends[1] == 1446);
	REQUIRE(out.cigars[1] == kAB271211_CIGAR);
	REQUIRE(out.scores[1] == 2430);
	REQUIRE(out.edit_distances[1] == 94);
	REQUIRE(round_1dp(out.identities[1]) == 935);
	REQUIRE(round_1dp(out.coverages[1]) == 962);
}

TEST_CASE("SortMeRNAAligner rejects malformed query batch shapes", "[sortmerna][wrapper]") {
	miint::SortMeRNAConfig cfg;
	cfg.num_threads = 1;
	miint::SortMeRNAAligner aligner(cfg, {kTestRef});
	miint::SortMeRNAResultBatch out;

	// read_ids / sequences size mismatch.
	{
		miint::SortMeRNAQueryBatch q;
		q.read_ids = {"a", "b"};
		q.sequences = {kAB271211};
		REQUIRE_THROWS_AS(aligner.align(q, out), std::invalid_argument);
	}
	// sequences2 present but different size from sequences.
	{
		miint::SortMeRNAQueryBatch q;
		q.read_ids = {"a"};
		q.sequences = {kAB271211};
		q.sequences2 = {kAB271211_RC, kAB271211_RC};
		REQUIRE_THROWS_AS(aligner.align(q, out), std::invalid_argument);
	}
}

TEST_CASE("SortMeRNAAligner rejects single-end batch on paired aligner", "[sortmerna][wrapper]") {
	miint::SortMeRNAConfig cfg;
	cfg.num_threads = 1;
	cfg.paired = true;
	miint::SortMeRNAAligner aligner(cfg, {kTestRef});

	miint::SortMeRNAQueryBatch queries;
	queries.read_ids.emplace_back("AB271211");
	queries.sequences.emplace_back(kAB271211);
	// No sequences2 — queries.is_paired() == false, mismatches aligner config.
	miint::SortMeRNAResultBatch out;
	REQUIRE_THROWS_AS(aligner.align(queries, out), std::invalid_argument);
}

TEST_CASE("SortMeRNAAligner rejects paired batch on single-end aligner", "[sortmerna][wrapper]") {
	miint::SortMeRNAConfig cfg;
	cfg.num_threads = 1;
	cfg.paired = false;
	miint::SortMeRNAAligner aligner(cfg, {kTestRef});

	miint::SortMeRNAQueryBatch queries;
	queries.read_ids.emplace_back("AB271211");
	queries.sequences.emplace_back(kAB271211);
	queries.sequences2.emplace_back(kAB271211_RC);
	miint::SortMeRNAResultBatch out;
	REQUIRE_THROWS_AS(aligner.align(queries, out), std::invalid_argument);
}

TEST_CASE("SortMeRNAAligner aligns reverse-complement read as strand=0", "[sortmerna][wrapper]") {
	miint::SortMeRNAConfig cfg;
	cfg.num_threads = 1;
	miint::SortMeRNAAligner aligner(cfg, {kTestRef});

	miint::SortMeRNAQueryBatch queries;
	queries.read_ids.emplace_back("AB271211_rc");
	queries.sequences.emplace_back(kAB271211_RC);

	miint::SortMeRNAResultBatch out;
	aligner.align(queries, out);

	REQUIRE(out.size() == 1);
	REQUIRE(out.aligned[0] == 1);
	REQUIRE(out.strands[0] == 0); // reverse-complement (SAM FLAG=16 → 0x10 bit set)
	REQUIRE(out.ref_names[0] == "Unc49508");
	REQUIRE(out.ref_starts[0] == 1);
	REQUIRE(out.ref_ends[0] == 1446);
	REQUIRE(out.cigars[0] == kAB271211_CIGAR);
	REQUIRE(out.scores[0] == 2430);
	REQUIRE(out.edit_distances[0] == 94);
	REQUIRE(round_1dp(out.identities[0]) == 935);
	REQUIRE(round_1dp(out.coverages[0]) == 962);
}

TEST_CASE("SortMeRNAAligner empty query batch is rejected", "[sortmerna][wrapper]") {
	miint::SortMeRNAConfig cfg;
	cfg.num_threads = 1;
	miint::SortMeRNAAligner aligner(cfg, {kTestRef});

	miint::SortMeRNAQueryBatch queries;
	miint::SortMeRNAResultBatch out;
	REQUIRE_THROWS_AS(aligner.align(queries, out), std::invalid_argument);
}

TEST_CASE("SortMeRNAAligner batch-split invariance: 1x6 == 3x2", "[sortmerna][wrapper]") {
	miint::SortMeRNAConfig cfg;
	cfg.num_threads = 1;
	miint::SortMeRNAAligner aligner(cfg, {kTestRef});

	// 6 copies of AB271211 — all strong hits, no threshold ambiguity.
	std::vector<std::string> ids = {"r0", "r1", "r2", "r3", "r4", "r5"};
	std::vector<std::string> seqs(6, kAB271211);

	// Monolithic: one batch of 6.
	miint::SortMeRNAResultBatch mono;
	{
		miint::SortMeRNAQueryBatch q;
		q.read_ids = ids;
		q.sequences = seqs;
		aligner.align(q, mono);
	}

	// Split: three batches of 2.
	miint::SortMeRNAResultBatch split;
	for (int i = 0; i < 3; ++i) {
		miint::SortMeRNAQueryBatch q;
		q.read_ids = {ids[2 * i], ids[2 * i + 1]};
		q.sequences = {seqs[2 * i], seqs[2 * i + 1]};
		aligner.align(q, split);
	}

	REQUIRE(mono.size() == 6);
	REQUIRE(split.size() == 6);

	// Every field byte-identical across monolithic vs split. The streaming
	// API is documented to guarantee batch-split invariance for identity,
	// coverage, score, edit_distance, CIGAR, ref coords — and by extension
	// for e_value (per-query Karlin-Altschul, n_read-local).
	for (size_t i = 0; i < 6; ++i) {
		REQUIRE(mono.read_ids[i] == split.read_ids[i]);
		REQUIRE(mono.aligned[i] == split.aligned[i]);
		REQUIRE(mono.strands[i] == split.strands[i]);
		REQUIRE(mono.ref_names[i] == split.ref_names[i]);
		REQUIRE(mono.ref_starts[i] == split.ref_starts[i]);
		REQUIRE(mono.ref_ends[i] == split.ref_ends[i]);
		REQUIRE(mono.cigars[i] == split.cigars[i]);
		REQUIRE(mono.scores[i] == split.scores[i]);
		REQUIRE(mono.edit_distances[i] == split.edit_distances[i]);
		REQUIRE(mono.identities[i] == split.identities[i]);
		REQUIRE(mono.coverages[i] == split.coverages[i]);
		REQUIRE(mono.e_values[i] == split.e_values[i]);
	}
}

TEST_CASE("SortMeRNAAligner produces identical results at num_threads 1 and 8", "[sortmerna][wrapper]") {
	// 12 AB271211 copies — enough work to exercise sortmerna's internal pool.
	std::vector<std::string> ids;
	std::vector<std::string> seqs;
	for (int i = 0; i < 12; ++i) {
		ids.push_back("r" + std::to_string(i));
		seqs.push_back(kAB271211);
	}

	auto run = [&](int32_t nthreads) {
		miint::SortMeRNAConfig cfg;
		cfg.num_threads = nthreads;
		miint::SortMeRNAAligner aligner(cfg, {kTestRef});
		miint::SortMeRNAQueryBatch q;
		q.read_ids = ids;
		q.sequences = seqs;
		miint::SortMeRNAResultBatch out;
		aligner.align(q, out);
		return out;
	};

	auto t1 = run(1);
	auto t8 = run(8);

	REQUIRE(t1.size() == 12);
	REQUIRE(t8.size() == 12);
	for (size_t i = 0; i < 12; ++i) {
		REQUIRE(t1.read_ids[i] == t8.read_ids[i]);
		REQUIRE(t1.aligned[i] == t8.aligned[i]);
		REQUIRE(t1.strands[i] == t8.strands[i]);
		REQUIRE(t1.ref_names[i] == t8.ref_names[i]);
		REQUIRE(t1.ref_starts[i] == t8.ref_starts[i]);
		REQUIRE(t1.ref_ends[i] == t8.ref_ends[i]);
		REQUIRE(t1.cigars[i] == t8.cigars[i]);
		REQUIRE(t1.scores[i] == t8.scores[i]);
		REQUIRE(t1.edit_distances[i] == t8.edit_distances[i]);
		REQUIRE(t1.identities[i] == t8.identities[i]);
		REQUIRE(t1.coverages[i] == t8.coverages[i]);
		// Per-query Karlin-Altschul form is deterministic at read level;
		// thread count must not perturb it.
		REQUIRE(t1.e_values[i] == t8.e_values[i]);
	}
}

TEST_CASE("SortMeRNAAligner appends to existing output batch", "[sortmerna][wrapper]") {
	miint::SortMeRNAConfig cfg;
	cfg.num_threads = 1;
	miint::SortMeRNAAligner aligner(cfg, {kTestRef});

	miint::SortMeRNAResultBatch out;
	for (int i = 0; i < 3; ++i) {
		miint::SortMeRNAQueryBatch q;
		q.read_ids.emplace_back("AB271211_" + std::to_string(i));
		q.sequences.emplace_back(kAB271211);
		aligner.align(q, out);
	}
	REQUIRE(out.size() == 3);
	REQUIRE(out.read_ids[0] == "AB271211_0");
	REQUIRE(out.read_ids[1] == "AB271211_1");
	REQUIRE(out.read_ids[2] == "AB271211_2");
}
