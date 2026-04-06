#include "../../src/include/deblur.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Helper: parse FASTA with ;size=N; annotation into DeblurSequence vector
// ---------------------------------------------------------------------------
static std::vector<miint::DeblurSequence> load_test_fasta(const std::string &path) {
	std::ifstream in(path);
	REQUIRE(in.good());

	std::vector<miint::DeblurSequence> seqs;
	std::string label;
	std::string seq_buf;

	auto flush = [&]() {
		if (label.empty()) {
			return;
		}
		// Parse size=N from label
		std::regex size_re(R"(size=(\d+))");
		std::smatch m;
		REQUIRE(std::regex_search(label, m, size_re));
		double freq = std::stod(m[1].str());

		// Strip label to just the ID part (before first ;)
		std::string id = label.substr(0, label.find(';'));

		seqs.push_back({id, seq_buf, freq});
		label.clear();
		seq_buf.clear();
	};

	std::string line;
	while (std::getline(in, line)) {
		if (line.empty()) {
			continue;
		}
		if (line[0] == '>') {
			flush();
			label = line.substr(1); // strip '>'
		} else {
			seq_buf += line;
		}
	}
	flush();
	return seqs;
}

// ---------------------------------------------------------------------------
// Cycle 1: compute_distance_hamming
// ---------------------------------------------------------------------------

TEST_CASE("compute_distance_hamming - identical sequences", "[deblur]") {
	std::string a = "ACGTACGT";
	auto result = miint::compute_distance_hamming(a.c_str(), a.c_str(), a.size());
	REQUIRE(result.hamming_distance == 0);
	REQUIRE(result.num_indels == 0);
	REQUIRE(result.num_substitutions == 0);
}

TEST_CASE("compute_distance_hamming - one substitution", "[deblur]") {
	std::string a = "ACGTACGT";
	std::string b = "ACGAACGT";
	auto result = miint::compute_distance_hamming(a.c_str(), b.c_str(), a.size());
	REQUIRE(result.hamming_distance == 1);
	REQUIRE(result.num_indels == 0);
	REQUIRE(result.num_substitutions == 1);
}

TEST_CASE("compute_distance_hamming - multiple substitutions", "[deblur]") {
	std::string a = "ACGTACGT";
	std::string b = "TCGAACGA";
	// diffs at pos 0 (A/T), 3 (T/A), 7 (T/A) = 3
	auto result = miint::compute_distance_hamming(a.c_str(), b.c_str(), a.size());
	REQUIRE(result.hamming_distance == 3);
	REQUIRE(result.num_indels == 0);
	REQUIRE(result.num_substitutions == 3);
}

TEST_CASE("compute_distance_hamming - internal gap (indel)", "[deblur]") {
	// seq_a has a gap at position 3, seq_b has a base there
	std::string a = "ACG-ACGT";
	std::string b = "ACGTACGT";
	// Full hamming: pos 3 differs ('-' vs 'T') = 1
	auto result = miint::compute_distance_hamming(a.c_str(), b.c_str(), a.size());
	REQUIRE(result.hamming_distance == 1);
	REQUIRE(result.num_indels == 1);
	REQUIRE(result.num_substitutions == 0);
}

TEST_CASE("compute_distance_hamming - trailing gap trimming", "[deblur]") {
	// seq_a is shorter (unaligned 4bp), seq_b is full 8bp
	// Trailing gaps in seq_a should be trimmed for indel/sub counting
	std::string a = "ACGT----";
	std::string b = "ACGTACGT";
	// Full hamming: positions 4-7 differ = 4
	// Trimmed region (min(4,8)=4): "ACGT" vs "ACGT" = 0 diffs, 0 indels
	// Since num_indels==0, Python uses full hamming: num_substitutions = 4 - 0 = 4
	// (In practice this case is prevented by deblur's equal-unaligned-length validation)
	auto result = miint::compute_distance_hamming(a.c_str(), b.c_str(), a.size());
	REQUIRE(result.hamming_distance == 4);
	REQUIRE(result.num_indels == 0);
	REQUIRE(result.num_substitutions == 4);
}

TEST_CASE("compute_distance_hamming - gap + substitution", "[deblur]") {
	std::string a = "AC-TACGT";
	std::string b = "ACGTACGA";
	// Full hamming: pos 2 ('-' vs 'G'), pos 7 ('T' vs 'A') = 2
	// Both have no trailing gaps, trim_len = 8
	// Indels: pos 2 (one is gap) = 1
	// num_subs = 2 - 1 = 1
	auto result = miint::compute_distance_hamming(a.c_str(), b.c_str(), a.size());
	REQUIRE(result.hamming_distance == 2);
	REQUIRE(result.num_indels == 1);
	REQUIRE(result.num_substitutions == 1);
}

TEST_CASE("compute_distance_hamming - both have trailing gaps", "[deblur]") {
	std::string a = "ACGT----";
	std::string b = "ACGA----";
	// Full hamming: pos 3 (T vs A) = 1
	// Trimmed to min(4,4) = 4: "ACGT" vs "ACGA" = 1 diff, no gaps
	auto result = miint::compute_distance_hamming(a.c_str(), b.c_str(), a.size());
	REQUIRE(result.hamming_distance == 1);
	REQUIRE(result.num_indels == 0);
	REQUIRE(result.num_substitutions == 1);
}

TEST_CASE("compute_distance_hamming - TEST_SEQS_1 E.Coli vs Error", "[deblur]") {
	// These are uppercased 149bp sequences from TEST_SEQS_1, differing at pos 0 (T vs A)
	std::string ecoli = "TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGT"
	                    "TAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTT"
	                    "GAGTCTCGTAGAGGGGGGCAGAATTCCAG";
	std::string error_seq = "AACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGT"
	                        "TAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTT"
	                        "GAGTCTCGTAGAGGGGGGCAGAATTCCAG";
	REQUIRE(ecoli.size() == 149);
	REQUIRE(error_seq.size() == 149);
	auto result = miint::compute_distance_hamming(ecoli.c_str(), error_seq.c_str(), ecoli.size());
	REQUIRE(result.hamming_distance == 1);
	REQUIRE(result.num_indels == 0);
	REQUIRE(result.num_substitutions == 1);
}

// ---------------------------------------------------------------------------
// Cycle 2: get_default_error_profile + deblur core
// ---------------------------------------------------------------------------

TEST_CASE("get_default_error_profile", "[deblur]") {
	auto profile = miint::get_default_error_profile();
	REQUIRE(profile.size() == 12);
	REQUIRE(profile[0] == 1.0);
	REQUIRE(profile[1] == 0.06);
	REQUIRE(profile[2] == 0.02);
	REQUIRE(profile[11] == 0.0005);
}

TEST_CASE("deblur - empty input", "[deblur]") {
	auto results = miint::deblur({});
	REQUIRE(results.empty());
}

TEST_CASE("deblur - single sequence passes through", "[deblur]") {
	std::vector<miint::DeblurSequence> seqs = {
	    {"only", "ACGTACGT", 100.0},
	};
	auto results = miint::deblur(std::move(seqs));
	REQUIRE(results.size() == 1);
	REQUIRE(results[0].label == "only");
	REQUIRE(results[0].abundance == 100);
	REQUIRE(results[0].sequence == "ACGTACGT");
}

TEST_CASE("deblur - TEST_SEQS_1 toy example", "[deblur]") {
	// E.Coli (size=1000) + Error (size=3), 149bp no gaps
	// Python deblur returns only E.Coli with abundance 1000
	std::vector<miint::DeblurSequence> seqs = {
	    {"E.Coli",
	     "TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGT"
	     "TAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTT"
	     "GAGTCTCGTAGAGGGGGGCAGAATTCCAG",
	     1000.0},
	    {"Error",
	     "AACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGT"
	     "TAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTT"
	     "GAGTCTCGTAGAGGGGGGCAGAATTCCAG",
	     3.0},
	};
	auto results = miint::deblur(std::move(seqs));
	REQUIRE(results.size() == 1);
	REQUIRE(results[0].label == "E.Coli");
	REQUIRE(results[0].abundance == 1000);
	// Output should have no gaps
	REQUIRE(results[0].sequence.find('-') == std::string::npos);
}

TEST_CASE("deblur - uppercases input (TEST_SEQS_1 lowercase)", "[deblur]") {
	// Same as TEST_SEQS_1 but lowercase, matching the original FASTA
	std::vector<miint::DeblurSequence> seqs = {
	    {"E.Coli",
	     "tacggagggtgcaagcgttaatcggaattactgggcgtaaagcgcacgcaggcggtttgt"
	     "taagtcagatgtgaaatccccgggctcaacctgggaactgcatctgatactggcaagctt"
	     "gagtctcgtagaggggggcagaattccag",
	     1000.0},
	    {"Error",
	     "aacggagggtgcaagcgttaatcggaattactgggcgtaaagcgcacgcaggcggtttgt"
	     "taagtcagatgtgaaatccccgggctcaacctgggaactgcatctgatactggcaagctt"
	     "gagtctcgtagaggggggcagaattccag",
	     3.0},
	};
	auto results = miint::deblur(std::move(seqs));
	REQUIRE(results.size() == 1);
	REQUIRE(results[0].label == "E.Coli");
	REQUIRE(results[0].abundance == 1000);
	// Output should be uppercased
	for (char c : results[0].sequence) {
		REQUIRE(c == std::toupper(c));
	}
}

// ---------------------------------------------------------------------------
// Cycle 3: TEST_SEQS_2, indels, custom params
// ---------------------------------------------------------------------------

TEST_CASE("deblur - TEST_SEQS_2 full dataset from file", "[deblur]") {
	auto seqs = load_test_fasta("data/deblur/test_seqs_2.fa");
	REQUIRE(seqs.size() == 61);

	auto results = miint::deblur(std::move(seqs));
	REQUIRE(results.size() == 1);
	REQUIRE(results[0].label == "E.Coli-999");
	REQUIRE(results[0].abundance == 720);
}

TEST_CASE("deblur - indel handling", "[deblur]") {
	// Matches test_deblur_indel from Python test suite.
	// Take TEST_SEQS_2, insert gap at position 10 of every sequence,
	// then add two indel sequences:
	//   indel1 (size=30): 'A' inserted at pos 10, trailing gap — gets removed
	//   indel2 (size=31): same modification — survives (above correction threshold)
	auto seqs = load_test_fasta("data/deblur/test_seqs_2.fa");
	REQUIRE(seqs.size() == 61);

	// Insert gap at position 10 of all sequences (matching Python test)
	for (auto &s : seqs) {
		s.aligned_sequence = s.aligned_sequence.substr(0, 10) + "-" + s.aligned_sequence.substr(10);
	}

	// Create indel1: 'A' at pos 10, trailing gap compensates length
	std::string base_seq = seqs[0].aligned_sequence;
	std::string indel_seq = base_seq.substr(0, 10) + "A" + base_seq.substr(11, base_seq.size() - 12) + "-";
	seqs.push_back({"indel1-read", indel_seq, 30.0});

	// Create indel2: same modification but higher frequency
	seqs.push_back({"indel2-read", indel_seq, 31.0});

	auto results = miint::deblur(std::move(seqs));
	REQUIRE(results.size() == 2);
	REQUIRE(results[0].label == "E.Coli-999");
	REQUIRE(results[1].label == "indel2-read");
}

TEST_CASE("deblur - custom error profile", "[deblur]") {
	auto seqs = load_test_fasta("data/deblur/test_seqs_2.fa");

	miint::DeblurParams params;
	params.error_dist = {1,         0.05,      0.000005,  0.000005,  0.000005,  0.000005,  0.0000025, 0.0000025,
	                     0.0000025, 0.0000025, 0.0000025, 0.0000005, 0.0000005, 0.0000005, 0.0000005};

	auto results = miint::deblur(std::move(seqs), params);
	REQUIRE(results.size() == 1);
	REQUIRE(results[0].label == "E.Coli-999");
}

TEST_CASE("deblur - default profile as explicit param matches default", "[deblur]") {
	auto seqs = load_test_fasta("data/deblur/test_seqs_2.fa");

	miint::DeblurParams params;
	params.error_dist = miint::get_default_error_profile();

	auto results = miint::deblur(std::move(seqs), params);
	REQUIRE(results.size() == 1);
	REQUIRE(results[0].label == "E.Coli-999");
	REQUIRE(results[0].abundance == 720);
}

// ---------------------------------------------------------------------------
// Additional tests from code review (#7)
// ---------------------------------------------------------------------------

TEST_CASE("deblur - mean_error affects normalization", "[deblur]") {
	// With very low mean_error, mod_factor is close to 1, so effective error
	// rates are lower. With mean_error=0.0001 the correction for 1-distance
	// at freq=1000 is num_err[1] = 0.06/mod_factor * 1000. For 149bp:
	//   mod_factor = (1-0.0001)^149 ≈ 0.985, num_err[1] ≈ 60.9 — still > 3
	// Error (freq=3) still gets removed.
	std::vector<miint::DeblurSequence> seqs = {
	    {"main",
	     "TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGT"
	     "TAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTT"
	     "GAGTCTCGTAGAGGGGGGCAGAATTCCAG",
	     1000.0},
	    {"error",
	     "AACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGT"
	     "TAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTT"
	     "GAGTCTCGTAGAGGGGGGCAGAATTCCAG",
	     3.0},
	};
	miint::DeblurParams params;
	params.mean_error = 0.0001;
	auto results = miint::deblur(std::move(seqs), params);
	REQUIRE(results.size() == 1);
	REQUIRE(results[0].label == "main");
}

TEST_CASE("deblur - indel_prob=0 suppresses all indel corrections", "[deblur]") {
	// With indel_prob=0, indel-based corrections are multiplied by 0,
	// so indel sequences are never removed by indel correction.
	auto seqs = load_test_fasta("data/deblur/test_seqs_2.fa");

	// Insert gap at pos 10
	for (auto &s : seqs) {
		s.aligned_sequence = s.aligned_sequence.substr(0, 10) + "-" + s.aligned_sequence.substr(10);
	}
	std::string base_seq = seqs[0].aligned_sequence;
	std::string indel_seq = base_seq.substr(0, 10) + "A" + base_seq.substr(11, base_seq.size() - 12) + "-";
	// indel1 at freq=30 would normally be removed with default indel_prob=0.01
	seqs.push_back({"indel1-read", indel_seq, 30.0});

	miint::DeblurParams params;
	params.indel_prob = 0.0;
	auto results = miint::deblur(std::move(seqs), params);
	// indel1 should survive since indel_prob=0 zeroes out the correction
	bool found_indel1 = false;
	for (auto &r : results) {
		if (r.label == "indel1-read") {
			found_indel1 = true;
		}
	}
	REQUIRE(found_indel1);
}

TEST_CASE("compute_distance_hamming - internal gaps with equal unaligned lengths", "[deblur]") {
	// Both sequences have 6 non-gap characters, but gaps at different positions.
	// This tests the case that deblur() will actually encounter.
	std::string a = "AC-GTACGT";
	std::string b = "ACG-TACGT";
	// Full hamming: pos 2 ('-' vs 'G'), pos 3 ('G' vs '-') = 2
	// Both have no trailing gaps, trim_len = 9
	// Indels: pos 2 (one gap), pos 3 (one gap) = 2
	// Trimmed hamming = 2 (indels present, so use trimmed)
	// num_subs = 2 - 2 = 0
	auto result = miint::compute_distance_hamming(a.c_str(), b.c_str(), a.size());
	REQUIRE(result.hamming_distance == 2);
	REQUIRE(result.num_indels == 2);
	REQUIRE(result.num_substitutions == 0);
}
