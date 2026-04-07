#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

struct DeblurParams {
	double mean_error = 0.005;
	std::vector<double> error_dist; // empty = use default
	double indel_prob = 0.01;
	int indel_max = 3;
};

struct DeblurSequence {
	std::string label;
	std::string aligned_sequence; // with gaps
	double frequency;
};

struct DeblurResult {
	std::string label;
	std::string sequence; // gaps stripped
	int64_t abundance;    // round(corrected_frequency) via banker's rounding
};

// Hamming distance result for aligned DNA sequences with gap awareness.
//
// Computes total mismatches and separately tracks gap-vs-base mismatches
// (indels) vs base-vs-base mismatches (substitutions). Operates on raw
// ASCII chars; inner loops are branchless to enable autovectorization.
// Gap character is '-'.
//
// hamming_distance: full aligned length, used for early-exit (> max_h_dist) check
// num_indels: in the trimmed region (trailing gaps stripped)
// num_substitutions: from trimmed hamming when indels present, else from full hamming
struct HammingResult {
	int hamming_distance;
	int num_indels;
	int num_substitutions;
};

HammingResult compute_distance_hamming(const char *seq_a, const char *seq_b, int len);

// Overload accepting pre-computed trailing-gap-stripped prefix lengths,
// avoiding redundant scanning in the inner loop of deblur().
HammingResult compute_distance_hamming(const char *seq_a, const char *seq_b, int len, int prefix_len_a,
                                       int prefix_len_b);

std::vector<double> get_default_error_profile();

std::vector<DeblurResult> deblur(std::vector<DeblurSequence> sequences, const DeblurParams &params = DeblurParams {});

} // namespace miint
