#include "deblur.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace miint {

// Hamming distance for aligned DNA sequences. Separately tracks gap-vs-base
// mismatches (indels) vs base-vs-base mismatches (substitutions). Operates on
// raw ASCII chars; inner loops are branchless for autovectorization.
// Gap character is '-'.
HammingResult compute_distance_hamming(const char *seq_a, const char *seq_b, int len) {
	// Full-length hamming distance (for > max_h_dist early exit check)
	int full_hamming = 0;
	for (int i = 0; i < len; i++) {
		full_hamming += (seq_a[i] != seq_b[i]);
	}

	// Compute unaligned lengths (strip trailing gaps)
	int len_a = len;
	while (len_a > 0 && seq_a[len_a - 1] == '-') {
		len_a--;
	}
	int len_b = len;
	while (len_b > 0 && seq_b[len_b - 1] == '-') {
		len_b--;
	}
	int trim_len = std::min(len_a, len_b);

	// Trimmed pass: count hamming and indels in the trimmed region
	int trimmed_hamming = 0;
	int num_indels = 0;
	for (int i = 0; i < trim_len; i++) {
		bool diff = (seq_a[i] != seq_b[i]);
		bool a_gap = (seq_a[i] == '-');
		bool b_gap = (seq_b[i] == '-');
		trimmed_hamming += diff;
		num_indels += (diff & (a_gap | b_gap));
	}

	// When indels are present, use trimmed hamming for substitution count
	// (matching Python deblur behavior: h_dist recomputed on trimmed region)
	int effective_hamming = (num_indels > 0) ? trimmed_hamming : full_hamming;
	int num_substitutions = effective_hamming - num_indels;

	return {full_hamming, num_indels, num_substitutions};
}

// Overload accepting pre-computed trailing-gap-stripped lengths, avoiding
// redundant scanning in the inner loop of deblur().
HammingResult compute_distance_hamming(const char *seq_a, const char *seq_b, int len, int prefix_len_a,
                                       int prefix_len_b) {
	int full_hamming = 0;
	for (int i = 0; i < len; i++) {
		full_hamming += (seq_a[i] != seq_b[i]);
	}

	int trim_len = std::min(prefix_len_a, prefix_len_b);
	int trimmed_hamming = 0;
	int num_indels = 0;
	for (int i = 0; i < trim_len; i++) {
		bool diff = (seq_a[i] != seq_b[i]);
		bool a_gap = (seq_a[i] == '-');
		bool b_gap = (seq_b[i] == '-');
		trimmed_hamming += diff;
		num_indels += (diff & (a_gap | b_gap));
	}

	int effective_hamming = (num_indels > 0) ? trimmed_hamming : full_hamming;
	int num_substitutions = effective_hamming - num_indels;

	return {full_hamming, num_indels, num_substitutions};
}

std::vector<double> get_default_error_profile() {
	return {1, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, 0.005, 0.001, 0.001, 0.001, 0.0005};
}

std::vector<DeblurResult> deblur(std::vector<DeblurSequence> sequences, const DeblurParams &params) {
	if (sequences.empty()) {
		return {};
	}

	const size_t n = sequences.size();

	// Uppercase all sequences (Python's Sequence.__init__ calls .upper())
	for (auto &s : sequences) {
		for (auto &c : s.aligned_sequence) {
			c = std::toupper(static_cast<unsigned char>(c));
		}
	}

	// Sort by frequency descending (stable sort for deterministic tie-breaking)
	std::stable_sort(sequences.begin(), sequences.end(),
	                 [](const DeblurSequence &a, const DeblurSequence &b) { return a.frequency > b.frequency; });

	// Validate: all aligned lengths must be equal
	const int aligned_len = static_cast<int>(sequences[0].aligned_sequence.size());
	for (size_t i = 1; i < n; i++) {
		if (static_cast<int>(sequences[i].aligned_sequence.size()) != aligned_len) {
			throw std::invalid_argument("Not all sequences have the same aligned length");
		}
	}

	// Two distinct "length" concepts for aligned sequences:
	//
	// prefix_lens[i]: index of last non-gap character + 1. Used for trimming
	//   trailing gaps in pairwise distance computation. Example: "ACG---" has
	//   prefix_len=3. This is what Python's rstrip('-') produces.
	//
	// unaligned_len: total count of non-gap characters. Used for error profile
	//   normalization (mod_factor). Example: "A-G---" has unaligned_len=2.
	//
	// For well-formed MSA output (same input lengths, only internal gap insertion),
	// these differ only when internal gaps are present.
	std::vector<int> prefix_lens(n);
	int unaligned_len = 0;
	for (size_t i = 0; i < n; i++) {
		const auto &seq = sequences[i].aligned_sequence;
		int last = aligned_len;
		while (last > 0 && seq[last - 1] == '-') {
			last--;
		}
		prefix_lens[i] = last;
		int non_gap = aligned_len - static_cast<int>(std::count(seq.begin(), seq.end(), '-'));
		if (i == 0) {
			unaligned_len = non_gap;
		} else if (non_gap != unaligned_len) {
			throw std::invalid_argument("Not all sequences have the same unaligned length");
		}
	}

	// Resolve and normalize error profile
	std::vector<double> error_dist = params.error_dist.empty() ? get_default_error_profile() : params.error_dist;
	double mod_factor = std::pow(1.0 - params.mean_error, unaligned_len);
	for (auto &e : error_dist) {
		e /= mod_factor;
	}

	const int max_h_dist = static_cast<int>(error_dist.size()) - 1;

	// Pre-allocate num_err buffer (reused per seq_i iteration)
	std::vector<double> num_err(error_dist.size());

	// Main deblur loop
	for (size_t i = 0; i < n; i++) {
		if (sequences[i].frequency <= 0) {
			continue;
		}

		for (size_t k = 0; k < error_dist.size(); k++) {
			num_err[k] = error_dist[k] * sequences[i].frequency;
		}

		if (num_err[1] < 0.1) {
			continue;
		}

		const char *seq_i_data = sequences[i].aligned_sequence.c_str();

		for (size_t j = 0; j < n; j++) {
			if (i == j) {
				continue;
			}

			auto dist = compute_distance_hamming(seq_i_data, sequences[j].aligned_sequence.c_str(), aligned_len,
			                                     prefix_lens[i], prefix_lens[j]);

			if (dist.hamming_distance > max_h_dist) {
				continue;
			}

			// Index is safe: num_substitutions <= hamming_distance <= max_h_dist.
			// When num_indels > 0 and num_substitutions == 0, this indexes
			// num_err[0] (= error_dist[0] * freq) — matching Python's behavior
			// where indel-only neighbors get a large base correction, then
			// reduced by indel_prob.
			double correction = num_err[dist.num_substitutions];

			if (dist.num_indels > params.indel_max) {
				correction = 0;
			} else if (dist.num_indels > 0) {
				correction *= params.indel_prob;
			}

			sequences[j].frequency -= correction;
		}
	}

	// Filter: keep sequences where round(frequency) > 0
	// std::nearbyint uses IEEE 754 round-to-nearest-even (banker's rounding),
	// matching Python 3's built-in round()
	std::vector<DeblurResult> results;
	results.reserve(n); // upper bound
	for (auto &s : sequences) {
		auto rounded = static_cast<int64_t>(std::nearbyint(s.frequency));
		if (rounded > 0) {
			std::string ungapped;
			ungapped.reserve(unaligned_len);
			for (char c : s.aligned_sequence) {
				if (c != '-') {
					ungapped += c;
				}
			}
			results.push_back({s.label, std::move(ungapped), rounded});
		}
	}

	return results;
}

} // namespace miint
