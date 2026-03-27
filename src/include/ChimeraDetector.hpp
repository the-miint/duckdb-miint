#pragma once

#include "KmerIndex.hpp"
#include "WFA2Aligner.hpp"

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace miint {

// Diff classification for a single column in a 3-way star alignment.
// Based on Edgar et al. 2011, Bioinformatics 27:2194-2200.
enum class DiffType : char {
	IGNORE = ' ',  // All identical, gap/gap-adjacent, or ambiguous base
	MATCH_A = 'A', // Query matches parent A only (A != B, Q == A)
	MATCH_B = 'B', // Query matches parent B only (A != B, Q == B)
	NO_VOTE = 'N', // Both parents agree, query differs (A == B, Q != A)
	ABSTAIN = '?'  // Query matches neither parent (A != B, Q != A, Q != B)
};

// A 3-way star alignment with Q as the hub, A and B as arms.
struct StarAlignment {
	std::string query_row;
	std::string parent_a_row;
	std::string parent_b_row;
	std::vector<DiffType> diffs;
};

StarAlignment build_star_alignment(const std::string &query_aligned_a, const std::string &subject_aligned_a,
                                   const std::string &query_aligned_b, const std::string &subject_aligned_b);

void classify_diffs(StarAlignment &star);

// WFA2 penalty parameters equivalent to vsearch's NW scoring.
// See ChimeraDetector.cpp for detailed derivation.
static constexpr int UCHIME_WFA2_MISMATCH = 6;
static constexpr int UCHIME_WFA2_GAP_OPEN = 20;
static constexpr int UCHIME_WFA2_GAP_EXTEND = 4;
static constexpr int SMOOTHING_WINDOW = 32;

// UCHIME scoring parameters with defaults from Edgar et al. 2011.
struct UchimeParams {
	double minh = 0.28;
	double xn = 8.0;
	double dn = 1.4;
	double mindiv = 0.8;
	int mindiffs = 3;
	double abskew = 2.0; // only used in de novo mode
};

// Result of breakpoint sweep.
struct BreakpointResult {
	double best_h = 0.0;
	int best_pos = -1; // alignment column index of best breakpoint
	bool reversed = false;
	int left_yes = 0, left_no = 0, left_abstain = 0;
	int right_yes = 0, right_no = 0, right_abstain = 0;
};

// Full UCHIME result for a single query (mirrors vsearch --uchimeout 18 columns).
struct UchimeResult {
	double score = 0.0;
	std::string query_label;
	std::string parent_a_label;
	std::string parent_b_label;
	std::string closest_parent_label;
	double id_query_model = 0.0;
	double id_query_a = 0.0;
	double id_query_b = 0.0;
	double id_a_b = 0.0;
	double id_query_top = 0.0;
	int left_yes = 0, left_no = 0, left_abstain = 0;
	int right_yes = 0, right_no = 0, right_abstain = 0;
	double divergence = 0.0;
	std::string flag = "N"; // Y, N, or ?
};

// Selected parent pair with cached alignment results.
struct ParentPair {
	uint32_t parent_a_idx;
	uint32_t parent_b_idx;
	WFA2FullResult align_a; // query aligned to parent A
	WFA2FullResult align_b; // query aligned to parent B
};

// Compute per-position match profile from a pairwise alignment.
// Returns vector of 0/1: 1 if bases match and neither is a gap, else 0.
std::vector<int> compute_match_profile(const std::string &query_aligned, const std::string &subject_aligned);

// Compute smoothed identity using a sliding window.
// Returns vector of running sums of match_profile over window_size positions.
// Output[i] is defined for i >= window_size-1 (earlier positions are partial).
std::vector<int> compute_smoothed(const std::vector<int> &match_profile, int window_size = SMOOTHING_WINDOW);

// Select the two best parents from candidate indices.
// Aligns query to each candidate, selects Parent A (most "wins" by smoothed identity),
// then wipes A's winning positions and selects Parent B.
// Returns nullopt if fewer than 2 valid candidates exist.
std::optional<ParentPair> select_parents(const std::string &query, const std::vector<uint32_t> &candidate_indices,
                                         const std::vector<std::string> &ref_sequences, WFA2Aligner &aligner);

// Sweep all breakpoints left-to-right over a diff vector.
BreakpointResult sweep_breakpoints(const std::vector<DiffType> &diffs, const UchimeParams &params);

// Full UCHIME chimera detection pipeline.
class ChimeraDetector {
public:
	explicit ChimeraDetector(const UchimeParams &params = UchimeParams {});

	// Load reference sequences. Builds k-mer index.
	void set_reference(const std::vector<std::string> &labels, const std::vector<std::string> &sequences);

	// Add a single sequence to the reference (for de novo incremental mode).
	// abundance is optional (default 0, only used when detect_denovo filters by abskew).
	void add_to_reference(const std::string &label, const std::string &sequence, int64_t abundance = 0);

	// Detect chimera for a single query. Thread-safe (aligner is per-thread).
	UchimeResult detect(const std::string &query_label, const std::string &query_sequence, WFA2Aligner &aligner) const;

	// Detect chimera with abundance skew filtering (de novo mode).
	// Candidate parents must have abundance >= abskew * query_abundance.
	UchimeResult detect_denovo(const std::string &query_label, const std::string &query_sequence,
	                           int64_t query_abundance, WFA2Aligner &aligner) const;

	const std::vector<std::string> &ref_labels() const {
		return ref_labels_;
	}
	const std::vector<std::string> &ref_sequences() const {
		return ref_sequences_;
	}

private:
	UchimeParams params_;
	KmerIndex kmer_index_;
	std::vector<std::string> ref_labels_;
	std::vector<std::string> ref_sequences_;
	std::vector<int64_t> ref_abundances_; // For de novo abundance skew filtering

	// Shared pipeline: given pre-filtered candidates, run select_parents through classification.
	UchimeResult detect_impl(const std::string &query_label, const std::string &query_sequence,
	                         const std::vector<uint32_t> &candidates, WFA2Aligner &aligner) const;

	// Compute identity % between two aligned sequences over non-ignored columns.
	static double compute_identity(const std::string &aligned_a, const std::string &aligned_b);

	// Compute model identity: query matches model (parent A left of breakpoint, B right).
	static double compute_model_identity(const StarAlignment &star, int breakpoint, bool reversed);
};

} // namespace miint
