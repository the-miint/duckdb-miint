#pragma once

#include "KmerIndex.hpp"
#include "WFA2Aligner.hpp"
#include "uchime_common.hpp"

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
// WFA2 penalty parameters for UCHIME chimera detection alignment.
//
// vsearch uses: match=+2, mismatch=-4, gap_open=-20, gap_extend=-2
// (interior), gap_extend=-1 (terminal).
//
// For the mismatch penalty: relative cost of mismatch vs match =
// match + |mismatch| = 2 + 4 = 6.
//
// For gap penalties: vsearch's gap_extend=-2 means each gap extension base
// costs 2 in absolute terms. The match-reward conversion would give 2+2=4,
// but vsearch also uses cheaper terminal gaps (extend=-1). Since 16S sequences
// have significant length variation (1300-1550bp), terminal gaps are common
// and using gap_extend=4 over-penalizes them compared to vsearch. Using
// gap_extend=2 (matching vsearch's raw interior gap_extend) produces
// alignments with more similar gap placement to vsearch's NW, especially
// for length-mismatched sequences. This is empirically validated against
// vsearch's chimera detection output on real 16S data.
static constexpr int UCHIME_WFA2_MISMATCH = 6;
static constexpr int UCHIME_WFA2_GAP_OPEN = 20;
static constexpr int UCHIME_WFA2_GAP_EXTEND = 2;
static constexpr int SMOOTHING_WINDOW = 32;

// Result of breakpoint sweep.
struct BreakpointResult {
	double best_h = 0.0;
	int best_pos = -1; // alignment column index of best breakpoint
	bool reversed = false;
	int left_yes = 0, left_no = 0, left_abstain = 0;
	int right_yes = 0, right_no = 0, right_abstain = 0;
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

	// Per-segment candidate search with alignment filtering (matching vsearch's flow).
	// For each of 4 segments: k-mer search → align top candidates to segment →
	// accept only those with >= min_identity (55%). Dedup accepted hits across segments.
	static constexpr double CHIMERA_MIN_ID = 0.55; // 55% identity threshold for candidate acceptance
	static constexpr int MAX_ACCEPTS_PER_SEGMENT = 4;
	static constexpr int MAX_REJECTS_PER_SEGMENT = 16;
	std::vector<uint32_t> find_chimera_candidates(const std::string &query, const std::string &masked_query,
	                                              WFA2Aligner &aligner) const;

	// Compute identity % between two aligned sequences over non-ignored columns.
	static double compute_identity(const std::string &aligned_a, const std::string &aligned_b);

	// Compute model identity: query matches model (parent A left of breakpoint, B right).
	static double compute_model_identity(const StarAlignment &star, int breakpoint, bool reversed);
};

} // namespace miint
