#pragma once

#include <cstdint>
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

// Build a 3-way star alignment from two pairwise alignments.
// qa/sa are query-aligned and subject-aligned strings from align_full(query, parentA).
// qb/sb are query-aligned and subject-aligned strings from align_full(query, parentB).
// The query acts as the hub: both pairwise alignments are merged by walking query
// positions in lockstep, inserting gap columns where one pairwise has gaps the other doesn't.
StarAlignment build_star_alignment(const std::string &query_aligned_a, const std::string &subject_aligned_a,
                                   const std::string &query_aligned_b, const std::string &subject_aligned_b);

// Classify each column in a star alignment into diff types.
// Fills star.diffs with one entry per column.
void classify_diffs(StarAlignment &star);

// WFA2 penalty parameters equivalent to vsearch's NW scoring.
//
// vsearch uses a score-maximization model: match=+2, mismatch=-4,
// gap_open=-20, gap_extend=-2. A gap of length k costs: gap_open + k*gap_extend.
//
// WFA2 uses a penalty-minimization model: mismatch, gap_open, gap_extend are
// positive penalties. A gap of length k costs: gap_open + k*gap_extend.
// Both use the same affine gap convention.
//
// Conversion: in vsearch, choosing mismatch over match costs (match - mismatch) = 6.
// Choosing gap_extend over match costs (match + |gap_extend|) = 4.
// gap_open is a one-time penalty: |gap_open| = 20.
//
// Note: our test sequences are substitution-only (no indels), so these gap
// penalties haven't been validated against vsearch's alignment for gapped cases.
// If gapped alignments diverge, this conversion should be revisited.
static constexpr int UCHIME_WFA2_MISMATCH = 6;
static constexpr int UCHIME_WFA2_GAP_OPEN = 20;
static constexpr int UCHIME_WFA2_GAP_EXTEND = 4;

} // namespace miint
