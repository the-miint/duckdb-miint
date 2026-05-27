#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace miint {

// A single per-read observation at one MSA column.
//   base in {A,C,G,T,-} (lower-case is upper-cased on insert by the aggregate)
//   qual = Phred score with no ASCII offset; 0 for gap observations
struct Observation {
	char base;
	std::uint8_t qual;
};

// The reduced consensus call at one column: a base + posterior-derived Q.
//   base == '-' means the column is suppressed from the output consensus.
struct ConsensusBase {
	char base;
	std::uint8_t qual;
};

// 5-state log-likelihood vote over {A,C,G,T,-} for a single MSA column.
//
// Model: P(obs | true) under a uniform error model.
//   For base observations: p_err = 10^(-q/10);
//   For gap observations:  p_err = GAP_ERR (constant);
//   P(obs | true == obs)   = 1 - p_err
//   P(obs | true == other) = p_err / 4
//
// Argmax over {A,C,G,T,-}, ties broken first by sum-of-Q for matching
// observations, then alphabetically ({A < C < G < T < -}). Posterior Q
// clamped to [0, 60] (UTINYINT) on output.
//
// Throws miint::InvalidInputException on empty observation list.
ConsensusBase VoteColumn(const std::vector<Observation> &obs);

// HP-aware post-correction: for each homopolymer run of length >= 2 in
// `consensus_seq`, look up the per-read HP length at the corresponding
// ungapped locus (derived from `aligned_seqs` + the MSA column index of
// each consensus position) and replace the run with `median(per-read HP
// lengths)` copies of the base. Reads without an HP of the consensus
// base at the corresponding ungapped locus are skipped from the median.
//
// Returns (corrected_seq, corrected_qual). Qual at every position in the
// emitted run is taken from `consensus_qual` at the run's first position
// (no new posterior synthesised here).
std::pair<std::string, std::vector<std::uint8_t>> HpCorrect(const std::string &consensus_seq,
                                                            const std::vector<std::uint8_t> &consensus_qual,
                                                            const std::vector<std::size_t> &msa_col_of_consensus_pos,
                                                            const std::vector<std::string> &aligned_seqs);

// End-to-end consensus pipeline: VoteColumn per MSA column + HpCorrect on
// the resulting naive consensus. Single-row groups bypass MSA and return
// (gap-stripped seq, qual) directly. Throws on width mismatches.
//
// Pulled into the kernel so the parallel-Combine semantic invariant
// ("partition-then-merge yields the same consensus as the single batch")
// is testable in isolation; see test/cpp/test_consensus_msa.cpp.
std::pair<std::string, std::vector<std::uint8_t>> BuildConsensus(const std::vector<std::string> &aligned_seqs,
                                                                 const std::vector<std::vector<std::uint8_t>> &quals);

} // namespace miint
