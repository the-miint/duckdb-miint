#pragma once

#include <memory>
#include <optional>
#include <string>

namespace miint {

struct WFA2CigarResult {
	int score;
	std::string cigar;
};

struct WFA2FullResult {
	int score;
	std::string cigar;
	std::string query_aligned;
	std::string subject_aligned;
};

// Thread-safety: WFA2Aligner instances are NOT thread-safe for concurrent calls.
// In DuckDB scalar functions, each thread gets its own instance via FunctionLocalState.
// Reusing an instance across rows within a single thread is safe and recommended.
class WFA2Aligner {
public:
	// Default penalties: mismatch=4, gap_open=6, gap_extend=2
	// These match BLASTN defaults and are standard for DNA sequence alignment
	// (Altschul et al., "Gapped BLAST and PSI-BLAST", NAR 1997).
	WFA2Aligner();
	WFA2Aligner(int mismatch, int gap_open, int gap_extend);
	~WFA2Aligner();

	// Non-copyable, movable
	WFA2Aligner(const WFA2Aligner &) = delete;
	WFA2Aligner &operator=(const WFA2Aligner &) = delete;
	WFA2Aligner(WFA2Aligner &&) noexcept;
	WFA2Aligner &operator=(WFA2Aligner &&) noexcept;

	// Returns nullopt on alignment failure (StatusMaxStepsReached, StatusOOM)
	std::optional<int> align_score(const std::string &query, const std::string &subject);
	std::optional<WFA2CigarResult> align_cigar(const std::string &query, const std::string &subject);
	std::optional<WFA2FullResult> align_full(const std::string &query, const std::string &subject);

private:
	struct Impl;
	std::unique_ptr<Impl> impl_;

	// Reconstruct aligned sequences from CIGAR ops
	static void reconstruct_aligned(const std::string &query, const std::string &subject, const std::string &cigar,
	                                std::string &query_aligned, std::string &subject_aligned);
};

} // namespace miint
