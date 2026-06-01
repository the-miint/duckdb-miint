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

	// Semi-global alignment: query anchored end-to-end, subject (reference) has
	// free end gaps (can overhang at both ends). This matches vsearch's chimera
	// detection alignment where terminal gaps on the target are penalized less.
	std::optional<WFA2FullResult> align_full_semiglobal(const std::string &query, const std::string &subject);

	// Semi-global with IUPAC-aware matching. Uses WFA2's lambda callback so
	// degenerate bases (N, R, Y, ...) match their compatible bases at zero cost.
	// Returns CIGAR+score only (no reconstructed aligned sequences).
	//
	// text_begin_free / text_end_free cap how many bases of the query (text) may
	// be trimmed at no cost from each end. Default 0 keeps the query anchored
	// end-to-end (legacy behavior). Pattern (subject) ends are always free up to
	// subject.size(); these knobs add free trim on the query side, enabling
	// cutadapt-style partial-overlap matches where the query's prefix/suffix
	// hangs off the subject's edge.
	std::optional<WFA2CigarResult> align_cigar_semiglobal_iupac(const std::string &query, const std::string &subject,
	                                                            int text_begin_free = 0, int text_end_free = 0);

private:
	struct Impl;
	std::unique_ptr<Impl> impl_;
};

} // namespace miint
