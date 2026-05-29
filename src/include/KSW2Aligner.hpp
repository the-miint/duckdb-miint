#pragma once

#include <memory>
#include <optional>
#include <string>

namespace miint {

struct KSW2CigarResult {
	int score;
	std::string cigar;
};

struct KSW2FullResult {
	int score;
	std::string cigar;
	std::string query_aligned;
	std::string subject_aligned;
};

// Thread-safety: KSW2Aligner instances are NOT thread-safe for concurrent calls.
// In DuckDB scalar functions, each thread gets its own instance via FunctionLocalState.
// Reusing an instance across rows within a single thread is safe and recommended.
//
// Scoring semantics: KSW2 returns native positive alignment scores (match contributes
// positively, mismatch + gap subtract). This differs from WFA2Aligner, which returns
// gap-affine penalties (identical = 0). Function names (_wfa2 vs _ksw2) telegraph the
// difference; do not compare scores across the two backends.
//
// CIGAR: KSW2 produces standard ops (M = match-or-mismatch lumped; I, D, N_SKIP).
// WFA2 produces extended ops (= and X distinguish match vs mismatch).
class KSW2Aligner {
public:
	// Defaults: match=2, mismatch=4, gap_open=6, gap_extend=2, w=-1 (no band), zdrop=-1 (disabled).
	// These mirror BLASTN-ish defaults symmetric with WFA2Aligner; w=-1/zdrop=-1 give exact end-to-end
	// alignment, matching WFA2 semantics.
	KSW2Aligner();
	KSW2Aligner(int match, int mismatch, int gap_open, int gap_extend, int w = -1, int zdrop = -1);

	// Dual-affine ctor: stores a second gap-penalty pair (gap_open2, gap_extend2) used by the
	// align_extd_* methods. All 8 args required to avoid arity ambiguity with the extz ctor above;
	// callers wanting "default" w/zdrop should pass -1, -1 explicitly.
	KSW2Aligner(int match, int mismatch, int gap_open, int gap_extend, int gap_open2, int gap_extend2, int w,
	            int zdrop);

	// Splice ctor: stores intron-open penalty (gap_open2) + non-canonical splice penalty (noncan)
	// used by the align_exts_* methods. All 7 args required to avoid arity ambiguity with the
	// other ctors. There is no `w` parameter — ksw_exts2_sse does not take a bandwidth.
	KSW2Aligner(int match, int mismatch, int gap_open, int gap_extend, int gap_open2, int noncan, int zdrop);

	~KSW2Aligner();

	KSW2Aligner(const KSW2Aligner &) = delete;
	KSW2Aligner &operator=(const KSW2Aligner &) = delete;
	KSW2Aligner(KSW2Aligner &&) noexcept;
	KSW2Aligner &operator=(KSW2Aligner &&) noexcept;

	// extz mode (ksw_extz2_sse): standard affine extension, end-to-end scored via ez->score.
	// Returns nullopt if the alignment fails to reach both ends (ez->score == KSW_NEG_INF).
	std::optional<int> align_extz_score(const std::string &query, const std::string &subject);
	std::optional<KSW2CigarResult> align_extz_cigar(const std::string &query, const std::string &subject);
	std::optional<KSW2FullResult> align_extz_full(const std::string &query, const std::string &subject);

	// extd mode (ksw_extd2_sse): dual-affine extension. For each gap of length L, KSW2 picks the
	// cheaper of (gap_open + L*gap_extend) and (gap_open2 + L*gap_extend2). MUST be called on an
	// instance constructed with the dual-affine (8-arg) ctor — on any other ctor, gap_open2 and
	// gap_extend2 hold internal sentinels that do NOT yield extz-equivalent scoring (the cheap-
	// extend sentinel makes the second pair attractive for any gap >= 1bp), so behavior is
	// unspecified. Enforcement is by convention; the caller is responsible.
	std::optional<int> align_extd_score(const std::string &query, const std::string &subject);
	std::optional<KSW2CigarResult> align_extd_cigar(const std::string &query, const std::string &subject);
	std::optional<KSW2FullResult> align_extd_full(const std::string &query, const std::string &subject);

	// exts mode (ksw_exts2_sse): splice-aware extension. Uses gap_open/gap_extend for non-splice
	// gaps; gap_open2 acts as the intron-open penalty (extension still uses gap_extend). `noncan`
	// adds a penalty when the chosen intron boundary is not canonical (GT-AG). Junction guidance
	// (`junc` array) is NULL in v1 — the alignment relies only on the score model. Flag is
	// KSW_EZ_SPLICE_FOR (forward strand). MUST be called on an instance constructed with the
	// splice (7-arg) ctor; on other ctors `noncan` is a sentinel and behavior is unspecified.
	std::optional<int> align_exts_score(const std::string &query, const std::string &subject);
	std::optional<KSW2CigarResult> align_exts_cigar(const std::string &query, const std::string &subject);
	std::optional<KSW2FullResult> align_exts_full(const std::string &query, const std::string &subject);

private:
	struct Impl;
	std::unique_ptr<Impl> impl_;
};

} // namespace miint
