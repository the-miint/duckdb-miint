#include "WFA2Aligner.hpp"

#include <WFA2-lib/bindings/cpp/WFAligner.hpp>
#include <stdexcept>

namespace miint {

struct WFA2Aligner::Impl {
	// Two aligner instances for different scopes, both using BiWFA (MemoryUltralow)
	// for O(s) memory instead of O(s^2):
	// - Score-scope: faster, only returns alignment score
	// - Alignment-scope: computes CIGAR traceback (but NOT used for score — see below)
	//
	// WFA2-lib v2.3.5 bug workaround: the BiWFA alignment-scope path does not
	// populate cigar->score when it falls back to unidirectional alignment for
	// sequences <= WF_BIALIGN_FALLBACK_MIN_LENGTH (100bp) or scores <=
	// WF_BIALIGN_FALLBACK_MIN_SCORE (250). The fallback in wavefront_bialign_base()
	// correctly appends CIGAR operations but returns WF_STATUS_OK without propagating
	// the base aligner's score to the parent cigar. The score remains at INT32_MIN.
	// See ext/WFA2-lib/bialign_score_bug.c for a standalone reproduction.
	//
	// Workaround: always obtain the score from the score-scope aligner (which does
	// not have this bug), and only use the alignment-scope aligner for CIGAR retrieval.
	std::unique_ptr<wfa::WFAlignerGapAffine> score_aligner;
	std::unique_ptr<wfa::WFAlignerGapAffine> alignment_aligner;

	Impl(int mismatch, int gap_open, int gap_extend) {
		score_aligner = std::make_unique<wfa::WFAlignerGapAffine>(mismatch, gap_open, gap_extend, wfa::WFAligner::Score,
		                                                          wfa::WFAligner::MemoryUltralow);
		alignment_aligner = std::make_unique<wfa::WFAlignerGapAffine>(
		    mismatch, gap_open, gap_extend, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryUltralow);
	}
};

WFA2Aligner::WFA2Aligner() : WFA2Aligner(4, 6, 2) {
}

WFA2Aligner::WFA2Aligner(int mismatch, int gap_open, int gap_extend) {
	if (mismatch <= 0) {
		throw std::invalid_argument("mismatch must be > 0");
	}
	if (gap_open < 0) {
		throw std::invalid_argument("gap_open must be >= 0");
	}
	if (gap_extend <= 0) {
		throw std::invalid_argument("gap_extend must be > 0");
	}
	impl_ = std::make_unique<Impl>(mismatch, gap_open, gap_extend);
}

WFA2Aligner::~WFA2Aligner() = default;
WFA2Aligner::WFA2Aligner(WFA2Aligner &&) noexcept = default;
WFA2Aligner &WFA2Aligner::operator=(WFA2Aligner &&) noexcept = default;

std::optional<int> WFA2Aligner::align_score(const std::string &query, const std::string &subject) {
	if (query.empty() && subject.empty()) {
		return 0;
	}
	// One empty, one non-empty: WFA2 handles correctly as all-gap alignment
	// (score = gap_open + gap_extend * len, CIGAR = <len>I or <len>D).
	// WFA2 convention: pattern=subject, text=query. This matches SAM convention
	// where subject is the reference and query is the read.
	auto status = impl_->score_aligner->alignEnd2End(subject, query);
	if (status != wfa::WFAligner::StatusAlgCompleted) {
		return std::nullopt;
	}
	// WFA2 returns negative scores (penalties as negative costs). We negate to
	// return positive values where lower score = more similar.
	return -(impl_->score_aligner->getAlignmentScore());
}

std::optional<WFA2CigarResult> WFA2Aligner::align_cigar(const std::string &query, const std::string &subject) {
	if (query.empty() && subject.empty()) {
		return WFA2CigarResult {0, ""};
	}
	// Run both aligners: score-scope for reliable score, alignment-scope for CIGAR.
	// See Impl comment for why we don't trust alignment_aligner->getAlignmentScore().
	auto score_status = impl_->score_aligner->alignEnd2End(subject, query);
	if (score_status != wfa::WFAligner::StatusAlgCompleted) {
		return std::nullopt;
	}
	auto align_status = impl_->alignment_aligner->alignEnd2End(subject, query);
	if (align_status != wfa::WFAligner::StatusAlgCompleted) {
		return std::nullopt;
	}
	WFA2CigarResult result;
	result.score = -(impl_->score_aligner->getAlignmentScore());
	// getCIGAR(true) produces extended CIGAR with =/X ops (not M).
	// This is more informative and SAM-spec compliant.
	result.cigar = impl_->alignment_aligner->getCIGAR(true);
	return result;
}

std::optional<WFA2FullResult> WFA2Aligner::align_full(const std::string &query, const std::string &subject) {
	if (query.empty() && subject.empty()) {
		return WFA2FullResult {0, "", "", ""};
	}
	// Run both aligners: score-scope for reliable score, alignment-scope for CIGAR.
	// See Impl comment for why we don't trust alignment_aligner->getAlignmentScore().
	auto score_status = impl_->score_aligner->alignEnd2End(subject, query);
	if (score_status != wfa::WFAligner::StatusAlgCompleted) {
		return std::nullopt;
	}
	auto align_status = impl_->alignment_aligner->alignEnd2End(subject, query);
	if (align_status != wfa::WFAligner::StatusAlgCompleted) {
		return std::nullopt;
	}
	WFA2FullResult result;
	result.score = -(impl_->score_aligner->getAlignmentScore());
	result.cigar = impl_->alignment_aligner->getCIGAR(true);
	reconstruct_aligned(query, subject, result.cigar, result.query_aligned, result.subject_aligned);
	return result;
}

// Reconstruct gapped alignment strings from an extended CIGAR and the original sequences.
//
// This is intentionally a separate implementation from ParseCigar() in
// alignment_functions_internal.hpp. ParseCigar() accumulates aggregate statistics
// (total match count, gap count, etc.) and never touches the underlying sequences.
// Here we need a positional walk: we consume characters from query and subject in
// lockstep with CIGAR ops, inserting gap characters ('-') to produce the aligned
// representation. The two loops have fundamentally different structure and output,
// so sharing code would add coupling without reducing complexity.
//
// We also intentionally avoid HTSlib's CIGAR utilities here. HTSlib operates on
// binary uint32_t arrays from bam1_t records; WFA2 produces string CIGARs. Converting
// string -> HTSlib binary -> positional walk would add an unnecessary dependency and
// extra conversion step for a straightforward ~30-line function.
void WFA2Aligner::reconstruct_aligned(const std::string &query, const std::string &subject, const std::string &cigar,
                                      std::string &query_aligned, std::string &subject_aligned) {
	query_aligned.clear();
	subject_aligned.clear();
	query_aligned.reserve(query.size() + subject.size());
	subject_aligned.reserve(query.size() + subject.size());

	size_t qi = 0; // query position
	size_t si = 0; // subject position
	size_t ci = 0; // cigar string position

	while (ci < cigar.size()) {
		// Parse run length prefix (e.g. "3" in "3=")
		size_t num_start = ci;
		while (ci < cigar.size() && cigar[ci] >= '0' && cigar[ci] <= '9') {
			ci++;
		}
		if (ci >= cigar.size()) {
			throw std::runtime_error("CIGAR string ends with digits but no operation character");
		}
		int count = 1;
		if (ci > num_start) {
			try {
				count = std::stoi(cigar.substr(num_start, ci - num_start));
			} catch (const std::out_of_range &) {
				throw std::runtime_error("CIGAR operation length overflows integer");
			}
			if (count <= 0) {
				throw std::runtime_error("CIGAR operation length must be positive");
			}
		}

		char op = cigar[ci++];
		for (int k = 0; k < count; k++) {
			switch (op) {
			case '=': // sequence match
			case 'X': // sequence mismatch
			case 'M': // match or mismatch
				if (qi >= query.size() || si >= subject.size()) {
					throw std::runtime_error("CIGAR consumes more bases than available in sequences");
				}
				query_aligned += query[qi++];
				subject_aligned += subject[si++];
				break;
			case 'I': // insertion in query (query has extra bases vs subject)
				if (qi >= query.size()) {
					throw std::runtime_error("CIGAR consumes more query bases than available");
				}
				query_aligned += query[qi++];
				subject_aligned += '-';
				break;
			case 'D': // deletion from query (subject has extra bases vs query)
				if (si >= subject.size()) {
					throw std::runtime_error("CIGAR consumes more subject bases than available");
				}
				query_aligned += '-';
				subject_aligned += subject[si++];
				break;
			default:
				throw std::runtime_error(std::string("Unknown CIGAR operation: ") + op);
			}
		}
	}

	if (qi != query.size()) {
		throw std::runtime_error("CIGAR did not consume all query bases: consumed " + std::to_string(qi) + " of " +
		                         std::to_string(query.size()));
	}
	if (si != subject.size()) {
		throw std::runtime_error("CIGAR did not consume all subject bases: consumed " + std::to_string(si) + " of " +
		                         std::to_string(subject.size()));
	}
}

} // namespace miint
