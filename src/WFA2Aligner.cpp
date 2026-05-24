#include "WFA2Aligner.hpp"
#include "sequence_utils.hpp"

#include <WFA2-lib/bindings/cpp/WFAligner.hpp>
#include <stdexcept>

namespace miint {

struct WFA2Aligner::Impl {
	// Score-scope aligner: faster path for align_score() when CIGAR is not needed.
	std::unique_ptr<wfa::WFAlignerGapAffine> score_aligner;

	// Alignment-scope aligner: computes CIGAR traceback + score. Uses MemoryMed
	// (piggyback backtrace) — near-MemoryHigh speed with bounded memory.
	std::unique_ptr<wfa::WFAlignerGapAffine> alignment_aligner;

	Impl(int mismatch, int gap_open, int gap_extend) {
		score_aligner = std::make_unique<wfa::WFAlignerGapAffine>(mismatch, gap_open, gap_extend, wfa::WFAligner::Score,
		                                                          wfa::WFAligner::MemoryMed);
		alignment_aligner = std::make_unique<wfa::WFAlignerGapAffine>(
		    mismatch, gap_open, gap_extend, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryMed);
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
	auto status = impl_->score_aligner->alignEnd2End(subject, query);
	if (status != wfa::WFAligner::StatusAlgCompleted) {
		return std::nullopt;
	}
	return -(impl_->score_aligner->getAlignmentScore());
}

std::optional<WFA2CigarResult> WFA2Aligner::align_cigar(const std::string &query, const std::string &subject) {
	if (query.empty() && subject.empty()) {
		return WFA2CigarResult {0, ""};
	}

	auto align_status = impl_->alignment_aligner->alignEnd2End(subject, query);
	if (align_status != wfa::WFAligner::StatusAlgCompleted) {
		return std::nullopt;
	}

	WFA2CigarResult result;
	result.cigar = impl_->alignment_aligner->getCIGAR(true);
	result.score = -(impl_->alignment_aligner->getAlignmentScore());
	return result;
}

std::optional<WFA2FullResult> WFA2Aligner::align_full(const std::string &query, const std::string &subject) {
	if (query.empty() && subject.empty()) {
		return WFA2FullResult {0, "", "", ""};
	}

	auto align_status = impl_->alignment_aligner->alignEnd2End(subject, query);
	if (align_status != wfa::WFAligner::StatusAlgCompleted) {
		return std::nullopt;
	}

	WFA2FullResult result;
	result.cigar = impl_->alignment_aligner->getCIGAR(true);
	result.score = -(impl_->alignment_aligner->getAlignmentScore());
	reconstruct_aligned(query, subject, result.cigar, result.query_aligned, result.subject_aligned);
	return result;
}

std::optional<WFA2FullResult> WFA2Aligner::align_full_semiglobal(const std::string &query, const std::string &subject) {
	if (query.empty() && subject.empty()) {
		return WFA2FullResult {0, "", "", ""};
	}

	// Semi-global: query (text) anchored end-to-end, subject (pattern) has free end gaps.
	// WFA2 convention: pattern=subject, text=query.
	// patternBeginFree/patternEndFree = subject length (all free)
	// textBeginFree/textEndFree = 0 (query anchored)
	auto align_status = impl_->alignment_aligner->alignEndsFree(subject, static_cast<int>(subject.size()),
	                                                            static_cast<int>(subject.size()), query, 0, 0);
	if (align_status != wfa::WFAligner::StatusAlgCompleted) {
		return std::nullopt;
	}

	WFA2FullResult result;
	result.cigar = impl_->alignment_aligner->getCIGAR(true);
	result.score = -(impl_->alignment_aligner->getAlignmentScore());
	reconstruct_aligned(query, subject, result.cigar, result.query_aligned, result.subject_aligned);
	return result;
}

namespace {
struct IupacMatchArgs {
	const char *pattern_data; // WFA2 pattern = subject
	int pattern_len;
	const char *text_data; // WFA2 text = query
	int text_len;
};

// Defensive OOB guard per WFA2 documentation (wavefront_sequences.h:46-52).
// WFA2 should never pass OOB indices for valid inputs, so this is a safety
// net — returning 0 (mismatch) is conservative and avoids UB on a library bug.
int iupac_match_callback(int pattern_pos, int text_pos, void *args) {
	auto *ctx = static_cast<IupacMatchArgs *>(args);
	if (pattern_pos >= ctx->pattern_len || text_pos >= ctx->text_len) {
		return 0;
	}
	return IupacMatch(ctx->pattern_data[pattern_pos], ctx->text_data[text_pos]) ? 1 : 0;
}
} // namespace

std::optional<WFA2CigarResult> WFA2Aligner::align_cigar_semiglobal_iupac(const std::string &query,
                                                                         const std::string &subject) {
	if (query.empty() && subject.empty()) {
		return WFA2CigarResult {0, ""};
	}

	IupacMatchArgs args {subject.data(), static_cast<int>(subject.size()), query.data(),
	                     static_cast<int>(query.size())};
	auto align_status = impl_->alignment_aligner->alignEndsFree(
	    iupac_match_callback, &args, static_cast<int>(subject.size()), static_cast<int>(subject.size()),
	    static_cast<int>(subject.size()), static_cast<int>(query.size()), 0, 0);
	if (align_status != wfa::WFAligner::StatusAlgCompleted) {
		return std::nullopt;
	}

	WFA2CigarResult result;
	result.cigar = impl_->alignment_aligner->getCIGAR(true);
	result.score = -(impl_->alignment_aligner->getAlignmentScore());
	return result;
}

// Reconstruct gapped alignment strings from an extended CIGAR and the original sequences.
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
