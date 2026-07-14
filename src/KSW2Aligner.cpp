#include "KSW2Aligner.hpp"

#include "cigar_reconstruction.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>

extern "C" {
#include <minimap2/ksw2.h>
extern unsigned char seq_nt4_table[256];
}

namespace miint {

namespace {

// Decode a packed KSW2 CIGAR array (uint32_t with op in low 4 bits, length in high 28 bits)
// into the human-readable run-length string used elsewhere (e.g., "4M2I3D").
// KSW2 only emits op codes 0-3 (M/I/D/N_SKIP); the BAM op character set is "MIDN".
std::string decode_ksw2_cigar(const uint32_t *cigar, int n_cigar) {
	static constexpr const char *kOpCharsBam = "MIDN";
	std::string out;
	for (int i = 0; i < n_cigar; i++) {
		int op = cigar[i] & 0xf;
		int len = cigar[i] >> 4;
		// op is masked with 0xf so the only reachable failure is the upper bound.
		if (op > 3) {
			throw std::runtime_error("KSW2 produced unexpected CIGAR op code " + std::to_string(op));
		}
		out += std::to_string(len);
		out += kOpCharsBam[op];
	}
	return out;
}

// Score for an all-gap alignment of length N: -(gap_open + N * gap_extend).
int all_gap_score(int gap_open, int gap_extend, size_t n) {
	return -(gap_open + static_cast<int>(n) * gap_extend);
}

} // namespace

// Growable byte buffer that grows capacity without value-initializing — used by encode_nt4
// to skip the memset that std::vector<uint8_t>::resize() incurs on growth. Capacity is held
// for the lifetime of the owning KSW2Aligner::Impl, so per-row encoding pays no allocation
// after the first oversize sequence.
struct EncBuf {
	std::unique_ptr<uint8_t[]> data;
	size_t capacity = 0;

	void reserve(size_t n) {
		if (n > capacity) {
			data.reset(new uint8_t[n]);
			capacity = n;
		}
	}
	uint8_t *get() {
		return data.get();
	}
};

struct KSW2Aligner::Impl {
	int gap_open;
	int gap_extend;
	int gap_open2;
	int gap_extend2;
	int noncan;
	int w;
	int zdrop;
	int8_t mat[25]; // 5x5 substitution matrix (ACGTN)

	// Reused encoding buffers — keep their capacity across align_* calls so per-row alignment
	// does not pay a fresh heap allocation. Pre-reserved for typical Illumina read length (256bp).
	// KSW2Aligner is documented non-thread-safe so a per-instance buffer is correct.
	EncBuf q_enc_buf;
	EncBuf s_enc_buf;

	// Equivalent of minimap2's static `ksw_gen_simple_mat` (in align.c), which is not exported
	// from libminimap2.a — reproduced here so the KSW2 substitution matrix has a sensible default
	// without modifying the vendored minimap2 source.
	// Some fields are mode-specific: gap_open2/gap_extend2 are read by align_extd_*; noncan is
	// read by align_exts_*. For modes that do not use a given field, the ctor caller passes a
	// documented sentinel and that field is ignored by the called kernel.
	Impl(int match, int mismatch, int gap_open_p, int gap_extend_p, int gap_open2_p, int gap_extend2_p, int noncan_p,
	     int w_p, int zdrop_p)
	    : gap_open(gap_open_p), gap_extend(gap_extend_p), gap_open2(gap_open2_p), gap_extend2(gap_extend2_p),
	      noncan(noncan_p), w(w_p), zdrop(zdrop_p) {
		q_enc_buf.reserve(256);
		s_enc_buf.reserve(256);
		// Diagonal among ACGT (0..3): +match. Off-diagonal among ACGT: -mismatch.
		// Any row or column for N (4): 0, making N a wildcard against any base.
		for (int i = 0; i < 5; i++) {
			for (int j = 0; j < 5; j++) {
				if (i == 4 || j == 4) {
					mat[i * 5 + j] = 0;
				} else if (i == j) {
					mat[i * 5 + j] = static_cast<int8_t>(match);
				} else {
					mat[i * 5 + j] = static_cast<int8_t>(-mismatch);
				}
			}
		}
	}

	// Encode an ASCII sequence into nt4 codes (A->0, C->1, G->2, T->3, anything else -> 4).
	// Returns the encoded buffer pointer, valid until the next call on the same EncBuf.
	static uint8_t *encode_nt4(const std::string &s, EncBuf &out) {
		out.reserve(s.size());
		uint8_t *p = out.get();
		for (size_t i = 0; i < s.size(); i++) {
			p[i] = seq_nt4_table[static_cast<uint8_t>(s[i])];
		}
		return p;
	}
};

namespace {

// KSW2 stores penalty / scoring params as int8_t (and the 5x5 substitution matrix is int8_t too).
// Reject values that would silently truncate when cast to int8_t. The non-negativity checks already
// restrict to [0, INT_MAX]; the upper bound here is INT8_MAX (127).
static constexpr int KSW2_INT8_MAX = 127;

void validate_common(int match, int mismatch, int gap_open, int gap_extend) {
	if (match <= 0 || match > KSW2_INT8_MAX) {
		throw std::invalid_argument("match must be in (0, 127]");
	}
	if (mismatch <= 0 || mismatch > KSW2_INT8_MAX) {
		throw std::invalid_argument("mismatch must be in (0, 127]");
	}
	if (gap_open < 0 || gap_open > KSW2_INT8_MAX) {
		throw std::invalid_argument("gap_open must be in [0, 127]");
	}
	if (gap_extend <= 0 || gap_extend > KSW2_INT8_MAX) {
		throw std::invalid_argument("gap_extend must be in (0, 127]");
	}
}

void validate_dual_affine(int gap_open2, int gap_extend2) {
	if (gap_open2 < 0 || gap_open2 > KSW2_INT8_MAX) {
		throw std::invalid_argument("gap_open2 must be in [0, 127]");
	}
	if (gap_extend2 <= 0 || gap_extend2 > KSW2_INT8_MAX) {
		throw std::invalid_argument("gap_extend2 must be in (0, 127]");
	}
}

void validate_splice(int gap_open2, int noncan) {
	if (gap_open2 < 0 || gap_open2 > KSW2_INT8_MAX) {
		throw std::invalid_argument("gap_open2 must be in [0, 127]");
	}
	if (noncan < 0 || noncan > KSW2_INT8_MAX) {
		throw std::invalid_argument("noncan must be in [0, 127]");
	}
}

} // namespace

KSW2Aligner::KSW2Aligner() : KSW2Aligner(2, 4, 6, 2, -1, -1) {
}

KSW2Aligner::KSW2Aligner(int match, int mismatch, int gap_open, int gap_extend, int w, int zdrop) {
	validate_common(match, mismatch, gap_open, gap_extend);
	// Placeholder values for fields not used by extz (gap_open2=0, gap_extend2=1, noncan=0).
	// Calling align_extd_*/align_exts_* on an extz-ctor instance is unspecified — see the
	// contract in KSW2Aligner.hpp. These values are not load-bearing.
	impl_ = std::make_unique<Impl>(match, mismatch, gap_open, gap_extend, /*gap_open2=*/0, /*gap_extend2=*/1,
	                               /*noncan=*/0, w, zdrop);
}

KSW2Aligner::KSW2Aligner(int match, int mismatch, int gap_open, int gap_extend, int gap_open2, int gap_extend2, int w,
                         int zdrop) {
	validate_common(match, mismatch, gap_open, gap_extend);
	validate_dual_affine(gap_open2, gap_extend2);
	impl_ =
	    std::make_unique<Impl>(match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, /*noncan=*/0, w, zdrop);
}

KSW2Aligner::KSW2Aligner(int match, int mismatch, int gap_open, int gap_extend, int gap_open2, int noncan, int zdrop) {
	validate_common(match, mismatch, gap_open, gap_extend);
	validate_splice(gap_open2, noncan);
	// Placeholder values for fields not used by exts: gap_extend2=1 (exts uses gap_extend for
	// intron extension), w=-1 (ksw_exts2_sse has no bandwidth parameter). Calling align_extd_*
	// on a splice-ctor instance is unspecified — see the contract in KSW2Aligner.hpp.
	impl_ = std::make_unique<Impl>(match, mismatch, gap_open, gap_extend, gap_open2, /*gap_extend2=*/1, noncan,
	                               /*w=*/-1, zdrop);
}

KSW2Aligner::~KSW2Aligner() = default;
KSW2Aligner::KSW2Aligner(KSW2Aligner &&) noexcept = default;
KSW2Aligner &KSW2Aligner::operator=(KSW2Aligner &&) noexcept = default;

std::optional<int> KSW2Aligner::align_extz_score(const std::string &query, const std::string &subject) {
	if (query.empty() && subject.empty()) {
		return 0;
	}
	if (query.empty() || subject.empty()) {
		size_t n = std::max(query.size(), subject.size());
		return all_gap_score(impl_->gap_open, impl_->gap_extend, n);
	}

	uint8_t *q_enc = Impl::encode_nt4(query, impl_->q_enc_buf);
	uint8_t *s_enc = Impl::encode_nt4(subject, impl_->s_enc_buf);

	ksw_extz_t ez {};
	ksw_extz2_sse(nullptr, static_cast<int>(query.size()), q_enc, static_cast<int>(subject.size()), s_enc, 5,
	              impl_->mat, static_cast<int8_t>(impl_->gap_open), static_cast<int8_t>(impl_->gap_extend), impl_->w,
	              impl_->zdrop, 0 /* end_bonus */, KSW_EZ_SCORE_ONLY, &ez);

	std::optional<int> result;
	if (ez.score != KSW_NEG_INF) {
		result = ez.score;
	}
	kfree(nullptr, ez.cigar);
	return result;
}

std::optional<KSW2CigarResult> KSW2Aligner::align_extz_cigar(const std::string &query, const std::string &subject) {
	if (query.empty() && subject.empty()) {
		return KSW2CigarResult {0, ""};
	}
	if (query.empty()) {
		return KSW2CigarResult {all_gap_score(impl_->gap_open, impl_->gap_extend, subject.size()),
		                        std::to_string(subject.size()) + "D"};
	}
	if (subject.empty()) {
		return KSW2CigarResult {all_gap_score(impl_->gap_open, impl_->gap_extend, query.size()),
		                        std::to_string(query.size()) + "I"};
	}

	uint8_t *q_enc = Impl::encode_nt4(query, impl_->q_enc_buf);
	uint8_t *s_enc = Impl::encode_nt4(subject, impl_->s_enc_buf);

	ksw_extz_t ez {};
	ksw_extz2_sse(nullptr, static_cast<int>(query.size()), q_enc, static_cast<int>(subject.size()), s_enc, 5,
	              impl_->mat, static_cast<int8_t>(impl_->gap_open), static_cast<int8_t>(impl_->gap_extend), impl_->w,
	              impl_->zdrop, 0, 0 /* flag: produce CIGAR */, &ez);

	std::optional<KSW2CigarResult> result;
	if (ez.score != KSW_NEG_INF) {
		// eqx post-pass: KSW2 emits only 'M'; split into '='/'X' so identity is readable
		// from the CIGAR alone (see eqx_split_cigar in cigar_reconstruction.hpp).
		result = KSW2CigarResult {ez.score, eqx_split_cigar(decode_ksw2_cigar(ez.cigar, ez.n_cigar), query, subject)};
	}
	kfree(nullptr, ez.cigar);
	return result;
}

std::optional<KSW2FullResult> KSW2Aligner::align_extz_full(const std::string &query, const std::string &subject) {
	auto cigar_result = align_extz_cigar(query, subject);
	if (!cigar_result.has_value()) {
		return std::nullopt;
	}

	KSW2FullResult result;
	result.score = cigar_result->score;
	result.cigar = std::move(cigar_result->cigar);
	reconstruct_aligned_from_cigar(query, subject, result.cigar, result.query_aligned, result.subject_aligned);
	return result;
}

// ---------------------------------------------------------------------------
// Dual-affine mode (extd): ksw_extd2_sse
// ---------------------------------------------------------------------------

namespace {

// All-gap score for dual-affine: picks the cheaper of the two affine pairs.
int dual_affine_all_gap_score(int gap_open, int gap_extend, int gap_open2, int gap_extend2, size_t n) {
	int g1 = gap_open + static_cast<int>(n) * gap_extend;
	int g2 = gap_open2 + static_cast<int>(n) * gap_extend2;
	return -std::min(g1, g2);
}

} // namespace

std::optional<int> KSW2Aligner::align_extd_score(const std::string &query, const std::string &subject) {
	if (query.empty() && subject.empty()) {
		return 0;
	}
	if (query.empty() || subject.empty()) {
		size_t n = std::max(query.size(), subject.size());
		return dual_affine_all_gap_score(impl_->gap_open, impl_->gap_extend, impl_->gap_open2, impl_->gap_extend2, n);
	}

	uint8_t *q_enc = Impl::encode_nt4(query, impl_->q_enc_buf);
	uint8_t *s_enc = Impl::encode_nt4(subject, impl_->s_enc_buf);

	ksw_extz_t ez {};
	ksw_extd2_sse(nullptr, static_cast<int>(query.size()), q_enc, static_cast<int>(subject.size()), s_enc, 5,
	              impl_->mat, static_cast<int8_t>(impl_->gap_open), static_cast<int8_t>(impl_->gap_extend),
	              static_cast<int8_t>(impl_->gap_open2), static_cast<int8_t>(impl_->gap_extend2), impl_->w,
	              impl_->zdrop, 0 /* end_bonus */, KSW_EZ_SCORE_ONLY, &ez);

	std::optional<int> result;
	if (ez.score != KSW_NEG_INF) {
		result = ez.score;
	}
	kfree(nullptr, ez.cigar);
	return result;
}

std::optional<KSW2CigarResult> KSW2Aligner::align_extd_cigar(const std::string &query, const std::string &subject) {
	if (query.empty() && subject.empty()) {
		return KSW2CigarResult {0, ""};
	}
	if (query.empty()) {
		return KSW2CigarResult {dual_affine_all_gap_score(impl_->gap_open, impl_->gap_extend, impl_->gap_open2,
		                                                  impl_->gap_extend2, subject.size()),
		                        std::to_string(subject.size()) + "D"};
	}
	if (subject.empty()) {
		return KSW2CigarResult {dual_affine_all_gap_score(impl_->gap_open, impl_->gap_extend, impl_->gap_open2,
		                                                  impl_->gap_extend2, query.size()),
		                        std::to_string(query.size()) + "I"};
	}

	uint8_t *q_enc = Impl::encode_nt4(query, impl_->q_enc_buf);
	uint8_t *s_enc = Impl::encode_nt4(subject, impl_->s_enc_buf);

	ksw_extz_t ez {};
	ksw_extd2_sse(nullptr, static_cast<int>(query.size()), q_enc, static_cast<int>(subject.size()), s_enc, 5,
	              impl_->mat, static_cast<int8_t>(impl_->gap_open), static_cast<int8_t>(impl_->gap_extend),
	              static_cast<int8_t>(impl_->gap_open2), static_cast<int8_t>(impl_->gap_extend2), impl_->w,
	              impl_->zdrop, 0, 0 /* flag: produce CIGAR */, &ez);

	std::optional<KSW2CigarResult> result;
	if (ez.score != KSW_NEG_INF) {
		// eqx post-pass: KSW2 emits only 'M'; split into '='/'X' so identity is readable
		// from the CIGAR alone (see eqx_split_cigar in cigar_reconstruction.hpp).
		result = KSW2CigarResult {ez.score, eqx_split_cigar(decode_ksw2_cigar(ez.cigar, ez.n_cigar), query, subject)};
	}
	kfree(nullptr, ez.cigar);
	return result;
}

std::optional<KSW2FullResult> KSW2Aligner::align_extd_full(const std::string &query, const std::string &subject) {
	auto cigar_result = align_extd_cigar(query, subject);
	if (!cigar_result.has_value()) {
		return std::nullopt;
	}

	KSW2FullResult result;
	result.score = cigar_result->score;
	result.cigar = std::move(cigar_result->cigar);
	reconstruct_aligned_from_cigar(query, subject, result.cigar, result.query_aligned, result.subject_aligned);
	return result;
}

// ---------------------------------------------------------------------------
// Splice-aware mode (exts): ksw_exts2_sse
// ---------------------------------------------------------------------------

std::optional<int> KSW2Aligner::align_exts_score(const std::string &query, const std::string &subject) {
	if (query.empty() && subject.empty()) {
		return 0;
	}
	if (query.empty() || subject.empty()) {
		// Splice mode treats one-empty as a plain all-gap alignment using the first affine pair
		// (gap_open2 is for intron-open, not whole-sequence absence).
		size_t n = std::max(query.size(), subject.size());
		return all_gap_score(impl_->gap_open, impl_->gap_extend, n);
	}

	uint8_t *q_enc = Impl::encode_nt4(query, impl_->q_enc_buf);
	uint8_t *s_enc = Impl::encode_nt4(subject, impl_->s_enc_buf);

	ksw_extz_t ez {};
	ksw_exts2_sse(nullptr, static_cast<int>(query.size()), q_enc, static_cast<int>(subject.size()), s_enc, 5,
	              impl_->mat, static_cast<int8_t>(impl_->gap_open), static_cast<int8_t>(impl_->gap_extend),
	              static_cast<int8_t>(impl_->gap_open2), static_cast<int8_t>(impl_->noncan), impl_->zdrop,
	              0 /* end_bonus */, 0 /* junc_bonus */, 0 /* junc_pen */, KSW_EZ_SPLICE_FOR, nullptr /* junc */, &ez);

	std::optional<int> result;
	if (ez.score != KSW_NEG_INF) {
		result = ez.score;
	}
	kfree(nullptr, ez.cigar);
	return result;
}

std::optional<KSW2CigarResult> KSW2Aligner::align_exts_cigar(const std::string &query, const std::string &subject) {
	if (query.empty() && subject.empty()) {
		return KSW2CigarResult {0, ""};
	}
	if (query.empty()) {
		return KSW2CigarResult {all_gap_score(impl_->gap_open, impl_->gap_extend, subject.size()),
		                        std::to_string(subject.size()) + "D"};
	}
	if (subject.empty()) {
		return KSW2CigarResult {all_gap_score(impl_->gap_open, impl_->gap_extend, query.size()),
		                        std::to_string(query.size()) + "I"};
	}

	uint8_t *q_enc = Impl::encode_nt4(query, impl_->q_enc_buf);
	uint8_t *s_enc = Impl::encode_nt4(subject, impl_->s_enc_buf);

	ksw_extz_t ez {};
	ksw_exts2_sse(nullptr, static_cast<int>(query.size()), q_enc, static_cast<int>(subject.size()), s_enc, 5,
	              impl_->mat, static_cast<int8_t>(impl_->gap_open), static_cast<int8_t>(impl_->gap_extend),
	              static_cast<int8_t>(impl_->gap_open2), static_cast<int8_t>(impl_->noncan), impl_->zdrop, 0, 0, 0,
	              KSW_EZ_SPLICE_FOR, nullptr, &ez);

	std::optional<KSW2CigarResult> result;
	if (ez.score != KSW_NEG_INF) {
		// eqx post-pass: split 'M' into '='/'X'. 'N' (intron skip) passes through unchanged
		// (see eqx_split_cigar in cigar_reconstruction.hpp).
		result = KSW2CigarResult {ez.score, eqx_split_cigar(decode_ksw2_cigar(ez.cigar, ez.n_cigar), query, subject)};
	}
	kfree(nullptr, ez.cigar);
	return result;
}

std::optional<KSW2FullResult> KSW2Aligner::align_exts_full(const std::string &query, const std::string &subject) {
	auto cigar_result = align_exts_cigar(query, subject);
	if (!cigar_result.has_value()) {
		return std::nullopt;
	}

	KSW2FullResult result;
	result.score = cigar_result->score;
	result.cigar = std::move(cigar_result->cigar);
	reconstruct_aligned_from_cigar(query, subject, result.cigar, result.query_aligned, result.subject_aligned);
	return result;
}

} // namespace miint
