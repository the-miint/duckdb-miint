#include "MsaConsensus.hpp"

#include "alignment_functions_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>

namespace miint {

namespace {

// 5 candidate true states: A, C, G, T, '-'. Index order matters: ties in
// log-likelihood break on this order (alphabetical-then-gap), which matches
// the deterministic tie-break the design requires.
constexpr std::array<char, 5> kCandidates = {'A', 'C', 'G', 'T', '-'};

// Gap-observation error probability (the analogue of 10^(-q/10) for non-gap
// observations). MAFFT places gaps imperfectly, so a fixed 5% accounts for
// the fact that an observed gap is *some* evidence for a true gap but not
// overwhelming. Tuned to keep "1 base / 4 gaps" suppressing the column and
// "4 base / 1 gap" keeping the base — see test/cpp/test_consensus_msa.cpp.
constexpr double kGapErrProb = 0.05;

// Posterior Q ceiling. Spec'd in FINDINGS Part 3 (the "Q capped at 60" cap).
constexpr std::uint8_t kPosteriorQMax = 60;

// Index of `b` in kCandidates, or -1 if not a recognised base. Lower-case is
// up-cased by the caller before reaching us.
int CandidateIndex(char b) {
	for (int i = 0; i < 5; ++i) {
		if (kCandidates[i] == b) {
			return i;
		}
	}
	return -1;
}

// Per-observation log P(obs | true == candidate[k]).
double LogProb(const Observation &obs, int true_idx) {
	const double p_err = (obs.base == '-') ? kGapErrProb : std::pow(10.0, -static_cast<double>(obs.qual) / 10.0);
	const int obs_idx = CandidateIndex(obs.base);
	if (obs_idx < 0) {
		// Unknown observation base (e.g. N, IUPAC ambiguity). Treat as no
		// information: log(1/5) for every candidate. The aggregate uppercases
		// + maps N to '-' upstream, so this is a defensive fallback.
		return std::log(1.0 / 5.0);
	}
	if (obs_idx == true_idx) {
		// Clamp away from exactly 1.0 to keep log finite when q is so high
		// that p_err underflows to 0.
		const double p_correct = 1.0 - p_err;
		return std::log(std::max(p_correct, 1e-300));
	}
	// One of the 4 alternative candidates.
	return std::log(std::max(p_err / 4.0, 1e-300));
}

// Numerically-stable log(sum_i exp(v[i])) via the log-sum-exp identity. The
// span is templated so we can call it for both the full 5-candidate ll array
// and the 4-candidate "all but best" subset used for posterior Q.
template <std::size_t N>
double LogSumExp(const std::array<double, N> &v) {
	const double m = *std::max_element(v.begin(), v.end());
	double s = 0.0;
	for (double x : v) {
		s += std::exp(x - m);
	}
	return m + std::log(s);
}

} // namespace

ConsensusBase VoteColumn(const std::vector<Observation> &obs) {
	if (obs.empty()) {
		throw miint::InvalidInputException("VoteColumn: empty observation list");
	}

	// Per-candidate log-likelihood and Q-sum for tie-break.
	std::array<double, 5> ll {};
	std::array<int, 5> qsum {};
	for (auto &o : obs) {
		const int obs_idx = CandidateIndex(o.base);
		if (obs_idx >= 0) {
			qsum[obs_idx] += o.qual;
		}
		for (int k = 0; k < 5; ++k) {
			ll[k] += LogProb(o, k);
		}
	}

	// Argmax with deterministic tie-break: higher ll wins; on near-tie
	// (within 1e-9 in log-likelihood units, well below floating noise) the
	// higher Q-sum for observations of that candidate wins; final fallback
	// is the kCandidates declaration order (A < C < G < T < -).
	int best = 0;
	for (int k = 1; k < 5; ++k) {
		const double diff = ll[k] - ll[best];
		if (diff > 1e-9) {
			best = k;
		} else if (diff > -1e-9) {
			if (qsum[k] > qsum[best]) {
				best = k;
			}
		}
	}

	// Posterior Q. Compute log(err) = log(sum_{k!=best} exp(ll[k]) / sum_all)
	// directly, NOT via log(1 - posterior), because at high confidence the
	// 1.0 - exp(log_post) cancellation loses all precision.
	const double log_z = LogSumExp(ll);
	std::array<double, 4> ll_others {};
	std::size_t oi = 0;
	for (int k = 0; k < 5; ++k) {
		if (k != best) {
			ll_others[oi++] = ll[k];
		}
	}
	const double log_err = LogSumExp(ll_others) - log_z;
	// Convert nats → Phred:  Q = -10 * log10(err) = (-10 / ln 10) * log_err
	double q = -10.0 / std::log(10.0) * log_err;
	if (!std::isfinite(q) || q > static_cast<double>(kPosteriorQMax)) {
		q = static_cast<double>(kPosteriorQMax);
	}
	if (q < 0.0) {
		q = 0.0;
	}

	ConsensusBase out;
	out.base = kCandidates[best];
	out.qual = static_cast<std::uint8_t>(std::round(q));
	return out;
}

namespace {

// Per-row precomputation: a prefix-count of non-gap chars per MSA column, and
// the ungapped sequence itself. We compute these once per call to HpCorrect
// because the same rows are visited once per HP run.
struct RowIndex {
	std::vector<std::size_t> prefix_nongap; // prefix_nongap[j] = count of non-gap chars in row[0..j)
	std::string ungapped;
};

std::vector<RowIndex> BuildRowIndex(const std::vector<std::string> &aligned_seqs) {
	std::vector<RowIndex> idx;
	idx.reserve(aligned_seqs.size());
	for (const auto &row : aligned_seqs) {
		RowIndex ri;
		ri.prefix_nongap.assign(row.size() + 1, 0);
		ri.ungapped.reserve(row.size());
		for (std::size_t j = 0; j < row.size(); ++j) {
			ri.prefix_nongap[j + 1] = ri.prefix_nongap[j] + (row[j] == '-' ? 0 : 1);
			if (row[j] != '-') {
				ri.ungapped.push_back(row[j]);
			}
		}
		idx.push_back(std::move(ri));
	}
	return idx;
}

// For a row, find the per-row homopolymer length around the ungapped position
// derived from MSA column `m_start`. The qpos lands on the first non-gap char
// at or after m_start. Returns 0 if the row has no HP of `base` there.
//
// V1 limitation: a row whose MSA column `m_start` is a gap, but which has an
// HP run of `base` further downstream within the same MSA-level HP block, is
// dropped from the median (returns 0). For a row with a leading-side gap
// inside a staggered HP this under-samples the read's HP length. In Karst
// UMI amplicons with well-aligned 5–30-read bins this is a small effect
// (most reads have the HP base at the run's start column); the spec-grade
// fix would search for the nearest-downstream HP of `base` within a window,
// at a cost of more parameters. Deferred to v2 if real data shows bias.
std::size_t RowHpLength(const RowIndex &ri, std::size_t m_start, char base) {
	if (m_start >= ri.prefix_nongap.size()) {
		return 0;
	}
	const std::size_t qpos = ri.prefix_nongap[m_start];
	if (qpos >= ri.ungapped.size()) {
		return 0;
	}
	if (ri.ungapped[qpos] != base) {
		return 0;
	}
	std::size_t lo = qpos;
	while (lo > 0 && ri.ungapped[lo - 1] == base) {
		--lo;
	}
	std::size_t hi = qpos;
	while (hi + 1 < ri.ungapped.size() && ri.ungapped[hi + 1] == base) {
		++hi;
	}
	return hi - lo + 1;
}

} // namespace

namespace {

// Build per-column observation lists from the per-row MSA + qual data. Each
// non-gap aligned-seq position consumes one qual element from the row's
// ungapped qual list; gap positions emit ('-', 0).
std::vector<std::vector<Observation>> BuildColumns(const std::vector<std::string> &aligned_seqs,
                                                   const std::vector<std::vector<std::uint8_t>> &quals) {
	const std::size_t width = aligned_seqs.front().size();
	std::vector<std::vector<Observation>> cols(width);
	for (std::size_t i = 0; i < aligned_seqs.size(); ++i) {
		const auto &row = aligned_seqs[i];
		const auto &q = quals[i];
		std::size_t qpos = 0;
		for (std::size_t j = 0; j < width; ++j) {
			if (row[j] == '-') {
				cols[j].push_back({'-', 0});
			} else {
				cols[j].push_back({row[j], q[qpos]});
				++qpos;
			}
		}
	}
	return cols;
}

// FINDINGS spec: "Bin of size 1: skip MSA entirely, emit the lone read's seq
// and qual directly." Strips gap characters but does not touch qual.
std::pair<std::string, std::vector<std::uint8_t>> SingleReadPassthrough(const std::string &aligned,
                                                                        const std::vector<std::uint8_t> &qual) {
	std::string seq;
	seq.reserve(aligned.size());
	for (char c : aligned) {
		if (c != '-') {
			seq.push_back(c);
		}
	}
	return {seq, qual};
}

} // namespace

std::pair<std::string, std::vector<std::uint8_t>> BuildConsensus(const std::vector<std::string> &aligned_seqs,
                                                                 const std::vector<std::vector<std::uint8_t>> &quals) {
	if (aligned_seqs.empty()) {
		throw miint::InvalidInputException("BuildConsensus: no rows supplied");
	}
	if (aligned_seqs.size() != quals.size()) {
		throw miint::InvalidInputException("BuildConsensus: aligned_seqs / quals length mismatch");
	}
	const std::size_t width = aligned_seqs.front().size();
	for (const auto &row : aligned_seqs) {
		if (row.size() != width) {
			throw miint::InvalidInputException("BuildConsensus: inconsistent MSA row widths in group");
		}
	}
	if (aligned_seqs.size() == 1) {
		return SingleReadPassthrough(aligned_seqs.front(), quals.front());
	}

	const auto cols = BuildColumns(aligned_seqs, quals);

	std::string naive_seq;
	std::vector<std::uint8_t> naive_qual;
	std::vector<std::size_t> msa_col_of_consensus_pos;
	naive_seq.reserve(cols.size());
	naive_qual.reserve(cols.size());
	msa_col_of_consensus_pos.reserve(cols.size());

	for (std::size_t j = 0; j < cols.size(); ++j) {
		const auto v = VoteColumn(cols[j]);
		if (v.base == '-') {
			continue; // column suppressed
		}
		naive_seq.push_back(v.base);
		naive_qual.push_back(v.qual);
		msa_col_of_consensus_pos.push_back(j);
	}

	return HpCorrect(naive_seq, naive_qual, msa_col_of_consensus_pos, aligned_seqs);
}

std::pair<std::string, std::vector<std::uint8_t>> HpCorrect(const std::string &consensus_seq,
                                                            const std::vector<std::uint8_t> &consensus_qual,
                                                            const std::vector<std::size_t> &msa_col_of_consensus_pos,
                                                            const std::vector<std::string> &aligned_seqs) {
	if (consensus_seq.size() != consensus_qual.size() || consensus_seq.size() != msa_col_of_consensus_pos.size()) {
		throw miint::InvalidInputException("HpCorrect: consensus / qual / msa-col-map size mismatch");
	}
	const auto row_index = BuildRowIndex(aligned_seqs);

	std::string out_seq;
	std::vector<std::uint8_t> out_qual;
	out_seq.reserve(consensus_seq.size());
	out_qual.reserve(consensus_qual.size());

	std::size_t i = 0;
	while (i < consensus_seq.size()) {
		std::size_t j = i;
		while (j < consensus_seq.size() && consensus_seq[j] == consensus_seq[i]) {
			++j;
		}
		const char base = consensus_seq[i];
		const std::size_t run_len = j - i;

		std::size_t emit_len = run_len;
		if (run_len >= 2 && !row_index.empty()) {
			const std::size_t m_start = msa_col_of_consensus_pos[i];
			std::vector<std::size_t> hp_lengths;
			hp_lengths.reserve(row_index.size());
			for (const auto &ri : row_index) {
				const std::size_t len = RowHpLength(ri, m_start, base);
				if (len > 0) {
					hp_lengths.push_back(len);
				}
			}
			if (!hp_lengths.empty()) {
				// Lower median: for even count, take the lower of the two
				// central values. Deterministic; avoids producing a length
				// that no read actually observed.
				std::sort(hp_lengths.begin(), hp_lengths.end());
				emit_len = hp_lengths[(hp_lengths.size() - 1) / 2];
			}
		}

		const std::uint8_t q = consensus_qual[i];
		for (std::size_t k = 0; k < emit_len; ++k) {
			out_seq.push_back(base);
			out_qual.push_back(q);
		}
		i = j;
	}

	return {out_seq, out_qual};
}

} // namespace miint
