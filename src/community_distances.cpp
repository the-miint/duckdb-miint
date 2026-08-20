#include "community_distances.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <system_error>
#include <thread>

namespace miint {

size_t CondensedRowBase(uint32_t n_samples, uint32_t i) {
	return static_cast<size_t>(static_cast<uint64_t>(i) * (2ull * n_samples - i - 1) / 2);
}

namespace {

// Accepted metric names, in a stable order that MUST match enum Metric below
// (the string is resolved to the enum once by index, before the O(n^2) loop).
constexpr std::array<const char *, 8> kMetrics = {"bray_curtis",   "euclidean", "jaccard", "soergel",
                                                  "morisita_horn", "pearson",   "chisq",   "gower"};

enum class Metric { BrayCurtis, Euclidean, Jaccard, Soergel, MorisitaHorn, Pearson, Chisq, Gower };

// Compile-time guard against the two lists drifting out of length: MetricIndex
// casts a kMetrics position straight to a Metric, so a future metric added to
// one list but not the other must fail to build, not silently misdispatch.
static_assert(kMetrics.size() == static_cast<size_t>(Metric::Gower) + 1,
              "kMetrics and enum class Metric must have the same number of entries");

// Whether each kMetrics entry is PAIRWISE-LOCAL, in the same order.
//
// Pairwise-local means d(i,j) is a function of samples i and j alone: it reads
// no statistic taken over the other rows, and it is unchanged by dropping
// feature columns that are zero in both i and j. Only such a metric can be
// computed one block of samples at a time (progressive_pcoa_from_features),
// because a block carries an arbitrary subset of the samples and only the
// features that subset happens to use.
//
// The three false entries each fail for a different reason, and all three fail
// SILENTLY — they return a plausible number that simply is not the same metric
// from block to block:
//   pearson  the per-sample mean is rowsum/f, and features zero in both samples
//            still contribute mi*mj to the covariance and mi^2/mj^2 to the
//            variances, so both the mean and the sums move with f
//   chisq    weights every term by grand_total/colsum[k]
//   gower    divides every term by colrange[k] = max_i M[i][k] - min_i M[i][k]
//
// Parallel-array rather than a field on the enum so that adding a metric to
// kMetrics without classifying it FAILS TO BUILD on the static_assert below,
// rather than defaulting to some answer nobody chose.
constexpr std::array<bool, 8> kPairwiseLocal = {
    true,  // bray_curtis    Sum|x-y| / Sum(x+y)          -- per-pair sums only
    true,  // euclidean      sqrt(Sum(x-y)^2)             -- per-pair sums only
    true,  // jaccard        (b+c)/(a+b+c)                -- per-pair counts only
    true,  // soergel        Sum|x-y| / Sum max(x,y)      -- per-pair sums only
    true,  // morisita_horn  normalizes by PER-SAMPLE totals, not column ones
    false, // pearson        feature-space dependent (see above)
    false, // chisq          needs column sums + grand total
    false, // gower          needs per-feature ranges
};
static_assert(kPairwiseLocal.size() == kMetrics.size(),
              "kPairwiseLocal must classify every kMetrics entry -- a new metric has to be "
              "declared pairwise-local or not before it can be used anywhere");

// Index of `metric` in kMetrics, or -1 if unknown.
int MetricIndex(const std::string &metric) {
	for (size_t i = 0; i < kMetrics.size(); ++i) {
		if (kMetrics[i] == metric) {
			return static_cast<int>(i);
		}
	}
	return -1;
}

// Number of unordered pairs for n samples, guarded against overflow of the
// reservation on tiny 32-bit builds (n is uint32; n*(n-1) fits in uint64).
size_t PairCount(uint32_t n) {
	return static_cast<size_t>(static_cast<uint64_t>(n) * static_cast<uint64_t>(n - 1) / 2);
}

// Base offset of row i's block in the condensed upper-triangle output: rows 0..i-1
// contribute sum_{t<i}(n-1-t) pairs. Every pair (i,j) therefore has a FIXED
// destination slot, independent of which thread computes it — which is what makes
// both pair loops bit-identical for any thread count.
//
// i*(2n-i-1) is always even (i and 2n-i-1 have opposite parity), so the halving is
// exact, and the whole product is formed in uint64 so it cannot wrap.
size_t RowBase(uint32_t n, uint32_t i) {
	return CondensedRowBase(n, i);
}

// Comma-separated metric names, either all of them or only the pairwise-local ones.
std::string JoinMetricNames(bool pairwise_local_only) {
	std::string out;
	for (size_t i = 0; i < kMetrics.size(); ++i) {
		if (pairwise_local_only && !kPairwiseLocal[i]) {
			continue;
		}
		if (!out.empty()) {
			out += ", ";
		}
		out += kMetrics[i];
	}
	return out;
}

// Validate the sample count and resolve a metric name, optionally requiring that it
// be pairwise-local. Shared so the dense and sparse entry points cannot drift on the
// message a caller sees for the same bad input.
Metric ResolveMetricOrThrow(uint32_t n_samples, const std::string &metric, bool require_pairwise_local) {
	if (n_samples < 2) {
		throw std::invalid_argument("community_distances requires at least 2 samples (got " +
		                            std::to_string(n_samples) + ")");
	}
	const int i = MetricIndex(metric);
	if (i < 0) {
		throw std::invalid_argument("community_distances: unknown metric '" + metric + "' (must be one of " +
		                            JoinMetricNames(/*pairwise_local_only=*/false) + ")");
	}
	if (require_pairwise_local && !kPairwiseLocal[static_cast<size_t>(i)]) {
		throw std::invalid_argument("community_distances: metric '" + metric +
		                            "' is not pairwise-local and cannot be computed from a single block of samples; "
		                            "it depends on statistics over the whole matrix. Pairwise-local metrics: " +
		                            JoinMetricNames(/*pairwise_local_only=*/true));
	}
	return static_cast<Metric>(i);
}

// Worker fan-out shared by both condensed pair loops. `compute_row(i)` must fill
// only row i's slots; every pair has a fixed destination (see RowBase), so which
// thread computes it cannot change a bit of the result and `n_threads` is purely a
// performance knob.
//
// Row i does (n-1-i) pairs, so a static contiguous split would overload the
// low-index threads; an atomic row cursor keeps the load balanced instead. The
// thread count is capped two ways — never more than there are rows to compute,
// and never more than the hardware can run concurrently, since extra threads
// on a CPU-bound loop only add context switches and spawning thousands risks a
// pids/ulimit hit. hardware_concurrency() reports 0 when it cannot detect; fall back
// to the requested count, leaving only the row cap.
//
// `row_end` is the number of rows to visit, normally n-1 (every row that has a
// pair). A caller that already holds the mutual distances of a trailing block of
// samples passes a smaller bound: because j > i, the pairs entirely inside that
// trailing block are exactly the rows at or after it, so skipping them is a loop
// bound rather than a per-pair test. Which thread reaches a row still cannot
// change a value — every pair has a fixed destination (see RowBase).
void RunPairLoop(uint32_t row_end, unsigned n_threads, const std::function<void(uint32_t)> &compute_row) {
	if (row_end == 0) {
		return;
	}
	if (n_threads <= 1) {
		for (uint32_t i = 0; i < row_end; ++i) {
			compute_row(i);
		}
		return;
	}
	unsigned hw = std::thread::hardware_concurrency();
	if (hw == 0) {
		hw = n_threads;
	}
	unsigned nt = std::min<unsigned>(n_threads, row_end);
	nt = std::min(nt, hw);
	std::atomic<uint32_t> next_row {0};
	auto worker = [&]() {
		for (;;) {
			const uint32_t i = next_row.fetch_add(1);
			if (i >= row_end) {
				break;
			}
			compute_row(i);
		}
	};
	std::vector<std::thread> pool;
	pool.reserve(nt - 1);
	try {
		for (unsigned t = 1; t < nt; ++t) {
			pool.emplace_back(worker);
		}
	} catch (const std::system_error &) {
		// A thread constructor failed (e.g. EAGAIN under a pids/ulimit cap). Do NOT
		// let it unwind past the joinable threads already in `pool` — a joinable
		// std::thread's destructor calls std::terminate() and would take the whole
		// process down. Swallow it and degrade: the threads already spawned keep
		// draining rows through the atomic cursor and the calling thread finishes the
		// remainder below, so the full result is still computed, just with fewer
		// workers. (`pool` is pre-reserved, so emplace_back itself cannot throw.)
	}
	worker();
	for (auto &th : pool) {
		th.join();
	}
}

// ── Dense Gram (GEMM) fast path ──────────────────────────────────────────────
//
// Three of the five pairwise-local metrics are a function of ONE inner product
// plus per-sample scalars, so all n(n-1)/2 of them can come out of a single
// matrix product instead of n(n-1)/2 independent merges:
//
//   euclidean      d^2   = Sx2_i + Sx2_j - 2*<x_i,x_j>
//   jaccard        a     = <1_i, 1_j> over the binary indicator, and
//                  d     = (|X|+|Y|-2a) / (|X|+|Y|-a)
//   morisita_horn  the numerator IS <x_i,x_j>
//
// bray_curtis and soergel are NOT expressible this way -- Sum|x-y| and
// Sum max(x,y) are not bilinear -- so they always take the merge. That is worth
// knowing rather than rediscovering: bray_curtis is the most-used microbiome
// metric and gets nothing from this path.
//
// WHY GATE IT: the merge visits the union of two samples' nonzeros; the GEMM
// visits the whole feature space. On real microbiome blocks that is a disaster,
// not a trade -- measured on a 2000-sample block, sub10k spans 56,142 features at
// 0.19% density, where the GEMM would do 4.5e11 flops against a merge that measures
// 0.8 s, and would need 898 MB for the dense operand against 2.6 MB of CSR. On a
// 38.8%-dense image block (673 features) the same comparison inverts and the GEMM is
// ~9.7x cheaper. So the dispatch is on DENSITY -- the work ratio is m*f/nnz =
// 1/density, with no n and no f in it -- and the crossover is MEASURED rather than
// argued from those two endpoints; see kGramDensityThreshold below for the ladder it
// was fitted on and the 0.96% it came out at.
//
// Eigen's GEMM is deterministic here because EIGEN_DONT_PARALLELIZE is set
// globally (see CMakeLists.txt) -- so this is bit-reproducible run to run for a
// given binary, and independent of `n_threads`, which only fans out the fill loop
// below. Across MACHINES it is not guaranteed: Eigen picks its blocking from
// runtime-queried cache sizes, so a different cache hierarchy can reassociate the
// sums. No test asserts an absolute distance to the bit (they use Catch2 Approx),
// and the tests that DO assert bits compare two paths inside one binary.
bool GemmEligible(Metric m) {
	return m == Metric::Euclidean || m == Metric::Jaccard || m == Metric::MorisitaHorn;
}

// Dispatch threshold, as a fraction of the block's cells that are nonzero.
//
// Fitted on a synthetic ladder that holds n, f, cluster structure and value
// distribution fixed and moves ONLY density, then validated on tables it was not
// fitted on -- so it is not tuned to whichever dataset was convenient. Measured at
// 20,000 samples over 2,000 features, one thread, as (merge wall) / (Gram wall):
//
//     density   0.8%   1.0%   1.2%   1.5%   2%    5%    15%    40%
//     speedup   0.83x  1.04x  1.25x  1.56x  2.1x  5.1x  13.0x  25.3x
//
// so the two kernels cost the same at about 0.96% and the Gram path is only worth
// taking some way above that. This is set at 1.5%, the lowest rung where the win is
// unambiguous rather than inside the noise, which leaves a 1.6x margin on the
// crossover -- because near it the trade is all cost and no benefit: equal time for
// a dense n x f operand the merge never allocates.
//
// Density is the whole rule because the work ratio is m*f/nnz = 1/density, with no n
// and no f in it. That is a theoretical claim, so it is measured rather than assumed:
// the same threshold was checked against 673 features (EMNIST, 9.7x at 38.8%) and
// 56,142 (a real microbiome table) and the crossover does not move with f.
constexpr double kGramDensityThreshold = 0.015;
// A hard ceiling on the dense OPERAND regardless of density: an optimization must
// never be the reason a query runs out of memory. It does not cover the n x n Gram,
// which is deliberately uncapped -- see CommunityDistancesUsesGramPath for why (it is
// 2x the result the caller already asked for, so it is a bounded factor, not a
// separate risk).
constexpr size_t kGramMaxOperandBytes = 512ull << 20;
// Below this the GEMM cannot pay for its own allocation, and staying on the merge
// keeps every small test fixture on the exactly-summed path.
constexpr uint32_t kGramMinSamples = 64;

bool UseGramPath(Metric m, uint32_t n, uint32_t f, size_t nnz, size_t max_operand_bytes) {
	if (!GemmEligible(m) || n < kGramMinSamples || f == 0) {
		return false;
	}
	if (max_operand_bytes == 0) {
		max_operand_bytes = kGramMaxOperandBytes;
	}
	// uint64 deliberately, not size_t: size_t is 32 bits on wasm32, which is a
	// first-class target here, and n * f overflows it for a whole-table feature
	// dictionary (millions of features). A wrapped product could land UNDER the cap and
	// then be used to size the dense operand — admitting the path and undersizing the
	// buffer at once. Same reason PairCount/RowBase promote to uint64 above.
	const uint64_t cells = static_cast<uint64_t>(n) * f;
	if (cells > max_operand_bytes / sizeof(double) || cells > std::numeric_limits<size_t>::max() / sizeof(double)) {
		return false;
	}
	return static_cast<double>(nnz) > kGramDensityThreshold * static_cast<double>(cells);
}

// `Sx2_i + Sx2_j - 2<x_i,x_j>` is a difference of large, nearly equal numbers when
// two rows are close, so it keeps almost none of its significant digits: the
// absolute error is about eps * (Sx2_i + Sx2_j) regardless of how small the true
// distance is. EXACT duplicates are safe by construction -- all three terms are the
// same sum over the same values, so they cancel to exactly 0 -- but NEAR duplicates
// with large magnitudes are destroyed. Measured on rows of ~1e6 differing in one
// feature by 1e-3: the true distance is 1e-3 and the product gives ~0.3.
//
// So any pair whose squared distance falls below this fraction of the row norms is
// recomputed exactly. At 1e-8 the residual relative error in d^2 is about
// eps/1e-8 ~ 2e-8, which is far below the fp32 narrowing every block gets anyway,
// and the guard only fires when d/||x|| < 1e-4 -- rare enough on real data to cost
// nothing, and exactly the regime where it must fire.
constexpr double kGramCancelRel = 1e-8;

// All pairwise distances for a GEMM-eligible metric, from a dense row-major
// operand. `n_cached_tail` has the same meaning as in the sparse entry point, and
// is honoured in the PRODUCT as well as the fill loop: the tail x tail triangle is
// never formed, which is the whole point of caching it.
std::vector<double> CondensedViaGram(const double *M, uint32_t n, uint32_t f, Metric m, unsigned n_threads,
                                     uint32_t n_cached_tail, bool operand_prebinarized = false) {
	using RowMajor = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

	// Jaccard's operand is the binary indicator, and presence is `> 0.0` -- the same
	// test the merge uses, so a negative cell counts as ABSENT in both. A caller that
	// materialized the operand itself can binarize while doing so and say so, which
	// saves holding a SECOND n x f matrix; the jaccard branch below reads only the Gram
	// and its diagonal, never the raw values, so a binarized operand loses it nothing.
	RowMajor binary;
	const double *operand = M;
	if (m == Metric::Jaccard && !operand_prebinarized) {
		binary.resize(n, f);
		for (uint32_t i = 0; i < n; ++i) {
			const double *xi = M + static_cast<size_t>(i) * f;
			for (uint32_t k = 0; k < f; ++k) {
				binary(i, k) = xi[k] > 0.0 ? 1.0 : 0.0;
			}
		}
		operand = binary.data();
	}
	const Eigen::Map<const RowMajor> X(operand, n, f);

	// Only the lower triangle is ever read (pair (i,j) with j > i reads S(j,i)), so
	// only the lower triangle is formed. With a cached tail the product splits into
	// the head triangle and the tail x head rectangle, leaving the tail's own
	// triangle uncomputed.
	Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
	const uint32_t head = n - n_cached_tail;
	if (n_cached_tail == 0) {
		S.selfadjointView<Eigen::Lower>().rankUpdate(X);
	} else if (head > 0) {
		// Reached only when n_cached_tail != 0 (the branch above takes 0), so the tail
		// rectangle is mandatory here, not optional. `head > 0` is defensive: the sparse
		// entry point returns early when the whole block is cached, so a zero head does
		// not arrive from there.
		S.topLeftCorner(head, head).selfadjointView<Eigen::Lower>().rankUpdate(X.topRows(head));
		S.bottomLeftCorner(n_cached_tail, head).noalias() =
		    X.bottomRows(n_cached_tail) * X.topRows(head).transpose();
	}

	// The Gram diagonal already holds Sum_k x_k^2 for euclidean and the presence
	// count for jaccard, so neither needs its own pass. Only morisita's Sum_k x_k
	// does -- and only when the tail was not cached, since a cached tail's own
	// diagonal was never formed.
	std::vector<double> diag(n, 0.0);
	for (uint32_t i = 0; i < head; ++i) {
		diag[i] = S(i, i);
	}
	for (uint32_t i = head; i < n; ++i) {
		const double *xi = operand + static_cast<size_t>(i) * f;
		double ss = 0.0;
		for (uint32_t k = 0; k < f; ++k) {
			ss += xi[k] * xi[k];
		}
		diag[i] = ss;
	}
	std::vector<double> rowsum;
	if (m == Metric::MorisitaHorn) {
		rowsum.assign(n, 0.0);
		for (uint32_t i = 0; i < n; ++i) {
			const double *xi = M + static_cast<size_t>(i) * f;
			double sv = 0.0;
			for (uint32_t k = 0; k < f; ++k) {
				sv += xi[k];
			}
			rowsum[i] = sv;
		}
	}

	std::vector<double> out(PairCount(n));
	auto compute_row = [&](uint32_t i) {
		size_t o = RowBase(n, i);
		const double *xi = M + static_cast<size_t>(i) * f;
		for (uint32_t j = i + 1; j < n; ++j, ++o) {
			const double dot = S(j, i);
			double dist = 0.0;
			switch (m) {
			case Metric::Euclidean: {
				double sq = diag[i] + diag[j] - 2.0 * dot;
				if (sq <= kGramCancelRel * (diag[i] + diag[j])) {
					// Cancellation territory -- recompute this pair exactly, in the same
					// ascending feature order the merge would use.
					const double *xj = M + static_cast<size_t>(j) * f;
					sq = 0.0;
					for (uint32_t k = 0; k < f; ++k) {
						const double dk = xi[k] - xj[k];
						sq += dk * dk;
					}
				}
				dist = sq > 0.0 ? std::sqrt(sq) : 0.0;
				break;
			}
			case Metric::Jaccard: {
				// diag holds |X| and |Y| (presence counts) off the binary Gram.
				const double union_sz = diag[i] + diag[j] - dot;
				dist = union_sz > 0.0 ? (union_sz - dot) / union_sz : 0.0;
				break;
			}
			case Metric::MorisitaHorn: {
				const double X_ = rowsum[i], Y_ = rowsum[j];
				if (X_ <= 0.0 && Y_ <= 0.0) {
					dist = 0.0;
				} else if (X_ <= 0.0 || Y_ <= 0.0) {
					dist = 1.0;
				} else {
					// No denominator guard, because the merge branch has none either: the
					// two paths have to agree on every input, including the ones where a
					// negative abundance makes this ill-conditioned.
					dist = 1.0 - 2 * dot / ((diag[i] / (X_ * X_) + diag[j] / (Y_ * Y_)) * X_ * Y_);
				}
				break;
			}
			case Metric::BrayCurtis:
			case Metric::Soergel:
			case Metric::Pearson:
			case Metric::Chisq:
			case Metric::Gower:
				// Unreachable: UseGramPath admits only the three above. Listed rather
				// than defaulted so a metric newly shown to be Gram-expressible has to
				// be handled here instead of silently yielding zeros.
				break;
			}
			out[o] = dist;
		}
	};
	RunPairLoop(n - std::max<uint32_t>(n_cached_tail, 1), n_threads, compute_row);
	return out;
}

} // namespace

bool IsValidCommunityMetric(const std::string &metric) {
	return MetricIndex(metric) >= 0;
}

std::string PairwiseLocalCommunityMetricList() {
	return JoinMetricNames(/*pairwise_local_only=*/true);
}

bool IsPairwiseLocalCommunityMetric(const std::string &metric) {
	const int i = MetricIndex(metric);
	// An unknown name is not admissible: refusing it here is what makes a typo
	// an error at bind rather than a metric silently computed per block.
	return i >= 0 && kPairwiseLocal[static_cast<size_t>(i)];
}

bool CommunityDistancesUsesGramPath(const std::string &metric, uint32_t n_samples, uint32_t n_features, size_t nnz,
                                    size_t max_operand_bytes) {
	const int i = MetricIndex(metric);
	if (i < 0) {
		return false;
	}
	return UseGramPath(static_cast<Metric>(i), n_samples, n_features, nnz, max_operand_bytes);
}

std::string CommunityMetricList() {
	return JoinMetricNames(/*pairwise_local_only=*/false);
}

std::vector<double> CommunityDistancesCondensed(const std::vector<double> &matrix, uint32_t n_samples,
                                                uint32_t n_features, const std::string &metric, unsigned n_threads) {
	const Metric m = ResolveMetricOrThrow(n_samples, metric, /*require_pairwise_local=*/false);
	if (matrix.size() != static_cast<size_t>(n_samples) * static_cast<size_t>(n_features)) {
		throw std::invalid_argument("community_distances: matrix size (" + std::to_string(matrix.size()) +
		                            ") does not match n_samples*n_features (" + std::to_string(n_samples) + "*" +
		                            std::to_string(n_features) + ")");
	}

	const uint32_t n = n_samples;
	const uint32_t f = n_features;
	const double *M = matrix.data();
	auto row = [&](uint32_t i) {
		return M + static_cast<size_t>(i) * f;
	};

	// One matrix product instead of n(n-1)/2 pair scans, when the block is dense
	// enough for that to be cheaper (see UseGramPath). Dispatched on the SAME
	// (metric, n, f, nnz) rule as the sparse entry point below and on the same dense
	// values, which is what keeps the two bit-identical on this path as they are on
	// the merge path. The nonzero count costs one O(n*f) pass, the same order as the
	// per-sample pre-pass just below it.
	if (GemmEligible(m)) {
		size_t nnz = 0;
		for (size_t c = 0; c < matrix.size(); ++c) {
			nnz += (matrix[c] != 0.0) ? 1 : 0;
		}
		// This entry point does not materialize an operand -- `matrix` is the caller's --
		// so it spends only the library default, never a caller's per-block budget.
		if (UseGramPath(m, n, f, nnz, /*max_operand_bytes=*/0)) {
			return CondensedViaGram(M, n, f, m, n_threads, /*n_cached_tail=*/0);
		}
	}

	// Per-sample aggregates reused across pairs.
	std::vector<double> rowsum(n, 0.0);   // Sum_k M[i][k]     (X, Y; bray denom; morisita/chisq)
	std::vector<double> rowsumsq(n, 0.0); // Sum_k M[i][k]^2   (morisita lambda)
	for (uint32_t i = 0; i < n; ++i) {
		const double *xi = row(i);
		double s = 0.0, ss = 0.0;
		for (uint32_t k = 0; k < f; ++k) {
			s += xi[k];
			ss += xi[k] * xi[k];
		}
		rowsum[i] = s;
		rowsumsq[i] = ss;
	}

	// Global column statistics for the correspondence-analysis chi-square and the
	// Gower distance (both are matrix-wide, not purely per-pair).
	std::vector<double> colsum;
	std::vector<double> colrange;
	double grand = 0.0;
	if (m == Metric::Chisq) {
		colsum.assign(f, 0.0);
		for (uint32_t i = 0; i < n; ++i) {
			const double *xi = row(i);
			for (uint32_t k = 0; k < f; ++k) {
				colsum[k] += xi[k];
			}
		}
		for (double c : colsum) {
			grand += c;
		}
	} else if (m == Metric::Gower) {
		colrange.assign(f, 0.0);
		for (uint32_t k = 0; k < f; ++k) {
			double lo = row(0)[k];
			double hi = lo;
			for (uint32_t i = 1; i < n; ++i) {
				const double v = row(i)[k];
				lo = std::min(lo, v);
				hi = std::max(hi, v);
			}
			colrange[k] = hi - lo;
		}
	}

	std::vector<double> out(PairCount(n));

	// All pairs (i, j>i) of one outer row, each written to its fixed slot.
	auto compute_row = [&](uint32_t i) {
		const double *xi = row(i);
		size_t o = RowBase(n, i);
		for (uint32_t j = i + 1; j < n; ++j) {
			const double *yj = row(j);
			double d = 0.0;

			switch (m) {
			case Metric::BrayCurtis: {
				double num = 0.0;
				for (uint32_t k = 0; k < f; ++k) {
					num += std::fabs(xi[k] - yj[k]);
				}
				const double den = rowsum[i] + rowsum[j];
				d = den > 0.0 ? num / den : 0.0;
				break;
			}
			case Metric::Euclidean: {
				double s = 0.0;
				for (uint32_t k = 0; k < f; ++k) {
					const double diff = xi[k] - yj[k];
					s += diff * diff;
				}
				d = std::sqrt(s);
				break;
			}
			case Metric::Jaccard: {
				// Binary presence/absence: (b+c)/(a+b+c).
				double bc = 0.0, abc = 0.0;
				for (uint32_t k = 0; k < f; ++k) {
					const bool px = xi[k] > 0.0;
					const bool py = yj[k] > 0.0;
					if (px || py) {
						abc += 1.0;
						if (px != py) {
							bc += 1.0;
						}
					}
				}
				d = abc > 0.0 ? bc / abc : 0.0;
				break;
			}
			case Metric::Soergel: {
				double num = 0.0, den = 0.0;
				for (uint32_t k = 0; k < f; ++k) {
					num += std::fabs(xi[k] - yj[k]);
					den += std::max(xi[k], yj[k]);
				}
				d = den > 0.0 ? num / den : 0.0;
				break;
			}
			case Metric::MorisitaHorn: {
				const double X = rowsum[i];
				const double Y = rowsum[j];
				// Empty-community handling. Reachable only from a dense caller: the
				// SQL wrapper reads a sparse feature table that drops zero cells, so
				// an all-zero sample never reaches here (it is absent from the
				// dictionary — see community_distances_function.cpp). Kept as defensive,
				// well-defined behavior for any future in-process dense caller.
				if (X <= 0.0 && Y <= 0.0) {
					d = 0.0; // two empty communities are identical
				} else if (X <= 0.0 || Y <= 0.0) {
					d = 1.0; // one empty: no shared abundance
				} else {
					double dot = 0.0;
					for (uint32_t k = 0; k < f; ++k) {
						dot += xi[k] * yj[k];
					}
					const double lx = rowsumsq[i] / (X * X);
					const double ly = rowsumsq[j] / (Y * Y);
					const double c_h = 2.0 * dot / ((lx + ly) * X * Y);
					d = 1.0 - c_h;
				}
				break;
			}
			case Metric::Pearson: {
				// 1 - Pearson correlation over features. Centered (two-pass) rather
				// than the Sum(x^2) - f*mean^2 shortcut: the shortcut suffers
				// catastrophic cancellation on high-mean/low-variance profiles
				// (common with a few dominant high-count features), which can both
				// diverge from the reference and spuriously report zero variance for
				// a genuinely non-constant row. De-means before the dot products.
				const double mi = rowsum[i] / static_cast<double>(f);
				const double mj = rowsum[j] / static_cast<double>(f);
				double cov = 0.0, varx = 0.0, vary = 0.0;
				for (uint32_t k = 0; k < f; ++k) {
					const double dxi = xi[k] - mi;
					const double dyj = yj[k] - mj;
					cov += dxi * dyj;
					varx += dxi * dxi;
					vary += dyj * dyj;
				}
				// Flat (constant) profiles have zero variance -> correlation is
				// undefined. We adopt the cogent3 `dist_pearson` convention (the
				// reference implementation of the metric Kuczynski 2010 used), NOT
				// scipy's: two flat rows are identical (r=1 -> distance 0); a flat vs
				// a non-flat row has no correlation (r=0 -> distance 1). This keeps a
				// constant-profile sample as a well-defined, finite distance rather
				// than a NaN that a downstream reader would reject.
				if (varx == 0.0 && vary == 0.0) {
					d = 0.0;
				} else if (varx == 0.0 || vary == 0.0) {
					d = 1.0;
				} else {
					d = 1.0 - cov / std::sqrt(varx * vary);
				}
				break;
			}
			case Metric::Chisq: {
				const double ri = rowsum[i];
				const double rj = rowsum[j];
				if (ri <= 0.0 && rj <= 0.0) {
					// Both empty: no row profiles exist, but the cogent3 `dist_chisq`
					// convention defines this as distance 0 (identical). Adopted for
					// faithful reproduction.
					d = 0.0;
				} else if (ri <= 0.0 || rj <= 0.0) {
					d = 1.0; // one empty vs non-empty -> maximal (cogent3 dist_chisq)
				} else {
					double s = 0.0;
					for (uint32_t k = 0; k < f; ++k) {
						if (colsum[k] <= 0.0) {
							continue; // globally absent feature contributes nothing
						}
						const double a = xi[k] / ri;
						const double b = yj[k] / rj;
						const double diff = a - b;
						s += (grand / colsum[k]) * diff * diff;
					}
					d = std::sqrt(s);
				}
				break;
			}
			case Metric::Gower: {
				double s = 0.0;
				for (uint32_t k = 0; k < f; ++k) {
					if (colrange[k] <= 0.0) {
						continue; // constant feature across all samples
					}
					s += std::fabs(xi[k] - yj[k]) / colrange[k];
				}
				d = s;
				break;
			}
			}

			out[o++] = d;
		}
	};

	RunPairLoop(n - 1, n_threads, compute_row);
	return out;
}

std::vector<double> CommunityDistancesCondensedSparse(const std::vector<uint32_t> &indptr,
                                                      const std::vector<uint32_t> &indices,
                                                      const std::vector<double> &values, uint32_t n_samples,
                                                      uint32_t n_features, const std::string &metric,
                                                      unsigned n_threads, uint32_t n_cached_tail,
                                                      size_t max_operand_bytes) {
	const Metric m = ResolveMetricOrThrow(n_samples, metric, /*require_pairwise_local=*/true);

	if (n_cached_tail > n_samples) {
		throw std::invalid_argument("community_distances: n_cached_tail (" + std::to_string(n_cached_tail) +
		                            ") exceeds n_samples (" + std::to_string(n_samples) + ")");
	}

	if (indptr.size() != static_cast<size_t>(n_samples) + 1) {
		throw std::invalid_argument("community_distances: indptr must have n_samples+1 (" +
		                            std::to_string(n_samples + 1) + ") entries, got " + std::to_string(indptr.size()));
	}
	if (indices.size() != values.size()) {
		throw std::invalid_argument("community_distances: indices (" + std::to_string(indices.size()) +
		                            ") and values (" + std::to_string(values.size()) + ") must be the same length");
	}
	if (indptr.front() != 0 || indptr.back() != indices.size()) {
		throw std::invalid_argument("community_distances: indptr must start at 0 and end at indices.size()");
	}
	for (uint32_t i = 0; i < n_samples; ++i) {
		if (indptr[i] > indptr[i + 1]) {
			throw std::invalid_argument("community_distances: indptr must be non-decreasing (row " + std::to_string(i) +
			                            ")");
		}
		for (uint32_t p = indptr[i]; p < indptr[i + 1]; ++p) {
			if (indices[p] >= n_features) {
				throw std::invalid_argument("community_distances: feature index " + std::to_string(indices[p]) +
				                            " is out of range for n_features=" + std::to_string(n_features));
			}
			if (p > indptr[i] && indices[p] <= indices[p - 1]) {
				throw std::invalid_argument("community_distances: feature indices must be strictly ascending within "
				                            "each row (row " +
				                            std::to_string(i) + ")");
			}
		}
	}

	const uint32_t n = n_samples;

	// Every pair is in the caller's cache, so there is nothing to compute. Worth an
	// early return rather than falling through: progressive_pcoa_from_features asks for
	// exactly this when it fetches the anchor-only reference block, and without it that
	// block densifies n*f doubles, forms an n x n Gram and walks every pre-pass -- and
	// then has all PairCount(n) of its results overwritten from the cache.
	if (n_cached_tail == n) {
		return std::vector<double>(PairCount(n), 0.0);
	}

	// Same dispatch as the dense entry point, on the same rule, so both agree.
	// Densifying costs n*f doubles -- which is exactly what UseGramPath's operand cap
	// bounds, and why that cap is applied on both sides even though the dense entry
	// point does not allocate.
	// VALUE nonzeros, not indices.size(). BuildBlockCsr coalesces duplicate cells by
	// SUMMING them, so a CSR can legitimately store an explicit 0.0 — and the dense
	// entry point counts `matrix[c] != 0.0`. Counting structural entries here would let
	// the same logical matrix take different kernels through the two entry points,
	// which is exactly the bit-identity the header promises.
	size_t nnz = 0;
	if (GemmEligible(m)) {
		for (size_t k = 0; k < values.size(); ++k) {
			nnz += (values[k] != 0.0) ? 1 : 0;
		}
	}
	if (UseGramPath(m, n, n_features, nnz, max_operand_bytes)) { // its own first test is !GemmEligible
		// Binarize WHILE densifying for jaccard, so the indicator does not cost a second
		// n x f matrix on top of this one -- which the operand cap does not budget for.
		const bool binarize = (m == Metric::Jaccard);
		std::vector<double> dense(static_cast<size_t>(n) * n_features, 0.0);
		for (uint32_t i = 0; i < n; ++i) {
			for (uint32_t pos = indptr[i]; pos < indptr[i + 1]; ++pos) {
				const double v = values[pos];
				dense[static_cast<size_t>(i) * n_features + indices[pos]] = binarize ? (v > 0.0 ? 1.0 : 0.0) : v;
			}
		}
		return CondensedViaGram(dense.data(), n, n_features, m, n_threads, n_cached_tail,
		                        /*operand_prebinarized=*/binarize);
	}

	// Per-sample aggregates, accumulated in ascending feature order exactly as the
	// dense pre-pass does (its extra zero terms add nothing).
	std::vector<double> rowsum(n, 0.0);
	std::vector<double> rowsumsq(n, 0.0);
	std::vector<uint32_t> present(n, 0); // nonzero count, for jaccard's presence sets
	for (uint32_t i = 0; i < n; ++i) {
		double s = 0.0, ss = 0.0;
		uint32_t np = 0;
		for (uint32_t p = indptr[i]; p < indptr[i + 1]; ++p) {
			const double v = values[p];
			s += v;
			ss += v * v;
			if (v > 0.0) {
				++np;
			}
		}
		rowsum[i] = s;
		rowsumsq[i] = ss;
		present[i] = np;
	}

	std::vector<double> out(PairCount(n));

	auto compute_row = [&](uint32_t i) {
		size_t o = RowBase(n, i);
		for (uint32_t j = i + 1; j < n; ++j, ++o) {
			// Every guard below compares against 0.0 with the SAME operator the dense
			// loop uses (`> 0.0` for a denominator or a presence test, `<= 0.0` for an
			// empty community) rather than testing equality. Nothing rejects negative
			// abundances on the way in — the feature-table filter drops only NULL, zero
			// and NaN — so `== 0.0` would let a negative denominator through and emit a
			// NEGATIVE distance, and would count a negative cell as "present". Both make
			// this function disagree with community_distances on the same table.
			//
			// Merge the two rows' nonzeros in ascending feature order. Every term the
			// dense loop would add for a feature absent from both is exactly 0.0, so
			// visiting only the union reproduces its sum bit for bit.
			uint32_t pi = indptr[i], pj = indptr[j];
			const uint32_t ei = indptr[i + 1], ej = indptr[j + 1];
			double abs_diff = 0.0; // Sum|x-y|
			double sq_diff = 0.0;  // Sum(x-y)^2
			double max_sum = 0.0;  // Sum max(x,y)
			double dot = 0.0;      // Sum x*y
			uint32_t shared = 0;   // features nonzero in BOTH (jaccard's a)
			while (pi < ei || pj < ej) {
				double x = 0.0, y = 0.0;
				if (pj >= ej || (pi < ei && indices[pi] < indices[pj])) {
					x = values[pi++];
				} else if (pi >= ei || indices[pj] < indices[pi]) {
					y = values[pj++];
				} else {
					x = values[pi++];
					y = values[pj++];
				}
				const double d = x - y;
				abs_diff += std::fabs(d);
				sq_diff += d * d;
				max_sum += std::max(x, y);
				dot += x * y;
				if (x > 0.0 && y > 0.0) {
					++shared;
				}
			}

			double dist = 0.0;
			switch (m) {
			case Metric::BrayCurtis: {
				const double den = rowsum[i] + rowsum[j];
				dist = den > 0.0 ? abs_diff / den : 0.0;
				break;
			}
			case Metric::Euclidean:
				dist = std::sqrt(sq_diff);
				break;
			case Metric::Jaccard: {
				// a = shared, b + c = the rest of the union.
				const uint32_t uni = present[i] + present[j] - shared;
				dist = uni == 0 ? 0.0 : static_cast<double>(uni - shared) / static_cast<double>(uni);
				break;
			}
			case Metric::Soergel:
				dist = max_sum > 0.0 ? abs_diff / max_sum : 0.0;
				break;
			case Metric::MorisitaHorn: {
				const double X = rowsum[i], Y = rowsum[j];
				if (X <= 0.0 && Y <= 0.0) {
					dist = 0.0;
				} else if (X <= 0.0 || Y <= 0.0) {
					dist = 1.0;
				} else {
					dist = 1.0 - 2 * dot / ((rowsumsq[i] / (X * X) + rowsumsq[j] / (Y * Y)) * X * Y);
				}
				break;
			}
			case Metric::Pearson:
			case Metric::Chisq:
			case Metric::Gower:
				// Unreachable: ResolveMetricOrThrow refused these above. Listed rather
				// than folded into a `default` so that a metric newly classified
				// pairwise-local has to be handled here instead of silently yielding a
				// block of zeros.
				break;
			}
			out[o] = dist;
		}
	};

	// Rows i >= n - n_cached_tail hold only pairs whose partner is also in the tail
	// (j > i), so stopping here is exactly the cached square and nothing else. The
	// per-sample pre-pass above deliberately still covers those samples — they are
	// still one endpoint of every cross pair. max(n_cached_tail, 1) keeps the
	// no-cache case at the usual n-1.
	RunPairLoop(n - std::max<uint32_t>(n_cached_tail, 1), n_threads, compute_row);
	return out;
}

} // namespace miint
