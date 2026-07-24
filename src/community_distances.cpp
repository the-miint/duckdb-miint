#include "community_distances.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <system_error>
#include <thread>

namespace miint {

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

} // namespace

bool IsValidCommunityMetric(const std::string &metric) {
	return MetricIndex(metric) >= 0;
}

std::string CommunityMetricList() {
	std::string out;
	for (size_t i = 0; i < kMetrics.size(); ++i) {
		if (i != 0) {
			out += ", ";
		}
		out += kMetrics[i];
	}
	return out;
}

std::vector<double> CommunityDistancesCondensed(const std::vector<double> &matrix, uint32_t n_samples,
                                                uint32_t n_features, const std::string &metric, unsigned n_threads) {
	if (n_samples < 2) {
		throw std::invalid_argument("community_distances requires at least 2 samples (got " +
		                            std::to_string(n_samples) + ")");
	}
	if (matrix.size() != static_cast<size_t>(n_samples) * static_cast<size_t>(n_features)) {
		throw std::invalid_argument("community_distances: matrix size (" + std::to_string(matrix.size()) +
		                            ") does not match n_samples*n_features (" + std::to_string(n_samples) + "*" +
		                            std::to_string(n_features) + ")");
	}
	const int metric_index = MetricIndex(metric);
	if (metric_index < 0) {
		throw std::invalid_argument("community_distances: unknown metric '" + metric + "' (must be one of " +
		                            CommunityMetricList() + ")");
	}
	const Metric m = static_cast<Metric>(metric_index);

	const uint32_t n = n_samples;
	const uint32_t f = n_features;
	const double *M = matrix.data();
	auto row = [&](uint32_t i) {
		return M + static_cast<size_t>(i) * f;
	};

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
	// PyCogent Gower distance (both are matrix-wide, not purely per-pair).
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

	// Base offset of row i's block in the condensed upper-triangle output: rows
	// 0..i-1 contribute sum_{t<i}(n-1-t) = i*(n-1) - i*(i-1)/2 pairs. Every pair
	// (i,j) thus has a FIXED destination slot, independent of which thread
	// computes it, so the result is bit-identical for any n_threads.
	auto row_base = [n](uint32_t i) -> size_t {
		// For i == 0 the (i - 1) subterm underflows to UINT32_MAX in uint32_t, but
		// its leading static_cast<size_t>(i) factor is 0, so the whole term is
		// exactly 0. Do NOT "simplify" the (i - 1) away as an overflow fix.
		return static_cast<size_t>(i) * (n - 1) - static_cast<size_t>(i) * (i - 1) / 2;
	};

	// All pairs (i, j>i) of one outer row, each written to its fixed slot.
	auto compute_row = [&](uint32_t i) {
		const double *xi = row(i);
		size_t o = row_base(i);
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
				// undefined. Follow PyCogent dist_pearson (the metric Kuczynski 2010
				// used, verified against its source), NOT scipy: two flat rows are
				// identical (r=1 -> distance 0); a flat vs a non-flat row has no
				// correlation (r=0 -> distance 1). This keeps a constant-profile
				// sample as a well-defined, finite distance rather than a NaN that a
				// downstream reader would reject.
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
					// Both empty: no row profiles, but PyCogent dist_chisq defines
					// this as distance 0 (identical). Follow it for faithful
					// reproduction (verified against the PyCogent source).
					d = 0.0;
				} else if (ri <= 0.0 || rj <= 0.0) {
					d = 1.0; // one empty vs non-empty -> maximal (PyCogent dist_chisq)
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

	if (n_threads <= 1) {
		for (uint32_t i = 0; i + 1 < n; ++i) {
			compute_row(i);
		}
	} else {
		// Dynamic row-stealing: row i does (n-1-i) pairs, so a static contiguous
		// split would overload the low-index threads. An atomic row cursor keeps
		// the load balanced while the fixed output slots keep the result
		// deterministic. Cap the OS-thread count two ways: never more than there
		// are rows with pairs (n-1), and never more than the hardware can run
		// concurrently — extra threads on this CPU-bound loop only add
		// context-switch overhead, and spawning thousands risks a pids/ulimit
		// hit. hardware_concurrency() reports 0 when it cannot detect; fall back
		// to the requested count (leaving only the n-1 cap) in that case.
		unsigned hw = std::thread::hardware_concurrency();
		if (hw == 0) {
			hw = n_threads;
		}
		unsigned nt = std::min<unsigned>(n_threads, n - 1);
		nt = std::min(nt, hw);
		std::atomic<uint32_t> next_row {0};
		auto worker = [&]() {
			for (;;) {
				const uint32_t i = next_row.fetch_add(1);
				if (i + 1 >= n) {
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
			// A thread constructor failed (e.g. EAGAIN under a pids/ulimit cap).
			// Do NOT let it unwind past the joinable threads already in `pool` —
			// a joinable std::thread's destructor calls std::terminate() and
			// would crash the entire process. Swallow it and degrade gracefully:
			// the threads already spawned keep draining rows through the atomic
			// cursor and the calling thread finishes the remainder below, so the
			// full result is still computed, just with fewer workers. (The vector
			// is pre-reserved, so emplace_back itself cannot throw bad_alloc; the
			// only exception here is the thread constructor's.)
		}
		worker(); // the calling thread participates as one worker
		for (auto &th : pool) {
			th.join();
		}
	}
	return out;
}

} // namespace miint
