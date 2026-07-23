#include "community_distances.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

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
                                                uint32_t n_features, const std::string &metric) {
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

	const double NaN = std::numeric_limits<double>::quiet_NaN();
	std::vector<double> out;
	out.reserve(PairCount(n));

	for (uint32_t i = 0; i + 1 < n; ++i) {
		const double *xi = row(i);
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
				// diverge from scipy's `correlation` and spuriously report zero
				// variance for a genuinely non-constant row. Matches scipy, which
				// also de-means before the dot products.
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
				// A constant profile has zero variance -> correlation is undefined
				// (scipy returns nan). The guard is explicit intent, not strictly
				// load-bearing: for an exactly-constant row every dxi is 0 so cov is
				// also exactly 0, and 0/0 would already yield NaN. We keep it so the
				// NaN is deliberate rather than an implicit IEEE 0/0. A downstream
				// reader (ReadDistanceTable) rejects the non-finite distance as
				// "not provided" rather than inventing a value.
				d = (varx > 0.0 && vary > 0.0) ? 1.0 - cov / std::sqrt(varx * vary) : NaN;
				break;
			}
			case Metric::Chisq: {
				const double ri = rowsum[i];
				const double rj = rowsum[j];
				if (ri <= 0.0 || rj <= 0.0) {
					d = NaN; // an all-zero sample has no row profile (see note above)
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

			out.push_back(d);
		}
	}
	return out;
}

} // namespace miint
