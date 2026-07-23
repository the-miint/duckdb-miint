#include "cluster_kmeans.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace miint {

namespace {

double SqDist(const double *a, const double *b, uint32_t d) {
	double s = 0.0;
	for (uint32_t i = 0; i < d; ++i) {
		const double diff = a[i] - b[i];
		s += diff * diff;
	}
	return s;
}

// Assign every point to its nearest centroid (ties -> lowest index, deterministic)
// and record the squared distance to that centroid in `dist2`. That snapshot is
// what the empty-cluster reseed consults, so its "farthest point" choice reflects
// one consistent set of centroids rather than the half-rewritten update array.
void AssignNearest(const double *X, uint32_t n, uint32_t d, const std::vector<double> &cent, int32_t k,
                   std::vector<int32_t> &assign, std::vector<double> &dist2) {
	for (uint32_t i = 0; i < n; ++i) {
		double best = std::numeric_limits<double>::infinity();
		int32_t bc = 0;
		for (int32_t c = 0; c < k; ++c) {
			const double dd = SqDist(&X[static_cast<size_t>(i) * d], &cent[static_cast<size_t>(c) * d], d);
			if (dd < best) {
				best = dd;
				bc = c;
			}
		}
		assign[i] = bc;
		dist2[i] = best;
	}
}

// One k-means++ seeded Lloyd run. Fills `assign` (size n) and returns inertia.
// `iters` receives the number of Lloyd iterations executed.
double KMeansOnce(const double *X, uint32_t n, uint32_t d, int32_t k, std::mt19937_64 &rng, int32_t max_iter,
                  std::vector<int32_t> &assign, int32_t &iters) {
	std::vector<double> cent(static_cast<size_t>(k) * d, 0.0);
	std::uniform_real_distribution<double> unit(0.0, 1.0);

	// k-means++ seeding: first centroid uniform, each subsequent centroid drawn
	// with probability proportional to squared distance to the nearest chosen one.
	{
		std::uniform_int_distribution<uint32_t> pick(0, n - 1);
		const uint32_t c0 = pick(rng);
		for (uint32_t j = 0; j < d; ++j) {
			cent[j] = X[static_cast<size_t>(c0) * d + j];
		}
	}
	std::vector<double> dmin(n, std::numeric_limits<double>::infinity());
	for (int32_t c = 1; c < k; ++c) {
		double total = 0.0;
		for (uint32_t i = 0; i < n; ++i) {
			const double dd = SqDist(&X[static_cast<size_t>(i) * d], &cent[static_cast<size_t>(c - 1) * d], d);
			if (dd < dmin[i]) {
				dmin[i] = dd;
			}
			total += dmin[i];
		}
		uint32_t chosen = 0;
		if (total <= 0.0) {
			// All points coincide with already-chosen centroids; pick arbitrarily.
			std::uniform_int_distribution<uint32_t> pick(0, n - 1);
			chosen = pick(rng);
		} else {
			const double target = unit(rng) * total;
			double acc = 0.0;
			chosen = n - 1;
			for (uint32_t i = 0; i < n; ++i) {
				acc += dmin[i];
				if (acc >= target) {
					chosen = i;
					break;
				}
			}
		}
		for (uint32_t j = 0; j < d; ++j) {
			cent[static_cast<size_t>(c) * d + j] = X[static_cast<size_t>(chosen) * d + j];
		}
	}

	assign.assign(n, 0);
	std::vector<int32_t> prev(n, -1);
	std::vector<double> sums(static_cast<size_t>(k) * d);
	std::vector<uint32_t> counts(k);
	std::vector<double> dist2(n, 0.0); // squared distance of each point to its assigned centroid
	std::vector<char> claimed(n, 0);   // points already taken as a reseed target this update pass
	iters = 0;
	for (int32_t it = 0; it < max_iter; ++it) {
		iters = it + 1;
		AssignNearest(X, n, d, cent, k, assign, dist2);
		if (it > 0 && assign == prev) {
			break; // converged
		}
		prev = assign;

		// Update step: centroid = mean of its members.
		std::fill(sums.begin(), sums.end(), 0.0);
		std::fill(counts.begin(), counts.end(), 0u);
		for (uint32_t i = 0; i < n; ++i) {
			const int32_t c = assign[i];
			counts[c]++;
			for (uint32_t j = 0; j < d; ++j) {
				sums[static_cast<size_t>(c) * d + j] += X[static_cast<size_t>(i) * d + j];
			}
		}
		for (int32_t c = 0; c < k; ++c) {
			if (counts[c] > 0) {
				for (uint32_t j = 0; j < d; ++j) {
					cent[static_cast<size_t>(c) * d + j] = sums[static_cast<size_t>(c) * d + j] / counts[c];
				}
			}
		}
		// Reseed empty clusters to the still-unclaimed point farthest from its own
		// centroid, using the frozen assignment-step distances (dist2) rather than
		// the partly-rewritten `cent`; claiming each chosen point keeps two empty
		// clusters in the same pass from collapsing onto one point. The final
		// assignment pass below then populates the reseeded centroids.
		std::fill(claimed.begin(), claimed.end(), 0);
		for (int32_t c = 0; c < k; ++c) {
			if (counts[c] > 0) {
				continue;
			}
			double worst = -1.0;
			uint32_t worst_i = 0;
			bool found = false;
			for (uint32_t i = 0; i < n; ++i) {
				if (!claimed[i] && dist2[i] > worst) {
					worst = dist2[i];
					worst_i = i;
					found = true;
				}
			}
			if (found) {
				claimed[worst_i] = 1;
				for (uint32_t j = 0; j < d; ++j) {
					cent[static_cast<size_t>(c) * d + j] = X[static_cast<size_t>(worst_i) * d + j];
				}
			}
		}
	}

	// Assign once more against the final centroids so the returned assignment
	// reflects the last update (including any reseed) and is consistent with the
	// inertia computed below. On a converged run this reproduces the break-state
	// assignment; it is the "assignment-only pass" excluded from `iters`.
	AssignNearest(X, n, d, cent, k, assign, dist2);

	double inertia = 0.0;
	for (uint32_t i = 0; i < n; ++i) {
		inertia += SqDist(&X[static_cast<size_t>(i) * d], &cent[static_cast<size_t>(assign[i]) * d], d);
	}
	return inertia;
}

// Relabel clusters so ids appear in order of first-seen point (point 0 -> 0).
void Canonicalize(std::vector<int32_t> &assign, int32_t k) {
	std::vector<int32_t> remap(k, -1);
	int32_t next = 0;
	for (int32_t a : assign) {
		if (remap[a] < 0) {
			remap[a] = next++;
		}
	}
	for (int32_t &a : assign) {
		a = remap[a];
	}
}

} // namespace

KMeansResult KMeans(const std::vector<double> &points, uint32_t n_points, uint32_t n_dims, int32_t k, int64_t seed,
                    int32_t max_iter, int32_t n_init) {
	if (n_points == 0) {
		throw std::invalid_argument("cluster_kmeans: no points");
	}
	if (n_dims == 0) {
		throw std::invalid_argument("cluster_kmeans: n_dims must be >= 1");
	}
	if (k < 1) {
		throw std::invalid_argument("cluster_kmeans: k must be >= 1 (got " + std::to_string(k) + ")");
	}
	if (static_cast<uint32_t>(k) > n_points) {
		throw std::invalid_argument("cluster_kmeans: k (" + std::to_string(k) + ") must be <= number of points (" +
		                            std::to_string(n_points) + ")");
	}
	if (points.size() != static_cast<size_t>(n_points) * static_cast<size_t>(n_dims)) {
		throw std::invalid_argument("cluster_kmeans: points size (" + std::to_string(points.size()) +
		                            ") does not match n_points*n_dims");
	}
	if (max_iter < 1) {
		throw std::invalid_argument("cluster_kmeans: max_iter must be >= 1");
	}
	if (n_init < 1) {
		throw std::invalid_argument("cluster_kmeans: n_init must be >= 1");
	}

	std::mt19937_64 rng(static_cast<uint64_t>(seed));
	const double *X = points.data();

	KMeansResult best;
	best.inertia = std::numeric_limits<double>::infinity();
	best.restarts = n_init;
	std::vector<int32_t> assign;
	for (int32_t r = 0; r < n_init; ++r) {
		int32_t iters = 0;
		const double inertia = KMeansOnce(X, n_points, n_dims, k, rng, max_iter, assign, iters);
		if (inertia < best.inertia) {
			best.inertia = inertia;
			best.assignments = assign;
			best.iterations = iters;
		}
	}
	Canonicalize(best.assignments, k);
	return best;
}

} // namespace miint
