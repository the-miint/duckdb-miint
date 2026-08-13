#pragma once
//
// Pure (DuckDB-free) k-means clustering of points in Euclidean space, used to
// score ordination recovery in the Kuczynski-2010 reproduction (k-means on PCoA
// coordinates -> fraction of sample-pairs whose co-membership matches truth).
//
// Header/impl split with NO DuckDB dependency so the Catch2 unit-test target
// links it without libduckdb (mirrors community_distances / simulate_resemblance).
// The DuckDB table-function wrapper (cluster_kmeans_function.cpp) reads a
// coordinate table -- the (sample_id, axis, coordinate) long form emitted by
// pcoa -- pivots it to a dense point matrix, calls in here, and emits
// (sample_id, cluster).

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

//! Result of a k-means run.
struct KMeansResult {
	std::vector<int32_t> assignments; //!< cluster id in [0, k) for each of the n_points, in input order
	double inertia = 0.0;             //!< sum over points of squared Euclidean distance to its centroid
	int32_t iterations = 0;           //!< Lloyd iterations of the best restart (excludes the assignment-only pass)
	int32_t restarts = 0;             //!< number of k-means++ restarts actually run
};

//! Lloyd's k-means (Lloyd 1982) with k-means++ seeding (Arthur & Vassilvitskii
//! 2007), `n_init` random restarts, and the
//! lowest-inertia restart returned. Deterministic for a given `seed` on a fixed
//! toolchain (all randomness flows from one std::mt19937_64(seed)); note that the
//! std::uniform_int/real_distribution mappings are not standardized across C++
//! standard libraries, so the exact seed draws may differ between e.g. libstdc++
//! and libc++. Empty clusters are re-seeded to the farthest still-unclaimed point.
//! When the points span at least `k` distinct locations this returns `k`
//! non-empty clusters; when there are fewer than `k` distinct locations the
//! result has fewer than `k` non-empty clusters -- matching scikit-learn's
//! KMeans, which likewise returns fewer than `k` labels (and warns) on duplicate
//! points, since `k` non-empty clusters are then impossible. (scikit-learn is the
//! parity-test oracle; see THIRD_PARTY_LICENSES.md.)
//!
//! `points` is row-major, `n_points` rows x `n_dims` columns. Cluster ids are
//! canonicalized so that clusters appear in order of first-seen point (id 0 is
//! the first point's cluster), making the labelling deterministic and
//! permutation-stable across equivalent solutions.
//!
//! Throws std::invalid_argument on: n_points == 0; n_dims == 0; k < 1; k >
//! n_points; points.size() != n_points*n_dims; n_init < 1; max_iter < 1.
KMeansResult KMeans(const std::vector<double> &points, uint32_t n_points, uint32_t n_dims, int32_t k, int64_t seed,
                    int32_t max_iter, int32_t n_init);

} // namespace miint

namespace duckdb {
class ExtensionLoader;
//! Registers the cluster_kmeans table function into the extension catalog.
void RegisterClusterKmeans(ExtensionLoader &loader);
} // namespace duckdb
