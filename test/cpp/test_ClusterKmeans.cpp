#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <set>
#include <stdexcept>
#include <vector>

#include "cluster_kmeans.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

// Three well-separated 2-D blobs (inter-blob ~10, intra ~1). k=3 must recover
// exactly these groups. Cluster ids are canonicalized by first-seen point, so
// the assignment vector is deterministic: blob A -> 0, blob B -> 1, blob C -> 2.
namespace {
const std::vector<double> BLOBS = {
    0,  0,  1,  0,  0,  1, // blob A (points 0-2)  near (0,0)
    10, 0,  11, 0,  10, 1, // blob B (points 3-5)  near (10,0)
    0,  10, 1,  10, 0,  11 // blob C (points 6-8)  near (0,10)
};
} // namespace

TEST_CASE("kmeans recovers well-separated blobs with canonical labels", "[cluster_kmeans]") {
	auto r = miint::KMeans(BLOBS, 9, 2, 3, /*seed*/ 42, /*max_iter*/ 100, /*n_init*/ 10);
	REQUIRE(r.assignments.size() == 9);
	// Exact canonical labelling: first-seen blob is cluster 0, etc.
	CHECK(r.assignments == std::vector<int32_t> {0, 0, 0, 1, 1, 1, 2, 2, 2});
	// Inertia = sum of within-blob SS. Each blob is 3 points {(0,0),(1,0),(0,1)}
	// translated; centroid (1/3,1/3); SS per blob = 3*((2/3)^2+(1/3)^2)*... ->
	// compute: dist^2 to centroid summed = (0-1/3)^2+(0-1/3)^2 + (1-1/3)^2+(0-1/3)^2
	// + (0-1/3)^2+(1-1/3)^2 = (2/9)+(5/9)+(5/9) = 12/9 = 4/3 per blob; 3 blobs -> 4.
	CHECK(r.inertia == Approx(4.0));
}

TEST_CASE("kmeans is deterministic for a fixed seed", "[cluster_kmeans]") {
	auto a = miint::KMeans(BLOBS, 9, 2, 3, 7, 100, 10);
	auto b = miint::KMeans(BLOBS, 9, 2, 3, 7, 100, 10);
	CHECK(a.assignments == b.assignments);
	CHECK(a.inertia == Approx(b.inertia));
}

TEST_CASE("kmeans k == n_points -> every point its own cluster, zero inertia", "[cluster_kmeans]") {
	auto r = miint::KMeans(BLOBS, 9, 2, 9, 1, 100, 5);
	CHECK(r.inertia == Approx(0.0).margin(1e-12));
	// Canonical labels: each point distinct, in first-seen order.
	CHECK(r.assignments == std::vector<int32_t> {0, 1, 2, 3, 4, 5, 6, 7, 8});
}

TEST_CASE("kmeans k == 1 -> one cluster, inertia = total scatter about the mean", "[cluster_kmeans]") {
	auto r = miint::KMeans(BLOBS, 9, 2, 1, 1, 100, 3);
	for (int32_t a : r.assignments) {
		CHECK(a == 0);
	}
	// Grand mean of the 9 points = (33/9, 33/9) = (3.6667, 3.6667). Total SS is
	// large (blobs are ~10 apart); just require it dwarfs the k=3 inertia (4.0).
	CHECK(r.inertia > 100.0);
}

TEST_CASE("kmeans cluster ids are contiguous 0..k-1", "[cluster_kmeans]") {
	auto r = miint::KMeans(BLOBS, 9, 2, 3, 3, 100, 10);
	int32_t mx = -1;
	for (int32_t a : r.assignments) {
		CHECK(a >= 0);
		mx = a > mx ? a : mx;
	}
	CHECK(mx == 2); // exactly 3 clusters, labelled 0,1,2
}

TEST_CASE("kmeans on fewer distinct locations than k returns fewer clusters (matches sklearn)", "[cluster_kmeans]") {
	// Five coincident points at 0 and one at 100: only two distinct locations, so
	// three non-empty clusters are impossible. scikit-learn returns two labels here
	// (with a ConvergenceWarning); we match that rather than fabricating a third
	// cluster the data cannot support. This pins the corrected header contract.
	const std::vector<double> pts = {0, 0, 0, 0, 0, 100};
	auto r = miint::KMeans(pts, 6, 1, 3, /*seed*/ 42, /*max_iter*/ 100, /*n_init*/ 10);
	const std::set<int32_t> distinct(r.assignments.begin(), r.assignments.end());
	CHECK(distinct.size() == 2);
	// Co-membership matches the true partition (permutation-invariant): the five
	// zeros share one cluster, the outlier is on its own.
	for (size_t i = 1; i < 5; ++i) {
		CHECK(r.assignments[i] == r.assignments[0]);
	}
	CHECK(r.assignments[5] != r.assignments[0]);
}

TEST_CASE("kmeans with a duplicate point but k distinct locations returns exactly k clusters", "[cluster_kmeans]") {
	// Distinct locations {0, 5, 10} with the point at 0 duplicated. k-means++ seeds
	// the three distinct locations, so all three clusters are non-empty -- a
	// duplicate only drops the cluster count below k when the number of DISTINCT
	// locations itself drops below k. Checked across seeds so it can't pass by luck.
	const std::vector<double> pts = {0, 0, 5, 10};
	for (int64_t seed : {int64_t(0), int64_t(1), int64_t(2), int64_t(7), int64_t(42)}) {
		auto r = miint::KMeans(pts, 4, 1, 3, seed, 100, 10);
		const std::set<int32_t> distinct(r.assignments.begin(), r.assignments.end());
		CHECK(distinct.size() == 3);
		CHECK(r.assignments[0] == r.assignments[1]); // the two coincident points co-cluster
	}
}

TEST_CASE("kmeans guards", "[cluster_kmeans]") {
	CHECK_THROWS_AS(miint::KMeans(BLOBS, 9, 2, 0, 1, 100, 10), std::invalid_argument);  // k < 1
	CHECK_THROWS_AS(miint::KMeans(BLOBS, 9, 2, 10, 1, 100, 10), std::invalid_argument); // k > n_points
	CHECK_THROWS_AS(miint::KMeans(BLOBS, 0, 2, 1, 1, 100, 10), std::invalid_argument);  // n_points 0
	CHECK_THROWS_AS(miint::KMeans(BLOBS, 9, 0, 1, 1, 100, 10), std::invalid_argument);  // n_dims 0
	CHECK_THROWS_AS(miint::KMeans(BLOBS, 9, 2, 3, 1, 0, 10), std::invalid_argument);    // max_iter < 1
	CHECK_THROWS_AS(miint::KMeans(BLOBS, 9, 2, 3, 1, 100, 0), std::invalid_argument);   // n_init < 1
	// size mismatch (n_points*n_dims != points.size())
	CHECK_THROWS_AS(miint::KMeans(BLOBS, 5, 2, 2, 1, 100, 10), std::invalid_argument);
	CHECK_THROWS_WITH(miint::KMeans(BLOBS, 5, 2, 2, 1, 100, 10), ContainsSubstring("size"));
}
