#pragma once
//
// Subset selection from ordination coordinates: choose `k` of `n` samples by a
// stated rule, returning their indices in selection order.
//
// Header/impl split with NO DuckDB dependency so the Catch2 unit-test target
// links it without libduckdb (mirrors cluster_kmeans / community_distances).
// The DuckDB table-function wrapper (pick_anchors_function.cpp) reads a
// coordinate table -- the (sample_id, axis, coordinate) long form emitted by
// pcoa -- pivots it to a dense point matrix, calls in here, and emits
// (anchor_rank, sample_id).
//
// WHY COORDINATES rather than a distance matrix: both rules below need only the
// N x d point cloud, so both are linear in N and run at millions of samples. A
// distance-matrix formulation would need a dense N x N fill and cap out around
// 25k samples for no gain.

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

//! Proportional stratified sampling (Neyman 1934; Cochran 1977, ch. 5) over the
//! point cloud: bin each of the `n_dims` axes into `n_bins` equal-frequency bins,
//! take the resulting grid cell as the stratum, and draw from every stratum in
//! proportion to its size.
//!
//! WHY THIS RULE IS THE DEFAULT for progressive-PCoA anchors. Four
//! literature-backed rules were measured against a full ordination of the
//! rarefied EMP 90 bp table (23,814 samples, unweighted UniFrac, 1,000 anchors),
//! scored by procrustes M2 at d=3: farthest-point/k-center 0.1745, statistical
//! leverage / CUR 0.0195, stratum medoids / k-medoids 0.0176, plain seeded random
//! 0.0079, proportional stratified 0.0045. Over five matched draws stratified beat
//! random 1.74x on the mean with non-overlapping ranges (stratified max 0.00604 <
//! random min 0.00707).
//!
//! The mechanism matters, because it is what picks this rule out of the four:
//! leverage, medoids and farthest-point each systematically prefer a KIND of
//! point -- high-influence, central, extreme -- and all three lose to an unbiased
//! draw. Stratified sampling is the only one that stays unbiased WITHIN a stratum
//! and merely equalizes coverage ACROSS strata. So the selection order below is
//! deliberately driven by a salted hash of the sample id, not by any geometric
//! quantity: within a stratum the choice must carry no preference at all.
//!
//! Ordering: within a stratum, samples are ranked by that salted hash; globally,
//! samples are ordered by `(within-stratum rank + 1) / (stratum size)`. Taking
//! the first `k` of that order draws floor(f * size) from each stratum for the
//! single cutoff fraction f that totals exactly `k` -- i.e. proportional
//! allocation. Because the order does not depend on `k`, selection is
//! prefix-stable: the first `m` of a `k`-selection are exactly the `m`-selection.
//!
//! `points` is row-major, `n_points` rows x `n_dims` columns; `ids` is parallel to
//! its rows and supplies the hash input, so the result depends on sample identity
//! rather than on row position. Deterministic for a given `seed` on any
//! toolchain: the hash is hand-rolled (FNV-1a + splitmix64) precisely so it does
//! not inherit std::hash's implementation-defined behavior.
//!
//! Note that the strata are a `n_bins ^ n_dims` grid. Keeping the number of
//! non-empty strata well below `k` is what makes the allocation proportional;
//! with far more strata than `k` the draw degenerates gracefully toward simple
//! random sampling over the most-populated cells.
//!
//! Throws std::invalid_argument on: n_points == 0; n_dims == 0; k < 1; k >
//! n_points; n_bins < 1; points.size() != n_points*n_dims; ids.size() !=
//! n_points; a `n_bins ^ n_dims` grid too large to index.
std::vector<uint32_t> SelectStratified(const std::vector<double> &points, const std::vector<std::string> &ids,
                                       uint32_t n_points, uint32_t n_dims, uint32_t k, uint32_t n_bins, int64_t seed);

//! Greedy farthest-point (k-center) selection (Gonzalez 1985): repeatedly take
//! the sample whose nearest already-selected neighbour is farthest away, which
//! is a 2-approximation to minimizing the maximum distance from any sample to
//! its nearest selected one.
//!
//! Rank 0 is the most peripheral sample, defined as the one farthest from the
//! centroid. In Euclidean space that is exactly the sample with the largest total
//! squared distance to all others, since sum_j ||x_i - x_j||^2 = n*||x_i - c||^2
//! + sum_j ||x_j - c||^2 -- the same "most peripheral" rule a distance-matrix
//! implementation would get from row sums, but computed in O(n*d) instead of
//! O(n^2*d).
//!
//! Deterministic, with no seed: every tie -- including rank 0's -- breaks to the
//! lowest row index, so a selection is a reproducible property of the data and
//! needs no seed recorded alongside it. Prefix-stable, since the greedy sequence
//! does not depend on `k`. Comparisons use squared Euclidean distance, which is
//! order-equivalent to Euclidean distance for both rules.
//!
//! NOT for progressive-PCoA anchors -- measured 15x worse than a random draw
//! there (see SelectStratified). It is a legitimate rule for diverse subset
//! selection generally: a maximally spread reference panel, a coverage-based
//! subsample, a diverse review set.
//!
//! Cost: O(n*d) for the centroid pass, then O(n*k*d), holding one double per
//! sample. `points` is row-major, `n_points` rows x `n_dims` columns.
//!
//! Throws std::invalid_argument on: n_points == 0; n_dims == 0; k < 1; k >
//! n_points; points.size() != n_points*n_dims.
std::vector<uint32_t> SelectFarthestPoint(const std::vector<double> &points, uint32_t n_points, uint32_t n_dims,
                                          uint32_t k);

} // namespace miint

namespace duckdb {
class ExtensionLoader;
//! Registers the pick_anchors table function into the extension catalog.
void RegisterPickAnchors(ExtensionLoader &loader);
} // namespace duckdb
