#pragma once
//
// Pure (DuckDB-free) non-phylogenetic community β-diversity distances.
//
// Computes the classic taxon-based (aphylogenetic) sample×sample distances used
// throughout Kuczynski et al. 2010 (Nature Methods 7:813-819) — the metrics the
// paper's central argument turns on (Jaccard vs Morisita-Horn vs χ² estimate
// fundamentally different things). Kept header/impl split with NO DuckDB
// dependency so the Catch2 unit-test target links it without libduckdb
// (mirrors simulate_resemblance.{hpp,cpp}); the DuckDB table-function wrapper
// (community_distances_function.cpp) reads the feature table, calls in here, and
// emits the condensed (sample_a, sample_b, distance) triple that pcoa /
// permanova / beta_* already consume.

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

//! True iff `metric` names a supported community distance (case-sensitive; the
//! wrapper lower-cases before calling). The accepted set is exactly the eight
//! figure metrics of Kuczynski 2010.
bool IsValidCommunityMetric(const std::string &metric);

//! True iff `metric` can be computed one BLOCK of samples at a time: d(i,j)
//! depends only on samples i and j — never on which other samples share the
//! matrix, and never on which features that matrix happens to carry. That is
//! what lets progressive_pcoa_from_features compute a block over
//! (anchors + batch) and get exactly the distances the full matrix would give.
//!
//! False for `pearson` (the per-sample mean is rowsum/n_features, and the
//! features that are zero in both samples still contribute to the covariance and
//! the variances, so the value moves with the block's feature space), `chisq`
//! (column sums and grand total) and `gower` (per-feature ranges). Computed per
//! block, those three would silently measure a DIFFERENT distance in every
//! block — no error, just a wrong ordination. Also false for an unknown metric.
bool IsPairwiseLocalCommunityMetric(const std::string &metric);

//! Human-readable comma-separated list of the PAIRWISE-LOCAL metrics only (for
//! errors that have to tell a caller which metrics a block-wise path accepts).
std::string PairwiseLocalCommunityMetricList();

//! Human-readable comma-separated list of accepted metric names (for errors).
std::string CommunityMetricList();

//! All pairwise community distances over a dense sample×feature abundance
//! matrix, returned condensed in row-major upper-triangle order: with i the
//! outer (0-based) sample index and j > i the inner, out[k] is the distance
//! between samples i and j as k runs 0 .. n*(n-1)/2. `matrix` is row-major,
//! `n_samples` rows × `n_features` columns; abundances are used AS GIVEN — the
//! function does NOT pre-normalize. Each metric applies whatever internal
//! normalization its own definition requires; feeding raw counts vs relative
//! abundance is the caller's modeling choice.
//!
//! Metrics (x, y = two sample rows; sums run over the k features; X=Σx, Y=Σy):
//!  - "bray_curtis"   Σ|xk-yk| / Σ(xk+yk)                        [0,1]; empty pair -> 0
//!  - "euclidean"     sqrt(Σ(xk-yk)^2)                            [0,inf)
//!  - "jaccard"       binary presence/absence: (b+c)/(a+b+c)      [0,1]; empty pair -> 0
//!  - "soergel"       Σ|xk-yk| / Σ max(xk,yk)                     [0,1]; empty pair -> 0
//!  - "morisita_horn" 1 - 2Σ(xk*yk) / ((Σxk^2/X^2 + Σyk^2/Y^2)*X*Y)
//!                                                                [0,1]; both-empty -> 0,
//!                                                                one-empty -> 1
//!  - "pearson"       1 - r, r = Pearson correlation over features [0,2];
//!                                                                both rows constant
//!                                                                -> 0, one constant
//!                                                                -> 1 (never NaN)
//!  - "chisq"         correspondence-analysis χ² distance:
//!                    sqrt(Σk (GT/colk)*(xk/rowx - yk/rowy)^2), colk>0 only,
//!                    GT = grand total; both row sums 0 -> 0, one 0 -> 1    [0,inf)
//!  - "gower"         Σk |xk-yk| / rangek (un-normalized;
//!                    rangek = max_i M[i][k] - min_i M[i][k]), rangek>0 only    [0,inf)
//!
//! χ² and Gower depend on GLOBAL column statistics (column sums / column ranges
//! over ALL samples), so they are matrix-wide, not purely per-pair.
//!
//! Primary sources for the metric definitions: Bray & Curtis 1957 (bray_curtis);
//! Jaccard 1912 (jaccard); Faith, Minchin & Belbin 1987 (chisq, gower in an
//! ecological setting); Gower 1971 (gower); Horn 1966 / Morisita 1959, as given
//! by Magurran 2004 p.246 (morisita_horn). The zero-variance / zero-row-sum
//! conventions above follow cogent3 `cogent3.maths.distance_transform` (BSD-3),
//! the reference implementation of the metrics Kuczynski et al. 2010 used; see
//! THIRD_PARTY_LICENSES.md.
//!
//! `n_threads` parallelizes the O(n^2 * f) pair loop. 0 or 1 runs fully serial
//! (no threads spawned; the caller resolves 0 = "follow DuckDB" before calling,
//! but 0 is also accepted here and treated as serial). >= 2 spawns up to
//! n_threads workers, itself capped at the number of rows with pairs (n-1) and
//! at std::thread::hardware_concurrency() — more threads than cores only add
//! context-switch overhead on this CPU-bound loop, and thousands risk a
//! pids/ulimit hit. The result is BIT-IDENTICAL for any thread count: each pair
//! writes to a fixed condensed slot, so threading only reorders WHEN each slot
//! is filled, never the value. The per-sample and global-column pre-passes are
//! always serial (they are O(n*f), not the hot loop).
//!
//! Throws std::invalid_argument on: n_samples < 2; matrix.size() !=
//! n_samples*n_features; or an unknown metric.
std::vector<double> CommunityDistancesCondensed(const std::vector<double> &matrix, uint32_t n_samples,
                                                uint32_t n_features, const std::string &metric, unsigned n_threads = 1);

//! The same distances from a CSR (sample-major) sparse matrix, for callers that
//! hold one block of a large, sparse feature table.
//!
//! WHY: the dense entry point above costs `pairs * n_features` because it walks
//! the whole feature space for every pair. On a real block that space is mostly
//! zeros — measured on a 1.2M-sample table, one 1100-sample block spans 11,018
//! features but averages 89 nonzeros per sample — so a merge over the union of
//! two samples' nonzeros does roughly 62x less arithmetic for the same answer.
//!
//! Layout: `indptr` has n_samples+1 entries, starts at 0, never decreases, and
//! ends at indices.size(); `indices` and `values` are parallel; within each row
//! the indices must be STRICTLY ASCENDING and less than n_features. Ascending is
//! not a convenience — the merge relies on it, and an unsorted row would yield a
//! wrong distance rather than an error, so it is validated.
//!
//! Results are BIT-IDENTICAL to CommunityDistancesCondensed on the same matrix.
//! Both accumulate in ascending feature order, and the terms the dense loop adds
//! for features absent from both samples are exactly 0.0. Also bit-identical for
//! any `n_threads`, on the same grounds as the dense loop.
//!
//! Only PAIRWISE-LOCAL metrics are accepted (see IsPairwiseLocalCommunityMetric).
//! chisq/gower/pearson are REFUSED rather than computed: this entry point exists
//! to serve one block at a time, and those read statistics over the whole table,
//! so a per-block value would silently be a different metric. Refusing in the
//! pure core puts that guarantee below every caller.
//!
//! Throws std::invalid_argument on: n_samples < 2; a malformed CSR (any of the
//! layout rules above); an unknown metric; or a metric that is not pairwise-local.
std::vector<double> CommunityDistancesCondensedSparse(const std::vector<uint32_t> &indptr,
                                                      const std::vector<uint32_t> &indices,
                                                      const std::vector<double> &values, uint32_t n_samples,
                                                      uint32_t n_features, const std::string &metric,
                                                      unsigned n_threads = 1);

} // namespace miint

namespace duckdb {
class ExtensionLoader;
//! Registers the community_distances table function into the extension catalog.
void RegisterCommunityDistances(ExtensionLoader &loader);
} // namespace duckdb
