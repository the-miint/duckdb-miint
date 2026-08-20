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

//! Index of pair (i, j), j > i, in the condensed upper-triangle layout every entry
//! point here returns: `CondensedRowBase(n, i) + (j - i - 1)`.
//!
//! WHY THIS IS PUBLIC: the layout is part of the return contract, so callers that
//! splice into or read out of a condensed vector have to reproduce
//! `i * (2n - i - 1) / 2`. Four independent copies of that expression existed before
//! this was exported, and an off-by-one in any of them silently mis-PAIRS samples
//! rather than failing — the worst possible failure mode for a distance matrix.
//!
//! Computed in 64-bit and returned as size_t: at n = 100,000 the last row base is
//! ~5e9, which overflows uint32 and, on a 32-bit size_t, the caller's own arithmetic
//! too — so a caller near that scale must check the result fits its index type.
size_t CondensedRowBase(uint32_t n_samples, uint32_t i);

//! True iff a block of this shape is evaluated by the dense Gram (GEMM) fast path
//! rather than the per-pair loop.
//!
//! WHY THIS IS PUBLIC: the dispatch is automatic and both kernels are internal, so
//! without this there is no way to assert that a given input actually took the
//! GEMM path — a test would pass just as happily against the pair loop, and a
//! mis-set threshold would be invisible. It is also the honest way for a caller to
//! reason about memory: a `true` here means an n_samples x n_features dense operand
//! is materialized.
//!
//! Only `euclidean`, `jaccard` and `morisita_horn` are ever eligible: each is a
//! function of one inner product plus per-sample scalars. `bray_curtis` and
//! `soergel` never are — Σ|x-y| and Σ max(x,y) are not bilinear — so this returns
//! false for them at every density.
//!
//! The rule is (metric eligible) AND n_samples >= 64 AND the dense operand fits a
//! fixed byte cap AND nnz / (n_samples * n_features) exceeds a density threshold.
//!
//! MEMORY, precisely: a `true` costs a dense n_samples x n_features operand, which
//! the byte cap bounds, PLUS an n_samples x n_samples Gram matrix, which it does not.
//! The Gram is not capped because it cannot be bounded independently of the request:
//! it is 2x the condensed result the caller is already asking for (n^2 doubles
//! against n(n-1)/2), so the path costs a bounded constant factor on an allocation
//! the caller has already chosen, never an unbounded one. Measured at a 2000-sample
//! block: +29 MB peak, against the 32 MB the n x n matrix predicts. Computing the
//! Gram in row panels would remove that term; it is not done because at block sizes
//! it does not matter, and at whole-table sizes the result array dominates anyway.
//! Density is the deciding term because the pair loop walks the union of two
//! samples' nonzeros while the GEMM walks the whole feature space, so the work ratio
//! is m*f/nnz = 1/density -- independent of both n and f. On a real microbiome block
//! (measured: 56,142 features at 0.19% density over 2000 samples) the GEMM is ~11x
//! slower and needs ~350x the memory; on a 38.8%-dense image block of 673 features it
//! is ~9.7x faster. The two kernels cost the same at ~0.96% density (measured on a
//! ladder that varies only density), and the threshold is set above that -- the
//! constant in community_distances.cpp carries the full table.
//! `max_operand_bytes` overrides the dense-operand ceiling for this one question;
//! 0 means the library default. A caller that runs several blocks AT ONCE must divide
//! its own budget by that concurrency, because this cap governs a single call and
//! knows nothing about how many are in flight.
bool CommunityDistancesUsesGramPath(const std::string &metric, uint32_t n_samples, uint32_t n_features, size_t nnz,
                                    size_t max_operand_bytes = 0);

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
//! On the merge path both accumulate in ascending feature order, and the terms the
//! dense loop adds for features absent from both samples are exactly 0.0. On the
//! dense-Gram path (see CommunityDistancesUsesGramPath) the guarantee survives for a
//! different reason: the dispatch rule is the same on both sides, so a block that
//! takes the Gram path here takes it there too, and the two then run the identical
//! kernel over the identical operand. Also bit-identical for any `n_threads`, on the
//! same grounds as the dense loop — the Gram path's product is serial regardless of
//! `n_threads`, which only fans out the fill.
//!
//! Bit-identity across the two ENTRY POINTS is not bit-identity across the two
//! KERNELS: for `euclidean` and `morisita_horn` a Gram-path distance can differ from
//! the merge's in its last bits, because `Sx2_i + Sx2_j - 2<x_i,x_j>` and
//! `Sum (x_k - y_k)^2` associate the same arithmetic differently. `jaccard` does not
//! — its inner product is a sum of 0/1 products, exact in double.
//!
//! Only PAIRWISE-LOCAL metrics are accepted (see IsPairwiseLocalCommunityMetric).
//! chisq/gower/pearson are REFUSED rather than computed: this entry point exists
//! to serve one block at a time, and those read statistics over the whole table,
//! so a per-block value would silently be a different metric. Refusing in the
//! pure core puts that guarantee below every caller.
//!
//! Throws std::invalid_argument on: n_samples < 2; a malformed CSR (any of the
//! layout rules above); an unknown metric; or a metric that is not pairwise-local.
//! `n_cached_tail` lets a caller that already HOLDS the mutual distances of the
//! last `n_cached_tail` samples skip re-deriving them: pairs with both indices in
//! that trailing square are left at 0.0 for the caller to fill, and every other
//! pair is computed exactly as it would have been (bit-identical, for any thread
//! count). 0 computes every pair.
//!
//! WHY: progressive_pcoa_from_features asks for one block per batch over (that
//! batch's samples ++ ALL anchors), so the anchor x anchor quadrant is identical in
//! every block of the run — see IsPairwiseLocalCommunityMetric for why it is
//! identical, which is exactly the property that makes caching it sound. That
//! quadrant is a(a-1)/2 of the block's (a+k)(a+k-1)/2 pairs: a quarter of the work
//! at the documented defaults (a = k = 1000), 44% at k = 500, and 64% at a = 200,
//! k = 50 — the small-batch regime the accuracy guidance actually recommends.
//!
//! The trailing position is not arbitrary: the core builds each request as the
//! batch's ids followed by the anchors, so the shared square is a SUFFIX, and
//! because j > i in the condensed layout the skipped pairs are exactly the rows
//! i >= n_samples - n_cached_tail. Skipping them is therefore a loop bound, not a
//! per-pair test.
//!
//! Throws std::invalid_argument if `n_cached_tail` exceeds `n_samples` — a tail
//! larger than the block means the caller's idea of what is cached disagrees with
//! the block it passed.
std::vector<double> CommunityDistancesCondensedSparse(const std::vector<uint32_t> &indptr,
                                                      const std::vector<uint32_t> &indices,
                                                      const std::vector<double> &values, uint32_t n_samples,
                                                      uint32_t n_features, const std::string &metric,
                                                      unsigned n_threads = 1, uint32_t n_cached_tail = 0,
                                                      size_t max_operand_bytes = 0);

} // namespace miint

namespace duckdb {
class ExtensionLoader;
//! Registers the community_distances table function into the extension catalog.
void RegisterCommunityDistances(ExtensionLoader &loader);
} // namespace duckdb
