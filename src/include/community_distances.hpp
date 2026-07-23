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
//!                                                                constant row -> NaN
//!  - "chisq"         correspondence-analysis χ² distance:
//!                    sqrt(Σk (GT/colk)*(xk/rowx - yk/rowy)^2), colk>0 only,
//!                    GT = grand total; needs positive row sums               [0,inf)
//!  - "gower"         Σk |xk-yk| / rangek (PyCogent dist_gower, un-normalized;
//!                    rangek = max_i M[i][k] - min_i M[i][k]), rangek>0 only    [0,inf)
//!
//! χ² and Gower depend on GLOBAL column statistics (column sums / column ranges
//! over ALL samples), so they are matrix-wide, not purely per-pair.
//!
//! Throws std::invalid_argument on: n_samples < 2; matrix.size() !=
//! n_samples*n_features; or an unknown metric.
std::vector<double> CommunityDistancesCondensed(const std::vector<double> &matrix, uint32_t n_samples,
                                                uint32_t n_features, const std::string &metric);

} // namespace miint

namespace duckdb {
class ExtensionLoader;
//! Registers the community_distances table function into the extension catalog.
void RegisterCommunityDistances(ExtensionLoader &loader);
} // namespace duckdb
