#pragma once

#include <cstdint>
#include <vector>

// Procrustes analysis core — a similarity transform (translation, uniform
// scaling, orthogonal rotation/reflection) that optimally superimposes one
// point cloud onto another.
//
// Ported/adapted from SciPy (BSD-3-Clause):
//   scipy.spatial._procrustes.procrustes
//   scipy.linalg._procrustes.orthogonal_procrustes
//   Copyright (c) 2001-2002 Enthought, Inc.; 2003- SciPy Developers.
//
// The "deconstructed" split (fit the transform on a paired subset, then apply
// it to all rows) generalizes the standard full-overlap case and supports the
// partial / reference-anchored use (author's own technique, q2-diversity#338).
// Full procrustes == the degenerate case where the paired subset is every row.
namespace miint::procrustes {

// Parameters of the fitted similarity transform. Enough to map either input
// into the shared standardized frame; SciPy's `procrustes` bakes fit+apply into
// one call, we keep them separable so a fit on anchors can be applied to a
// larger set.
struct ProcrustesFit {
	uint32_t d = 0;                      // dimensionality (columns)
	std::vector<double> ref_translate;   // d: column means of the centered reference block
	std::vector<double> other_translate; // d: column means of the centered other block
	double ref_norm = 0.0;               // Frobenius norm of the centered reference block
	double other_norm = 0.0;             // Frobenius norm of the centered other block
	std::vector<double> R;               // d*d row-major orthogonal matrix (U * V^T)
	double scale = 0.0;                  // s = sum of singular values of ref_std^T * other_std
};

// Fit the transform on paired points. `ref_paired` / `other_paired` are
// (n_paired x d) row-major. Throws std::invalid_argument when n_paired < d + 1
// (too few points to determine a d-dimensional transform after centering) or
// when either centered block has zero Frobenius norm (< 2 unique points —
// mirrors SciPy's "must contain >1 unique points").
ProcrustesFit FitProcrustes(const double *ref_paired, const double *other_paired, uint32_t n_paired, uint32_t d);

// Map reference rows into the standardized frame: (X - ref_translate) / ref_norm.
// `ref` and `out` are (n x d) row-major; `out` must not alias `ref`.
void ApplyToReference(const ProcrustesFit &fit, const double *ref, uint32_t n, double *out);

// Map other rows into the standardized frame:
//   ((X - other_translate) / other_norm) * R^T * scale
// `other` and `out` are (n x d) row-major; `out` must not alias `other`.
void ApplyToOther(const ProcrustesFit &fit, const double *other, uint32_t n, double *out);

// Disparity M^2 = sum of squared pointwise differences between the standardized
// reference and fitted other over `n` rows. Both are (n x d) row-major. For the
// full case this equals SciPy's `procrustes` disparity; for the partial case it
// is evaluated over the paired points the transform was fit on.
double Disparity(const double *ref_std, const double *other_fit, uint32_t n, uint32_t d);

// Convenience: the full-overlap disparity M^2 between two (n x d) row-major point
// clouds — fit the transform on all n points, apply, and score. Equals SciPy's
// `scipy.spatial.procrustes` disparity. Throws the same std::invalid_argument as
// FitProcrustes on degenerate/non-finite input.
double FullDisparity(const double *ref, const double *other, uint32_t n, uint32_t d);

// Monte Carlo PROTEST p-value for a full procrustes fit (q2-diversity
// `_procrustes_monte_carlo` parity). Permutes `other`'s rows `permutations`
// times, recomputes the full-overlap disparity each trial, counts trials whose
// disparity is strictly below `true_m2`, and returns (count + 1)/(permutations
// + 1) — the +1 avoids a zero p-value (scikit-bio convention). `ref`/`other` are
// the (n x d) row-major clouds the true fit used. Reproducible given `seed`
// within one C++ standard library (unlike q2, which uses an unseeded
// default_rng); the exact draw sequence is NOT portable across libstdc++/libc++/
// MSVC because std::uniform_int_distribution's algorithm is implementation-
// defined. `permutations == 0` returns NaN. NOTE: exact q2 parity is impossible
// (q2 is unseeded; NumPy's PRNG differs from this std::mt19937_64) — agreement is
// statistical, within Monte Carlo error.
double MonteCarloPValue(const double *ref, const double *other, uint32_t n, uint32_t d, double true_m2,
                        uint32_t permutations, uint64_t seed);

} // namespace miint::procrustes
