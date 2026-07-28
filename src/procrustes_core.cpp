#include "procrustes_core.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Dense>

// Ported/adapted from SciPy (BSD-3-Clause):
//   scipy.spatial._procrustes.procrustes
//   scipy.linalg._procrustes.orthogonal_procrustes
//   Copyright (c) 2001-2002 Enthought, Inc.; 2003- SciPy Developers.
// See procrustes_core.hpp for the licensing note. The fit/apply split
// generalizes SciPy's single-call procrustes to the partial (anchored) case.
namespace miint::procrustes {

namespace {
using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Eigen::Map;
} // namespace

ProcrustesFit FitProcrustes(const double *ref_paired, const double *other_paired, uint32_t n_paired, uint32_t d) {
	if (d == 0) {
		throw std::invalid_argument("procrustes: dimensionality d must be >= 1");
	}
	if (n_paired < d + 1) {
		throw std::invalid_argument("procrustes: need at least d+1 paired points to fit a d-dimensional transform");
	}

	// Reject non-finite inputs up front (SciPy orthogonal_procrustes uses
	// asarray_chkfinite). NaN/Inf would otherwise flow into the SVD and produce a
	// silently-garbage transform rather than a clear error — important on the SQL
	// path where coordinates come from user tables.
	const std::size_t count = static_cast<std::size_t>(n_paired) * d;
	for (std::size_t i = 0; i < count; ++i) {
		if (!std::isfinite(ref_paired[i]) || !std::isfinite(other_paired[i])) {
			throw std::invalid_argument("procrustes: input contains non-finite values (NaN or Inf)");
		}
	}

	Map<const RowMajorMatrix> ref(ref_paired, n_paired, d);
	Map<const RowMajorMatrix> other(other_paired, n_paired, d);

	// Center each block on its column means (SciPy: mtx -= mtx.mean(0)).
	Eigen::RowVectorXd ref_mean = ref.colwise().mean();
	Eigen::RowVectorXd other_mean = other.colwise().mean();
	RowMajorMatrix ref_c = ref.rowwise() - ref_mean;
	RowMajorMatrix other_c = other.rowwise() - other_mean;

	// Scale each to unit Frobenius norm (SciPy: tr(AA^T) = 1).
	const double ref_norm = ref_c.norm();
	const double other_norm = other_c.norm();
	if (ref_norm == 0.0 || other_norm == 0.0) {
		throw std::invalid_argument("procrustes: input matrices must contain >1 unique points");
	}
	ref_c /= ref_norm;
	other_c /= other_norm;

	// Optimal orthogonal map (SciPy orthogonal_procrustes(A, B)): SVD of A^T B,
	// R = U V^T, scale = sum of singular values. Reflections are admitted
	// (det R may be -1), matching SciPy.
	//
	// Degeneracy note: when A^T B has repeated (tied) singular values, the U/V
	// singular vectors are only determined up to a rotation within the degenerate
	// subspace (and a shared sign for simple values), so R itself can differ
	// between SVD backends (Eigen JacobiSVD here vs LAPACK in the SciPy oracle)
	// while both remain optimal. The disparity M^2 = 1 - s^2 depends only on the
	// singular values and is invariant; the per-coordinate output is not, so
	// parity fixtures with tied spectra can legitimately differ in orientation.
	const Eigen::MatrixXd cross = ref_c.transpose() * other_c; // d x d
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(cross, Eigen::ComputeFullU | Eigen::ComputeFullV);
	const Eigen::MatrixXd R = svd.matrixU() * svd.matrixV().transpose();

	ProcrustesFit fit;
	fit.d = d;
	fit.ref_translate.assign(ref_mean.data(), ref_mean.data() + d);
	fit.other_translate.assign(other_mean.data(), other_mean.data() + d);
	fit.ref_norm = ref_norm;
	fit.other_norm = other_norm;
	fit.R.resize(static_cast<size_t>(d) * d);
	Map<RowMajorMatrix>(fit.R.data(), d, d) = R;
	fit.scale = svd.singularValues().sum();
	return fit;
}

void ApplyToReference(const ProcrustesFit &fit, const double *ref, uint32_t n, double *out) {
	const uint32_t d = fit.d;
	Map<const RowMajorMatrix> X(ref, n, d);
	Map<const Eigen::RowVectorXd> t(fit.ref_translate.data(), d);
	Map<RowMajorMatrix>(out, n, d) = (X.rowwise() - t) / fit.ref_norm;
}

void ApplyToOther(const ProcrustesFit &fit, const double *other, uint32_t n, double *out) {
	const uint32_t d = fit.d;
	Map<const RowMajorMatrix> X(other, n, d);
	Map<const Eigen::RowVectorXd> t(fit.other_translate.data(), d);
	Map<const RowMajorMatrix> R(fit.R.data(), d, d);
	const RowMajorMatrix centered = (X.rowwise() - t) / fit.other_norm;
	// SciPy: mtx2 = mtx2_std @ R.T * s.
	Map<RowMajorMatrix>(out, n, d) = centered * R.transpose() * fit.scale;
}

double Disparity(const double *ref_std, const double *other_fit, uint32_t n, uint32_t d) {
	Map<const RowMajorMatrix> A(ref_std, n, d);
	Map<const RowMajorMatrix> B(other_fit, n, d);
	return (A - B).array().square().sum();
}

double FullDisparity(const double *ref, const double *other, uint32_t n, uint32_t d) {
	const ProcrustesFit fit = FitProcrustes(ref, other, n, d);
	std::vector<double> ref_std(static_cast<std::size_t>(n) * d);
	std::vector<double> other_fit(static_cast<std::size_t>(n) * d);
	ApplyToReference(fit, ref, n, ref_std.data());
	ApplyToOther(fit, other, n, other_fit.data());
	return Disparity(ref_std.data(), other_fit.data(), n, d);
}

double MonteCarloPValue(const double *ref, const double *other, uint32_t n, uint32_t d, double true_m2,
                        uint32_t permutations, uint64_t seed) {
	if (permutations == 0) {
		return std::numeric_limits<double>::quiet_NaN(); // mirrors q2's disable path
	}

	// Reproducible PRNG (q2 uses an unseeded default_rng, so parity is statistical
	// only — see the header note). Each trial applies a further Fisher-Yates
	// shuffle to the running row order, matching q2's cumulative in-place shuffle
	// of `other`, then rescores the full-overlap disparity.
	std::mt19937_64 rng(seed);
	std::vector<uint32_t> perm(n);
	for (uint32_t i = 0; i < n; ++i) {
		perm[i] = i;
	}
	std::vector<double> shuffled(static_cast<std::size_t>(n) * d);

	uint32_t below = 0;
	for (uint32_t trial = 0; trial < permutations; ++trial) {
		for (uint32_t i = n; i > 1; --i) {
			std::uniform_int_distribution<uint32_t> dist(0, i - 1);
			std::swap(perm[i - 1], perm[dist(rng)]);
		}
		for (uint32_t r = 0; r < n; ++r) {
			const double *src = other + static_cast<std::size_t>(perm[r]) * d;
			double *dst = shuffled.data() + static_cast<std::size_t>(r) * d;
			for (uint32_t k = 0; k < d; ++k) {
				dst[k] = src[k];
			}
		}
		if (FullDisparity(ref, shuffled.data(), n, d) < true_m2) {
			++below;
		}
	}
	return static_cast<double>(below + 1) / static_cast<double>(permutations + 1);
}

} // namespace miint::procrustes
