#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "procrustes_core.hpp"

using Catch::Approx;
using miint::procrustes::ApplyToOther;
using miint::procrustes::ApplyToReference;
using miint::procrustes::Disparity;
using miint::procrustes::FitProcrustes;
using miint::procrustes::FullDisparity;
using miint::procrustes::MonteCarloPValue;
using miint::procrustes::ProcrustesFit;

namespace {

// Row-major (n x 2) similarity image: out[i] = scale * (in[i] @ M) + t, where M
// is a 2x2 row-major matrix. A negative-determinant M encodes a reflection.
std::vector<double> image2d(const std::vector<double> &in, uint32_t n, const double M[4], double scale,
                            const double t[2]) {
	std::vector<double> out(static_cast<size_t>(n) * 2);
	for (uint32_t i = 0; i < n; ++i) {
		const double x = in[i * 2 + 0];
		const double y = in[i * 2 + 1];
		out[i * 2 + 0] = scale * (x * M[0] + y * M[2]) + t[0];
		out[i * 2 + 1] = scale * (x * M[1] + y * M[3]) + t[1];
	}
	return out;
}

// Row-major (n x d) similarity image for arbitrary d: out[i] = scale * (in[i] @ M)
// + t, where M is a d x d row-major matrix (a d-dimensional generalization of
// image2d, used to exercise the RowMajor Map round-trip at d >= 3).
std::vector<double> imageNd(const std::vector<double> &in, uint32_t n, uint32_t d, const double *M, double scale,
                            const double *t) {
	std::vector<double> out(static_cast<size_t>(n) * d);
	for (uint32_t i = 0; i < n; ++i) {
		for (uint32_t j = 0; j < d; ++j) {
			double acc = 0.0;
			for (uint32_t k = 0; k < d; ++k) {
				acc += in[i * d + k] * M[k * d + j];
			}
			out[i * d + j] = scale * acc + t[j];
		}
	}
	return out;
}

// 5 distinct, non-collinear points spanning 2-D (norm > 0 after centering).
const std::vector<double> kRef5 = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.5, 2.0};

} // namespace

TEST_CASE("procrustes recovers a pure similarity image (M^2 ~ 0)", "[procrustes]") {
	const double c = std::cos(0.6), s = std::sin(0.6);
	const double rot[4] = {c, -s, s, c}; // proper rotation, det = +1
	const double t[2] = {7.0, -3.0};
	auto other = image2d(kRef5, 5, rot, 2.5, t);

	auto fit = FitProcrustes(kRef5.data(), other.data(), 5, 2);
	std::vector<double> ref_std(10), other_fit(10);
	ApplyToReference(fit, kRef5.data(), 5, ref_std.data());
	ApplyToOther(fit, other.data(), 5, other_fit.data());

	// A similarity image standardizes to the same shape, so the optimal fit is exact.
	REQUIRE(Disparity(ref_std.data(), other_fit.data(), 5, 2) == Approx(0.0).margin(1e-12));
	for (int i = 0; i < 10; ++i) {
		REQUIRE(other_fit[i] == Approx(ref_std[i]).margin(1e-9));
	}
}

TEST_CASE("procrustes recovers an image that includes a reflection", "[procrustes]") {
	const double c = std::cos(0.9), s = std::sin(0.9);
	const double refl[4] = {c, s, s, -c}; // det = -1 (rotation + reflection)
	const double t[2] = {-4.0, 11.0};
	auto other = image2d(kRef5, 5, refl, 1.8, t);

	auto fit = FitProcrustes(kRef5.data(), other.data(), 5, 2);
	std::vector<double> ref_std(10), other_fit(10);
	ApplyToReference(fit, kRef5.data(), 5, ref_std.data());
	ApplyToOther(fit, other.data(), 5, other_fit.data());
	// Procrustes admits reflections, so a reflected image is still recovered exactly.
	REQUIRE(Disparity(ref_std.data(), other_fit.data(), 5, 2) == Approx(0.0).margin(1e-12));
}

TEST_CASE("procrustes recovers a similarity image in 3-D (M^2 ~ 0)", "[procrustes]") {
	// The 2-D cases exercise the RowMajor Map with d == 2; this one guards the
	// column-major/row-major round-trip and the d x d SVD at d == 3.
	// 5 non-coplanar points spanning 3-D.
	const std::vector<double> ref3 = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
	// Proper rotation M = Rz(a) * Rx(b) (product of two rotations -> det = +1),
	// written out row-major.
	const double ca = std::cos(0.5), sa = std::sin(0.5), cb = std::cos(0.4), sb = std::sin(0.4);
	const double M[9] = {ca, -sa * cb, sa * sb, sa, ca * cb, -ca * sb, 0.0, sb, cb};
	const double t[3] = {5.0, -2.0, 3.0};
	auto other = imageNd(ref3, 5, 3, M, 2.0, t);

	auto fit = FitProcrustes(ref3.data(), other.data(), 5, 3);
	REQUIRE(fit.d == 3u);
	std::vector<double> ref_std(15), other_fit(15);
	ApplyToReference(fit, ref3.data(), 5, ref_std.data());
	ApplyToOther(fit, other.data(), 5, other_fit.data());
	REQUIRE(Disparity(ref_std.data(), other_fit.data(), 5, 3) == Approx(0.0).margin(1e-12));
	for (int i = 0; i < 15; ++i) {
		REQUIRE(other_fit[i] == Approx(ref_std[i]).margin(1e-9));
	}
}

TEST_CASE("procrustes fits at the n == d + 1 solvability boundary", "[procrustes]") {
	// Exactly d + 1 = 3 non-collinear points is the minimum that determines a
	// 2-D transform after centering; it must succeed (one fewer throws — see the
	// rejection test).
	const std::vector<double> ref3 = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
	const double rot[4] = {std::cos(0.4), -std::sin(0.4), std::sin(0.4), std::cos(0.4)};
	const double t[2] = {2.0, -1.0};
	auto other = imageNd(ref3, 3, 2, rot, 1.5, t);

	auto fit = FitProcrustes(ref3.data(), other.data(), 3, 2);
	std::vector<double> ref_std(6), other_fit(6);
	ApplyToReference(fit, ref3.data(), 3, ref_std.data());
	ApplyToOther(fit, other.data(), 3, other_fit.data());
	REQUIRE(Disparity(ref_std.data(), other_fit.data(), 3, 2) == Approx(0.0).margin(1e-12));
}

TEST_CASE("procrustes rotation matrix R is orthogonal", "[procrustes]") {
	const double rot[4] = {std::cos(0.3), -std::sin(0.3), std::sin(0.3), std::cos(0.3)};
	const double t[2] = {1.0, 2.0};
	auto other = image2d(kRef5, 5, rot, 1.0, t);
	auto fit = FitProcrustes(kRef5.data(), other.data(), 5, 2);

	REQUIRE(fit.d == 2u);
	REQUIRE(fit.R.size() == 4u);
	// R^T R == I.
	for (uint32_t i = 0; i < 2; ++i) {
		for (uint32_t j = 0; j < 2; ++j) {
			double dot = 0.0;
			for (uint32_t k = 0; k < 2; ++k) {
				dot += fit.R[k * 2 + i] * fit.R[k * 2 + j];
			}
			REQUIRE(dot == Approx(i == j ? 1.0 : 0.0).margin(1e-12));
		}
	}
}

TEST_CASE("procrustes disparity is in [0, 1] for unrelated clouds", "[procrustes]") {
	// Why the bound is [0, 1] (not [0, 2]): after standardization ref_std has unit
	// Frobenius norm and other_fit = (other_std @ R^T) * s with R orthogonal, so
	// ||other_fit||^2 = s^2. The cross term is <ref_std, other_std @ R^T> = tr(Sigma)
	// = s (the SVD trace identity for R = U V^T), hence
	//   M^2 = ||ref_std||^2 - 2s*s + s^2 = 1 - 2s^2 + s^2 = 1 - s^2.
	// s = sum of singular values of the unit-norm cross product is <= 1 by
	// Cauchy-Schwarz, so M^2 = 1 - s^2 lies in [0, 1]. This matches SciPy's
	// `procrustes` disparity exactly for the full-overlap case.
	// `other` is a genuinely different cloud, not an image of `ref`.
	const std::vector<double> other = {0.2, 1.7, 2.0, 0.1, -1.0, 0.5, 0.3, 0.3, 1.4, -0.8};
	auto fit = FitProcrustes(kRef5.data(), other.data(), 5, 2);
	std::vector<double> ref_std(10), other_fit(10);
	ApplyToReference(fit, kRef5.data(), 5, ref_std.data());
	ApplyToOther(fit, other.data(), 5, other_fit.data());
	const double m2 = Disparity(ref_std.data(), other_fit.data(), 5, 2);
	REQUIRE(m2 >= -1e-12);
	REQUIRE(m2 <= 1.0 + 1e-9);
}

TEST_CASE("procrustes partial: fit on anchors, apply to all (extras land correctly)", "[procrustes]") {
	// ref-space truth for 5 anchors + 3 extra points.
	std::vector<double> ref_all = {0.0, 0.0, 1.0,  0.0, 0.0, 1.0, 1.0, 1.0, 0.5, 2.0, // anchors
	                               2.0, 0.5, -1.0, 1.5, 0.7, -0.9};                   // extras
	const double c = std::cos(0.7), s = std::sin(0.7);
	const double rot[4] = {c, -s, s, c};
	const double t[2] = {3.0, -6.0};
	auto other_all = image2d(ref_all, 8, rot, 4.0, t); // whole cloud transformed identically

	// Reference ordination knows only the 5 anchors; other knows all 8.
	auto fit = FitProcrustes(ref_all.data(), other_all.data(), 5, 2); // fit on anchors
	std::vector<double> ref_std(10), other_fit(16);
	ApplyToReference(fit, ref_all.data(), 5, ref_std.data());
	ApplyToOther(fit, other_all.data(), 8, other_fit.data());

	// Anchors align exactly.
	REQUIRE(Disparity(ref_std.data(), other_fit.data(), 5, 2) == Approx(0.0).margin(1e-12));

	// The extras (rows 5..7 of `other`) must land at the anchor-standardized
	// ref-space positions — this is the crux of the partial technique.
	for (uint32_t r = 5; r < 8; ++r) {
		for (uint32_t j = 0; j < 2; ++j) {
			const double expected = (ref_all[r * 2 + j] - fit.ref_translate[j]) / fit.ref_norm;
			REQUIRE(other_fit[r * 2 + j] == Approx(expected).margin(1e-9));
		}
	}
}

TEST_CASE("procrustes rejects too-few paired points and degenerate input", "[procrustes]") {
	// Need >= d + 1 = 3 paired points for a 2-D transform.
	const std::vector<double> two_pts = {0.0, 0.0, 1.0, 1.0};
	REQUIRE_THROWS_AS(FitProcrustes(two_pts.data(), two_pts.data(), 2, 2), std::invalid_argument);

	// Three identical points -> zero norm after centering (< 2 unique points).
	const std::vector<double> same = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	REQUIRE_THROWS_AS(FitProcrustes(same.data(), same.data(), 3, 2), std::invalid_argument);

	// Non-finite inputs are rejected up front (SciPy asarray_chkfinite parity);
	// otherwise NaN/Inf silently poisons the SVD and yields a garbage transform.
	std::vector<double> with_nan = kRef5;
	with_nan[4] = std::numeric_limits<double>::quiet_NaN();
	REQUIRE_THROWS_AS(FitProcrustes(kRef5.data(), with_nan.data(), 5, 2), std::invalid_argument);
	std::vector<double> with_inf = kRef5;
	with_inf[7] = std::numeric_limits<double>::infinity();
	REQUIRE_THROWS_AS(FitProcrustes(with_inf.data(), kRef5.data(), 5, 2), std::invalid_argument);
}

TEST_CASE("procrustes FullDisparity equals the fit/apply/score pipeline", "[procrustes]") {
	// FullDisparity is the one-call form of fit-on-all + apply + Disparity; it must
	// reproduce that pipeline exactly (it's what the Monte Carlo test recomputes
	// per permutation) and stay within SciPy's [0, 1] disparity bound.
	const std::vector<double> other = {0.2, 1.7, 2.0, 0.1, -1.0, 0.5, 0.3, 0.3, 1.4, -0.8};
	auto fit = FitProcrustes(kRef5.data(), other.data(), 5, 2);
	std::vector<double> ref_std(10), other_fit(10);
	ApplyToReference(fit, kRef5.data(), 5, ref_std.data());
	ApplyToOther(fit, other.data(), 5, other_fit.data());
	const double manual = Disparity(ref_std.data(), other_fit.data(), 5, 2);

	const double full = FullDisparity(kRef5.data(), other.data(), 5, 2);
	REQUIRE(full == Approx(manual).margin(1e-12));
	REQUIRE(full >= -1e-12);
	REQUIRE(full <= 1.0 + 1e-9);
}

TEST_CASE("procrustes Monte Carlo p-value: determinism, range, and perfect-fit floor", "[procrustes]") {
	const double c = std::cos(0.6), s = std::sin(0.6);
	const double rot[4] = {c, -s, s, c};
	const double t[2] = {7.0, -3.0};
	auto image = image2d(kRef5, 5, rot, 2.5, t); // a pure similarity image of kRef5

	const uint32_t perms = 49;
	const double true_m2 = FullDisparity(kRef5.data(), image.data(), 5, 2); // ~ 0

	// Deterministic given the seed (q2 is unseeded; we are not).
	const double p1 = MonteCarloPValue(kRef5.data(), image.data(), 5, 2, true_m2, perms, 123);
	const double p2 = MonteCarloPValue(kRef5.data(), image.data(), 5, 2, true_m2, perms, 123);
	REQUIRE(p1 == p2);

	// A perfect fit cannot be beaten by any row permutation, so zero trials fall
	// strictly below true_m2 and the p-value sits at its floor (0 + 1)/(perms + 1).
	REQUIRE(p1 == Approx(1.0 / (perms + 1)).margin(1e-12));

	// Range is (0, 1].
	REQUIRE(p1 > 0.0);
	REQUIRE(p1 <= 1.0);

	// permutations == 0 -> NaN (mirrors q2's disable path).
	REQUIRE(std::isnan(MonteCarloPValue(kRef5.data(), image.data(), 5, 2, true_m2, 0, 1)));
}
