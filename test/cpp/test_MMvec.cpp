#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "LBFGSpp/LBFGS.h"
#include "mmvec.hpp"
#include "mmvec_oracle.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using Eigen::VectorXd;
using miint::mmvec::ComputeSufficientStats;
using miint::mmvec::SparseCounts;

// ---------------------------------------------------------------------------
// Vendored LBFGS++ (ext/LBFGSpp/, upstream v0.4.0)
//
// These cases guard the *dependency*, not mmvec. They exist because three
// properties we rely on are invisible in our own source and would break
// silently on a vendor bump or a botched re-drop:
//
//   1. The flattened vendoring compiles at all. Upstream splits the headers
//      across include/ and include/LBFGSpp/; we vendor them flat and rely on
//      `include_directories(ext)` plus the headers' own relative includes to
//      resolve, with no patching and no new include directory. That claim is
//      only ever tested by compiling this file. Likewise the globally-set
//      EIGEN_MPL2_ONLY -- BFGSMat.h pulls <Eigen/LU>, and if a future Eigen
//      restricted that under MPL2 the failure shows up here first.
//   2. The solver actually solves. A truncated or corrupted header drop can
//      still compile.
//   3. The LBFGSParam defaults, which are NOT the values mmvec needs. mmvec
//      reproduces SciPy's L-BFGS-B stopping rule, and every default pinned
//      below has to be overridden to do it: m 6 -> 10 (SciPy maxcor),
//      epsilon_rel 1e-5 -> 0 (SciPy has no relative gradient test), and
//      past/delta 0/0 -> 1/1e-9 (SciPy's ftol). If upstream ever changes one
//      of these, our overrides silently start meaning something else, so the
//      defaults are asserted rather than trusted.
//
// Note that LBFGSParam::ftol and ::wolfe are the line-search Armijo and Wolfe
// constants -- they are NOT SciPy's convergence ftol, which maps to
// past/delta. See ext/LBFGSpp/PROVENANCE.md.
// ---------------------------------------------------------------------------
namespace {
// f(x) = sum_i (x_i - i)^2, minimized at x_i = i with f = 0. Strictly convex
// and exactly solvable, so any deviation is the solver's fault, not the
// problem's.
struct Quadratic {
	double operator()(const VectorXd &x, VectorXd &grad) {
		double f = 0.0;
		for (int i = 0; i < x.size(); ++i) {
			const double d = x[i] - static_cast<double>(i);
			f += d * d;
			grad[i] = 2.0 * d;
		}
		return f;
	}
};

// f(x) = sum(x^2), but the reported gradient carries the wrong sign, so every
// direction the solver believes is downhill actually goes uphill and the line
// search cannot satisfy Armijo at any step size.
struct SignFlippedGradient {
	double operator()(const VectorXd &x, VectorXd &grad) {
		grad = -2.0 * x;
		return x.squaredNorm();
	}
};

// ---------------------------------------------------------------------------
// The hand-worked pair the sufficient-statistics cases below are checked
// against. 3 samples; X has 2 features, Y has 3. Written dense here and given
// to the core as COO with the zeros dropped, which is how the SQL layer will
// hand it over.
//
//   X =  s0 [1, 2]        Y =  s0 [5, 1, 0]
//        s1 [3, 0]             s1 [0, 2, 7]
//        s2 [0, 4]             s2 [1, 0, 3]
//
//   X^T Y = [ 1*5+3*0+0*1,  1*1+3*2+0*0,  1*0+3*7+0*3 ]  =  [  5,  7, 21 ]
//           [ 2*5+0*0+4*1,  2*1+0*2+4*0,  2*0+0*7+4*3 ]     [ 14,  2, 12 ]
//
// No all-zero row or column in either table, so the pair is valid input.
// Every value is a small integer, so the expectations are exact, not
// approximate.
const std::vector<double> kHandExpectedYSums = {5, 7, 21, 14, 2, 12};
const std::vector<double> kHandExpectedNSums = {33, 28};

SparseCounts HandX() {
	SparseCounts x;
	x.n_rows = 3;
	x.n_cols = 2;
	x.row = {0, 0, 1, 2};
	x.col = {0, 1, 0, 1};
	x.val = {1, 2, 3, 4};
	return x;
}

SparseCounts HandY() {
	SparseCounts y;
	y.n_rows = 3;
	y.n_cols = 3;
	y.row = {0, 0, 1, 1, 2, 2};
	y.col = {0, 1, 1, 2, 0, 2};
	y.val = {5, 1, 2, 7, 1, 3};
	return y;
}

// Reverse the COO arrays in lockstep: same cells, opposite order.
SparseCounts Reversed(SparseCounts t) {
	std::reverse(t.row.begin(), t.row.end());
	std::reverse(t.col.begin(), t.col.end());
	std::reverse(t.val.begin(), t.val.end());
	return t;
}
} // namespace

TEST_CASE("vendored LBFGS++ minimizes a convex quadratic", "[mmvec][lbfgspp]") {
	LBFGSpp::LBFGSParam<double> param;
	param.m = 10;
	param.epsilon = 1e-5;
	param.epsilon_rel = 0.0;
	param.past = 1;
	param.delta = 1e-9;
	param.max_iterations = 200;

	LBFGSpp::LBFGSSolver<double, LBFGSpp::LineSearchMoreThuente> solver(param);
	Quadratic fun;
	VectorXd x = VectorXd::Zero(5);
	double fx = 0.0;
	const int niter = solver.minimize(fun, x, fx);

	REQUIRE(niter > 0);
	REQUIRE(fx == Approx(0.0).margin(1e-20));
	for (int i = 0; i < 5; ++i) {
		REQUIRE(x[i] == Approx(static_cast<double>(i)).margin(1e-10));
	}
}

TEST_CASE("LBFGSParam defaults are the ones mmvec must override", "[mmvec][lbfgspp]") {
	const LBFGSpp::LBFGSParam<double> d;

	// History size: SciPy's maxcor is 10, LBFGS++ defaults to 6.
	REQUIRE(d.m == 6);
	// SciPy tests only an absolute gradient tolerance; LBFGS++ adds a relative
	// one that must be switched off.
	REQUIRE(d.epsilon == Approx(1e-5));
	REQUIRE(d.epsilon_rel == Approx(1e-5));
	// The objective-change test is off by default; SciPy's ftol=1e-9 maps onto
	// past=1 with delta=1e-9.
	REQUIRE(d.past == 0);
	REQUIRE(d.delta == Approx(0.0));
	// Line-search constants. SciPy's dcsrch uses ftol=1e-3 and gtol=0.9, so
	// wolfe already agrees but ftol does not.
	REQUIRE(d.ftol == Approx(1e-4));
	REQUIRE(d.wolfe == Approx(0.9));
	REQUIRE(d.max_linesearch == 20);
}

TEST_CASE("LBFGS++ throws on line-search failure and abandons x", "[mmvec][lbfgspp]") {
	// SciPy's L-BFGS-B reports abnormal termination through a status code and
	// still hands back its last iterate; LBFGS++ throws instead. mmvec has to
	// survive non-convergence -- real data never reaches a stationary point --
	// so the fit keeps its own best-so-far (theta, f) snapshot rather than
	// reading x and fx back out of the solver. This case pins the behaviour
	// that makes the snapshot necessary: on failure the exception escapes and
	// x is left at the starting point, not at the best point seen.
	LBFGSpp::LBFGSParam<double> param;
	param.max_iterations = 50;
	LBFGSpp::LBFGSSolver<double, LBFGSpp::LineSearchMoreThuente> solver(param);
	SignFlippedGradient fun;
	VectorXd x = VectorXd::Constant(3, 1.0);
	double fx = 0.0;

	REQUIRE_THROWS_AS(solver.minimize(fun, x, fx), std::runtime_error);
	REQUIRE_THROWS_WITH(solver.minimize(fun, x, fx), ContainsSubstring("line search"));
	REQUIRE(x[0] == Approx(1.0));
}

// ---------------------------------------------------------------------------
// Sufficient statistics: y_sums = X^T Y and n_sums = rowsum(y_sums)
//
// These two arrays are the data's ENTIRE contribution to the objective -- the
// sample axis appears nowhere else in the algorithm. A wrong y_sums therefore
// does not produce an error or a failure to converge; it produces a fit that
// converges beautifully to the wrong answer. Hence exact hand-computed
// expectations rather than a self-consistency check.
// ---------------------------------------------------------------------------

TEST_CASE("sufficient statistics reproduce a hand-computed X^T Y", "[mmvec]") {
	const auto s = ComputeSufficientStats(HandX(), HandY());

	REQUIRE(s.n_features_x == 2);
	REQUIRE(s.n_features_y == 3);
	// Exact: every input is a small integer, so no rounding is involved.
	REQUIRE(s.y_sums == kHandExpectedYSums);
	REQUIRE(s.n_sums == kHandExpectedNSums);
}

TEST_CASE("n_sums agrees with its independent definition", "[mmvec]") {
	// n_sums is computed as rowsum(y_sums), but it is DEFINED in the model as
	// SUM_n X[n,i] * SUM_j Y[n,j] -- the total Y depth each X-feature is
	// conditioned on. The two are equal in exact arithmetic; this pins the
	// identity that lets one pass over the data serve both, so a future
	// "optimization" that breaks it fails here rather than silently biasing
	// every log-normalizer term.
	const auto s = ComputeSufficientStats(HandX(), HandY());

	// rowsum(Y) = [6, 9, 4]; n_sums[i] = SUM_n X[n,i] * rowsum(Y)[n].
	const std::vector<double> by_definition = {1 * 6 + 3 * 9 + 0 * 4, 2 * 6 + 0 * 9 + 4 * 4};
	REQUIRE(s.n_sums == by_definition);

	// And it really is the row sum of y_sums, cell for cell.
	for (int64_t i = 0; i < s.n_features_x; ++i) {
		const auto begin = s.y_sums.begin() + static_cast<std::ptrdiff_t>(i * s.n_features_y);
		const double row = std::accumulate(begin, begin + static_cast<std::ptrdiff_t>(s.n_features_y), 0.0);
		REQUIRE(s.n_sums[static_cast<size_t>(i)] == row);
	}
}

TEST_CASE("sufficient statistics are bit-identical under COO reordering", "[mmvec]") {
	// The SQL layer cannot promise a row order -- a parallel scan of a feature
	// table delivers rows in whatever order the threads finish. The core is
	// nonetheless required to be deterministic, which holds because each output
	// cell (i,j) receives exactly one term per sample, so the summation order is
	// fixed by ascending sample index alone.
	//
	// A naive implementation that accumulated in COO order would satisfy the
	// integer case above and still fail here, so the values are chosen so that
	// summation order is observable: cell (0,0) receives 1.0, 1e-16 and 1e-16
	// from samples 0, 1 and 2. Ascending gives (1 + 1e-16) + 1e-16 == 1.0
	// exactly, because 1e-16 is below half an ulp of 1.0 and is swallowed twice;
	// the reverse order sums the two small terms first and yields
	// 1.0000000000000002. Both are pinned below, so this case detects a
	// COO-order-dependent accumulation AND pins which order is canonical.
	SparseCounts x;
	x.n_rows = 3;
	x.n_cols = 2;
	x.row = {0, 0, 1, 1, 2, 2};
	x.col = {0, 1, 0, 1, 0, 1};
	x.val = {1.0, 2.0, 1e-16, 3.0, 1e-16, 4.0};

	SparseCounts y;
	y.n_rows = 3;
	y.n_cols = 2;
	y.row = {0, 0, 1, 1, 2, 2};
	y.col = {0, 1, 0, 1, 0, 1};
	y.val = {1.0, 5.0, 1.0, 6.0, 1.0, 7.0};

	const auto forward = ComputeSufficientStats(x, y);
	const auto reversed = ComputeSufficientStats(Reversed(x), Reversed(y));

	// Ascending-sample accumulation, exactly.
	REQUIRE(forward.y_sums[0] == 1.0);
	REQUIRE(forward.y_sums[0] != 1.0000000000000002);
	// Bit-identical, not merely close.
	REQUIRE(reversed.y_sums == forward.y_sums);
	REQUIRE(reversed.n_sums == forward.n_sums);
}

TEST_CASE("sufficient statistics reject structurally invalid tables", "[mmvec]") {
	SECTION("parallel arrays of unequal length") {
		auto x = HandX();
		x.val.pop_back();
		REQUIRE_THROWS_WITH(ComputeSufficientStats(x, HandY()), ContainsSubstring("same length"));
	}

	SECTION("sample counts differ between the tables") {
		auto y = HandY();
		y.n_rows = 4;
		REQUIRE_THROWS_WITH(ComputeSufficientStats(HandX(), y), ContainsSubstring("number of samples"));
	}

	SECTION("no samples at all") {
		SparseCounts x;
		x.n_rows = 0;
		x.n_cols = 2;
		SparseCounts y;
		y.n_rows = 0;
		y.n_cols = 3;
		REQUIRE_THROWS_WITH(ComputeSufficientStats(x, y), ContainsSubstring("at least one sample"));
	}

	SECTION("row index out of range") {
		auto x = HandX();
		x.row[0] = 3;
		REQUIRE_THROWS_WITH(ComputeSufficientStats(x, HandY()), ContainsSubstring("out of range"));
	}

	SECTION("column index out of range") {
		auto y = HandY();
		y.col[2] = 3;
		REQUIRE_THROWS_WITH(ComputeSufficientStats(HandX(), y), ContainsSubstring("out of range"));
	}

	SECTION("negative index") {
		auto x = HandX();
		x.col[1] = -1;
		REQUIRE_THROWS_WITH(ComputeSufficientStats(x, HandY()), ContainsSubstring("out of range"));
	}
}

TEST_CASE("sufficient statistics require a fittable model shape", "[mmvec]") {
	SECTION("a single Y-feature leaves nothing but the reference category") {
		// Y-feature 0 is pinned to a zero logit, so with d2 == 1 every
		// conditional probability is 1 by construction and there is no
		// likelihood to maximize. Rejecting this beats returning a "fit" whose
		// parameters are pure prior.
		SparseCounts y;
		y.n_rows = 3;
		y.n_cols = 1;
		y.row = {0, 1, 2};
		y.col = {0, 0, 0};
		y.val = {1, 1, 1};
		REQUIRE_THROWS_WITH(ComputeSufficientStats(HandX(), y), ContainsSubstring("at least two"));
	}

	SECTION("no X-features") {
		SparseCounts x;
		x.n_rows = 3;
		x.n_cols = 0;
		REQUIRE_THROWS_WITH(ComputeSufficientStats(x, HandY()), ContainsSubstring("at least one"));
	}
}

TEST_CASE("sufficient statistics reject values that are not counts", "[mmvec]") {
	SECTION("negative value") {
		// A CLR- or log-transformed table is a plausible mistake to make, and
		// MMvec's multinomial likelihood is meaningless on one. Refuse rather
		// than fit it.
		auto x = HandX();
		x.val[0] = -1.0;
		REQUIRE_THROWS_WITH(ComputeSufficientStats(x, HandY()), ContainsSubstring("negative"));
	}

	SECTION("NaN value") {
		auto y = HandY();
		y.val[3] = std::nan("");
		REQUIRE_THROWS_WITH(ComputeSufficientStats(HandX(), y), ContainsSubstring("finite"));
	}

	SECTION("infinite value") {
		auto x = HandX();
		x.val[2] = std::numeric_limits<double>::infinity();
		REQUIRE_THROWS_WITH(ComputeSufficientStats(x, HandY()), ContainsSubstring("finite"));
	}
}

TEST_CASE("sufficient statistics reject duplicate cells", "[mmvec]") {
	// A count split across two COO rows for the same (sample, feature) would
	// contribute x_a*y + x_b*y where the model wants (x_a+x_b)*y -- equal in
	// exact arithmetic, but it breaks the one-term-per-cell property that makes
	// the result independent of COO order. Reject rather than silently give up
	// determinism.
	auto x = HandX();
	x.row.push_back(0);
	x.col.push_back(0);
	x.val.push_back(1);
	REQUIRE_THROWS_WITH(ComputeSufficientStats(x, HandY()), ContainsSubstring("duplicate"));
}

TEST_CASE("sufficient statistics reject all-zero rows and columns", "[mmvec]") {
	// scikit-bio rejects these at fit time and so do we. An all-zero feature
	// column contributes nothing to the likelihood, so its embedding is
	// determined entirely by the prior: the fit converges and that feature's
	// output is meaningless. Failing loudly is the only honest option, and the
	// message has to tell the user what to do about it.
	SECTION("all-zero X column") {
		SparseCounts x;
		x.n_rows = 3;
		x.n_cols = 3; // feature 2 never appears
		x.row = {0, 0, 1, 2};
		x.col = {0, 1, 0, 1};
		x.val = {1, 2, 3, 4};
		REQUIRE_THROWS_WITH(ComputeSufficientStats(x, HandY()), ContainsSubstring("all-zero column"));
	}

	SECTION("all-zero Y column") {
		auto y = HandY();
		y.n_cols = 4; // feature 3 never appears
		REQUIRE_THROWS_WITH(ComputeSufficientStats(HandX(), y), ContainsSubstring("all-zero column"));
	}

	SECTION("all-zero X row") {
		SparseCounts x;
		x.n_rows = 3;
		x.n_cols = 2; // sample 1 never appears
		x.row = {0, 0, 2, 2};
		x.col = {0, 1, 0, 1};
		x.val = {1, 2, 3, 4};
		REQUIRE_THROWS_WITH(ComputeSufficientStats(x, HandY()), ContainsSubstring("all-zero row"));
	}

	SECTION("all-zero Y row") {
		SparseCounts y;
		y.n_rows = 3;
		y.n_cols = 3; // sample 2 never appears
		y.row = {0, 0, 1, 1};
		y.col = {0, 1, 1, 2};
		y.val = {5, 1, 2, 7};
		REQUIRE_THROWS_WITH(ComputeSufficientStats(HandX(), y), ContainsSubstring("all-zero row"));
	}

	SECTION("a stored zero does not rescue an otherwise-absent column") {
		// The check is on sums, not on index presence: a column that appears
		// only with explicit zeros is still all-zero.
		auto x = HandX();
		x.n_cols = 3;
		x.row.push_back(0);
		x.col.push_back(2);
		x.val.push_back(0.0);
		REQUIRE_THROWS_WITH(ComputeSufficientStats(x, HandY()), ContainsSubstring("all-zero column"));
	}
}

// ---------------------------------------------------------------------------
// Cross-check against the committed fixtures at real scale.
//
// The hand case above pins the formula on six cells. This one runs the sparse
// outer product against a straightforward dense triple loop over the actual
// oracle inputs -- up to 50 samples x 8 X-features x 10 Y-features -- so a bug
// that only shows up with many samples, sparse gaps, or a feature absent from
// some samples has somewhere to fail.
//
// Bit-identity is the right assertion, not approximate equality: the dense loop
// accumulates over samples in ascending order, exactly as the sparse one does,
// and skipping a zero cell is indistinguishable from adding 0.0 to a
// non-negative running sum.
// ---------------------------------------------------------------------------
namespace {
SparseCounts ToCoo(const std::vector<double> &dense, int64_t n_rows, int64_t n_cols) {
	SparseCounts t;
	t.n_rows = n_rows;
	t.n_cols = n_cols;
	for (int64_t r = 0; r < n_rows; ++r) {
		for (int64_t c = 0; c < n_cols; ++c) {
			const double v = dense[static_cast<size_t>(r * n_cols + c)];
			if (v != 0.0) {
				t.row.push_back(r);
				t.col.push_back(c);
				t.val.push_back(v);
			}
		}
	}
	return t;
}

// Textbook X^T Y over dense row-major matrices, summing samples in order.
std::vector<double> DenseXtY(const std::vector<double> &x, const std::vector<double> &y, int64_t n_samples, int64_t d1,
                             int64_t d2) {
	std::vector<double> out(static_cast<size_t>(d1 * d2), 0.0);
	for (int64_t i = 0; i < d1; ++i) {
		for (int64_t j = 0; j < d2; ++j) {
			double acc = 0.0;
			for (int64_t n = 0; n < n_samples; ++n) {
				acc += x[static_cast<size_t>(n * d1 + i)] * y[static_cast<size_t>(n * d2 + j)];
			}
			out[static_cast<size_t>(i * d2 + j)] = acc;
		}
	}
	return out;
}

void CheckAgainstDense(const std::vector<double> &xd, const std::vector<double> &yd, int64_t n_samples, int64_t d1,
                       int64_t d2) {
	REQUIRE(xd.size() == static_cast<size_t>(n_samples * d1));
	REQUIRE(yd.size() == static_cast<size_t>(n_samples * d2));

	const auto s = ComputeSufficientStats(ToCoo(xd, n_samples, d1), ToCoo(yd, n_samples, d2));
	REQUIRE(s.n_features_x == d1);
	REQUIRE(s.n_features_y == d2);
	REQUIRE(s.y_sums == DenseXtY(xd, yd, n_samples, d1, d2));

	// The multinomial depth each X-feature is conditioned on, computed the
	// other way round: SUM_n X[n,i] * (total Y counts in sample n).
	std::vector<double> depth(static_cast<size_t>(d1), 0.0);
	for (int64_t n = 0; n < n_samples; ++n) {
		double y_total = 0.0;
		for (int64_t j = 0; j < d2; ++j) {
			y_total += yd[static_cast<size_t>(n * d2 + j)];
		}
		for (int64_t i = 0; i < d1; ++i) {
			depth[static_cast<size_t>(i)] += xd[static_cast<size_t>(n * d1 + i)] * y_total;
		}
	}
	for (size_t i = 0; i < depth.size(); ++i) {
		REQUIRE(s.n_sums[i] == Approx(depth[i]).epsilon(1e-14));
	}
}
} // namespace

TEST_CASE("sufficient statistics match a dense X^T Y on the oracle fixtures", "[mmvec]") {
	SECTION("toy (5 x 8 / 5 x 6)") {
		using namespace miint::mmvec_oracle::toy;
		CheckAgainstDense(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY);
	}
	SECTION("synth_a (8 x 4 / 8 x 6)") {
		using namespace miint::mmvec_oracle::synth_a;
		CheckAgainstDense(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY);
	}
	SECTION("synth_b (50 x 8 / 50 x 10)") {
		using namespace miint::mmvec_oracle::synth_b;
		CheckAgainstDense(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY);
	}
}

// ---------------------------------------------------------------------------
// The objective: negative log-posterior and its gradient
//
// T1 is the tightest tier in the oracle because it is optimizer-free: given the
// committed theta, any correct implementation must land on the same loss,
// gradient and logits up to floating-point summation order. It is therefore the
// one place a sign error, a mishandled reference category, a missing prior term
// or a transposed gradient block cannot hide.
//
// The finite-difference cases below matter for a second reason: scikit-bio
// grad-checks only its Adam minibatch objective, so the full-batch gradient --
// the one L-BFGS actually descends, and the default path -- has no gradient test
// upstream at all. Both are checked here.
// ---------------------------------------------------------------------------
namespace {
using miint::mmvec::BatchNorm;
using miint::mmvec::BatchNormFactor;
using miint::mmvec::ComputeLogits;
using miint::mmvec::ComputeSufficientStats;
using miint::mmvec::FullBatchLossAndGradient;
using miint::mmvec::MinibatchInputs;
using miint::mmvec::MinibatchLossAndGradient;
using miint::mmvec::ModelShape;
using miint::mmvec::NumParams;
using miint::mmvec::Priors;
using miint::mmvec::SufficientStats;
using miint::mmvec::Workspace;

// The oracle's priors: scikit-bio's defaults, and the ones every carved value
// above was produced under.
Priors OraclePriors() {
	return Priors {0.0, 1.0, 0.0, 1.0};
}

struct FixtureCase {
	ModelShape shape;
	SufficientStats stats;
	MinibatchInputs mb;
};

// Build both the summed-away statistics and the minibatch inputs from a case's
// carved dense matrices, so the objective sees exactly the oracle's data.
FixtureCase MakeCase(const std::vector<double> &xd, const std::vector<double> &yd, int64_t n_samples, int64_t d1,
                     int64_t d2, int32_t p) {
	FixtureCase c;
	c.shape = ModelShape {d1, d2, p};
	const auto x_coo = ToCoo(xd, n_samples, d1);
	c.stats = ComputeSufficientStats(x_coo, ToCoo(yd, n_samples, d2));
	c.mb.n_samples = n_samples;
	c.mb.x_row = x_coo.row;
	c.mb.x_col = x_coo.col;
	c.mb.x_val = x_coo.val;
	c.mb.y_dense = yd;
	return c;
}

// Central differences: (f(t+h) - f(t-h)) / 2h, with h scaled to each parameter
// so the check behaves the same for a parameter near 0 as for one near 5.
void CheckGradientByFiniteDifference(const std::vector<double> &theta,
                                     const std::function<double(const std::vector<double> &)> &loss,
                                     const std::vector<double> &analytic, double rel_tol) {
	REQUIRE(analytic.size() == theta.size());
	double worst = 0.0;
	for (size_t i = 0; i < theta.size(); ++i) {
		const double h = 1e-6 * std::max(1.0, std::abs(theta[i]));
		auto up = theta;
		auto down = theta;
		up[i] += h;
		down[i] -= h;
		const double numeric = (loss(up) - loss(down)) / (2.0 * h);
		const double denom = std::max({std::abs(numeric), std::abs(analytic[i]), 1.0});
		worst = std::max(worst, std::abs(numeric - analytic[i]) / denom);
	}
	INFO("worst relative gradient error " << worst);
	REQUIRE(worst < rel_tol);
}
} // namespace

TEST_CASE("full-batch objective reproduces the T1 oracle", "[mmvec]") {
	using miint::mmvec_oracle::kT1LossRelTol;
	using miint::mmvec_oracle::kT1Tol;

	auto check = [&](const std::vector<double> &xd, const std::vector<double> &yd, int64_t n_samples, int64_t d1,
	                 int64_t d2, int32_t p, const std::vector<double> &theta, double want_loss,
	                 const std::vector<double> &want_grad, const std::vector<double> &want_logits) {
		auto c = MakeCase(xd, yd, n_samples, d1, d2, p);
		REQUIRE(theta.size() == static_cast<size_t>(NumParams(c.shape)));

		Workspace ws;
		std::vector<double> grad;
		const double loss = FullBatchLossAndGradient(c.shape, OraclePriors(), c.stats, theta, ws, grad);

		REQUIRE(loss == Approx(want_loss).epsilon(kT1LossRelTol));
		REQUIRE(grad.size() == want_grad.size());
		for (size_t i = 0; i < want_grad.size(); ++i) {
			INFO("gradient element " << i);
			REQUIRE(grad[i] == Approx(want_grad[i]).margin(kT1Tol));
		}

		const auto logits = ComputeLogits(c.shape, theta);
		REQUIRE(logits.size() == want_logits.size());
		for (size_t i = 0; i < want_logits.size(); ++i) {
			INFO("logit element " << i);
			REQUIRE(logits[i] == Approx(want_logits[i]).margin(kT1Tol));
		}
		// The reference category is pinned, not merely small.
		for (int64_t i = 0; i < d1; ++i) {
			REQUIRE(logits[static_cast<size_t>(i * d2)] == 0.0);
		}
	};

	SECTION("toy") {
		using namespace miint::mmvec_oracle::toy;
		check(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents, kTheta, kT1Loss, kT1Grad,
		      kT1Logits);
	}
	SECTION("synth_a") {
		using namespace miint::mmvec_oracle::synth_a;
		check(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents, kTheta, kT1Loss, kT1Grad,
		      kT1Logits);
	}
	SECTION("synth_b") {
		using namespace miint::mmvec_oracle::synth_b;
		check(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents, kTheta, kT1Loss, kT1Grad,
		      kT1Logits);
	}
}

TEST_CASE("full-batch gradient matches central finite differences", "[mmvec]") {
	// This is the gradient L-BFGS descends, and scikit-bio has no test for it.
	// Every element of theta is perturbed, so a block that is transposed,
	// mis-offset in the packed layout, or missing its prior term shows up here
	// even where T1's single carved theta might not distinguish it.
	using namespace miint::mmvec_oracle::toy;
	auto c = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
	const Priors priors {0.25, 2.0, -0.5, 0.75}; // deliberately not the defaults

	Workspace ws;
	std::vector<double> grad;
	FullBatchLossAndGradient(c.shape, priors, c.stats, kTheta, ws, grad);

	auto loss_only = [&](const std::vector<double> &t) {
		Workspace w;
		std::vector<double> g;
		return FullBatchLossAndGradient(c.shape, priors, c.stats, t, w, g);
	};
	CheckGradientByFiniteDifference(kTheta, loss_only, grad, 1e-6);
}

TEST_CASE("minibatch gradient matches central finite differences", "[mmvec]") {
	// The batch deliberately repeats X features (index 17 appears many times in
	// the committed oracle sequence too), so the scatter-accumulation path is
	// exercised rather than a one-hit-per-feature special case.
	using namespace miint::mmvec_oracle::toy;
	auto c = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
	const std::vector<int64_t> batch = {0, 3, 3, 7, 11, 3, 0, 19};
	const Priors priors {0.25, 2.0, -0.5, 0.75};
	const double norm = 2.5;

	Workspace ws;
	std::vector<double> grad;
	MinibatchLossAndGradient(c.shape, priors, c.mb, batch, norm, kTheta, ws, grad);

	auto loss_only = [&](const std::vector<double> &t) {
		Workspace w;
		std::vector<double> g;
		return MinibatchLossAndGradient(c.shape, priors, c.mb, batch, norm, t, w, g);
	};
	CheckGradientByFiniteDifference(kTheta, loss_only, grad, 1e-6);
}

TEST_CASE("minibatch objective reduces to the full-batch one", "[mmvec]") {
	// The two paths are the same likelihood counted differently: the full-batch
	// objective weights each (sample, X-feature) pair by its X count, while the
	// minibatch treats each sampled nonzero as one draw of the sample's whole Y
	// vector. So with every X count equal to 1 and the "batch" being every
	// nonzero exactly once at norm = 1, they must agree -- which ties the
	// summed-away statistics to the per-draw formulation and would catch a
	// weighting or reference-column error in either path.
	//
	// Not bit-exact: the full-batch path reduces over d1 rows via a matrix
	// product while the minibatch accumulates in batch order.
	const int64_t n_samples = 4;
	const int64_t d1 = 3;
	const int64_t d2 = 4;
	const int32_t p = 2;
	// X is 0/1 only. Every row and column has a 1, so the pair is valid.
	const std::vector<double> xd = {1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1};
	const std::vector<double> yd = {3, 1, 0, 2, 0, 4, 1, 1, 2, 0, 5, 0, 1, 2, 0, 3};

	auto c = MakeCase(xd, yd, n_samples, d1, d2, p);
	const Priors priors {0.1, 1.5, -0.2, 0.9};
	std::vector<double> theta(static_cast<size_t>(NumParams(c.shape)));
	for (size_t i = 0; i < theta.size(); ++i) {
		theta[i] = 0.3 * static_cast<double>(i % 7) - 0.8;
	}

	std::vector<int64_t> all_nonzeros(c.mb.x_row.size());
	std::iota(all_nonzeros.begin(), all_nonzeros.end(), 0);

	Workspace ws_fb;
	Workspace ws_mb;
	std::vector<double> grad_fb;
	std::vector<double> grad_mb;
	const double loss_fb = FullBatchLossAndGradient(c.shape, priors, c.stats, theta, ws_fb, grad_fb);
	const double loss_mb = MinibatchLossAndGradient(c.shape, priors, c.mb, all_nonzeros, 1.0, theta, ws_mb, grad_mb);

	REQUIRE(loss_mb == Approx(loss_fb).epsilon(1e-12));
	REQUIRE(grad_mb.size() == grad_fb.size());
	for (size_t i = 0; i < grad_fb.size(); ++i) {
		INFO("gradient element " << i);
		REQUIRE(grad_mb[i] == Approx(grad_fb[i]).epsilon(1e-11).margin(1e-11));
	}
}

TEST_CASE("batch_norm scales the likelihood but never the prior", "[mmvec]") {
	// 'legacy' is not a different step size for the same posterior -- it
	// multiplies the likelihood by n_samples/sum(X) while leaving the priors
	// alone, which is a differently-regularized model. This pins that: the two
	// losses must be an exact affine function of norm with a common intercept,
	// and that intercept is the unscaled prior.
	using namespace miint::mmvec_oracle::toy;
	auto c = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
	const std::vector<int64_t> batch = {1, 4, 9, 12};
	const Priors priors {0.0, 1.0, 0.0, 1.0};

	double x_total = 0.0;
	for (const double v : kXCounts) {
		x_total += v;
	}
	const auto batch_size = static_cast<int64_t>(batch.size());
	const double norm_u = BatchNormFactor(BatchNorm::Unbiased, x_total, kNSamples, batch_size);
	const double norm_l = BatchNormFactor(BatchNorm::Legacy, x_total, kNSamples, batch_size);

	// The documented formulas, not just "some number".
	REQUIRE(norm_u == Approx(x_total / static_cast<double>(batch_size)));
	REQUIRE(norm_l == Approx(static_cast<double>(kNSamples) / static_cast<double>(batch_size)));
	REQUIRE(norm_u != Approx(norm_l));

	Workspace ws;
	std::vector<double> grad;
	const double loss_u = MinibatchLossAndGradient(c.shape, priors, c.mb, batch, norm_u, kTheta, ws, grad);
	const double loss_l = MinibatchLossAndGradient(c.shape, priors, c.mb, batch, norm_l, kTheta, ws, grad);

	// loss(norm) = norm*data + prior, so data is recoverable from the two points
	// and both intercepts must agree.
	const double data = (loss_u - loss_l) / (norm_u - norm_l);
	const double intercept_u = loss_u - norm_u * data;
	const double intercept_l = loss_l - norm_l * data;
	REQUIRE(intercept_u == Approx(intercept_l).epsilon(1e-9));

	// And that intercept is exactly the Gaussian prior at these parameters.
	double prior = 0.0;
	for (const double t : kTheta) {
		prior += 0.5 * t * t;
	}
	REQUIRE(intercept_u == Approx(prior).epsilon(1e-9));
}

TEST_CASE("the Gaussian prior uses the right mean and scale per block", "[mmvec]") {
	// Four distinct prior values, so applying the X scale to a Y block (or
	// swapping mean and scale) cannot pass. The data term is independent of the
	// priors, so driving the scales to a huge value isolates it, and the
	// difference from the real priors must equal the closed form.
	using namespace miint::mmvec_oracle::toy;
	auto c = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
	const Priors priors {0.5, 2.0, -0.25, 0.5};
	const Priors negligible {0.5, 1e8, -0.25, 1e8};

	Workspace ws;
	std::vector<double> grad;
	const double with_prior = FullBatchLossAndGradient(c.shape, priors, c.stats, kTheta, ws, grad);
	const double data_only = FullBatchLossAndGradient(c.shape, negligible, c.stats, kTheta, ws, grad);

	// Blocks in packed order: x_main (d1*p), x_bias (d1), y_main (p*(d2-1)), y_bias (d2-1).
	const size_t n_x = static_cast<size_t>(kNFeaturesX * kNComponents + kNFeaturesX);
	double expected = 0.0;
	for (size_t i = 0; i < kTheta.size(); ++i) {
		const bool x_side = i < n_x;
		const double mean = x_side ? priors.x_prior_mean : priors.y_prior_mean;
		const double scale = x_side ? priors.x_prior_scale : priors.y_prior_scale;
		const double d = kTheta[i] - mean;
		expected += 0.5 * d * d / (scale * scale);
	}
	REQUIRE(with_prior - data_only == Approx(expected).epsilon(1e-9));
}

TEST_CASE("the objective rejects inconsistent shapes and parameters", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;
	auto c = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
	Workspace ws;
	std::vector<double> grad;

	SECTION("theta of the wrong length") {
		auto short_theta = kTheta;
		short_theta.pop_back();
		REQUIRE_THROWS_WITH(FullBatchLossAndGradient(c.shape, OraclePriors(), c.stats, short_theta, ws, grad),
		                    ContainsSubstring("parameters, expected"));
	}

	SECTION("n_components below one") {
		ModelShape bad = c.shape;
		bad.n_components = 0;
		REQUIRE_THROWS_WITH(FullBatchLossAndGradient(bad, OraclePriors(), c.stats, kTheta, ws, grad),
		                    ContainsSubstring("n_components must be >= 1"));
	}

	SECTION("non-positive prior scale") {
		// scikit-bio validates this at fit time; a zero scale is not a prior.
		REQUIRE_THROWS_WITH(FullBatchLossAndGradient(c.shape, Priors {0.0, 0.0, 0.0, 1.0}, c.stats, kTheta, ws, grad),
		                    ContainsSubstring("x_prior_scale must be > 0"));
		REQUIRE_THROWS_WITH(FullBatchLossAndGradient(c.shape, Priors {0.0, 1.0, 0.0, -1.0}, c.stats, kTheta, ws, grad),
		                    ContainsSubstring("y_prior_scale must be > 0"));
	}

	SECTION("statistics that do not match the shape") {
		ModelShape bad = c.shape;
		bad.n_features_x = kNFeaturesX + 1;
		std::vector<double> theta(static_cast<size_t>(NumParams(bad)), 0.1);
		REQUIRE_THROWS_WITH(FullBatchLossAndGradient(bad, OraclePriors(), c.stats, theta, ws, grad),
		                    ContainsSubstring("but the model shape is"));
	}

	SECTION("an empty minibatch") {
		REQUIRE_THROWS_WITH(MinibatchLossAndGradient(c.shape, OraclePriors(), c.mb, {}, 1.0, kTheta, ws, grad),
		                    ContainsSubstring("minibatch is empty"));
	}

	SECTION("a minibatch index past the last nonzero") {
		const auto nnz = static_cast<int64_t>(c.mb.x_row.size());
		REQUIRE_THROWS_WITH(MinibatchLossAndGradient(c.shape, OraclePriors(), c.mb, {nnz}, 1.0, kTheta, ws, grad),
		                    ContainsSubstring("out of range"));
	}

	SECTION("a batch size of zero has no normalization factor") {
		REQUIRE_THROWS_WITH(BatchNormFactor(BatchNorm::Unbiased, 100.0, 5, 0),
		                    ContainsSubstring("batch_size must be >= 1"));
	}
}

TEST_CASE("the log-normalizer survives extreme logits", "[mmvec]") {
	// The shift in log_norm = m + log(exp(-m) + SUM exp(logits - m)) is an
	// algebraic identity for any m, so clamping m at zero cannot be caught by
	// comparing values on ordinary data -- it buys numerical range only. That
	// makes overflow the only observable behaviour, and it is what this case
	// pins.
	//
	// theta is built so every non-reference logit equals a single constant c:
	// x_main and y_main are zero, y_bias is zero, and x_bias is c.
	using namespace miint::mmvec_oracle::toy;
	auto c = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
	const auto layout_n = static_cast<size_t>(NumParams(c.shape));
	const auto n_x_main = static_cast<size_t>(kNFeaturesX * kNComponents);
	// Priors would otherwise contribute 0.5*c^2 per bias and swamp the data term.
	const Priors no_prior {0.0, 1e12, 0.0, 1e12};
	const int64_t nref = kNFeaturesY - 1;

	// Total observed non-reference counts, and the grouped Y totals.
	double y_nonref_total = 0.0;
	for (int64_t i = 0; i < kNFeaturesX; ++i) {
		for (int64_t j = 1; j < kNFeaturesY; ++j) {
			y_nonref_total += c.stats.y_sums[static_cast<size_t>(i * kNFeaturesY + j)];
		}
	}
	double n_sums_total = 0.0;
	for (const double v : c.stats.n_sums) {
		n_sums_total += v;
	}

	auto loss_at = [&](double bias) {
		std::vector<double> theta(layout_n, 0.0);
		for (size_t i = 0; i < static_cast<size_t>(kNFeaturesX); ++i) {
			theta[n_x_main + i] = bias;
		}
		Workspace ws;
		std::vector<double> grad;
		const double loss = FullBatchLossAndGradient(c.shape, no_prior, c.stats, theta, ws, grad);
		for (const double g : grad) {
			REQUIRE(std::isfinite(g));
		}
		return loss;
	};

	SECTION("logits far below zero: the reference category takes all the mass") {
		// exp(-900) underflows to 0, so log_norm is exactly log(1) = 0 and the
		// loss reduces to -c * (total non-reference counts).
		const double loss = loss_at(-900.0);
		REQUIRE(std::isfinite(loss));
		REQUIRE(loss == Approx(900.0 * y_nonref_total).epsilon(1e-12));
	}

	SECTION("logits far above zero: the reference category is negligible") {
		// log_norm = 900 + log(nref) since every shifted exponential is exactly 1
		// and the reference's exp(-900) underflows away.
		const double loss = loss_at(900.0);
		REQUIRE(std::isfinite(loss));
		const double expected = n_sums_total * (900.0 + std::log(static_cast<double>(nref))) - 900.0 * y_nonref_total;
		REQUIRE(loss == Approx(expected).epsilon(1e-12));
	}
}

// ---------------------------------------------------------------------------
// L-BFGS
//
// T2 is the converged tier, and it is deliberately loose: duckdb-miint matches
// scikit-bio's hard-coded ftol=1e-9 / gtol=1e-5 rather than converging tighter,
// and under that rule the stopping point depends on the starting point. The
// optimum itself is unique to ~1e-7, so what is being checked is "the same
// optimum, found by a different L-BFGS implementation", not "the same iterates".
//
// Measured relative deviation of the final loss from the carved oracle:
// toy 1.1e-10, synth_a 3.0e-8, synth_b 2.6e-7. The band below is 1e-6.
// Evaluation counts are NOT compared -- ours were 70/125/206 against the
// oracle's 69/113/211, because the line searches differ.
// ---------------------------------------------------------------------------
namespace {
using miint::mmvec::FitLbfgsFromInit;
using miint::mmvec::kLbfgsGtol;
using miint::mmvec::Model;

constexpr double kT2LossRelBand = 1e-6;
} // namespace

TEST_CASE("L-BFGS reaches the T2 optimum from the committed init", "[mmvec]") {
	auto check = [&](const std::vector<double> &xd, const std::vector<double> &yd, int64_t n_samples, int64_t d1,
	                 int64_t d2, int32_t p, const std::vector<double> &theta, double t1_loss, int max_iter,
	                 double t2_loss) {
		auto c = MakeCase(xd, yd, n_samples, d1, d2, p);
		const Model m = FitLbfgsFromInit(c.shape, OraclePriors(), c.stats, theta, max_iter);

		// The first evaluation happens at theta0, which is the T1 oracle's theta --
		// so the fit and the fixed-theta tier must agree exactly at that point.
		// This ties the optimizer to the objective the oracle was built from.
		REQUIRE_FALSE(m.loss_curve.empty());
		REQUIRE(m.loss_curve.front() == Approx(t1_loss).epsilon(1e-12));

		REQUIRE(m.final_loss == Approx(t2_loss).epsilon(kT2LossRelBand));
		REQUIRE(m.final_loss < m.loss_curve.front()); // the fit actually improved
		REQUIRE(m.theta.size() == static_cast<size_t>(NumParams(c.shape)));
		REQUIRE(m.n_evals == static_cast<int64_t>(m.loss_curve.size()));

		// theta is the best point over all evaluations, so the reported loss must
		// be the minimum of the trace -- not merely the last entry, which for a
		// line-search probe can be higher.
		const double best = *std::min_element(m.loss_curve.begin(), m.loss_curve.end());
		REQUIRE(m.final_loss == Approx(best).epsilon(1e-12));
	};

	SECTION("toy") {
		using namespace miint::mmvec_oracle::toy;
		check(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents, kTheta, kT1Loss, kT2MaxIter,
		      kT2FinalLoss);
	}
	SECTION("synth_a") {
		using namespace miint::mmvec_oracle::synth_a;
		check(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents, kTheta, kT1Loss, kT2MaxIter,
		      kT2FinalLoss);
	}
	SECTION("synth_b") {
		using namespace miint::mmvec_oracle::synth_b;
		check(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents, kTheta, kT1Loss, kT2MaxIter,
		      kT2FinalLoss);
	}
}

TEST_CASE("a converged L-BFGS fit is not necessarily a stationary point", "[mmvec]") {
	// scikit-bio's ftol=1e-9 fires when the objective stops changing, which
	// happens well before the gradient vanishes. `toy` stops legitimately at
	// iteration 65 of 100 with max|gradient| ~ 1e-3 -- four orders of magnitude
	// above gtol. That is exactly why max_abs_grad is reported next to
	// `converged` instead of being folded into it, and this case pins the
	// distinction so nobody later "simplifies" converged to mean max|g| <= gtol.
	using namespace miint::mmvec_oracle::toy;
	auto c = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
	const Model m = FitLbfgsFromInit(c.shape, OraclePriors(), c.stats, kTheta, kT2MaxIter);

	REQUIRE(m.converged);
	REQUIRE(m.n_iter < kT2MaxIter);
	REQUIRE(m.max_abs_grad > kLbfgsGtol);
	REQUIRE_THAT(m.message, ContainsSubstring("converged"));
	REQUIRE_THAT(m.message, ContainsSubstring("max|gradient|"));
}

TEST_CASE("L-BFGS reports failure to converge instead of hiding it", "[mmvec]") {
	// scikit-bio neither raises nor warns when the optimizer gives up, so a caller
	// cannot tell a converged fit from an exhausted one. Both non-convergence
	// routes are checked here, and in each the model must still be usable: the
	// best point seen is returned rather than nothing.
	SECTION("a starved iteration budget") {
		using namespace miint::mmvec_oracle::toy;
		auto c = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
		const Model m = FitLbfgsFromInit(c.shape, OraclePriors(), c.stats, kTheta, 1);

		REQUIRE_FALSE(m.converged);
		REQUIRE(m.n_iter == 1);
		REQUIRE_THAT(m.message, ContainsSubstring("iteration limit"));
		REQUIRE(m.theta.size() == static_cast<size_t>(NumParams(c.shape)));
		const double best = *std::min_element(m.loss_curve.begin(), m.loss_curve.end());
		REQUIRE(m.final_loss == Approx(best).epsilon(1e-12));
	}

	SECTION("a case that exhausts its budget at scikit-bio's own max_iter") {
		// Measured, not chosen: at scikit-bio's own settings synth_b runs the full
		// 200 iterations without satisfying a convergence test, and SciPy likewise
		// terminates abnormally on it. If this ever starts converging, that is a
		// real behavioural change worth looking at rather than retuning.
		using namespace miint::mmvec_oracle::synth_b;
		auto c = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
		const Model m = FitLbfgsFromInit(c.shape, OraclePriors(), c.stats, kTheta, kT2MaxIter);

		REQUIRE_FALSE(m.converged);
		REQUIRE(m.n_iter == kT2MaxIter);
		REQUIRE_THAT(m.message, ContainsSubstring("iteration limit"));
		// Still lands on the right optimum despite not converging.
		REQUIRE(m.final_loss == Approx(kT2FinalLoss).epsilon(kT2LossRelBand));
	}
}

TEST_CASE("L-BFGS is deterministic", "[mmvec]") {
	// The fit takes its starting point as an argument and touches no RNG, so two
	// runs must agree bit for bit -- not merely to a tolerance. This is what lets
	// the SQL layer promise reproducible output for a given seed.
	using namespace miint::mmvec_oracle::synth_a;
	auto c = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
	const Model a = FitLbfgsFromInit(c.shape, OraclePriors(), c.stats, kTheta, kT2MaxIter);
	const Model b = FitLbfgsFromInit(c.shape, OraclePriors(), c.stats, kTheta, kT2MaxIter);

	REQUIRE(a.theta == b.theta);
	REQUIRE(a.loss_curve == b.loss_curve);
	REQUIRE(a.final_loss == b.final_loss);
	REQUIRE(a.n_iter == b.n_iter);
	REQUIRE(a.converged == b.converged);
}

TEST_CASE("L-BFGS rejects a nonsensical iteration budget", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;
	auto c = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
	REQUIRE_THROWS_WITH(FitLbfgsFromInit(c.shape, OraclePriors(), c.stats, kTheta, 0),
	                    ContainsSubstring("max_iter must be >= 1"));
	// And it inherits the objective's own validation rather than duplicating it.
	auto short_theta = kTheta;
	short_theta.pop_back();
	REQUIRE_THROWS_WITH(FitLbfgsFromInit(c.shape, OraclePriors(), c.stats, short_theta, 10),
	                    ContainsSubstring("parameters, expected"));
}

TEST_CASE("a failed line search still yields a usable model", "[mmvec]") {
	// This is the only situation in which the best-so-far snapshot changes the
	// answer, and it is why the snapshot exists. Measured: on every well-behaved
	// fixture run the snapshot and the solver's own final iterate are the SAME
	// point, so nothing is gained there -- but LBFGS++ signals line-search failure
	// by throwing, and (pinned separately above) leaves its iterate at the
	// starting point when it does. Without the snapshot the fit would propagate an
	// exception instead of returning the best parameters it found.
	//
	// A starting point at 1e150 reaches that failure deterministically: the
	// objective is enormous, and the line search shrinks its step below the
	// minimum allowed before making progress.
	using namespace miint::mmvec_oracle::toy;
	auto c = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
	const std::vector<double> wild(static_cast<size_t>(NumParams(c.shape)), 1e150);

	const Model m = FitLbfgsFromInit(c.shape, OraclePriors(), c.stats, wild, 100);

	REQUIRE_FALSE(m.converged);
	REQUIRE_THAT(m.message, ContainsSubstring("line search failed"));
	// The model is still complete and finite, not a half-written wreck.
	REQUIRE(m.theta.size() == static_cast<size_t>(NumParams(c.shape)));
	for (const double t : m.theta) {
		REQUIRE(std::isfinite(t));
	}
	REQUIRE(std::isfinite(m.final_loss));
	REQUIRE_FALSE(m.loss_curve.empty());
	const double best = *std::min_element(m.loss_curve.begin(), m.loss_curve.end());
	REQUIRE(m.final_loss == Approx(best).epsilon(1e-12));
}

TEST_CASE("the objective rejects a non-finite starting point", "[mmvec]") {
	// A NaN in theta would otherwise propagate silently: every loss and gradient
	// becomes NaN, no comparison fails, and the fit "succeeds" with a model made
	// entirely of NaN.
	using namespace miint::mmvec_oracle::toy;
	auto c = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
	auto bad = kTheta;
	bad[7] = std::nan("");
	REQUIRE_THROWS_WITH(FitLbfgsFromInit(c.shape, OraclePriors(), c.stats, bad, 10),
	                    ContainsSubstring("is not finite"));
	bad[7] = std::numeric_limits<double>::infinity();
	REQUIRE_THROWS_WITH(FitLbfgsFromInit(c.shape, OraclePriors(), c.stats, bad, 10),
	                    ContainsSubstring("is not finite"));
}

// ---------------------------------------------------------------------------
// Portable RNG
//
// The bit stream is std::mt19937_64, whose output the standard fixes exactly, but
// every transform from bits to doubles is written by hand. <random>'s
// distributions are implementation-defined, so relying on them would make a
// seeded fit produce different numbers on libstdc++, libc++ and Emscripten --
// which is exactly the caveat cluster_kmeans.hpp carries and which mmvec is
// meant not to have.
// ---------------------------------------------------------------------------
namespace {
using miint::mmvec::AdamParams;
using miint::mmvec::FitAdam;
using miint::mmvec::FitAdamWithIndices;
using miint::mmvec::InitTheta;
using miint::mmvec::Rng;
using miint::mmvec::SampleWeightedIndices;
} // namespace

TEST_CASE("the RNG bit stream is the standard mt19937_64", "[mmvec]") {
	// Non-circular: 9981545732273789042 is the standard's own published value for
	// the 10000th output of a default-seeded (5489) mt19937_64, so this pins both
	// that the engine is that engine and that Uniform() is the documented
	// top-53-bits mapping. If either changed, every seeded fit would change.
	Rng rng(5489);
	double last = 0.0;
	for (int i = 0; i < 10000; ++i) {
		last = rng.Uniform();
	}
	const auto expected = static_cast<double>(9981545732273789042ULL >> 11) * 0x1.0p-53;
	REQUIRE(last == expected);
}

TEST_CASE("the RNG draws are pinned for portability", "[mmvec]") {
	// These literals are golden values produced by this implementation, not by an
	// external oracle -- there is no external oracle, because scikit-bio draws from
	// numpy's PCG64 stream and cannot be matched. Their job is to fail if a
	// toolchain, a standard library, or a WASM build ever changes the sequence.
	const std::vector<double> want_uniform = {
	    0.15979336337046079, 0.99214520962982877, 0.039569025844865657,
	    0.59749466269467166, 0.54228496999260445, 0.057159791465323573,
	};
	const std::vector<double> want_normal = {
	    -0.48132337199836744, 0.10191855551453786, 0.064987953338865465,
	    -0.68060303256354293, 1.8863239328876753,  -1.0961189116175776,
	};

	Rng u(0);
	for (const double want : want_uniform) {
		REQUIRE(u.Uniform() == want);
	}
	Rng n(0);
	for (const double want : want_normal) {
		REQUIRE(n.Normal() == want);
	}
}

TEST_CASE("the RNG produces the distributions it claims", "[mmvec]") {
	// The golden literals above would pass just as happily if Normal() returned
	// uniforms, so the distributions are also checked on their own terms.
	SECTION("Uniform stays in range and has the right mean") {
		Rng rng(7);
		const int n = 200000;
		double sum = 0.0;
		bool all_in_range = true;
		for (int i = 0; i < n; ++i) {
			const double x = rng.Uniform();
			all_in_range = all_in_range && x >= 0.0 && x < 1.0;
			sum += x;
		}
		REQUIRE(all_in_range);
		REQUIRE(sum / n == Approx(0.5).margin(0.005));
	}

	SECTION("Normal has mean 0 and variance 1") {
		Rng rng(12345);
		const int n = 400000;
		double sum = 0.0;
		double sum_sq = 0.0;
		for (int i = 0; i < n; ++i) {
			const double z = rng.Normal();
			sum += z;
			sum_sq += z * z;
		}
		const double mean = sum / n;
		REQUIRE(mean == Approx(0.0).margin(0.01));
		REQUIRE(sum_sq / n - mean * mean == Approx(1.0).margin(0.01));
	}

	SECTION("Normal's cached second draw does not depend on how many are taken") {
		// The polar method generates pairs and caches the spare. If that cache were
		// mishandled, taking draws in different-sized chunks would diverge.
		Rng a(99);
		std::vector<double> in_one;
		for (int i = 0; i < 9; ++i) {
			in_one.push_back(a.Normal());
		}
		Rng b(99);
		std::vector<double> in_chunks;
		for (int i = 0; i < 4; ++i) {
			in_chunks.push_back(b.Normal());
		}
		for (int i = 0; i < 5; ++i) {
			in_chunks.push_back(b.Normal());
		}
		REQUIRE(in_chunks == in_one);
	}
}

TEST_CASE("weighted index sampling follows the weights", "[mmvec]") {
	SECTION("empirical frequencies converge on the weights") {
		const std::vector<double> weights = {1.0, 3.0, 0.0, 6.0};
		Rng rng(4);
		const auto draws = SampleWeightedIndices(weights, 200000, rng);
		REQUIRE(draws.size() == 200000);
		std::vector<int> hits(weights.size(), 0);
		bool all_in_range = true;
		for (const int64_t d : draws) {
			all_in_range = all_in_range && d >= 0 && d < static_cast<int64_t>(weights.size());
			hits[static_cast<size_t>(d)]++;
		}
		REQUIRE(all_in_range);
		// A zero weight must be unreachable, not merely unlikely -- otherwise an
		// absent cell could enter a minibatch.
		REQUIRE(hits[2] == 0);
		REQUIRE(static_cast<double>(hits[0]) / 200000.0 == Approx(0.1).margin(0.01));
		REQUIRE(static_cast<double>(hits[1]) / 200000.0 == Approx(0.3).margin(0.01));
		REQUIRE(static_cast<double>(hits[3]) / 200000.0 == Approx(0.6).margin(0.01));
	}

	SECTION("the same seed gives the same draws") {
		const std::vector<double> weights = {2.0, 1.0, 5.0};
		Rng a(11);
		Rng b(11);
		REQUIRE(SampleWeightedIndices(weights, 500, a) == SampleWeightedIndices(weights, 500, b));
	}

	SECTION("degenerate weights are refused") {
		Rng rng(1);
		REQUIRE_THROWS_WITH(SampleWeightedIndices({}, 5, rng), ContainsSubstring("empty weight vector"));
		REQUIRE_THROWS_WITH(SampleWeightedIndices({0.0, 0.0}, 5, rng), ContainsSubstring("sum to zero"));
		REQUIRE_THROWS_WITH(SampleWeightedIndices({1.0, -1.0}, 5, rng), ContainsSubstring("negative"));
		REQUIRE_THROWS_WITH(SampleWeightedIndices({1.0, std::nan("")}, 5, rng), ContainsSubstring("not finite"));
	}
}

TEST_CASE("InitTheta draws the right number of standard normals", "[mmvec]") {
	const ModelShape shape {6, 5, 2};
	const auto theta = InitTheta(shape, 0);
	REQUIRE(theta.size() == static_cast<size_t>(NumParams(shape)));
	for (const double t : theta) {
		REQUIRE(std::isfinite(t));
	}
	// Same seed, same parameters; different seed, different parameters.
	REQUIRE(InitTheta(shape, 0) == theta);
	REQUIRE_FALSE(InitTheta(shape, 1) == theta);
	// A degenerate shape is rejected before NumParams is used as a size.
	REQUIRE_THROWS_WITH(InitTheta(ModelShape {6, 1, 2}, 0), ContainsSubstring("n_features_y must be >= 2"));
	REQUIRE_THROWS_WITH(InitTheta(ModelShape {6, 5, 0}, 0), ContainsSubstring("n_components must be >= 1"));
}

// ---------------------------------------------------------------------------
// Adam
//
// T3 is exact, not banded, because the oracle supplies the minibatch index
// sequence -- which removes the RNG from the comparison entirely and leaves only
// arithmetic. Measured worst relative deviation across both batch_norm modes and
// all 20 recorded losses: 2.4e-16, with several values bit-identical.
//
// Gradient clipping is exercised on every one of those updates (measured global
// gradient norms of 390-670 against a clipnorm of 10), so T3 pins the clipping
// formula -- the global L2 norm across all four blocks, applied after the priors
// are already in the gradient -- as tightly as it pins everything else.
// ---------------------------------------------------------------------------
namespace {
// The Adam oracle's inputs and hyperparameters for the toy case.
struct AdamOracleCase {
	ModelShape shape;
	MinibatchInputs mb;
	AdamParams params;
	std::vector<int64_t> batches;
};

AdamOracleCase ToyAdamCase(BatchNorm mode) {
	using namespace miint::mmvec_oracle;
	using namespace miint::mmvec_oracle::toy;
	AdamOracleCase c;
	c.shape = ModelShape {kNFeaturesX, kNFeaturesY, kNComponents};
	const auto x_coo = ToCoo(kXCounts, kNSamples, kNFeaturesX);
	c.mb.n_samples = kNSamples;
	c.mb.x_row = x_coo.row;
	c.mb.x_col = x_coo.col;
	c.mb.x_val = x_coo.val;
	c.mb.y_dense = kYCounts;
	c.params.learning_rate = kAdamLearningRate;
	c.params.beta_1 = kAdamBeta1;
	c.params.beta_2 = kAdamBeta2;
	c.params.clipnorm = kAdamClipnorm;
	c.params.batch_size = kAdamBatchSize;
	c.params.batch_norm = mode;
	c.batches.assign(kAdamIndexBlocks.begin(), kAdamIndexBlocks.end());
	return c;
}
} // namespace

TEST_CASE("Adam reproduces the T3 oracle under injected minibatch indices", "[mmvec]") {
	using miint::mmvec_oracle::kT3RelTol;
	using namespace miint::mmvec_oracle::toy;

	auto check = [&](BatchNorm mode, const std::vector<double> &want_curve, double want_final) {
		auto c = ToyAdamCase(mode);
		const Model m = FitAdamWithIndices(c.shape, OraclePriors(), c.mb, c.params, kTheta, c.batches);

		REQUIRE(m.n_iter == kAdamNumUpdates);
		REQUIRE(m.loss_curve.size() == want_curve.size());
		for (size_t i = 0; i < want_curve.size(); ++i) {
			INFO("update " << i);
			REQUIRE(m.loss_curve[i] == Approx(want_curve[i]).epsilon(kT3RelTol));
		}
		// The oracle's "final loss" is the last recorded minibatch loss.
		REQUIRE(m.loss_curve.back() == Approx(want_final).epsilon(kT3RelTol));
		REQUIRE(m.theta.size() == static_cast<size_t>(NumParams(c.shape)));
	};

	SECTION("batch_norm = unbiased") {
		check(BatchNorm::Unbiased, kAdamUnbiasedLossCurve, kAdamUnbiasedFinalLoss);
	}
	SECTION("batch_norm = legacy") {
		check(BatchNorm::Legacy, kAdamLegacyLossCurve, kAdamLegacyFinalLoss);
	}
}

TEST_CASE("Adam's two batch_norm modes fit different models", "[mmvec]") {
	// Not merely different step sizes: 'legacy' rescales the likelihood relative to
	// unscaled priors, so it targets a more heavily regularized posterior. The
	// parameters must actually end up somewhere else.
	using namespace miint::mmvec_oracle::toy;
	auto u = ToyAdamCase(BatchNorm::Unbiased);
	auto l = ToyAdamCase(BatchNorm::Legacy);
	const Model mu = FitAdamWithIndices(u.shape, OraclePriors(), u.mb, u.params, kTheta, u.batches);
	const Model ml = FitAdamWithIndices(l.shape, OraclePriors(), l.mb, l.params, kTheta, l.batches);

	REQUIRE_FALSE(mu.theta == ml.theta);
	REQUIRE_FALSE(mu.loss_curve == ml.loss_curve);
}

TEST_CASE("Adam reports honestly that it has no convergence test", "[mmvec]") {
	// Adam runs its schedule and stops. Reporting `converged` would be a
	// fabrication, so it is always false and the message says why. `final_loss` is
	// the FULL-BATCH objective at the final parameters -- not the last minibatch
	// loss, which is one noisy draw and not comparable to anything.
	using namespace miint::mmvec_oracle::toy;
	auto c = ToyAdamCase(BatchNorm::Unbiased);
	const Model m = FitAdamWithIndices(c.shape, OraclePriors(), c.mb, c.params, kTheta, c.batches);

	REQUIRE_FALSE(m.converged);
	REQUIRE_THAT(m.message, ContainsSubstring("no convergence test"));
	REQUIRE(m.max_abs_grad > 0.0);

	// final_loss is the full-batch objective at theta, so it must equal what the
	// objective reports there -- and differ from the last minibatch estimate.
	Workspace ws;
	std::vector<double> grad;
	auto stats_case = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
	const double full = FullBatchLossAndGradient(c.shape, OraclePriors(), stats_case.stats, m.theta, ws, grad);
	REQUIRE(m.final_loss == Approx(full).epsilon(1e-12));
	REQUIRE_FALSE(m.final_loss == Approx(m.loss_curve.back()).epsilon(1e-6));

	// 10 Adam updates at a 1e-3 learning rate barely move the objective, but they
	// must move it the right way.
	REQUIRE(m.final_loss < kT1Loss);
}

TEST_CASE("seeded Adam is deterministic and sizes its own schedule", "[mmvec]") {
	// An epoch is a pass over X's nonzero cells, so the update count comes from the
	// data: max(1, nnz/batch_size) per epoch. toy has 27 nonzero X cells, so at a
	// batch size of 10 that is 2 updates per epoch.
	using namespace miint::mmvec_oracle::toy;
	auto c = ToyAdamCase(BatchNorm::Unbiased);
	const int64_t epochs = 5;
	const auto nnz = static_cast<int64_t>(c.mb.x_val.size());
	REQUIRE(nnz == 27);
	const int64_t expected_updates = epochs * std::max<int64_t>(1, nnz / c.params.batch_size);
	REQUIRE(expected_updates == 10);

	const Model a = FitAdam(c.shape, OraclePriors(), c.mb, c.params, kTheta, epochs, 0);
	const Model b = FitAdam(c.shape, OraclePriors(), c.mb, c.params, kTheta, epochs, 0);
	const Model d = FitAdam(c.shape, OraclePriors(), c.mb, c.params, kTheta, epochs, 1);

	REQUIRE(a.n_iter == expected_updates);
	REQUIRE(a.theta == b.theta);
	REQUIRE(a.loss_curve == b.loss_curve);
	// A different seed draws different minibatches, so it lands elsewhere.
	REQUIRE_FALSE(a.theta == d.theta);

	// A batch larger than the whole data set still performs one update per epoch
	// rather than none.
	AdamParams big = c.params;
	big.batch_size = 1000;
	const Model wide = FitAdam(c.shape, OraclePriors(), c.mb, big, kTheta, 3, 0);
	REQUIRE(wide.n_iter == 3);
}

TEST_CASE("gradient clipping is wired into the Adam update", "[mmvec]") {
	// T3 pins the clipping arithmetic exactly (it fires on all 10 updates there).
	// This adds the cruder guarantee that the threshold is actually consulted: a
	// clipnorm far below the observed gradient norms must change the trajectory,
	// and one far above must leave it untouched.
	using namespace miint::mmvec_oracle::toy;
	auto c = ToyAdamCase(BatchNorm::Unbiased);

	AdamParams tight = c.params;
	tight.clipnorm = 1e-6;
	AdamParams loose = c.params;
	loose.clipnorm = 1e12; // above every observed norm, so clipping never triggers

	const Model reference = FitAdamWithIndices(c.shape, OraclePriors(), c.mb, c.params, kTheta, c.batches);
	const Model clipped = FitAdamWithIndices(c.shape, OraclePriors(), c.mb, tight, kTheta, c.batches);
	const Model unclipped = FitAdamWithIndices(c.shape, OraclePriors(), c.mb, loose, kTheta, c.batches);

	REQUIRE_FALSE(clipped.theta == reference.theta);
	REQUIRE_FALSE(unclipped.theta == reference.theta);
}

TEST_CASE("Adam rejects invalid hyperparameters and batch sequences", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;
	auto c = ToyAdamCase(BatchNorm::Unbiased);

	SECTION("a batch sequence that is not a whole number of batches") {
		auto ragged = c.batches;
		ragged.pop_back();
		REQUIRE_THROWS_WITH(FitAdamWithIndices(c.shape, OraclePriors(), c.mb, c.params, kTheta, ragged),
		                    ContainsSubstring("not a whole number of batches"));
	}
	SECTION("no batches at all") {
		REQUIRE_THROWS_WITH(FitAdamWithIndices(c.shape, OraclePriors(), c.mb, c.params, kTheta, {}),
		                    ContainsSubstring("no minibatches supplied"));
	}
	SECTION("out-of-range hyperparameters") {
		AdamParams bad = c.params;
		bad.batch_size = 0;
		REQUIRE_THROWS_WITH(FitAdamWithIndices(c.shape, OraclePriors(), c.mb, bad, kTheta, c.batches),
		                    ContainsSubstring("batch_size must be >= 1"));
		bad = c.params;
		bad.learning_rate = 0.0;
		REQUIRE_THROWS_WITH(FitAdamWithIndices(c.shape, OraclePriors(), c.mb, bad, kTheta, c.batches),
		                    ContainsSubstring("learning_rate must be > 0"));
		bad = c.params;
		bad.beta_1 = 1.0; // would make the bias correction 1 - 1 = 0
		REQUIRE_THROWS_WITH(FitAdamWithIndices(c.shape, OraclePriors(), c.mb, bad, kTheta, c.batches),
		                    ContainsSubstring("beta_1 and beta_2"));
		bad = c.params;
		bad.clipnorm = 0.0;
		REQUIRE_THROWS_WITH(FitAdamWithIndices(c.shape, OraclePriors(), c.mb, bad, kTheta, c.batches),
		                    ContainsSubstring("clipnorm must be > 0"));
	}
	SECTION("epochs below one") {
		REQUIRE_THROWS_WITH(FitAdam(c.shape, OraclePriors(), c.mb, c.params, kTheta, 0, 0),
		                    ContainsSubstring("epochs must be >= 1"));
	}
	SECTION("X counts missing for the sampling weights") {
		auto no_vals = c.mb;
		no_vals.x_val.clear();
		REQUIRE_THROWS_WITH(FitAdam(c.shape, OraclePriors(), no_vals, c.params, kTheta, 5, 0),
		                    ContainsSubstring("needs X counts"));
	}
}

// ---------------------------------------------------------------------------
// Derived quantities: ranks, probs, predict, score
//
// These four are what a fitted model is read through, and they are the only
// outputs worth comparing across implementations at all: the embeddings are
// identified only up to an orthogonal rotation of the latent axes, so two fits
// that agree perfectly here can still have unrelated x_main/y_main. Everything
// below therefore checks a rotation-invariant quantity.
//
// The T1 cases pin all four against the oracle at the carved theta -- no fitting
// involved, so a disagreement is a formula error and nothing else. The invariant
// cases catch the errors T1 structurally cannot: T1 would pass just as happily if
// `Ranks` centred columns instead of rows on a square matrix, or if `Predict`
// consumed raw counts on a fixture whose rows happened to have equal depth.
// ---------------------------------------------------------------------------
namespace {
using miint::mmvec::Predict;
using miint::mmvec::Probs;
using miint::mmvec::Ranks;
using miint::mmvec::Score;

// Row-wise softmax written independently of the core's, so that comparing it
// against Probs() tests the shift-invariance identity rather than restating one
// implementation twice.
std::vector<double> SoftmaxRows(const std::vector<double> &m, int64_t rows, int64_t cols) {
	std::vector<double> out(m.size());
	for (int64_t i = 0; i < rows; ++i) {
		const auto begin = m.begin() + static_cast<std::ptrdiff_t>(i * cols);
		const double max = *std::max_element(begin, begin + static_cast<std::ptrdiff_t>(cols));
		double total = 0.0;
		for (int64_t j = 0; j < cols; ++j) {
			out[static_cast<size_t>(i * cols + j)] = std::exp(m[static_cast<size_t>(i * cols + j)] - max);
			total += out[static_cast<size_t>(i * cols + j)];
		}
		for (int64_t j = 0; j < cols; ++j) {
			out[static_cast<size_t>(i * cols + j)] /= total;
		}
	}
	return out;
}

double RowSum(const std::vector<double> &m, int64_t row, int64_t cols) {
	const auto begin = m.begin() + static_cast<std::ptrdiff_t>(row * cols);
	return std::accumulate(begin, begin + static_cast<std::ptrdiff_t>(cols), 0.0);
}
} // namespace

TEST_CASE("derived quantities reproduce the T1 oracle", "[mmvec]") {
	using miint::mmvec_oracle::kT1Tol;

	// `px`/`py` are the predict/score inputs: held-out samples for toy, and the
	// training matrices for synth_a/synth_b, exactly as the oracle recorded them.
	auto check = [&](int64_t d1, int64_t d2, int32_t p, const std::vector<double> &theta,
	                 const std::vector<double> &want_ranks, const std::vector<double> &want_probs,
	                 const std::vector<double> &want_predict, double want_score, const std::vector<double> &px,
	                 const std::vector<double> &py, int64_t n_predict) {
		const ModelShape shape {d1, d2, p};
		REQUIRE(theta.size() == static_cast<size_t>(NumParams(shape)));

		const auto ranks = Ranks(shape, theta);
		REQUIRE(ranks.size() == want_ranks.size());
		for (size_t i = 0; i < want_ranks.size(); ++i) {
			INFO("rank element " << i);
			REQUIRE(ranks[i] == Approx(want_ranks[i]).margin(kT1Tol));
		}

		const auto probs = Probs(shape, theta);
		REQUIRE(probs.size() == want_probs.size());
		for (size_t i = 0; i < want_probs.size(); ++i) {
			INFO("prob element " << i);
			REQUIRE(probs[i] == Approx(want_probs[i]).margin(kT1Tol));
		}

		const auto x_coo = ToCoo(px, n_predict, d1);
		const auto predict = Predict(shape, theta, x_coo);
		REQUIRE(predict.size() == want_predict.size());
		for (size_t i = 0; i < want_predict.size(); ++i) {
			INFO("predict element " << i);
			REQUIRE(predict[i] == Approx(want_predict[i]).margin(kT1Tol));
		}

		REQUIRE(Score(shape, theta, x_coo, ToCoo(py, n_predict, d2)) == Approx(want_score).margin(kT1Tol));
	};

	SECTION("toy") {
		using namespace miint::mmvec_oracle::toy;
		// The only fixture with genuinely held-out predict/score samples.
		check(kNFeaturesX, kNFeaturesY, kNComponents, kTheta, kT1Ranks, kT1Probs, kT1Predict, kT1Score, kPredictXCounts,
		      kScoreYCounts, static_cast<int64_t>(kPredictSampleIds.size()));
	}
	SECTION("synth_a") {
		using namespace miint::mmvec_oracle::synth_a;
		check(kNFeaturesX, kNFeaturesY, kNComponents, kTheta, kT1Ranks, kT1Probs, kT1Predict, kT1Score, kXCounts,
		      kYCounts, kNSamples);
	}
	SECTION("synth_b") {
		using namespace miint::mmvec_oracle::synth_b;
		check(kNFeaturesX, kNFeaturesY, kNComponents, kTheta, kT1Ranks, kT1Probs, kT1Predict, kT1Score, kXCounts,
		      kYCounts, kNSamples);
	}
}

TEST_CASE("ranks are the logits with each row centred", "[mmvec]") {
	using miint::mmvec_oracle::kRowCenteringTol;
	using namespace miint::mmvec_oracle::toy;
	const ModelShape shape {kNFeaturesX, kNFeaturesY, kNComponents};

	const auto logits = ComputeLogits(shape, kTheta);
	const auto ranks = Ranks(shape, kTheta);

	for (int64_t i = 0; i < kNFeaturesX; ++i) {
		INFO("X feature " << i);
		// Rows summing to zero is the defining property: it is what makes one X
		// feature's ranks comparable to another's, and what removes the influence
		// of which Y feature happened to be the pinned reference.
		REQUIRE(RowSum(ranks, i, kNFeaturesY) == Approx(0.0).margin(kRowCenteringTol));

		// Centring is a per-row shift, so every within-row DIFFERENCE must survive
		// it untouched. This is what distinguishes row-centring from column-
		// centring or from a global shift, none of which T1 could tell apart on a
		// fixture with equal row and column counts.
		for (int64_t j = 1; j < kNFeaturesY; ++j) {
			const size_t a = static_cast<size_t>(i * kNFeaturesY);
			const size_t b = static_cast<size_t>(i * kNFeaturesY + j);
			REQUIRE((ranks[b] - ranks[a]) == Approx(logits[b] - logits[a]).margin(1e-12));
		}
	}
}

TEST_CASE("probs are a proper distribution and agree with softmax(ranks)", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;
	const ModelShape shape {kNFeaturesX, kNFeaturesY, kNComponents};

	const auto probs = Probs(shape, kTheta);
	for (int64_t i = 0; i < kNFeaturesX; ++i) {
		INFO("X feature " << i);
		REQUIRE(RowSum(probs, i, kNFeaturesY) == Approx(1.0).margin(1e-12));
	}
	for (const double v : probs) {
		REQUIRE(v > 0.0);
		REQUIRE(v < 1.0);
	}

	// softmax is shift-invariant, so centring the logits first cannot change the
	// result. The core computes probs from the logits; this recomputes them from
	// the ranks through an independent softmax, so the two paths agreeing is a
	// real check on both rather than a tautology.
	const auto from_ranks = SoftmaxRows(Ranks(shape, kTheta), kNFeaturesX, kNFeaturesY);
	REQUIRE(from_ranks.size() == probs.size());
	for (size_t i = 0; i < probs.size(); ++i) {
		INFO("prob element " << i);
		REQUIRE(from_ranks[i] == Approx(probs[i]).margin(1e-14));
	}
}

TEST_CASE("probs survive logits far outside exp's range", "[mmvec]") {
	// The row-maximum subtraction in the softmax is an algebraic identity, so
	// every comparison against the oracle passes with or without it -- the same
	// blind spot the objective's log-normalizer had. Only a theta extreme enough
	// to overflow exp() can tell the two apart: here the non-reference logits
	// reach ~2760, where a naive exp() gives inf and inf/inf gives NaN for every
	// probability. A fit that wanders somewhere extreme must still produce
	// readable output rather than a matrix of NaNs.
	using namespace miint::mmvec_oracle::toy;
	const ModelShape shape {kNFeaturesX, kNFeaturesY, kNComponents};
	const std::vector<double> wild(static_cast<size_t>(NumParams(shape)), 30.0);

	const auto probs = Probs(shape, wild);
	for (const double v : probs) {
		REQUIRE(std::isfinite(v));
		REQUIRE(v >= 0.0);
	}
	for (int64_t i = 0; i < kNFeaturesX; ++i) {
		INFO("X feature " << i);
		REQUIRE(RowSum(probs, i, kNFeaturesY) == Approx(1.0).margin(1e-12));
	}
	// Ranks stay finite too -- they are a shift, so nothing can overflow, but the
	// row-centring still has to hold at this magnitude.
	const auto ranks = Ranks(shape, wild);
	for (int64_t i = 0; i < kNFeaturesX; ++i) {
		REQUIRE(RowSum(ranks, i, kNFeaturesY) == Approx(0.0).margin(1e-9));
	}
}

TEST_CASE("predict reads X as proportions, not as counts", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;
	const ModelShape shape {kNFeaturesX, kNFeaturesY, kNComponents};
	const int64_t n = static_cast<int64_t>(kPredictSampleIds.size());
	const auto base = ToCoo(kPredictXCounts, n, kNFeaturesX);
	const auto want = Predict(shape, kTheta, base);

	SECTION("rows are distributions") {
		for (int64_t s = 0; s < n; ++s) {
			INFO("sample " << s);
			REQUIRE(RowSum(want, s, kNFeaturesY) == Approx(1.0).margin(1e-12));
		}
	}

	SECTION("scaling one sample's depth leaves its prediction unchanged") {
		// Sequencing depth is an artefact of the assay, not biology. A sample
		// sequenced twice as deeply must get the same predicted profile, which
		// holds only if X is row-normalized before the product -- an
		// implementation that used raw counts would pass every T1 comparison and
		// fail here.
		auto scaled = base;
		for (size_t k = 0; k < scaled.val.size(); ++k) {
			if (scaled.row[k] == 1) {
				scaled.val[k] *= 7.0;
			}
		}
		const auto got = Predict(shape, kTheta, scaled);
		for (size_t i = 0; i < want.size(); ++i) {
			INFO("predict element " << i);
			REQUIRE(got[i] == Approx(want[i]).margin(1e-12));
		}
	}

	SECTION("the result does not depend on the order of the COO entries") {
		// The SQL layer cannot promise a scan order. Unlike the sufficient
		// statistics this is not automatic -- an output cell takes one term per X
		// nonzero in its row, so the summation order has to be imposed -- and
		// bit-identity is the assertion because imposing it is the whole point.
		const auto got = Predict(shape, kTheta, Reversed(base));
		REQUIRE(got == want);
	}
}

TEST_CASE("score is Q^2 against the mean-community baseline", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;
	const ModelShape shape {kNFeaturesX, kNFeaturesY, kNComponents};
	const int64_t n = static_cast<int64_t>(kPredictSampleIds.size());
	const auto x_coo = ToCoo(kPredictXCounts, n, kNFeaturesX);

	SECTION("it matches the definition computed independently") {
		// Recomputed here from Predict's output and the Y proportions, so this
		// pins the SS_res/SS_tot assembly -- in particular that the baseline is a
		// PER-FEATURE mean over samples, not one global mean, which is a plausible
		// misreading that would still give a number in the right range.
		const auto pred = Predict(shape, kTheta, x_coo);
		std::vector<double> y_props(static_cast<size_t>(n * kNFeaturesY), 0.0);
		for (int64_t s = 0; s < n; ++s) {
			const double total = RowSum(kScoreYCounts, s, kNFeaturesY);
			for (int64_t j = 0; j < kNFeaturesY; ++j) {
				y_props[static_cast<size_t>(s * kNFeaturesY + j)] =
				    kScoreYCounts[static_cast<size_t>(s * kNFeaturesY + j)] / total;
			}
		}
		std::vector<double> col_mean(static_cast<size_t>(kNFeaturesY), 0.0);
		for (int64_t s = 0; s < n; ++s) {
			for (int64_t j = 0; j < kNFeaturesY; ++j) {
				col_mean[static_cast<size_t>(j)] +=
				    y_props[static_cast<size_t>(s * kNFeaturesY + j)] / static_cast<double>(n);
			}
		}
		double ss_res = 0.0;
		double ss_tot = 0.0;
		for (int64_t s = 0; s < n; ++s) {
			for (int64_t j = 0; j < kNFeaturesY; ++j) {
				const size_t i = static_cast<size_t>(s * kNFeaturesY + j);
				ss_res += (y_props[i] - pred[i]) * (y_props[i] - pred[i]);
				ss_tot +=
				    (y_props[i] - col_mean[static_cast<size_t>(j)]) * (y_props[i] - col_mean[static_cast<size_t>(j)]);
			}
		}
		REQUIRE(Score(shape, kTheta, x_coo, ToCoo(kScoreYCounts, n, kNFeaturesY)) ==
		        Approx(1.0 - ss_res / ss_tot).margin(1e-12));
	}

	SECTION("it is exactly zero for a single sample") {
		// With one sample the mean community IS that sample, so SS_tot is exactly
		// zero and Q^2 is 0/0. scikit-bio returns 0.0 and so do we; the assertion
		// is on the exact value because the alternative is a NaN propagating
		// silently through anything downstream.
		//
		// One sample is the only case that reaches this branch reliably. Several
		// IDENTICAL samples do NOT: colmean sums n copies of a value and divides
		// by n, which is off by an ulp for most values (measured: 5/21 over three
		// samples), leaving SS_tot at ~1e-33 rather than 0. The guard tests exact
		// equality, matching scikit-bio, so an almost-constant Y is not caught and
		// scores hugely negative instead. Widening the guard to a tolerance would
		// change results relative to scikit-bio, so it is a deliberate M6 ledger
		// question rather than something to quietly patch here.
		const std::vector<double> one_x(kPredictXCounts.begin(),
		                                kPredictXCounts.begin() + static_cast<std::ptrdiff_t>(kNFeaturesX));
		const std::vector<double> one_y(kScoreYCounts.begin(),
		                                kScoreYCounts.begin() + static_cast<std::ptrdiff_t>(kNFeaturesY));
		REQUIRE(Score(shape, kTheta, ToCoo(one_x, 1, kNFeaturesX), ToCoo(one_y, 1, kNFeaturesY)) == 0.0);
	}
}

TEST_CASE("predict and score reject unusable inputs", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;
	const ModelShape shape {kNFeaturesX, kNFeaturesY, kNComponents};
	const int64_t n = static_cast<int64_t>(kPredictSampleIds.size());
	const auto x_coo = ToCoo(kPredictXCounts, n, kNFeaturesX);
	const auto y_coo = ToCoo(kScoreYCounts, n, kNFeaturesY);

	SECTION("a sample with no counts at all") {
		// Proportions are undefined for it, so there is nothing to predict.
		auto x = x_coo;
		x.n_rows += 1;
		REQUIRE_THROWS_WITH(Predict(shape, kTheta, x), ContainsSubstring("all-zero row"));
		REQUIRE_THROWS_WITH(Predict(shape, kTheta, x), ContainsSubstring("proportions"));
	}

	SECTION("an X whose width is not the model's") {
		// scikit-bio does not check this and consumes X positionally; we do. A
		// differently-filtered X table is a realistic mistake and silently
		// produces numbers otherwise.
		auto narrow = x_coo;
		narrow.n_cols = kNFeaturesX - 1;
		narrow.row.clear();
		narrow.col.clear();
		narrow.val.clear();
		for (int64_t s = 0; s < n; ++s) {
			for (int64_t c = 0; c < kNFeaturesX - 1; ++c) {
				narrow.row.push_back(s);
				narrow.col.push_back(c);
				narrow.val.push_back(1.0);
			}
		}
		REQUIRE_THROWS_WITH(Predict(shape, kTheta, narrow), ContainsSubstring("features but the model was fitted on"));
	}

	SECTION("a Y whose width is not the model's") {
		auto wide = y_coo;
		wide.n_cols = kNFeaturesY + 1;
		REQUIRE_THROWS_WITH(Score(shape, kTheta, x_coo, wide),
		                    ContainsSubstring("features but the model was fitted on"));
	}

	SECTION("X and Y disagreeing on the sample count") {
		auto y = y_coo;
		y.n_rows += 1;
		REQUIRE_THROWS_WITH(Score(shape, kTheta, x_coo, y), ContainsSubstring("same number of samples"));
	}

	SECTION("a non-finite theta") {
		auto bad = kTheta;
		bad[3] = std::numeric_limits<double>::quiet_NaN();
		REQUIRE_THROWS_WITH(Ranks(shape, bad), ContainsSubstring("not finite"));
		REQUIRE_THROWS_WITH(Probs(shape, bad), ContainsSubstring("not finite"));
		REQUIRE_THROWS_WITH(Predict(shape, bad, x_coo), ContainsSubstring("not finite"));
		REQUIRE_THROWS_WITH(Score(shape, bad, x_coo, y_coo), ContainsSubstring("not finite"));
	}

	SECTION("a feature absent from the supplied samples is accepted") {
		// Fitting rejects an all-zero feature column, because then the feature's
		// parameters are set by the prior alone. Predicting must NOT: the columns
		// are only weights on already-fitted rows of probs, and a microbe missing
		// from three held-out samples is ordinary data, not a broken table.
		std::vector<double> dense(kPredictXCounts);
		for (int64_t s = 0; s < n; ++s) {
			dense[static_cast<size_t>(s * kNFeaturesX)] = 0.0; // wipe feature 0
		}
		const auto sparse_x = ToCoo(dense, n, kNFeaturesX);
		REQUIRE_NOTHROW(Predict(shape, kTheta, sparse_x));
		REQUIRE_NOTHROW(Score(shape, kTheta, sparse_x, y_coo));
	}
}

TEST_CASE("a converged fit reproduces the T2 oracle's derived quantities", "[mmvec]") {
	// Deferred out of the L-BFGS phase, which could only check the objective: the
	// parameters themselves are not comparable across implementations (rotation),
	// so this is where a converged fit is actually held to the oracle.
	//
	// The tolerance is measured, not chosen. kT2RankTol is 0.01 against a spread
	// of 4.1e-04 in max|dranks| across four starting points under scikit-bio's own
	// stopping rule -- the rule stops well short of a stationary point, so two
	// correct implementations land at genuinely different nearby points. Tightening
	// this would be testing the line search, not the model.
	using namespace miint::mmvec_oracle::toy;

	auto c = MakeCase(kXCounts, kYCounts, kNSamples, kNFeaturesX, kNFeaturesY, kNComponents);
	const Model m = FitLbfgsFromInit(c.shape, OraclePriors(), c.stats, kTheta, kT2MaxIter);
	REQUIRE(m.converged);

	const auto ranks = Ranks(m.shape, m.theta);
	REQUIRE(ranks.size() == kT2Ranks.size());
	for (size_t i = 0; i < kT2Ranks.size(); ++i) {
		INFO("rank element " << i);
		REQUIRE(ranks[i] == Approx(kT2Ranks[i]).margin(kT2RankTol));
	}

	const auto probs = Probs(m.shape, m.theta);
	REQUIRE(probs.size() == kT2Probs.size());
	for (size_t i = 0; i < kT2Probs.size(); ++i) {
		INFO("prob element " << i);
		REQUIRE(probs[i] == Approx(kT2Probs[i]).margin(kT2RankTol));
	}

	const int64_t n = static_cast<int64_t>(kPredictSampleIds.size());
	REQUIRE(Score(m.shape, m.theta, ToCoo(kPredictXCounts, n, kNFeaturesX), ToCoo(kScoreYCounts, n, kNFeaturesY)) ==
	        Approx(kT2Score).margin(kT2RankTol));

	// The fit is worth something on held-out samples, where the unfitted T1 theta
	// scores about -2. Cheap, but it is the only assertion here that would catch a
	// fit that reproduced the oracle's ranks while predicting nothing.
	REQUIRE(kT2Score > 0.0);
}

TEST_CASE("Adam's T3 fit reproduces the oracle's ranks", "[mmvec]") {
	// Deferred out of the Adam phase for the same reason, but at a completely
	// different tolerance: with the minibatch indices injected there is no
	// sampling difference left, so the two implementations follow the same
	// trajectory and the ranks agree to rounding rather than to a measured band.
	using miint::mmvec_oracle::kT3RelTol;
	using namespace miint::mmvec_oracle::toy;

	auto check = [&](BatchNorm mode, const std::vector<double> &want) {
		auto c = ToyAdamCase(mode);
		const Model m = FitAdamWithIndices(c.shape, OraclePriors(), c.mb, c.params, kTheta, c.batches);
		const auto ranks = Ranks(m.shape, m.theta);
		REQUIRE(ranks.size() == want.size());
		for (size_t i = 0; i < want.size(); ++i) {
			INFO("rank element " << i);
			// Relative OR absolute: a rank that lands near zero has no meaningful
			// relative scale, and Catch2's Approx accepts either test.
			REQUIRE(ranks[i] == Approx(want[i]).epsilon(kT3RelTol).margin(kT3RelTol));
		}
	};

	SECTION("batch_norm = unbiased") {
		check(BatchNorm::Unbiased, kAdamUnbiasedRanks);
	}
	SECTION("batch_norm = legacy") {
		check(BatchNorm::Legacy, kAdamLegacyRanks);
	}
}
