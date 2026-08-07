//
// mmvec's numeric kernel. Split out of mmvec.cpp so it can be compiled more than
// once at different instruction sets -- Eigen fixes its packet width at compile
// time, so a wider ISA is a different translation unit, not a runtime choice.
// See mmvec_kernel.hpp for why the seam sits exactly here.
//
// Nothing in this file may depend on anything but its arguments: the same source
// is compiled several times, and every copy has to compute the same function of
// its inputs (bit-for-bit within one instruction set).

#include "mmvec_kernel.hpp"

#include <cmath>

namespace miint::mmvec {

namespace {

//! Sum of squared deviations of `n` values from `mean`, i.e. `||v - mean||^2`.
double SumSquaredDev(const double *v, size_t n, double mean) {
	double acc = 0.0;
	for (size_t i = 0; i < n; ++i) {
		const double d = v[i] - mean;
		acc += d * d;
	}
	return acc;
}

//! The Gaussian prior's contribution to the objective:
//! `0.5 * ||block - mean||^2 / scale^2`, summed over the four blocks.
double PriorLoss(const ParamLayout &l, const Priors &priors, const double *theta) {
	const double x_sq = priors.x_prior_scale * priors.x_prior_scale;
	const double y_sq = priors.y_prior_scale * priors.y_prior_scale;
	double loss = 0.0;
	loss += 0.5 * SumSquaredDev(theta + l.x_main, static_cast<size_t>(l.d1 * l.p), priors.x_prior_mean) / x_sq;
	loss += 0.5 * SumSquaredDev(theta + l.x_bias, static_cast<size_t>(l.d1), priors.x_prior_mean) / x_sq;
	loss += 0.5 * SumSquaredDev(theta + l.y_main, static_cast<size_t>(l.p * l.nref), priors.y_prior_mean) / y_sq;
	loss += 0.5 * SumSquaredDev(theta + l.y_bias, static_cast<size_t>(l.nref), priors.y_prior_mean) / y_sq;
	return loss;
}

//! Add `(block - mean) / scale^2` -- the prior's gradient -- to each block.
void AddPriorGradient(const ParamLayout &l, const Priors &priors, const double *theta, double *grad) {
	const double x_sq = priors.x_prior_scale * priors.x_prior_scale;
	const double y_sq = priors.y_prior_scale * priors.y_prior_scale;
	const auto add = [&](size_t begin, size_t n, double mean, double scale_sq) {
		for (size_t i = begin; i < begin + n; ++i) {
			grad[i] += (theta[i] - mean) / scale_sq;
		}
	};
	add(l.x_main, static_cast<size_t>(l.d1 * l.p), priors.x_prior_mean, x_sq);
	add(l.x_bias, static_cast<size_t>(l.d1), priors.x_prior_mean, x_sq);
	add(l.y_main, static_cast<size_t>(l.p * l.nref), priors.y_prior_mean, y_sq);
	add(l.y_bias, static_cast<size_t>(l.nref), priors.y_prior_mean, y_sq);
}

} // namespace

//! The model equation, in the one place it is written:
//!
//!     out(r, j) = x_main[i_r] . y_main[:,j] + x_bias[i_r] + y_bias[j]
//!
//! for the d2-1 NON-REFERENCE Y features, over `n_rows` rows whose X features are
//! given by `x_index` (nullptr = identity). `out` is row-major (n_rows x d2-1).
//!
//! The rows' X-side parameters are gathered into `ws.x_main_rows` /
//! `ws.x_bias_rows` on the way, because the gradient needs the same gathered
//! matrix afterwards. For an identity map that gather is a straight copy of
//! x_main -- d1*p elements, which is what buys a single code path for both
//! objectives instead of the same algebra written twice.
void FillNonRefLogits(const ParamLayout &l, const double *theta, const int64_t *x_index, int64_t n_rows, Workspace &ws,
                      double *out) {
	const ConstRowMatrixMap y_main(theta + l.y_main, l.p, l.nref);
	RowMatrixMap x_main_rows(ws.x_main_rows.data(), n_rows, l.p);
	for (int64_t r = 0; r < n_rows; ++r) {
		const int64_t i = MapRow(x_index, r);
		for (int64_t k = 0; k < l.p; ++k) {
			x_main_rows(r, k) = theta[l.x_main + static_cast<size_t>(i * l.p + k)];
		}
	}

	RowMatrixMap logits(out, n_rows, l.nref);
	logits.noalias() = x_main_rows * y_main;
	for (int64_t r = 0; r < n_rows; ++r) {
		// Read straight from theta. x_main has to be gathered because Eigen needs
		// it contiguous for the product above; the bias is a scalar per row and
		// gathering it bought nothing but a buffer.
		const double xb = theta[l.x_bias + static_cast<size_t>(MapRow(x_index, r))];
		for (int64_t j = 0; j < l.nref; ++j) {
			logits(r, j) += xb + theta[l.y_bias + static_cast<size_t>(j)];
		}
	}
}

//! One evaluation of the negative log-posterior and its gradient.
//!
//! Both the full-batch and the minibatch objective are this function; they are
//! genuinely the same computation over a different set of rows. A "row" is an
//! X-feature paired with an observed non-reference Y count vector and that
//! vector's total, and the two paths supply:
//!
//!                       | full batch                 | minibatch
//!   ------------------- | -------------------------- | -----------------------
//!   n_rows              | d1                         | batch size
//!   x_index             | identity (nullptr)         | X feature per draw
//!   obs                 | y_sums, no indirection     | y_dense via sample index
//!   totals              | n_sums                     | Y row sums of the draws
//!   norm                | 1                          | BatchNormFactor(...)
//!
//! `norm` scales the data term only. Priors are added at full strength either
//! way, because they are not part of the sampled sum.
//!
//! `grad` must already have `l.total` elements. Every element is written here --
//! the dy_* blocks by assignment, the dx_* blocks zeroed and then accumulated --
//! so callers need only size it, not clear it.
double EvaluateObjective(const ParamLayout &l, const Priors &priors, int64_t n_rows, const int64_t *x_index,
                         const ObsView &obs, const double *totals, double norm, const double *theta, Workspace &ws,
                         double *grad) {
	const ConstRowMatrixMap y_main(theta + l.y_main, l.p, l.nref);
	FillNonRefLogits(l, theta, x_index, n_rows, ws, ws.logits.data());
	const ConstRowMatrixMap x_main_rows(ws.x_main_rows.data(), n_rows, l.p);
	const RowMatrixMap logits(ws.logits.data(), n_rows, l.nref);

	// Log-normalizer with the reference category left implicit at logit 0:
	//   m        = max(0, rowmax(logits))
	//   log_norm = m + log(exp(-m) + SUM_j exp(logits - m))
	//
	// The `exp(-m)` term is the reference category: it is exp(0 - m), the
	// reference's own shifted exponential, which is how the zero column enters
	// without being materialized.
	//
	// The shift by m is an algebraic identity for ANY m, so clamping it at 0
	// changes no value -- it buys numerical range. m >= 0 keeps every shifted
	// exponential, the reference's included, at or below 1, so nothing overflows;
	// without the clamp a row of strongly negative logits would send exp(-m) to
	// infinity. resids carries exp(logits - m) forward for reuse.
	RowMatrixMap resids(ws.resids.data(), n_rows, l.nref);
	double data = 0.0;
	for (int64_t r = 0; r < n_rows; ++r) {
		// A plain compare rather than std::fmax, which GCC emits as a CALL into libm
		// (4.4% of a cystic-fibrosis fit, its PLT stub included) because fmax's
		// NaN-propagation semantics are not maxsd's. The two agree here by
		// construction: `m` starts at 0.0 and a NaN logit leaves it unchanged under
		// either form, so neither can ever make `m` NaN. Measured bit-identical on
		// every fixture, both optimizers, at 200 and 1000 iterations.
		double m = 0.0;
		for (int64_t j = 0; j < l.nref; ++j) {
			const double v = logits(r, j);
			if (v > m) {
				m = v;
			}
		}
		double shifted_sum = 0.0;
		for (int64_t j = 0; j < l.nref; ++j) {
			const double e = std::exp(logits(r, j) - m);
			resids(r, j) = e;
			shifted_sum += e;
		}
		// The shifted normalizer, kept as well as its log: it IS the row scale the
		// probabilities divide by, so recovering it as exp(log_norm - m) would be a
		// log-and-exp round trip of a number already in hand -- an extra
		// transcendental per row per evaluation, and a double rounding that bites
		// hardest when m is large, which is exactly where precision matters.
		const double row_scale = std::exp(-m) + shifted_sum;
		const double log_norm = m + std::log(row_scale);

		// Data term for this row, and the residual it contributes to the
		// gradient. resids goes exp-shifted -> probability -> residual in place.
		const double *obs_row = obs.base + MapRow(obs.row_map, r) * obs.stride;
		const double total = totals[static_cast<size_t>(r)];
		double dot_obs_logits = 0.0;
		for (int64_t j = 0; j < l.nref; ++j) {
			const double o = obs_row[j + 1];
			dot_obs_logits += o * logits(r, j);
			resids(r, j) = o - total * (resids(r, j) / row_scale);
		}
		data += total * log_norm - dot_obs_logits;
	}

	// Gradient of the data term. dy_* land directly in `grad`; the X-side ones go
	// through a per-row buffer first because a minibatch can name the same X
	// feature several times and must accumulate, in batch order, exactly as
	// numpy's add.at does.
	const ConstRowMatrixMap resids_const(ws.resids.data(), n_rows, l.nref);
	RowMatrixMap dy_main(grad + l.y_main, l.p, l.nref);
	dy_main.noalias() = -norm * (x_main_rows.transpose() * resids_const);
	// One row-major sweep rather than nref column reductions. `resids` is row-major,
	// so `.col(j).sum()` strides nref*8 bytes through every row and pays a cache miss
	// per element -- 8.2% of a cystic-fibrosis fit for the same 1.25M additions the
	// row-contiguous X-side scatter below does in 1.4%. Fusing cuts that phase by 75%.
	//
	// Bit-identical, not merely equivalent: each j still accumulates over ascending r,
	// which is the order Eigen's redux uses for a strided vector. Verified on every
	// fixture, both optimizers, at 200 and 1000 iterations.
	double *dy_bias = grad + l.y_bias;
	for (int64_t j = 0; j < l.nref; ++j) {
		dy_bias[j] = 0.0;
	}
	for (int64_t r = 0; r < n_rows; ++r) {
		const double *resid_row = ws.resids.data() + static_cast<size_t>(r) * static_cast<size_t>(l.nref);
		for (int64_t j = 0; j < l.nref; ++j) {
			dy_bias[j] += resid_row[j];
		}
	}
	for (int64_t j = 0; j < l.nref; ++j) {
		dy_bias[j] *= -norm;
	}

	RowMatrixMap dx_main_rows(ws.dx_main_rows.data(), n_rows, l.p);
	dx_main_rows.noalias() = resids_const * y_main.transpose();
	for (size_t i = 0; i < static_cast<size_t>(l.d1 * l.p); ++i) {
		grad[l.x_main + i] = 0.0;
	}
	for (size_t i = 0; i < static_cast<size_t>(l.d1); ++i) {
		grad[l.x_bias + i] = 0.0;
	}
	for (int64_t r = 0; r < n_rows; ++r) {
		const int64_t i = MapRow(x_index, r);
		for (int64_t k = 0; k < l.p; ++k) {
			grad[l.x_main + static_cast<size_t>(i * l.p + k)] += -norm * dx_main_rows(r, k);
		}
		// Same reduction, taken where it is consumed rather than staged in a
		// buffer first. resids_const is not touched by the scatter, so the value
		// and the summation order are unchanged.
		grad[l.x_bias + static_cast<size_t>(i)] += -norm * resids_const.row(r).sum();
	}

	AddPriorGradient(l, priors, theta, grad);
	return norm * data + PriorLoss(l, priors, theta);
}

} // namespace miint::mmvec
