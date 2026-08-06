#include "mmvec.hpp"

#include <Eigen/Core>

#include "LBFGSpp/LBFGS.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace miint::mmvec {

namespace {

//! Row-major, to match the packed parameter layout so Eigen can map the blocks
//! of `theta` and `grad` in place with no copying or transposing.
using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using RowMatrixMap = Eigen::Map<RowMatrix>;
using ConstRowMatrixMap = Eigen::Map<const RowMatrix>;

//! A validated count table regrouped by sample, so sample `s` owns entries
//! `[offsets[s], offsets[s+1])`. Indices are stored as size_t because they have
//! already been range-checked and are only ever used as offsets from here on.
struct GroupedCounts {
	std::vector<size_t> offsets; //!< n_rows + 1
	std::vector<size_t> col;
	std::vector<double> val;
};

//! Whether an all-zero feature column is an error.
//!
//! It is when FITTING: that feature contributes nothing to the likelihood, so its
//! parameters are driven entirely by the prior and the fit "succeeds" while the
//! feature's output is meaningless. It is not when PREDICTING or SCORING: there
//! the columns are not parameters, only weights on already-fitted rows of
//! `probs`, and a feature that happens to be absent from the supplied samples is
//! ordinary data. scikit-bio draws the same line.
enum class ColumnPolicy {
	kRejectAllZero,
	kAllowAllZero,
};

//! Validate one table's entries and regroup them by sample. `label` names the
//! table in error messages ("X" or "Y"). Assumes the caller has already checked
//! the parallel-array lengths and the declared shape, so `t.n_rows` and
//! `t.n_cols` can be trusted as sizes here.
//!
//! Grouping is a counting sort, so entries keep their input order within a
//! sample. That does not affect the result -- see the bit-identity argument on
//! ComputeSufficientStats -- but it does keep the duplicate-cell report
//! deterministic.
GroupedCounts GroupBySample(const SparseCounts &t, const std::string &label, ColumnPolicy columns) {
	const size_t nnz = t.val.size();
	const auto n_rows = static_cast<size_t>(t.n_rows);
	const auto n_cols = static_cast<size_t>(t.n_cols);

	GroupedCounts g;
	g.offsets.assign(n_rows + 1, 0);
	std::vector<double> row_sums(n_rows, 0.0);
	std::vector<double> col_sums(n_cols, 0.0);

	// Pass 1: bounds and values, accumulating the per-axis sums for the
	// all-zero checks and the per-sample counts for the regrouping.
	for (size_t k = 0; k < nnz; ++k) {
		const int64_t r = t.row[k];
		const int64_t c = t.col[k];
		if (r < 0 || r >= t.n_rows) {
			throw std::invalid_argument("mmvec: " + label + " row index " + std::to_string(r) + " out of range [0, " +
			                            std::to_string(t.n_rows) + ")");
		}
		if (c < 0 || c >= t.n_cols) {
			throw std::invalid_argument("mmvec: " + label + " column index " + std::to_string(c) +
			                            " out of range [0, " + std::to_string(t.n_cols) + ")");
		}
		const double v = t.val[k];
		// Finiteness first: -inf would satisfy the negativity test too, and
		// "not finite" is the more useful diagnosis.
		if (!std::isfinite(v)) {
			throw std::invalid_argument("mmvec: " + label + " value at entry " + std::to_string(k) +
			                            " is not finite (sample " + std::to_string(r) + ", feature " +
			                            std::to_string(c) + ")");
		}
		if (v < 0.0) {
			throw std::invalid_argument("mmvec: " + label + " value at entry " + std::to_string(k) + " is negative (" +
			                            std::to_string(v) + ", sample " + std::to_string(r) + ", feature " +
			                            std::to_string(c) + "); mmvec models counts");
		}
		g.offsets[static_cast<size_t>(r) + 1]++;
		row_sums[static_cast<size_t>(r)] += v;
		col_sums[static_cast<size_t>(c)] += v;
	}

	// Pass 2: turn the per-sample counts into offsets, then scatter.
	for (size_t s = 0; s < n_rows; ++s) {
		g.offsets[s + 1] += g.offsets[s];
	}
	g.col.resize(nnz);
	g.val.resize(nnz);
	std::vector<size_t> cursor(g.offsets.begin(), g.offsets.end() - 1);
	for (size_t k = 0; k < nnz; ++k) {
		const size_t pos = cursor[static_cast<size_t>(t.row[k])]++;
		g.col[pos] = static_cast<size_t>(t.col[k]);
		g.val[pos] = t.val[k];
	}

	// Pass 3: duplicate cells, via a per-feature stamp holding the last sample
	// that used it. O(nnz) with no hashing, and correct because the entries are
	// now visited in ascending sample order.
	std::vector<int64_t> last_sample(n_cols, -1);
	for (size_t s = 0; s < n_rows; ++s) {
		for (size_t k = g.offsets[s]; k < g.offsets[s + 1]; ++k) {
			if (last_sample[g.col[k]] == static_cast<int64_t>(s)) {
				throw std::invalid_argument("mmvec: " + label + " has a duplicate entry for sample " +
				                            std::to_string(s) + ", feature " + std::to_string(g.col[k]) +
				                            "; each (sample, feature) cell must appear at most once");
			}
			last_sample[g.col[k]] = static_cast<int64_t>(s);
		}
	}

	// All-zero rows and columns. The sums decide it, not index presence: a
	// feature that appears only with explicitly stored zeros is still all-zero.
	// An all-zero row is always fatal -- proportions are undefined for it either
	// way -- so only the column test is policy.
	const bool fitting = columns == ColumnPolicy::kRejectAllZero;
	const std::string row_remedy = fitting
	                                   ? "); remove all-zero rows and columns before fitting"
	                                   : "); every sample needs a nonzero total for its counts to become proportions";
	for (size_t s = 0; s < n_rows; ++s) {
		if (row_sums[s] == 0.0) {
			throw std::invalid_argument("mmvec: " + label + " contains an all-zero row (sample index " +
			                            std::to_string(s) + row_remedy);
		}
	}
	if (fitting) {
		for (size_t c = 0; c < n_cols; ++c) {
			if (col_sums[c] == 0.0) {
				throw std::invalid_argument("mmvec: " + label + " contains an all-zero column (feature index " +
				                            std::to_string(c) + "); remove all-zero rows and columns before fitting");
			}
		}
	}

	return g;
}

//! Index sample `s`'s entries in ascending feature order and return its total.
//!
//! The ordering is what makes `Predict` and `Score` independent of the order the
//! COO entries arrived in: both accumulate one term per stored cell, so unlike
//! `ComputeSufficientStats` their floating-point summation order is not pinned by
//! the sample index alone. The SQL layer's scan order must not change a result.
//!
//! The total is > 0 for any table that passed GroupBySample, which rejects
//! all-zero rows under either column policy -- so callers may divide by it.
double OrderRowByFeature(const GroupedCounts &g, int64_t s, std::vector<size_t> &order) {
	const size_t begin = g.offsets[static_cast<size_t>(s)];
	const size_t end = g.offsets[static_cast<size_t>(s) + 1];
	order.resize(end - begin);
	std::iota(order.begin(), order.end(), begin);
	std::sort(order.begin(), order.end(), [&g](size_t a, size_t b) { return g.col[a] < g.col[b]; });
	double total = 0.0;
	for (const size_t k : order) {
		total += g.val[k];
	}
	return total;
}

//! max|v|, the convergence diagnostic both optimizers report. Zero for an empty
//! vector, which cannot arise here -- a model always has parameters.
double MaxAbs(const std::vector<double> &v) {
	double worst = 0.0;
	for (const double x : v) {
		worst = std::fmax(worst, std::fabs(x));
	}
	return worst;
}

//! Subtract each row's mean from it, in place: logits -> ranks.
void CenterRowsInPlace(double *m, int64_t rows, int64_t cols) {
	for (int64_t i = 0; i < rows; ++i) {
		double *row = m + i * cols;
		double total = 0.0;
		for (int64_t j = 0; j < cols; ++j) {
			total += row[j];
		}
		const double mean = total / static_cast<double>(cols);
		for (int64_t j = 0; j < cols; ++j) {
			row[j] -= mean;
		}
	}
}

//! Replace each row by its softmax, in place: logits (or ranks) -> probs.
//!
//! Shifting by the row maximum is an algebraic identity for the softmax, and the
//! particular choice of the maximum is what keeps it safe: the largest exponent
//! becomes exactly 1, so no term can overflow and the denominator is >= 1, which
//! rules out a division by zero however extreme the logits are.
void SoftmaxRowsInPlace(double *m, int64_t rows, int64_t cols) {
	for (int64_t i = 0; i < rows; ++i) {
		double *row = m + i * cols;
		const double max = *std::max_element(row, row + cols);
		double total = 0.0;
		for (int64_t j = 0; j < cols; ++j) {
			row[j] = std::exp(row[j] - max);
			total += row[j];
		}
		for (int64_t j = 0; j < cols; ++j) {
			row[j] /= total;
		}
	}
}

//! Reject a table whose three parallel arrays disagree in length.
void CheckParallel(const SparseCounts &t, const std::string &label) {
	if (t.row.size() != t.col.size() || t.row.size() != t.val.size()) {
		throw std::invalid_argument("mmvec: " + label + " row, col and val arrays must all be the same length (got " +
		                            std::to_string(t.row.size()) + ", " + std::to_string(t.col.size()) + ", " +
		                            std::to_string(t.val.size()) + ")");
	}
}

// ---------------------------------------------------------------------------
// The objective
// ---------------------------------------------------------------------------

//! ELEMENT offsets of the four blocks inside a packed theta or gradient vector --
//! indices into a `double*`, not byte offsets.
struct ParamLayout {
	int64_t d1 = 0;   //!< X features
	int64_t d2 = 0;   //!< Y features, including the reference
	int64_t nref = 0; //!< d2 - 1, the number of parameterized Y features
	int64_t p = 0;    //!< latent dimension
	size_t x_main = 0;
	size_t x_bias = 0;
	size_t y_main = 0;
	size_t y_bias = 0;
	size_t total = 0;
};

ParamLayout Layout(const ModelShape &shape) {
	ParamLayout l;
	l.d1 = shape.n_features_x;
	l.d2 = shape.n_features_y;
	l.nref = shape.n_features_y - 1;
	l.p = shape.n_components;
	l.x_main = 0;
	l.x_bias = l.x_main + static_cast<size_t>(l.d1 * l.p);
	l.y_main = l.x_bias + static_cast<size_t>(l.d1);
	l.y_bias = l.y_main + static_cast<size_t>(l.p * l.nref);
	l.total = l.y_bias + static_cast<size_t>(l.nref);
	return l;
}

//! A row-major count matrix the objective reads observed Y counts out of, with
//! an optional row indirection. Element (r, j) of the NON-REFERENCE block is
//! `base[Row(r) * stride + j + 1]` -- the `+ 1` skips the reference column.
//!
//! The two paths differ only here: the full-batch objective reads `y_sums`
//! directly (one row per X feature, no indirection), while the minibatch reads
//! the dense Y table through the batch's sample indices.
struct ObsView {
	const double *base = nullptr;
	int64_t stride = 0;
	const int64_t *row_map = nullptr; //!< nullptr means the identity map
};

//! Resolve an optional index map, where nullptr means identity.
inline int64_t MapRow(const int64_t *map, int64_t r) {
	return map != nullptr ? map[r] : r;
}

//! Validate a shape on its own terms. Split out from the theta check because a
//! degenerate shape must be rejected BEFORE NumParams is used as an allocation
//! size (n_features_y = 0 would make it negative).
void CheckShape(const ModelShape &shape) {
	if (shape.n_features_x < 1) {
		throw std::invalid_argument("mmvec: n_features_x must be >= 1 (got " + std::to_string(shape.n_features_x) +
		                            ")");
	}
	if (shape.n_features_y < 2) {
		throw std::invalid_argument("mmvec: n_features_y must be >= 2 (got " + std::to_string(shape.n_features_y) +
		                            ")");
	}
	if (shape.n_components < 1) {
		throw std::invalid_argument("mmvec: n_components must be >= 1 (got " + std::to_string(shape.n_components) +
		                            ")");
	}
}

//! Validate a shape together with the theta it will be used with.
//!
//! Length before finiteness, deliberately: a theta of the wrong length is not
//! this model's theta at all, so reporting which of its elements is a NaN would
//! name an index that means nothing.
void CheckShapeAndTheta(const ModelShape &shape, const std::vector<double> &theta) {
	CheckShape(shape);
	const size_t want = static_cast<size_t>(NumParams(shape));
	if (theta.size() != want) {
		throw std::invalid_argument("mmvec: theta has " + std::to_string(theta.size()) + " parameters, expected " +
		                            std::to_string(want) + " for d1=" + std::to_string(shape.n_features_x) + ", d2=" +
		                            std::to_string(shape.n_features_y) + ", p=" + std::to_string(shape.n_components));
	}
	// A single non-finite parameter poisons every loss and gradient downstream
	// without any comparison failing, which reads as a successful all-NaN fit.
	for (size_t i = 0; i < theta.size(); ++i) {
		if (!std::isfinite(theta[i])) {
			throw std::invalid_argument("mmvec: theta[" + std::to_string(i) +
			                            "] is not finite; the objective cannot be evaluated there");
		}
	}
}

//! scikit-bio's fit-time check: a prior scale of zero or less is not a prior.
void CheckPriors(const Priors &priors) {
	if (!(priors.x_prior_scale > 0.0)) {
		throw std::invalid_argument("mmvec: x_prior_scale must be > 0 (got " + std::to_string(priors.x_prior_scale) +
		                            ")");
	}
	if (!(priors.y_prior_scale > 0.0)) {
		throw std::invalid_argument("mmvec: y_prior_scale must be > 0 (got " + std::to_string(priors.y_prior_scale) +
		                            ")");
	}
}

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

//! The Y row total for every sample: the sum over ALL d2 features, reference
//! included. Loop-invariant across Adam updates, so Adam computes this once.
std::vector<double> ComputeSampleTotals(const MinibatchInputs &inputs, int64_t n_features_y) {
	std::vector<double> totals(static_cast<size_t>(inputs.n_samples), 0.0);
	for (int64_t s = 0; s < inputs.n_samples; ++s) {
		const double *y_row = inputs.y_dense.data() + s * n_features_y;
		double total = 0.0;
		for (int64_t j = 0; j < n_features_y; ++j) {
			total += y_row[j];
		}
		totals[static_cast<size_t>(s)] = total;
	}
	return totals;
}

//! Resolve a minibatch into the two per-row index maps the objective needs and the
//! Y row totals it conditions on, validating every index on the way -- so the hot
//! loop inside EvaluateObjective does not have to.
//!
//! `sample_totals` is the precomputed total for EVERY sample when the caller has it.
//! Adam does, because the totals do not change between updates: re-deriving them per
//! update costs 4.5 billion additions over a default cystic-fibrosis fit to produce
//! only n_samples = 172 distinct values. Pass nullptr to sum each named row on the
//! spot, which is what a one-shot caller wants -- precomputing every sample would be
//! wasted work when the batch names only a few. Either way the sum runs over
//! ascending j, so the two paths agree bit-for-bit.
void ResolveBatch(const ModelShape &shape, const MinibatchInputs &inputs, const std::vector<int64_t> &batch,
                  const double *sample_totals, Workspace &ws) {
	const auto nnz = static_cast<int64_t>(inputs.x_row.size());
	for (size_t b = 0; b < batch.size(); ++b) {
		const int64_t e = batch[b];
		if (e < 0 || e >= nnz) {
			throw std::invalid_argument("mmvec: minibatch index " + std::to_string(e) + " out of range [0, " +
			                            std::to_string(nnz) + ")");
		}
		const int64_t s = inputs.x_row[static_cast<size_t>(e)];
		const int64_t i = inputs.x_col[static_cast<size_t>(e)];
		if (s < 0 || s >= inputs.n_samples) {
			throw std::invalid_argument("mmvec: minibatch entry " + std::to_string(e) + " names sample " +
			                            std::to_string(s) + ", out of range [0, " + std::to_string(inputs.n_samples) +
			                            ")");
		}
		if (i < 0 || i >= shape.n_features_x) {
			throw std::invalid_argument("mmvec: minibatch entry " + std::to_string(e) + " names X feature " +
			                            std::to_string(i) + ", out of range [0, " + std::to_string(shape.n_features_x) +
			                            ")");
		}
		ws.sample_index[b] = s;
		ws.x_index[b] = i;

		if (sample_totals != nullptr) {
			ws.row_totals[b] = sample_totals[static_cast<size_t>(s)];
			continue;
		}
		double total = 0.0;
		const double *y_row = inputs.y_dense.data() + s * shape.n_features_y;
		for (int64_t j = 0; j < shape.n_features_y; ++j) {
			total += y_row[j];
		}
		ws.row_totals[b] = total;
	}
}

//! The full-batch objective packaged for LBFGS++, which needs a callable taking
//! (x, grad) and returning the value.
//!
//! It carries two things the solver cannot give us:
//!
//!   * the evaluation trace, appended on every call because scikit-bio's
//!     `loss_curve_` counts evaluations, line-search probes included;
//!   * a best-so-far (theta, loss) snapshot, which is not defensive padding.
//!     LBFGS++ throws on line-search failure and may leave `x` mid-step, and the
//!     reference implementation's own runs terminate abnormally on real data, so
//!     without this a failed line search would lose the fit entirely.
class LbfgsObjective {
public:
	LbfgsObjective(const ParamLayout &l, const Priors &priors, const SufficientStats &stats, Workspace &ws)
	    : l_(l), priors_(priors), obs_ {stats.y_sums.data(), l.d2, nullptr}, totals_(stats.n_sums.data()), ws_(ws) {
	}

	double operator()(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
		const double loss =
		    EvaluateObjective(l_, priors_, l_.d1, nullptr, obs_, totals_, 1.0, x.data(), ws_, grad.data());
		loss_curve.push_back(loss);
		if (loss < best_loss) {
			best_loss = loss;
			best_theta = x;
		}
		return loss;
	}

	std::vector<double> loss_curve;
	double best_loss = std::numeric_limits<double>::infinity();
	Eigen::VectorXd best_theta;

private:
	const ParamLayout &l_;
	const Priors &priors_;
	ObsView obs_;
	const double *totals_;
	Workspace &ws_;
};

} // namespace

int64_t NumParams(const ModelShape &shape) {
	const int64_t nref = shape.n_features_y - 1;
	return shape.n_features_x * shape.n_components + shape.n_features_x + shape.n_components * nref + nref;
}

void Workspace::Resize(const ModelShape &shape, int64_t n_rows) {
	const auto rows = static_cast<size_t>(n_rows);
	const auto nref = static_cast<size_t>(shape.n_features_y - 1);
	const auto p = static_cast<size_t>(shape.n_components);
	logits.resize(rows * nref);
	resids.resize(rows * nref);
	row_totals.resize(rows);
	x_main_rows.resize(rows * p);
	dx_main_rows.resize(rows * p);
	x_index.resize(rows);
	sample_index.resize(rows);
}

std::vector<double> ComputeLogits(const ModelShape &shape, const std::vector<double> &theta) {
	CheckShapeAndTheta(shape, theta);
	const ParamLayout l = Layout(shape);

	// The same equation the objective uses -- one row per X feature, identity
	// index map -- then padded out with the pinned zero reference column. Routed
	// through a scratch buffer rather than restating the model: this runs once
	// per fit (for `ranks` and `probs`), not once per iteration, so the extra
	// (d1 x d2-1) copy costs nothing worth trading a second copy of the model
	// equation for.
	// Sized by hand rather than through Workspace::Resize, which allocates all of
	// the objective's buffers. FillNonRefLogits touches exactly these two, and
	// `resids` alone is d1 x (d2-1) -- 10 MB at cystic-fibrosis scale, zero-filled
	// and never read, on every ranks / probs / predict / score call.
	Workspace ws;
	ws.logits.resize(static_cast<size_t>(l.d1 * l.nref));
	ws.x_main_rows.resize(static_cast<size_t>(l.d1 * l.p));
	FillNonRefLogits(l, theta.data(), nullptr, l.d1, ws, ws.logits.data());

	std::vector<double> out(static_cast<size_t>(l.d1 * l.d2), 0.0);
	for (int64_t i = 0; i < l.d1; ++i) {
		out[static_cast<size_t>(i * l.d2)] = 0.0; // Y feature 0 is the reference
		for (int64_t j = 0; j < l.nref; ++j) {
			const double v = ws.logits[static_cast<size_t>(i * l.nref + j)];
			// A finite theta can still produce an infinite logit: the dot product
			// overflows once the parameters reach about 1e154. The softmax then
			// computes exp(inf - inf) and every probability in the row comes back
			// NaN -- silently, since nothing downstream compares against anything.
			// A fit cannot wander there (the Gaussian prior makes it astronomically
			// expensive), but the SQL layer rebuilds theta from a model RELATION
			// whose values a user can write by hand, so it is reachable and has to
			// fail loudly. Checked here rather than in CheckShapeAndTheta because
			// the parameters themselves are finite; it is the product that is not.
			// This is the read path -- the objective never calls ComputeLogits --
			// so the scan costs nothing per iteration.
			if (!std::isfinite(v)) {
				throw std::invalid_argument("mmvec: the logit for X feature " + std::to_string(i) + ", Y feature " +
				                            std::to_string(j + 1) +
				                            " is not finite; the parameters are finite but their product overflowed");
			}
			out[static_cast<size_t>(i * l.d2 + j + 1)] = v;
		}
	}
	return out;
}

double FullBatchLossAndGradient(const ModelShape &shape, const Priors &priors, const SufficientStats &stats,
                                const std::vector<double> &theta, Workspace &ws, std::vector<double> &grad) {
	CheckShapeAndTheta(shape, theta);
	CheckPriors(priors);
	if (stats.n_features_x != shape.n_features_x || stats.n_features_y != shape.n_features_y) {
		throw std::invalid_argument("mmvec: sufficient statistics are d1=" + std::to_string(stats.n_features_x) +
		                            ", d2=" + std::to_string(stats.n_features_y) + " but the model shape is d1=" +
		                            std::to_string(shape.n_features_x) + ", d2=" + std::to_string(shape.n_features_y));
	}

	const ParamLayout l = Layout(shape);
	ws.Resize(shape, l.d1);
	grad.resize(l.total);

	// One row per X feature, reading y_sums directly: no index indirection, and
	// the totals are already grouped as n_sums. norm is 1 -- this IS the
	// objective the minibatch version estimates.
	const ObsView obs {stats.y_sums.data(), l.d2, nullptr};
	return EvaluateObjective(l, priors, l.d1, nullptr, obs, stats.n_sums.data(), 1.0, theta.data(), ws, grad.data());
}

double BatchNormFactor(BatchNorm mode, double x_total, int64_t n_samples, int64_t batch_size) {
	if (batch_size < 1) {
		throw std::invalid_argument("mmvec: batch_size must be >= 1 (got " + std::to_string(batch_size) + ")");
	}
	if (n_samples < 1) {
		throw std::invalid_argument("mmvec: n_samples must be >= 1 (got " + std::to_string(n_samples) + ")");
	}
	if (!(x_total > 0.0)) {
		throw std::invalid_argument("mmvec: total X count must be > 0 (got " + std::to_string(x_total) + ")");
	}
	const double numerator = mode == BatchNorm::Unbiased ? x_total : static_cast<double>(n_samples);
	return numerator / static_cast<double>(batch_size);
}

double MinibatchLossAndGradient(const ModelShape &shape, const Priors &priors, const MinibatchInputs &inputs,
                                const std::vector<int64_t> &batch, double norm, const std::vector<double> &theta,
                                Workspace &ws, std::vector<double> &grad) {
	CheckShapeAndTheta(shape, theta);
	CheckPriors(priors);
	if (inputs.x_row.size() != inputs.x_col.size() || inputs.x_row.size() != inputs.x_val.size()) {
		throw std::invalid_argument("mmvec: minibatch x_row, x_col and x_val must all be the same length (got " +
		                            std::to_string(inputs.x_row.size()) + ", " + std::to_string(inputs.x_col.size()) +
		                            ", " + std::to_string(inputs.x_val.size()) + ")");
	}
	if (inputs.n_samples < 1) {
		throw std::invalid_argument("mmvec: minibatch requires at least one sample");
	}
	if (inputs.y_dense.size() != static_cast<size_t>(inputs.n_samples * shape.n_features_y)) {
		throw std::invalid_argument("mmvec: minibatch y_dense has " + std::to_string(inputs.y_dense.size()) +
		                            " values, expected " + std::to_string(inputs.n_samples * shape.n_features_y) +
		                            " for " + std::to_string(inputs.n_samples) + " samples x " +
		                            std::to_string(shape.n_features_y) + " features");
	}
	if (batch.empty()) {
		throw std::invalid_argument("mmvec: minibatch is empty");
	}

	const ParamLayout l = Layout(shape);
	const auto n_rows = static_cast<int64_t>(batch.size());
	ws.Resize(shape, n_rows);
	grad.resize(l.total);

	ResolveBatch(shape, inputs, batch, nullptr, ws);

	const ObsView obs {inputs.y_dense.data(), shape.n_features_y, ws.sample_index.data()};
	return EvaluateObjective(l, priors, n_rows, ws.x_index.data(), obs, ws.row_totals.data(), norm, theta.data(), ws,
	                         grad.data());
}

Model FitLbfgsFromInit(const ModelShape &shape, const Priors &priors, const SufficientStats &stats,
                       const std::vector<double> &theta0, int64_t max_iter) {
	if (max_iter < 1) {
		throw std::invalid_argument("mmvec: max_iter must be >= 1 (got " + std::to_string(max_iter) + ")");
	}
	// Validate through the objective so the rules live in one place; this also
	// means the checks happen once here rather than on every evaluation.
	Workspace ws;
	std::vector<double> probe_grad;
	const double initial_loss = FullBatchLossAndGradient(shape, priors, stats, theta0, ws, probe_grad);

	const ParamLayout l = Layout(shape);
	LBFGSpp::LBFGSParam<double> param;
	// The mapping onto scikit-bio's scipy.optimize.minimize(..., method='L-BFGS-B',
	// options={'maxiter': max_iter, 'ftol': 1e-9, 'gtol': 1e-5}).
	param.m = kLbfgsHistory;
	param.epsilon = kLbfgsGtol; // NOTE: LBFGS++ tests ||g||_2 here where SciPy
	                            // tests ||g||_inf. Since ||g||_2 >= ||g||_inf this
	                            // can only stop later, never earlier, than SciPy
	                            // would; `converged` below is decided on the
	                            // L-infinity norm regardless.
	param.epsilon_rel = 0.0;    // SciPy has no relative gradient test.
	param.past = 1;             // With delta below this is exactly SciPy's ftol
	param.delta = kLbfgsFtol;   // test: |f_{k-1} - f_k| <= delta*max(|f_k|,|f_{k-1}|,1).
	param.max_iterations = static_cast<int>(max_iter);
	// Line-search constants matching SciPy's dcsrch. param.ftol and param.wolfe
	// are the Armijo and Wolfe parameters -- NOT the convergence ftol above.
	param.max_linesearch = 20;
	param.ftol = 1e-3;
	param.wolfe = 0.9;

	LBFGSpp::LBFGSSolver<double, LBFGSpp::LineSearchMoreThuente> solver(param);
	LbfgsObjective objective(l, priors, stats, ws);
	Eigen::VectorXd x = Eigen::Map<const Eigen::VectorXd>(theta0.data(), static_cast<Eigen::Index>(theta0.size()));
	double fx = initial_loss;

	Model model;
	model.shape = shape;
	std::string outcome;
	int niter = 0;
	bool niter_known = true;
	try {
		niter = solver.minimize(objective, x, fx);
		// LBFGS++ returns the same iteration count whether a convergence test
		// fired or the cap was reached, so they are told apart by the count. On
		// the exact boundary the convergence tests do run first, so this
		// attributes a boundary stop to the cap -- under-reporting convergence
		// rather than over-reporting it.
		if (niter >= static_cast<int>(max_iter)) {
			model.converged = false;
			outcome = "reached the iteration limit of " + std::to_string(max_iter);
		} else {
			model.converged = true;
			outcome = "converged";
		}
	} catch (const std::exception &e) {
		// LBFGS++ reports line-search failure by throwing, where SciPy returns a
		// status. The snapshot inside the functor is what makes this recoverable;
		// the solver's own iterate is not trustworthy after a throw.
		//
		// The iteration count is genuinely LOST on this path: `minimize` returns it,
		// so throwing leaves `niter` at its initial 0 no matter how many iterations
		// actually ran, and LBFGS++ exposes no other accessor for it. Reporting a
		// confident "0 iterations" for a fit that managed forty-seven would be a
		// fabrication, so the message says the count is unavailable instead.
		// The evaluation count is unaffected -- the functor counts those itself.
		model.converged = false;
		niter_known = false;
		outcome = std::string("line search failed (") + e.what() + "); the best point seen was kept";
	}

	model.theta.assign(objective.best_theta.data(), objective.best_theta.data() + objective.best_theta.size());
	model.loss_curve = std::move(objective.loss_curve);
	model.n_iter = niter;

	// Recompute at the point actually being returned, so max|gradient| describes
	// `theta` and not some iterate that was discarded.
	std::vector<double> grad;
	model.final_loss = FullBatchLossAndGradient(shape, priors, stats, model.theta, ws, grad);
	model.max_abs_grad = MaxAbs(grad);

	model.message = outcome + "; " +
	                (niter_known ? std::to_string(niter) + " iterations" : std::string("iteration count unavailable")) +
	                ", " + std::to_string(model.loss_curve.size()) +
	                " objective evaluations, max|gradient| = " + std::to_string(model.max_abs_grad);
	return model;
}

Rng::Rng(uint64_t seed) : engine_(seed) {
}

double Rng::Uniform() {
	// Top 53 bits scaled by 2^-53: the canonical mt19937_64 -> [0,1) mapping, and
	// the only step here that <random> would have made implementation-defined.
	return static_cast<double>(engine_() >> 11) * 0x1.0p-53;
}

double Rng::Normal() {
	if (has_cached_normal_) {
		has_cached_normal_ = false;
		return cached_normal_;
	}
	// Marsaglia polar method: rejection-sample a point in the open unit disc, then
	// map it to two independent standard normals.
	double u = 0.0;
	double v = 0.0;
	double s = 0.0;
	do {
		u = 2.0 * Uniform() - 1.0;
		v = 2.0 * Uniform() - 1.0;
		s = u * u + v * v;
	} while (s >= 1.0 || s == 0.0);
	const double factor = std::sqrt(-2.0 * std::log(s) / s);
	cached_normal_ = v * factor;
	has_cached_normal_ = true;
	return u * factor;
}

std::vector<int64_t> SampleWeightedIndices(const std::vector<double> &weights, int64_t count, Rng &rng) {
	if (weights.empty()) {
		throw std::invalid_argument("mmvec: cannot sample from an empty weight vector");
	}
	if (count < 0) {
		throw std::invalid_argument("mmvec: sample count must be >= 0 (got " + std::to_string(count) + ")");
	}
	std::vector<double> cdf(weights.size());
	double total = 0.0;
	for (size_t i = 0; i < weights.size(); ++i) {
		const double w = weights[i];
		if (!std::isfinite(w)) {
			throw std::invalid_argument("mmvec: weight " + std::to_string(i) + " is not finite");
		}
		if (w < 0.0) {
			throw std::invalid_argument("mmvec: weight " + std::to_string(i) + " is negative (" + std::to_string(w) +
			                            ")");
		}
		total += w;
		cdf[i] = total;
	}
	if (!(total > 0.0)) {
		throw std::invalid_argument("mmvec: weights sum to zero, so nothing can be sampled");
	}

	std::vector<int64_t> out(static_cast<size_t>(count));
	for (int64_t k = 0; k < count; ++k) {
		const double target = rng.Uniform() * total;
		// upper_bound, not lower_bound: a zero-weight entry leaves the cumulative
		// sum unchanged and must never be selectable.
		const auto it = std::upper_bound(cdf.begin(), cdf.end(), target);
		auto idx = static_cast<size_t>(std::distance(cdf.begin(), it));
		if (idx >= cdf.size()) {
			idx = cdf.size() - 1; // rounding guard: target can reach `total`
		}
		out[static_cast<size_t>(k)] = static_cast<int64_t>(idx);
	}
	return out;
}

std::vector<double> InitTheta(const ModelShape &shape, uint64_t seed) {
	CheckShape(shape);
	std::vector<double> theta(static_cast<size_t>(NumParams(shape)));
	Rng rng(seed);
	for (double &t : theta) {
		t = rng.Normal();
	}
	return theta;
}

Model FitAdamWithIndices(const ModelShape &shape, const Priors &priors, const MinibatchInputs &inputs,
                         const AdamParams &params, const std::vector<double> &theta0,
                         const std::vector<int64_t> &batches) {
	if (params.batch_size < 1) {
		throw std::invalid_argument("mmvec: batch_size must be >= 1 (got " + std::to_string(params.batch_size) + ")");
	}
	if (batches.empty()) {
		throw std::invalid_argument("mmvec: no minibatches supplied");
	}
	if (static_cast<int64_t>(batches.size()) % params.batch_size != 0) {
		throw std::invalid_argument("mmvec: the minibatch index sequence holds " + std::to_string(batches.size()) +
		                            " entries, which is not a whole number of batches of " +
		                            std::to_string(params.batch_size));
	}
	if (!(params.learning_rate > 0.0)) {
		throw std::invalid_argument("mmvec: learning_rate must be > 0 (got " + std::to_string(params.learning_rate) +
		                            ")");
	}
	if (!(params.beta_1 >= 0.0 && params.beta_1 < 1.0) || !(params.beta_2 >= 0.0 && params.beta_2 < 1.0)) {
		throw std::invalid_argument("mmvec: beta_1 and beta_2 must lie in [0, 1)");
	}
	if (!(params.clipnorm > 0.0)) {
		throw std::invalid_argument("mmvec: clipnorm must be > 0 (got " + std::to_string(params.clipnorm) + ")");
	}
	double x_total = 0.0;
	for (const double v : inputs.x_val) {
		x_total += v;
	}
	const double norm = BatchNormFactor(params.batch_norm, x_total, inputs.n_samples, params.batch_size);

	const int64_t n_updates = static_cast<int64_t>(batches.size()) / params.batch_size;

	// The summed-away statistics, built BEFORE the training loop rather than after
	// it. They are needed only at the end, for the full-batch `final_loss`, but
	// they also carry the degenerate-data validation -- all-zero rows and columns,
	// duplicate cells -- that the minibatch objective never performs. Computing
	// them last meant a table L-BFGS rejects outright could run every Adam update
	// to completion and only then throw, destroying a fully-trained model on the
	// way out. Building them first makes the two optimizers fail on the same
	// inputs at the same point, at no extra cost: it is the same single call,
	// moved.
	SparseCounts x_coo;
	x_coo.n_rows = inputs.n_samples;
	x_coo.n_cols = shape.n_features_x;
	x_coo.row = inputs.x_row;
	x_coo.col = inputs.x_col;
	x_coo.val = inputs.x_val;
	SparseCounts y_coo;
	y_coo.n_rows = inputs.n_samples;
	y_coo.n_cols = shape.n_features_y;
	for (int64_t s = 0; s < inputs.n_samples; ++s) {
		for (int64_t j = 0; j < shape.n_features_y; ++j) {
			const double v = inputs.y_dense[static_cast<size_t>(s * shape.n_features_y + j)];
			if (v != 0.0) {
				y_coo.row.push_back(s);
				y_coo.col.push_back(j);
				y_coo.val.push_back(v);
			}
		}
	}
	const SufficientStats stats = ComputeSufficientStats(x_coo, y_coo);

	Model model;
	model.shape = shape;
	model.theta = theta0;
	model.loss_curve.reserve(static_cast<size_t>(n_updates));

	// Adam state. One flat pair of moment vectors rather than four: the update is
	// elementwise and every block shares the same hyperparameters, so a single
	// pass over the packed vector is exactly the four per-block updates the
	// reference implementation performs.
	const auto n_params = static_cast<size_t>(NumParams(shape));
	std::vector<double> moment1(n_params, 0.0);
	std::vector<double> moment2(n_params, 0.0);
	constexpr double kAdamEps = 1e-8;

	// Two things the minibatch objective re-derives on every call are invariant across
	// updates, and at cystic-fibrosis defaults each is billions of operations: the Y
	// row totals (4.5e9 additions to produce only n_samples = 172 distinct values) and
	// the validation sweep over theta (2.5e9 isfinite checks). Both are done once.
	//
	// Update 1 still goes through the fully-validating public entry point, so every
	// input error is raised in the same place with the same message as before; what
	// updates 2..N skip is only the REPEAT of work whose answer cannot have changed.
	// Theta itself does change, so the loss is checked for finiteness below -- an O(1)
	// test that catches a diverged fit where the O(n_params) sweep used to.
	const ParamLayout l = Layout(shape);
	const std::vector<double> sample_totals = ComputeSampleTotals(inputs, shape.n_features_y);
	const ObsView obs {inputs.y_dense.data(), shape.n_features_y, nullptr};

	Workspace ws;
	std::vector<double> grad;
	std::vector<int64_t> batch(static_cast<size_t>(params.batch_size));
	for (int64_t update = 1; update <= n_updates; ++update) {
		const auto begin = batches.begin() + static_cast<std::ptrdiff_t>((update - 1) * params.batch_size);
		batch.assign(begin, begin + static_cast<std::ptrdiff_t>(params.batch_size));

		double loss = 0.0;
		if (update == 1) {
			loss = MinibatchLossAndGradient(shape, priors, inputs, batch, norm, model.theta, ws, grad);
		} else {
			ResolveBatch(shape, inputs, batch, sample_totals.data(), ws);
			const ObsView batch_obs {obs.base, obs.stride, ws.sample_index.data()};
			loss = EvaluateObjective(l, priors, params.batch_size, ws.x_index.data(), batch_obs, ws.row_totals.data(),
			                         norm, model.theta.data(), ws, grad.data());
		}
		if (!std::isfinite(loss)) {
			throw std::invalid_argument("mmvec: Adam diverged -- the minibatch loss became non-finite at update " +
			                            std::to_string(update) +
			                            ". Lower learning_rate, or tighten clipnorm, and refit.");
		}
		model.loss_curve.push_back(loss);

		// Global L2 clipping, over all four blocks at once and AFTER the priors are
		// already in the gradient.
		double sq = 0.0;
		for (const double g : grad) {
			sq += g * g;
		}
		const double global_norm = std::sqrt(sq);
		if (global_norm > params.clipnorm) {
			const double scale = params.clipnorm / global_norm;
			for (double &g : grad) {
				g *= scale;
			}
		}

		const double bc1 = 1.0 - std::pow(params.beta_1, static_cast<double>(update));
		const double bc2 = 1.0 - std::pow(params.beta_2, static_cast<double>(update));
		for (size_t i = 0; i < n_params; ++i) {
			moment1[i] = params.beta_1 * moment1[i] + (1.0 - params.beta_1) * grad[i];
			moment2[i] = params.beta_2 * moment2[i] + (1.0 - params.beta_2) * grad[i] * grad[i];
			const double m_hat = moment1[i] / bc1;
			const double v_hat = moment2[i] / bc2;
			model.theta[i] -= params.learning_rate * m_hat / (std::sqrt(v_hat) + kAdamEps);
		}
	}

	model.n_iter = n_updates;

	// The full-batch objective at the final parameters, so final_loss means the
	// same thing as it does for L-BFGS. The last minibatch loss is a noisy
	// estimate over one draw and stays in loss_curve where it belongs. `stats` was
	// built before the loop; see the note there for why.
	std::vector<double> full_grad;
	model.final_loss = FullBatchLossAndGradient(shape, priors, stats, model.theta, ws, full_grad);
	model.max_abs_grad = MaxAbs(full_grad);

	// Adam has no convergence test at all -- it runs its schedule and stops -- so
	// claiming convergence would be a fabrication.
	model.converged = false;
	model.message = "ran " + std::to_string(n_updates) +
	                " Adam updates to completion (Adam has no convergence test); max|gradient| = " +
	                std::to_string(model.max_abs_grad);
	return model;
}

Model FitAdam(const ModelShape &shape, const Priors &priors, const MinibatchInputs &inputs, const AdamParams &params,
              const std::vector<double> &theta0, int64_t epochs, uint64_t seed) {
	if (epochs < 1) {
		throw std::invalid_argument("mmvec: epochs must be >= 1 (got " + std::to_string(epochs) + ")");
	}
	if (params.batch_size < 1) {
		throw std::invalid_argument("mmvec: batch_size must be >= 1 (got " + std::to_string(params.batch_size) + ")");
	}
	if (inputs.x_val.empty()) {
		throw std::invalid_argument("mmvec: Adam needs X counts to weight its sampling");
	}

	// One epoch is a pass over the nonzero cells, so the update count follows from
	// the data rather than from `epochs` alone.
	const auto nnz = static_cast<int64_t>(inputs.x_val.size());
	const int64_t per_epoch = std::max<int64_t>(1, nnz / params.batch_size);
	const int64_t n_updates = epochs * per_epoch;

	// The X counts are the sampling weights as they stand. SampleWeightedIndices
	// builds its own cumulative sums and draws Uniform() * total, so it is
	// scale-invariant by construction -- pre-dividing by the total would allocate
	// an nnz-sized copy to compute a distribution it recomputes anyway.
	Rng rng(seed);
	const std::vector<int64_t> batches = SampleWeightedIndices(inputs.x_val, n_updates * params.batch_size, rng);
	return FitAdamWithIndices(shape, priors, inputs, params, theta0, batches);
}

SufficientStats ComputeSufficientStats(const SparseCounts &x, const SparseCounts &y) {
	CheckParallel(x, "X");
	CheckParallel(y, "Y");

	// Shape reconciliation comes before the per-table passes so that a sample
	// count that disagrees between the tables is reported as such, rather than
	// as a spurious all-zero row in whichever table claims the extra sample.
	if (x.n_rows != y.n_rows) {
		throw std::invalid_argument("mmvec: X and Y must have the same number of samples (got " +
		                            std::to_string(x.n_rows) + " and " + std::to_string(y.n_rows) + ")");
	}
	if (x.n_rows < 1) {
		throw std::invalid_argument("mmvec: requires at least one sample");
	}
	if (x.n_cols < 1) {
		throw std::invalid_argument("mmvec: X must have at least one feature");
	}
	if (y.n_cols < 2) {
		throw std::invalid_argument("mmvec: Y must have at least two features (got " + std::to_string(y.n_cols) +
		                            "); feature 0 is the reference category and carries a fixed zero logit, so a "
		                            "single feature leaves no likelihood to fit");
	}

	const auto gx = GroupBySample(x, "X", ColumnPolicy::kRejectAllZero);
	const auto gy = GroupBySample(y, "Y", ColumnPolicy::kRejectAllZero);

	SufficientStats stats;
	stats.n_features_x = x.n_cols;
	stats.n_features_y = y.n_cols;
	const auto d1 = static_cast<size_t>(x.n_cols);
	const auto d2 = static_cast<size_t>(y.n_cols);
	stats.y_sums.assign(d1 * d2, 0.0);

	// Per-sample sparse outer product: neither table is ever densified, and the
	// cost is sum_n nnz_x(n) * nnz_y(n) rather than n_samples * d1 * d2.
	//
	// Samples are visited in ascending order and each (i, j) cell takes exactly
	// one term per sample, so the summation order -- and therefore the result,
	// bit for bit -- is fixed however the caller ordered the COO arrays.
	for (size_t s = 0; s < static_cast<size_t>(x.n_rows); ++s) {
		const size_t yb = gy.offsets[s];
		const size_t ye = gy.offsets[s + 1];
		for (size_t a = gx.offsets[s]; a < gx.offsets[s + 1]; ++a) {
			const double xv = gx.val[a];
			double *out = stats.y_sums.data() + gx.col[a] * d2;
			for (size_t b = yb; b < ye; ++b) {
				out[gy.col[b]] += xv * gy.val[b];
			}
		}
	}

	// n_sums as the row sums of y_sums rather than as a second aggregation over
	// the data, so the loss term dot(n_sums, log_norm) and the gradient term
	// y_sums - n_sums * probs stay exactly consistent. See the header.
	stats.n_sums.assign(d1, 0.0);
	for (size_t i = 0; i < d1; ++i) {
		double row = 0.0;
		for (size_t j = 0; j < d2; ++j) {
			row += stats.y_sums[i * d2 + j];
		}
		stats.n_sums[i] = row;
	}

	return stats;
}

// ---------------------------------------------------------------------------
// Derived quantities
// ---------------------------------------------------------------------------

std::vector<double> Ranks(const ModelShape &shape, const std::vector<double> &theta) {
	std::vector<double> out = ComputeLogits(shape, theta);
	CenterRowsInPlace(out.data(), shape.n_features_x, shape.n_features_y);
	return out;
}

std::vector<double> Probs(const ModelShape &shape, const std::vector<double> &theta) {
	std::vector<double> out = ComputeLogits(shape, theta);
	SoftmaxRowsInPlace(out.data(), shape.n_features_x, shape.n_features_y);
	return out;
}

std::vector<double> Predict(const ModelShape &shape, const std::vector<double> &theta, const SparseCounts &x) {
	CheckShapeAndTheta(shape, theta);
	CheckParallel(x, "X");
	if (x.n_rows < 1) {
		throw std::invalid_argument("mmvec: requires at least one sample to predict");
	}
	// Checked, where scikit-bio does not: X is consumed positionally, so a table
	// with a different feature count is silently a different model's X upstream.
	if (x.n_cols != shape.n_features_x) {
		throw std::invalid_argument("mmvec: X has " + std::to_string(x.n_cols) +
		                            " features but the model was fitted on " + std::to_string(shape.n_features_x));
	}
	const auto gx = GroupBySample(x, "X", ColumnPolicy::kAllowAllZero);
	const std::vector<double> probs = Probs(shape, theta);

	// rowprops(x) @ probs, accumulated over x's stored cells only: O(nnz*d2)
	// rather than the O(n*d1*d2) a dense product would cost, which matters
	// because microbiome X tables are mostly zeros.
	const int64_t d2 = shape.n_features_y;
	std::vector<double> out(static_cast<size_t>(x.n_rows * d2), 0.0);
	std::vector<size_t> order;
	for (int64_t s = 0; s < x.n_rows; ++s) {
		const double total = OrderRowByFeature(gx, s, order);
		double *out_row = &out[static_cast<size_t>(s * d2)];
		for (const size_t k : order) {
			const double weight = gx.val[k] / total;
			const double *probs_row = &probs[gx.col[k] * static_cast<size_t>(d2)];
			for (int64_t j = 0; j < d2; ++j) {
				out_row[j] += weight * probs_row[j];
			}
		}
	}
	return out;
}

double Score(const ModelShape &shape, const std::vector<double> &theta, const SparseCounts &x, const SparseCounts &y) {
	const std::vector<double> pred = Predict(shape, theta, x);
	CheckParallel(y, "Y");
	if (y.n_rows != x.n_rows) {
		throw std::invalid_argument("mmvec: X and Y must have the same number of samples (got " +
		                            std::to_string(x.n_rows) + " and " + std::to_string(y.n_rows) + ")");
	}
	if (y.n_cols != shape.n_features_y) {
		throw std::invalid_argument("mmvec: Y has " + std::to_string(y.n_cols) +
		                            " features but the model was fitted on " + std::to_string(shape.n_features_y));
	}
	const auto gy = GroupBySample(y, "Y", ColumnPolicy::kAllowAllZero);

	const int64_t n = y.n_rows;
	const int64_t d2 = shape.n_features_y;

	// Y as proportions, densely: both sums below run over every cell, structural
	// zeros included, so there is nothing for the sparse form to skip.
	std::vector<double> y_props(static_cast<size_t>(n * d2), 0.0);
	std::vector<size_t> order;
	for (int64_t s = 0; s < n; ++s) {
		const double total = OrderRowByFeature(gy, s, order);
		for (const size_t k : order) {
			y_props[static_cast<size_t>(s * d2) + gy.col[k]] = gy.val[k] / total;
		}
	}

	// The baseline Q^2 measures against is the mean abundance of each Y feature
	// across samples -- "the average community" -- not one global constant.
	std::vector<double> col_mean(static_cast<size_t>(d2), 0.0);
	for (int64_t s = 0; s < n; ++s) {
		for (int64_t j = 0; j < d2; ++j) {
			col_mean[static_cast<size_t>(j)] += y_props[static_cast<size_t>(s * d2 + j)];
		}
	}
	for (int64_t j = 0; j < d2; ++j) {
		col_mean[static_cast<size_t>(j)] /= static_cast<double>(n);
	}

	double ss_res = 0.0;
	double ss_tot = 0.0;
	for (int64_t s = 0; s < n; ++s) {
		for (int64_t j = 0; j < d2; ++j) {
			const size_t idx = static_cast<size_t>(s * d2 + j);
			const double res = y_props[idx] - pred[idx];
			const double tot = y_props[idx] - col_mean[static_cast<size_t>(j)];
			ss_res += res * res;
			ss_tot += tot * tot;
		}
	}
	if (ss_tot == 0.0) {
		// Every sample shares one profile, so there is no variance to explain.
		// scikit-bio returns 0 rather than dividing by zero, and 0 is the honest
		// answer: the mean baseline is already exact.
		return 0.0;
	}
	return 1.0 - ss_res / ss_tot;
}

} // namespace miint::mmvec
