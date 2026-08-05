#pragma once
//
// Pure (DuckDB-free) MMvec: joint embeddings of two paired count modalities via
// a conditional multinomial likelihood (Morton et al. 2019, Nat Methods
// 16:1306-1314). Given per-sample counts of X-features (canonically microbes)
// and Y-features (canonically metabolites), MMvec fits
//
//     P(Y_j | X_i) = softmax_j( x_main[i] . y_main[:,j] + x_bias[i] + y_bias[j] )
//
// by maximum a posteriori under Gaussian priors, and reports the conditional
// log-probabilities as microbe-metabolite co-occurrence scores.
//
// Header/impl split with NO DuckDB dependency so the Catch2 unit-test target
// links it without libduckdb (mirrors community_distances.{hpp,cpp} and
// simulate_resemblance.{hpp,cpp}); the DuckDB table-function wrapper reads the
// two feature tables, assigns the feature indices, calls in here, and emits the
// results as relations.
//
// Reimplemented from scikit-bio 0.7.3's `skbio.stats.ordination.mmvec` to
// reproduce its results; no scikit-bio code is vendored or executed. Expected
// values live in test/cpp/mmvec_oracle.hpp. See THIRD_PARTY_LICENSES.md.

#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace miint::mmvec {

//! One paired count table in COO (long) form: `row[k]`, `col[k]`, `val[k]`
//! describe a single (sample, feature) cell, with the three arrays parallel.
//! Cells absent from the arrays are zero -- this is the sparse feature-table
//! contract the rest of the codebase uses.
//!
//! Indices are assigned BY THE CALLER, and the assignment is load-bearing:
//! Y-feature 0 is MMvec's reference category (its logit is pinned to zero), and
//! because the Gaussian priors are not softmax-shift-invariant the fitted
//! optimum genuinely depends on which feature that is. The core consumes the
//! indices as given so the choice of ordering rule lives in exactly one place
//! (the SQL wrapper, which sorts feature ids lexicographically). See
//! data/mmvec/README.md.
struct SparseCounts {
	int64_t n_rows = 0; //!< samples
	int64_t n_cols = 0; //!< features
	std::vector<int64_t> row;
	std::vector<int64_t> col;
	std::vector<double> val;
};

//! The data's ENTIRE contribution to the MMvec objective, with the sample axis
//! summed away:
//!
//!     y_sums[i,j] = SUM_n X[n,i] * Y[n,j]          (d1 x d2, row-major)
//!     n_sums[i]   = SUM_j y_sums[i,j]              (d1)
//!
//! Every later iteration reads only these, so per-iteration cost is O(d1*d2)
//! regardless of how many samples were aggregated.
struct SufficientStats {
	int64_t n_features_x = 0; //!< d1
	int64_t n_features_y = 0; //!< d2
	std::vector<double> y_sums;
	std::vector<double> n_sums;
};

//! Sum the sample axis out of a pair of count tables: `y_sums = X^T Y` and its
//! row sums, accumulated as a per-sample sparse outer product so neither table
//! is ever densified.
//!
//! `n_sums` is deliberately derived as `rowsum(y_sums)` rather than
//! independently as `SUM_n X[n,i] * SUM_j Y[n,j]`. The two are equal in exact
//! arithmetic (and were verified bit-identical against scikit-bio on every
//! committed fixture), but deriving one from the other makes the loss term
//! `dot(n_sums, log_norm)` and the gradient term `y_sums - n_sums * probs`
//! exactly consistent by construction instead of consistent to within a
//! rounding error.
//!
//! The result is BIT-IDENTICAL for any ordering of the COO arrays. Each output
//! cell (i,j) receives exactly one contribution per sample -- `X[n,i]*Y[n,j]`
//! -- so the summation order is fixed by ascending sample index alone, never by
//! the input row order. That is what lets the SQL layer stream rows in whatever
//! order a parallel scan produces. Duplicate (row, col) entries are rejected
//! precisely because they would break this: a count split across two entries
//! would contribute `x_a*y + x_b*y` in place of `(x_a+x_b)*y`.
//!
//! Throws std::invalid_argument on:
//!  - either table's row/col/val arrays not being the same length;
//!  - `n_rows` differing between the tables, or being < 1;
//!  - `x.n_cols` < 1, or `y.n_cols` < 2 (with a lone Y-feature there is nothing
//!    but the reference category, so every logit is pinned to zero and the model
//!    has no content to fit);
//!  - any index outside [0, n_rows) x [0, n_cols);
//!  - a duplicate (row, col) cell;
//!  - a value that is negative or not finite (MMvec is a multinomial model on
//!    counts; a transformed table with negative entries is a modeling error, not
//!    an input to silently fit);
//!  - an all-zero row or column in either table. scikit-bio rejects these too,
//!    and for a reason worth restating: an all-zero feature column contributes
//!    nothing to the likelihood, so its parameters are driven entirely by the
//!    prior and the fit succeeds while that feature's output is meaningless.
SufficientStats ComputeSufficientStats(const SparseCounts &x, const SparseCounts &y);

//! Model dimensions. `n_components` is the latent dimension p.
struct ModelShape {
	int64_t n_features_x = 0; //!< d1
	int64_t n_features_y = 0; //!< d2
	int32_t n_components = 0; //!< p
};

//! Length of the packed parameter vector: `d1*p + d1 + p*(d2-1) + (d2-1)`.
//!
//! Parameters are packed as `[x_main (d1 x p), x_bias (d1), y_main (p x d2-1),
//! y_bias (d2-1)]`, every block C-order row-major -- the same layout scikit-bio
//! uses, so a theta vector is interchangeable between the two implementations.
//! The Y-side blocks are one column short because Y-feature 0 is the reference
//! category and holds no parameters at all.
int64_t NumParams(const ModelShape &shape);

//! Gaussian prior on every parameter, contributing
//! `0.5 * ||block - mean||^2 / scale^2` to the objective. The X-side prior
//! applies to `x_main` and `x_bias`, the Y-side to `y_main` and `y_bias`.
//!
//! These are NOT softmax-shift-invariant, which is exactly why the reference
//! category matters (see SparseCounts). Scales must be strictly positive.
struct Priors {
	double x_prior_mean = 0.0;
	double x_prior_scale = 1.0;
	double y_prior_mean = 0.0;
	double y_prior_scale = 1.0;
};

//! Preallocated scratch for one objective evaluation. A fit calls the objective
//! on the order of a thousand times and the intermediates are (rows x d2-1) --
//! 10 MB at cystic-fibrosis scale -- so they are hoisted out of the call.
//! Reused across evaluations; every field is fully overwritten on each call.
//!
//! A "row" is one unit of the likelihood: an X feature for the full-batch
//! objective, one sampled (sample, X-feature) draw for the minibatch.
struct Workspace {
	std::vector<double> logits;       //!< rows x (d2-1), non-reference only
	std::vector<double> resids;       //!< rows x (d2-1); shifted exps, then probabilities, then residuals, in place
	std::vector<double> row_totals;   //!< rows; the Y total each row is conditioned on
	std::vector<double> x_main_rows;  //!< rows x p, the rows' gathered x_main
	std::vector<double> x_bias_rows;  //!< rows
	std::vector<double> dx_main_rows; //!< rows x p, before scattering back to d1 x p
	std::vector<double> dx_bias_rows; //!< rows

	//! Size the buffers for `n_rows` rows under `shape`. Idempotent.
	void Resize(const ModelShape &shape, int64_t n_rows);
};

//! The full (d1 x d2) logit matrix, row-major, including the pinned zero column
//! for Y-feature 0:
//!
//!     logits[i,0]   = 0
//!     logits[i,j+1] = x_main[i] . y_main[:,j] + x_bias[i] + y_bias[j]
//!
//! The objective never materializes this -- it works on the (d1 x d2-1)
//! non-reference block -- but `ranks` and `probs` are defined over the full
//! matrix, so it is built explicitly for them.
//!
//! Throws std::invalid_argument if `theta.size() != NumParams(shape)` or the
//! shape is degenerate.
std::vector<double> ComputeLogits(const ModelShape &shape, const std::vector<double> &theta);

//! Negative log-posterior and its gradient at `theta`, over the whole data set.
//! This is the objective L-BFGS minimizes.
//!
//!     loss = SUM_i n_sums[i]*log_norm[i] - SUM_ij y_sums[i,j+1]*logits[i,j]
//!            + 0.5*||x_main - x_mean||^2 / x_scale^2 + (same for x_bias,
//!              y_main, y_bias)
//!
//! `grad` is resized to `NumParams(shape)` and written in the same packed layout
//! as `theta`. Returns the loss.
//!
//! Cost is O(d1*d2*p), independent of the sample count -- the samples were
//! summed away into `stats`.
//!
//! Throws std::invalid_argument if `theta` is the wrong length, the shape
//! disagrees with `stats`, `n_components < 1`, or a prior scale is <= 0.
double FullBatchLossAndGradient(const ModelShape &shape, const Priors &priors, const SufficientStats &stats,
                                const std::vector<double> &theta, Workspace &ws, std::vector<double> &grad);

//! Which likelihood scale the minibatch objective uses. The two are NOT two
//! step sizes for one posterior -- see BatchNormFactor.
enum class BatchNorm {
	//! `sum(X) / batch_size`: makes the minibatch loss an unbiased estimator of
	//! the full-batch objective, so Adam and L-BFGS target the same posterior.
	Unbiased,
	//! `n_samples / batch_size`: upstream mmvec's scaling. Kept for comparison
	//! against published mmvec output.
	Legacy,
};

//! The likelihood scale for a minibatch of `batch_size` draws.
//!
//! `Legacy` multiplies the likelihood by `n_samples / sum(X)` relative to
//! `Unbiased` while leaving the priors untouched, so it is a materially
//! different posterior -- typically far more heavily regularized -- rather than
//! a rescaled version of the same one. That is why `Unbiased` is the default and
//! `Legacy` exists only for reproducing upstream numbers.
//!
//! Throws std::invalid_argument if `batch_size < 1`, `n_samples < 1`, or
//! `x_total <= 0`.
double BatchNormFactor(BatchNorm mode, double x_total, int64_t n_samples, int64_t batch_size);

//! What the minibatch objective needs, which is deliberately NOT the sufficient
//! statistics. Adam's sampled unit is a single nonzero (sample, X-feature) cell
//! of X, and the likelihood for that draw is the sample's WHOLE Y row, so the
//! sample axis cannot be summed away first.
//!
//! `x_row` / `x_col` / `x_val` are X's COO triples, one entry per nonzero; a
//! batch is a list of indices into them. All three are always required, even
//! though the objective itself never reads `x_val` -- an X count enters the
//! minibatch likelihood only by making that cell more likely to be drawn. The
//! Adam driver does need it, for the sampling weights and the normalization
//! factor, and a field that is required only sometimes is worse than one the
//! caller always fills in from the COO it already has.
//!
//! `y_dense` is the (n_samples x d2) Y table, row-major and including the
//! reference column. It is dense because `delta` spans every non-reference
//! category regardless of how sparse the row is.
struct MinibatchInputs {
	int64_t n_samples = 0;
	std::vector<int64_t> x_row;
	std::vector<int64_t> x_col;
	std::vector<double> x_val;
	std::vector<double> y_dense;
};

//! Negative log-posterior and its gradient over one minibatch -- the objective
//! Adam descends.
//!
//!     loss = norm * (SUM_b total_b*log_norm_b - SUM_bj Y[n_b,j+1]*logits[b,j])
//!            + the same prior terms as the full-batch objective
//!
//! `norm` comes from BatchNormFactor. Priors enter at FULL strength, unscaled by
//! `norm`: they are not part of the sampled sum, so scaling them would change
//! the posterior rather than estimate it.
//!
//! An X-feature appearing several times in `batch` accumulates several times, in
//! batch order -- matching numpy's `add.at`, which the reference implementation
//! relies on and whose summation order the committed Adam oracle encodes.
//!
//! `grad` is resized to `NumParams(shape)` and uses the packed layout. Returns
//! the loss.
//!
//! Throws std::invalid_argument on a wrong-length `theta`, a shape inconsistent
//! with `inputs`, a non-positive prior scale, an out-of-range batch index, or an
//! empty batch.
double MinibatchLossAndGradient(const ModelShape &shape, const Priors &priors, const MinibatchInputs &inputs,
                                const std::vector<int64_t> &batch, double norm, const std::vector<double> &theta,
                                Workspace &ws, std::vector<double> &grad);

//! scikit-bio's stopping rule, hard-coded there and matched here on purpose.
//!
//! Under this rule the converged parameters remain starting-point dependent at
//! roughly 1e-4 -- the optimum itself is unique to about 1e-7, so that spread is
//! an artifact of stopping early, not of multiple optima. Converging tighter
//! would be easy and would break parity, so it is not done.
inline constexpr double kLbfgsFtol = 1e-9; //!< relative objective-change tolerance
inline constexpr double kLbfgsGtol = 1e-5; //!< tolerance on max|gradient|
inline constexpr int kLbfgsHistory = 10;   //!< L-BFGS correction pairs (SciPy's maxcor)

//! A fitted model: the packed parameters plus an honest account of how the fit
//! went.
//!
//! Only the packed `theta` is stored, not four unpacked copies of it, so there is
//! one source of truth for the parameters; `ComputeLogits` and the derived
//! quantities read it directly.
struct Model {
	ModelShape shape;
	std::vector<double> theta; //!< NumParams(shape) values, packed layout

	//! One entry per objective EVALUATION, line-search probes included -- not per
	//! optimizer iteration. This matches scikit-bio, whose `loss_curve_` is
	//! likewise an evaluation trace despite being documented as per-iteration, so
	//! the length is a property of the line search and must never be asserted for
	//! equality across implementations.
	std::vector<double> loss_curve;

	//! The FULL-BATCH objective at `theta`, always -- so it is comparable between
	//! optimizers and against the L-BFGS tier. For L-BFGS that is also the minimum
	//! over all evaluations. For Adam it is deliberately NOT `loss_curve.back()`:
	//! minibatch losses are noisy estimates over different draws and are not
	//! comparable to each other, let alone across optimizers.
	double final_loss = 0.0;
	int64_t n_evals = 0; //!< == loss_curve.size()
	int64_t n_iter = 0;  //!< L-BFGS iterations, or Adam parameter updates

	//! max|gradient| at `theta`. Reported separately from `converged` because the
	//! two answer different questions, and on this stopping rule they often
	//! disagree: the objective-change test can fire while the gradient is still
	//! far from zero (measured on the committed fixtures: max|g| of 1e-3 to 1e-1
	//! at a legitimately-converged stop). Anyone who needs a stationary point
	//! rather than a stalled one should read this, not `converged`.
	double max_abs_grad = 0.0;
	//! True when the optimizer stopped on one of its own convergence tests, rather
	//! than exhausting `max_iter` or failing its line search -- the same meaning
	//! as SciPy's `success`, which is what scikit-bio surfaces.
	//!
	//! Adam therefore always reports false: it has no convergence test at all and
	//! simply runs its fixed epoch schedule to the end. That is not a failure, and
	//! `message` says so.
	//!
	//! It does NOT mean max|gradient| is small; see `max_abs_grad`. LBFGS++
	//! returns the same iteration count whether a convergence test fired or the
	//! cap was hit, so the two are told apart by `n_iter < max_iter`. That is
	//! conservative on the exact boundary -- a convergence test satisfied on the
	//! final permitted iteration is reported as hitting the cap -- so this flag
	//! can under-report convergence but never over-report it.
	//!
	//! Real data routinely fails to converge. scikit-bio neither raises nor warns
	//! in that case, which is why this is surfaced explicitly.
	bool converged = false;

	//! Human-readable outcome: which stopping condition fired, the iteration and
	//! evaluation counts, and max|gradient|.
	std::string message;
};

//! Minimize the full-batch objective with L-BFGS, starting from `theta0`.
//!
//! The starting point is an explicit argument rather than drawn internally, both
//! because the committed oracle pins its own init and because it keeps this
//! function free of any RNG.
//!
//! `theta` is the best point seen across ALL objective evaluations rather than
//! the solver's last iterate. Measured on every committed fixture, including the
//! runs that hit the iteration cap, those are the SAME point -- so this is not a
//! divergence from scikit-bio in practice, and it is never worse than the last
//! iterate. It earns its place in one situation: LBFGS++ signals line-search
//! failure by throwing and leaves its iterate at the starting point, so without
//! the snapshot such a fit would raise instead of returning the best parameters
//! it had already found.
//!
//! Throws std::invalid_argument if `max_iter < 1`, or for any reason
//! FullBatchLossAndGradient would -- including a non-finite `theta0`, which would
//! otherwise poison every loss and gradient without any comparison failing.
Model FitLbfgsFromInit(const ModelShape &shape, const Priors &priors, const SufficientStats &stats,
                       const std::vector<double> &theta0, int64_t max_iter);

//! Deterministic pseudo-random draws that are byte-identical on every platform.
//!
//! This is NOT a hand-rolled generator. `std::mt19937_64` supplies every bit, and
//! the standard fixes its output exactly -- period, equidistribution and seeding
//! are all upstream and untouched (the engine pin in the tests checks the
//! standard's own published 10000th-output vector). Only the transforms from bits
//! to doubles are written out here instead of taken from `<random>`, because that
//! is the sole part the standard leaves open. `std::uniform_real_distribution`,
//! `std::normal_distribution` and `std::discrete_distribution` are all
//! implementation-defined, so the same seed legitimately produces different
//! numbers on libstdc++, libc++ and Emscripten.
//!
//! `cluster_kmeans` documents that hazard and accepts it -- its seeded draws may
//! differ between standard libraries. mmvec deliberately does not, because a
//! seeded fit is expected to reproduce everywhere including WASM, and because
//! accepting it would leave the Adam path comparable only statistically. The cost
//! is about thirty lines.
class Rng {
public:
	explicit Rng(uint64_t seed);

	//! Uniform on [0, 1), from the top 53 bits -- the canonical construction.
	double Uniform();

	//! Standard normal by the Marsaglia polar method, which needs only sqrt and
	//! log (no erfinv, no lookup tables) and so is exactly reproducible. Draws
	//! come in pairs; the spare is cached and returned by the next call.
	double Normal();

private:
	std::mt19937_64 engine_;
	double cached_normal_ = 0.0;
	bool has_cached_normal_ = false;
};

//! `count` indices into `weights`, drawn with replacement with probability
//! proportional to the weights. Implemented as a binary search over the
//! cumulative sums, which is deterministic given the uniform stream.
//!
//! Throws std::invalid_argument if `weights` is empty, `count < 0`, any weight is
//! negative or non-finite, or the weights sum to zero.
std::vector<int64_t> SampleWeightedIndices(const std::vector<double> &weights, int64_t count, Rng &rng);

//! Standard-normal starting parameters, drawn in scikit-bio's order: all of
//! `x_main`, then `x_bias`, then `y_main`, then `y_bias` -- i.e. simply in packed
//! order, which is why one sequential pass suffices.
//!
//! The VALUES cannot match scikit-bio, which draws from numpy's PCG64 stream;
//! only the shape, the order and the distribution do. That is why every parity
//! tier in the oracle injects its starting point explicitly instead of relying on
//! a seed.
std::vector<double> InitTheta(const ModelShape &shape, uint64_t seed);

//! Adam hyperparameters, defaulted to scikit-bio's own values. Note `beta_2` is
//! 0.95, not the 0.999 usual elsewhere.
struct AdamParams {
	double learning_rate = 1e-3;
	int64_t batch_size = 50;
	double beta_1 = 0.9;
	double beta_2 = 0.95;
	double clipnorm = 10.0;
	BatchNorm batch_norm = BatchNorm::Unbiased;
};

//! Adam over the minibatch objective, with the draw sequence supplied by the
//! caller: `batches` is the per-update index lists concatenated, so its length
//! must be a multiple of `params.batch_size` and it defines the update count.
//!
//! This exists because it removes the RNG from any comparison, which is the only
//! way Adam can be checked exactly against the reference implementation. It is
//! the entry point the committed Adam oracle drives.
//!
//! Per update, in this order: evaluate the minibatch loss and gradient (priors
//! included), clip the gradient to `clipnorm` by its GLOBAL L2 norm across all
//! four blocks together, then apply Adam with bias corrections
//! `1 - beta^update`. The loss recorded is the one measured BEFORE that update,
//! matching scikit-bio.
//!
//! Throws std::invalid_argument if `batches` is empty or not a whole number of
//! batches, if a hyperparameter is out of range, or for any reason
//! MinibatchLossAndGradient would.
Model FitAdamWithIndices(const ModelShape &shape, const Priors &priors, const MinibatchInputs &inputs,
                         const AdamParams &params, const std::vector<double> &theta0,
                         const std::vector<int64_t> &batches);

//! Adam drawing its own minibatches from `seed`, over `epochs` passes.
//!
//! Each epoch performs `max(1, nnz / batch_size)` updates, where nnz counts X's
//! nonzero cells -- so `epochs` is not the number of parameter updates.
//! Cells are drawn with replacement, weighted by their X count.
//!
//! Throws std::invalid_argument if `epochs < 1`, or for any reason
//! FitAdamWithIndices would.
Model FitAdam(const ModelShape &shape, const Priors &priors, const MinibatchInputs &inputs, const AdamParams &params,
              const std::vector<double> &theta0, int64_t epochs, uint64_t seed);

//! The four quantities a fitted model is actually read through.
//!
//! All of them take `(shape, theta)` rather than a `Model`, matching
//! `ComputeLogits`, so the parity tiers can evaluate them at the oracle's carved
//! `theta` without fabricating a fit. A caller holding a `Model` passes
//! `model.shape, model.theta`.
//!
//! The embeddings themselves are NOT here, and deliberately: they are identified
//! only up to an orthogonal rotation of the latent axes, so two runs that agree
//! perfectly on every quantity below can still have unrelated `x_main`/`y_main`.
//! These four are the rotation-invariant ones, which is why they are what parity
//! is measured on and what the SQL layer will expose.

//! Row-centred logits, (d1 x d2) row-major: `logits - rowmean(logits)`.
//!
//! This is MMvec's headline output -- the log-conditional-probability of each Y
//! feature given each X feature, relative to that X feature's average. Centring
//! makes rows comparable to each other and removes the arbitrariness of which Y
//! feature was pinned as the reference; rows sum to zero by construction.
//!
//! Throws std::invalid_argument if `theta` is the wrong length or non-finite, or
//! the shape is degenerate.
std::vector<double> Ranks(const ModelShape &shape, const std::vector<double> &theta);

//! Conditional probabilities P(Y_j | X_i), (d1 x d2) row-major, rows summing to 1.
//!
//! Computed as `softmax(logits)`, which is algebraically `softmax(ranks)` -- the
//! softmax is shift-invariant, so centring cannot change it. Both forms subtract
//! the row maximum before exponentiating, so both are numerically stable and they
//! agree to rounding; the tests assert that agreement across the two paths rather
//! than assuming it.
//!
//! Throws std::invalid_argument for the same reasons as Ranks.
std::vector<double> Probs(const ModelShape &shape, const std::vector<double> &theta);

//! Expected Y proportions for each sample of `x`: `rowprops(x) @ Probs()`,
//! (x.n_rows x d2) row-major, rows summing to 1.
//!
//! Reading a sample's X counts as proportions is what makes the prediction
//! invariant to sequencing depth -- scaling a row of `x` cannot change its
//! prediction, and there is a test for exactly that.
//!
//! Accumulated sparsely (O(nnz*d2), not the O(n*d1*d2) of a dense product) and in
//! ascending feature order within each sample, so the result is BIT-IDENTICAL for
//! any ordering of the COO arrays. The sparse form is not a micro-optimization:
//! the cystic-fibrosis fixture's X is 9810 cells of a possible 172*2720, i.e. 2%
//! dense, which is 4.5M multiply-adds here against 216M densely. Unlike `ComputeSufficientStats` that is not automatic
//! here -- an output cell receives one contribution per X nonzero in its row, so without the ordering the
//! floating-point summation order would follow whatever order a parallel scan produced.
//!
//! Throws std::invalid_argument if `theta` is bad, `x`'s parallel arrays
//! disagree, `x.n_rows < 1`, `x.n_cols != shape.n_features_x`, an index is out of
//! range, a value is negative or non-finite, a cell is duplicated, or a sample's
//! counts are all zero (there are no proportions to form). An all-zero FEATURE is
//! accepted here, unlike when fitting: a feature absent from these samples merely
//! contributes no weight, whereas when fitting its parameters would be determined
//! entirely by the prior.
//!
//! DIVERGENCE (for the M6 ledger): scikit-bio's `predict` validates neither X's
//! width nor its feature identities, and consumes X positionally. We check the
//! width, which is cheap and catches the realistic error of passing a
//! differently-filtered X table. Aligning feature IDs against the fitted model is
//! the SQL layer's job, where the IDs actually exist.
std::vector<double> Predict(const ModelShape &shape, const std::vector<double> &theta, const SparseCounts &x);

//! Q^2, the fraction of variance in Y's proportions explained by the prediction:
//!
//!     SS_res = SUM_nj (y_props[n,j] - predict(x)[n,j])^2
//!     SS_tot = SUM_nj (y_props[n,j] - colmean(y_props)[j])^2
//!     Q^2    = 1 - SS_res/SS_tot
//!
//! `colmean` is per-Y-feature over samples, so the baseline being beaten is "the
//! mean community", not a global constant. Unlike R^2 this is evaluated on data
//! the model did not fit, so it is unbounded below -- a badly fitting model
//! scores arbitrarily negative, and the toy fixture's own T1 value is about -2.0.
//! Returns exactly 0.0 when `SS_tot == 0`, matching scikit-bio rather than
//! dividing by zero. That test is exact equality, which in practice means it
//! fires only for a SINGLE sample -- where the mean community is that sample, so
//! the answer is unambiguous. Several identical samples do NOT reach it:
//! `colmean` sums n copies and divides by n, which is an ulp off for most values,
//! leaving `SS_tot` around 1e-33 and the score at some enormous negative number.
//! An almost-degenerate Y is therefore not protected. Widening the guard to a
//! tolerance would change results relative to scikit-bio, so it is left alone
//! here and raised in the M6 divergence ledger instead.
//!
//! Throws std::invalid_argument for any reason Predict would, or if `y`'s shape
//! disagrees with `x.n_rows` or `shape.n_features_y`.
double Score(const ModelShape &shape, const std::vector<double> &theta, const SparseCounts &x, const SparseCounts &y);

} // namespace miint::mmvec

namespace duckdb {
class ExtensionLoader;
//! Registers the mmvec_fit table function into the extension catalog. Declared
//! here behind a forward declaration, so this header stays free of DuckDB
//! includes and the Catch2 target still links the core without libduckdb.
void RegisterMmvecFit(ExtensionLoader &loader);
} // namespace duckdb
