#pragma once
//
// mmvec's numeric kernel -- the two functions that dominate a fit, split into
// their own translation unit.
//
// This is an INTERNAL header. It is not part of mmvec's public surface (that is
// `mmvec.hpp`); it exists because Eigen chooses its packet width at COMPILE time
// from `__AVX__` / `__AVX512F__` and offers no runtime dispatch. Getting a wider
// instruction set into these two functions therefore means compiling them again
// with different flags, which means they cannot live in `mmvec.cpp` alongside
// code that must stay on the baseline instruction set.
//
// The seam is deliberately narrow: everything below is either a plain data
// layout or one of the two entry points. See docs/internals/embedded-tools.md.

#include "mmvec.hpp"

#include <Eigen/Core>

#include <cstdint>

namespace miint::mmvec {

//! Row-major, to match the packed parameter layout so Eigen can map the blocks
//! of `theta` and `grad` in place with no copying or transposing.
using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using RowMatrixMap = Eigen::Map<RowMatrix>;
using ConstRowMatrixMap = Eigen::Map<const RowMatrix>;

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

inline ParamLayout Layout(const ModelShape &shape) {
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

//! The two entry points, declared once and instantiated into one namespace per
//! instruction set. mmvec_kernel.cpp is compiled once per variant with
//! `-DMIINT_ISA_VARIANT=<name>`, so each build lands in its own namespace and the
//! three copies coexist in one binary without a symbol clash.
//!
//! A macro rather than three hand-written copies because the signatures must not
//! drift between tiers -- a variant declared with a stale signature would compile
//! and then dispatch to the wrong thing.
//!
//! FillNonRefLogits: the model equation, evaluated into `out` (row-major,
//! n_rows x d2-1). EvaluateObjective: one evaluation of the negative
//! log-posterior and its gradient, the hot function -- two thirds of an L-BFGS
//! fit and two fifths of an Adam one. Full contracts at the definitions.
#define MIINT_MMVEC_KERNEL_DECLARATIONS                                                                                \
	void FillNonRefLogits(const ParamLayout &l, const double *theta, const int64_t *x_index, int64_t n_rows,           \
	                      Workspace &ws, double *out);                                                                 \
	double EvaluateObjective(const ParamLayout &l, const Priors &priors, int64_t n_rows, const int64_t *x_index,       \
	                         const ObsView &obs, const double *totals, double norm, const double *theta,               \
	                         Workspace &ws, double *grad);

//! Always built, and the only variant guaranteed to exist.
//!
//! It is also where the READ path stays. `ComputeLogits` -- and therefore `Ranks`
//! / `Probs` / `Predict` / `Score` -- calls FillNonRefLogits once per fit rather
//! than once per iteration, so it has no performance stake in a wider register.
//! Pinning it here keeps every derived carved value bit-exact on every CPU and
//! confines instruction-set dependence to the fitted theta alone. Measured: at a
//! fixed theta the wide variants move `ranks` by 4e-15 and change no top-ranked
//! feature, so this costs nothing and buys an exact oracle.
namespace baseline {
MIINT_MMVEC_KERNEL_DECLARATIONS
}

#ifdef MIINT_HAS_AVX2
namespace avx2 {
MIINT_MMVEC_KERNEL_DECLARATIONS
}
#endif

#ifdef MIINT_HAS_AVX512
namespace avx512 {
MIINT_MMVEC_KERNEL_DECLARATIONS
}
#endif

} // namespace miint::mmvec
