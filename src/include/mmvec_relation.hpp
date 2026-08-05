#pragma once
//
// Translation in BOTH directions between the SQL long-form relations and the
// arrays the MMvec core works in: on the way in, id dictionaries, the sample join
// and the ORDERING RULE; on the way out, the fitted parameter vector unpacked
// into (modality, feature, axis, value) rows. One responsibility, stated once:
// relation shape <-> core shape.
//
// Split out of mmvec_fit_function.cpp, and DuckDB-free, because both halves are
// the wrapper's logic rather than the model's AND are the two places where a
// mistake is silent rather than loud:
//
//   * Which Y feature ends up at index 0 decides the reference category. The
//     Gaussian priors are not softmax-shift-invariant, so a different ordering
//     reaches a genuinely different optimum -- it does not fail, it just answers
//     differently.
//   * theta packs four blocks of different shapes back to back. Reading y_main
//     (p x d2-1, row-major) with the indices transposed, or getting a block
//     boundary off by one, still yields in-range values of the right count -- so
//     row counts, NULL patterns and the pinned-zero reference all still look
//     right while every emitted coordinate is mislabelled.
//
// Both therefore need assertions against an independently known quantity, which
// means they cannot live behind libduckdb (mirrors the
// community_distances.{hpp,cpp} split). BuildModelRows is checked by
// reconstructing the logits from its output and comparing against ComputeLogits
// at the oracle's carved theta.
//
// Consumes exactly what unifrac_internal::ReadFeatureTable returns, with no
// conversion: CooRow lives in feature_table_row.hpp precisely so generic
// consumers need not pull in the UniFrac feature.

#include <string>
#include <vector>

#include "feature_table_row.hpp"
#include "mmvec.hpp"

namespace miint::mmvec {

//! The two count tables with indices assigned, plus the dictionaries that map
//! those indices back to the user's ids for output.
struct IngestedTables {
	std::vector<std::string> sample_ids;    //!< lexicographically sorted; shared by both tables
	std::vector<std::string> x_feature_ids; //!< lexicographically sorted; index i is X feature i
	std::vector<std::string> y_feature_ids; //!< lexicographically sorted; index 0 is the REFERENCE
	SparseCounts x;
	SparseCounts y;
};

//! Assign sample and feature indices to a pair of long-form count tables.
//!
//! Feature ids are sorted LEXICOGRAPHICALLY, as byte strings. This is the
//! ordering rule, and it is load-bearing rather than cosmetic: `y_feature_ids[0]`
//! becomes MMvec's reference category, whose logit is pinned to zero, and since
//! the priors are not shift-invariant that choice moves the fitted ranks (up to
//! 31% of their magnitude on the toy fixture) and reaches a different optimum.
//! Sorting as bytes, not numerically, means "f10" sorts before "f9" -- matching
//! how the committed fixtures were generated, so the carved oracle applies.
//!
//! Samples are matched BY ID, never by position: the two relations are read
//! independently and a parallel scan fixes neither one's row order. Sample ids
//! are sorted too, which the model does not require -- the sample axis is summed
//! away into the sufficient statistics -- but which makes the assignment
//! reproducible for anyone inspecting it.
//!
//! Rows are consumed as given otherwise. In particular duplicate (sample,
//! feature) cells are REJECTED rather than summed: `X^T Y` would receive
//! `x_a*y + x_b*y` where the user meant `(x_a+x_b)*y`, so silently summing would
//! quietly change the fit. The core rejects them too, but only knows indices, so
//! the check happens here where the offending ids can be named.
//!
//! Throws std::invalid_argument if:
//!  - either table is empty;
//!  - the two tables' sample sets differ (naming the count and some examples);
//!  - a (sample, feature) cell is duplicated within either table;
//!  - fewer than 1 X feature, or fewer than 2 Y features (with only the
//!    reference category present every logit is pinned to zero and there is no
//!    likelihood to fit).
//!
//! Values are NOT validated here -- negatives, non-finite counts and all-zero
//! rows or columns are the core's contract, and ComputeSufficientStats reports
//! them.
IngestedTables IngestPairedTables(const std::vector<miint::unifrac::CooRow> &x_rows,
                                  const std::vector<miint::unifrac::CooRow> &y_rows);

//! One row of the model relation, before the ids and the broadcast diagnostics
//! are attached. `id_index` indexes `x_feature_ids` or `y_feature_ids` according
//! to `kind`, and is unused for `Loss`.
struct ModelRow {
	enum class Kind : uint8_t { X, Y, Loss };
	Kind kind = Kind::X;
	int64_t id_index = 0;
	int32_t axis = 0;
	double value = 0.0;
};

//! Unpack a fitted model into relation rows: every X feature, then every Y
//! feature, then the loss curve.
//!
//! `axis` runs 1..p for the embedding coordinates and then 0 for the bias, which
//! is exactly the column order of scikit-bio's
//! `x_embeddings = hstack([x_main, x_bias])`. For `Loss` rows it is the 1-based
//! evaluation ordinal -- evaluations, not optimizer iterations, matching
//! scikit-bio's `loss_curve_`, whose length is a property of the line search.
//!
//! Y feature 0 is emitted with every coordinate AND its bias at exactly zero. It
//! holds no parameters in `theta` at all -- `y_main` is `p x (d2-1)` -- but
//! scikit-bio publishes `y_embeddings` as `(d2, p+1)` with row 0 all zeros. Doing
//! the same keeps the relation rectangular and stops a user's lexicographically
//! first Y feature vanishing from its own model.
//!
//! Emits `(d1 + d2) * (p + 1) + loss_curve.size()` rows. Throws
//! std::invalid_argument if the model's shape and theta disagree.
std::vector<ModelRow> BuildModelRows(const Model &model);

} // namespace miint::mmvec
