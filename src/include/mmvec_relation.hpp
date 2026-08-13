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

#include <optional>
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

//! The `modality` string a ModelRow::Kind is written as, and read back as.
//!
//! The one place the enum and the wire vocabulary meet. Without it the mapping
//! gets rewritten at every site that crosses the boundary -- the SQL emitter, the
//! parser's validation, and the round-trip test, which all had their own copy --
//! and a renamed modality would drift the encode and decode paths apart while the
//! round-trip test, carrying its own third copy, stayed green. That is the same
//! silent-agreement hazard the block-offset comment at the top of
//! mmvec_relation.cpp warns about.
inline const char *ModalityName(ModelRow::Kind kind) {
	switch (kind) {
	case ModelRow::Kind::X:
		return "x";
	case ModelRow::Kind::Y:
		return "y";
	case ModelRow::Kind::Loss:
		return "loss";
	}
	return "loss"; // unreachable; enumerated above, and silences -Wreturn-type
}

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

//! One cell of a model relation as SQL presents it: an id, not an index.
struct ModelCell {
	//! Exactly "x", "y" or "loss", lowercase. The caller lowercases; anything else
	//! is rejected by name rather than guessed at.
	std::string modality;
	//! The id `modality` SELECTS -- `x_feature_id` for "x", `y_feature_id` for "y"
	//! -- and `nullopt` for "loss". Deliberately not "whichever id is non-NULL":
	//! that would read a Y id as an X feature name on a mislabelled row, inventing
	//! a feature instead of complaining. `std::optional` rather than an empty-string
	//! sentinel because "" is a legal VARCHAR feature id.
	std::optional<std::string> feature_id;
	int32_t axis = 0;
	double value = 0.0;
};

//! A model relation read back into the arrays the core works in.
struct ParsedModel {
	ModelShape shape;
	std::vector<double> theta;
	std::vector<std::string> x_feature_ids; //!< index i is X feature i
	std::vector<std::string> y_feature_ids; //!< index 0 is the REFERENCE, then the rest sorted
};

//! The inverse of BuildModelRows: rebuild `(shape, theta)` and both id
//! dictionaries from a model relation, so a fitted model can be read back and
//! evaluated. `modality = 'loss'` cells are ignored, so a model table can be
//! passed through unfiltered.
//!
//! ORDER IS NOT RECOVERED, AND DOES NOT NEED TO BE. Lexicographic feature order
//! is load-bearing when FITTING, because it decides which Y feature becomes the
//! reference category and the Gaussian priors are not shift-invariant. On the way
//! back in, the only structural requirement is that the reference sits at index 0:
//! `ComputeLogits` pins Y feature 0's logit to zero and derives the rest from
//! `y_main[:, j-1]`. Any order for the others works, because this function packs
//! `theta` and labels `y_feature_ids` in the SAME order, so permuting a
//! non-reference feature moves its column in both places and its value follows.
//!
//! The reference is therefore identified by the property that actually defines it
//! -- it is the Y feature holding no parameters, which BuildModelRows emits as
//! all zeros -- and NOT by lexicographic position. That matters in practice:
//! joining a model to metadata and renaming features to readable names is an
//! ordinary thing to do, and it must not silently reassign the reference category.
//! Exactly one all-zero Y feature is required; zero or several is ambiguous and
//! rejected rather than guessed at.
//!
//! Throws std::invalid_argument if:
//!  - a `modality` is not one of "x", "y", "loss";
//!  - an "x" or "y" cell has a NULL feature id;
//!  - either modality has no cells at all;
//!  - a value is not finite, or an `axis` is negative;
//!  - a (modality, feature, axis) cell is duplicated;
//!  - the grid is incomplete -- some feature is missing an axis (naming both);
//!  - the two modalities disagree on the number of axes;
//!  - `p < 1`, or `d2 < 2` (with only a reference category there is no model);
//!  - the number of all-zero Y features is not exactly one.
ParsedModel ParseModelCells(const std::vector<ModelCell> &cells);

//! A long-form table indexed against a FITTED MODEL's feature dictionary, so its
//! columns line up with that model's theta.
struct AlignedTable {
	std::vector<std::string> sample_ids; //!< lexicographically sorted
	SparseCounts counts;                 //!< `n_cols` is the MODEL's width, always
};

//! Index a held-out X table against the model's X features, ready for Predict.
//!
//! This is the counterpart to IngestPairedTables, and the difference is the whole
//! point: ingest builds a dictionary FROM the data, whereas here the dictionary is
//! the model's and the data is fitted into it. Building a fresh dictionary over a
//! held-out table instead would yield in-range column indices of plausibly the
//! right count that name entirely different features -- a silent wrong answer, not
//! a failure. `counts.n_cols` is therefore `x_feature_ids.size()` regardless of
//! which features actually appear.
//!
//! UNKNOWN vs MISSING, following scikit-bio:
//!  - A feature the model never saw is an ERROR, naming the count and examples.
//!    scikit-bio's `predict` takes a dense array with no ids, so an extra column is
//!    a shape mismatch there; and there is no `P(Y | unknown microbe)` to
//!    contribute. Restricting the table to the model's features is the caller's
//!    decision to make explicitly, in SQL, not ours to make silently.
//!  - A model feature ABSENT from the table is fine: in scikit-bio it would be a
//!    zero column, and in long form "no rows" is exactly that. `Predict` documents
//!    accepting an all-zero feature -- unlike when fitting, where its parameters
//!    would come from the prior alone -- and this is the case it accepts it for.
//!
//! Sample ids are collected from the table itself (the model does not constrain
//! them; predicting for new samples is the point) and sorted, so the output order
//! is reproducible.
//!
//! Throws std::invalid_argument if the table has no usable cells, if any feature is
//! unknown to the model, or if a (sample, feature) cell is duplicated. Values are
//! NOT validated here -- negatives, non-finite counts and all-zero samples are
//! `Predict`'s contract, and it reports them.
AlignedTable AlignXToModel(const std::vector<miint::unifrac::CooRow> &x_rows,
                           const std::vector<std::string> &x_feature_ids);

//! A held-out X/Y pair indexed against a fitted model, over ONE sample dictionary.
struct AlignedPair {
	std::vector<std::string> sample_ids; //!< lexicographically sorted; shared by both
	SparseCounts x;                      //!< `n_cols` is the model's X width
	SparseCounts y;                      //!< `n_cols` is the model's Y width
};

//! Index a held-out X and Y against the model's two feature dictionaries, ready
//! for Score.
//!
//! Every rule AlignXToModel states applies to both modalities here -- the model's
//! dictionary rather than a fresh one, unknown features rejected, absent features
//! left as zero columns -- plus the one thing scoring adds: X and Y must describe
//! the SAME samples, matched BY ID, and they are indexed over ONE sample dictionary
//! so that row n of `x` and row n of `y` are the same sample. Score compares the
//! two tables cell by cell, so two independently built sample orders would pair one
//! sample's microbes against another's metabolites and report a plausible number
//! for it. A sample set mismatch is rejected with the same message
//! IngestPairedTables gives, because it is the same mistake.
//!
//! Throws std::invalid_argument if either table has no usable cells, if the two
//! sample sets differ, if any feature is unknown to the model, or if a (sample,
//! feature) cell is duplicated. Values are Score's contract, and it reports them.
AlignedPair AlignPairedToModel(const std::vector<miint::unifrac::CooRow> &x_rows,
                               const std::vector<miint::unifrac::CooRow> &y_rows,
                               const std::vector<std::string> &x_feature_ids,
                               const std::vector<std::string> &y_feature_ids);

} // namespace miint::mmvec
