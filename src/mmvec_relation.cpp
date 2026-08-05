#include "mmvec_relation.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

namespace miint::mmvec {

namespace {

//! Where each of theta's four blocks starts: x_main (d1 x p), x_bias (d1),
//! y_main (p x d2-1), y_bias (d2-1), every block C-order row-major.
//!
//! Spelled out ONCE, and shared by both directions. Note what that costs: an error
//! made consistently in BuildModelRows and ParseModelCells CANCELS, so the
//! round-trip test cannot see it -- verified by mutation, not assumed. Two tests in
//! test_MMvecRelation.cpp are the anchors that can, and neither is redundant:
//! "model rows reconstruct the oracle's logits" reads the emitted ROWS through an
//! independent reconstruction, and "parsing a relation written from theta recovers
//! theta" feeds ParseModelCells cells built from an independent spelling of this
//! layout. Don't drop either on the grounds that the round-trip covers it.
struct Offsets {
	size_t x_main;
	size_t x_bias;
	size_t y_main;
	size_t y_bias;
};

Offsets BlockOffsets(int64_t d1, int64_t d2, int64_t p) {
	Offsets o {};
	o.x_main = 0;
	o.x_bias = static_cast<size_t>(d1 * p);
	o.y_main = o.x_bias + static_cast<size_t>(d1);
	o.y_bias = o.y_main + static_cast<size_t>(p * (d2 - 1));
	return o;
}

//! Distinct ids from one axis of a long-form table, lexicographically sorted,
//! paired with the id -> index map that sorting implies.
struct Dictionary {
	std::vector<std::string> ids;
	std::unordered_map<std::string, int64_t> index;
};

//! `select` picks which of a row's two ids to collect.
template <typename Select>
Dictionary BuildDictionary(const std::vector<miint::unifrac::CooRow> &rows, Select select) {
	Dictionary d;
	d.ids.reserve(rows.size());
	for (const auto &r : rows) {
		d.ids.push_back(select(r));
	}
	// Sort-then-unique rather than a hash set, so the order is the sort order and
	// not an insertion order that happened to be sorted on the test fixtures.
	std::sort(d.ids.begin(), d.ids.end());
	d.ids.erase(std::unique(d.ids.begin(), d.ids.end()), d.ids.end());

	d.index.reserve(d.ids.size());
	for (size_t i = 0; i < d.ids.size(); ++i) {
		d.index.emplace(d.ids[i], static_cast<int64_t>(i));
	}
	return d;
}

//! Ids in `a` that are absent from `b`, at most `limit` of them, for error text.
std::vector<std::string> MissingFrom(const std::vector<std::string> &a,
                                     const std::unordered_map<std::string, int64_t> &b, size_t limit) {
	std::vector<std::string> out;
	for (const auto &id : a) {
		if (b.find(id) == b.end()) {
			out.push_back(id);
			if (out.size() == limit) {
				break;
			}
		}
	}
	return out;
}

std::string Join(const std::vector<std::string> &ids) {
	std::string out;
	for (size_t i = 0; i < ids.size(); ++i) {
		if (i != 0) {
			out += ", ";
		}
		out += ids[i];
	}
	return out;
}

//! A Dictionary over ids that are already fixed -- a fitted model's feature list.
//! Unlike BuildDictionary this does NOT sort: the order is the model's, and for Y
//! features that order is load-bearing (index 0 is the reference category).
Dictionary DictionaryFromIds(const std::vector<std::string> &ids) {
	Dictionary d;
	d.ids = ids;
	d.index.reserve(d.ids.size());
	// Indexed off d.ids, not the parameter, so the two fields cannot disagree --
	// the same way BuildDictionary does it. Reading the parameter instead would let
	// a future reorder of d.ids leave the mapping silently stale.
	for (size_t i = 0; i < d.ids.size(); ++i) {
		d.index.emplace(d.ids[i], static_cast<int64_t>(i));
	}
	return d;
}

//! The two tables must describe the same samples, matched BY ID. Shared by the
//! fitting and the scoring paths: it is the same mistake in both, and a user who
//! has read this message once should not have to learn a second wording of it.
void RequireSameSamples(const Dictionary &x_samples, const Dictionary &y_samples) {
	if (x_samples.ids == y_samples.ids) {
		return;
	}
	// Intersecting instead would silently answer a question about fewer samples
	// than the user asked about.
	const auto x_only = MissingFrom(x_samples.ids, y_samples.index, 3);
	const auto y_only = MissingFrom(y_samples.ids, x_samples.index, 3);
	std::string detail;
	if (!x_only.empty()) {
		detail += "; only in X: " + Join(x_only);
	}
	if (!y_only.empty()) {
		detail += "; only in Y: " + Join(y_only);
	}
	throw std::invalid_argument("mmvec: the X and Y feature-tables must describe the same samples (X has " +
	                            std::to_string(x_samples.ids.size()) + ", Y has " +
	                            std::to_string(y_samples.ids.size()) + ")" + detail);
}

//! Reject cells naming a feature the model never saw, before Index dereferences
//! them. Index uses .at(), so without this the failure would be an out_of_range
//! with no message a user could act on.
void RejectUnknownFeatures(const std::vector<miint::unifrac::CooRow> &rows, const Dictionary &features,
                           const std::string &label) {
	std::vector<std::string> examples;
	size_t n_unknown = 0;
	std::unordered_map<std::string, size_t> seen;
	for (const auto &r : rows) {
		if (features.index.find(r.feature_id) != features.index.end()) {
			continue;
		}
		if (seen.emplace(r.feature_id, 0).second) {
			++n_unknown;
			if (examples.size() < 3) {
				examples.push_back(r.feature_id);
			}
		}
	}
	if (n_unknown != 0) {
		throw std::invalid_argument(
		    "mmvec: the " + label + " table has " + std::to_string(n_unknown) +
		    " feature(s) the model was never fitted on (" + Join(examples) +
		    (n_unknown > examples.size() ? ", ..." : "") +
		    "); there is no conditional probability for a feature the model has not seen, so restrict the table to the "
		    "model's own features first -- e.g. WHERE feature_id IN (SELECT DISTINCT " +
		    (label == "X" ? std::string("x_feature_id") : std::string("y_feature_id")) +
		    " FROM <model> WHERE modality = '" + (label == "X" ? "x" : "y") + "') -- or refit including them");
	}
}

//! Index one table against the shared sample dictionary and its own feature
//! dictionary, rejecting duplicated cells by id.
SparseCounts Index(const std::vector<miint::unifrac::CooRow> &rows, const Dictionary &samples,
                   const Dictionary &features, const std::string &label) {
	SparseCounts t;
	t.n_rows = static_cast<int64_t>(samples.ids.size());
	t.n_cols = static_cast<int64_t>(features.ids.size());
	t.row.reserve(rows.size());
	t.col.reserve(rows.size());
	t.val.reserve(rows.size());

	// Duplicate detection over the linearized (sample, feature) cell. The core
	// rejects duplicates too, but only in terms of indices; catching them here
	// lets the message name the ids the user actually wrote.
	std::unordered_map<int64_t, size_t> seen;
	seen.reserve(rows.size());

	for (const auto &r : rows) {
		const int64_t si = samples.index.at(r.sample_id);
		const int64_t fi = features.index.at(r.feature_id);
		const int64_t cell = si * t.n_cols + fi;
		if (!seen.emplace(cell, 0).second) {
			throw std::invalid_argument("mmvec: " + label + " has a duplicate entry for sample '" + r.sample_id +
			                            "', feature '" + r.feature_id +
			                            "'; each (sample, feature) cell must appear at most once (summing them would "
			                            "change the fit, so they are not merged)");
		}
		t.row.push_back(si);
		t.col.push_back(fi);
		t.val.push_back(r.count);
	}
	return t;
}

//! Index one held-out table into a FITTED MODEL's feature dictionary, over an
//! already-built sample dictionary. Unchanged from the fitting path except for
//! where the feature dictionary comes from, which is exactly the intended
//! difference: same unknown-feature rejection, same duplicate-cell messages.
SparseCounts AlignOne(const std::vector<miint::unifrac::CooRow> &rows, const Dictionary &samples,
                      const std::vector<std::string> &feature_ids, const std::string &label) {
	const auto features = DictionaryFromIds(feature_ids);
	RejectUnknownFeatures(rows, features, label);
	return Index(rows, samples, features, label);
}

} // namespace

IngestedTables IngestPairedTables(const std::vector<miint::unifrac::CooRow> &x_rows,
                                  const std::vector<miint::unifrac::CooRow> &y_rows) {
	if (x_rows.empty()) {
		throw std::invalid_argument("mmvec: the X feature-table has no usable cells");
	}
	if (y_rows.empty()) {
		throw std::invalid_argument("mmvec: the Y feature-table has no usable cells");
	}

	const auto x_samples = BuildDictionary(x_rows, [](const miint::unifrac::CooRow &r) { return r.sample_id; });
	const auto y_samples = BuildDictionary(y_rows, [](const miint::unifrac::CooRow &r) { return r.sample_id; });

	// Matched by id, so a mismatch is an error rather than an intersection:
	// dropping a sample silently would change the sufficient statistics, and a
	// user whose two tables disagree wants to hear about it.
	RequireSameSamples(x_samples, y_samples);

	const auto x_features = BuildDictionary(x_rows, [](const miint::unifrac::CooRow &r) { return r.feature_id; });
	const auto y_features = BuildDictionary(y_rows, [](const miint::unifrac::CooRow &r) { return r.feature_id; });

	// Checked here rather than left to the core so the message can talk about Y
	// features instead of a column count.
	if (y_features.ids.size() < 2) {
		throw std::invalid_argument("mmvec: the Y feature-table needs at least two features (got " +
		                            std::to_string(y_features.ids.size()) + "); feature '" +
		                            (y_features.ids.empty() ? std::string("?") : y_features.ids[0]) +
		                            "' would be the reference category, whose logit is pinned to zero, leaving no "
		                            "likelihood to fit");
	}

	IngestedTables out;
	out.sample_ids = x_samples.ids;
	out.x_feature_ids = x_features.ids;
	out.y_feature_ids = y_features.ids;
	out.x = Index(x_rows, x_samples, x_features, "X");
	out.y = Index(y_rows, y_samples, y_features, "Y");
	return out;
}

AlignedTable AlignXToModel(const std::vector<miint::unifrac::CooRow> &x_rows,
                           const std::vector<std::string> &x_feature_ids) {
	if (x_rows.empty()) {
		throw std::invalid_argument("mmvec: the X feature-table has no usable cells");
	}
	const auto samples = BuildDictionary(x_rows, [](const miint::unifrac::CooRow &r) { return r.sample_id; });

	AlignedTable out;
	out.sample_ids = samples.ids;
	out.counts = AlignOne(x_rows, samples, x_feature_ids, "X");
	return out;
}

AlignedPair AlignPairedToModel(const std::vector<miint::unifrac::CooRow> &x_rows,
                               const std::vector<miint::unifrac::CooRow> &y_rows,
                               const std::vector<std::string> &x_feature_ids,
                               const std::vector<std::string> &y_feature_ids) {
	if (x_rows.empty()) {
		throw std::invalid_argument("mmvec: the X feature-table has no usable cells");
	}
	if (y_rows.empty()) {
		throw std::invalid_argument("mmvec: the Y feature-table has no usable cells");
	}
	const auto x_samples = BuildDictionary(x_rows, [](const miint::unifrac::CooRow &r) { return r.sample_id; });
	const auto y_samples = BuildDictionary(y_rows, [](const miint::unifrac::CooRow &r) { return r.sample_id; });
	RequireSameSamples(x_samples, y_samples);

	AlignedPair out;
	out.sample_ids = x_samples.ids;
	// One dictionary indexes both tables. The check above has just proved the two
	// are equal, so this is not what makes the pairing correct -- it is what keeps
	// it correct if that check is ever relaxed.
	out.x = AlignOne(x_rows, x_samples, x_feature_ids, "X");
	out.y = AlignOne(y_rows, x_samples, y_feature_ids, "Y");
	return out;
}

std::vector<ModelRow> BuildModelRows(const Model &model) {
	const int64_t d1 = model.shape.n_features_x;
	const int64_t d2 = model.shape.n_features_y;
	const int64_t p = model.shape.n_components;
	if (model.theta.size() != static_cast<size_t>(NumParams(model.shape))) {
		throw std::invalid_argument("mmvec: the model's theta has " + std::to_string(model.theta.size()) +
		                            " parameters, expected " + std::to_string(NumParams(model.shape)) +
		                            " for its own shape");
	}

	const auto [x_main, x_bias, y_main, y_bias] = BlockOffsets(d1, d2, p);

	std::vector<ModelRow> rows;
	rows.reserve(static_cast<size_t>((d1 + d2) * (p + 1)) + model.loss_curve.size());
	const auto push = [&rows](ModelRow::Kind kind, int64_t id, int64_t axis, double value) {
		rows.push_back({kind, id, static_cast<int32_t>(axis), value});
	};

	for (int64_t i = 0; i < d1; ++i) {
		for (int64_t k = 0; k < p; ++k) {
			push(ModelRow::Kind::X, i, k + 1, model.theta[x_main + static_cast<size_t>(i * p + k)]);
		}
		push(ModelRow::Kind::X, i, 0, model.theta[x_bias + static_cast<size_t>(i)]);
	}

	for (int64_t j = 0; j < d2; ++j) {
		for (int64_t k = 0; k < p; ++k) {
			// y_main is (p x d2-1) ROW-MAJOR, so component k of non-reference
			// feature j lives at k*(d2-1) + (j-1) -- not (j-1)*p + k. Both are in
			// range, which is why the transposed read is a silent error.
			const double v = j == 0 ? 0.0 : model.theta[y_main + static_cast<size_t>(k * (d2 - 1) + (j - 1))];
			push(ModelRow::Kind::Y, j, k + 1, v);
		}
		push(ModelRow::Kind::Y, j, 0, j == 0 ? 0.0 : model.theta[y_bias + static_cast<size_t>(j - 1)]);
	}

	for (size_t e = 0; e < model.loss_curve.size(); ++e) {
		push(ModelRow::Kind::Loss, 0, static_cast<int64_t>(e) + 1, model.loss_curve[e]);
	}
	return rows;
}

namespace {

//! One modality's cells, verified complete and laid out as an (n_features x (p+1))
//! grid with `axis` as the column, so reading a parameter back is a plain lookup.
struct Grid {
	std::vector<std::string> ids; //!< lexicographically sorted, deduplicated
	int32_t p = 0;                //!< the number of embedding axes; axes run 0..p
	std::vector<double> value;    //!< ids.size() x (p+1), row-major
};

//! The parameter for feature `i` on `axis`. Only valid once GatherGrid has proved
//! the grid complete -- which is why that check is not optional.
double GridAt(const Grid &g, int64_t i, int32_t axis) {
	return g.value[static_cast<size_t>(i * (g.p + 1) + axis)];
}

Grid GatherGrid(const std::vector<ModelCell> &cells, const char *modality) {
	// (feature index, axis) -> value, gathered first and validated on the sorted
	// order. Deliberately NOT a dense (n_features x (p+1)) array filled as we go:
	// `p` comes from the user's own `axis` column, so allocating from it before
	// proving it plausible would turn one bogus axis value into a huge allocation.
	struct Entry {
		int64_t feature;
		int32_t axis;
		double value;
	};
	std::vector<Entry> entries;
	std::vector<std::string> ids;
	int32_t max_axis = -1;
	for (const auto &c : cells) {
		if (c.modality != modality) {
			continue;
		}
		if (!c.feature_id.has_value()) {
			throw std::invalid_argument(std::string("mmvec: a modality '") + modality +
			                            "' row of the model relation has a NULL feature id");
		}
		if (c.axis < 0) {
			throw std::invalid_argument(std::string("mmvec: modality '") + modality + "' feature '" + *c.feature_id +
			                            "' has a negative axis " + std::to_string(c.axis));
		}
		if (!std::isfinite(c.value)) {
			throw std::invalid_argument(std::string("mmvec: modality '") + modality + "' feature '" + *c.feature_id +
			                            "' axis " + std::to_string(c.axis) + " is not finite");
		}
		ids.push_back(*c.feature_id);
		max_axis = std::max(max_axis, c.axis);
	}
	if (ids.empty()) {
		throw std::invalid_argument(std::string("mmvec: the model relation has no modality '") + modality +
		                            "' rows; it must carry both modalities of a fitted model");
	}

	Grid g;
	g.ids = std::move(ids);
	std::sort(g.ids.begin(), g.ids.end());
	g.ids.erase(std::unique(g.ids.begin(), g.ids.end()), g.ids.end());
	g.p = max_axis;

	// Every feature needs axes 0..p, so p+1 rows per feature and at least one
	// feature: an axis beyond the row count cannot possibly be filled. Checked
	// before any product is formed, so the arithmetic below cannot overflow.
	const size_t n_cells = static_cast<size_t>(
	    std::count_if(cells.begin(), cells.end(), [modality](const ModelCell &c) { return c.modality == modality; }));
	if (static_cast<size_t>(g.p) + 1 > n_cells) {
		throw std::invalid_argument(std::string("mmvec: modality '") + modality + "' has axes 0.." +
		                            std::to_string(g.p) + " but only " + std::to_string(n_cells) +
		                            " rows; the model relation cannot be complete");
	}

	std::unordered_map<std::string, int64_t> index;
	index.reserve(g.ids.size());
	for (size_t i = 0; i < g.ids.size(); ++i) {
		index.emplace(g.ids[i], static_cast<int64_t>(i));
	}
	entries.reserve(n_cells);
	for (const auto &c : cells) {
		if (c.modality == modality) {
			entries.push_back({index.at(*c.feature_id), c.axis, c.value});
		}
	}
	std::sort(entries.begin(), entries.end(), [](const Entry &a, const Entry &b) {
		return a.feature != b.feature ? a.feature < b.feature : a.axis < b.axis;
	});

	for (size_t k = 1; k < entries.size(); ++k) {
		if (entries[k].feature == entries[k - 1].feature && entries[k].axis == entries[k - 1].axis) {
			throw std::invalid_argument(std::string("mmvec: modality '") + modality + "' feature '" +
			                            g.ids[static_cast<size_t>(entries[k].feature)] +
			                            "' has more than one row for axis " + std::to_string(entries[k].axis));
		}
	}

	// With no duplicates and every axis in 0..p, having exactly ids.size()*(p+1)
	// entries forces the grid to be complete. Short of that it is not, and the walk
	// names the first hole -- a filtered or lossily joined model relation must not
	// be silently completed with a parameter the user never fitted.
	const size_t stride = static_cast<size_t>(g.p) + 1;
	if (entries.size() != g.ids.size() * stride) {
		for (size_t i = 0; i < g.ids.size(); ++i) {
			for (int32_t axis = 0; axis <= g.p; ++axis) {
				const size_t k = i * stride + static_cast<size_t>(axis);
				if (k >= entries.size() || entries[k].feature != static_cast<int64_t>(i) || entries[k].axis != axis) {
					throw std::invalid_argument(std::string("mmvec: modality '") + modality + "' feature '" + g.ids[i] +
					                            "' is missing axis " + std::to_string(axis) + " (axes 0.." +
					                            std::to_string(g.p) + " must all be present)");
				}
			}
		}
		throw std::invalid_argument(std::string("mmvec: modality '") + modality +
		                            "' has an inconsistent number of rows for its features");
	}

	g.value.reserve(entries.size());
	for (const auto &e : entries) {
		g.value.push_back(e.value);
	}
	return g;
}

} // namespace

ParsedModel ParseModelCells(const std::vector<ModelCell> &cells) {
	for (const auto &c : cells) {
		if (c.modality != "x" && c.modality != "y" && c.modality != "loss") {
			throw std::invalid_argument("mmvec: the model relation has modality '" + c.modality +
			                            "'; expected 'x', 'y' or 'loss'");
		}
	}

	const Grid x = GatherGrid(cells, "x");
	const Grid y = GatherGrid(cells, "y");
	if (x.p != y.p) {
		throw std::invalid_argument("mmvec: the model relation's two modalities disagree on the number of axes (x has "
		                            "0.." +
		                            std::to_string(x.p) + ", y has 0.." + std::to_string(y.p) + ")");
	}
	const int32_t p = x.p;
	const int64_t d1 = static_cast<int64_t>(x.ids.size());
	const int64_t d2 = static_cast<int64_t>(y.ids.size());
	if (p < 1) {
		throw std::invalid_argument("mmvec: the model relation has only axis 0 (the bias); a model needs at least one "
		                            "embedding axis");
	}
	if (d2 < 2) {
		throw std::invalid_argument("mmvec: the model relation has " + std::to_string(d2) +
		                            " Y feature(s); at least two are needed, since one is the reference category");
	}

	// The reference category holds no parameters, which BuildModelRows emits as an
	// all-zero row. Identified by that property rather than by id order, so
	// renaming features cannot silently move it -- see the header contract.
	std::vector<int64_t> all_zero;
	for (int64_t j = 0; j < d2; ++j) {
		bool zero = true;
		for (int32_t axis = 0; axis <= p && zero; ++axis) {
			zero = GridAt(y, j, axis) == 0.0;
		}
		if (zero) {
			all_zero.push_back(j);
		}
	}
	if (all_zero.size() != 1) {
		std::string detail;
		for (size_t k = 0; k < all_zero.size() && k < 3; ++k) {
			detail += (k == 0 ? " (" : ", ") + y.ids[static_cast<size_t>(all_zero[k])];
		}
		if (!detail.empty()) {
			detail += ")";
		}
		throw std::invalid_argument(
		    "mmvec: the model relation has " + std::to_string(all_zero.size()) + " all-zero Y features" + detail +
		    "; exactly one is required -- it is the reference category, whose logit is pinned to zero and which holds "
		    "no parameters");
	}

	// Reference first, then the rest in sorted order. Any order works as long as
	// theta and y_feature_ids agree, which they do by construction here.
	std::vector<int64_t> y_order;
	y_order.reserve(static_cast<size_t>(d2));
	y_order.push_back(all_zero[0]);
	for (int64_t j = 0; j < d2; ++j) {
		if (j != all_zero[0]) {
			y_order.push_back(j);
		}
	}

	ParsedModel out;
	out.shape = {d1, d2, p};
	out.x_feature_ids = x.ids;
	out.y_feature_ids.reserve(static_cast<size_t>(d2));
	for (const int64_t j : y_order) {
		out.y_feature_ids.push_back(y.ids[static_cast<size_t>(j)]);
	}

	const auto [x_main, x_bias, y_main, y_bias] = BlockOffsets(d1, d2, p);
	out.theta.assign(static_cast<size_t>(NumParams(out.shape)), 0.0);
	for (int64_t i = 0; i < d1; ++i) {
		for (int32_t k = 0; k < p; ++k) {
			out.theta[x_main + static_cast<size_t>(i * p + k)] = GridAt(x, i, k + 1);
		}
		out.theta[x_bias + static_cast<size_t>(i)] = GridAt(x, i, 0);
	}
	for (int64_t j = 1; j < d2; ++j) {
		const int64_t src = y_order[static_cast<size_t>(j)];
		for (int32_t k = 0; k < p; ++k) {
			// Mirrors BuildModelRows: y_main is (p x d2-1) row-major, so component k
			// of non-reference feature j lives at k*(d2-1) + (j-1).
			out.theta[y_main + static_cast<size_t>(k * (d2 - 1) + (j - 1))] = GridAt(y, src, k + 1);
		}
		out.theta[y_bias + static_cast<size_t>(j - 1)] = GridAt(y, src, 0);
	}
	return out;
}

} // namespace miint::mmvec
