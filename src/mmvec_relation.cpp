#include "mmvec_relation.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_map>

namespace miint::mmvec {

namespace {

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
	if (x_samples.ids != y_samples.ids) {
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

std::vector<ModelRow> BuildModelRows(const Model &model) {
	const int64_t d1 = model.shape.n_features_x;
	const int64_t d2 = model.shape.n_features_y;
	const int64_t p = model.shape.n_components;
	if (model.theta.size() != static_cast<size_t>(NumParams(model.shape))) {
		throw std::invalid_argument("mmvec: the model's theta has " + std::to_string(model.theta.size()) +
		                            " parameters, expected " + std::to_string(NumParams(model.shape)) +
		                            " for its own shape");
	}

	// The packed layout, spelled out once: x_main (d1 x p), x_bias (d1),
	// y_main (p x d2-1), y_bias (d2-1), every block C-order row-major.
	const size_t x_main = 0;
	const size_t x_bias = static_cast<size_t>(d1 * p);
	const size_t y_main = x_bias + static_cast<size_t>(d1);
	const size_t y_bias = y_main + static_cast<size_t>(p * (d2 - 1));

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

} // namespace miint::mmvec
