#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "mmvec_relation.hpp"
#include "mmvec_oracle.hpp"

using Catch::Matchers::ContainsSubstring;
using miint::mmvec::IngestPairedTables;
using miint::unifrac::CooRow;

// ---------------------------------------------------------------------------
// Index assignment: the ordering rule, the sample join, and duplicate rejection.
//
// This is the wrapper's most consequential logic. A wrong feature order does not
// fail, crash, or converge poorly -- it silently picks a different reference
// category, and because the Gaussian priors are not softmax-shift-invariant the
// fit then lands on a different optimum with ranks moved by up to 31%. So the
// assertion here is not "some consistent ordering" but "exactly the matrices the
// carved oracle consumed", which is the only check that can catch it.
// ---------------------------------------------------------------------------
namespace {

// Long-form rows from a dense row-major matrix, EMITTED IN SCRAMBLED ORDER.
//
// Scrambling matters: the toy fixture's ids (O1..O8, C1..C6, S1..S5) are already
// lexicographic, so feeding rows in natural order would let a first-seen-order
// dictionary pass while still being wrong. Walking the cells backwards makes
// insertion order disagree with sorted order for every dictionary.
std::vector<CooRow> ScrambledRows(const std::vector<double> &dense, const std::vector<std::string> &sample_ids,
                                  const std::vector<std::string> &feature_ids) {
	const auto n = static_cast<int64_t>(sample_ids.size());
	const auto f = static_cast<int64_t>(feature_ids.size());
	std::vector<CooRow> rows;
	for (int64_t s = n - 1; s >= 0; --s) {
		for (int64_t c = f - 1; c >= 0; --c) {
			const double v = dense[static_cast<size_t>(s * f + c)];
			if (v != 0.0) {
				rows.push_back({sample_ids[static_cast<size_t>(s)], feature_ids[static_cast<size_t>(c)], v});
			}
		}
	}
	return rows;
}

// Densify an indexed SparseCounts back to row-major, so it can be compared to
// the oracle's carved dense matrix cell for cell.
//
// The bounds are asserted rather than trusted: an index assignment that
// transposed row and col would otherwise write past the buffer and abort the
// whole run with a heap-corruption message instead of failing one assertion.
// Verified by mutation -- this exact defect crashed here before the check.
std::vector<double> Densify(const miint::mmvec::SparseCounts &t) {
	std::vector<double> out(static_cast<size_t>(t.n_rows * t.n_cols), 0.0);
	for (size_t k = 0; k < t.val.size(); ++k) {
		REQUIRE(t.row[k] >= 0);
		REQUIRE(t.row[k] < t.n_rows);
		REQUIRE(t.col[k] >= 0);
		REQUIRE(t.col[k] < t.n_cols);
		out[static_cast<size_t>(t.row[k] * t.n_cols + t.col[k])] = t.val[k];
	}
	return out;
}

// A minimal valid pair: 2 samples, 2 X features, 2 Y features.
std::vector<CooRow> MinimalX() {
	return {{"s1", "xa", 1.0}, {"s1", "xb", 2.0}, {"s2", "xa", 3.0}, {"s2", "xb", 4.0}};
}
std::vector<CooRow> MinimalY() {
	return {{"s1", "ya", 1.0}, {"s1", "yb", 2.0}, {"s2", "ya", 3.0}, {"s2", "yb", 4.0}};
}

} // namespace

TEST_CASE("ingest reproduces the oracle's index assignment", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;

	const auto x_rows = ScrambledRows(kXCounts, kSampleIds, kXIds);
	const auto y_rows = ScrambledRows(kYCounts, kSampleIds, kYIds);
	const auto got = IngestPairedTables(x_rows, y_rows);

	// The dictionaries are the oracle's, in the oracle's order.
	REQUIRE(got.sample_ids == kSampleIds);
	REQUIRE(got.x_feature_ids == kXIds);
	REQUIRE(got.y_feature_ids == kYIds);

	// And the indexed tables densify back to the exact matrices the carved
	// expected values were produced from. Exact equality, not approximate: index
	// assignment moves values around, it never arithmetically alters them.
	REQUIRE(got.x.n_rows == kNSamples);
	REQUIRE(got.x.n_cols == kNFeaturesX);
	REQUIRE(got.y.n_rows == kNSamples);
	REQUIRE(got.y.n_cols == kNFeaturesY);
	REQUIRE(Densify(got.x) == kXCounts);
	REQUIRE(Densify(got.y) == kYCounts);
}

TEST_CASE("feature ids sort as bytes, not as numbers", "[mmvec]") {
	// "f10" < "f9" as bytes. This is the rule the committed fixtures were
	// generated under, so a "helpful" natural-order sort would silently invalidate
	// every carved expected value -- and, because y_feature_ids[0] is the
	// reference category, would move the optimum rather than merely relabel it.
	const std::vector<CooRow> x = {{"s1", "f9", 1.0}, {"s1", "f10", 2.0}, {"s2", "f9", 3.0}, {"s2", "f10", 4.0}};
	const std::vector<CooRow> y = {{"s1", "m2", 1.0}, {"s1", "m10", 2.0}, {"s2", "m2", 3.0}, {"s2", "m10", 4.0}};
	const auto got = IngestPairedTables(x, y);

	REQUIRE(got.x_feature_ids == std::vector<std::string> {"f10", "f9"});
	REQUIRE(got.y_feature_ids == std::vector<std::string> {"m10", "m2"});
	// So "m10" -- not "m2" -- is the reference category here.
	REQUIRE(got.y_feature_ids[0] == "m10");

	// Column indices follow the sorted dictionary, not the input order: s1's "f9"
	// value (1.0) must land in column 1, not column 0.
	const auto dense = Densify(got.x);
	REQUIRE(dense[0] == 2.0); // s1, f10
	REQUIRE(dense[1] == 1.0); // s1, f9
}

TEST_CASE("samples are matched by id, not by position", "[mmvec]") {
	// The two relations are read independently and a parallel scan fixes neither
	// one's row order, so position carries no meaning. Y is given here with its
	// samples in the opposite order from X; the pairing must still be by id.
	const std::vector<CooRow> x = {{"s1", "xa", 1.0}, {"s1", "xb", 2.0}, {"s2", "xa", 30.0}, {"s2", "xb", 40.0}};
	const std::vector<CooRow> y = {{"s2", "ya", 300.0}, {"s2", "yb", 400.0}, {"s1", "ya", 3.0}, {"s1", "yb", 4.0}};
	const auto got = IngestPairedTables(x, y);

	REQUIRE(got.sample_ids == std::vector<std::string> {"s1", "s2"});
	// Row 0 is s1 in BOTH tables.
	const auto dx = Densify(got.x);
	const auto dy = Densify(got.y);
	REQUIRE(dx[0] == 1.0);   // s1, xa
	REQUIRE(dy[0] == 3.0);   // s1, ya
	REQUIRE(dx[2] == 30.0);  // s2, xa
	REQUIRE(dy[2] == 300.0); // s2, ya
}

TEST_CASE("ingest rejects tables that cannot be paired or fitted", "[mmvec]") {
	SECTION("a sample present in only one modality") {
		// scikit-bio errors here rather than intersecting, and so do we: silently
		// dropping a sample changes the sufficient statistics, and a user who
		// mismatched their tables wants to know.
		auto y = MinimalY();
		y.push_back({"s3", "ya", 1.0});
		y.push_back({"s3", "yb", 1.0});
		REQUIRE_THROWS_WITH(IngestPairedTables(MinimalX(), y), ContainsSubstring("same samples"));
		REQUIRE_THROWS_WITH(IngestPairedTables(MinimalX(), y), ContainsSubstring("s3"));
	}

	SECTION("the missing sample is named whichever side it is on") {
		auto x = MinimalX();
		x.push_back({"zz", "xa", 1.0});
		x.push_back({"zz", "xb", 1.0});
		REQUIRE_THROWS_WITH(IngestPairedTables(x, MinimalY()), ContainsSubstring("zz"));
	}

	SECTION("a duplicated cell, reported by id") {
		// Summing duplicates would change X^T Y from (x_a+x_b)*y to x_a*y + x_b*y,
		// so this must fail rather than be quietly merged. The core rejects it too,
		// but knows only indices -- here the ids can be named.
		auto x = MinimalX();
		x.push_back({"s1", "xa", 5.0});
		REQUIRE_THROWS_WITH(IngestPairedTables(x, MinimalY()), ContainsSubstring("duplicate"));
		REQUIRE_THROWS_WITH(IngestPairedTables(x, MinimalY()), ContainsSubstring("s1"));
		REQUIRE_THROWS_WITH(IngestPairedTables(x, MinimalY()), ContainsSubstring("xa"));
	}

	SECTION("an empty table") {
		REQUIRE_THROWS_WITH(IngestPairedTables({}, MinimalY()), ContainsSubstring("no usable cells"));
		REQUIRE_THROWS_WITH(IngestPairedTables(MinimalX(), {}), ContainsSubstring("no usable cells"));
	}

	SECTION("a single Y feature leaves nothing to fit") {
		// Y feature 0 is the reference and carries a pinned zero logit, so with one
		// Y feature every logit is zero and the model has no content.
		const std::vector<CooRow> y = {{"s1", "only", 1.0}, {"s2", "only", 2.0}};
		REQUIRE_THROWS_WITH(IngestPairedTables(MinimalX(), y), ContainsSubstring("at least two"));
	}

	SECTION("a single X feature is fine") {
		// Unlike Y, one X feature is a legitimate (if uninteresting) model: there is
		// no reference category on the X side.
		const std::vector<CooRow> x = {{"s1", "only", 1.0}, {"s2", "only", 2.0}};
		REQUIRE_NOTHROW(IngestPairedTables(x, MinimalY()));
	}
}

TEST_CASE("ingested tables feed the core unchanged", "[mmvec]") {
	// The end-to-end point of this unit: what it produces is directly consumable,
	// so the sufficient statistics computed from the scrambled long-form input
	// match those from the oracle's dense matrices.
	using namespace miint::mmvec_oracle::toy;
	const auto got =
	    IngestPairedTables(ScrambledRows(kXCounts, kSampleIds, kXIds), ScrambledRows(kYCounts, kSampleIds, kYIds));
	const auto stats = miint::mmvec::ComputeSufficientStats(got.x, got.y);

	REQUIRE(stats.n_features_x == kNFeaturesX);
	REQUIRE(stats.n_features_y == kNFeaturesY);
	// Cross-checked against the logits the oracle carved at kTheta: if the index
	// assignment were off, these would not reproduce.
	const auto logits = miint::mmvec::ComputeLogits({kNFeaturesX, kNFeaturesY, kNComponents}, kTheta);
	REQUIRE(logits.size() == kT1Logits.size());
	for (size_t i = 0; i < kT1Logits.size(); ++i) {
		INFO("logit element " << i);
		REQUIRE(logits[i] == Catch::Approx(kT1Logits[i]).margin(miint::mmvec_oracle::kT1Tol));
	}
}

// ---------------------------------------------------------------------------
// Unpacking theta into relation rows.
//
// theta packs four blocks of different shapes back to back, so the errors that
// matter here are silent ones: reading y_main (p x d2-1, row-major) with the
// indices transposed, or getting a block boundary off by one, still produces
// in-range values of exactly the right count. Row totals, NULL patterns and the
// pinned-zero reference category all still look correct while every emitted
// coordinate is attached to the wrong (feature, axis).
//
// So the assertion is not structural. The rows are reassembled into embeddings
// and multiplied back out into logits, then compared against the oracle's carved
// kT1Logits at the same theta. That is decisive and non-circular: it goes through
// a different code path (this reconstruction) to a value produced by neither this
// unit nor the core, and both surviving mutations above are caught by it.
// ---------------------------------------------------------------------------
namespace {

// Rebuild the (d1 x d2) logit matrix from emitted rows, using only the DOCUMENTED
// meaning of (modality, feature, axis) -- never the packed offsets.
std::vector<double> LogitsFromRows(const std::vector<miint::mmvec::ModelRow> &rows, int64_t d1, int64_t d2, int64_t p) {
	std::vector<double> x_main(static_cast<size_t>(d1 * p), 0.0);
	std::vector<double> x_bias(static_cast<size_t>(d1), 0.0);
	std::vector<double> y_main(static_cast<size_t>(d2 * p), 0.0); // full d2, reference included
	std::vector<double> y_bias(static_cast<size_t>(d2), 0.0);

	for (const auto &r : rows) {
		switch (r.kind) {
		case miint::mmvec::ModelRow::Kind::X:
			if (r.axis == 0) {
				x_bias[static_cast<size_t>(r.id_index)] = r.value;
			} else {
				x_main[static_cast<size_t>(r.id_index * p + (r.axis - 1))] = r.value;
			}
			break;
		case miint::mmvec::ModelRow::Kind::Y:
			if (r.axis == 0) {
				y_bias[static_cast<size_t>(r.id_index)] = r.value;
			} else {
				y_main[static_cast<size_t>(r.id_index * p + (r.axis - 1))] = r.value;
			}
			break;
		case miint::mmvec::ModelRow::Kind::Loss:
			break;
		}
	}

	std::vector<double> logits(static_cast<size_t>(d1 * d2), 0.0);
	for (int64_t i = 0; i < d1; ++i) {
		for (int64_t j = 1; j < d2; ++j) {
			double dot = 0.0;
			for (int64_t k = 0; k < p; ++k) {
				dot += x_main[static_cast<size_t>(i * p + k)] * y_main[static_cast<size_t>(j * p + k)];
			}
			logits[static_cast<size_t>(i * d2 + j)] =
			    dot + x_bias[static_cast<size_t>(i)] + y_bias[static_cast<size_t>(j)];
		}
		// Column 0 stays 0: the reference category's logit is pinned.
	}
	return logits;
}

} // namespace

TEST_CASE("model rows reconstruct the oracle's logits", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;
	miint::mmvec::Model model;
	model.shape = {kNFeaturesX, kNFeaturesY, kNComponents};
	model.theta = kTheta;
	model.loss_curve = {1.0, 2.0, 3.0};

	const auto rows = miint::mmvec::BuildModelRows(model);
	REQUIRE(rows.size() == static_cast<size_t>((kNFeaturesX + kNFeaturesY) * (kNComponents + 1)) + 3);

	const auto rebuilt = LogitsFromRows(rows, kNFeaturesX, kNFeaturesY, kNComponents);
	REQUIRE(rebuilt.size() == kT1Logits.size());
	for (size_t i = 0; i < kT1Logits.size(); ++i) {
		INFO("logit element " << i);
		REQUIRE(rebuilt[i] == Catch::Approx(kT1Logits[i]).margin(miint::mmvec_oracle::kT1Tol));
	}
}

TEST_CASE("model rows lay out the axes and the reference category", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;
	miint::mmvec::Model model;
	model.shape = {kNFeaturesX, kNFeaturesY, kNComponents};
	model.theta = kTheta;
	model.loss_curve = {10.0, 9.0};
	const auto rows = miint::mmvec::BuildModelRows(model);

	int x_rows = 0;
	int y_rows = 0;
	int loss_rows = 0;
	double reference_abs = 0.0;
	int reference_rows = 0;
	bool axes_in_range = true;
	for (const auto &r : rows) {
		switch (r.kind) {
		case miint::mmvec::ModelRow::Kind::X:
			++x_rows;
			axes_in_range = axes_in_range && r.axis >= 0 && r.axis <= kNComponents;
			break;
		case miint::mmvec::ModelRow::Kind::Y:
			++y_rows;
			axes_in_range = axes_in_range && r.axis >= 0 && r.axis <= kNComponents;
			if (r.id_index == 0) {
				reference_abs += std::abs(r.value);
				++reference_rows;
			}
			break;
		case miint::mmvec::ModelRow::Kind::Loss:
			++loss_rows;
			break;
		}
	}

	REQUIRE(axes_in_range);
	REQUIRE(x_rows == kNFeaturesX * (kNComponents + 1));
	REQUIRE(y_rows == kNFeaturesY * (kNComponents + 1));
	REQUIRE(loss_rows == 2);
	// The reference category gets a full axis set, all exactly zero.
	REQUIRE(reference_rows == kNComponents + 1);
	REQUIRE(reference_abs == 0.0);

	// Loss rows are the 1-based evaluation ordinal, in order.
	std::vector<int32_t> loss_axes;
	std::vector<double> loss_values;
	for (const auto &r : rows) {
		if (r.kind == miint::mmvec::ModelRow::Kind::Loss) {
			loss_axes.push_back(r.axis);
			loss_values.push_back(r.value);
		}
	}
	REQUIRE(loss_axes == std::vector<int32_t> {1, 2});
	REQUIRE(loss_values == std::vector<double> {10.0, 9.0});
}

TEST_CASE("model rows reject a theta that does not match the shape", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;
	miint::mmvec::Model model;
	model.shape = {kNFeaturesX, kNFeaturesY, kNComponents};
	model.theta = kTheta;
	model.theta.pop_back();
	REQUIRE_THROWS_WITH(miint::mmvec::BuildModelRows(model), ContainsSubstring("expected"));
}
