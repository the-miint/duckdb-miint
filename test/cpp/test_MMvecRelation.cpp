#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
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

// ---------------------------------------------------------------------------
// Reading a model relation back: ParseModelCells, the inverse of BuildModelRows.
//
// Two independent assertions are needed, and neither replaces the other:
//
//  1. ROUND-TRIP IDENTITY on theta. Catches an inverse that disagrees with the
//     forward direction -- a transposed y_main read, a shifted axis.
//  2. The oracle's carved kT1Logits, reached from the PARSED theta. Catches what
//     the round-trip structurally cannot: the two directions now share one
//     spelling of theta's block offsets (BlockOffsets), so a wrong offset is wrong
//     consistently and CANCELS in the round-trip. Only an external value sees it.
//
// Neither is a self-consistency check on its own terms, which is the whole point:
// theta packs four differently-shaped blocks, so a mislabelled coordinate stays in
// range and produces a plausible-looking model of something else.
// ---------------------------------------------------------------------------
namespace {

using miint::mmvec::ModelCell;

// The model relation as SQL would present it: BuildModelRows' indices resolved
// through the id dictionaries, exactly as the wrapper will do it.
std::vector<ModelCell> CellsFromRows(const std::vector<miint::mmvec::ModelRow> &rows,
                                     const std::vector<std::string> &x_ids, const std::vector<std::string> &y_ids) {
	std::vector<ModelCell> cells;
	cells.reserve(rows.size());
	for (const auto &r : rows) {
		// Deliberately routed through ModalityName rather than spelling the three
		// strings again here: a private copy would let a renamed modality drift the
		// emitter and the parser apart with this very round trip still passing.
		const std::string modality = miint::mmvec::ModalityName(r.kind);
		switch (r.kind) {
		case miint::mmvec::ModelRow::Kind::X:
			cells.push_back({modality, x_ids[static_cast<size_t>(r.id_index)], r.axis, r.value});
			break;
		case miint::mmvec::ModelRow::Kind::Y:
			cells.push_back({modality, y_ids[static_cast<size_t>(r.id_index)], r.axis, r.value});
			break;
		case miint::mmvec::ModelRow::Kind::Loss:
			cells.push_back({modality, std::nullopt, r.axis, r.value});
			break;
		}
	}
	return cells;
}

// The oracle's toy model, as a model relation.
std::vector<ModelCell> ToyModelCells() {
	using namespace miint::mmvec_oracle::toy;
	miint::mmvec::Model model;
	model.shape = {kNFeaturesX, kNFeaturesY, kNComponents};
	model.theta = kTheta;
	model.loss_curve = {5.0, 4.0};
	return CellsFromRows(miint::mmvec::BuildModelRows(model), kXIds, kYIds);
}

// The same relation, but written from kTheta using the DOCUMENTED block layout
// spelled out HERE, independently of mmvec_relation.cpp. That independence is the
// whole value: the round-trip test above cannot see an error that is present in
// both directions and therefore cancels, because the two directions share one
// spelling of the layout. This does not go through BuildModelRows at all.
std::vector<ModelCell> ToyModelCellsFromTheta() {
	using namespace miint::mmvec_oracle::toy;
	const int64_t d1 = kNFeaturesX;
	const int64_t d2 = kNFeaturesY;
	const int64_t p = kNComponents;
	const size_t x_main = 0;
	const size_t x_bias = static_cast<size_t>(d1 * p);
	const size_t y_main = x_bias + static_cast<size_t>(d1);
	const size_t y_bias = y_main + static_cast<size_t>(p * (d2 - 1));

	std::vector<ModelCell> cells;
	for (int64_t i = 0; i < d1; ++i) {
		for (int64_t k = 0; k < p; ++k) {
			cells.push_back({"x", kXIds[static_cast<size_t>(i)], static_cast<int32_t>(k + 1),
			                 kTheta[x_main + static_cast<size_t>(i * p + k)]});
		}
		cells.push_back({"x", kXIds[static_cast<size_t>(i)], 0, kTheta[x_bias + static_cast<size_t>(i)]});
	}
	for (int64_t j = 0; j < d2; ++j) {
		for (int64_t k = 0; k < p; ++k) {
			const double v = j == 0 ? 0.0 : kTheta[y_main + static_cast<size_t>(k * (d2 - 1) + (j - 1))];
			cells.push_back({"y", kYIds[static_cast<size_t>(j)], static_cast<int32_t>(k + 1), v});
		}
		cells.push_back(
		    {"y", kYIds[static_cast<size_t>(j)], 0, j == 0 ? 0.0 : kTheta[y_bias + static_cast<size_t>(j - 1)]});
	}
	return cells;
}

} // namespace

TEST_CASE("a model relation round-trips through theta exactly", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;
	const auto parsed = miint::mmvec::ParseModelCells(ToyModelCells());

	REQUIRE(parsed.shape.n_features_x == kNFeaturesX);
	REQUIRE(parsed.shape.n_features_y == kNFeaturesY);
	REQUIRE(parsed.shape.n_components == kNComponents);
	REQUIRE(parsed.x_feature_ids == kXIds);
	// The reference lands at index 0; here it is also the lexicographically first
	// id, because that is what mmvec_fit's ordering rule made it.
	REQUIRE(parsed.y_feature_ids == kYIds);

	// EXACT, not approximate: nothing in this path does arithmetic on a parameter,
	// so anything other than bit equality is a bug rather than rounding.
	REQUIRE(parsed.theta.size() == kTheta.size());
	for (size_t i = 0; i < kTheta.size(); ++i) {
		INFO("theta element " << i);
		REQUIRE(parsed.theta[i] == kTheta[i]);
	}
}

TEST_CASE("parsing a relation written from theta recovers theta", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;
	// The anchor. Cells come from the independent layout above, so BuildModelRows is
	// not in the loop and cannot cancel an error in ParseModelCells. A transposed
	// y_main read applied to BOTH directions passes the round-trip test and fails
	// here -- verified by mutation, not assumed.
	const auto parsed = miint::mmvec::ParseModelCells(ToyModelCellsFromTheta());
	REQUIRE(parsed.shape.n_features_x == kNFeaturesX);
	REQUIRE(parsed.shape.n_features_y == kNFeaturesY);
	REQUIRE(parsed.shape.n_components == kNComponents);
	REQUIRE(parsed.theta.size() == kTheta.size());
	for (size_t i = 0; i < kTheta.size(); ++i) {
		INFO("theta element " << i);
		REQUIRE(parsed.theta[i] == kTheta[i]);
	}
	// And the recovered model reproduces the carved logits, which is what the SQL
	// layer will actually read it through.
	const auto logits = miint::mmvec::ComputeLogits(parsed.shape, parsed.theta);
	REQUIRE(logits.size() == kT1Logits.size());
	for (size_t i = 0; i < kT1Logits.size(); ++i) {
		INFO("logit element " << i);
		REQUIRE(logits[i] == Catch::Approx(kT1Logits[i]).margin(miint::mmvec_oracle::kT1Tol));
	}
}

TEST_CASE("both ways of writing the toy model relation agree", "[mmvec]") {
	// If these two ever disagree, one of the two layout spellings is wrong and the
	// tests above will be measuring the wrong thing. Cheap, and it pins the
	// independence rather than leaving it as a claim in a comment.
	const auto from_rows = ToyModelCells();
	const auto from_theta = ToyModelCellsFromTheta();
	std::vector<std::string> a;
	std::vector<std::string> b;
	for (const auto &c : from_rows) {
		if (c.modality != "loss") {
			a.push_back(c.modality + "|" + *c.feature_id + "|" + std::to_string(c.axis) + "|" +
			            std::to_string(c.value));
		}
	}
	for (const auto &c : from_theta) {
		b.push_back(c.modality + "|" + *c.feature_id + "|" + std::to_string(c.axis) + "|" + std::to_string(c.value));
	}
	std::sort(a.begin(), a.end());
	std::sort(b.begin(), b.end());
	REQUIRE(a == b);
}

TEST_CASE("parsing a model ignores the loss curve and tolerates its absence", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;
	auto with_loss = ToyModelCells();
	std::vector<ModelCell> without_loss;
	for (const auto &c : with_loss) {
		if (c.modality != "loss") {
			without_loss.push_back(c);
		}
	}
	REQUIRE(without_loss.size() + 2 == with_loss.size());

	const auto a = miint::mmvec::ParseModelCells(with_loss);
	const auto b = miint::mmvec::ParseModelCells(without_loss);
	REQUIRE(a.theta == b.theta);
	REQUIRE(a.y_feature_ids == b.y_feature_ids);
}

TEST_CASE("parsing a model does not depend on the order of the cells", "[mmvec]") {
	auto cells = ToyModelCells();
	const auto forward = miint::mmvec::ParseModelCells(cells);
	std::reverse(cells.begin(), cells.end());
	const auto reversed = miint::mmvec::ParseModelCells(cells);
	// A relation has no row order -- NO_ORDER is declared on the producing function
	// -- so a parse that depended on it would be wrong for reasons no fixture would
	// reveal, since mmvec_fit happens to emit in dictionary order.
	REQUIRE(forward.theta == reversed.theta);
	REQUIRE(forward.x_feature_ids == reversed.x_feature_ids);
	REQUIRE(forward.y_feature_ids == reversed.y_feature_ids);
}

TEST_CASE("the reference category is found by its zeros, not by its id", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;
	// Rename the Y features so the reference ('C1') is no longer lexicographically
	// first: C1 -> zz_ref sorts LAST. Joining a model to metadata and renaming
	// features to readable names is ordinary, and must not move the reference.
	auto cells = ToyModelCells();
	for (auto &c : cells) {
		if (c.modality == "y" && c.feature_id == std::string("C1")) {
			c.feature_id = "zz_ref";
		}
	}
	const auto parsed = miint::mmvec::ParseModelCells(cells);
	REQUIRE(parsed.y_feature_ids[0] == "zz_ref");
	// Renaming is presentational: the model itself is untouched.
	const auto original = miint::mmvec::ParseModelCells(ToyModelCells());
	REQUIRE(parsed.theta == original.theta);
	for (size_t j = 1; j < parsed.y_feature_ids.size(); ++j) {
		REQUIRE(parsed.y_feature_ids[j] == original.y_feature_ids[j]);
	}
}

TEST_CASE("parsing a model rejects a relation it cannot trust", "[mmvec]") {
	using namespace miint::mmvec_oracle::toy;

	SECTION("an unknown modality") {
		auto cells = ToyModelCells();
		cells[0].modality = "fit";
		REQUIRE_THROWS_WITH(miint::mmvec::ParseModelCells(cells), ContainsSubstring("modality 'fit'"));
	}

	SECTION("a NULL feature id on a parameter row") {
		auto cells = ToyModelCells();
		cells[0].feature_id = std::nullopt;
		REQUIRE_THROWS_WITH(miint::mmvec::ParseModelCells(cells), ContainsSubstring("NULL feature id"));
	}

	SECTION("a missing modality") {
		std::vector<ModelCell> only_x;
		for (const auto &c : ToyModelCells()) {
			if (c.modality == "x") {
				only_x.push_back(c);
			}
		}
		REQUIRE_THROWS_WITH(miint::mmvec::ParseModelCells(only_x), ContainsSubstring("no modality 'y' rows"));
	}

	SECTION("a duplicated cell") {
		auto cells = ToyModelCells();
		cells.push_back(cells[0]);
		REQUIRE_THROWS_WITH(miint::mmvec::ParseModelCells(cells), ContainsSubstring("more than one row for axis"));
	}

	SECTION("an incomplete grid names the feature and the axis") {
		auto cells = ToyModelCells();
		// Drop O3's axis 2 only. The count is still plausible and every other
		// feature is whole, which is exactly why this has to be checked per slot.
		cells.erase(std::remove_if(cells.begin(), cells.end(),
		                           [](const ModelCell &c) {
			                           return c.modality == "x" && c.feature_id == std::string("O3") && c.axis == 2;
		                           }),
		            cells.end());
		REQUIRE_THROWS_WITH(miint::mmvec::ParseModelCells(cells), ContainsSubstring("'O3' is missing axis 2"));
	}

	SECTION("the modalities disagree on the number of axes") {
		auto cells = ToyModelCells();
		for (const auto &c : ToyModelCells()) {
			if (c.modality == "y" && c.axis == 1) {
				ModelCell extra = c;
				extra.axis = kNComponents + 1;
				cells.push_back(extra);
			}
		}
		REQUIRE_THROWS_WITH(miint::mmvec::ParseModelCells(cells), ContainsSubstring("disagree on the number of axes"));
	}

	SECTION("a non-finite parameter") {
		auto cells = ToyModelCells();
		cells[0].value = std::numeric_limits<double>::infinity();
		REQUIRE_THROWS_WITH(miint::mmvec::ParseModelCells(cells), ContainsSubstring("is not finite"));
	}

	SECTION("a negative axis") {
		auto cells = ToyModelCells();
		cells[0].axis = -1;
		REQUIRE_THROWS_WITH(miint::mmvec::ParseModelCells(cells), ContainsSubstring("negative axis"));
	}

	SECTION("no reference category") {
		auto cells = ToyModelCells();
		for (auto &c : cells) {
			if (c.modality == "y" && c.feature_id == std::string("C1") && c.axis == 0) {
				c.value = 0.5; // C1 now holds a parameter, so nothing is pinned
			}
		}
		REQUIRE_THROWS_WITH(miint::mmvec::ParseModelCells(cells), ContainsSubstring("has 0 all-zero Y features"));
	}

	SECTION("two reference categories is ambiguous, not a guess") {
		auto cells = ToyModelCells();
		for (auto &c : cells) {
			if (c.modality == "y" && c.feature_id == std::string("C4")) {
				c.value = 0.0;
			}
		}
		REQUIRE_THROWS_WITH(miint::mmvec::ParseModelCells(cells), ContainsSubstring("has 2 all-zero Y features"));
		REQUIRE_THROWS_WITH(miint::mmvec::ParseModelCells(cells), ContainsSubstring("C4"));
	}

	SECTION("only a bias axis is not a model") {
		std::vector<ModelCell> bias_only;
		for (const auto &c : ToyModelCells()) {
			if (c.modality != "loss" && c.axis == 0) {
				bias_only.push_back(c);
			}
		}
		REQUIRE_THROWS_WITH(miint::mmvec::ParseModelCells(bias_only), ContainsSubstring("at least one embedding axis"));
	}

	SECTION("a single Y feature leaves nothing to model") {
		std::vector<ModelCell> one_y;
		for (const auto &c : ToyModelCells()) {
			if (c.modality == "x" || (c.modality == "y" && c.feature_id == std::string("C1"))) {
				one_y.push_back(c);
			}
		}
		REQUIRE_THROWS_WITH(miint::mmvec::ParseModelCells(one_y), ContainsSubstring("at least two are needed"));
	}
}

// ---------------------------------------------------------------------------
// Aligning a held-out table to a fitted model's dictionary.
//
// The counterpart to IngestPairedTables, and the opposite direction: ingest builds
// a dictionary FROM the data, here the dictionary is the model's and the data is
// fitted into it. Getting that backwards is a silent wrong answer -- in-range
// column indices of plausibly the right count, naming different features -- and
// the SQL-level guard for it is thin (measured: a fresh dictionary passes the
// entire exact-parity block in mmvec_predict.test, because that fixture's own
// dictionary happens to coincide with the model's). So it is pinned here directly.
// ---------------------------------------------------------------------------

TEST_CASE("aligning to a model keeps the model's width and column order", "[mmvec]") {
	// The model knows four features; the table carries only two of them, and the
	// table's own sorted dictionary would be ("xb", "xd") -- two columns at indices
	// 0 and 1. Aligned to the model they must be columns 1 and 3 of four.
	const std::vector<std::string> model_features = {"xa", "xb", "xc", "xd"};
	const std::vector<CooRow> rows = {{"s1", "xd", 4.0}, {"s1", "xb", 2.0}, {"s2", "xb", 20.0}};

	const auto got = miint::mmvec::AlignXToModel(rows, model_features);
	REQUIRE(got.sample_ids == std::vector<std::string> {"s1", "s2"});
	REQUIRE(got.counts.n_rows == 2);
	REQUIRE(got.counts.n_cols == 4); // the MODEL's width, not the table's two features

	const auto dense = Densify(got.counts);
	REQUIRE(dense == std::vector<double> {0.0, 2.0, 0.0, 4.0, 0.0, 20.0, 0.0, 0.0});
}

TEST_CASE("aligning follows the model's id order, without re-sorting it", "[mmvec]") {
	// The dictionary is taken AS GIVEN. That matters for the Y side, where the model
	// puts the reference category at index 0 regardless of its id -- a helpful
	// re-sort here would silently move the reference and change what every emitted
	// column means. X is exercised because it shares the code path.
	const std::vector<std::string> model_features = {"zz", "aa", "mm"};
	const std::vector<CooRow> rows = {{"s1", "aa", 1.0}, {"s1", "zz", 2.0}, {"s1", "mm", 3.0}};

	const auto got = miint::mmvec::AlignXToModel(rows, model_features);
	REQUIRE(got.counts.n_cols == 3);
	// Column order is zz, aa, mm -- the model's -- not the sorted aa, mm, zz.
	REQUIRE(Densify(got.counts) == std::vector<double> {2.0, 1.0, 3.0});
}

TEST_CASE("aligning rejects a table it cannot fit into the model", "[mmvec]") {
	const std::vector<std::string> model_features = {"xa", "xb"};

	SECTION("a feature the model never saw, named and counted") {
		const std::vector<CooRow> rows = {{"s1", "xa", 1.0}, {"s1", "surprise", 2.0}};
		REQUIRE_THROWS_WITH(miint::mmvec::AlignXToModel(rows, model_features),
		                    ContainsSubstring("1 feature(s) the model was never fitted on"));
		REQUIRE_THROWS_WITH(miint::mmvec::AlignXToModel(rows, model_features), ContainsSubstring("surprise"));
		// The message tells the user what to do about it, because a sample-wise
		// train/test split produces this case routinely.
		REQUIRE_THROWS_WITH(miint::mmvec::AlignXToModel(rows, model_features),
		                    ContainsSubstring("restrict the table to the model's own features"));
	}

	SECTION("several unknown features are counted once each, not per cell") {
		const std::vector<CooRow> rows = {{"s1", "p", 1.0}, {"s2", "p", 1.0}, {"s1", "q", 1.0}, {"s1", "xa", 1.0}};
		REQUIRE_THROWS_WITH(miint::mmvec::AlignXToModel(rows, model_features),
		                    ContainsSubstring("2 feature(s) the model was never fitted on"));
	}

	SECTION("an empty table") {
		REQUIRE_THROWS_WITH(miint::mmvec::AlignXToModel({}, model_features), ContainsSubstring("no usable cells"));
	}

	SECTION("a duplicated cell, with the same message the fitting path gives") {
		const std::vector<CooRow> rows = {{"s1", "xa", 1.0}, {"s1", "xa", 2.0}};
		REQUIRE_THROWS_WITH(miint::mmvec::AlignXToModel(rows, model_features), ContainsSubstring("duplicate entry"));
	}
}

// ---------------------------------------------------------------------------
// The paired form, for Score: the same alignment on both modalities, plus the one
// thing scoring adds -- X and Y must describe the SAME samples, over ONE sample
// dictionary. Score compares the two tables cell by cell, so two independently
// built sample orders would pair one sample's microbes against another's
// metabolites and report a plausible number for it.
// ---------------------------------------------------------------------------

TEST_CASE("pairing aligns both modalities into the model's own dictionaries", "[mmvec]") {
	// Neither table's own dictionary is the model's: X's sorted features are
	// (xa, xb, xc) but only two appear, and Y's model order (zz, aa, mm) is
	// deliberately unsorted so a re-sort would show up as moved columns. The Y case
	// is the one that matters -- index 0 is the reference category.
	const std::vector<std::string> x_model = {"xa", "xb", "xc"};
	const std::vector<std::string> y_model = {"zz", "aa", "mm"};
	const std::vector<CooRow> x_rows = {{"s1", "xc", 4.0}, {"s1", "xa", 1.0}, {"s2", "xb", 2.0}};
	const std::vector<CooRow> y_rows = {{"s1", "aa", 3.0}, {"s2", "zz", 5.0}, {"s2", "mm", 7.0}};

	const auto got = miint::mmvec::AlignPairedToModel(x_rows, y_rows, x_model, y_model);
	REQUIRE(got.sample_ids == std::vector<std::string> {"s1", "s2"});
	REQUIRE(got.x.n_rows == 2);
	REQUIRE(got.y.n_rows == 2);
	REQUIRE(got.x.n_cols == 3);
	REQUIRE(got.y.n_cols == 3);

	REQUIRE(Densify(got.x) == std::vector<double> {1.0, 0.0, 4.0, 0.0, 2.0, 0.0});
	// (zz, aa, mm), the model's order -- not the sorted (aa, mm, zz).
	REQUIRE(Densify(got.y) == std::vector<double> {0.0, 3.0, 0.0, 5.0, 0.0, 7.0});
}

TEST_CASE("pairing puts both modalities on ONE sample dictionary", "[mmvec]") {
	// Row n of x and row n of y must be the same sample. Written so that a per-table
	// dictionary would still produce two valid-looking tables of the right shape: Y's
	// rows arrive in the opposite order to X's, so pairing by arrival rather than by
	// id would transpose the two samples and quietly score the wrong pairs.
	const std::vector<std::string> x_model = {"xa"};
	const std::vector<std::string> y_model = {"ya", "yb"};
	const std::vector<CooRow> x_rows = {{"s1", "xa", 1.0}, {"s2", "xa", 2.0}};
	const std::vector<CooRow> y_rows = {{"s2", "yb", 20.0}, {"s1", "ya", 10.0}};

	const auto got = miint::mmvec::AlignPairedToModel(x_rows, y_rows, x_model, y_model);
	REQUIRE(got.sample_ids == std::vector<std::string> {"s1", "s2"});
	// s1 (row 0) holds ya, s2 (row 1) holds yb -- by id, not by arrival.
	REQUIRE(Densify(got.y) == std::vector<double> {10.0, 0.0, 0.0, 20.0});
	REQUIRE(Densify(got.x) == std::vector<double> {1.0, 2.0});
}

TEST_CASE("pairing rejects tables it cannot align", "[mmvec]") {
	const std::vector<std::string> x_model = {"xa", "xb"};
	const std::vector<std::string> y_model = {"ya", "yb"};
	const std::vector<CooRow> x_ok = {{"s1", "xa", 1.0}, {"s2", "xb", 2.0}};
	const std::vector<CooRow> y_ok = {{"s1", "ya", 1.0}, {"s2", "yb", 2.0}};

	SECTION("a sample present in only one modality, with the fitting path's message") {
		// Deliberately the SAME text IngestPairedTables produces: it is the same
		// mistake, and a user who has seen it once should not have to learn it twice.
		auto x = x_ok;
		x.push_back({"s3", "xa", 1.0});
		REQUIRE_THROWS_WITH(miint::mmvec::AlignPairedToModel(x, y_ok, x_model, y_model),
		                    ContainsSubstring("must describe the same samples"));
		REQUIRE_THROWS_WITH(miint::mmvec::AlignPairedToModel(x, y_ok, x_model, y_model),
		                    ContainsSubstring("only in X: s3"));

		auto y = y_ok;
		y.push_back({"s4", "ya", 1.0});
		REQUIRE_THROWS_WITH(miint::mmvec::AlignPairedToModel(x_ok, y, x_model, y_model),
		                    ContainsSubstring("only in Y: s4"));
	}

	SECTION("a mismatch of the same SIZE, which is the case this check exists for") {
		// Where the two counts differ, Score would notice by itself -- x.n_rows and
		// y.n_rows disagree. Here they do not, so nothing downstream can tell that s2's
		// microbes are about to be scored against s9's metabolites. Both sides are
		// named, not just the count.
		const std::vector<CooRow> y = {{"s1", "ya", 1.0}, {"s9", "yb", 2.0}};
		REQUIRE_THROWS_WITH(miint::mmvec::AlignPairedToModel(x_ok, y, x_model, y_model),
		                    ContainsSubstring("(X has 2, Y has 2)"));
		REQUIRE_THROWS_WITH(miint::mmvec::AlignPairedToModel(x_ok, y, x_model, y_model),
		                    ContainsSubstring("only in X: s2; only in Y: s9"));
	}

	SECTION("an unknown Y feature is rejected exactly as an unknown X feature is") {
		auto y = y_ok;
		y.push_back({"s1", "mystery", 1.0});
		REQUIRE_THROWS_WITH(miint::mmvec::AlignPairedToModel(x_ok, y, x_model, y_model),
		                    ContainsSubstring("1 feature(s) the model was never fitted on"));
		REQUIRE_THROWS_WITH(miint::mmvec::AlignPairedToModel(x_ok, y, x_model, y_model), ContainsSubstring("mystery"));
		// The remedy names the Y column and the Y modality, not X's. Getting this
		// wrong would send the user to filter the wrong table.
		REQUIRE_THROWS_WITH(miint::mmvec::AlignPairedToModel(x_ok, y, x_model, y_model),
		                    ContainsSubstring("SELECT DISTINCT y_feature_id"));
		REQUIRE_THROWS_WITH(miint::mmvec::AlignPairedToModel(x_ok, y, x_model, y_model),
		                    ContainsSubstring("modality = 'y'"));
	}

	SECTION("an unknown X feature, still") {
		auto x = x_ok;
		x.push_back({"s1", "surprise", 1.0});
		REQUIRE_THROWS_WITH(miint::mmvec::AlignPairedToModel(x, y_ok, x_model, y_model),
		                    ContainsSubstring("SELECT DISTINCT x_feature_id"));
	}

	SECTION("an empty table, either side") {
		REQUIRE_THROWS_WITH(miint::mmvec::AlignPairedToModel({}, y_ok, x_model, y_model),
		                    ContainsSubstring("the X feature-table has no usable cells"));
		REQUIRE_THROWS_WITH(miint::mmvec::AlignPairedToModel(x_ok, {}, x_model, y_model),
		                    ContainsSubstring("the Y feature-table has no usable cells"));
	}

	SECTION("a duplicated cell, labelled with the modality it is in") {
		auto y = y_ok;
		y.push_back({"s1", "ya", 5.0});
		REQUIRE_THROWS_WITH(miint::mmvec::AlignPairedToModel(x_ok, y, x_model, y_model),
		                    ContainsSubstring("Y has a duplicate entry for sample 's1', feature 'ya'"));
	}
}
