#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "NewickTree.hpp"
#include "unifrac_bptree.hpp"
#include "api.hpp"
#include "unifrac_subsample_bridge.hpp"
#include "unifrac_support_biom.hpp"

#include <stdexcept>
#include <unordered_map>
#include <vector>

using miint::unifrac::BridgeSubsample;
using miint::unifrac::CooRow;
using miint::unifrac::UnifracBptreeView;
using miint::unifrac::UnifracSupportBiomView;

namespace {

// Sample totals — used to predict subsample drop behaviour.
//   S1: 6, S2: 6, S3: 8, S4: 4
// At subsample_depth=4 all four samples survive.
// At subsample_depth=7 only S3 survives (others below 7).
std::vector<CooRow> MakeRichRows() {
	return {
	    {"S1", "A", 3.0}, {"S1", "B", 3.0},                                     // total 6
	    {"S2", "A", 2.0}, {"S2", "B", 2.0}, {"S2", "C", 2.0},                   // total 6
	    {"S3", "A", 2.0}, {"S3", "B", 2.0}, {"S3", "C", 2.0}, {"S3", "D", 2.0}, // total 8
	    {"S4", "B", 1.0}, {"S4", "C", 1.0}, {"S4", "D", 1.0}, {"S4", "E", 1.0}, // total 4
	};
}

} // namespace

TEST_CASE("BridgeSubsample preserves all samples when their totals meet the depth", "[unifrac][subsample_bridge]") {
	auto biom = UnifracSupportBiomView::FromCoo(MakeRichRows());
	auto sub = BridgeSubsample(biom, /*subsample_depth*/ 4, /*with_replacement*/ false, /*seed*/ 42);

	const auto *b = sub.support_biom();
	REQUIRE(b->n_samples == 4); // S1..S4 all survive at depth 4
}

TEST_CASE("BridgeSubsample (no replacement): column sums equal subsample_depth", "[unifrac][subsample_bridge]") {
	// Without replacement, every surviving sample is rarefied to exactly
	// `subsample_depth` total counts — the defining invariant of the
	// subsample bridge that downstream Faith PD relies on.
	auto biom = UnifracSupportBiomView::FromCoo(MakeRichRows());
	const uint32_t depth = 4;
	auto sub = BridgeSubsample(biom, depth, /*with_replacement*/ false, /*seed*/ 42);

	const auto *b = sub.support_biom();
	std::unordered_map<int, double> column_sums;
	for (int o = 0; o < b->n_obs; ++o) {
		const uint32_t start = b->indptr[o];
		const uint32_t end = b->indptr[o + 1];
		for (uint32_t k = start; k < end; ++k) {
			column_sums[static_cast<int>(b->indices[k])] += b->data[k];
		}
	}
	for (int s = 0; s < b->n_samples; ++s) {
		REQUIRE(column_sums[s] == Catch::Approx(static_cast<double>(depth)));
	}
}

TEST_CASE("BridgeSubsample drops samples whose total counts fall below the depth", "[unifrac][subsample_bridge]") {
	auto biom = UnifracSupportBiomView::FromCoo(MakeRichRows());
	// At depth 7, only S3 (total 8) survives; S1/S2 (totals 6) and S4
	// (total 4) get dropped.
	auto sub = BridgeSubsample(biom, /*subsample_depth*/ 7, /*with_replacement*/ false, /*seed*/ 42);

	const auto *b = sub.support_biom();
	REQUIRE(b->n_samples == 1);
	REQUIRE(std::string(b->sample_ids[0]) == "S3");
}

TEST_CASE("BridgeSubsample is reproducible under the same seed", "[unifrac][subsample_bridge]") {
	// Two subsamples with the same seed must yield byte-identical CSR
	// arrays — subsample_table_inmem_seeded is per-call seeded so this
	// holds without any global RNG coordination.
	auto biom1 = UnifracSupportBiomView::FromCoo(MakeRichRows());
	auto biom2 = UnifracSupportBiomView::FromCoo(MakeRichRows());

	auto sub_a = BridgeSubsample(biom1, /*subsample_depth*/ 4, /*with_replacement*/ false, /*seed*/ 42);
	auto sub_b = BridgeSubsample(biom2, /*subsample_depth*/ 4, /*with_replacement*/ false, /*seed*/ 42);
	const auto *a = sub_a.support_biom();
	const auto *b = sub_b.support_biom();

	REQUIRE(a->n_samples == b->n_samples);
	REQUIRE(a->n_obs == b->n_obs);
	REQUIRE(a->nnz == b->nnz);
	for (int i = 0; i < a->nnz; ++i) {
		REQUIRE(a->indices[i] == b->indices[i]);
		REQUIRE(a->data[i] == b->data[i]);
	}
	for (int i = 0; i <= a->n_obs; ++i) {
		REQUIRE(a->indptr[i] == b->indptr[i]);
	}
}

TEST_CASE("BridgeSubsample with different seeds produces different output", "[unifrac][subsample_bridge]") {
	// Guard against a no-op bridge passing the same-seed test: two different
	// seeds at the same depth must produce CSR data that differs in at least
	// one cell. With depth=4 and 4–8 source counts per sample, the multinomial
	// draw has plenty of room for variation across seeds.
	auto biom = UnifracSupportBiomView::FromCoo(MakeRichRows());
	auto sub_a = BridgeSubsample(biom, /*subsample_depth*/ 4, /*with_replacement*/ false, /*seed*/ 42);
	auto sub_b = BridgeSubsample(biom, /*subsample_depth*/ 4, /*with_replacement*/ false, /*seed*/ 1337);
	const auto *a = sub_a.support_biom();
	const auto *b = sub_b.support_biom();

	REQUIRE(a->n_samples == b->n_samples);
	REQUIRE(a->n_obs == b->n_obs);

	bool differs = false;
	if (a->nnz != b->nnz) {
		differs = true;
	} else {
		for (int i = 0; i < a->nnz && !differs; ++i) {
			if (a->indices[i] != b->indices[i] || a->data[i] != b->data[i]) {
				differs = true;
			}
		}
	}
	REQUIRE(differs);
}

TEST_CASE("BridgeSubsample with_replacement=true: per-sample totals still equal depth", "[unifrac][subsample_bridge]") {
	// Multinomial sampling (with_replacement=true) draws exactly `depth` reads
	// per sample, so the column-sum invariant holds for the with-replacement
	// path too. This exercises biom_subsampled's WeightedSampleWithReplacement
	// branch, which the no-replacement test never reaches.
	auto biom = UnifracSupportBiomView::FromCoo(MakeRichRows());
	const uint32_t depth = 4;
	auto sub = BridgeSubsample(biom, depth, /*with_replacement*/ true, /*seed*/ 42);

	const auto *b = sub.support_biom();
	REQUIRE(b->n_samples == 4);
	std::unordered_map<int, double> column_sums;
	for (int o = 0; o < b->n_obs; ++o) {
		for (uint32_t k = b->indptr[o]; k < b->indptr[o + 1]; ++k) {
			column_sums[static_cast<int>(b->indices[k])] += b->data[k];
		}
	}
	for (int s = 0; s < b->n_samples; ++s) {
		REQUIRE(column_sums[s] == Catch::Approx(static_cast<double>(depth)));
	}
}

TEST_CASE("BridgeSubsample output is consumable by faith_pd_inmem", "[unifrac][subsample_bridge][libssu]") {
	// End-to-end smoke: the whole point of the bridge is that downstream
	// libssu APIs can consume the reconstructed support_biom_t. If the
	// canonical-ordering invariant is wrong, faith_pd_inmem returns
	// nonsense; this catches that immediately.
	auto biom = UnifracSupportBiomView::FromCoo(MakeRichRows());
	auto sub = BridgeSubsample(biom, /*subsample_depth*/ 4, /*with_replacement*/ false, /*seed*/ 42);

	auto tree = miint::NewickTree::parse("(((A:0.1,B:0.2)Int1:0.3,(C:0.4,D:0.5)Int2:0.6)Int3:0.7,E:0.8);");
	auto bptree = UnifracBptreeView::FromNewickTree(tree);

	r_vec *result = nullptr;
	ComputeStatus s = faith_pd_inmem(sub.support_biom(), bptree.support_bptree(), &result);
	REQUIRE(s == okay);
	REQUIRE(result != nullptr);
	REQUIRE(result->n_samples == sub.support_biom()->n_samples);
	// All Faith PD values are positive (every sample retains at least one
	// non-zero tip after rarefaction; their root-paths sum to > 0).
	for (unsigned int i = 0; i < result->n_samples; ++i) {
		REQUIRE(result->values[i] > 0.0);
	}
	destroy_results_vec(&result);
}
