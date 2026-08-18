#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "unifrac_bptree.hpp"
#include "api.hpp" // faith_pd_inmem, r_vec*, destroy_results_vec
#include "unifrac_support_biom.hpp"

#include "NewickTree.hpp"

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <vector>

using miint::unifrac::CooRow;
using miint::unifrac::UnifracBptreeView;
using miint::unifrac::UnifracSupportBiomView;

// Libssu's support_biom_t is CSR (obs-major) per api.hpp:223 + biom_inmem.cpp:144-156.
// Layout invariants:
//   indptr  length n_obs+1; indptr[i]..indptr[i+1] is obs i's slice
//   indices length nnz;     sample column index within each obs row
//   data    length nnz;     count at (obs, sample)
// Any deviation produces undefined behavior on the first libssu call.

TEST_CASE("UnifracSupportBiomView throws on empty input", "[unifrac][support_biom]") {
	REQUIRE_THROWS_AS(UnifracSupportBiomView::FromCoo({}), std::invalid_argument);
}

TEST_CASE("UnifracSupportBiomView canonical CSR shape", "[unifrac][support_biom]") {
	// 3 samples x 4 features, sparse, deliberately out of input order.
	//   F1: S1=1, S3=4
	//   F2: S2=2, S3=5
	//   F3: S1=3
	//   F4: S3=6
	std::vector<CooRow> rows = {
	    {"S2", "F2", 2.0}, {"S1", "F1", 1.0}, {"S3", "F4", 6.0},
	    {"S3", "F1", 4.0}, {"S1", "F3", 3.0}, {"S3", "F2", 5.0},
	};
	auto view = UnifracSupportBiomView::FromCoo(rows);
	const support_biom_t *biom = view.support_biom();

	REQUIRE(biom->n_samples == 3);
	REQUIRE(biom->n_obs == 4);
	REQUIRE(biom->nnz == 6);

	REQUIRE(std::string(biom->sample_ids[0]) == "S1");
	REQUIRE(std::string(biom->sample_ids[1]) == "S2");
	REQUIRE(std::string(biom->sample_ids[2]) == "S3");

	REQUIRE(std::string(biom->obs_ids[0]) == "F1");
	REQUIRE(std::string(biom->obs_ids[1]) == "F2");
	REQUIRE(std::string(biom->obs_ids[2]) == "F3");
	REQUIRE(std::string(biom->obs_ids[3]) == "F4");

	// indptr is obs-indexed: length n_obs+1, monotone, starts at 0, ends at nnz.
	std::vector<uint32_t> indptr(biom->indptr, biom->indptr + biom->n_obs + 1);
	REQUIRE(indptr == std::vector<uint32_t> {0, 2, 4, 5, 6});

	// indices and data, obs-major, sample-sorted within each obs row.
	// F1 row: S1=1.0, S3=4.0  -> sample idx 0, 2
	// F2 row: S2=2.0, S3=5.0  -> sample idx 1, 2
	// F3 row: S1=3.0           -> sample idx 0
	// F4 row: S3=6.0           -> sample idx 2
	std::vector<uint32_t> indices(biom->indices, biom->indices + biom->nnz);
	std::vector<double> data(biom->data, biom->data + biom->nnz);
	REQUIRE(indices == std::vector<uint32_t> {0, 2, 1, 2, 0, 2});
	REQUIRE(data == std::vector<double> {1.0, 4.0, 2.0, 5.0, 3.0, 6.0});
}

TEST_CASE("UnifracSupportBiomView indices strictly increasing within each obs row", "[unifrac][support_biom]") {
	std::vector<CooRow> rows = {
	    {"S2", "F1", 1.0}, {"S1", "F1", 2.0}, {"S3", "F1", 3.0}, {"S2", "F2", 4.0}, {"S1", "F2", 5.0},
	};
	auto view = UnifracSupportBiomView::FromCoo(rows);
	const support_biom_t *biom = view.support_biom();

	for (int o = 0; o < biom->n_obs; ++o) {
		const uint32_t start = biom->indptr[o];
		const uint32_t end = biom->indptr[o + 1];
		for (uint32_t k = start + 1; k < end; ++k) {
			REQUIRE(biom->indices[k] > biom->indices[k - 1]);
		}
	}
}

TEST_CASE("UnifracSupportBiomView sums duplicate (sample,feature) rows", "[unifrac][support_biom]") {
	// Microbiome convention: duplicates are summed. A silent throw would
	// footgun users who join counts from multiple inputs before computing.
	std::vector<CooRow> rows = {
	    {"S1", "F1", 1.5},
	    {"S1", "F1", 2.5},
	    {"S2", "F2", 0.5},
	    {"S2", "F2", 0.5},
	};
	auto view = UnifracSupportBiomView::FromCoo(rows);
	const support_biom_t *biom = view.support_biom();

	REQUIRE(biom->nnz == 2);
	std::vector<double> data(biom->data, biom->data + biom->nnz);
	REQUIRE(data == std::vector<double> {4.0, 1.0});
}

TEST_CASE("UnifracSupportBiomView drops zero-valued cells", "[unifrac][support_biom]") {
	// Sparse storage invariant: zero-aggregates are dropped entirely.
	std::vector<CooRow> rows = {
	    {"S1", "F1", 1.0},
	    {"S1", "F2", 1.0},
	    {"S1", "F2", -1.0},
	    {"S2", "F1", 2.0},
	};
	auto view = UnifracSupportBiomView::FromCoo(rows);
	const support_biom_t *biom = view.support_biom();

	REQUIRE(biom->nnz == 2);
	// F2 (obs idx 1) should have an empty slice: indptr[1] == indptr[2].
	REQUIRE(biom->indptr[1] == biom->indptr[2]);
}

TEST_CASE("UnifracSupportBiomView is deterministic across input row order", "[unifrac][support_biom]") {
	// Seeded reproducibility downstream depends on byte-identical CSR
	// arrays across input shuffles.
	std::vector<CooRow> base = {
	    {"S1", "F1", 1.0}, {"S1", "F3", 3.0}, {"S2", "F2", 2.0},
	    {"S3", "F1", 4.0}, {"S3", "F2", 5.0}, {"S3", "F4", 6.0},
	};

	auto extract = [](const UnifracSupportBiomView &v) {
		const support_biom_t *b = v.support_biom();
		std::vector<uint32_t> indices(b->indices, b->indices + b->nnz);
		std::vector<uint32_t> indptr(b->indptr, b->indptr + b->n_obs + 1);
		std::vector<double> data(b->data, b->data + b->nnz);
		return std::make_tuple(indices, indptr, data);
	};

	auto view_base = UnifracSupportBiomView::FromCoo(base);
	auto [indices_base, indptr_base, data_base] = extract(view_base);

	auto rotated = base;
	std::rotate(rotated.begin(), rotated.begin() + 3, rotated.end());
	auto [indices_r, indptr_r, data_r] = extract(UnifracSupportBiomView::FromCoo(rotated));

	auto reversed = base;
	std::reverse(reversed.begin(), reversed.end());
	auto [indices_rev, indptr_rev, data_rev] = extract(UnifracSupportBiomView::FromCoo(reversed));

	REQUIRE(indices_base == indices_r);
	REQUIRE(indices_base == indices_rev);
	REQUIRE(indptr_base == indptr_r);
	REQUIRE(indptr_base == indptr_rev);
	REQUIRE(data_base == data_r);
	REQUIRE(data_base == data_rev);
}

TEST_CASE("UnifracSupportBiomView feeds faith_pd_inmem with hand-checked values", "[unifrac][support_biom][libssu]") {
	// LIVE libssu integration: confirms the CSR encoding actually drives
	// libssu correctly. Without this, shape tests can pass while the
	// orientation is silently wrong and produces undefined behavior on
	// every distance/Faith-PD call.
	//
	// Tree: ((A:0.1,B:0.2)Int1:0.3,(C:0.4,D:0.5)Int2:0.6);
	// Faith PD for a sample = sum of branch lengths on the minimal subtree
	// covering the tips present (with non-zero count) in that sample.
	//   S1 has only A -> branch path A + Int1 = 0.1 + 0.3 = 0.4
	//   S2 has A,B    -> A + B + Int1            = 0.1 + 0.2 + 0.3 = 0.6
	//   S3 has C,D    -> C + D + Int2            = 0.4 + 0.5 + 0.6 = 1.5
	//   S4 has A,C    -> A + Int1 + C + Int2     = 0.1 + 0.3 + 0.4 + 0.6 = 1.4
	std::vector<CooRow> rows = {
	    {"S1", "A", 1.0}, {"S2", "A", 1.0}, {"S2", "B", 1.0}, {"S3", "C", 1.0},
	    {"S3", "D", 1.0}, {"S4", "A", 1.0}, {"S4", "C", 1.0},
	};
	auto biom_view = UnifracSupportBiomView::FromCoo(rows);

	auto tree = miint::NewickTree::parse("((A:0.1,B:0.2)Int1:0.3,(C:0.4,D:0.5)Int2:0.6);");
	auto bptree_view = UnifracBptreeView::FromNewickTree(tree);

	r_vec *result = nullptr;
	ComputeStatus status = faith_pd_inmem(biom_view.support_biom(), bptree_view.support_bptree(), &result);
	REQUIRE(status == okay);
	REQUIRE(result != nullptr);
	REQUIRE(result->n_samples == 4);

	// Sample IDs are sorted lexicographically by the adapter: S1, S2, S3, S4.
	REQUIRE(std::string(result->sample_ids[0]) == "S1");
	REQUIRE(std::string(result->sample_ids[1]) == "S2");
	REQUIRE(std::string(result->sample_ids[2]) == "S3");
	REQUIRE(std::string(result->sample_ids[3]) == "S4");

	REQUIRE(result->values[0] == Catch::Approx(0.4));
	REQUIRE(result->values[1] == Catch::Approx(0.6));
	REQUIRE(result->values[2] == Catch::Approx(1.5));
	REQUIRE(result->values[3] == Catch::Approx(1.4));

	destroy_results_vec(&result);
}

TEST_CASE("UnifracSupportBiomView exposes its features as the canonical dictionary", "[unifrac][support_biom]") {
	// progressive_pcoa_from_unifrac shears the tree to feature_ids() before every
	// block, so this accessor has to be the SAME set the CSR was built against —
	// deduplicated and lexicographically sorted, matching obs_ids one for one. If
	// it ever drifted (an id present in the matrix but missing here), the shear
	// would prune a tip some sample in that block actually uses and the block's
	// distances would silently change.
	std::vector<CooRow> rows = {
	    {"S2", "F2", 2.0}, {"S1", "F1", 1.0}, {"S3", "F4", 6.0}, {"S3", "F1", 4.0},
	    {"S1", "F3", 3.0}, {"S3", "F2", 5.0}, {"S1", "F1", 7.0}, // duplicate id, must not duplicate here
	};
	auto view = UnifracSupportBiomView::FromCoo(rows);
	const support_biom_t *biom = view.support_biom();
	const auto &fids = view.feature_ids();

	REQUIRE(fids == std::vector<std::string> {"F1", "F2", "F3", "F4"});
	REQUIRE(fids.size() == static_cast<size_t>(biom->n_obs));
	for (int i = 0; i < biom->n_obs; ++i) {
		REQUIRE(fids[static_cast<size_t>(i)] == std::string(biom->obs_ids[i]));
	}
}

TEST_CASE("UnifracSupportBiomView canonicalizes ids the dictionary never sees in order", "[unifrac][support_biom]") {
	// The dictionary is built from a hash set, whose iteration order is arbitrary
	// and differs from insertion order — the sort is what makes the result
	// canonical. Ids chosen so that first-appearance order, reverse order and
	// lexicographic order all disagree.
	std::vector<CooRow> rows = {
	    {"b", "z", 1.0}, {"c", "y", 2.0}, {"a", "x", 3.0}, {"c", "z", 4.0}, {"a", "y", 5.0},
	};
	auto view = UnifracSupportBiomView::FromCoo(rows);
	const support_biom_t *biom = view.support_biom();

	REQUIRE(std::string(biom->sample_ids[0]) == "a");
	REQUIRE(std::string(biom->sample_ids[1]) == "b");
	REQUIRE(std::string(biom->sample_ids[2]) == "c");
	REQUIRE(view.feature_ids() == std::vector<std::string> {"x", "y", "z"});
}
