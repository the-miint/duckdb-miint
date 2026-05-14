#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "NewickTree.hpp"
#include "unifrac_bptree.hpp"
#include "unifrac_support_biom.hpp" // for the live faith_pd test fixture

#include "unifrac_libssu.hpp" // load_bptree_opaque, convert_bptree_opaque, faith_pd_inmem

#include <string>
#include <vector>

using miint::unifrac::CooRow;
using miint::unifrac::UnifracSupportBiomView;

using miint::unifrac::UnifracBptreeView;
using miint::unifrac::ValidateTreeCoversFeatures;

// Libssu reads `support_bptree_t` as the balanced-parens (BP) encoding of a
// rooted tree: pre-order open-paren on entering a node, post-order close-paren
// on exiting. names[k] / lengths[k] are populated at the open-paren index and
// blank/0 at the close. Misencoding silently produces wrong UniFrac distances,
// so the tests below pin the encoding and exercise libssu's own loader to
// confirm the round-trip.

TEST_CASE("UnifracBptreeView encodes 3-tip Newick into 5-node BP", "[unifrac][bptree]") {
	// ((A:0.1,B:0.2):0.3,C:0.4); -> 3 tips + 2 internal = 5 nodes -> 10 parens.
	auto tree = miint::NewickTree::parse("((A:0.1,B:0.2):0.3,C:0.4);");
	auto view = UnifracBptreeView::FromNewickTree(tree);
	const support_bptree_t *bp = view.support_bptree();

	REQUIRE(bp->n_parens == 10);

	// Open at 0, close at 9 (the root).
	REQUIRE(bp->structure[0] == true);
	REQUIRE(bp->structure[bp->n_parens - 1] == false);

	// Parens balance: cumulative open-minus-close stays >= 0 and ends at 0.
	int depth = 0;
	for (int i = 0; i < bp->n_parens; ++i) {
		depth += bp->structure[i] ? 1 : -1;
		REQUIRE(depth >= 0);
	}
	REQUIRE(depth == 0);

	// Each tip name appears exactly once at an open-paren position; close-paren
	// positions have empty names. Internal nodes here are unnamed.
	int count_A = 0, count_B = 0, count_C = 0, count_other = 0;
	for (int i = 0; i < bp->n_parens; ++i) {
		std::string name = bp->names[i];
		if (bp->structure[i]) {
			if (name == "A") {
				++count_A;
			} else if (name == "B") {
				++count_B;
			} else if (name == "C") {
				++count_C;
			} else if (!name.empty()) {
				++count_other;
			}
		} else {
			REQUIRE(name.empty());
		}
	}
	REQUIRE(count_A == 1);
	REQUIRE(count_B == 1);
	REQUIRE(count_C == 1);
	REQUIRE(count_other == 0);
}

TEST_CASE("UnifracBptreeView BP element count matches libssu's Newick parser", "[unifrac][bptree]") {
	// Sanity check on the cardinality of the BP encoding (n_parens == 2*nodes).
	// This catches gross encoding bugs but does NOT prove semantic equivalence
	// — two distinct topologies with the same node count produce identical
	// n_parens. The live faith_pd test below is the actual semantic check.
	const char *newick = "((A:0.1,B:0.2):0.3,C:0.4);";

	auto tree = miint::NewickTree::parse(newick);
	auto view = UnifracBptreeView::FromNewickTree(tree);

	opaque_bptree_t *via_newick = nullptr;
	load_bptree_opaque(newick, &via_newick);
	REQUIRE(via_newick != nullptr);
	int els_via_newick = get_bptree_opaque_els(via_newick);
	destroy_bptree_opaque(&via_newick);

	opaque_bptree_t *via_ours = nullptr;
	convert_bptree_opaque(view.support_bptree(), &via_ours);
	REQUIRE(via_ours != nullptr);
	int els_via_ours = get_bptree_opaque_els(via_ours);
	destroy_bptree_opaque(&via_ours);

	REQUIRE(els_via_ours == els_via_newick);
}

TEST_CASE("UnifracBptreeView drives faith_pd_inmem on an asymmetric tree", "[unifrac][bptree][libssu]") {
	// LIVE libssu test on an asymmetric, deeper-nested tree to complement the
	// symmetric 4-tip fixture in test_UnifracSupportBiom. If BP encoding is
	// wrong on any tree shape, Faith PD diverges from hand-checked values.
	//
	// Tree: (((A:0.1,B:0.2):0.3,C:0.4):0.5,D:0.6);
	//   A and B share parent (depth 3), then merge with C at depth 2, finally
	//   meet D at the root. Asymmetric: A is 3 internal hops from D.
	//
	// Faith PD = sum of branch lengths on the minimal subtree covering tips:
	//   S1 has only A     -> A + (A,B-parent) + ((A,B),C-parent) = 0.1+0.3+0.5 = 0.9
	//   S2 has only D     -> D                                   = 0.6
	//   S3 has A,D        -> A + (A,B-parent) + ((A,B),C-parent) + D = 0.1+0.3+0.5+0.6 = 1.5
	//   S4 has B,C        -> B + (A,B-parent) + C + ((A,B),C-parent) = 0.2+0.3+0.4+0.5 = 1.4
	std::vector<CooRow> rows = {
	    {"S1", "A", 1.0}, {"S2", "D", 1.0}, {"S3", "A", 1.0}, {"S3", "D", 1.0}, {"S4", "B", 1.0}, {"S4", "C", 1.0},
	};
	auto biom_view = UnifracSupportBiomView::FromCoo(rows);

	auto tree = miint::NewickTree::parse("(((A:0.1,B:0.2):0.3,C:0.4):0.5,D:0.6);");
	auto bptree_view = UnifracBptreeView::FromNewickTree(tree);

	r_vec *result = nullptr;
	ComputeStatus status = faith_pd_inmem(biom_view.support_biom(), bptree_view.support_bptree(), &result);
	REQUIRE(status == okay);
	REQUIRE(result != nullptr);
	REQUIRE(result->n_samples == 4);

	REQUIRE(result->values[0] == Catch::Approx(0.9));
	REQUIRE(result->values[1] == Catch::Approx(0.6));
	REQUIRE(result->values[2] == Catch::Approx(1.5));
	REQUIRE(result->values[3] == Catch::Approx(1.4));

	destroy_results_vec(&result);
}

TEST_CASE("UnifracBptreeView preserves quoted/special-char tip names", "[unifrac][bptree]") {
	// NewickTree::parse already supports quoted labels — guard against
	// regression in our BP layer.
	auto tree = miint::NewickTree::parse("(('GG OTU 1':0.1,B:0.2):0.3,C:0.4);");
	auto view = UnifracBptreeView::FromNewickTree(tree);
	const support_bptree_t *bp = view.support_bptree();

	bool found = false;
	for (int i = 0; i < bp->n_parens; ++i) {
		if (std::string(bp->names[i]) == "GG OTU 1") {
			found = true;
			break;
		}
	}
	REQUIRE(found);
}

TEST_CASE("UnifracBptreeView places branch length at open-paren only", "[unifrac][bptree]") {
	// Convention pinned: branch_length[open_pos] = node's branch length,
	// branch_length[close_pos] = 0. Confirms our encoding agrees with the
	// libssu test fixture format.
	auto tree = miint::NewickTree::parse("(A:0.1,B:0.2);");
	auto view = UnifracBptreeView::FromNewickTree(tree);
	const support_bptree_t *bp = view.support_bptree();

	double total_at_open = 0.0;
	double total_at_close = 0.0;
	for (int i = 0; i < bp->n_parens; ++i) {
		if (bp->structure[i]) {
			total_at_open += bp->lengths[i];
		} else {
			total_at_close += bp->lengths[i];
		}
	}
	REQUIRE(total_at_close == 0.0);
	REQUIRE(total_at_open == Catch::Approx(0.3)); // 0.1 + 0.2 (root has no branch length)
}

TEST_CASE("ValidateTreeCoversFeatures passes when features subset tree tips", "[unifrac][bptree]") {
	auto tree = miint::NewickTree::parse("((A:0.1,B:0.2):0.3,C:0.4);");
	REQUIRE_NOTHROW(ValidateTreeCoversFeatures(tree, {"A", "B", "C"}));
	REQUIRE_NOTHROW(ValidateTreeCoversFeatures(tree, {"A"}));
	REQUIRE_NOTHROW(ValidateTreeCoversFeatures(tree, {}));
}

TEST_CASE("ValidateTreeCoversFeatures throws naming the missing feature", "[unifrac][bptree]") {
	// WHY: libssu returns table_and_tree_do_not_overlap with no diagnostic of
	// WHICH feature was missing. Pre-validating in our layer is the only way
	// to surface that to users so they can fix upstream filtering.
	auto tree = miint::NewickTree::parse("((A:0.1,B:0.2):0.3,C:0.4);");
	REQUIRE_THROWS_WITH(ValidateTreeCoversFeatures(tree, {"A", "MISSING_F", "B"}),
	                    Catch::Matchers::ContainsSubstring("MISSING_F"));
}

TEST_CASE("ValidateTreeCoversFeatures does not match internal node names", "[unifrac][bptree]") {
	// Internal nodes can carry labels (bootstrap support, named clades). They
	// are NOT tips and must not satisfy a feature-id lookup. WHY: a feature
	// table maps observation IDs to tips, never to internals.
	auto tree = miint::NewickTree::parse("((A:0.1,B:0.2)Internal:0.3,C:0.4);");
	REQUIRE_NOTHROW(ValidateTreeCoversFeatures(tree, {"A", "B", "C"}));
	REQUIRE_THROWS_WITH(ValidateTreeCoversFeatures(tree, {"Internal"}), Catch::Matchers::ContainsSubstring("Internal"));
}

TEST_CASE("ValidateTreeCoversFeatures ignores extra tips in tree", "[unifrac][bptree]") {
	// Libssu silently ignores tree tips that don't appear in the feature
	// table — so we must mirror that. Only the converse (feature not in tree)
	// is an error.
	auto tree = miint::NewickTree::parse("(((A:0.1,B:0.2):0.3,C:0.4):0.5,D:0.6);");
	REQUIRE_NOTHROW(ValidateTreeCoversFeatures(tree, {"A", "B"}));
}
