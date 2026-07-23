#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <stdexcept>
#include <vector>

#include "cluster_upgma.hpp"

using Catch::Approx;

// A 4-taxon distance matrix (A,B,C,D = indices 0,1,2,3) with a hand-worked UPGMA:
//   A-B=2, C-D=4, all cross pairs =6.
// Merge 1: A,B (dist 2)  -> node 4, height 1
// Merge 2: C,D (dist 4)  -> node 5, height 2   (avg(AB,C)=avg(AB,D)=6 > 4)
// Merge 3: AB,CD (dist 6)-> node 6 (root), height 3
// Branch lengths (height difference to parent): A=B=1, C=D=2, AB=2, CD=1, root=0.
namespace {
const std::vector<double> D4 = {
    0, 2, 6, 6, // A
    2, 0, 6, 6, // B
    6, 6, 0, 4, // C
    6, 6, 4, 0  // D
};
} // namespace

TEST_CASE("upgma 4-taxon hand-computed topology, heights, branch lengths", "[cluster_upgma]") {
	auto t = miint::UpgmaAverageLinkage(D4, 4);
	REQUIRE(t.size() == 7); // 2n-1

	// Tips 0..3 in position == leaf index, height 0.
	for (int i = 0; i < 4; ++i) {
		CHECK(t[i].leaf_index == i);
		CHECK(t[i].height == Approx(0.0));
	}
	// Internal nodes have leaf_index -1.
	CHECK(t[4].leaf_index == -1);
	CHECK(t[5].leaf_index == -1);
	CHECK(t[6].leaf_index == -1);

	// Parent structure: A,B -> 4 ; C,D -> 5 ; 4,5 -> 6 ; 6 = root.
	CHECK(t[0].parent == 4);
	CHECK(t[1].parent == 4);
	CHECK(t[2].parent == 5);
	CHECK(t[3].parent == 5);
	CHECK(t[4].parent == 6);
	CHECK(t[5].parent == 6);
	CHECK(t[6].parent == -1);

	// Heights.
	CHECK(t[4].height == Approx(1.0));
	CHECK(t[5].height == Approx(2.0));
	CHECK(t[6].height == Approx(3.0));

	// Branch lengths = parent height - own height.
	CHECK(t[0].branch_length == Approx(1.0));
	CHECK(t[1].branch_length == Approx(1.0));
	CHECK(t[2].branch_length == Approx(2.0));
	CHECK(t[3].branch_length == Approx(2.0));
	CHECK(t[4].branch_length == Approx(2.0)); // 3 - 1
	CHECK(t[5].branch_length == Approx(1.0)); // 3 - 2
	CHECK(t[6].branch_length == Approx(0.0)); // root
}

TEST_CASE("upgma two samples -> single merge at half the distance", "[cluster_upgma]") {
	const std::vector<double> d2 = {0, 5, 5, 0};
	auto t = miint::UpgmaAverageLinkage(d2, 2);
	REQUIRE(t.size() == 3);
	CHECK(t[0].leaf_index == 0);
	CHECK(t[1].leaf_index == 1);
	CHECK(t[2].leaf_index == -1);
	CHECK(t[0].parent == 2);
	CHECK(t[1].parent == 2);
	CHECK(t[2].parent == -1);
	CHECK(t[2].height == Approx(2.5)); // 5/2
	CHECK(t[0].branch_length == Approx(2.5));
	CHECK(t[1].branch_length == Approx(2.5));
}

TEST_CASE("upgma recovers two well-separated groups as monophyletic subtrees", "[cluster_upgma]") {
	// 4 samples: {0,1} close, {2,3} close, groups far apart. The two group MRCAs
	// must each cover exactly their group (monophyly).
	const std::vector<double> d = {
	    0, 1, 9, 9, //
	    1, 0, 9, 9, //
	    9, 9, 0, 1, //
	    9, 9, 1, 0  //
	};
	auto t = miint::UpgmaAverageLinkage(d, 4);
	REQUIRE(t.size() == 7);
	// 0 and 1 share their immediate parent; 2 and 3 share theirs; those two
	// parents are distinct and both children of the root.
	CHECK(t[0].parent == t[1].parent);
	CHECK(t[2].parent == t[3].parent);
	CHECK(t[0].parent != t[2].parent);
}

TEST_CASE("upgma guards", "[cluster_upgma]") {
	CHECK_THROWS_AS(miint::UpgmaAverageLinkage({0}, 1), std::invalid_argument);       // n < 2
	CHECK_THROWS_AS(miint::UpgmaAverageLinkage({0, 1, 1}, 2), std::invalid_argument); // size != n*n
}
