#include <catch2/catch_test_macros.hpp>
#include "taxdump_parser.hpp"

using namespace miint;

namespace {
// Locate a parsed node by taxid (order-independent assertions).
const TaxdumpNode *FindNode(const std::vector<TaxdumpNode> &nodes, int64_t taxid) {
	for (const auto &n : nodes) {
		if (n.taxid == taxid) {
			return &n;
		}
	}
	return nullptr;
}

// A minimal but structurally-real nodes.dmp. taxid 1 is self-parented (the root);
// taxid 90 (strain) is a leaf. Field count varies on purpose — only the first
// three fields (taxid, parent, rank) are meaningful, and the parser must tolerate
// lines with fewer than the real 13 fields.
const char *kNodes = "1\t|\t1\t|\tno rank\t|\t\t|\n"
                     "10\t|\t1\t|\tno rank\t|\t\t|\n"
                     "20\t|\t10\t|\tsuperkingdom\t|\t\t|\n"
                     "70\t|\t20\t|\tgenus\t|\t\t|\n"
                     "80\t|\t70\t|\tspecies\t|\t\t|\n"
                     "90\t|\t80\t|\tstrain\t|\t\t|\n";

// names.dmp: taxid 80 carries three name classes with the scientific name written
// LAST, so a parser that naively takes the first name for a taxid would fail.
const char *kNames = "1\t|\troot\t|\t\t|\tscientific name\t|\n"
                     "10\t|\tcellular organisms\t|\t\t|\tscientific name\t|\n"
                     "20\t|\tBacteria\t|\t\t|\tscientific name\t|\n"
                     "70\t|\tEscherichia\t|\t\t|\tscientific name\t|\n"
                     "80\t|\tBacillus coli\t|\t\t|\tsynonym\t|\n"
                     "80\t|\tE. coli\t|\t\t|\tcommon name\t|\n"
                     "80\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n"
                     "90\t|\tEscherichia coli K-12\t|\t\t|\tscientific name\t|\n";
} // namespace

// The .dmp delimiter (fields joined by "\t|\t", line terminated by "\t|") is quirky;
// getting it wrong silently shifts every field, so pin it down directly.
TEST_CASE("TaxdumpParser::SplitFields parses the dmp delimiters", "[taxdump]") {
	SECTION("full nodes.dmp line: leading fields plus a preserved trailing empty field") {
		std::string line = "80\t|\t70\t|\tspecies\t|\t\t|\t0\t|\t0\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|\n";
		auto f = TaxdumpParser::SplitFields(line);
		REQUIRE(f.size() == 13);
		CHECK(f[0] == "80");
		CHECK(f[1] == "70");
		CHECK(f[2] == "species");
		CHECK(f[3].empty());  // embl code (empty)
		CHECK(f[12].empty()); // comments (trailing empty field must survive)
	}
	SECTION("names.dmp line") {
		auto f = TaxdumpParser::SplitFields("80\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n");
		REQUIRE(f.size() == 4);
		CHECK(f[0] == "80");
		CHECK(f[1] == "Escherichia coli");
		CHECK(f[3] == "scientific name");
	}
	SECTION("delnodes.dmp single-field line") {
		auto f = TaxdumpParser::SplitFields("888\t|\n");
		REQUIRE(f.size() == 1);
		CHECK(f[0] == "888");
	}
	SECTION("line without a trailing newline still parses") {
		auto f = TaxdumpParser::SplitFields("999\t|\t80\t|");
		REQUIRE(f.size() == 2);
		CHECK(f[0] == "999");
		CHECK(f[1] == "80");
	}
}

TEST_CASE("TaxdumpParser::ParseNodes builds live tree nodes", "[taxdump]") {
	auto nodes = TaxdumpParser::ParseNodes(kNodes, kNames);

	SECTION("every node in nodes.dmp is emitted") {
		CHECK(nodes.size() == 6);
	}
	SECTION("the self-parented root gets a NULL parent") {
		// NCBI encodes the root (taxid 1) as its own parent; the tree-table
		// convention (matching read_newick) is parent_index = NULL for the root.
		const auto *root = FindNode(nodes, 1);
		REQUIRE(root != nullptr);
		CHECK_FALSE(root->parent_taxid.has_value());
		CHECK(root->name == "root");
		CHECK(root->rank == "no rank");
	}
	SECTION("a non-root node keeps its parent and rank") {
		const auto *n = FindNode(nodes, 80);
		REQUIRE(n != nullptr);
		REQUIRE(n->parent_taxid.has_value());
		CHECK(n->parent_taxid.value() == 70);
		CHECK(n->rank == "species");
	}
	SECTION("the scientific name wins over synonym/common regardless of file order") {
		const auto *n = FindNode(nodes, 80);
		REQUIRE(n != nullptr);
		CHECK(n->name == "Escherichia coli");
	}
	SECTION("is_tip is true only for leaves") {
		// A wrong is_tip breaks any downstream leaf-only logic (e.g. tip counts,
		// diversity over the taxonomy), so verify both a leaf and interior nodes.
		REQUIRE(FindNode(nodes, 90) != nullptr);
		CHECK(FindNode(nodes, 90)->is_tip);       // strain leaf
		CHECK_FALSE(FindNode(nodes, 80)->is_tip); // has child 90
		CHECK_FALSE(FindNode(nodes, 1)->is_tip);  // root has children
	}
}

// merged.dmp remaps retired taxids to current ones; without it, older taxids in
// user data silently fail to join the live tree.
TEST_CASE("TaxdumpParser::ParseMerged reads retired->current pairs", "[taxdump]") {
	auto merged = TaxdumpParser::ParseMerged("999\t|\t80\t|\n1000\t|\t70\t|\n");
	REQUIRE(merged.size() == 2);
	CHECK(merged[0].old_taxid == 999);
	CHECK(merged[0].new_taxid == 80);
	CHECK(merged[1].old_taxid == 1000);
	CHECK(merged[1].new_taxid == 70);
}

TEST_CASE("TaxdumpParser::ParseDeleted reads deleted taxids", "[taxdump]") {
	auto deleted = TaxdumpParser::ParseDeleted("888\t|\n889\t|\n");
	REQUIRE(deleted.size() == 2);
	CHECK(deleted[0] == 888);
	CHECK(deleted[1] == 889);
}
