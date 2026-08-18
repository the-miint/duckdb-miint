#include <catch2/catch_test_macros.hpp>
#include "taxdump_parser.hpp"

#include "taxdump_test_helpers.hpp"

#include <string>

using namespace miint;
using miint_test::ThrowsMentioning;

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

// ParseNodes deliberately keeps only the scientific name, which leaves every other
// name class in names.dmp unreachable. ParseNames is the verbatim view that makes
// them reachable — consumers that need e.g. the genbank common name pivot over it,
// so the contract is "every row, unfiltered, in file order".
TEST_CASE("TaxdumpParser::ParseNames emits every names.dmp row verbatim", "[taxdump]") {
	auto names = TaxdumpParser::ParseNames(kNames);

	SECTION("no row is filtered out and file order is preserved") {
		REQUIRE(names.size() == 8);
		CHECK(names[0].taxid == 1);
		CHECK(names[0].name == "root");
		CHECK(names[0].name_class == "scientific name");
		CHECK(names[7].taxid == 90);
		CHECK(names[7].name == "Escherichia coli K-12");
	}

	SECTION("a taxon's non-scientific classes survive, unlike in ParseNodes") {
		// The point of this function: taxid 80's synonym and common name are
		// dropped by ParseNodes and must be recoverable here. If this ever
		// regresses to scientific-name-only, the downstream common-name lookup
		// silently returns NULL for every taxon.
		std::vector<std::string> classes_for_80;
		for (const auto &n : names) {
			if (n.taxid == 80) {
				classes_for_80.push_back(n.name_class);
			}
		}
		REQUIRE(classes_for_80.size() == 3);
		CHECK(classes_for_80[0] == "synonym");
		CHECK(classes_for_80[1] == "common name");
		CHECK(classes_for_80[2] == "scientific name");
	}

	SECTION("unique_name is carried through when NCBI populates it") {
		// unique_name is the disambiguator NCBI supplies when a name string is
		// shared across taxa; it is empty on most rows but must not be dropped.
		auto disambiguated = TaxdumpParser::ParseNames("32199\t|\tBacteria\t|\tBacteria <bacteria>\t|\tsynonym\t|\n");
		REQUIRE(disambiguated.size() == 1);
		CHECK(disambiguated[0].unique_name == "Bacteria <bacteria>");
		CHECK(disambiguated[0].name == "Bacteria");
	}
}

// NCBI's taxdump_readme.txt documents names.dmp, merged.dmp and delnodes.dmp with a
// fixed width (4, 2 and 1 fields), and we consume every one of those fields. So a row
// of the wrong width means the format moved, and reading it positionally would hand
// back values from the wrong columns -- a `name_class` that is really a `unique name`
// looks entirely well-formed downstream. Fail loudly instead, naming the member and
// line so the report is actionable.
//
// nodes.dmp is deliberately NOT in this group: the readme documents 13 fields for
// taxdump and 18 for new_taxdump, so its width legitimately varies and only the first
// three are consumed. See the ParseNodes tolerance test below.
TEST_CASE("TaxdumpParser rejects fixed-width members whose width moved", "[taxdump]") {
	SECTION("names.dmp must be exactly 4 fields") {
		// Too few: no name class to read at all.
		CHECK(ThrowsMentioning([] { TaxdumpParser::ParseNames("80\t|\tEscherichia coli\t|\n"); },
		                       {"names.dmp", "line 1", "4"}));
		// Too many: the old guard was a lower bound only, so this read positionally
		// and silently mislabelled the name class. This is the case the review found.
		CHECK(ThrowsMentioning(
		    [] { TaxdumpParser::ParseNames("80\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\textra\t|\n"); },
		    {"names.dmp", "line 1", "5"}));
		// The line number must point at the offending row, not the first row.
		CHECK(ThrowsMentioning(
		    [] {
			    TaxdumpParser::ParseNames("1\t|\troot\t|\t\t|\tscientific name\t|\n"
			                              "10\t|\tcellular organisms\t|\t\t|\tscientific name\t|\n"
			                              "20\t|\tBacteria\t|\n");
		    },
		    {"names.dmp", "line 3"}));
	}

	SECTION("merged.dmp must be exactly 2 fields") {
		CHECK(ThrowsMentioning([] { TaxdumpParser::ParseMerged("999\t|\n"); }, {"merged.dmp", "line 1"}));
		CHECK(ThrowsMentioning([] { TaxdumpParser::ParseMerged("999\t|\t80\t|\t70\t|\n"); },
		                       {"merged.dmp", "line 1", "3"}));
	}

	SECTION("delnodes.dmp must be exactly 1 field") {
		// The pre-existing guard here (`f.empty() || f[0].empty()`) could never reject
		// a wrong width, because SplitFields always emplaces at least one element.
		CHECK(
		    ThrowsMentioning([] { TaxdumpParser::ParseDeleted("888\t|\t999\t|\n"); }, {"delnodes.dmp", "line 1", "2"}));
	}

	SECTION("well-formed members of the documented width still parse") {
		CHECK(TaxdumpParser::ParseNames("80\t|\tEscherichia coli\t|\t\t|\tscientific name\t|\n").size() == 1);
		CHECK(TaxdumpParser::ParseMerged("999\t|\t80\t|\n").size() == 1);
		CHECK(TaxdumpParser::ParseDeleted("888\t|\n").size() == 1);
	}
}

// nodes.dmp is the one member whose width is not fixed across dump flavours, so it
// keeps a lower-bound guard. Only taxid/parent/rank are consumed, and appended
// columns must not break reading -- otherwise a new_taxdump-format nodes.dmp (18
// fields) would stop parsing.
TEST_CASE("TaxdumpParser::ParseNodes tolerates varying nodes.dmp widths", "[taxdump]") {
	const char *names = "1\t|\troot\t|\t\t|\tscientific name\t|\n";

	SECTION("the 13-field taxdump and 18-field new_taxdump widths both parse") {
		std::string thirteen = "1\t|\t1\t|\tno rank\t|\t\t|\t0\t|\t0\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|\n";
		std::string eighteen = "1\t|\t1\t|\tno rank\t|\t\t|\t0\t|\t0\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|"
		                       "\t0\t|\t0\t|\t1\t|\t0\t|\t0\t|\n";
		CHECK(TaxdumpParser::ParseNodes(thirteen, names).size() == 1);
		CHECK(TaxdumpParser::ParseNodes(eighteen, names).size() == 1);
	}

	SECTION("fewer than the three consumed fields is still an error, not a silent skip") {
		// A truncated dump must fail loudly rather than yield a short tree that looks
		// complete -- consistent with Gunzip rejecting a truncated gzip stream.
		CHECK(ThrowsMentioning([&] { TaxdumpParser::ParseNodes("1\t|\t1\t|\n", names); }, {"nodes.dmp", "line 1"}));
	}
}

// strtoll reported neither a partial parse nor a range error, so three classes of bad
// input became ordinary-looking BIGINTs. Only the first was detectable downstream (0
// is not a valid taxid); truncation and clamping produced in-range values a consumer
// could not tell from real data.
TEST_CASE("TaxdumpParser rejects taxids that are not integers", "[taxdump]") {
	SECTION("a non-numeric taxid is rejected instead of becoming 0") {
		CHECK(ThrowsMentioning([] { TaxdumpParser::ParseDeleted("notataxid\t|\n"); },
		                       {"delnodes.dmp", "line 1", "notataxid"}));
	}

	SECTION("a partially-numeric taxid is rejected instead of being truncated") {
		// "123abc" previously yielded 123 -- in range, and indistinguishable from a
		// real taxid 123.
		CHECK(ThrowsMentioning([] { TaxdumpParser::ParseDeleted("123abc\t|\n"); }, {"delnodes.dmp", "123abc"}));
	}

	SECTION("an out-of-range taxid is rejected instead of being clamped") {
		// Previously clamped to INT64_MAX, again in range for a BIGINT column.
		CHECK(ThrowsMentioning([] { TaxdumpParser::ParseDeleted("99999999999999999999\t|\n"); },
		                       {"delnodes.dmp", "99999999999999999999"}));
	}

	SECTION("an empty taxid field is rejected") {
		CHECK(ThrowsMentioning([] { TaxdumpParser::ParseMerged("\t|\t80\t|\n"); }, {"merged.dmp"}));
	}

	SECTION("every numeric field is checked, not only the first") {
		// merged.dmp's new_tax_id is just as load-bearing as old_tax_id: a bad value
		// there would silently remap a retired taxid onto the wrong live node.
		CHECK(ThrowsMentioning([] { TaxdumpParser::ParseMerged("999\t|\tnotanumber\t|\n"); },
		                       {"merged.dmp", "notanumber"}));
	}

	SECTION("real taxid magnitudes are unaffected") {
		// The largest live NCBI taxid is ~3.8M, so nothing legitimate is near the
		// int64 boundary; this pins that the new check does not reject real data.
		auto deleted = TaxdumpParser::ParseDeleted("3799730\t|\n1\t|\n");
		REQUIRE(deleted.size() == 2);
		CHECK(deleted[0] == 3799730);
		CHECK(deleted[1] == 1);
	}
}
