#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "TaxonomyStore.hpp"

#include <stdexcept>
#include <type_traits>
#include <utility>

using namespace miint;
using Catch::Matchers::ContainsSubstring;

namespace {
enum Rank : size_t { kDomain, kPhylum, kClass, kOrder, kFamily, kGenus, kSpecies, kStrain };

// A row from its feature id and its names, shallowest first; unlisted ranks are NULL.
TaxonomyRow Row(std::string id, std::initializer_list<std::optional<std::string>> names) {
	TaxonomyRow row {std::move(id), {}};
	size_t i = 0;
	for (const auto &name : names) {
		row.ranks[i++] = name;
	}
	return row;
}

// The row's cleaned name at `rank`, or "<missing>".
std::string Text(const TaxonomyStore &store, size_t row, size_t rank) {
	NameId name = store.Name(row, rank);
	return name == kMissingName ? "<missing>" : store.NameText(rank, name);
}

// The message a TaxonomyStore over `rows` throws, or "" if it builds.
std::string ErrorOf(const std::vector<TaxonomyRow> &rows,
                    const std::vector<std::string> &tokens = DefaultBadNameTokens()) {
	try {
		TaxonomyStore store(rows, tokens);
	} catch (const std::invalid_argument &e) {
		return e.what();
	}
	return "";
}
} // namespace

// Names are compared as NameIds: equal text at one rank must be one NameId, and the same
// text at two ranks two different taxa.
TEST_CASE("TaxonomyStore interns names per rank", "[taxonomy_store]") {
	TaxonomyStore store({Row("G1", {"A", "B", "C", "D", "E", "F", "F g"}),
	                     Row("G2", {"A", "B", "C", "D", "E", "F", "F h"}), Row("G3", {"A", "X"}),
	                     Row("G4", {"A", "B", "C", "D", "Same", "Same"})});
	REQUIRE(store.RowCount() == 4);
	CHECK(store.FeatureId(1) == "G2");
	CHECK(store.FindRow("G3") == std::optional<size_t>(2));
	CHECK_FALSE(store.FindRow("G9").has_value());

	CHECK(store.Name(0, kGenus) == store.Name(1, kGenus));
	CHECK(store.Name(0, kSpecies) != store.Name(1, kSpecies));
	CHECK(store.NameCount(kDomain) == 1);
	CHECK(store.NameCount(kPhylum) == 2);
	CHECK(store.NameCount(kSpecies) == 2);
	CHECK(store.NameCount(kStrain) == 0);
	CHECK(Text(store, 1, kSpecies) == "F h");

	// Family "Same" and genus "Same" are two taxa, not one name with two parents.
	CHECK(Text(store, 3, kFamily) == "Same");
	CHECK(Text(store, 3, kGenus) == "Same");

	// A lineage ends where its names end; the ranks after are missing, not a problem.
	CHECK(Text(store, 2, kClass) == "<missing>");
	CHECK(store.IsInformative(2));
}

// FindRow's index views the store's own id strings, so the store refuses copies. A move
// must keep those strings where they are, even once the moved-from store is gone.
static_assert(!std::is_copy_constructible_v<TaxonomyStore> && !std::is_copy_assignable_v<TaxonomyStore>);
TEST_CASE("TaxonomyStore finds rows by feature_id after a move", "[taxonomy_store]") {
	const std::string long_id = "a feature id longer than any small-string buffer";
	std::optional<TaxonomyStore> source(std::in_place,
	                                    std::vector<TaxonomyRow> {Row("G1", {"A", "B"}), Row(long_id, {"A", "C"})});
	TaxonomyStore moved(std::move(*source));
	source.reset();
	CHECK(moved.FindRow("G1") == std::optional<size_t>(0));
	CHECK(moved.FindRow(long_id) == std::optional<size_t>(1));
}

// Each name has one parent, so walking Parent() from a name gives its one lineage.
TEST_CASE("TaxonomyStore records each name's parent", "[taxonomy_store]") {
	TaxonomyStore store({Row("G1", {"A", "B", "C", "D", "E", "F"}), Row("G2", {"A", "B", "C", "D", "E2", "F2"})});
	auto parent_text = [&](size_t row, size_t rank) {
		NameId name = store.Name(row, rank);
		REQUIRE(name != kMissingName);
		NameId parent = store.Parent(rank, name);
		REQUIRE(parent != kMissingName);
		return store.NameText(rank - 1, parent);
	};
	CHECK(parent_text(0, kGenus) == "E");
	CHECK(parent_text(1, kGenus) == "E2");
	CHECK(parent_text(0, kPhylum) == "A");
	REQUIRE(store.Name(0, kDomain) != kMissingName);
	CHECK(store.Parent(kDomain, store.Name(0, kDomain)) == kMissingName);
}

// Rank values must be bare names. Each of these would be a name no other row spells the
// same way, or, for ';', one that reads as two ranks once lineages are joined.
TEST_CASE("TaxonomyStore rejects malformed values, naming the id and column", "[taxonomy_store]") {
	struct Case {
		std::string value;
		std::string reason;
	};
	for (const auto &c : std::vector<Case> {{"", "empty"},
	                                        {" Foo", "whitespace"},
	                                        {"Foo ", "whitespace"},
	                                        {"Foo\t", "whitespace"},
	                                        {"g__Foo", "rank prefix"},
	                                        {"Foo;Bar", "';'"}}) {
		INFO("value '" << c.value << "'");
		auto message = ErrorOf({Row("ok", {"A", "B"}), Row("G1", {"A", "B", "C", "D", "E", c.value})});
		CHECK_THAT(message, ContainsSubstring("feature_id 'G1', column 'genus'"));
		CHECK_THAT(message, ContainsSubstring(c.reason));
	}
}

// A NULL above a name is a hole this store will not fill: the caller supplies a gap-free
// taxonomy. (tax2tree keeps a gap below the phylum, and empties a lineage with no phylum.)
TEST_CASE("TaxonomyStore rejects gaps", "[taxonomy_store]") {
	auto message = ErrorOf({Row("G1", {"A", std::nullopt, "C"})});
	CHECK_THAT(message, ContainsSubstring("feature_id 'G1'"));
	CHECK_THAT(message, ContainsSubstring("'phylum' is NULL"));
	// A NULL domain under a named phylum is a gap too.
	CHECK_THAT(ErrorOf({Row("G2", {std::nullopt, "B"})}), ContainsSubstring("'domain' is NULL"));
}

// FindRow looks rows up by feature_id: a repeated id has no one row, and an empty one names
// nothing.
TEST_CASE("TaxonomyStore rejects empty and duplicate feature_ids", "[taxonomy_store]") {
	// One problem per repeated id, however often it repeats.
	auto repeated =
	    ErrorOf({Row("G1", {"A", "B"}), Row("G2", {"A", "B"}), Row("G1", {"A", "C"}), Row("G1", {"A", "D"})});
	CHECK_THAT(repeated, ContainsSubstring("taxonomy: 1 problem:"));
	CHECK_THAT(repeated, ContainsSubstring("feature_id 'G1' appears 3 times"));
	CHECK_THAT(ErrorOf({Row("", {"A", "B"})}), ContainsSubstring("empty feature_id"));
	// A row with no id is named by its position, so its other problems can be found.
	CHECK_THAT(ErrorOf({Row("G1", {"A", "B"}), Row("", {"A", "B", " C"})}),
	           ContainsSubstring("row 2 (empty feature_id), column 'class'"));
}

// tax2tree's own token list, in its order, then incertae sedis in two spellings, which
// tax2tree lacks. How cleaning applies the tokens differs too; see TaxonomyStore.hpp.
TEST_CASE("DefaultBadNameTokens is tax2tree's BAD_NAMES, then incertae sedis", "[taxonomy_store]") {
	CHECK(DefaultBadNameTokens() ==
	      std::vector<std::string> {"environmental sample", "uncultured", "UNNAMEABLE", "unclassified", "unidentified",
	                                "cluster", "isolate", "environmental samples", "incertae sedis", "incertae_sedis"});
}

// An incertae sedis taxon holds members whose placement is uncertain; it is not a clade, so
// it is a placeholder like "uncultured" and is cut. The cost is every name below it: NCBI
// files the genus Sulfobacillus under "Clostridiales Family XVII. Incertae Sedis", and the
// cut drops the genus with the family.
TEST_CASE("TaxonomyStore reads incertae sedis as a placeholder", "[taxonomy_store]") {
	SECTION("a uniquely named one truncates its lineage, the names below it included") {
		TaxonomyStore store({Row("G1", {"A", "B", "C", "D", "Clostridiales Family XVII. Incertae Sedis",
		                                "Sulfobacillus", "Sulfobacillus thermosulfidooxidans"})});
		CHECK(Text(store, 0, kOrder) == "D");
		CHECK(Text(store, 0, kFamily) == "<missing>");
		CHECK(Text(store, 0, kGenus) == "<missing>");
		REQUIRE(store.Truncations().size() == 1);
		CHECK(store.Truncations()[0].token == "incertae sedis");
	}
	SECTION("the underscore spelling is cut too, in any case") {
		TaxonomyStore store({Row("G1", {"A", "B", "C", "D", "D_fam_Incertae_sedis", "Foo"})});
		CHECK(Text(store, 0, kFamily) == "<missing>");
		REQUIRE(store.Truncations().size() == 1);
		CHECK(store.Truncations()[0].token == "incertae_sedis");
	}
	SECTION("one bare placeholder under several parents is cut, not a multi-parent error") {
		// Cleaning runs before the name checks, so the placeholder never reaches them.
		std::vector<TaxonomyRow> rows;
		for (int i = 1; i <= 3; i++) {
			rows.push_back(Row("G" + std::to_string(i), {"A", "B", "C", "D" + std::to_string(i), "Incertae Sedis"}));
		}
		CHECK(ErrorOf(rows).empty());
	}
}

// Emptying only the bad rank would leave the names below it under a missing rank: a gap.
// Truncating at the first bad name never leaves one.
TEST_CASE("TaxonomyStore truncates a lineage at its first bad name", "[taxonomy_store]") {
	SECTION("the bad rank and every rank below it become missing, and the cut is recorded") {
		TaxonomyStore store({Row("G1", {"A", "B", "C", "D", "uncultured family", "Foo", "Foo bar"})});
		CHECK(Text(store, 0, kOrder) == "D");
		CHECK(Text(store, 0, kFamily) == "<missing>");
		CHECK(Text(store, 0, kGenus) == "<missing>");
		CHECK(Text(store, 0, kSpecies) == "<missing>");
		REQUIRE(store.Truncations().size() == 1);
		CHECK(store.Truncations()[0].row == 0);
		CHECK(store.Truncations()[0].rank == kFamily);
		CHECK(store.Truncations()[0].token == "uncultured");
	}
	SECTION("matching ignores ASCII case, in the value and in the token") {
		TaxonomyStore store({Row("G1", {"A", "B", "C", "Uncultured order"}), Row("G2", {"A", "B", "unCULTured"})},
		                    {"UNCULTURED"});
		CHECK(Text(store, 0, kOrder) == "<missing>");
		CHECK(Text(store, 1, kClass) == "<missing>");
	}
	SECTION("only the first bad rank is recorded") {
		TaxonomyStore store({Row("G1", {"A", "B", "C", "D", "uncultured X", "isolate Y"})});
		REQUIRE(store.Truncations().size() == 1);
		CHECK(store.Truncations()[0].rank == kFamily);
	}
	SECTION("an empty token list disables the filter") {
		TaxonomyStore store({Row("G1", {"A", "B", "C", "D", "uncultured family"})}, {});
		CHECK(Text(store, 0, kFamily) == "uncultured family");
		CHECK(store.Truncations().empty());
	}
	SECTION("an empty token is an error: it would match every name") {
		CHECK_THAT(ErrorOf({Row("G1", {"A", "B"})}, {""}), ContainsSubstring("bad-name token is empty"));
	}
	SECTION("the caller can supply its own tokens") {
		TaxonomyStore store({Row("G1", {"A", "B", "C", "Candidatus Foo"})}, {"candidatus"});
		CHECK(Text(store, 0, kOrder) == "<missing>");
		REQUIRE(store.Truncations().size() == 1);
		CHECK(store.Truncations()[0].token == "candidatus");
	}
}

// A lineage missing its domain or its phylum is emptied entirely.
TEST_CASE("TaxonomyStore wipes a lineage missing its domain or phylum", "[taxonomy_store]") {
	SECTION("a bad phylum truncates to the domain, which is then wiped") {
		// tax2tree tests for a phylum before its bad-name filter, so there a bad phylum
		// leaves a lineage with a domain and no phylum.
		TaxonomyStore store({Row("G1", {"A", "unclassified Bacteria", "C"})});
		CHECK_FALSE(store.IsInformative(0));
		CHECK(Text(store, 0, kDomain) == "<missing>");
		CHECK(store.WipedRows() == std::vector<size_t> {0});
		REQUIRE(store.Truncations().size() == 1);
		CHECK(store.Truncations()[0].rank == kPhylum);
	}
	SECTION("domain 'Unclassified' is wiped by the bad-name rule, with no special case") {
		TaxonomyStore store({Row("G1", {"Unclassified", "B"})});
		CHECK_FALSE(store.IsInformative(0));
		CHECK(store.WipedRows() == std::vector<size_t> {0});
	}
	SECTION("a domain alone is wiped") {
		TaxonomyStore store({Row("G1", {"A"})});
		CHECK_FALSE(store.IsInformative(0));
		CHECK(store.WipedRows() == std::vector<size_t> {0});
	}
	SECTION("an all-NULL row is non-informative, but nothing was wiped") {
		TaxonomyStore store({Row("G1", {})});
		CHECK_FALSE(store.IsInformative(0));
		CHECK(store.WipedRows().empty());
	}
	SECTION("a wiped lineage contributes no names") {
		TaxonomyStore store({Row("G1", {"A"}), Row("G2", {"Z", "Y"})});
		CHECK(store.NameCount(kDomain) == 1);
	}
}

// Every name must have one parent, so walking Parent() from it is unambiguous.
TEST_CASE("TaxonomyStore rejects a name with more than one parent", "[taxonomy_store]") {
	SECTION("the error names the name, its rank and its parents, each in a row's lineage") {
		auto message =
		    ErrorOf({Row("G1", {"A", "B", "C", "D", "E1", "Foo"}), Row("G2", {"A", "B", "C", "D", "E2", "Foo"})});
		CHECK_THAT(message, ContainsSubstring("'Foo' (genus) is under 2 different family names"));
		CHECK_THAT(message, ContainsSubstring("[d__A; p__B; c__C; o__D; f__E1]"));
		CHECK_THAT(message, ContainsSubstring("[d__A; p__B; c__C; o__D; f__E2]"));
	}
	SECTION("each lineage shown is one a row actually has") {
		// C is under B and B2, so walking up from E2 would give p__B; G2 says p__B2.
		auto message =
		    ErrorOf({Row("G1", {"A", "B", "C", "D", "E", "F"}), Row("G2", {"A", "B2", "C", "D", "E2", "F"})});
		CHECK_THAT(message, ContainsSubstring("[d__A; p__B2; c__C; o__D; f__E2]"));
	}
	SECTION("a parent seen again is counted once") {
		// E2 repeats, and so does the first parent, E1.
		auto message =
		    ErrorOf({Row("G1", {"A", "B", "C", "D", "E1", "Foo"}), Row("G2", {"A", "B", "C", "D", "E2", "Foo"}),
		             Row("G3", {"A", "B", "C", "D", "E2", "Foo"}), Row("G4", {"A", "B", "C", "D", "E1", "Foo"})});
		CHECK_THAT(message, ContainsSubstring("is under 2 different family names"));
	}
	SECTION("a name under many parents lists the first 3") {
		std::vector<TaxonomyRow> rows;
		for (int i = 1; i <= 5; i++) {
			rows.push_back(Row("G" + std::to_string(i), {"A", "B", "C", "D" + std::to_string(i), "Foo"}));
		}
		auto message = ErrorOf(rows);
		CHECK_THAT(message, ContainsSubstring("'Foo' (family) is under 5 different order names"));
		CHECK_THAT(message, ContainsSubstring("o__D3], and 2 more"));
		CHECK_THAT(message, !ContainsSubstring("o__D4"));
	}
	SECTION("the hierarchy is built from cleaned lineages") {
		// G2's family is a placeholder, so cleaning removes it and the genus below it, and
		// Foo keeps the one parent G1 gives it.
		std::vector<TaxonomyRow> rows {Row("G1", {"A", "B", "C", "D", "E1", "Foo"}),
		                               Row("G2", {"A", "B", "C", "D", "uncultured E2", "Foo"})};
		REQUIRE(ErrorOf(rows).empty());
		TaxonomyStore store(rows);
		CHECK(store.NameCount(kFamily) == 1);
		const NameId foo = store.Name(0, kGenus);
		REQUIRE(foo != kMissingName);
		const NameId parent = store.Parent(kGenus, foo);
		REQUIRE(parent != kMissingName);
		CHECK(store.NameText(kFamily, parent) == "E1");
		CHECK(Text(store, 1, kGenus) == "<missing>");
	}
	SECTION("names are not checked while any row has a problem") {
		// G2's family has a stray space, so it is not E1. Checking parents over it would add
		// a conflict that fixing the value removes.
		auto message =
		    ErrorOf({Row("G1", {"A", "B", "C", "D", "E1", "Foo"}), Row("G2", {"A", "B", "C", "D", " E1", "Foo"})});
		CHECK_THAT(message, ContainsSubstring("taxonomy: 1 problem:"));
		CHECK_THAT(message, !ContainsSubstring("different family names"));
	}
}

// One run shows the shape of a broken taxonomy without an unbounded message.
TEST_CASE("TaxonomyStore reports the problem count and the first 10", "[taxonomy_store]") {
	// Ids descend, so input order and id order disagree.
	std::vector<TaxonomyRow> rows;
	for (int i = 12; i >= 1; i--) {
		std::string id = (i < 10 ? "B0" : "B") + std::to_string(i);
		rows.push_back(Row(id, {"A", "B", "C", "D", "E", " leading space"}));
	}
	auto message = ErrorOf(rows);
	CHECK_THAT(message, ContainsSubstring("taxonomy: 12 problems; the first 10:"));
	CHECK_THAT(message, ContainsSubstring("feature_id 'B03'"));
	CHECK_THAT(message, !ContainsSubstring("feature_id 'B02'"));
	const auto b12 = message.find("'B12'");
	const auto b11 = message.find("'B11'");
	REQUIRE(b12 != std::string::npos);
	REQUIRE(b11 != std::string::npos);
	CHECK(b12 < b11);
}

// --- FindRowViolations and FindMultipleParents --------------------------------------------
//
// The checks behind both TaxonomyStore's errors and taxonomy_violations' rows. A caller that
// keeps the violations (miintax stores them, reports totals and can write its files anyway)
// needs every one of them, each pointing at its input row, rather than the first 10 in a
// message.

namespace {
const char *KindText(TaxonomyViolationKind kind) {
	switch (kind) {
	case TaxonomyViolationKind::InvalidId:
		return "invalid_id";
	case TaxonomyViolationKind::InvalidValue:
		return "invalid_value";
	case TaxonomyViolationKind::DuplicateId:
		return "duplicate_id";
	case TaxonomyViolationKind::Gap:
		return "gap";
	}
	return "?";
}

// Each violation as one line, "row 2 gap phylum", so that a failure shows the whole list.
std::vector<std::string> Describe(const std::vector<RowViolation> &violations) {
	std::vector<std::string> out;
	for (const auto &v : violations) {
		std::string line = "row " + std::to_string(v.row) + " " + KindText(v.kind);
		switch (v.kind) {
		case TaxonomyViolationKind::InvalidValue:
			line += std::string(" ") + kLineageRanks[v.rank].column + ": " + (v.reason ? v.reason : "<no reason>");
			break;
		case TaxonomyViolationKind::DuplicateId:
			line += " x" + std::to_string(v.count);
			break;
		case TaxonomyViolationKind::Gap:
			line += std::string(" ") + kLineageRanks[v.rank].column;
			break;
		case TaxonomyViolationKind::InvalidId:
			break;
		}
		out.push_back(std::move(line));
	}
	return out;
}

// Each label as one line, "family 'fA': 'oA' [0 1]; NULL [2]", its parents in the order given,
// with the text read from the rows as MultipleParents documents.
std::vector<std::string> Describe(const std::vector<TaxonomyRow> &rows, const std::vector<MultipleParents> &found) {
	std::vector<std::string> out;
	for (const auto &label : found) {
		REQUIRE(label.rank >= 1);
		REQUIRE(!label.rows_by_parent.empty());
		REQUIRE(!label.rows_by_parent[0].empty());
		const auto &name = rows[label.rows_by_parent[0][0]].ranks[label.rank];
		REQUIRE(name);
		std::string line = std::string(kLineageRanks[label.rank].column) + " '" + *name + "':";
		for (size_t p = 0; p < label.rows_by_parent.size(); p++) {
			const auto &under = label.rows_by_parent[p];
			REQUIRE(!under.empty());
			const auto &parent = rows[under[0]].ranks[label.rank - 1];
			line += (p == 0 ? " " : "; ") + (parent ? "'" + *parent + "'" : std::string("NULL")) + " [";
			for (size_t i = 0; i < under.size(); i++) {
				line += (i == 0 ? "" : " ") + std::to_string(under[i]);
			}
			line += "]";
		}
		out.push_back(std::move(line));
	}
	return out;
}

using Lines = std::vector<std::string>;
} // namespace

// Each row is checked in full, so a row with several problems shows them all, and they come in
// the order TaxonomyStore's message lists them: input order, and within a row its id, its
// values by rank, then its gap.
TEST_CASE("FindRowViolations reports every row problem in input order", "[taxonomy_store]") {
	std::vector<TaxonomyRow> rows {Row("G1", {"A", "B", "C"}), Row("", {"A", " B"}),
	                               Row("G2", {"A", std::nullopt, "C", std::nullopt, "E"}),
	                               Row("G3", {"A", "B", "C", "", "g__F"})};
	CHECK(Describe(FindRowViolations(rows)) == Lines {"row 1 invalid_id",
	                                                  "row 1 invalid_value phylum: starts or ends "
	                                                  "with whitespace",
	                                                  "row 2 gap phylum", "row 3 invalid_value order: is empty",
	                                                  "row 3 invalid_value family: starts with a rank prefix such "
	                                                  "as 'p__'"});
	// A clean taxonomy has none, so no rows can be told apart from rows not checked.
	CHECK(FindRowViolations({Row("G1", {"A", "B"}), Row("G2", {"A", "C"})}).empty());
	CHECK(FindRowViolations({}).empty());
}

// The reasons are the text taxonomy_violations returns in `detail`; each is a value no other
// row spells the same way, or, for ';', one that reads as two ranks once lineages are joined.
TEST_CASE("FindRowViolations gives the reason a value is not a name", "[taxonomy_store]") {
	struct Case {
		std::string value;
		std::string reason;
	};
	for (const auto &c : std::vector<Case> {{"", "is empty"},
	                                        {" Foo", "starts or ends with whitespace"},
	                                        {"Foo\t", "starts or ends with whitespace"},
	                                        {"g__Foo", "starts with a rank prefix such as 'p__'"},
	                                        {"Foo;Bar", "contains ';'"}}) {
		INFO("value '" << c.value << "'");
		CHECK(Describe(FindRowViolations({Row("G1", {"A", "B", "C", "D", "E", c.value})})) ==
		      Lines {"row 0 invalid_value genus: " + c.reason});
	}
}

// miintax counts how many rows carry a repeated id (its I2), and a caller looking for the rows
// at fault needs each of them, so every occurrence is reported with the total. An empty id is
// reported as such, and is not a repeat: it names no feature.
TEST_CASE("FindRowViolations reports every occurrence of a repeated id", "[taxonomy_store]") {
	std::vector<TaxonomyRow> rows {Row("G1", {"A", "B"}), Row("G2", {"A", "B"}), Row("G1", {"A", "C"}),
	                               Row("G1", {"A", "D"}), Row("", {"A"}),        Row("", {"A"})};
	CHECK(Describe(FindRowViolations(rows)) == Lines {"row 0 duplicate_id x3", "row 2 duplicate_id x3",
	                                                  "row 3 duplicate_id x3", "row 4 invalid_id", "row 5 invalid_id"});
}

// A gap is a NULL rank above a set one; the first NULL rank is the one reported, once per row
// (miintax's I4 names it). NULLs at the end of a lineage are not gaps: a lineage may stop at
// any rank, or be empty.
TEST_CASE("FindRowViolations reports a gap at its first NULL rank", "[taxonomy_store]") {
	CHECK(Describe(FindRowViolations({Row("G1", {"A", std::nullopt, "C"})})) == Lines {"row 0 gap phylum"});
	// d__dA; p__; c__cZ; o__; f__fZ: two NULLs above set ranks, one gap, at the phylum.
	CHECK(Describe(FindRowViolations({Row("L2", {"dA", std::nullopt, "cZ", std::nullopt, "fZ"})})) ==
	      Lines {"row 0 gap phylum"});
	CHECK(Describe(FindRowViolations({Row("L2", {"dA", "pA", std::nullopt, "oZ", "fZ", "gZ", "gZ sZ"})})) ==
	      Lines {"row 0 gap class"});
	// A NULL domain under a set phylum is a gap too.
	CHECK(Describe(FindRowViolations({Row("G2", {std::nullopt, "B"})})) == Lines {"row 0 gap domain"});
	CHECK(FindRowViolations({Row("G1", {"A", "B"}), Row("G2", {}), Row("G3", {"A", "B", "C", "D", "E", "F", "G", "H"})})
	          .empty());
}

// The parent rules miintax's I1 checks: a label is a rank and a whole name, and its parent is
// the value directly above it, NULL included, so an empty family is a parent like any other.
TEST_CASE("FindMultipleParents counts a NULL rank as a parent", "[taxonomy_store]") {
	std::vector<TaxonomyRow> rows {Row("G1", {"dA", "pA", "cA", "oA", "fA", "gA", "gA sA"}),
	                               Row("L1", {"dA", "pA", "cA", "oA", "fA", "gA"}),
	                               Row("L1", {"dA", "pA", "cA", "oA", std::nullopt, "gA", "gA sB"})};
	CHECK(Describe(rows, FindMultipleParents(rows)) == Lines {"genus 'gA': 'fA' [0 1]; NULL [2]"});
}

TEST_CASE("FindMultipleParents checks a species under its genus", "[taxonomy_store]") {
	std::vector<TaxonomyRow> rows {Row("G1", {"dA", "pA", "cA", "oA", "fA", "gA", "gA sA"}),
	                               Row("L1", {"dA", "pA", "cA", "oA", "fA", "gB", "gA sA"})};
	CHECK(Describe(rows, FindMultipleParents(rows)) == Lines {"species 'gA sA': 'gA' [0]; 'gB' [1]"});
}

// Names compare whole: GTDB's suffixes make Bacillus and Bacillus_A two taxa, which may sit in
// different families. The same text at two ranks is two labels too. A domain has no parent.
TEST_CASE("FindMultipleParents compares labels whole and per rank", "[taxonomy_store]") {
	std::vector<TaxonomyRow> rows {Row("G1", {"dA", "pA", "cA", "oA", "f1", "Bacillus"}),
	                               Row("G2", {"dA", "pA", "cA", "oA", "f2", "Bacillus_A"}),
	                               Row("G3", {"dA", "pA", "cA", "oA", "X", "gX"}),
	                               Row("G4", {"dA", "pA", "cA", "oA", "f3", "X"}), Row("G5", {"dB", "pB"})};
	CHECK(FindMultipleParents(rows).empty());
	// A phylum under two domains is a label under two parents.
	rows.push_back(Row("G6", {"dB", "pA"}));
	CHECK(Describe(rows, FindMultipleParents(rows)) == Lines {"phylum 'pA': 'dA' [0 1 2 3]; 'dB' [5]"});
}

// The parent check runs over every row, a repeated id's included: miintax checks one parent per
// label over its primary and secondary files together, and an id may be in both. Every row under
// each parent is listed, so the rows at fault can be found.
TEST_CASE("FindMultipleParents lists every row under each parent", "[taxonomy_store]") {
	std::vector<TaxonomyRow> rows {Row("G1", {"dA", "pA", "cA", "oA", "fA"}), Row("G1", {"dA", "pA", "cA", "oA", "fA"}),
	                               Row("G2", {"dA", "pA", "cA", "oB", "fA"}),
	                               Row("G3", {"dA", "pA", "cA", "oA", "fA"})};
	CHECK(Describe(rows, FindMultipleParents(rows)) == Lines {"family 'fA': 'oA' [0 1 3]; 'oB' [2]"});
}

// The order is the input's, so a run is reproducible: labels in the order their second parent
// first appears (as TaxonomyStore's message lists them), and each label's parents in the order
// they were first seen.
TEST_CASE("FindMultipleParents keeps input order", "[taxonomy_store]") {
	std::vector<TaxonomyRow> rows {
	    Row("R0", {"d", "p", "c", "oA", "fX"}), Row("R1", {"d", "p", "c", "oA", "fY"}),
	    Row("R2", {"d", "p", "c", "oB", "fY"}), Row("R3", {"d", "p", "c", "oB", "fX"}),
	    Row("R4", {"d", "p", "c", "oC", "fY"}),
	};
	CHECK(Describe(rows, FindMultipleParents(rows)) ==
	      Lines {"family 'fY': 'oA' [1]; 'oB' [2]; 'oC' [4]", "family 'fX': 'oA' [0]; 'oB' [3]"});
}

// TaxonomyStore checks parents over its cleaned lineages, which it describes by depth: the
// ranks a cut removed must not count, or a placeholder cleaning drops would still conflict.
TEST_CASE("FindMultipleParents reads only the ranks above each row's depth", "[taxonomy_store]") {
	std::vector<TaxonomyRow> rows {Row("G1", {"dA", "pA", "cA", "oA", "fA"}),
	                               Row("G2", {"dA", "pA", "cA", "oB", "fA"})};
	CHECK(Describe(rows, FindMultipleParents(rows)) == Lines {"family 'fA': 'oA' [0]; 'oB' [1]"});
	const std::vector<size_t> cut_at_family {5, 4};
	CHECK(FindMultipleParents(rows, &cut_at_family).empty());
	const std::vector<size_t> wiped {5, 0};
	CHECK(FindMultipleParents(rows, &wiped).empty());
	const std::vector<size_t> whole {8, 8};
	CHECK(Describe(rows, FindMultipleParents(rows, &whole)) == Lines {"family 'fA': 'oA' [0]; 'oB' [1]"});
}
