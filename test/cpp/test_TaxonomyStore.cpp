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

// Decoration compares names as integers; equal text at one rank must be one NameId, and
// the same text at two ranks two different taxa.
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

// Backfill walks up this map to supply the ranks between a placed name and its nearest
// named ancestor.
// FindRow's index views the store's own id strings. A copy would view the original's, so
// copying is refused; a move keeps the strings where they are, so the index still works.
static_assert(!std::is_copy_constructible_v<TaxonomyStore> && !std::is_copy_assignable_v<TaxonomyStore>);
TEST_CASE("TaxonomyStore finds rows by feature_id after a move", "[taxonomy_store]") {
	TaxonomyStore store({Row("G1", {"A", "B"}), Row("a much longer feature id than any SSO buffer", {"A", "C"})});
	TaxonomyStore moved(std::move(store));
	CHECK(moved.FindRow("G1") == std::optional<size_t>(0));
	CHECK(moved.FindRow("a much longer feature id than any SSO buffer") == std::optional<size_t>(1));
}

TEST_CASE("TaxonomyStore records each name's parent", "[taxonomy_store]") {
	TaxonomyStore store({Row("G1", {"A", "B", "C", "D", "E", "F"}), Row("G2", {"A", "B", "C", "D", "E2", "F2"})});
	auto parent_text = [&](size_t row, size_t rank) {
		return store.NameText(rank - 1, store.Parent(rank, store.Name(row, rank)));
	};
	CHECK(parent_text(0, kGenus) == "E");
	CHECK(parent_text(1, kGenus) == "E2");
	CHECK(parent_text(0, kPhylum) == "A");
	CHECK(store.Name(0, kDomain) != kMissingName);
	CHECK(store.Parent(kDomain, store.Name(0, kDomain)) == kMissingName);
}

// Every converted or hand-built table must hold bare names. Each of these would become a
// name that no other row spells the same way, or break the "; "-joined tree labels.
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
		CHECK_THAT(message, ContainsSubstring("feature_id 'G1'"));
		CHECK_THAT(message, ContainsSubstring("genus"));
		CHECK_THAT(message, ContainsSubstring(c.reason));
	}
}

// A NULL above a name is a hole the caller must fill; tax2tree does not guess what goes there.
TEST_CASE("TaxonomyStore rejects gaps", "[taxonomy_store]") {
	auto message = ErrorOf({Row("G1", {"A", std::nullopt, "C"})});
	CHECK_THAT(message, ContainsSubstring("feature_id 'G1'"));
	CHECK_THAT(message, ContainsSubstring("'phylum' is NULL"));
	// A NULL domain under a named phylum is a gap too.
	CHECK_THAT(ErrorOf({Row("G2", {std::nullopt, "B"})}), ContainsSubstring("'domain' is NULL"));
}

// Lineages are matched to tips by feature_id. A repeated id would count one tip twice;
// an empty one could match an unnamed tip.
TEST_CASE("TaxonomyStore rejects empty and duplicate feature_ids", "[taxonomy_store]") {
	CHECK_THAT(ErrorOf({Row("G1", {"A", "B"}), Row("G2", {"A", "B"}), Row("G1", {"A", "C"})}),
	           ContainsSubstring("feature_id 'G1' appears 2 times"));
	CHECK_THAT(ErrorOf({Row("", {"A", "B"})}), ContainsSubstring("empty feature_id"));
	// A row with no id is named by its position, so its other problems can be found.
	CHECK_THAT(ErrorOf({Row("G1", {"A", "B"}), Row("", {"A", "B", " C"})}),
	           ContainsSubstring("row 2 (empty feature_id), column 'class'"));
}

// The reference nulls only the bad field, leaving the names below it under a missing
// rank: a gap. Truncating at the first bad name never makes one.
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
		TaxonomyStore store({Row("G1", {"A", "B", "C", "Uncultured order"}), Row("G2", {"A", "B", "unnameable"})});
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

// A lineage without a domain and phylum says too little to place anywhere.
TEST_CASE("TaxonomyStore wipes a lineage missing its domain or phylum", "[taxonomy_store]") {
	SECTION("a bad phylum truncates to the domain, which is then wiped (A.3)") {
		TaxonomyStore store({Row("G1", {"A", "unclassified Bacteria", "C"})});
		CHECK_FALSE(store.IsInformative(0));
		CHECK(Text(store, 0, kDomain) == "<missing>");
		CHECK(store.WipedRows() == std::vector<size_t> {0});
		REQUIRE(store.Truncations().size() == 1);
		CHECK(store.Truncations()[0].rank == kPhylum);
	}
	SECTION("domain 'Unclassified' is wiped with no special case (A.12)") {
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

// Backfill and correction look a name's ancestors up in the hierarchy. A name under two
// parents leaves them nothing to look up but a guess.
TEST_CASE("TaxonomyStore rejects a name with more than one parent", "[taxonomy_store]") {
	SECTION("the error names the name, its rank and every parent lineage") {
		auto message =
		    ErrorOf({Row("G1", {"A", "B", "C", "D", "E1", "Foo"}), Row("G2", {"A", "B", "C", "D", "E2", "Foo"})});
		CHECK_THAT(message, ContainsSubstring("'Foo'"));
		CHECK_THAT(message, ContainsSubstring("genus"));
		CHECK_THAT(message, ContainsSubstring("d__A; p__B; c__C; o__D; f__E1"));
		CHECK_THAT(message, ContainsSubstring("d__A; p__B; c__C; o__D; f__E2"));
	}
	SECTION("the usual incertae sedis placeholder, under two orders") {
		auto message = ErrorOf({Row("G1", {"A", "B", "C", "D1", "Incertae Sedis", "G1"}),
		                        Row("G2", {"A", "B", "C", "D2", "Incertae Sedis", "G2"})});
		CHECK_THAT(message, ContainsSubstring("'Incertae Sedis'"));
	}
	SECTION("the hierarchy is built from cleaned lineages") {
		// G2's family is a placeholder, so cleaning removes it and the genus below it.
		TaxonomyStore store(
		    {Row("G1", {"A", "B", "C", "D", "E1", "Foo"}), Row("G2", {"A", "B", "C", "D", "uncultured E2", "Foo"})});
		CHECK(Text(store, 1, kGenus) == "<missing>");
	}
}

// One run shows the shape of a broken taxonomy without an unbounded message.
TEST_CASE("TaxonomyStore reports the problem count and the first 10", "[taxonomy_store]") {
	std::vector<TaxonomyRow> rows;
	for (int i = 1; i <= 12; i++) {
		std::string id = (i < 10 ? "B0" : "B") + std::to_string(i);
		rows.push_back(Row(id, {"A", "B", "C", "D", "E", " leading space"}));
	}
	auto message = ErrorOf(rows);
	CHECK_THAT(message, ContainsSubstring("12 problems"));
	CHECK_THAT(message, ContainsSubstring("feature_id 'B10'"));
	CHECK_THAT(message, !ContainsSubstring("feature_id 'B11'"));
	CHECK(message.find("'B01'") < message.find("'B02'"));
}
