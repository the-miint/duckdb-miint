// Parity with tax2tree: AssignNames, run on the committed fixture in data/tax2tree/assign_names/,
// must leave every internal node exactly the names tax2tree's own pick_names and
// name_node_score_fold left it. test/scripts/gen_tax2tree_assign_names_oracle.py ran tax2tree to
// make the fixture; provenance.txt records the commit and the command. A missing or unreadable
// fixture fails this test; nothing skips it.
//
// Where AssignNames differs from tax2tree, and how the fixture keeps clear of each difference:
//
// 1. F1 comparison. AssignNames compares F1 exactly, in integers (F1Score::Compare), so equal F1
//    from different counts is a tie; tax2tree compares doubles, which can order such a pair by
//    rounding. The generator rejects, and draws again, any case where a name's candidates include
//    such a pair.
// 2. Strain's min_count. AssignNames uses 1 at strain (MinCountAt); tax2tree uses one min_count at
//    every rank. The fixture has seven ranks, so no strain.
// 3. Cleaning. AssignNames reads lineages as TaxonomyStore cleaned them, with its bad-name tokens
//    and rules; tax2tree cleans with its own. No fixture name holds a bad-name token, and every
//    lineage has a domain and a phylum, so both keep every lineage whole.
// 4. Tip names. miint matches tip names as given; tax2tree strips every apostrophe from them. The
//    fixture's tips are named T<index>, or by genome ids in the golden case.
//
// Decorate and EmitLineages run on a second fixture, data/tax2tree/decorate/, made by the same
// script from tax2tree's next stages (set_preliminary_name_and_rank, backfill_names_gap,
// commonname_promotion, pull_consensus_strings), and must give every internal node tax2tree's
// label and every tip tax2tree's lineage. Differences 1 to 4 apply there too, and five more; the
// generator leaves out each case that hits one, and provenance.txt says which and why:
//
// 5. A root named below the domain. Decorate fills the root's label up to the domain, as for any
//    node with no named ancestor; tax2tree keeps the root's deepest name alone.
// 6. A node below the root with no named ancestor. Decorate fills its label up to the domain;
//    tax2tree keeps its deepest name alone, losing the shallower names the fold gave it.
// 7. Inverted ranks: a node named at a shallower rank than its nearest named ancestor. Decorate
//    keeps its name and records the inversion; tax2tree empties its label.
// 8. Promotion of a name no deeper than the node's names. Decorate stops the node's promotion;
//    tax2tree appends the name, leaving the label out of rank order or with two names at one
//    rank. Only an inverted-rank node leads to it.
// 9. Polyphyletic merges. Decorate merges only a name and its own variant (Bacillus with
//    Bacillus_A); tax2tree merges any name with one that contains it and has a capital-letter
//    suffix somewhere (Baci with Bacillus_A).

#include <catch2/catch_test_macros.hpp>

#include "TaxonomyDecorator.hpp"
#include "TaxonomyStore.hpp"
#include "TreeIndex.hpp"
#include "taxonomy_decorator_fixtures.hpp"
#include "taxonomy_lineage.hpp"

#include <zlib.h>

#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

using namespace miint;

namespace {

const std::string kDir = "data/tax2tree/assign_names/";
const std::string kDecorateDir = "data/tax2tree/decorate/";

using Rows = std::vector<std::vector<std::string>>;

// A gzip-compressed tab-separated file's rows after its header, every field kept, empty ones
// included. zlib reports a file cut short only when it is closed, so the close is checked too.
Rows ReadTsvGz(const std::string &path) {
	INFO(path);
	gzFile in = gzopen(path.c_str(), "rb");
	REQUIRE(in != nullptr);
	std::string text;
	char buffer[1 << 14];
	int n;
	while ((n = gzread(in, buffer, sizeof(buffer))) > 0) {
		text.append(buffer, static_cast<size_t>(n));
	}
	const int closed = gzclose(in);
	REQUIRE(n == 0);
	REQUIRE(closed == Z_OK);
	REQUIRE(!text.empty());
	REQUIRE(text.back() == '\n');

	Rows rows;
	std::vector<std::string> fields(1);
	for (char ch : text) {
		if (ch == '\t') {
			fields.emplace_back();
		} else if (ch == '\n') {
			rows.push_back(std::move(fields));
			fields.assign(1, std::string());
		} else {
			fields.back() += ch;
		}
	}
	rows.erase(rows.begin());
	return rows;
}

// A table's rows by their first field, the case, with that field dropped; each case's rows in
// file order.
std::map<std::string, Rows> ByCase(const Rows &rows, size_t width) {
	std::map<std::string, Rows> out;
	for (const auto &fields : rows) {
		REQUIRE(fields.size() == width);
		out[fields[0]].emplace_back(fields.begin() + 1, fields.end());
	}
	return out;
}

uint32_t ToIndex(const std::string &field) {
	return static_cast<uint32_t>(std::stoul(field));
}

// "node_index<TAB>rank<TAB>name", the expected file's line for one assigned name.
std::string Line(uint32_t node, size_t rank, const std::string &name) {
	return std::to_string(node) + "\t" + kLineageRanks[rank].column + "\t" + name;
}

// A case's taxonomy store, and the counting index over its tree with each named tip's row.
struct FixtureCase {
	TaxonomyStore store;
	CladeCountIndex index;
};

FixtureCase LoadCase(const Rows &tree_rows, const Rows &taxonomy_rows) {
	std::vector<uint32_t> parent(tree_rows.size());
	for (uint32_t v = 0; v < tree_rows.size(); v++) {
		const auto &fields = tree_rows[v];
		REQUIRE(ToIndex(fields[0]) == v);
		parent[v] = fields[1].empty() ? TreeIndex::NO_PARENT : ToIndex(fields[1]);
	}

	std::vector<TaxonomyRow> rows;
	for (const auto &fields : taxonomy_rows) {
		TaxonomyRow taxonomy_row {fields[0], {}};
		for (size_t r = 0; r < kLineageRanks.size(); r++) {
			if (!fields[1 + r].empty()) {
				taxonomy_row.ranks[r] = fields[1 + r];
			}
		}
		rows.push_back(std::move(taxonomy_row));
	}
	TaxonomyStore store(rows);
	std::vector<std::optional<size_t>> row_of_node(tree_rows.size());
	size_t rows_on_tips = 0;
	for (uint32_t v = 0; v < tree_rows.size(); v++) {
		if (!tree_rows[v][2].empty()) {
			row_of_node[v] = store.FindRow(tree_rows[v][2]);
			rows_on_tips += row_of_node[v].has_value();
		}
	}
	// Every taxonomy row is a tip's, and cleaning keeps every lineage whole (difference 3).
	// Otherwise a fixture that drifted from its tree, or a cleaning rule that differs from
	// tax2tree's, would only lose tips, and could still agree on what was left.
	REQUIRE(rows_on_tips == store.RowCount());
	REQUIRE(store.Truncations().empty());
	for (size_t row = 0; row < store.RowCount(); row++) {
		REQUIRE(store.IsInformative(row));
	}
	CladeCountIndex index(TreeIndex::Build(parent, std::vector<double>(parent.size(), 1.0)), store, row_of_node);
	return {std::move(store), std::move(index)};
}

// The name the label of `holder` gives at `rank`, or "" for no holder.
} // namespace

TEST_CASE("AssignNames matches tax2tree on the committed fixture", "[taxonomy_decorator]") {
	const auto cases = ReadTsvGz(kDir + "cases.tsv.gz");
	REQUIRE(cases.size() == 20);
	const auto trees = ByCase(ReadTsvGz(kDir + "trees.tsv.gz"), 4);
	const auto taxonomies = ByCase(ReadTsvGz(kDir + "taxonomy.tsv.gz"), 2 + kLineageRanks.size());
	const auto expectations = ByCase(ReadTsvGz(kDir + "expected.tsv.gz"), 4);
	// Every table holds rows for exactly the listed cases.
	for (const auto *table : {&trees, &taxonomies, &expectations}) {
		REQUIRE(table->size() == cases.size());
	}
	size_t compared = 0, below_root = 0;
	for (const auto &row : cases) {
		REQUIRE(row.size() == 2);
		const std::string &name = row[0];
		const uint32_t min_count = ToIndex(row[1]);
		INFO(name << ", min_count " << min_count);
		REQUIRE(trees.count(name) == 1);
		REQUIRE(taxonomies.count(name) == 1);
		REQUIRE(expectations.count(name) == 1);

		const auto &tree_rows = trees.at(name);
		const auto fixture = LoadCase(tree_rows, taxonomies.at(name));
		const TaxonomyStore &store = fixture.store;
		const CladeCountIndex &index = fixture.index;
		const uint32_t root = index.Tree().root();

		std::vector<std::string> expected;
		for (const auto &fields : expectations.at(name)) {
			expected.push_back(fields[0] + "\t" + fields[1] + "\t" + fields[2]);
			below_root += ToIndex(fields[0]) != root;
		}
		// Both lists by node, then rank, as the generator writes them and NamesAt orders them.
		const auto assignment = AssignNames(index, min_count);
		std::vector<std::string> observed;
		for (uint32_t v = 0; v < tree_rows.size(); v++) {
			for (const auto &score : assignment.NamesAt(v)) {
				observed.push_back(Line(v, score.rank, store.NameText(score.rank, score.name)));
			}
		}
		CHECK(observed == expected);
		compared += expected.size();
	}
	// The fixture must assign names, and below the root, or agreement shows little.
	CHECK(compared > 100);
	CHECK(below_root > 50);
}

TEST_CASE("Decorate and EmitLineages match tax2tree on the committed fixture", "[taxonomy_decorator]") {
	const auto cases = ReadTsvGz(kDecorateDir + "cases.tsv.gz");
	REQUIRE(cases.size() == 21);
	const auto trees = ByCase(ReadTsvGz(kDecorateDir + "trees.tsv.gz"), 4);
	const auto taxonomies = ByCase(ReadTsvGz(kDecorateDir + "taxonomy.tsv.gz"), 2 + kLineageRanks.size());
	const auto labels = ByCase(ReadTsvGz(kDecorateDir + "expected_labels.tsv.gz"), 3);
	const auto lineages = ByCase(ReadTsvGz(kDecorateDir + "expected_lineages.tsv.gz"), 2 + kLineageRanks.size());
	for (const auto *table : {&trees, &taxonomies, &labels, &lineages}) {
		REQUIRE(table->size() == cases.size());
	}
	size_t labels_compared = 0, multi_rank_labels = 0, lineages_compared = 0;
	bool golden = false;
	for (const auto &row : cases) {
		REQUIRE(row.size() == 2);
		const std::string &name = row[0];
		const uint32_t min_count = ToIndex(row[1]);
		INFO(name << ", min_count " << min_count);
		golden = golden || name == "golden";
		REQUIRE(trees.count(name) == 1);
		REQUIRE(taxonomies.count(name) == 1);
		REQUIRE(labels.count(name) == 1);
		REQUIRE(lineages.count(name) == 1);

		const auto &tree_rows = trees.at(name);
		const uint32_t n = static_cast<uint32_t>(tree_rows.size());
		const auto fixture = LoadCase(tree_rows, taxonomies.at(name));
		const TaxonomyStore &store = fixture.store;
		const CladeCountIndex &index = fixture.index;

		const auto result = Decorate(index, store, min_count);
		REQUIRE(result.labels.size() == n);
		// No case hits an inversion (difference 7).
		CHECK(result.inverted.empty());
		std::vector<std::string> observed_labels;
		std::vector<std::vector<LabelPart>> parsed(n);
		for (uint32_t v = 0; v < n; v++) {
			const std::string text = FormatLabel(result.labels[v], store);
			if (!text.empty()) {
				observed_labels.push_back(std::to_string(v) + "\t" + text);
			}
			parsed[v] = ParseNodeLabel(text);
		}
		std::vector<std::string> expected_labels;
		for (const auto &fields : labels.at(name)) {
			expected_labels.push_back(fields[0] + "\t" + fields[1]);
			multi_rank_labels += fields[1].find(';') != std::string::npos;
		}
		CHECK(observed_labels == expected_labels);
		labels_compared += expected_labels.size();

		// Every tip, by node, as the generator writes them and EmitLineages lists them.
		const auto emitted = EmitLineages(index.Tree(), parsed);
		REQUIRE(emitted.tips.size() == emitted.holders.size());
		CHECK(emitted.repeats.empty());
		std::vector<std::string> observed_lineages;
		for (size_t i = 0; i < emitted.tips.size(); i++) {
			std::string line = tree_rows[emitted.tips[i]][2];
			for (size_t r = 0; r < kLineageRanks.size(); r++) {
				line += "\t" + test::HolderName(parsed, emitted.holders[i][r], r);
			}
			observed_lineages.push_back(line);
		}
		std::vector<std::string> expected_lineages;
		for (const auto &fields : lineages.at(name)) {
			std::string line = fields[0];
			for (size_t r = 1; r < fields.size(); r++) {
				line += "\t" + fields[r];
			}
			expected_lineages.push_back(line);
		}
		CHECK(observed_lineages == expected_lineages);
		lineages_compared += expected_lineages.size();
	}
	// The fixture must label nodes, many with several ranks, or agreement shows little.
	CHECK(golden);
	CHECK(labels_compared > 100);
	CHECK(multi_rank_labels > 30);
	CHECK(lineages_compared > 500);
}
