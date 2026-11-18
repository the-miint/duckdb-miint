// Parity with tax2tree: AssignNames, run on the committed fixture in data/tax2tree/placement/,
// must leave every internal node exactly the names tax2tree's own placement stages left it
// (pick_names, then name_node_score_fold). test/scripts/gen_tax2tree_placement_oracle.py ran
// tax2tree to make the fixture; provenance.txt records the commit and the command. A missing
// or unreadable fixture fails this test; nothing skips it.
//
// Divergence ledger: where miint deliberately differs from tax2tree, and how the fixture keeps
// clear of each difference. A new divergence is reviewed before it is added here.
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
//    fixture's tips are named T<index>.

#include <catch2/catch_test_macros.hpp>

#include "TaxonomyDecorator.hpp"
#include "TaxonomyStore.hpp"
#include "TreeIndex.hpp"
#include "taxonomy_lineage.hpp"

#include <cstdint>
#include <fstream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

using namespace miint;

namespace {

const std::string kDir = "data/tax2tree/placement/";

// A tab-separated file's rows after its header, every field kept, empty ones included.
std::vector<std::vector<std::string>> ReadTsv(const std::string &path) {
	std::ifstream in(path);
	INFO(path);
	REQUIRE(in.is_open());
	std::string line;
	REQUIRE(std::getline(in, line));
	std::vector<std::vector<std::string>> rows;
	while (std::getline(in, line)) {
		std::vector<std::string> fields(1);
		for (char ch : line) {
			if (ch == '\t') {
				fields.emplace_back();
			} else {
				fields.back() += ch;
			}
		}
		rows.push_back(std::move(fields));
	}
	return rows;
}

uint32_t ToIndex(const std::string &field) {
	return static_cast<uint32_t>(std::stoul(field));
}

// "node_index<TAB>rank<TAB>name", the expected file's line for one assigned name.
std::string Line(uint32_t node, size_t rank, const std::string &name) {
	return std::to_string(node) + "\t" + kLineageRanks[rank].column + "\t" + name;
}

} // namespace

TEST_CASE("AssignNames matches tax2tree's placement on the committed fixture", "[taxonomy_decorator]") {
	const auto cases = ReadTsv(kDir + "cases.tsv");
	REQUIRE(cases.size() == 20);
	size_t compared = 0, below_root = 0;
	for (const auto &row : cases) {
		REQUIRE(row.size() == 2);
		const std::string &name = row[0];
		const uint32_t min_count = ToIndex(row[1]);
		INFO(name << ", min_count " << min_count);

		const auto tree_rows = ReadTsv(kDir + name + ".tree.tsv");
		std::vector<uint32_t> parent(tree_rows.size());
		uint32_t root = TreeIndex::NO_PARENT;
		for (uint32_t v = 0; v < tree_rows.size(); v++) {
			const auto &fields = tree_rows[v];
			REQUIRE(fields.size() == 3);
			REQUIRE(ToIndex(fields[0]) == v);
			parent[v] = fields[1].empty() ? TreeIndex::NO_PARENT : ToIndex(fields[1]);
			if (fields[1].empty()) {
				root = v;
			}
		}

		std::vector<TaxonomyRow> rows;
		for (const auto &fields : ReadTsv(kDir + name + ".taxonomy.tsv")) {
			REQUIRE(fields.size() == 1 + kLineageRanks.size());
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
		// Every taxonomy row is a tip's, and cleaning keeps every lineage whole (divergence 3).
		// Otherwise a fixture that drifted from its tree, or a cleaning rule that differs from
		// tax2tree's, would only lose tips, and could still agree on what was left.
		REQUIRE(rows_on_tips == store.RowCount());
		REQUIRE(store.Truncations().empty());
		for (size_t row = 0; row < store.RowCount(); row++) {
			REQUIRE(store.IsInformative(row));
		}
		CladeCountIndex index(TreeIndex::Build(parent, std::vector<double>(parent.size(), 1.0)), store, row_of_node);

		std::vector<std::string> expected;
		for (const auto &fields : ReadTsv(kDir + name + ".expected.tsv")) {
			REQUIRE(fields.size() == 3);
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
