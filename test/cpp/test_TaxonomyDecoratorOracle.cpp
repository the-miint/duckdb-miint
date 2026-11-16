// A second, independent check of CladeCountIndex: the same counts computed by a recursive CTE
// that walks parent links over (node, parent) and (tip, rank, name) tables in an in-process
// DuckDB, then groups. The index finds a subtree by preorder ranges and binary searches, so an
// indexing off-by-one cannot hide in both.

#include <catch2/catch_test_macros.hpp>

#include "TaxonomyDecorator.hpp"
#include "TaxonomyStore.hpp"
#include "taxonomy_decorator_fixtures.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <cstdint>
#include <map>
#include <random>
#include <string>
#include <tuple>

using namespace miint;
using namespace miint::test;

namespace {

void Exec(duckdb::Connection &con, const std::string &sql) {
	auto result = con.Query(sql);
	INFO(sql << ": " << result->GetError());
	REQUIRE_FALSE(result->HasError());
}

// Each row of a query's result, as its columns read as uint32.
template <class F>
void ForEachRow(duckdb::Connection &con, const std::string &sql, F f) {
	auto result = con.Query(sql);
	INFO(sql << ": " << result->GetError());
	REQUIRE_FALSE(result->HasError());
	while (auto chunk = result->Fetch()) {
		for (duckdb::idx_t row = 0; row < chunk->size(); row++) {
			f([&](duckdb::idx_t col) { return chunk->GetValue(col, row).GetValue<uint32_t>(); });
		}
	}
}

// Replaces the tables' contents with one case: every node and its parent, and one row per
// named rank of every informative tip.
void Load(duckdb::Connection &con, const TaxonomyCase &c, const TaxonomyStore &store) {
	Exec(con, "DELETE FROM tree");
	Exec(con, "DELETE FROM tip_name");
	std::string nodes = "INSERT INTO tree VALUES ";
	std::string names;
	for (uint32_t v = 0; v < c.tree.parent.size(); v++) {
		const uint32_t parent = c.tree.parent[v];
		nodes += (v ? ", (" : "(") + std::to_string(v) + ", " +
		         (parent == kNoParent ? std::string("NULL") : std::to_string(parent)) + ")";
		const auto &row = c.row_of_node[v];
		if (!row || !store.IsInformative(*row)) {
			continue;
		}
		for (size_t r = 0; r < kRankCount; r++) {
			const NameId name = store.Name(*row, r);
			if (name != kMissingName) {
				names += (names.empty() ? "(" : ", (") + std::to_string(v) + ", " + std::to_string(r) + ", " +
				         std::to_string(name) + ")";
			}
		}
	}
	Exec(con, nodes);
	if (!names.empty()) {
		Exec(con, "INSERT INTO tip_name VALUES " + names);
	}
}

// Every informative tip (one with a tip_name row) paired with itself and each of its
// ancestors.
const char *kUp = R"(
WITH RECURSIVE up(tip, node) AS (
    SELECT DISTINCT node, node FROM tip_name
    UNION ALL
    SELECT up.tip, tree.parent FROM up JOIN tree ON tree.node = up.node WHERE tree.parent IS NOT NULL
)
)";

// The same for every tip, informative or not: a node that is no node's parent.
const char *kUpFromEveryTip = R"(
WITH RECURSIVE up(tip, node) AS (
    SELECT node, node FROM tree WHERE node NOT IN (SELECT parent FROM tree WHERE parent IS NOT NULL)
    UNION ALL
    SELECT up.tip, tree.parent FROM up JOIN tree ON tree.node = up.node WHERE tree.parent IS NOT NULL
)
)";

} // namespace

TEST_CASE("CladeCountIndex agrees with recursive-CTE counts", "[taxonomy_decorator]") {
	duckdb::DuckDB db(nullptr);
	duckdb::Connection con(db);
	Exec(con, "CREATE TABLE tree (node INTEGER, parent INTEGER)");
	Exec(con, "CREATE TABLE tip_name (node INTEGER, rank INTEGER, name INTEGER)");

	std::mt19937 rng(20261119);
	size_t compared = 0;
	for (int trial = 0; trial < 40; trial++) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 50)(rng);
		const auto c = random_taxonomy_case(rng, n);
		TaxonomyStore store(c.rows);
		const auto index = build_index(c, store);
		Load(con, c, store);
		INFO("trial " << trial << ", " << n << " nodes");

		std::map<std::tuple<uint32_t, size_t, NameId>, uint32_t> clade;
		ForEachRow(con, std::string(kUp) + R"(
            SELECT up.node, t.rank, t.name, count(*) FROM up JOIN tip_name t ON t.node = up.tip
            GROUP BY up.node, t.rank, t.name)",
		           [&](auto col) {
			           clade[{col(0), col(1), col(2)}] = col(3);
		           });
		std::map<uint32_t, uint32_t> informative_tip_count;
		ForEachRow(con, std::string(kUp) + "SELECT node, count(*) FROM up GROUP BY node",
		           [&](auto col) { informative_tip_count[col(0)] = col(1); });
		std::map<uint32_t, uint32_t> tip_count;
		ForEachRow(con, std::string(kUpFromEveryTip) + "SELECT node, count(*) FROM up GROUP BY node",
		           [&](auto col) { tip_count[col(0)] = col(1); });
		std::map<std::pair<size_t, NameId>, uint32_t> total;
		ForEachRow(con, "SELECT rank, name, count(*) FROM tip_name GROUP BY rank, name", [&](auto col) {
			total[{col(0), col(1)}] = col(2);
		});
		compared += clade.size();

		// Every (node, rank, name): the CTE's count, or 0 where it has no row.
		for (uint32_t v = 0; v < n; v++) {
			INFO("node " << v);
			REQUIRE(index.InformativeTipCount(v) == (informative_tip_count.count(v) ? informative_tip_count[v] : 0));
			// Every node has a tip below it, so every node has a row.
			REQUIRE(index.TipCount(v) == tip_count.at(v));
			for (size_t r = 0; r < kRankCount; r++) {
				for (NameId name = 0; name < store.NameCount(r); name++) {
					auto found = clade.find({v, r, name});
					REQUIRE(index.CladeCount(v, r, name) == (found == clade.end() ? 0 : found->second));
				}
			}
		}
		for (size_t r = 0; r < kRankCount; r++) {
			for (NameId name = 0; name < store.NameCount(r); name++) {
				auto found = total.find({r, name});
				REQUIRE(index.Total(r, name) == (found == total.end() ? 0 : found->second));
			}
		}
	}
	// The cases must have produced counts to compare, or agreement proves nothing.
	CHECK(compared > 1000);
}
