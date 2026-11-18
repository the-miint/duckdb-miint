// A second, independent check of CladeCountIndex and AssignNames: the same counts computed by
// a recursive CTE that walks parent links over (node, parent) and (tip, rank, name) tables in
// an in-process DuckDB, then groups, and from them rank safety, the contiguity rule and the
// fold as joins and window functions. The C++ finds a subtree by preorder ranges and binary
// searches, and a name's nodes by walking up from its lowest ones, so an indexing off-by-one
// or a missed candidate cannot hide in both.

#include <catch2/catch_test_macros.hpp>

#include "TaxonomyDecorator.hpp"
#include "TaxonomyStore.hpp"
#include "taxonomy_decorator_fixtures.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <vector>

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

// The counts both tests read, as CTEs: every informative tip (one with a tip_name row) paired
// with itself and each of its ancestors (up), the same for every tip, informative or not, a tip
// being a node that is no node's parent (every_up), and from them each node's count of each name
// (counts), of informative tips (informative) and of all tips (tips), and each name's total
// (totals).
const char *kCounts = R"(
WITH RECURSIVE up(tip, node) AS (
    SELECT DISTINCT node, node FROM tip_name
    UNION ALL
    SELECT up.tip, tree.parent FROM up JOIN tree ON tree.node = up.node WHERE tree.parent IS NOT NULL
),
every_up(tip, node) AS (
    SELECT node, node FROM tree WHERE node NOT IN (SELECT parent FROM tree WHERE parent IS NOT NULL)
    UNION ALL
    SELECT every_up.tip, tree.parent FROM every_up JOIN tree ON tree.node = every_up.node
    WHERE tree.parent IS NOT NULL
),
counts AS (
    SELECT up.node, t.rank, t.name, count(*) AS c FROM up JOIN tip_name t ON t.node = up.tip
    GROUP BY up.node, t.rank, t.name
),
informative AS (SELECT node, count(*) AS m FROM up GROUP BY node),
tips AS (SELECT node, count(*) AS tip_count FROM every_up GROUP BY node),
totals AS (SELECT rank, name, count(*) AS total FROM tip_name GROUP BY rank, name)
)";

// After kCounts: rank safety, the contiguity rule and the fold. MIN_COUNT is replaced by the
// case's min_count. Strain is rank 7.
//
// The fold orders by F1 as a double, which is exact here: 2c and total + m are integers of at
// most about 100, division is correctly rounded, so equal fractions give equal doubles, and two
// different fractions with denominators that small differ by at least 1/10^4, far more than a
// rounding error.
const char *kAssign = R"(
-- Postorder, children in ascending index: a node's path of zero-padded indices from the root,
-- ended by '~', sorts after every descendant's path and before a later sibling's.
, paths(node, path) AS (
    SELECT node, lpad(node::VARCHAR, 6, '0') FROM tree WHERE parent IS NULL
    UNION ALL
    SELECT tree.node, paths.path || '/' || lpad(tree.node::VARCHAR, 6, '0')
    FROM paths JOIN tree ON tree.parent = paths.node
),
post AS (SELECT node, row_number() OVER (ORDER BY path || '~') AS post_rank FROM paths),
internal AS (SELECT DISTINCT parent AS node FROM tree WHERE parent IS NOT NULL),
eligible AS (
    SELECT counts.node, counts.rank, counts.name, counts.c, totals.total
    FROM counts JOIN totals USING (rank, name)
    WHERE counts.c >= greatest(CASE WHEN counts.rank = 7 THEN 1 ELSE MIN_COUNT END, (totals.total + 1) // 2)
),
safe AS (
    SELECT node, rank, any_value(name) AS name FROM eligible
    WHERE node IN (SELECT node FROM internal)
    GROUP BY node, rank HAVING count(*) = 1
),
-- A safe rank is picked when every rank from the node's first safe rank to it is safe.
picked AS (
    SELECT node, rank, name FROM (
        SELECT *, min(rank) OVER (PARTITION BY node) AS first_rank,
               row_number() OVER (PARTITION BY node ORDER BY rank) AS k
        FROM safe)
    WHERE rank - first_rank + 1 = k
),
winners AS (
    SELECT picked.rank, picked.name, picked.node, eligible.c, informative.m, eligible.total
    FROM picked
    JOIN eligible USING (node, rank, name)
    JOIN informative USING (node)
    JOIN tips USING (node)
    JOIN post USING (node)
    QUALIFY row_number() OVER (
        PARTITION BY picked.rank, picked.name
        ORDER BY (2 * eligible.c)::DOUBLE / (eligible.total + informative.m) DESC, tips.tip_count ASC,
                 post.post_rank DESC) = 1
)
)";

// The winners, then every name with tips that has none, with its reason as NotPlacedReason
// numbers it.
const char *kWinners = "SELECT rank, name, node, c, m, total FROM winners ORDER BY rank, name";
const char *kNotPlaced = R"(
SELECT t.rank, t.name, t.total,
       CASE WHEN NOT EXISTS (SELECT 1 FROM eligible e WHERE e.rank = t.rank AND e.name = t.name) THEN 0
            WHEN NOT EXISTS (SELECT 1 FROM safe s WHERE s.rank = t.rank AND s.name = t.name) THEN 1
            ELSE 2 END
FROM totals t
WHERE NOT EXISTS (SELECT 1 FROM winners w WHERE w.rank = t.rank AND w.name = t.name)
ORDER BY t.rank, t.name
)";

// kAssign with its one MIN_COUNT replaced.
std::string WithMinCount(std::string sql, uint32_t min_count) {
	const std::string token = "MIN_COUNT";
	const auto at = sql.find(token);
	REQUIRE(at != std::string::npos);
	REQUIRE(sql.find(token, at + 1) == std::string::npos);
	sql.replace(at, token.size(), std::to_string(min_count));
	return sql;
}

void CreateTables(duckdb::Connection &con) {
	Exec(con, "CREATE TABLE tree (node INTEGER, parent INTEGER)");
	Exec(con, "CREATE TABLE tip_name (node INTEGER, rank INTEGER, name INTEGER)");
}

} // namespace

TEST_CASE("CladeCountIndex agrees with recursive-CTE counts", "[taxonomy_decorator]") {
	duckdb::DuckDB db(nullptr);
	duckdb::Connection con(db);
	CreateTables(con);

	std::mt19937 rng(20261119);
	size_t compared = 0;
	for (int trial = 0; trial < 40; trial++) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 50)(rng);
		const auto c = random_taxonomy_case(rng, n);
		TaxonomyStore store(c.rows);
		const auto index = build_index(c, store);
		Load(con, c, store);
		INFO("trial " << trial << ", " << n << " nodes");

		const std::string counts = kCounts;
		std::map<std::tuple<uint32_t, size_t, NameId>, uint32_t> clade;
		ForEachRow(con, counts + "SELECT node, rank, name, c FROM counts", [&](auto col) {
			clade[{col(0), col(1), col(2)}] = col(3);
		});
		std::map<uint32_t, uint32_t> informative_tip_count;
		ForEachRow(con, counts + "SELECT node, m FROM informative",
		           [&](auto col) { informative_tip_count[col(0)] = col(1); });
		std::map<uint32_t, uint32_t> tip_count;
		ForEachRow(con, counts + "SELECT node, tip_count FROM tips", [&](auto col) { tip_count[col(0)] = col(1); });
		std::map<std::pair<size_t, NameId>, uint32_t> total;
		ForEachRow(con, counts + "SELECT rank, name, total FROM totals", [&](auto col) {
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

// Rank safety, the contiguity rule and the fold, computed relationally, must assign every name
// to the same node with the same counts, and report the same names, with the same reasons, as
// assigned nowhere. Every relational result is first checked against the brute force
// (naive_assign), so a wrong query cannot pass for a wrong AssignNames, and that check runs over
// every case even while AssignNames fails.
TEST_CASE("AssignNames agrees with a relational assignment", "[taxonomy_decorator]") {
	duckdb::DuckDB db(nullptr);
	duckdb::Connection con(db);
	CreateTables(con);

	struct Trial {
		TaxonomyCase c;
		uint32_t min_count;
		std::vector<NameScore> winners;
		std::vector<NotPlacedName> unassigned;
	};
	std::vector<Trial> trials;
	std::mt19937 rng(20261123);
	size_t compared = 0, not_placed = 0, won_on_tip_count = 0, won_on_postorder = 0;
	std::array<size_t, 3> reasons {};
	for (int k = 0; k < 40; k++) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 50)(rng);
		const auto c = random_taxonomy_case(rng, n);
		TaxonomyStore store(c.rows);
		const auto naive = naive_counts(c, store);
		Load(con, c, store);
		for (uint32_t min_count = 1; min_count <= 3; min_count++) {
			INFO("trial " << k << ", " << n << " nodes, min_count " << min_count);
			const std::string ctes = std::string(kCounts) + WithMinCount(kAssign, min_count);
			Trial trial {c, min_count, {}, {}};
			ForEachRow(con, ctes + kWinners, [&](auto col) {
				trial.winners.push_back({col(0), col(1), col(2), col(3), col(4), col(5)});
			});
			ForEachRow(con, ctes + kNotPlaced, [&](auto col) {
				trial.unassigned.push_back({col(0), col(1), col(2), static_cast<NotPlacedReason>(col(3))});
			});
			const auto expected = naive_assign(c, store, naive, min_count);
			REQUIRE(trial.winners == expected.scores);
			REQUIRE(trial.unassigned == expected.not_placed);
			compared += trial.winners.size();
			not_placed += trial.unassigned.size();
			won_on_tip_count += expected.won_on_tip_count;
			won_on_postorder += expected.won_on_postorder;
			for (const auto &name : expected.not_placed) {
				reasons[static_cast<size_t>(name.reason)]++;
			}
			trials.push_back(std::move(trial));
		}
	}
	// The cases must have assigned names, left some unassigned for each reason, and needed each
	// tie-break, or agreement proves little about them.
	CHECK(compared > 500);
	CHECK(not_placed > 100);
	CHECK(won_on_tip_count > 0);
	CHECK(won_on_postorder > 0);
	CHECK(reasons[static_cast<size_t>(NotPlacedReason::kNeverEligible)] > 0);
	CHECK(reasons[static_cast<size_t>(NotPlacedReason::kNeverRankSafe)] > 0);
	CHECK(reasons[static_cast<size_t>(NotPlacedReason::kBlockedByContiguity)] > 0);

	for (size_t k = 0; k < trials.size(); k++) {
		const auto &trial = trials[k];
		INFO("trial " << k / 3 << ", min_count " << trial.min_count);
		TaxonomyStore store(trial.c.rows);
		const auto assignment = AssignNames(build_index(trial.c, store), trial.min_count);
		REQUIRE(assignment.Scores() == trial.winners);
		REQUIRE(assignment.NotPlaced() == trial.unassigned);
	}
}
