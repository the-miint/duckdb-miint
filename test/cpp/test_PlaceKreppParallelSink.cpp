// place_krepp's consumer runs on every thread.
//
// Under the default preserve_insertion_order, a table function that declares neither an order nor
// a batch index feeds a sink that keeps the order by running on one thread, and the whole pipeline
// with it (docs/internals/duckdb-engine-notes.md, "Output order costs threads and memory").
// place_krepp's rows come out in no fixed order, so its SELECT and INSERT must get parallel sinks.
//
// The SQL tests cannot see this: the rows, and their count, are the same on one thread as on many.
// The test asks the planner instead, which needs a bind and no krepp index.

#include <catch2/catch_test_macros.hpp>

#include "duckdb/execution/operator/helper/physical_result_collector.hpp"
#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/main/client_data.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/prepared_statement_data.hpp"

#include <algorithm>
#include <string>

namespace {

using namespace duckdb;

void Exec(Connection &con, const std::string &sql) {
	auto result = con.Query(sql);
	INFO(sql);
	REQUIRE_FALSE(result->HasError());
}

//! The plan PREPARE built for `name`. Planning binds place_krepp without running InitGlobal, so the
//! index path is never opened.
PreparedStatementData &Prepared(Connection &con, const std::string &name, const std::string &sql) {
	Exec(con, "PREPARE " + name + " AS " + sql);
	auto &statements = ClientData::Get(*con.context).prepared_statements;
	auto entry = std::find_if(statements.begin(), statements.end(),
	                          [&](const auto &statement) { return statement.first == name; });
	REQUIRE(entry != statements.end());
	REQUIRE(entry->second->physical_plan);
	return *entry->second;
}

} // namespace

TEST_CASE("place_krepp feeds a parallel sink under preserve_insertion_order", "[krepp]") {
	DuckDB db(nullptr);
	Connection con(db);
	Exec(con, "SET threads = 4");
	// preserve_insertion_order is on by default; set it so the test does not lean on the default.
	Exec(con, "SET preserve_insertion_order = true");
	Exec(con, "CREATE TABLE reads AS SELECT 'r1' AS read_id, 'ACGT' AS sequence1");
	const std::string scan = "SELECT * FROM place_krepp(query_table := 'reads', index_path := 'no_such_index')";

	SECTION("a SELECT's result sink") {
		auto &prepared = Prepared(con, "q", scan);
		auto collector = PhysicalResultCollector::GetResultCollector(*con.context, prepared);
		REQUIRE(collector->ParallelSink());
	}

	SECTION("an INSERT's sink") {
		// A CTAS is planned the same way, but DuckDB cannot prepare one.
		Exec(con, "CREATE TABLE placements (fragment VARCHAR, edge_num BIGINT, likelihood DOUBLE, like_weight_ratio "
		          "DOUBLE, distal_length DOUBLE, pendant_length DOUBLE, distance DOUBLE)");
		auto &prepared = Prepared(con, "i", "INSERT INTO placements " + scan);
		auto &root = prepared.physical_plan->Root();
		REQUIRE(root.type == PhysicalOperatorType::INSERT);
		REQUIRE(root.ParallelSink());
	}
}
