// QuerySequenceStream's sub-batches, closed by rows or by bases (src/include/sequence_table_reader.hpp).
//
// A sub-batch is the unit of work a table function's threads claim. place_krepp's threads each take
// one sub-batch at a time, so a call over a few thousand long contigs that fits in one row-sized
// sub-batch runs on one thread while the rest sit idle (#474). Closing a sub-batch once it holds
// enough bases spreads those contigs across threads without making short-read calls pay for tiny
// batches.
//
// What has to hold, and why:
//  - Every row comes out exactly once and in scan order. A row lost or repeated at a batch boundary
//    is a query silently unplaced or placed twice.
//  - A batch closes at whichever limit it reaches first, rows or bases.
//  - The row that reaches the bases limit is in the batch it closes, so a query longer than the limit
//    is a batch of its own rather than an empty batch, which callers read as end-of-stream.
//  - The bases count runs across DuckDB chunk boundaries and is carried into the remainder after a
//    split, since a batch of long rows is routinely assembled from more than one chunk.
//  - Both mates count, so a paired relation's batches hold the bases they say they do.
//  - No bases limit (0) is the row-only behaviour every other caller already has.

#include <catch2/catch_test_macros.hpp>

#include "sequence_table_reader.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <string>
#include <vector>

namespace {

using namespace duckdb;

void Exec(Connection &con, const std::string &sql) {
	auto result = con.Query(sql);
	INFO(sql);
	REQUIRE_FALSE(result->HasError());
}

//! Every sub-batch the stream hands out, as its read_ids, until it returns an empty one.
//! A table function builds the stream inside the caller's query, where the caller's context has an
//! active transaction (the schema check and the TEMP-object inheritance read its catalog), so the
//! test opens one around it.
std::vector<std::vector<std::string>> Batches(Connection &con, const std::string &table, idx_t rows, idx_t bases) {
	std::vector<std::vector<std::string>> out;
	con.BeginTransaction();
	{
		auto schema = ValidateSequenceTableSchema(*con.context, table, /*allow_bigint=*/false);
		QuerySequenceStream stream(*con.context, table, schema, rows, bases);
		while (true) {
			auto batch = stream.FetchSubBatch();
			if (batch.size() == 0) {
				break;
			}
			out.push_back(batch.read_ids);
		}
	}
	con.Commit();
	return out;
}

std::vector<size_t> Sizes(const std::vector<std::vector<std::string>> &batches) {
	std::vector<size_t> sizes;
	for (auto &b : batches) {
		sizes.push_back(b.size());
	}
	return sizes;
}

//! The batches' ids run together, to compare against the relation's own order.
std::vector<std::string> Flatten(const std::vector<std::vector<std::string>> &batches) {
	std::vector<std::string> ids;
	for (auto &b : batches) {
		ids.insert(ids.end(), b.begin(), b.end());
	}
	return ids;
}

std::vector<std::string> Ids(idx_t n) {
	std::vector<std::string> ids;
	for (idx_t i = 0; i < n; i++) {
		ids.push_back("r" + std::to_string(i));
	}
	return ids;
}

//! `n` rows r0..r{n-1}, each `len` bases, in that order.
void MakeReads(Connection &con, const std::string &table, idx_t n, idx_t len) {
	Exec(con, "CREATE TABLE " + table + " AS SELECT 'r' || i AS read_id, repeat('A', " + std::to_string(len) +
	              ") AS sequence1 FROM range(" + std::to_string(n) + ") t(i) ORDER BY i");
}

} // namespace

TEST_CASE("QuerySequenceStream closes a sub-batch by bases", "[sequence_table_reader]") {
	DuckDB db(nullptr);
	Connection con(db);

	SECTION("the row reaching the limit closes the batch it is in") {
		MakeReads(con, "reads", 10, 100);
		auto batches = Batches(con, "reads", STANDARD_VECTOR_SIZE, 250);
		// 100, 200, 300: the third row reaches 250, so it is the batch's last.
		REQUIRE(Sizes(batches) == std::vector<size_t> {3, 3, 3, 1});
		REQUIRE(Flatten(batches) == Ids(10));
	}

	SECTION("a row longer than the limit is a batch of its own") {
		MakeReads(con, "reads", 3, 1000);
		auto batches = Batches(con, "reads", STANDARD_VECTOR_SIZE, 250);
		REQUIRE(Sizes(batches) == std::vector<size_t> {1, 1, 1});
		REQUIRE(Flatten(batches) == Ids(3));
	}

	SECTION("the row limit still applies when it is reached first") {
		MakeReads(con, "reads", 5, 100);
		auto batches = Batches(con, "reads", 2, 1000);
		REQUIRE(Sizes(batches) == std::vector<size_t> {2, 2, 1});
		REQUIRE(Flatten(batches) == Ids(5));
	}

	SECTION("the count runs across chunk boundaries and into the remainder") {
		// 5,000 rows arrive as three chunks (2,048 + 2,048 + 904). 30,000 bases is 3,000 rows, so the
		// first batch spans two chunks and is split mid-chunk; the second is the 2,000 left over, of
		// which 1,096 were already held from that split.
		MakeReads(con, "reads", 5000, 10);
		auto batches = Batches(con, "reads", 1000000, 30000);
		REQUIRE(Sizes(batches) == std::vector<size_t> {3000, 2000});
		REQUIRE(Flatten(batches) == Ids(5000));
	}

	SECTION("both mates count") {
		Exec(con, "CREATE TABLE pairs AS SELECT 'r' || i AS read_id, repeat('A', 50) AS sequence1, "
		          "repeat('C', 50) AS sequence2 FROM range(4) t(i) ORDER BY i");
		auto batches = Batches(con, "pairs", STANDARD_VECTOR_SIZE, 200);
		REQUIRE(Sizes(batches) == std::vector<size_t> {2, 2});
		REQUIRE(Flatten(batches) == Ids(4));
	}

	SECTION("no bases limit is the row-only behaviour") {
		MakeReads(con, "reads", 10, 100000);
		auto batches = Batches(con, "reads", 4, 0);
		REQUIRE(Sizes(batches) == std::vector<size_t> {4, 4, 2});
		REQUIRE(Flatten(batches) == Ids(10));
	}
}
