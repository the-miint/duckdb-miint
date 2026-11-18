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
//  - A batch holds at most the row limit and at most the bases limit, whichever binds first. A row
//    that would take it past the bases limit starts the next batch (predict_genes' kBatchBases rule).
//  - A row longer than the bases limit is a batch of its own, not an empty batch, which callers read
//    as end-of-stream.
//  - The bases count runs across DuckDB chunk boundaries and on into the rows left after a cut, since
//    a batch of long rows is routinely assembled from more than one chunk.
//  - Only sequence1 counts: place_krepp, the caller, places sequence1 alone, so counting a mate would
//    close its batches at half the bases it actually places.
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

	SECTION("a row that would pass the limit starts the next batch") {
		MakeReads(con, "reads", 10, 100);
		auto batches = Batches(con, "reads", STANDARD_VECTOR_SIZE, 250);
		// 100, 200 fit; a third row would make 300.
		REQUIRE(Sizes(batches) == std::vector<size_t> {2, 2, 2, 2, 2});
		REQUIRE(Flatten(batches) == Ids(10));
	}

	SECTION("a batch may fill the limit exactly") {
		MakeReads(con, "reads", 10, 100);
		auto batches = Batches(con, "reads", STANDARD_VECTOR_SIZE, 300);
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

	SECTION("the count runs across chunk boundaries and into the rows left after a cut") {
		// No chunk holds more than STANDARD_VECTOR_SIZE rows, so a first batch of one and a half
		// chunks' worth must be assembled from at least two chunks and cut partway through one; the
		// second batch then starts with the rows left over from that cut.
		const idx_t v = STANDARD_VECTOR_SIZE;
		const idx_t first = v + v / 2;
		MakeReads(con, "reads", first + v, 10);
		auto batches = Batches(con, "reads", 1000000000, first * 10);
		REQUIRE(Sizes(batches) == std::vector<size_t> {first, v});
		REQUIRE(Flatten(batches) == Ids(first + v));
	}

	SECTION("only sequence1 counts") {
		// 50 bases of sequence1 per row: two rows fill 100. Counting the 50-base mate as well would
		// put one row in each batch.
		Exec(con, "CREATE TABLE pairs AS SELECT 'r' || i AS read_id, repeat('A', 50) AS sequence1, "
		          "repeat('C', 50) AS sequence2 FROM range(4) t(i) ORDER BY i");
		auto batches = Batches(con, "pairs", STANDARD_VECTOR_SIZE, 100);
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
