// QuerySequenceStream's sub-batches, closed by rows or by bases (src/include/sequence_table_reader.hpp).
//
// A sub-batch is the unit of work a table function's threads claim, one at a time. A few thousand long
// contigs fit in one row-sized sub-batch, so a bases limit is what spreads them across threads.
//
// What has to hold, and why:
//  - Every row comes out exactly once, in scan order, with all of its columns. A row lost, repeated
//    or misaligned at a batch boundary is a query silently unplaced, placed twice, or placed with
//    another query's sequence.
//  - A batch holds at most the row limit and at most the bases limit, whichever binds first. A row
//    that would take it past the bases limit starts the next batch.
//  - A row longer than the bases limit is a batch of its own, not an empty batch, which callers read
//    as end-of-stream.
//  - Batches are cut correctly however the rows arrived: across DuckDB chunk boundaries, and after
//    the rows already handed out are dropped to make room for the next chunk.
//  - Only sequence1 counts toward the bases limit.
//  - No bases limit (0) is the row-only behaviour, which every caller without one relies on.
//  - Once a chunk fails to convert, the stream raises on every later call rather than handing out
//    rows from columns the failure left uneven.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "sequence_table_reader.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <string>
#include <vector>

namespace {

using namespace duckdb;
using Catch::Matchers::ContainsSubstring;

void Exec(Connection &con, const std::string &sql) {
	auto result = con.Query(sql);
	INFO(sql);
	REQUIRE_FALSE(result->HasError());
}

//! Every sub-batch the stream hands out, until it returns an empty one.
//!
//! A table function builds the stream inside the caller's query, where the caller's context has an
//! active transaction (the schema check and the TEMP-object inheritance read its catalog), so the
//! test opens one around it.
std::vector<miint::SequenceRecordBatch> Batches(Connection &con, const std::string &table, idx_t rows, idx_t bases) {
	std::vector<miint::SequenceRecordBatch> out;
	con.BeginTransaction();
	{
		auto schema = ValidateSequenceTableSchema(*con.context, table, /*allow_bigint=*/false);
		QuerySequenceStream stream(*con.context, table, schema, rows, bases);
		while (true) {
			auto batch = stream.FetchSubBatch();
			if (batch.size() == 0) {
				break;
			}
			out.push_back(std::move(batch));
		}
	}
	con.Commit();
	return out;
}

std::vector<size_t> Sizes(const std::vector<miint::SequenceRecordBatch> &batches) {
	std::vector<size_t> sizes;
	for (auto &b : batches) {
		sizes.push_back(b.size());
	}
	return sizes;
}

//! The batches' ids run together, to compare against the relation's own order.
std::vector<std::string> Flatten(const std::vector<miint::SequenceRecordBatch> &batches) {
	std::vector<std::string> ids;
	for (auto &b : batches) {
		ids.insert(ids.end(), b.read_ids.begin(), b.read_ids.end());
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

//! `count` batches of `size` rows, then one of `rest` if it is non-zero.
std::vector<size_t> Expect(size_t count, size_t size, size_t rest = 0) {
	std::vector<size_t> sizes(count, size);
	if (rest > 0) {
		sizes.push_back(rest);
	}
	return sizes;
}

} // namespace

TEST_CASE("QuerySequenceStream closes a sub-batch by bases", "[sequence_table_reader]") {
	DuckDB db(nullptr);
	Connection con(db);
	const idx_t v = STANDARD_VECTOR_SIZE;

	SECTION("a row that would pass the limit starts the next batch") {
		MakeReads(con, "reads", 10, 100);
		auto batches = Batches(con, "reads", v, 250);
		// 100, 200 fit; a third row would make 300.
		REQUIRE(Sizes(batches) == Expect(5, 2));
		REQUIRE(Flatten(batches) == Ids(10));
	}

	SECTION("a batch may fill the limit exactly") {
		MakeReads(con, "reads", 10, 100);
		auto batches = Batches(con, "reads", v, 300);
		REQUIRE(Sizes(batches) == Expect(3, 3, 1));
		REQUIRE(Flatten(batches) == Ids(10));
	}

	SECTION("a row longer than the limit is a batch of its own") {
		MakeReads(con, "reads", 3, 1000);
		auto batches = Batches(con, "reads", v, 250);
		REQUIRE(Sizes(batches) == Expect(3, 1));
		REQUIRE(Flatten(batches) == Ids(3));
	}

	SECTION("the row limit still applies when it is reached first") {
		MakeReads(con, "reads", 5, 100);
		auto batches = Batches(con, "reads", 2, 1000);
		REQUIRE(Sizes(batches) == Expect(2, 2, 1));
		REQUIRE(Flatten(batches) == Ids(5));
	}

	SECTION("bases cuts across chunk boundaries, before and after rows are dropped") {
		// No chunk holds more than v rows, so a batch of one and a half chunks' worth is assembled
		// from at least two chunks and cut partway through one. The second batch starts with the
		// rows left over from that cut and needs a third chunk, which arrives only after the first
		// batch's rows are dropped, so its cut is counted from the new front.
		const idx_t first = v + v / 2;
		MakeReads(con, "reads", 3 * first, 10);
		auto batches = Batches(con, "reads", 1000000000, first * 10);
		REQUIRE(Sizes(batches) == Expect(3, first));
		REQUIRE(Flatten(batches) == Ids(3 * first));
	}

	SECTION("rows only, across chunks, at a size that does not divide a chunk") {
		const idx_t n = 2 * v + 3;
		MakeReads(con, "reads", n, 10);
		auto batches = Batches(con, "reads", 1000, 0);
		REQUIRE(Sizes(batches) == Expect(n / 1000, 1000, n % 1000));
		REQUIRE(Flatten(batches) == Ids(n));
	}

	SECTION("only sequence1 counts") {
		// 50 bases of sequence1 per row: two rows fill 100. Counting the 50-base mate as well would
		// put one row in each batch.
		Exec(con, "CREATE TABLE pairs AS SELECT 'r' || i AS read_id, repeat('A', 50) AS sequence1, "
		          "repeat('C', 50) AS sequence2 FROM range(4) t(i) ORDER BY i");
		auto batches = Batches(con, "pairs", v, 100);
		REQUIRE(Sizes(batches) == Expect(2, 2));
		REQUIRE(Flatten(batches) == Ids(4));
	}

	SECTION("every column stays with its row") {
		// Paired, with qualities, and every column a function of the row number, so a column moved
		// out of step at a cut or a drop shows as a mismatch on the row it landed on. Cut by both
		// limits over more than two chunks.
		const idx_t n = 2 * v + 5;
		Exec(con, "CREATE TABLE full_rows AS SELECT 'r' || i AS read_id, repeat('A', 1 + i % 7) AS sequence1, "
		          "repeat('C', 1 + i % 5) AS sequence2, [(i % 90)::UTINYINT] AS qual1, "
		          "[(i % 60)::UTINYINT] AS qual2 FROM range(" +
		              std::to_string(n) + ") t(i) ORDER BY i");
		auto batches = Batches(con, "full_rows", 700, 1000);
		REQUIRE(Flatten(batches) == Ids(n));
		idx_t i = 0;
		for (auto &b : batches) {
			REQUIRE(b.is_paired);
			REQUIRE(b.sequences1.size() == b.size());
			REQUIRE(b.sequences2.size() == b.size());
			REQUIRE(b.quals1.size() == b.size());
			REQUIRE(b.quals2.size() == b.size());
			for (idx_t row = 0; row < b.size(); row++, i++) {
				INFO("row " << i);
				REQUIRE(b.sequences1[row] == std::string(1 + i % 7, 'A'));
				REQUIRE(b.sequences2[row] == std::string(1 + i % 5, 'C'));
				REQUIRE(b.quals1[row].as_string() == std::string(1, static_cast<char>(33 + i % 90)));
				REQUIRE(b.quals2[row].as_string() == std::string(1, static_cast<char>(33 + i % 60)));
			}
		}
	}
}

TEST_CASE("QuerySequenceStream raises on every call after a chunk fails", "[sequence_table_reader]") {
	DuckDB db(nullptr);
	Connection con(db);
	const idx_t v = STANDARD_VECTOR_SIZE;
	// A quality of 94 is out of range and throws while its chunk is converted, after that row's
	// read_id and sequence are already held. The bad row is in the second chunk, so the first
	// batches succeed.
	Exec(con, "CREATE TABLE reads AS SELECT 'r' || i AS read_id, 'ACGT' AS sequence1, "
	          "[(CASE WHEN i = " +
	              std::to_string(v + 5) + " THEN 94 ELSE 30 END)::UTINYINT] AS qual1 FROM range(" +
	              std::to_string(2 * v) + ") t(i) ORDER BY i");
	con.BeginTransaction();
	{
		auto schema = ValidateSequenceTableSchema(*con.context, "reads", /*allow_bigint=*/false);
		QuerySequenceStream stream(*con.context, "reads", schema, 100, 0);
		bool raised = false;
		idx_t delivered = 0;
		while (!raised) {
			try {
				auto batch = stream.FetchSubBatch();
				REQUIRE(batch.size() > 0);
				delivered += batch.size();
			} catch (const std::exception &) {
				raised = true;
			}
		}
		// Only rows from the first chunk were handed out.
		REQUIRE(delivered <= v);
		// And nothing after: neither rows nor an empty batch, which would read as a clean end.
		REQUIRE_THROWS_WITH(stream.FetchSubBatch(), ContainsSubstring("earlier batch"));
		REQUIRE_THROWS_WITH(stream.FetchSubBatch(), ContainsSubstring("earlier batch"));
	}
	con.Rollback();
}
