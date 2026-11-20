// COPY FORMAT SAM/BAM keep a batch's records (SAM lines, built BAM records) in prepare_batch, in a buffer outside
// DuckDB's memory_limit, up to 64 MiB: the record that would pass it, and every row after it, prepare_batch keeps as
// rows, and flush_batch builds them again as it writes them. DuckDB caps a batch by its input rows, which
// SEQUENCE_DATA's reads are not part of, so a batch of small rows can make records of any size; the SQL tests compare
// the two paths' output but cannot see which one a batch took. These tests drive RecordPrepareBatch with a stand-in
// encoder whose records have a set size, over batches that span several chunks, which pins the decision on both sides
// of 64 MiB without 64 MiB of real records, and that every row is still checked while the batch is prepared. Then
// RecordFlushBatch, which writes the kept records and builds the deferred rows again, skipping the chunks that hold
// only kept rows: each row must be written once, in order, whether deferral starts at a chunk's first row or inside
// a chunk. The SQL tests' deferred batches are each one chunk.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "copy_format_common.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <cstring>
#include <string>

namespace {

using namespace duckdb;
using Catch::Matchers::ContainsSubstring;

constexpr idx_t MIB = 1024 * 1024;
constexpr idx_t THRESHOLD = 64 * MIB;
// Rows that make exactly THRESHOLD bytes at RECORD_BYTES each; two chunks' worth, so a batch spans chunks.
constexpr idx_t ROWS = 2 * STANDARD_VECTOR_SIZE;
constexpr idx_t RECORD_BYTES = THRESHOLD / ROWS;

// Each row is one record of `record_bytes`, which starts with the row's index in the batch (rows from Rows below hold
// it). Counts the rows it encodes and the rows it only checks, and fails the check of row `bad_row`, if set.
class SizedRecordEncoder {
public:
	explicit SizedRecordEncoder(idx_t record_bytes_p = RECORD_BYTES, idx_t bad_row_p = DConstants::INVALID_INDEX)
	    : record_bytes(record_bytes_p), bad_row(bad_row_p) {
	}

	void SetChunk(DataChunk &chunk_p) {
		chunk = &chunk_p;
	}

	void Encode(idx_t row, RecordBuffer &out) {
		auto record = out.Reserve(record_bytes);
		std::memset(record, 'A', record_bytes);
		Store<uint64_t>(RowIndex(row), record);
		encoded++;
		Fail(row);
	}

	void Check(idx_t row) {
		checked++;
		Fail(row);
	}

	idx_t encoded = 0;
	idx_t checked = 0;

private:
	idx_t RowIndex(idx_t row) const {
		return NumericCast<idx_t>(chunk->data[0].GetValue(row).GetValue<int64_t>());
	}

	void Fail(idx_t row) {
		if (RowIndex(row) == bad_row) {
			throw InvalidInputException("bad row %llu", static_cast<unsigned long long>(bad_row));
		}
	}

	idx_t record_bytes;
	idx_t bad_row;
	DataChunk *chunk = nullptr;
};

// `count` rows, row i holding i, in a collection of more than one chunk.
unique_ptr<ColumnDataCollection> Rows(Connection &con, idx_t count) {
	auto result = con.Query("SELECT i FROM range(" + std::to_string(count) + ") t(i)");
	INFO(result->GetError());
	REQUIRE_FALSE(result->HasError());
	auto rows = result->TakeCollection();
	REQUIRE(rows->ChunkCount() > 1);
	return rows;
}

} // namespace

TEST_CASE("SAM/BAM records over 64 MiB are built again in flush_batch, and checked in prepare_batch", "[copy_sam]") {
	DuckDB db(nullptr);
	Connection con(db);
	auto &context = *con.context;

	SECTION("64 MiB exactly is kept") {
		SizedRecordEncoder encoder;
		auto prepared = RecordPrepareBatch(context, encoder, Rows(con, ROWS));
		auto &batch = prepared->Cast<RecordBatchData>();
		REQUIRE_FALSE(batch.deferred);
		REQUIRE(batch.records.Size() == THRESHOLD);
		REQUIRE(batch.records.Count() == ROWS);
		REQUIRE(encoder.encoded == ROWS);
		REQUIRE(encoder.checked == 0);
	}

	SECTION("past 64 MiB, the records that fit are kept, the rows from there on deferred, and every later row "
	        "checked") {
		SizedRecordEncoder encoder;
		auto prepared = RecordPrepareBatch(context, encoder, Rows(con, 2 * ROWS));
		auto &batch = prepared->Cast<RecordBatchData>();
		REQUIRE(batch.deferred);
		REQUIRE(batch.deferred->Count() == 2 * ROWS);
		// The first ROWS records make exactly 64 MiB and are kept; record ROWS passed it, and is built again with the
		// rows after it.
		REQUIRE(batch.records.Count() == ROWS);
		REQUIRE(batch.records.Size() == THRESHOLD);
		REQUIRE(batch.deferred_from == ROWS);
		// That record was the last one encoded; the rows after it were checked, not encoded.
		REQUIRE(encoder.encoded == ROWS + 1);
		REQUIRE(encoder.checked == ROWS - 1);
	}

	SECTION("a row after the records pass 64 MiB that fails its check fails the prepare") {
		SizedRecordEncoder encoder(RECORD_BYTES, 2 * ROWS - 1);
		REQUIRE_THROWS_WITH(RecordPrepareBatch(context, encoder, Rows(con, 2 * ROWS)),
		                    ContainsSubstring("bad row " + std::to_string(2 * ROWS - 1)));
	}
}

TEST_CASE("flush_batch writes a deferred batch's kept records, then its deferred rows, each once and in order",
          "[copy_sam]") {
	DuckDB db(nullptr);
	Connection con(db);
	auto &context = *con.context;

	// Prepares 2 * ROWS rows of `record_bytes` records, deferred from `deferred_from`, flushes them, and returns the
	// index of the first written record that is not the next row (2 * ROWS if every row is written once, in order).
	const auto first_out_of_place = [&](idx_t record_bytes, idx_t deferred_from) {
		SizedRecordEncoder encoder(record_bytes);
		auto prepared = RecordPrepareBatch(context, encoder, Rows(con, 2 * ROWS));
		auto &batch = prepared->Cast<RecordBatchData>();
		REQUIRE(batch.deferred);
		REQUIRE(batch.deferred_from == deferred_from);
		vector<idx_t> written;
		RecordFlushBatch(
		    context, batch, [&]() { return SizedRecordEncoder(record_bytes); },
		    [&](const_data_ptr_t data, idx_t size) {
			    REQUIRE(size == record_bytes);
			    written.push_back(NumericCast<idx_t>(Load<uint64_t>(data)));
		    });
		idx_t i = 0;
		while (i < written.size() && written[i] == i) {
			i++;
		}
		INFO("records written: " << written.size());
		REQUIRE(written.size() == 2 * ROWS);
		return i;
	};

	SECTION("deferred from the first row of a later chunk: the chunks before it are skipped whole") {
		// ROWS records of RECORD_BYTES make exactly 64 MiB, so row ROWS, the first of the third chunk, passes it.
		REQUIRE(ROWS % STANDARD_VECTOR_SIZE == 0);
		REQUIRE(first_out_of_place(RECORD_BYTES, ROWS) == 2 * ROWS);
	}

	SECTION("deferred from inside a chunk") {
		// 3,355 records of 20,000 bytes fit in 64 MiB (67,100,000 bytes); row 3,355 is the 1,308th of the second chunk.
		constexpr idx_t SIZE = 20000;
		constexpr idx_t FIT = THRESHOLD / SIZE;
		REQUIRE(FIT % STANDARD_VECTOR_SIZE != 0);
		REQUIRE(first_out_of_place(SIZE, FIT) == 2 * ROWS);
	}
}
