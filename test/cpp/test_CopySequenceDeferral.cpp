// COPY FORMAT FASTA/FASTQ format a batch's text in prepare_batch, into buffers outside DuckDB's memory_limit, unless
// that text would pass 64 MiB: then prepare_batch keeps the rows, and flush_batch formats them, writing out each time
// 1 MiB has gathered (a record larger than that is held whole). The SQL tests compare the two paths' output but cannot
// see which one a batch took, so a batch of many records that stopped being deferred would hold all its text in memory
// and fail no SQL test. These tests drive SequencePrepareBatch with a stand-in FORMAT whose records claim a size, over
// batches that span several chunks, which pins the decision on both sides of 64 MiB without 64 MiB of real records.

#include <catch2/catch_test_macros.hpp>

#include "copy_format_common.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <string>

namespace {

using namespace duckdb;

constexpr idx_t MIB = 1024 * 1024;
constexpr idx_t THRESHOLD = 64 * MIB;
// Rows that make exactly THRESHOLD bytes at RECORD_BYTES each; two chunks' worth, so a batch spans chunks.
constexpr idx_t ROWS = 2 * STANDARD_VECTOR_SIZE;
constexpr idx_t RECORD_BYTES = THRESHOLD / ROWS;

// Each row is one record of `record_bytes` (filler text) in R1, and with `paired` one more in R2. Counts Check calls.
class SizedRecordFormat {
public:
	struct Record {};

	SizedRecordFormat(idx_t record_bytes_p, bool paired_p)
	    : record_bytes(record_bytes_p), paired(paired_p), filler(MIB, 'A') {
	}

	template <class EMIT>
	void Visit(DataChunk &chunk, SequenceBatchData &, EMIT &&emit) {
		for (idx_t row = 0; row < chunk.size(); row++) {
			emit(SequenceStream::R1, Record {});
			if (paired) {
				emit(SequenceStream::R2, Record {});
			}
		}
	}

	idx_t Size(const Record &) const {
		return record_bytes;
	}

	void Write(MemoryStream &stream, const Record &) {
		for (idx_t written = 0; written < record_bytes; written += MIB) {
			stream.WriteData(const_data_ptr_cast(filler.data()), MinValue(MIB, record_bytes - written));
		}
	}

	void Check(const Record &) {
		checked++;
	}

	idx_t checked = 0;

private:
	idx_t record_bytes;
	bool paired;
	std::string filler;
};

// `count` rows in a collection of more than one chunk.
unique_ptr<ColumnDataCollection> Rows(Connection &con, idx_t count) {
	auto result = con.Query("SELECT i FROM range(" + std::to_string(count) + ") t(i)");
	INFO(result->GetError());
	REQUIRE_FALSE(result->HasError());
	auto rows = result->TakeCollection();
	REQUIRE(rows->ChunkCount() > 1);
	return rows;
}

} // namespace

TEST_CASE("FASTA/FASTQ text over 64 MiB is formatted in flush_batch, and checked in prepare_batch", "[copy_fastx]") {
	DuckDB db(nullptr);
	Connection con(db);
	auto &context = *con.context;

	SECTION("64 MiB exactly is formatted in prepare_batch") {
		SizedRecordFormat format(RECORD_BYTES, false);
		auto prepared = SequencePrepareBatch(context, format, Rows(con, ROWS));
		auto &batch = prepared->Cast<SequenceBatchData>();
		REQUIRE_FALSE(batch.deferred);
		REQUIRE(batch.r1);
		REQUIRE(batch.r1->GetPosition() == THRESHOLD);
		REQUIRE(format.checked == 0);
	}

	SECTION("one record past 64 MiB, the rows are kept for flush_batch, and every record is checked first") {
		SizedRecordFormat format(RECORD_BYTES, false);
		auto prepared = SequencePrepareBatch(context, format, Rows(con, ROWS + 1));
		auto &batch = prepared->Cast<SequenceBatchData>();
		REQUIRE(batch.deferred);
		REQUIRE(batch.deferred->Count() == ROWS + 1);
		REQUIRE(format.checked == ROWS + 1);
	}

	SECTION("R1 and R2 count together, and a kept batch has no buffer for either") {
		// Each file's text is just over 32 MiB; only their sum passes 64 MiB.
		SizedRecordFormat format(RECORD_BYTES / 2, true);
		auto prepared = SequencePrepareBatch(context, format, Rows(con, ROWS + 1));
		auto &batch = prepared->Cast<SequenceBatchData>();
		REQUIRE(batch.deferred);
		REQUIRE_FALSE(batch.r1);
		REQUIRE_FALSE(batch.r2);
		REQUIRE(format.checked == 2 * (ROWS + 1));
	}
}
