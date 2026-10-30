// Unit tests for RypeInputStream's filter streams (src/include/rype_input_stream.hpp): the input
// rype_filter hands RYpe one group at a time, holding each group's rows so that the ones RYpe does
// not match can be returned whole.
//
// What has to hold, and why:
//  - Put together, the groups are the relation's rows, each exactly once and in scan order, and the
//    Arrow ids run on without a gap. rype_filter maps RYpe's hits (Arrow ids) back to rows through
//    the held runs, so a row missing, repeated or off by one is a read silently kept or dropped.
//  - A group ends with end-of-stream after its rows, and the next group carries on at the next row,
//    mid-chunk if need be. One RYpe call must see exactly one group, whatever its own pull timing.
//  - A held row is the row whose sequence went to Arrow under that id.
//  - The held chunk carries every column, typed as the relation types it; Arrow gets only the
//    sequences.
//
// The tests play RYpe's part over the C stream interface: they move the struct out of the stream
// as arrow-rs's from_raw does, pull to end-of-stream, and release the copy.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "rype_common.hpp"
#include "rype_input_stream.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <algorithm>
#include <condition_variable>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using namespace duckdb;
using Catch::Matchers::ContainsSubstring;

constexpr idx_t ROWS = 5000;

// What one group sent through the C stream, as RYpe sees it.
struct SentGroup {
	std::vector<std::string> formats;
	std::vector<int64_t> ids;
	std::vector<std::string> sequences;
	//! Empty unless the stream sent a pair_sequence column; a NULL mate reads as "<null>".
	std::vector<std::string> pair_sequences;
};

bool RowValid(const ArrowArray &array, idx_t row) {
	if (array.null_count == 0 || !array.buffers[0]) {
		return true;
	}
	const auto bits = static_cast<const uint8_t *>(array.buffers[0]);
	const auto i = static_cast<idx_t>(array.offset) + row;
	return (bits[i / 8] >> (i % 8)) & 1;
}

std::string ReadLargeBinary(const ArrowArray &array, idx_t row) {
	if (!RowValid(array, row)) {
		return "<null>";
	}
	const auto offsets = static_cast<const int64_t *>(array.buffers[1]);
	const auto data = static_cast<const char *>(array.buffers[2]);
	const auto i = static_cast<idx_t>(array.offset) + row;
	return std::string(data + offsets[i], static_cast<size_t>(offsets[i + 1] - offsets[i]));
}

// Take the stream's struct the way rype_classify_arrow_ex does (arrow-rs's from_raw moves it out
// and leaves an empty one behind), so this returns the copy that RYpe would pull and release.
ArrowArrayStream TakeStream(RypeInputStream &input) {
	ArrowArrayStream taken = input.stream;
	input.stream = ArrowArrayStream {};
	return taken;
}

// Consume one group as RYpe does: read the schema, pull to end-of-stream, release.
SentGroup ConsumeGroup(RypeInputStream &input) {
	auto taken = TakeStream(input);
	SentGroup sent;

	ArrowSchema schema;
	REQUIRE(taken.get_schema(&taken, &schema) == 0);
	for (int64_t c = 0; c < schema.n_children; c++) {
		sent.formats.emplace_back(schema.children[c]->format);
	}
	schema.release(&schema);

	for (;;) {
		ArrowArray batch;
		const int rc = taken.get_next(&taken, &batch);
		INFO("get_next: " << (rc == 0 ? "" : taken.get_last_error(&taken)));
		REQUIRE(rc == 0);
		if (!batch.release) {
			break;
		}
		REQUIRE(static_cast<size_t>(batch.n_children) == sent.formats.size());
		REQUIRE(batch.n_children >= 2);
		const auto &ids = *batch.children[0];
		const auto id_data = static_cast<const int64_t *>(ids.buffers[1]);
		for (idx_t r = 0; r < static_cast<idx_t>(batch.length); r++) {
			sent.ids.push_back(id_data[ids.offset + r]);
			sent.sequences.push_back(ReadLargeBinary(*batch.children[1], r));
			if (batch.n_children > 2) {
				sent.pair_sequences.push_back(ReadLargeBinary(*batch.children[2], r));
			}
		}
		batch.release(&batch);
	}

	// End-of-stream stays put: a further pull still ends rather than starting the next group.
	ArrowArray again;
	REQUIRE(taken.get_next(&taken, &again) == 0);
	CHECK(again.release == nullptr);

	taken.release(&taken);
	CHECK(taken.release == nullptr);
	return sent;
}

// Read rows [from, to) of `run`'s column `col` as strings, NULL as "<null>".
std::vector<std::string> HeldStrings(const RypeHeldRows &run, idx_t col) {
	std::vector<std::string> out;
	for (idx_t i = run.from; i < run.to; i++) {
		auto v = run.chunk->GetValue(col, i);
		out.push_back(v.IsNull() ? "<null>" : v.ToString());
	}
	return out;
}

struct FilterStreamFixture {
	DuckDB db;
	//! The caller's connection, which makes the relations.
	Connection con;
	//! The connection the stream scans on, as rype_filter's InitGlobal makes one.
	Connection input;

	FilterStreamFixture() : db(nullptr), con(db), input(db) {
		// One thread, so a table's scan order is its insertion order and can be checked against
		// an ordinal column.
		Exec("SET threads = 1");
		ConfigureRypeArrowExport(input);
	}

	void Exec(const std::string &sql) {
		auto result = con.Query(sql);
		INFO(sql << ": " << (result->HasError() ? result->GetError() : ""));
		REQUIRE_FALSE(result->HasError());
	}

	unique_ptr<RypeInputStream> Open(const std::string &relation, const vector<string> &columns, idx_t group_rows,
	                                 idx_t batch_rows = STANDARD_VECTOR_SIZE, idx_t batch_bytes = 0) {
		RypeInputStreamOptions options;
		options.relation_quoted = SQLIdentifier::ToString(relation);
		options.source_name = relation;
		options.batch_rows = batch_rows;
		options.batch_bytes = batch_bytes;
		auto stream = BuildRypeFilterStream(input, columns, std::move(options));
		stream->SetGroupRows(group_rows);
		return stream;
	}
};

// Pull every group of `stream`, a scan of `SELECT i AS ord, 'S' || i AS sequence1 FROM range(rows)`, and
// check that together they are its rows, each once and in scan order, under ids that run on.
void ExpectWholeGroups(RypeInputStream &stream, idx_t group_rows, idx_t rows) {
	idx_t next = 0;
	idx_t groups = 0;
	while (stream.BeginGroup()) {
		INFO("group " << groups);
		// Once every row has gone out there is no group to begin; an empty one would loop forever.
		REQUIRE(next < rows);
		const auto sent = ConsumeGroup(stream);
		const auto held = stream.TakeGroup();

		const idx_t expect = MinValue<idx_t>(group_rows, rows - next);
		std::vector<int64_t> expect_ids;
		std::vector<std::string> expect_sequences;
		for (idx_t i = 0; i < expect; i++) {
			expect_ids.push_back(static_cast<int64_t>(next + i));
			expect_sequences.push_back("S" + std::to_string(next + i));
		}
		CHECK(sent.formats == std::vector<std::string> {"l", "Z"});
		CHECK(sent.ids == expect_ids);
		CHECK(sent.sequences == expect_sequences);

		// The held runs, read in order, are the same rows under the same ids.
		std::vector<int64_t> held_ids;
		std::vector<int64_t> held_ords;
		std::vector<std::string> held_sequences;
		for (idx_t r = 0; r < held.size(); r++) {
			const auto &run = held[r];
			REQUIRE(run.from < run.to);
			REQUIRE(run.to <= run.chunk->size());
			// One run per chunk, however many Arrow batches it took: each run becomes one unit of
			// output work, so a byte ceiling of one row per batch must not make one unit per row.
			if (r > 0) {
				CHECK(held[r - 1].chunk.get() != run.chunk.get());
			}
			for (idx_t i = run.from; i < run.to; i++) {
				held_ids.push_back(static_cast<int64_t>(run.first_id + (i - run.from)));
				held_ords.push_back(run.chunk->GetValue(0, i).GetValue<int64_t>());
			}
			const auto strings = HeldStrings(run, 1);
			held_sequences.insert(held_sequences.end(), strings.begin(), strings.end());
		}
		CHECK(held_ids == expect_ids);
		CHECK(held_ords == expect_ids);
		CHECK(held_sequences == expect_sequences);

		next += expect;
		groups++;
	}
	CHECK(next == rows);
	CHECK(groups == (rows + group_rows - 1) / group_rows);
}

// Pull groups as rype_filter does until the stream ends or fails, counting the rows sent. Returns the
// error, or "" for a clean end.
std::string PullUntilError(RypeInputStream &stream, idx_t &rows) {
	std::string error;
	try {
		while (error.empty() && stream.BeginGroup()) {
			auto taken = TakeStream(stream);
			for (;;) {
				ArrowArray batch;
				if (taken.get_next(&taken, &batch) != 0) {
					error = taken.get_last_error(&taken);
					break;
				}
				if (!batch.release) {
					break;
				}
				rows += static_cast<idx_t>(batch.length);
				batch.release(&batch);
			}
			taken.release(&taken);
			stream.TakeGroup();
		}
	} catch (std::exception &e) {
		// The error may also surface while BeginGroup looks for the next group's first row.
		error = e.what();
	}
	return error;
}

// Runs ReadAhead on a thread of its own whenever the stream asks, as rype_filter's job on DuckDB's
// async pool does. Must go before the stream does.
class ReadAheadThread {
public:
	ReadAheadThread(RypeInputStream &stream_p, idx_t bytes) : stream(stream_p) {
		stream.SetReadAhead(bytes, [this]() {
			std::lock_guard<std::mutex> guard(lock);
			pending++;
			most_pending = std::max(most_pending, pending);
			requests++;
			wake.notify_all();
		});
		worker = std::thread([this]() { Run(); });
	}

	~ReadAheadThread() {
		{
			std::lock_guard<std::mutex> guard(lock);
			stop = true;
		}
		wake.notify_all();
		worker.join();
	}

	//! Requests not yet served, at most at once, and in all. Read on the test's thread.
	idx_t MostPending() {
		std::lock_guard<std::mutex> guard(lock);
		return most_pending;
	}
	idx_t Requests() {
		std::lock_guard<std::mutex> guard(lock);
		return requests;
	}

private:
	void Run() {
		std::unique_lock<std::mutex> guard(lock);
		for (;;) {
			wake.wait(guard, [this]() { return pending > 0 || stop; });
			if (pending == 0) {
				return;
			}
			pending--;
			guard.unlock();
			stream.ReadAhead();
			guard.lock();
		}
	}

	RypeInputStream &stream;
	std::mutex lock;
	std::condition_variable wake;
	idx_t pending = 0;
	idx_t most_pending = 0;
	idx_t requests = 0;
	bool stop = false;
	std::thread worker;
};

} // namespace

TEST_CASE("Filter stream groups are the relation's rows, each once and in order, with ids that run on",
          "[RypeInputStream]") {
	FilterStreamFixture f;
	f.Exec("CREATE TABLE reads AS SELECT i AS ord, 'S' || i AS sequence1 FROM range(" + std::to_string(ROWS) +
	       ") t(i)");

	// Group sizes either side of a 2048-row scan chunk, and ones that do not divide it, so groups
	// end mid-chunk and chunks straddle groups.
	const idx_t group_rows = GENERATE(as<idx_t> {}, 1, 7, 2047, 2048, 2049, 3000, 5000, 6000);
	// The Arrow batch's own ceilings close batches inside a group; the group must still be whole.
	const auto ceiling = GENERATE(std::pair<idx_t, idx_t>(STANDARD_VECTOR_SIZE, 0), std::pair<idx_t, idx_t>(5, 0),
	                              std::pair<idx_t, idx_t>(STANDARD_VECTOR_SIZE, 1));
	INFO("group_rows=" << group_rows << " batch_rows=" << ceiling.first << " batch_bytes=" << ceiling.second);

	auto stream = f.Open("reads", {"ord", "sequence1"}, group_rows, ceiling.first, ceiling.second);
	ExpectWholeGroups(*stream, group_rows, ROWS);
}

TEST_CASE("A group that ends mid-chunk leaves the rest of the chunk to the next group, shared, not copied",
          "[RypeInputStream]") {
	FilterStreamFixture f;
	// Twenty rows come back as one scan chunk, so groups of seven end inside it twice.
	f.Exec("CREATE TABLE reads AS SELECT i AS ord, 'S' || i AS sequence1 FROM range(20) t(i)");
	auto stream = f.Open("reads", {"ord", "sequence1"}, 7);

	std::vector<std::vector<RypeHeldRows>> groups;
	// Bounded, so a stream that kept beginning groups fails the size check below instead of hanging.
	while (groups.size() <= 3 && stream->BeginGroup()) {
		ConsumeGroup(*stream);
		groups.push_back(stream->TakeGroup());
	}
	REQUIRE(groups.size() == 3);
	for (const auto &g : groups) {
		REQUIRE(g.size() == 1);
	}
	const auto &first = groups[0][0];
	REQUIRE(first.chunk->size() == 20);
	CHECK((first.from == 0 && first.to == 7 && first.first_id == 0));
	CHECK((groups[1][0].from == 7 && groups[1][0].to == 14 && groups[1][0].first_id == 7));
	CHECK((groups[2][0].from == 14 && groups[2][0].to == 20 && groups[2][0].first_id == 14));
	CHECK(groups[1][0].chunk.get() == first.chunk.get());
	CHECK(groups[2][0].chunk.get() == first.chunk.get());
}

TEST_CASE("A filter stream holds every column as the relation types it, and sends Arrow only the sequences",
          "[RypeInputStream]") {
	FilterStreamFixture f;
	// sequence1 and sequence2 sit away from positions 1 and 2, and sequence2 comes first, so a stream
	// that read them by position would send the wrong column.
	f.Exec("CREATE TABLE reads (extra INTEGER, sequence2 VARCHAR, read_id VARCHAR, sequence1 VARCHAR, "
	       "qual1 UTINYINT[])");
	f.Exec("INSERT INTO reads VALUES (10, 'TTTT', 'a', 'ACGTAC', [30, 31]), (11, NULL, 'b', 'GG', [32]), "
	       "(12, 'CC', 'c', 'TTTTTTTTTT', [])");
	const vector<string> columns {"extra", "sequence2", "read_id", "sequence1", "qual1"};
	auto stream = f.Open("reads", columns, 100);

	const auto sizing = stream->SampleSizing();
	CHECK(sizing.sampled);
	CHECK(sizing.avg_read_length == (6 + 2 + 10) / 3);
	CHECK(sizing.is_paired);

	REQUIRE(stream->BeginGroup());
	const auto sent = ConsumeGroup(*stream);
	// LargeBinary for both, sequence1 included though it is VARCHAR: it is reinterpreted, not
	// exported as text.
	CHECK(sent.formats == std::vector<std::string> {"l", "Z", "Z"});
	// Sampling consumed nothing: the first group still starts at the first row.
	CHECK(sent.ids == std::vector<int64_t> {0, 1, 2});
	CHECK(sent.sequences == std::vector<std::string> {"ACGTAC", "GG", "TTTTTTTTTT"});
	CHECK(sent.pair_sequences == std::vector<std::string> {"TTTT", "<null>", "CC"});

	const auto held = stream->TakeGroup();
	REQUIRE(held.size() == 1);
	const auto &chunk = *held[0].chunk;
	REQUIRE(chunk.ColumnCount() == 5);
	CHECK(chunk.data[0].GetType() == LogicalType::INTEGER);
	CHECK(chunk.data[1].GetType() == LogicalType::VARCHAR);
	CHECK(chunk.data[2].GetType() == LogicalType::VARCHAR);
	CHECK(chunk.data[3].GetType() == LogicalType::VARCHAR);
	CHECK(chunk.data[4].GetType() == LogicalType::LIST(LogicalType::UTINYINT));
	CHECK(HeldStrings(held[0], 2) == std::vector<std::string> {"a", "b", "c"});
	CHECK(HeldStrings(held[0], 4) == std::vector<std::string> {"[30, 31]", "[32]", "[]"});

	CHECK_FALSE(stream->BeginGroup());
}

TEST_CASE("SampleSizing reports the bytes a held row costs, every column counted", "[RypeInputStream]") {
	// rype_filter sizes its groups by what holding a row costs, and a held row carries every column,
	// not just the sequences RYpe is sent. Two relations with the same sequences, one with a 1000-byte
	// extra column, must differ by at least that much per row.
	FilterStreamFixture f;
	f.Exec("CREATE TABLE lean AS SELECT repeat('ACGT', 25) AS sequence1 FROM range(3000) t(i)");
	f.Exec("CREATE TABLE wide AS SELECT repeat('ACGT', 25) AS sequence1, repeat('x', 1000) || i AS note "
	       "FROM range(3000) t(i)");
	const auto lean = f.Open("lean", {"sequence1"}, 7)->SampleSizing();
	const auto wide = f.Open("wide", {"sequence1", "note"}, 7)->SampleSizing();
	CHECK(lean.avg_read_length == 100);
	CHECK(wide.avg_read_length == 100);
	CHECK(lean.held_bytes_per_row >= 100);
	CHECK(wide.held_bytes_per_row >= lean.held_bytes_per_row + 1000);
}

TEST_CASE("A BLOB sequence1 with no sequence2 goes to Arrow as a single-end batch", "[RypeInputStream]") {
	FilterStreamFixture f;
	f.Exec("CREATE TABLE reads AS SELECT 'r' || i AS read_id, ('S' || i)::BLOB AS sequence1 FROM range(3) t(i)");
	auto stream = f.Open("reads", {"read_id", "sequence1"}, 100);
	CHECK_FALSE(stream->SampleSizing().is_paired);

	REQUIRE(stream->BeginGroup());
	const auto sent = ConsumeGroup(*stream);
	CHECK(sent.formats == std::vector<std::string> {"l", "Z"});
	CHECK(sent.sequences == std::vector<std::string> {"S0", "S1", "S2"});
	CHECK(sent.pair_sequences.empty());
	const auto held = stream->TakeGroup();
	REQUIRE(held.size() == 1);
	CHECK(held[0].chunk->data[1].GetType() == LogicalType::BLOB);
}

TEST_CASE("An empty relation begins no group", "[RypeInputStream]") {
	FilterStreamFixture f;
	f.Exec("CREATE TABLE reads (sequence1 VARCHAR)");
	auto stream = f.Open("reads", {"sequence1"}, 7);
	CHECK_FALSE(stream->SampleSizing().sampled);
	CHECK_FALSE(stream->BeginGroup());
}

TEST_CASE("A filter stream refuses a relation whose sequence1 is missing or not text", "[RypeInputStream]") {
	FilterStreamFixture f;
	f.Exec("CREATE TABLE no_sequence AS SELECT 1 AS read_id");
	f.Exec("CREATE TABLE int_sequence AS SELECT 1 AS sequence1");
	CHECK_THROWS_WITH(f.Open("no_sequence", {"read_id"}, 7),
	                  ContainsSubstring("no_sequence") && ContainsSubstring("sequence1"));
	CHECK_THROWS_WITH(f.Open("int_sequence", {"sequence1"}, 7),
	                  ContainsSubstring("int_sequence") && ContainsSubstring("INTEGER"));
}

TEST_CASE("A NULL sequence1 fails the pull, naming the row", "[RypeInputStream]") {
	FilterStreamFixture f;
	f.Exec("CREATE TABLE reads AS SELECT CASE WHEN i = 3 THEN NULL ELSE 'S' || i END AS sequence1 FROM range(6) t(i)");
	auto stream = f.Open("reads", {"sequence1"}, 100);
	REQUIRE(stream->BeginGroup());
	auto taken = TakeStream(*stream);
	ArrowArray batch;
	REQUIRE(taken.get_next(&taken, &batch) != 0);
	const std::string error = taken.get_last_error(&taken);
	CHECK_THAT(error,
	           ContainsSubstring("sequence1 is NULL") && ContainsSubstring("'reads'") && ContainsSubstring("row 3"));
	taken.release(&taken);
}

TEST_CASE("A scan error in a later group is reported, never taken for the end of the input", "[RypeInputStream]") {
	FilterStreamFixture f;
	// The error fires on row 4000, two scan chunks in, after earlier groups have gone out cleanly.
	f.Exec("CREATE VIEW reads AS SELECT CASE WHEN i = 4000 THEN error('scan failed at 4000') ELSE 'S' || i END "
	       "AS sequence1 FROM range(" +
	       std::to_string(ROWS) + ") t(i)");
	auto stream = f.Open("reads", {"sequence1"}, 7);

	idx_t rows = 0;
	const auto error = PullUntilError(*stream, rows);
	INFO("rows before the error: " << rows);
	CHECK(rows < 4000);
	CHECK_THAT(error, ContainsSubstring("scan failed at 4000"));
}

TEST_CASE("Read-ahead fetches the scan while nothing is pulled, and stops at its bytes", "[RypeInputStream]") {
	// A streaming result produces only as it is fetched from, so a scan fetched from get_next alone
	// stands still while RYpe classifies. ReadAhead has to move the scan on with nothing pulled, and
	// stop at its bytes, which rype_filter's memory budget counts.
	FilterStreamFixture f;
	const idx_t rows = 10 * STANDARD_VECTOR_SIZE;
	f.Exec("CREATE TABLE reads AS SELECT i AS ord, 'S' || i AS sequence1 FROM range(" + std::to_string(rows) +
	       ") t(i)");

	// What one scan chunk of this table costs: a read-ahead of one byte fetches exactly one.
	idx_t chunk_bytes = 0;
	{
		auto probe = f.Open("reads", {"ord", "sequence1"}, 7);
		probe->SetReadAhead(1, []() {});
		probe->ReadAhead();
		chunk_bytes = probe->ReadAheadBytes();
	}
	REQUIRE(chunk_bytes > 0);

	auto stream = f.Open("reads", {"ord", "sequence1"}, 7);
	idx_t starts = 0;
	stream->SetReadAhead(3 * chunk_bytes, [&starts]() { starts++; });
	stream->ReadAhead();
	// Three chunks ahead, with nothing pulled, and not a fourth.
	CHECK(stream->ReadAheadBytes() >= 3 * chunk_bytes);
	CHECK(stream->ReadAheadBytes() < 4 * chunk_bytes);
	// Called directly, not asked for.
	CHECK(starts == 0);

	// Nothing serves the request below, so the pulls fetch every later chunk themselves: a pull never
	// waits for a ReadAhead that has not started, which on a busy pool may not start for a while.
	ExpectWholeGroups(*stream, 7, rows);
	// Asked for once, when the chunks read ahead fell to half the bytes, and not again while unserved.
	CHECK(starts == 1);
}

TEST_CASE("Groups stay whole while read-ahead runs on another thread", "[RypeInputStream]") {
	FilterStreamFixture f;
	const idx_t rows = 10 * STANDARD_VECTOR_SIZE + 99;
	f.Exec("CREATE TABLE reads AS SELECT i AS ord, 'S' || i AS sequence1 FROM range(" + std::to_string(rows) +
	       ") t(i)");
	const idx_t group_rows = GENERATE(as<idx_t> {}, 7, 2047, 3000);
	// One byte keeps a single chunk ahead, so read-ahead is asked for at nearly every chunk and races
	// the pulls for the scan; the most there is reads the whole scan ahead.
	const idx_t bytes = GENERATE(as<idx_t> {}, 1, NumericLimits<idx_t>::Maximum());
	INFO("group_rows=" << group_rows << " read-ahead bytes=" << bytes);

	auto stream = f.Open("reads", {"ord", "sequence1"}, group_rows);
	ReadAheadThread reader(*stream, bytes);
	ExpectWholeGroups(*stream, group_rows, rows);
	CHECK(reader.Requests() > 0);
	// A request is not repeated before it is served: rype_filter would otherwise queue a job per chunk.
	CHECK(reader.MostPending() == 1);
}

TEST_CASE("A scan error met reading ahead is reported to the pull that reaches it", "[RypeInputStream]") {
	FilterStreamFixture f;
	f.Exec("CREATE VIEW reads AS SELECT CASE WHEN i = 4000 THEN error('scan failed at 4000') ELSE 'S' || i END "
	       "AS sequence1 FROM range(" +
	       std::to_string(ROWS) + ") t(i)");
	auto stream = f.Open("reads", {"sequence1"}, 7);
	stream->SetReadAhead(NumericLimits<idx_t>::Maximum(), []() {});
	// Met on the read-ahead's thread and kept there: thrown from ReadAhead, it would be lost, and the
	// pulls would take the end of the chunks for the end of the input.
	REQUIRE_NOTHROW(stream->ReadAhead());

	idx_t rows = 0;
	const auto error = PullUntilError(*stream, rows);
	INFO("rows before the error: " << rows);
	CHECK(rows < 4000);
	CHECK_THAT(error, ContainsSubstring("scan failed at 4000"));
}

TEST_CASE("Abort stops read-ahead and fails the next pull", "[RypeInputStream]") {
	// rype_filter aborts the stream when a query ends early (a LIMIT, an error), and a read-ahead job
	// must then stop rather than keep reading the relation.
	FilterStreamFixture f;
	f.Exec("CREATE TABLE reads AS SELECT 'S' || i AS sequence1 FROM range(" + std::to_string(ROWS) + ") t(i)");
	auto stream = f.Open("reads", {"sequence1"}, 7);
	stream->SetReadAhead(NumericLimits<idx_t>::Maximum(), []() {});
	REQUIRE(stream->BeginGroup());
	stream->Abort();
	stream->ReadAhead();
	CHECK(stream->ReadAheadBytes() == 0);

	auto taken = TakeStream(*stream);
	ArrowArray batch;
	CHECK(taken.get_next(&taken, &batch) != 0);
	taken.release(&taken);
}
