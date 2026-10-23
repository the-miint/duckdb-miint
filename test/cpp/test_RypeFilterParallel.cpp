// rype_filter on many threads, driven through DuckDB's streaming result.
//
// The SQL tests cannot reach this path: sqllogictest materializes every result. A streamed result
// has a small buffer, so the threads producing it are held back until the reader takes a chunk,
// and the reader may stop, or interrupt the query, with groups still being classified and mapped.
// What has to hold, and why:
//  - The rows come out whole: every surviving row once, in whatever order the threads emit them.
//    A parked thread that is never woken shows up as a hang, and a unit dropped or emitted twice as
//    the wrong rows.
//  - Interrupting the query, or dropping its result half read, ends it promptly, and the
//    connection runs the next query correctly. The background work must stop and be waited for,
//    not left running against freed state.
//
// A hang is the failure looked for, so each query runs against a deadline, past which the process
// aborts with a message instead of leaving the suite to hang.

#include <catch2/catch_test_macros.hpp>

#include "miint_streaming_query.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <future>
#include <string>
#include <thread>
#include <vector>

namespace {

using namespace duckdb;

constexpr auto DEADLINE = std::chrono::seconds(120);

template <class BODY>
void WithDeadline(const char *what, BODY &&body) {
	std::promise<void> done;
	auto finished = done.get_future();
	std::thread watchdog([&]() {
		if (finished.wait_for(DEADLINE) != std::future_status::ready) {
			std::fprintf(stderr, "%s did not finish within %lld s: a thread was parked and never woken?\n", what,
			             static_cast<long long>(DEADLINE.count()));
			std::abort();
		}
	});
	try {
		body();
	} catch (...) {
		done.set_value();
		watchdog.join();
		throw;
	}
	done.set_value();
	watchdog.join();
}

struct ParallelFilterFixture {
	DuckDB db;
	Connection con;
	std::string mmi;
	std::string filter_sql;

	ParallelFilterFixture() : db(nullptr), con(db) {
		mmi = (std::filesystem::temp_directory_path() /
		       ("miint_rype_filter_parallel_" +
		        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".mmi"))
		          .string();
		// The fixture of test/sql/rype_filter_parallel.test: random host sequences no RYpe bucket
		// covers, and reads that are host pieces, random, or ACGT repeats RYpe always matches.
		Exec("CREATE MACRO random_dna(seed, len) AS (SELECT string_agg(substr('ACGT', (hash(seed * 1000003 + j) "
		     "% 4)::INTEGER + 1, 1), '' ORDER BY j) FROM range(len) t(j))");
		Exec("CREATE TABLE host AS SELECT 'host' || h AS read_id, random_dna(h, 3000) AS sequence1 FROM range(3) "
		     "t(h)");
		Exec("CREATE MACRO host_piece(seed, len) AS substr((SELECT sequence1 FROM host WHERE read_id = 'host' || "
		     "(seed % 3)), (seed % (3000 - len))::INTEGER + 1, len)");
		Exec("CREATE MACRO read_of(kind, seed) AS CASE kind WHEN 0 THEN host_piece(seed, 150) WHEN 1 THEN "
		     "reverse(translate(host_piece(seed, 150), 'ACGT', 'TGCA')) WHEN 2 THEN random_dna(seed + 100000, 150) "
		     "ELSE 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT' END");
		Exec("SELECT success FROM save_minimap2_index('host', '" + mmi + "')");
		Exec("CREATE TABLE pairs AS SELECT i AS ord, 'r' || i AS read_id, read_of(hash(i) % 4, i) AS sequence1, "
		     "CASE WHEN i % 7 = 0 THEN NULL ELSE read_of(hash(i + 7919) % 4, i + 50000) END AS sequence2 "
		     "FROM range(3000) t(i) ORDER BY i");
		filter_sql = "SELECT ord FROM rype_filter('pairs', 'data/rype/test.ryxdi', '" + mmi + "', threshold := 0.05)";
		Exec("CREATE TABLE rype_kept AS SELECT * FROM pairs ANTI JOIN rype_classify('data/rype/test.ryxdi', 'pairs', "
		     "threshold := 0.05) USING (read_id)");
		Exec(
		    "CREATE TABLE expected AS SELECT ord FROM rype_kept ANTI JOIN align_minimap2('rype_kept', index_path := '" +
		    mmi + "') USING (read_id)");

		Exec("SET threads = 8");
		// Many small groups, so several are in flight whenever the reader stops.
		Exec("SET miint_rype_filter_group_rows = 7");
		Exec("SET max_streaming_buffer_size = '1KB'");
	}

	~ParallelFilterFixture() {
		std::error_code ignored;
		std::filesystem::remove(mmi, ignored);
	}

	void Exec(const std::string &sql) {
		auto result = con.Query(sql);
		if (result->HasError()) {
			FAIL(sql << ": " << result->GetError());
		}
	}

	//! The first column of a query's rows, which must be BIGINT, sorted.
	std::vector<int64_t> Ords(const std::string &sql) {
		auto stream = SubmitStream(con, sql);
		std::vector<int64_t> ords;
		while (auto chunk = stream->Fetch()) {
			Collect(*chunk, ords);
		}
		INFO(sql << ": " << stream->GetError());
		REQUIRE_FALSE(stream->HasError());
		std::sort(ords.begin(), ords.end());
		return ords;
	}

	static void Collect(DataChunk &chunk, std::vector<int64_t> &ords) {
		for (idx_t row = 0; row < chunk.size(); row++) {
			ords.push_back(chunk.GetValue(0, row).GetValue<int64_t>());
		}
	}
};

} // namespace

TEST_CASE("rype_filter streamed on eight threads comes out whole", "[rype_filter]") {
	ParallelFilterFixture f;
	const auto expected = f.Ords("SELECT ord FROM expected");
	REQUIRE(expected.size() > 100);

	WithDeadline("a streamed rype_filter read slowly", [&]() {
		auto stream = SubmitStream(f.con, f.filter_sql);
		REQUIRE_FALSE(stream->HasError());
		std::vector<int64_t> ords;
		while (auto chunk = stream->Fetch()) {
			ParallelFilterFixture::Collect(*chunk, ords);
			// A slow reader keeps the buffer full, so the threads producing it are held back.
			std::this_thread::sleep_for(std::chrono::milliseconds(1));
		}
		INFO(stream->GetError());
		REQUIRE_FALSE(stream->HasError());
		std::sort(ords.begin(), ords.end());
		REQUIRE(ords == expected);
	});
}

TEST_CASE("rype_filter stopped half read leaves the connection working", "[rype_filter]") {
	ParallelFilterFixture f;
	const auto expected = f.Ords("SELECT ord FROM expected");
	REQUIRE(expected.size() > 100);

	SECTION("interrupted") {
		WithDeadline("an interrupted rype_filter", [&]() {
			auto stream = SubmitStream(f.con, f.filter_sql);
			REQUIRE_FALSE(stream->HasError());
			REQUIRE(stream->Fetch());
			f.con.Interrupt();
			// The rest of the result may already be buffered; past it the query ends, by error or by
			// finishing, and does not hang.
			while (stream->Fetch()) {
			}
		});
	}
	SECTION("dropped") {
		WithDeadline("a dropped rype_filter result", [&]() {
			auto stream = SubmitStream(f.con, f.filter_sql);
			REQUIRE_FALSE(stream->HasError());
			REQUIRE(stream->Fetch());
			stream.reset();
		});
	}

	WithDeadline("the query after it", [&]() { REQUIRE(f.Ords(f.filter_sql) == expected); });
}
