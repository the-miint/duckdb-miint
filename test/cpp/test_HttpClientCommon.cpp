// What the HTTP clients share (src/include/http_client_common.hpp): the URL settings,
// to which they append each request's path, and the wait that Ctrl-C ends. Also how
// NCBIClient spaces its requests, which needs a running DuckDB as these do.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "concurrent_test_util.hpp"
#include "http_client_common.hpp"
#include "ncbi_client.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <chrono>
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

} // namespace

TEST_CASE("NCBIClient sends E-utilities requests to miint_ncbi_eutils_url", "[http_client]") {
	DuckDB db(nullptr);
	Connection con(db);

	SECTION("NCBI's by default") {
		CHECK(miint::NCBIClient(*con.context).EutilsBase() == miint::NCBIClient::EUTILS_BASE);
	}

	SECTION("a trailing slash is dropped") {
		// A path appended after it would otherwise start "//".
		Exec(con, "SET miint_ncbi_eutils_url = 'http://mirror.example/entrez/eutils//'");
		CHECK(miint::NCBIClient(*con.context).EutilsBase() == "http://mirror.example/entrez/eutils");
	}

	SECTION("an empty URL is refused, naming the setting") {
		// Accepted, every request would be a relative URL, failing without naming the setting.
		Exec(con, "SET miint_ncbi_eutils_url = ''");
		REQUIRE_THROWS_WITH(miint::NCBIClient(*con.context),
		                    ContainsSubstring("miint_ncbi_eutils_url cannot be empty"));
	}
}

TEST_CASE("SleepUnlessInterrupted ends at Ctrl-C", "[http_client]") {
	// The SQL tests stand the deadline in for Ctrl-C; this is the Ctrl-C half of the
	// check, which Connection::Interrupt trips as the shell's Ctrl-C handler does.
	DuckDB db(nullptr);
	Connection con(db);
	con.Interrupt();
	const auto start = std::chrono::steady_clock::now();
	REQUIRE_THROWS_AS(miint::SleepUnlessInterrupted(*con.context, std::chrono::seconds(60)), InterruptException);
	// One slice, not the minute; the bound is loose so a loaded machine cannot fail it.
	CHECK(std::chrono::steady_clock::now() - start < std::chrono::seconds(10));
}

// Requests to one server with one API key are spaced as one stream, whichever client, connection
// or database sends them. Each test sends to its own ".invalid" host, never contacted, so no
// other test's requests share its pacing. The intervals are 1 s / 3 requests without a key and
// 1 s / 10 requests with one, times the 1.2 margin.
const auto NO_KEY_INTERVAL = std::chrono::milliseconds(400);
const auto KEYED_INTERVAL = std::chrono::milliseconds(120);

TEST_CASE("NCBIClient paces requests to one server and API key across clients and databases", "[http_client]") {
	DuckDB db1(nullptr);
	DuckDB db2(nullptr);
	Connection con1(db1);
	Connection con2(db2);
	// What read_ncbi sends, epost then efetch with the WebEnv it returned, then a later query's
	// efetch: three URLs on one server.
	const std::string base = "http://pacing-shared.invalid/entrez/eutils";
	auto start = std::chrono::steady_clock::now();
	miint::NCBIClient(*con1.context).RespectRateLimit(base + "/epost.fcgi");
	miint::NCBIClient(*con2.context).RespectRateLimit(base + "/efetch.fcgi?db=nuccore&query_key=1&WebEnv=MCID_a");
	miint::NCBIClient(*con1.context).RespectRateLimit(base + "/efetch.fcgi?db=taxonomy&query_key=1&WebEnv=MCID_b");
	// The first request may go at once and each later one an interval after the one before, so
	// the third goes no sooner than two intervals after the start. A loaded machine can only
	// make it later.
	CHECK(std::chrono::steady_clock::now() - start >= 2 * NO_KEY_INTERVAL);

	// With a key, the same at the keyed interval.
	start = std::chrono::steady_clock::now();
	miint::NCBIClient(*con1.context, "pacing-shared-key").RespectRateLimit(base + "/epost.fcgi");
	miint::NCBIClient(*con2.context, "pacing-shared-key").RespectRateLimit(base + "/efetch.fcgi?db=nuccore");
	CHECK(std::chrono::steady_clock::now() - start >= KEYED_INTERVAL);
}

TEST_CASE("NCBIClient lets one waiting query at a time take an opening", "[http_client]") {
	// Queries waiting together wake for the same opening; only one may send in it.
	DuckDB db(nullptr);
	constexpr size_t n_queries = 3;
	std::vector<std::unique_ptr<Connection>> cons;
	for (size_t i = 0; i < n_queries; i++) {
		cons.push_back(std::make_unique<Connection>(db));
	}
	const std::string url = "http://pacing-concurrent.invalid/entrez/eutils/efetch.fcgi";
	const auto start = std::chrono::steady_clock::now();
	miint::test::RunConcurrently(n_queries,
	                             [&](size_t i) { miint::NCBIClient(*cons[i]->context).RespectRateLimit(url); });
	// In whatever order they go, the last goes no sooner than two intervals after the start.
	CHECK(std::chrono::steady_clock::now() - start >= (n_queries - 1) * NO_KEY_INTERVAL);
}

TEST_CASE("NCBIClient paces each server and API key on its own", "[http_client]") {
	// Two servers, or two keys, are two quotas: requests to one never wait for the other's.
	// Each request below is the first on its pacing, so none waits.
	DuckDB db(nullptr);
	Connection con(db);
	auto start = std::chrono::steady_clock::now();
	miint::NCBIClient(*con.context).RespectRateLimit("http://pacing-a.invalid/entrez/eutils/epost.fcgi");
	miint::NCBIClient(*con.context).RespectRateLimit("http://pacing-b.invalid/entrez/eutils/epost.fcgi");
	CHECK(std::chrono::steady_clock::now() - start < NO_KEY_INTERVAL);

	const std::string efetch = "http://pacing-a.invalid/entrez/eutils/efetch.fcgi";
	start = std::chrono::steady_clock::now();
	miint::NCBIClient(*con.context, "pacing-key-1").RespectRateLimit(efetch);
	miint::NCBIClient(*con.context, "pacing-key-2").RespectRateLimit(efetch);
	CHECK(std::chrono::steady_clock::now() - start < KEYED_INTERVAL);
}

TEST_CASE("NCBIClient leaves no wait behind when Ctrl-C ends one", "[http_client]") {
	// A query cancelled while it waits for its turn sends nothing, so the next request must not
	// wait for it as well.
	DuckDB db(nullptr);
	Connection con1(db);
	Connection con2(db);
	const std::string url = "http://pacing-interrupted.invalid/entrez/eutils/efetch.fcgi";
	const auto start = std::chrono::steady_clock::now();
	miint::NCBIClient(*con1.context).RespectRateLimit(url);
	con2.Interrupt();
	REQUIRE_THROWS_AS(miint::NCBIClient(*con2.context).RespectRateLimit(url), InterruptException);
	miint::NCBIClient(*con1.context).RespectRateLimit(url);
	// con1's second request goes one interval after its first, as if con2 had never asked.
	CHECK(std::chrono::steady_clock::now() - start < 2 * NO_KEY_INTERVAL);
}

// read_ncbi_taxdump fetches an http(s) source with an HTTP GET and opens any other through
// DuckDB's file systems. The scheme is DuckDB's parse of it, and schemes ignore case
// (RFC 3986, section 3.1).
TEST_CASE("IsHttpUrl tells http(s) URLs by their scheme", "[http_client]") {
	for (const char *url : {"http://h/x", "https://h/taxdump.tar.gz", "HTTPS://h/x", "Http://h:8080/a/b?q=1"}) {
		INFO(url);
		CHECK(miint::IsHttpUrl(url));
	}
	for (const char *other : {"s3://bucket/x", "file:///a/b", "/data/http://x", "data/taxdump", "httpx://h/x",
	                          "http:/h/x", "C:\\taxdump", ""}) {
		INFO(other);
		CHECK_FALSE(miint::IsHttpUrl(other));
	}
}
