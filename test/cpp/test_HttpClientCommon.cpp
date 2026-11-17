// What the HTTP clients share (src/include/http_client_common.hpp): the URL settings,
// to which they append each request's path, and the wait that Ctrl-C ends.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "http_client_common.hpp"
#include "ncbi_client.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <chrono>
#include <string>

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

// NCBI's quota is per IP address without an API key and per key with one, and it does not
// reset between queries. Each statement builds its own NCBIClient, so pacing kept per client
// let back-to-back statements exceed the quota. NCBI then answers 429. Natively that is
// retried, but in a browser the 429 carries no CORS header, so duckdb-wasm sees a failed
// request and the statement fails (the-miint/duckdb-miint#459).
// Each test sends to its own server (".invalid" hosts, never contacted), so no other test's
// requests count against its quota. Without a key the quota is 3 requests a second, and requests
// go PACING_MARGIN further apart than that.
const auto NO_KEY_INTERVAL = std::chrono::milliseconds(
    static_cast<int>(1000.0 * miint::NCBIClient::PACING_MARGIN / miint::NCBIClient::RATE_LIMIT_NO_KEY));

TEST_CASE("NCBIClient paces requests across clients that share a server and API key", "[http_client]") {
	DuckDB db(nullptr);
	Connection con1(db);
	Connection con2(db);
	const std::string url = "http://pacing-shared.invalid/entrez/eutils/efetch.fcgi";
	const auto start = std::chrono::steady_clock::now();
	miint::NCBIClient(*con1.context).RespectRateLimit(url);
	miint::NCBIClient(*con2.context).RespectRateLimit(url);
	miint::NCBIClient(*con1.context).RespectRateLimit(url);
	// The first request may go at once and each later one an interval after the one before, so
	// the third goes no sooner than two intervals after the start. A loaded machine can only
	// make it later.
	CHECK(std::chrono::steady_clock::now() - start >= 2 * NO_KEY_INTERVAL);
}

TEST_CASE("NCBIClient paces each server and API key on its own", "[http_client]") {
	// Two servers, or two keys, are two quotas: requests to one never wait for the other's.
	DuckDB db(nullptr);
	Connection con(db);
	const auto start = std::chrono::steady_clock::now();
	miint::NCBIClient(*con.context).RespectRateLimit("http://pacing-a.invalid/entrez/eutils/epost.fcgi");
	miint::NCBIClient(*con.context).RespectRateLimit("http://pacing-b.invalid/entrez/eutils/epost.fcgi");
	miint::NCBIClient(*con.context, "pacing-key").RespectRateLimit("http://pacing-a.invalid/entrez/eutils/efetch.fcgi");
	// Each is the first request on its pacer, so none waits.
	CHECK(std::chrono::steady_clock::now() - start < NO_KEY_INTERVAL);
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
