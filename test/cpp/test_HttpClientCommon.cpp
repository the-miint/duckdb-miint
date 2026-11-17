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
TEST_CASE("NCBIClient paces requests across clients that share an API key", "[http_client]") {
	DuckDB db(nullptr);
	Connection con1(db);
	Connection con2(db);
	// A key no other test uses, so no earlier request counts against it. With a key the quota
	// is 10 requests a second: one every 100 ms.
	const std::string key = "test-http-client-common-shared-pacer";
	const auto interval = std::chrono::milliseconds(static_cast<int>(1000.0 / miint::NCBIClient::RATE_LIMIT_WITH_KEY));
	const auto start = std::chrono::steady_clock::now();
	miint::NCBIClient(*con1.context, key).RespectRateLimit();
	miint::NCBIClient(*con2.context, key).RespectRateLimit();
	miint::NCBIClient(*con1.context, key).RespectRateLimit();
	// The first request may go at once and each later one an interval after the one before, so
	// the third goes no sooner than two intervals after the start. A loaded machine can only
	// make it later.
	CHECK(std::chrono::steady_clock::now() - start >= 2 * interval);
}

TEST_CASE("NCBIClient paces each API key on its own", "[http_client]") {
	// Two keys are two quotas: one key's requests never wait for another's.
	DuckDB db(nullptr);
	Connection con(db);
	const auto interval = std::chrono::milliseconds(static_cast<int>(1000.0 / miint::NCBIClient::RATE_LIMIT_WITH_KEY));
	const auto start = std::chrono::steady_clock::now();
	miint::NCBIClient(*con.context, "test-http-client-common-key-a").RespectRateLimit();
	miint::NCBIClient(*con.context, "test-http-client-common-key-b").RespectRateLimit();
	CHECK(std::chrono::steady_clock::now() - start < interval);
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
