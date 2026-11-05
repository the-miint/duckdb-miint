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
