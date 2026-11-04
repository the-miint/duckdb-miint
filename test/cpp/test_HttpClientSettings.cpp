// The NCBI and BLAST clients take their URLs from settings, which a user can point at a
// mirror and the tests point at the mock (test/sql/http_client_waits_mock.test). The
// clients append each request's path to the URL, so a trailing slash would double up,
// and an empty one would send a relative URL that fails without naming the setting.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "ncbi_client.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

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
		Exec(con, "SET miint_ncbi_eutils_url = 'http://mirror.example/entrez/eutils//'");
		CHECK(miint::NCBIClient(*con.context).EutilsBase() == "http://mirror.example/entrez/eutils");
	}

	SECTION("an empty URL is refused, naming the setting") {
		Exec(con, "SET miint_ncbi_eutils_url = ''");
		REQUIRE_THROWS_WITH(miint::NCBIClient(*con.context),
		                    ContainsSubstring("miint_ncbi_eutils_url cannot be empty"));
	}
}
