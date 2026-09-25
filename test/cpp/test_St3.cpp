#ifdef MIINT_HAS_ST3

#include <catch2/catch_test_macros.hpp>

// st3.h refers to the Arrow C Data Interface structs by pointer and deliberately
// does not declare them; DuckDB ships the same struct names in this header.
#include "duckdb/common/arrow/arrow.hpp"

#include "st3.h"

// The one thing this test proves is that the st3 C ABI is compiled into the
// Rust umbrella archive and reachable from a miint translation unit. The
// behavioural coverage lives in the SQL tests, which run the real pipeline
// against st3's committed fixtures.
TEST_CASE("st3 C ABI links and reports the v1 ABI version", "[st3]") {
	REQUIRE(st3_abi_version() == 0);
}

#endif // MIINT_HAS_ST3
