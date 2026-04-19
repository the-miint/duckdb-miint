#include <catch2/catch_test_macros.hpp>

#include <rocksdb/version.h>
#include <smr_api.h>
#include <cstring>
#include <string>

TEST_CASE("vcpkg provides rocksdb headers", "[sortmerna][build]") {
	REQUIRE(ROCKSDB_MAJOR >= 7);
}

TEST_CASE("sortmerna static lib is linkable", "[sortmerna][build]") {
	const char *v = smr_version();
	REQUIRE(v != nullptr);
	REQUIRE(std::strcmp(v, "4.4.0") == 0);
}

TEST_CASE("sortmerna config ABI matches header", "[sortmerna][build]") {
	smr_config_t cfg;
	smr_config_init(&cfg);
	REQUIRE(cfg.struct_size == static_cast<uint32_t>(sizeof(smr_config_t)));
}

TEST_CASE("MIINT_HAS_SORTMERNA is defined when enabled", "[sortmerna][build]") {
#ifndef MIINT_HAS_SORTMERNA
	FAIL("MIINT_HAS_SORTMERNA not defined");
#endif
	REQUIRE(true);
}
