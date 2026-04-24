#include <catch2/catch_test_macros.hpp>

#include "ena_search_field_registry.hpp"

TEST_CASE("ENASearchFieldRegistry recognizes host_body_site", "[ena][registry]") {
	CHECK(miint::ENASearchFieldRegistry::IsSearchableSampleField("host_body_site"));
	CHECK(miint::ENASearchFieldRegistry::IsSearchableSampleField("collection_date"));
	CHECK(miint::ENASearchFieldRegistry::IsSearchableSampleField("scientific_name"));
	CHECK(miint::ENASearchFieldRegistry::IsSearchableSampleField("sample_accession"));
}

TEST_CASE("ENASearchFieldRegistry rejects unknown tag", "[ena][registry]") {
	CHECK_FALSE(miint::ENASearchFieldRegistry::IsSearchableSampleField("some_custom_tag"));
	CHECK_FALSE(miint::ENASearchFieldRegistry::IsSearchableSampleField("ENV_FEATURE"));
	// Typo in a real field: still not in the set.
	CHECK_FALSE(miint::ENASearchFieldRegistry::IsSearchableSampleField("host_bodysite"));
}

TEST_CASE("ENASearchFieldRegistry is case-insensitive", "[ena][registry]") {
	// Real filters may come from users typing WHERE tag='HOST_BODY_SITE'.
	CHECK(miint::ENASearchFieldRegistry::IsSearchableSampleField("HOST_BODY_SITE"));
	CHECK(miint::ENASearchFieldRegistry::IsSearchableSampleField("Host_Body_Site"));
	CHECK(miint::ENASearchFieldRegistry::IsSearchableSampleField("COLLECTION_DATE"));
}

TEST_CASE("ENASearchFieldRegistry trims whitespace", "[ena][registry]") {
	CHECK(miint::ENASearchFieldRegistry::IsSearchableSampleField("  host_body_site  "));
	CHECK(miint::ENASearchFieldRegistry::IsSearchableSampleField("\tcollection_date\n"));
	CHECK_FALSE(miint::ENASearchFieldRegistry::IsSearchableSampleField(""));
	CHECK_FALSE(miint::ENASearchFieldRegistry::IsSearchableSampleField("   "));
}

TEST_CASE("ENASearchFieldRegistry set is non-empty and sorted", "[ena][registry]") {
	const auto &fields = miint::ENASearchFieldRegistry::SearchableSampleFields();
	// Sanity check: the curated list should have at least a few dozen entries.
	// Drops from regeneration will trip this and force an update to the test.
	REQUIRE(fields.size() >= 40);
	REQUIRE(fields.count("host_body_site") == 1);
	REQUIRE(fields.count("environment_biome") == 1);
}
