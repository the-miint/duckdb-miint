#include <catch2/catch_test_macros.hpp>

#include "ena_attributes_filter.hpp"

#include <memory>
#include <string>
#include <utility>
#include <vector>

using miint::ENAAttributePushdown;
using miint::ENAFilterNode;
using miint::ExtractPushdownPredicates;

// Helper: wrap a single ENAFilterNode as the conjunct-list that
// ExtractPushdownPredicates expects (DuckDB already splits top-level AND).
static std::vector<std::unique_ptr<ENAFilterNode>> SingleConjunct(std::unique_ptr<ENAFilterNode> node) {
	std::vector<std::unique_ptr<ENAFilterNode>> out;
	out.push_back(std::move(node));
	return out;
}

TEST_CASE("ExtractPushdownPredicates detects tag = 'X'", "[ena][pushdown]") {
	auto conjuncts = SingleConjunct(ENAFilterNode::MakeEqual("tag", "host_body_site"));

	auto result = ExtractPushdownPredicates(conjuncts);

	REQUIRE(result.tags.size() == 1);
	CHECK(result.tags[0] == "host_body_site");
	CHECK(result.tag_value_pairs.empty());
}

TEST_CASE("ExtractPushdownPredicates detects tag IN (...) when all in allowlist", "[ena][pushdown]") {
	auto conjuncts = SingleConjunct(ENAFilterNode::MakeIn("tag", {"host_body_site", "collection_date"}));

	auto result = ExtractPushdownPredicates(conjuncts);

	REQUIRE(result.tags.size() == 2);
	CHECK(result.tags[0] == "host_body_site");
	CHECK(result.tags[1] == "collection_date");
	CHECK(result.tag_value_pairs.empty());
}

TEST_CASE("ExtractPushdownPredicates detects tag = X AND value = Y", "[ena][pushdown]") {
	// DuckDB splits top-level AND into two conjuncts.
	std::vector<std::unique_ptr<ENAFilterNode>> conjuncts;
	conjuncts.push_back(ENAFilterNode::MakeEqual("tag", "host_body_site"));
	conjuncts.push_back(ENAFilterNode::MakeEqual("value", "UBERON:feces"));

	auto result = ExtractPushdownPredicates(conjuncts);

	REQUIRE(result.tags.size() == 1);
	CHECK(result.tags[0] == "host_body_site");
	REQUIRE(result.tag_value_pairs.size() == 1);
	CHECK(result.tag_value_pairs[0].first == "host_body_site");
	CHECK(result.tag_value_pairs[0].second == "UBERON:feces");
}

TEST_CASE("ExtractPushdownPredicates falls back when any tag is unknown (IN list)", "[ena][pushdown]") {
	auto conjuncts = SingleConjunct(ENAFilterNode::MakeIn("tag", {"host_body_site", "some_custom_tag"}));

	auto result = ExtractPushdownPredicates(conjuncts);

	CHECK(result.tags.empty());
	CHECK(result.tag_value_pairs.empty());
}

TEST_CASE("ExtractPushdownPredicates falls back on OR across valid+unknown tags", "[ena][pushdown]") {
	std::vector<std::unique_ptr<ENAFilterNode>> or_children;
	or_children.push_back(ENAFilterNode::MakeEqual("tag", "host_body_site"));
	or_children.push_back(ENAFilterNode::MakeEqual("tag", "some_custom_tag"));

	auto conjuncts = SingleConjunct(ENAFilterNode::MakeOr(std::move(or_children)));

	auto result = ExtractPushdownPredicates(conjuncts);

	CHECK(result.tags.empty());
	CHECK(result.tag_value_pairs.empty());
}

TEST_CASE("ExtractPushdownPredicates falls back on non-equality operator on tag", "[ena][pushdown]") {
	// tag LIKE 'host_%' — ENAFilterNode::MakeUnsupported() represents an
	// expression the DuckDB→AST translator couldn't map onto EQ/IN/AND/OR.
	auto conjuncts = SingleConjunct(ENAFilterNode::MakeUnsupported());

	auto result = ExtractPushdownPredicates(conjuncts);

	CHECK(result.tags.empty());
	CHECK(result.tag_value_pairs.empty());
}

TEST_CASE("ExtractPushdownPredicates falls back on non-equality operator on value", "[ena][pushdown]") {
	// tag='host_body_site' AND value LIKE 'UBERON:%' — the tag conjunct maps
	// cleanly, but the value conjunct is unsupported. Any unsupported conjunct
	// forces fallback (not a partial push).
	std::vector<std::unique_ptr<ENAFilterNode>> conjuncts;
	conjuncts.push_back(ENAFilterNode::MakeEqual("tag", "host_body_site"));
	conjuncts.push_back(ENAFilterNode::MakeUnsupported());

	auto result = ExtractPushdownPredicates(conjuncts);

	CHECK(result.tags.empty());
	CHECK(result.tag_value_pairs.empty());
}

TEST_CASE("ExtractPushdownPredicates falls back when value is constrained without a pinned tag", "[ena][pushdown]") {
	// value='UBERON:feces' alone is ambiguous (which tag?). Fallback required.
	auto conjuncts = SingleConjunct(ENAFilterNode::MakeEqual("value", "UBERON:feces"));

	auto result = ExtractPushdownPredicates(conjuncts);

	CHECK(result.tags.empty());
	CHECK(result.tag_value_pairs.empty());
}

TEST_CASE("ExtractPushdownPredicates falls back when sample_accession is constrained", "[ena][pushdown]") {
	// tag='host_body_site' AND sample_accession='SAMN...' — we can't translate
	// the sample_accession predicate through the /search endpoint with the same
	// semantics as the XML path's output, so we reject the whole thing for
	// safety.
	std::vector<std::unique_ptr<ENAFilterNode>> conjuncts;
	conjuncts.push_back(ENAFilterNode::MakeEqual("tag", "host_body_site"));
	conjuncts.push_back(ENAFilterNode::MakeEqual("sample_accession", "SAMN12345"));

	auto result = ExtractPushdownPredicates(conjuncts);

	CHECK(result.tags.empty());
	CHECK(result.tag_value_pairs.empty());
}

TEST_CASE("ExtractPushdownPredicates is case-insensitive for tag values", "[ena][pushdown]") {
	// WHERE tag='HOST_BODY_SITE' — user-typed uppercase must still match the
	// registry and be stored lowercase in the output.
	auto conjuncts = SingleConjunct(ENAFilterNode::MakeEqual("tag", "HOST_BODY_SITE"));

	auto result = ExtractPushdownPredicates(conjuncts);

	REQUIRE(result.tags.size() == 1);
	CHECK(result.tags[0] == "host_body_site");
}

TEST_CASE("ExtractPushdownPredicates deduplicates repeated tags", "[ena][pushdown]") {
	// WHERE tag IN ('host_body_site','HOST_BODY_SITE') — two spellings of the
	// same canonical tag should collapse.
	auto conjuncts = SingleConjunct(ENAFilterNode::MakeIn("tag", {"host_body_site", "HOST_BODY_SITE"}));

	auto result = ExtractPushdownPredicates(conjuncts);

	REQUIRE(result.tags.size() == 1);
	CHECK(result.tags[0] == "host_body_site");
}

TEST_CASE("ExtractPushdownPredicates falls back on multiple tag IN conjuncts", "[ena][pushdown]") {
	// `tag IN ('a','b') AND tag IN ('c','d')` has SQL INTERSECTION semantics
	// (tag must be in both sets). Emitting the union to ENA would broaden the
	// wire query; reject for strict correctness.
	std::vector<std::unique_ptr<ENAFilterNode>> conjuncts;
	conjuncts.push_back(ENAFilterNode::MakeIn("tag", {"host_body_site", "collection_date"}));
	conjuncts.push_back(ENAFilterNode::MakeIn("tag", {"collection_date", "scientific_name"}));

	auto result = ExtractPushdownPredicates(conjuncts);

	CHECK(result.tags.empty());
	CHECK(result.tag_value_pairs.empty());
}

TEST_CASE("ExtractPushdownPredicates falls back on tag IN (...) AND value=Y", "[ena][pushdown]") {
	// value= pin requires a single-tag EQ pin (one row → one tag → unambiguous
	// column to constrain). With tag IN, the column to constrain is ambiguous.
	std::vector<std::unique_ptr<ENAFilterNode>> conjuncts;
	conjuncts.push_back(ENAFilterNode::MakeIn("tag", {"host_body_site", "collection_date"}));
	conjuncts.push_back(ENAFilterNode::MakeEqual("value", "UBERON:feces"));

	auto result = ExtractPushdownPredicates(conjuncts);

	CHECK(result.tags.empty());
	CHECK(result.tag_value_pairs.empty());
}

TEST_CASE("ExtractPushdownPredicates on empty input returns empty (no pushdown)", "[ena][pushdown]") {
	std::vector<std::unique_ptr<ENAFilterNode>> conjuncts;

	auto result = ExtractPushdownPredicates(conjuncts);

	CHECK(result.tags.empty());
	CHECK(result.tag_value_pairs.empty());
}
