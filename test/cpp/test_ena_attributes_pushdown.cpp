#include <catch2/catch_test_macros.hpp>

#include "ena_attributes_filter.hpp"
#include "ena_parser.hpp"

#include <string>
#include <vector>

using miint::BuildStructuredSearchURL;
using miint::ENAParser;
using miint::ENATSVResult;
using miint::SampleAttribute;
using miint::UnpivotStructuredTSV;

TEST_CASE("BuildStructuredSearchURL emits tag-only query when no pinned values", "[ena][pushdown]") {
	std::vector<std::string> sample_accs = {"SAMN111", "SAMN222"};
	std::vector<std::string> tags = {"host_body_site", "collection_date"};
	std::vector<std::pair<std::string, std::string>> pairs;

	auto url = BuildStructuredSearchURL(sample_accs, tags, pairs);

	CHECK(url.find("https://www.ebi.ac.uk/ena/portal/api/search?") == 0);
	CHECK(url.find("result=sample") != std::string::npos);
	CHECK(url.find("fields=sample_accession,host_body_site,collection_date") != std::string::npos);
	// Accession IN-filter, URL-encoded.
	CHECK(url.find("sample_accession%20IN%20%28%22SAMN111%22%2C%22SAMN222%22%29") != std::string::npos);
	// No AND clauses — no tag_value_pairs were supplied.
	CHECK(url.find("%20AND%20") == std::string::npos);
	CHECK(url.find("format=tsv") != std::string::npos);
}

TEST_CASE("BuildStructuredSearchURL appends AND tag=value for each pinned pair", "[ena][pushdown]") {
	std::vector<std::string> sample_accs = {"SAMN111"};
	std::vector<std::string> tags = {"host_body_site"};
	std::vector<std::pair<std::string, std::string>> pairs = {{"host_body_site", "UBERON:feces"}};

	auto url = BuildStructuredSearchURL(sample_accs, tags, pairs);

	// The `:` in "UBERON:feces" must be percent-encoded as %3A.
	CHECK(url.find("%20AND%20host_body_site%3D%22UBERON%3Afeces%22") != std::string::npos);
	CHECK(url.find("fields=sample_accession,host_body_site") != std::string::npos);
}

TEST_CASE("BuildStructuredSearchURL deduplicates fields while preserving order", "[ena][pushdown]") {
	// Caller may pass the same tag twice (e.g., from a weirdly-canonicalized
	// filter). `sample_accession` must stay first; duplicates drop.
	std::vector<std::string> sample_accs = {"SAMN111"};
	std::vector<std::string> tags = {"host_body_site", "host_body_site"};
	std::vector<std::pair<std::string, std::string>> pairs;

	auto url = BuildStructuredSearchURL(sample_accs, tags, pairs);

	CHECK(url.find("fields=sample_accession,host_body_site") != std::string::npos);
	// No double host_body_site.
	CHECK(url.find("host_body_site,host_body_site") == std::string::npos);
}

TEST_CASE("BuildStructuredSearchURL rejects empty sample_accessions", "[ena][pushdown]") {
	std::vector<std::string> tags = {"host_body_site"};
	std::vector<std::pair<std::string, std::string>> pairs;

	CHECK_THROWS(BuildStructuredSearchURL({}, tags, pairs));
}

TEST_CASE("BuildStructuredSearchURL rejects empty tags", "[ena][pushdown]") {
	std::vector<std::string> sample_accs = {"SAMN111"};
	std::vector<std::pair<std::string, std::string>> pairs;

	CHECK_THROWS(BuildStructuredSearchURL(sample_accs, {}, pairs));
}

TEST_CASE("BuildStructuredSearchURL rejects pinned tag not in tags", "[ena][pushdown]") {
	// Callers that hand-build a pushdown (rather than going through
	// ExtractPushdownPredicates) could supply a tag_value_pair whose tag is
	// absent from `tags`. The resulting URL would carry a query clause for a
	// field we didn't ask ENA to return. Reject to prevent this mismatch.
	std::vector<std::string> sample_accs = {"SAMN111"};
	std::vector<std::string> tags = {"host_body_site"};
	std::vector<std::pair<std::string, std::string>> pairs = {{"collection_date", "2020-03-14"}};

	CHECK_THROWS(BuildStructuredSearchURL(sample_accs, tags, pairs));
}

TEST_CASE("BuildStructuredSearchURL percent-encodes control chars in values", "[ena][pushdown]") {
	// Values come from user-supplied SQL literals; a malformed or hostile
	// value shouldn't be able to break out of the quoted `field="..."`
	// context. Exercise CR / LF / NUL / tab / quote / ampersand / space
	// specifically.
	std::vector<std::string> sample_accs = {"SAMN111"};
	std::vector<std::string> tags = {"host_body_site"};

	std::string value_with_controls;
	value_with_controls.push_back('\r');
	value_with_controls.push_back('\n');
	value_with_controls.push_back('\t');
	value_with_controls.push_back('\0');
	value_with_controls.push_back(' ');
	value_with_controls.push_back('&');
	value_with_controls.push_back('"');
	value_with_controls.push_back('a');

	std::vector<std::pair<std::string, std::string>> pairs = {{"host_body_site", value_with_controls}};

	auto url = BuildStructuredSearchURL(sample_accs, tags, pairs);

	// Each control char must appear as its upper-case hex %XX escape.
	CHECK(url.find("%0D") != std::string::npos); // CR
	CHECK(url.find("%0A") != std::string::npos); // LF
	CHECK(url.find("%09") != std::string::npos); // TAB
	CHECK(url.find("%00") != std::string::npos); // NUL
	CHECK(url.find("%20") != std::string::npos); // space
	CHECK(url.find("%26") != std::string::npos); // &
	CHECK(url.find("%22") != std::string::npos); // " (also closes the field="" pair, so appears 2x)
	// Literal plaintext char stays unescaped.
	CHECK(url.find("a%22") != std::string::npos);
	// No raw CR/LF/NUL survived.
	CHECK(url.find('\r') == std::string::npos);
	CHECK(url.find('\n') == std::string::npos);
	CHECK(url.find('\0') == std::string::npos);
}

TEST_CASE("BuildStructuredSearchURL rejects injection via accession", "[ena][pushdown]") {
	// Accession validation is delegated to ENAParser::ValidateAccession,
	// which rejects non-alphanumeric (modulo _ - .). Confirm upstream
	// validation fires.
	std::vector<std::string> sample_accs = {"SAMN\"; DROP TABLE"};
	std::vector<std::string> tags = {"host_body_site"};
	std::vector<std::pair<std::string, std::string>> pairs;

	CHECK_THROWS(BuildStructuredSearchURL(sample_accs, tags, pairs));
}

TEST_CASE("UnpivotStructuredTSV emits one row per non-empty tag column", "[ena][pushdown]") {
	ENATSVResult parsed;
	parsed.column_names = {"sample_accession", "host_body_site", "collection_date"};
	parsed.rows = {
	    {"SAMN111", "UBERON:feces", "2020-03-14"},
	    {"SAMN222", "UBERON:saliva", "2021-06-02"},
	};

	auto rows = UnpivotStructuredTSV(parsed, {"host_body_site", "collection_date"});

	REQUIRE(rows.size() == 4);
	CHECK(rows[0].sample_accession == "SAMN111");
	CHECK(rows[0].tag == "host_body_site");
	CHECK(rows[0].value == "UBERON:feces");
	CHECK(rows[1].sample_accession == "SAMN111");
	CHECK(rows[1].tag == "collection_date");
	CHECK(rows[1].value == "2020-03-14");
	CHECK(rows[2].sample_accession == "SAMN222");
	CHECK(rows[2].tag == "host_body_site");
	CHECK(rows[2].value == "UBERON:saliva");
	CHECK(rows[3].sample_accession == "SAMN222");
	CHECK(rows[3].tag == "collection_date");
	CHECK(rows[3].value == "2021-06-02");
}

TEST_CASE("UnpivotStructuredTSV skips empty cells", "[ena][pushdown]") {
	// Matches XML-path behavior: samples without a given attribute produce no
	// row for that (sample, tag) pair.
	ENATSVResult parsed;
	parsed.column_names = {"sample_accession", "host_body_site", "collection_date"};
	parsed.rows = {
	    {"SAMN111", "", "2020-03-14"},    // missing host_body_site
	    {"SAMN222", "UBERON:saliva", ""}, // missing collection_date
	};

	auto rows = UnpivotStructuredTSV(parsed, {"host_body_site", "collection_date"});

	REQUIRE(rows.size() == 2);
	CHECK(rows[0].sample_accession == "SAMN111");
	CHECK(rows[0].tag == "collection_date");
	CHECK(rows[1].sample_accession == "SAMN222");
	CHECK(rows[1].tag == "host_body_site");
}

TEST_CASE("UnpivotStructuredTSV handles short rows (missing trailing columns)", "[ena][pushdown]") {
	// ENA sometimes returns rows shorter than the header row when all trailing
	// columns are empty. Treat out-of-bounds cells as empty.
	ENATSVResult parsed;
	parsed.column_names = {"sample_accession", "host_body_site", "collection_date"};
	parsed.rows = {
	    {"SAMN111", "UBERON:feces"}, // collection_date column absent in row
	};

	auto rows = UnpivotStructuredTSV(parsed, {"host_body_site", "collection_date"});

	REQUIRE(rows.size() == 1);
	CHECK(rows[0].tag == "host_body_site");
}

TEST_CASE("UnpivotStructuredTSV returns empty when sample_accession column is missing", "[ena][pushdown]") {
	ENATSVResult parsed;
	parsed.column_names = {"host_body_site"};
	parsed.rows = {{"UBERON:feces"}};

	auto rows = UnpivotStructuredTSV(parsed, {"host_body_site"});

	CHECK(rows.empty());
}

TEST_CASE("UnpivotStructuredTSV case-insensitive column matching", "[ena][pushdown]") {
	ENATSVResult parsed;
	parsed.column_names = {"SAMPLE_ACCESSION", "HOST_BODY_SITE"};
	parsed.rows = {{"SAMN111", "UBERON:feces"}};

	auto rows = UnpivotStructuredTSV(parsed, {"host_body_site"});

	REQUIRE(rows.size() == 1);
	CHECK(rows[0].sample_accession == "SAMN111");
	CHECK(rows[0].tag == "host_body_site");
	CHECK(rows[0].value == "UBERON:feces");
}
