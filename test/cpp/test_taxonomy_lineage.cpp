#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "taxonomy_lineage.hpp"

#include <stdexcept>

using namespace miint;
using Catch::Matchers::ContainsSubstring;

namespace {
// A two-taxon taxonomy efetch response. 9606 uses the legacy "superkingdom" rank
// for its domain; 83333 uses the newer "domain" rank AND is a strain, exercising
// both the superkingdom/domain alias and the strain column. The nested <Taxon>
// elements inside <LineageEx> must NOT be mistaken for top-level query taxa.
const char *kTaxonomyXml = R"(<?xml version="1.0"?>
<TaxaSet>
<Taxon>
<TaxId>9606</TaxId>
<ScientificName>Homo sapiens</ScientificName>
<ParentTaxId>9605</ParentTaxId>
<Rank>species</Rank>
<LineageEx>
<Taxon><TaxId>131567</TaxId><ScientificName>cellular organisms</ScientificName><Rank>no rank</Rank></Taxon>
<Taxon><TaxId>2759</TaxId><ScientificName>Eukaryota</ScientificName><Rank>superkingdom</Rank></Taxon>
<Taxon><TaxId>33208</TaxId><ScientificName>Metazoa</ScientificName><Rank>kingdom</Rank></Taxon>
<Taxon><TaxId>7711</TaxId><ScientificName>Chordata</ScientificName><Rank>phylum</Rank></Taxon>
<Taxon><TaxId>40674</TaxId><ScientificName>Mammalia</ScientificName><Rank>class</Rank></Taxon>
<Taxon><TaxId>9443</TaxId><ScientificName>Primates</ScientificName><Rank>order</Rank></Taxon>
<Taxon><TaxId>9604</TaxId><ScientificName>Hominidae</ScientificName><Rank>family</Rank></Taxon>
<Taxon><TaxId>9605</TaxId><ScientificName>Homo</ScientificName><Rank>genus</Rank></Taxon>
</LineageEx>
</Taxon>
<Taxon>
<TaxId>83333</TaxId>
<ScientificName>Escherichia coli K-12</ScientificName>
<ParentTaxId>562</ParentTaxId>
<Rank>strain</Rank>
<LineageEx>
<Taxon><TaxId>2</TaxId><ScientificName>Bacteria</ScientificName><Rank>domain</Rank></Taxon>
<Taxon><TaxId>1224</TaxId><ScientificName>Pseudomonadota</ScientificName><Rank>phylum</Rank></Taxon>
<Taxon><TaxId>561</TaxId><ScientificName>Escherichia</ScientificName><Rank>genus</Rank></Taxon>
<Taxon><TaxId>562</TaxId><ScientificName>Escherichia coli</ScientificName><Rank>species</Rank></Taxon>
</LineageEx>
</Taxon>
</TaxaSet>
)";
} // namespace

TEST_CASE("RankColumn maps NCBI ranks to collapse columns", "[taxdump][lineage]") {
	// The superkingdom/domain alias is the whole reason this is centralized: NCBI is
	// mid-transition and both strings must land in the same column.
	CHECK(RankColumn("superkingdom") == "domain");
	CHECK(RankColumn("domain") == "domain");
	CHECK(RankColumn("phylum") == "phylum");
	CHECK(RankColumn("species") == "species");
	CHECK(RankColumn("strain") == "strain");
	CHECK(RankColumn("clade").empty()); // non-standard ranks are dropped
	CHECK(RankColumn("subspecies").empty());
	CHECK(RankColumn("no rank").empty());
}

// A queried taxid that NCBI has merged resolves to a taxon whose TaxId is the
// current one; the old id(s) appear under <AkaTaxIds>. Capturing them lets
// read_ncbi_lineage recognise the requested id as resolved instead of "omitted".
TEST_CASE("ParseTaxonomyXML captures AkaTaxIds for merged taxids", "[taxdump][lineage]") {
	const char *xml = R"(<TaxaSet>
<Taxon>
<TaxId>67890</TaxId>
<ScientificName>Example taxon</ScientificName>
<ParentTaxId>2</ParentTaxId>
<Rank>species</Rank>
<AkaTaxIds>
<TaxId>12345</TaxId>
<TaxId>12346</TaxId>
</AkaTaxIds>
<LineageEx>
<Taxon><TaxId>2</TaxId><ScientificName>Bacteria</ScientificName><Rank>domain</Rank></Taxon>
</LineageEx>
</Taxon>
</TaxaSet>)";
	auto lineages = ParseTaxonomyXML(xml);
	REQUIRE(lineages.size() == 1);
	CHECK(lineages[0].taxid == 67890); // the FIRST <TaxId> is the current taxon, not an aka id
	CHECK(lineages[0].domain == "Bacteria");
	REQUIRE(lineages[0].aka_taxids.size() == 2);
	CHECK(lineages[0].aka_taxids[0] == 12345);
	CHECK(lineages[0].aka_taxids[1] == 12346);
}

// PARITY PIN: the online path (read_ncbi_lineage -> ParseTaxonomyXML) and the offline
// path (the taxonomy_lineage SQL macro) MUST emit identical values for the same lineage
// -- that is the whole point of the shared schema. This XML reproduces the synthetic
// taxdump tree's Escherichia coli K-12 lineage (taxid 90), and the expected values below
// are byte-for-byte the golden asserted for taxid 90 in test/sql/taxonomy_lineage.test.
// If the C++ rank-collapse (RankColumn/CollapseLineage/FormatLineage) or the SQL macro
// drifts, one of the two tests breaks -- keep them in lockstep.
TEST_CASE("ParseTaxonomyXML matches the taxonomy_lineage macro golden (online==offline)", "[taxdump][lineage]") {
	const char *xml = R"(<TaxaSet>
<Taxon>
<TaxId>90</TaxId>
<ScientificName>Escherichia coli K-12</ScientificName>
<ParentTaxId>80</ParentTaxId>
<Rank>strain</Rank>
<LineageEx>
<Taxon><TaxId>1</TaxId><ScientificName>root</ScientificName><Rank>no rank</Rank></Taxon>
<Taxon><TaxId>10</TaxId><ScientificName>cellular organisms</ScientificName><Rank>no rank</Rank></Taxon>
<Taxon><TaxId>20</TaxId><ScientificName>Bacteria</ScientificName><Rank>superkingdom</Rank></Taxon>
<Taxon><TaxId>30</TaxId><ScientificName>Pseudomonadota</ScientificName><Rank>phylum</Rank></Taxon>
<Taxon><TaxId>40</TaxId><ScientificName>Gammaproteobacteria</ScientificName><Rank>class</Rank></Taxon>
<Taxon><TaxId>50</TaxId><ScientificName>Enterobacterales</ScientificName><Rank>order</Rank></Taxon>
<Taxon><TaxId>60</TaxId><ScientificName>Enterobacteriaceae</ScientificName><Rank>family</Rank></Taxon>
<Taxon><TaxId>70</TaxId><ScientificName>Escherichia</ScientificName><Rank>genus</Rank></Taxon>
<Taxon><TaxId>80</TaxId><ScientificName>Escherichia coli</ScientificName><Rank>species</Rank></Taxon>
</LineageEx>
</Taxon>
</TaxaSet>)";
	auto lineages = ParseTaxonomyXML(xml);
	REQUIRE(lineages.size() == 1);
	const auto &l = lineages[0];
	CHECK(l.taxid == 90);
	CHECK(l.name == "Escherichia coli K-12");
	CHECK(l.rank == "strain");
	CHECK(l.domain == "Bacteria"); // 'superkingdom' -> domain, as in the synthetic tree
	CHECK(l.phylum == "Pseudomonadota");
	CHECK(l.tax_class == "Gammaproteobacteria");
	CHECK(l.tax_order == "Enterobacterales");
	CHECK(l.family == "Enterobacteriaceae");
	CHECK(l.genus == "Escherichia");
	CHECK(l.species == "Escherichia coli");
	CHECK(l.strain == "Escherichia coli K-12");
	CHECK(l.lineage == "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__"
	                   "Enterobacteriaceae;g__Escherichia;s__Escherichia coli;t__Escherichia coli K-12");
}

TEST_CASE("ParseTaxonomyXML collapses lineages and splits top-level taxa", "[taxdump][lineage]") {
	auto lineages = ParseTaxonomyXML(kTaxonomyXml);
	REQUIRE(lineages.size() == 2); // nested <Taxon> in <LineageEx> must not add spurious rows

	SECTION("eukaryote via legacy 'superkingdom' rank, species-level, no strain") {
		const auto &h = lineages[0];
		CHECK(h.taxid == 9606);
		CHECK(h.name == "Homo sapiens");
		CHECK(h.rank == "species");
		CHECK(h.domain == "Eukaryota"); // superkingdom -> domain
		CHECK(h.phylum == "Chordata");
		CHECK(h.tax_class == "Mammalia");
		CHECK(h.tax_order == "Primates");
		CHECK(h.family == "Hominidae");
		CHECK(h.genus == "Homo");
		CHECK(h.species == "Homo sapiens"); // the query taxon itself
		CHECK(h.strain.empty());
		CHECK(h.lineage == "d__Eukaryota;p__Chordata;c__Mammalia;o__Primates;f__Hominidae;g__Homo;s__Homo sapiens");
	}

	SECTION("bacterium via 'domain' rank, strain-level, gaps left empty") {
		const auto &e = lineages[1];
		CHECK(e.taxid == 83333);
		CHECK(e.rank == "strain");
		CHECK(e.domain == "Bacteria"); // 'domain' rank string
		CHECK(e.phylum == "Pseudomonadota");
		CHECK(e.tax_class.empty()); // absent ranks stay empty
		CHECK(e.tax_order.empty());
		CHECK(e.family.empty());
		CHECK(e.genus == "Escherichia");
		CHECK(e.species == "Escherichia coli");
		CHECK(e.strain == "Escherichia coli K-12"); // the query taxon itself
		CHECK(e.lineage ==
		      "d__Bacteria;p__Pseudomonadota;c__;o__;f__;g__Escherichia;s__Escherichia coli;t__Escherichia coli K-12");
	}
}

namespace {
// The eight shared-schema columns, shallowest first; missing ranks pad with nullopt.
LineageRanks Ranks(std::initializer_list<std::optional<std::string>> values) {
	LineageRanks out;
	size_t i = 0;
	for (const auto &v : values) {
		out[i++] = v;
	}
	return out;
}

// What a shared-schema column holds for one of Lineage's string fields: NULL when empty.
std::optional<std::string> Column(const std::string &value) {
	return value.empty() ? std::nullopt : std::optional<std::string>(value);
}
} // namespace

// read_ncbi_lineage's output columns and ParseLineageString's field positions both
// come from kLineageRanks, and prefixed lineages use the letters dpcofgst. Reordering
// it changes the shared lineage schema.
TEST_CASE("kLineageRanks is the shared schema's column order and rank letters", "[lineage][parse_lineage]") {
	std::string columns;
	std::string letters;
	for (const auto &rank : kLineageRanks) {
		columns += std::string(columns.empty() ? "" : ",") + rank.column;
		letters += rank.letter;
	}
	CHECK(columns == "domain,phylum,class,order,family,genus,species,strain");
	CHECK(letters == "dpcofgst");
}

// CollapseLineage, FormatLineage and read_ncbi_lineage all reach a rank's value through
// kLineageRanks' field pointer, so a rank paired with the wrong Lineage member would file
// every value under the wrong column.
TEST_CASE("kLineageRanks pairs each column with its Lineage field", "[lineage]") {
	Lineage l;
	l.domain = "domain";
	l.phylum = "phylum";
	l.tax_class = "class";
	l.tax_order = "order";
	l.family = "family";
	l.genus = "genus";
	l.species = "species";
	l.strain = "strain";
	for (const auto &rank : kLineageRanks) {
		CHECK(l.*rank.field == rank.column);
	}
}

// tax2tree node labels and GG2 trees join ranks with "; "; miint's own lineage
// strings (FormatLineage) use ";". The converter has to read both.
TEST_CASE("ParseLineageString reads both separator styles", "[lineage][parse_lineage]") {
	auto expected = Ranks({"Bacteria", "Bacillota", "Bacilli"});
	CHECK(ParseLineageString("d__Bacteria; p__Bacillota; c__Bacilli") == expected);
	CHECK(ParseLineageString("d__Bacteria;p__Bacillota;c__Bacilli") == expected);
	// Whitespace around a field is not part of the name, including a CRLF file's '\r'.
	CHECK(ParseLineageString("  d__Bacteria ;\tp__Bacillota;c__Bacilli\r") == expected);
}

// A lineage string written by taxonomy_lineage / read_ncbi_lineage parses back into the
// columns it was written from. (A name containing ';' cannot; it errors instead.)
TEST_CASE("ParseLineageString reads back what FormatLineage writes", "[lineage][parse_lineage]") {
	auto lineages = ParseTaxonomyXML(kTaxonomyXml);
	REQUIRE(lineages.size() == 2);
	for (const auto &l : lineages) {
		INFO(l.lineage);
		auto expected = Ranks({Column(l.domain), Column(l.phylum), Column(l.tax_class), Column(l.tax_order),
		                       Column(l.family), Column(l.genus), Column(l.species), Column(l.strain)});
		CHECK(ParseLineageString(l.lineage) == expected);
	}
}

// NCBI has species names containing ';' ("Rhodococcus sp. SJ-2;3"). FormatLineage writes
// them unescaped, so the name splits into an extra field. That must be an error, never a
// lineage with the name cut short and its tail filed under the next rank.
TEST_CASE("ParseLineageString rejects a lineage whose name contained ';'", "[lineage][parse_lineage]") {
	auto lin = CollapseLineage(1, "Rhodococcus sp. SJ-2;3", "species",
	                           {{2, "Bacteria", "domain"}, {3, "Rhodococcus", "genus"}});
	REQUIRE(lin.lineage == "d__Bacteria;p__;c__;o__;f__;g__Rhodococcus;s__Rhodococcus sp. SJ-2;3");
	CHECK_THROWS_WITH(ParseLineageString(lin.lineage), ContainsSubstring("must start with 't__'"));
}

TEST_CASE("ParseLineageString: empty bodies and short lineages are missing ranks", "[lineage][parse_lineage]") {
	SECTION("all eight ranks") {
		CHECK(ParseLineageString("d__A; p__B; c__C; o__D; f__E; g__F; s__F g; t__F g H") ==
		      Ranks({"A", "B", "C", "D", "E", "F", "F g", "F g H"}));
	}
	SECTION("fewer than eight fields leaves the deeper ranks missing") {
		CHECK(ParseLineageString("d__A; p__B") == Ranks({"A", "B"}));
	}
	SECTION("an empty body is a missing rank, Greengenes style") {
		CHECK(ParseLineageString("d__A; p__B; c__; o__; f__; g__; s__") == Ranks({"A", "B"}));
	}
	SECTION("a lineage with nothing assigned is all-missing, not an error") {
		CHECK(ParseLineageString("d__; p__") == LineageRanks {});
	}
	SECTION("an internal empty body parses: gaps are not rejected here") {
		CHECK(ParseLineageString("d__A; p__; c__C") == Ranks({"A", std::nullopt, "C"}));
	}
	SECTION("a body may contain '__' after its prefix") {
		CHECK(ParseLineageString("d__A; p__B__x") == Ranks({"A", "B__x"}));
	}
}

// A field's rank comes from its position, and its letter must agree. Anything else
// means the string does not follow d/p/c/o/f/g/s/t, and guessing would misfile names.
TEST_CASE("ParseLineageString rejects a rank letter out of position", "[lineage][parse_lineage]") {
	CHECK_THROWS_AS(ParseLineageString("p__A; d__B"), std::invalid_argument);
	// A skipped rank shifts every later field.
	CHECK_THROWS_WITH(ParseLineageString("d__A; c__B"), ContainsSubstring("p__"));
	CHECK_THROWS_WITH(ParseLineageString("d__A; c__B"), ContainsSubstring("d__A; c__B"));
	CHECK_THROWS_AS(ParseLineageString("d__A; p__B; c__C; o__D; f__E; g__F; s__G; s__H"), std::invalid_argument);
	// A ninth field has no rank to go to.
	CHECK_THROWS_WITH(ParseLineageString("d__A; p__B; c__C; o__D; f__E; g__F; s__G; t__H; t__I"),
	                  ContainsSubstring("at most 8 fields"));
}

// Greengenes 13_8 writes k__ for the domain. Accepting it only when the caller says
// so means one call cannot mix d__ and k__ rows.
TEST_CASE("ParseLineageString accepts k__ only when declared", "[lineage][parse_lineage]") {
	CHECK_THROWS_WITH(ParseLineageString("k__Bacteria; p__Bacillota"), ContainsSubstring("pass domain_letter := 'k'"));
	CHECK(ParseLineageString("k__Bacteria; p__Bacillota", 'k') == Ranks({"Bacteria", "Bacillota"}));
	CHECK_THROWS_WITH(ParseLineageString("d__Bacteria; p__Bacillota", 'k'),
	                  ContainsSubstring("pass domain_letter := 'd'"));
	// k is a domain synonym only; it is never accepted below position 0.
	CHECK_THROWS_AS(ParseLineageString("k__Bacteria; k__Bacillota", 'k'), std::invalid_argument);
	CHECK_THROWS_WITH(ParseLineageString("d__Bacteria", 'x'), ContainsSubstring("domain_letter"));
}

TEST_CASE("ParseLineageString rejects malformed fields", "[lineage][parse_lineage]") {
	CHECK_THROWS_AS(ParseLineageString(""), std::invalid_argument);
	CHECK_THROWS_AS(ParseLineageString("   "), std::invalid_argument);
	// A trailing or doubled ';' is an empty field, which carries no rank letter.
	CHECK_THROWS_AS(ParseLineageString("d__A; p__B;"), std::invalid_argument);
	CHECK_THROWS_AS(ParseLineageString("d__A;; p__B"), std::invalid_argument);
	// Unprefixed strings (raw SILVA) are unsupported: a field's rank would be unknown.
	CHECK_THROWS_AS(ParseLineageString("Bacteria; Bacillota"), std::invalid_argument);
	// So are QIIME-formatted SILVA's level prefixes, which name a depth, not a rank.
	CHECK_THROWS_AS(ParseLineageString("D_0__Bacteria; D_1__Bacillota"), std::invalid_argument);
	CHECK_THROWS_AS(ParseLineageString("d_Bacteria"), std::invalid_argument);
	CHECK_THROWS_AS(ParseLineageString("D__Bacteria"), std::invalid_argument);
}

// A name that starts with whitespace, or with a second rank prefix (usually one the caller
// forgot to strip), would differ from the same name spelled correctly elsewhere.
TEST_CASE("ParseLineageString rejects a name starting with whitespace or a second prefix", "[lineage][parse_lineage]") {
	// Field trimming cannot reach whitespace between the prefix and the name.
	CHECK_THROWS_WITH(ParseLineageString("d__ Bacteria"), ContainsSubstring("whitespace"));
	CHECK_THROWS_WITH(ParseLineageString("d__Bacteria; p__p__Bacillota"), ContainsSubstring("second rank prefix"));
	// Any lowercase letter counts, not only the field's own.
	CHECK_THROWS_WITH(ParseLineageString("d__Bacteria; p__d__Bacillota"), ContainsSubstring("second rank prefix"));
}
