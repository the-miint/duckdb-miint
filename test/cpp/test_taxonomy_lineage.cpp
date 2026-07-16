#include <catch2/catch_test_macros.hpp>
#include "taxonomy_lineage.hpp"

using namespace miint;

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
