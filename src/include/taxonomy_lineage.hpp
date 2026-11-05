#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace miint {

// The eight rank columns of the shared lineage schema, shallowest first, with the
// letter each carries in a prefixed lineage string ("d__", "p__", ...).
struct LineageRank {
	const char *column;
	char letter;
};
inline constexpr std::array<LineageRank, 8> kLineageRanks = {{{"domain", 'd'},
                                                              {"phylum", 'p'},
                                                              {"class", 'c'},
                                                              {"order", 'o'},
                                                              {"family", 'f'},
                                                              {"genus", 'g'},
                                                              {"species", 's'},
                                                              {"strain", 't'}}};

// One value per rank column, in kLineageRanks order; nullopt is a missing rank.
using LineageRanks = std::array<std::optional<std::string>, kLineageRanks.size()>;

// Parse a prefixed lineage string ("d__A; p__B" or miint's "d__A;p__B") into the rank
// columns. Each ';'-separated field, trimmed, is <letter>__<name> with the letter
// matching its position; an empty name is a missing rank and absent trailing fields are
// missing. `domain_letter` ('d' or 'k') is the only letter accepted at position 0.
// Throws std::invalid_argument on anything else, including a name with leading
// whitespace or one that itself looks prefixed. Gaps are not checked here.
LineageRanks ParseLineageString(std::string_view lineage, char domain_letter = 'd');

// One node in a taxon's lineage (an ancestor, or the query taxon itself).
struct TaxonNode {
	int64_t taxid;
	std::string name;
	std::string rank;
};

// A rank-collapsed lineage for one query taxon. This is the shared schema emitted
// by both read_ncbi_lineage (online) and the taxonomy_lineage macro (offline over
// a taxdump tree) so the two are drop-in interchangeable. `tax_class` / `tax_order`
// back the SQL columns "class" / "order" (reserved words).
struct Lineage {
	int64_t taxid;    // the query taxon
	std::string name; // its scientific name
	std::string rank; // its rank
	std::string domain;
	std::string phylum;
	std::string tax_class;
	std::string tax_order;
	std::string family;
	std::string genus;
	std::string species;
	std::string strain;
	std::string lineage; // formatted: d__..;p__..;c__..;o__..;f__..;g__..;s__..[;t__..]

	// Retired taxids that NCBI has merged into this (current) taxon, from the
	// <AkaTaxIds> block of a taxonomy efetch. NOT part of the output schema — used
	// only so read_ncbi_lineage can recognise a requested-but-merged taxid as
	// resolved rather than reporting it as omitted. Empty for the offline path.
	std::vector<int64_t> aka_taxids;
};

// Map an NCBI rank string to its collapse column. Both the legacy "superkingdom"
// and the newer "domain" map to "domain"; the other seven standard ranks map to
// their like-named column; everything else (clade, subspecies, no rank, ...)
// returns "" and is not collapsed.
std::string RankColumn(const std::string &rank);

// Collapse a set of ancestors plus the query taxon into the standard ranks and the
// formatted lineage string.
Lineage CollapseLineage(int64_t taxid, const std::string &name, const std::string &rank,
                        const std::vector<TaxonNode> &ancestors);

// Parse NCBI taxonomy efetch XML (a <TaxaSet> of top-level <Taxon> elements, each
// with a nested <LineageEx>) into one collapsed Lineage per query taxon, in order.
std::vector<Lineage> ParseTaxonomyXML(const std::string &xml);

} // namespace miint
