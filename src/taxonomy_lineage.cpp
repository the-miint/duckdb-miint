#include "taxonomy_lineage.hpp"

#include "ncbi_parser.hpp"

#include <cstdlib>

namespace miint {

namespace {

int64_t ParseInt(const std::string &s) {
	return static_cast<int64_t>(std::strtoll(s.c_str(), nullptr, 10));
}

// Split `xml` into the inner content of each TOP-LEVEL <Taxon>...</Taxon> element,
// using depth counting so the nested <Taxon> elements inside <LineageEx> do not
// cause a split. NCBI emits both levels as bare "<Taxon>" (no attributes).
std::vector<std::string> TopLevelTaxonBlocks(const std::string &xml) {
	static const std::string kOpen = "<Taxon>";
	static const std::string kClose = "</Taxon>";
	std::vector<std::string> blocks;
	size_t pos = 0;
	int depth = 0;
	size_t block_start = std::string::npos;
	while (pos < xml.size()) {
		size_t next_open = xml.find(kOpen, pos);
		size_t next_close = xml.find(kClose, pos);
		if (next_close == std::string::npos) {
			break; // no more complete elements
		}
		if (next_open != std::string::npos && next_open < next_close) {
			if (depth == 0) {
				block_start = next_open + kOpen.size();
			}
			++depth;
			pos = next_open + kOpen.size();
		} else {
			--depth;
			if (depth == 0 && block_start != std::string::npos) {
				blocks.push_back(xml.substr(block_start, next_close - block_start));
				block_start = std::string::npos;
			}
			pos = next_close + kClose.size();
		}
	}
	return blocks;
}

void AssignRank(Lineage &lin, const std::string &rank, const std::string &name) {
	std::string col = RankColumn(rank);
	if (col == "domain") {
		lin.domain = name;
	} else if (col == "phylum") {
		lin.phylum = name;
	} else if (col == "class") {
		lin.tax_class = name;
	} else if (col == "order") {
		lin.tax_order = name;
	} else if (col == "family") {
		lin.family = name;
	} else if (col == "genus") {
		lin.genus = name;
	} else if (col == "species") {
		lin.species = name;
	} else if (col == "strain") {
		lin.strain = name;
	}
}

// Extract every <TaxId>N</TaxId> value in `xml` (used for the <AkaTaxIds> block,
// which may list more than one retired taxid).
std::vector<int64_t> ExtractAllTaxIds(const std::string &xml) {
	static const std::string kOpen = "<TaxId>";
	static const std::string kClose = "</TaxId>";
	std::vector<int64_t> ids;
	size_t pos = 0;
	while (true) {
		size_t s = xml.find(kOpen, pos);
		if (s == std::string::npos) {
			break;
		}
		size_t inner = s + kOpen.size();
		size_t e = xml.find(kClose, inner);
		if (e == std::string::npos) {
			break;
		}
		ids.push_back(ParseInt(xml.substr(inner, e - inner)));
		pos = e + kClose.size();
	}
	return ids;
}

std::string FormatLineage(const Lineage &l) {
	std::string s = "d__" + l.domain + ";p__" + l.phylum + ";c__" + l.tax_class + ";o__" + l.tax_order + ";f__" +
	                l.family + ";g__" + l.genus + ";s__" + l.species;
	if (!l.strain.empty()) {
		s += ";t__" + l.strain;
	}
	return s;
}

} // namespace

std::string RankColumn(const std::string &rank) {
	if (rank == "superkingdom" || rank == "domain") {
		return "domain";
	}
	if (rank == "phylum" || rank == "class" || rank == "order" || rank == "family" || rank == "genus" ||
	    rank == "species" || rank == "strain") {
		return rank;
	}
	return {};
}

Lineage CollapseLineage(int64_t taxid, const std::string &name, const std::string &rank,
                        const std::vector<TaxonNode> &ancestors) {
	Lineage lin;
	lin.taxid = taxid;
	lin.name = name;
	lin.rank = rank;
	for (const auto &a : ancestors) {
		AssignRank(lin, a.rank, a.name);
	}
	// The query taxon itself last, so its own rank (e.g. species/strain) wins over
	// any ancestor at the same level.
	AssignRank(lin, rank, name);
	lin.lineage = FormatLineage(lin);
	return lin;
}

std::vector<Lineage> ParseTaxonomyXML(const std::string &xml) {
	std::vector<Lineage> out;
	for (const auto &taxon : TopLevelTaxonBlocks(xml)) {
		// The query taxon's own TaxId/ScientificName/Rank are the first occurrences,
		// appearing before the nested <LineageEx>.
		int64_t taxid = ParseInt(NCBIParser::ExtractXMLTagValue(taxon, "TaxId"));
		std::string name = NCBIParser::ExtractXMLTagValue(taxon, "ScientificName");
		std::string rank = NCBIParser::ExtractXMLTagValue(taxon, "Rank");

		std::vector<TaxonNode> ancestors;
		std::string lineage_ex = NCBIParser::ExtractXMLTagValue(taxon, "LineageEx");
		for (const auto &anc : TopLevelTaxonBlocks(lineage_ex)) {
			ancestors.push_back(TaxonNode {ParseInt(NCBIParser::ExtractXMLTagValue(anc, "TaxId")),
			                               NCBIParser::ExtractXMLTagValue(anc, "ScientificName"),
			                               NCBIParser::ExtractXMLTagValue(anc, "Rank")});
		}
		Lineage lin = CollapseLineage(taxid, name, rank, ancestors);
		// Retired taxids merged into this taxon (present only when a queried taxid
		// was itself a merged id). The query taxid above is the FIRST <TaxId>, which
		// is the current one, so this only adds the aka ids.
		lin.aka_taxids = ExtractAllTaxIds(NCBIParser::ExtractXMLTagValue(taxon, "AkaTaxIds"));
		out.push_back(std::move(lin));
	}
	return out;
}

} // namespace miint
