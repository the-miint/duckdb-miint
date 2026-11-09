#include "taxonomy_lineage.hpp"

#include "ncbi_parser.hpp"

#include <cstdlib>
#include <stdexcept>

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
	if (const auto *column = FindRankColumn(rank)) {
		lin.*column->field = name;
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

// ASCII whitespace, trimmed from both ends of a lineage field. Names may contain it
// ("Escherichia coli") but may not start with it. Spelled out because std::isspace
// depends on the locale, and this file is DuckDB-free, so StringUtil::CharacterIsSpace
// is not available.
bool IsLineageSpace(char c) {
	return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
}

// A rank prefix is a letter and two underscores: "d__", "p__", ...
constexpr size_t kRankPrefixWidth = 3;

bool HasRankPrefix(std::string_view field, char letter) {
	return field.size() >= kRankPrefixWidth && field[0] == letter && field[1] == '_' && field[2] == '_';
}

// The prefix letter of `rank` in a lineage string: its own, except the domain's, which the
// caller declares ('d' or 'k').
char RankLetter(size_t rank, char domain_letter) {
	return rank == 0 ? domain_letter : kLineageRanks[rank].letter;
}

// The rank whose prefix `field` starts with, or kLineageRanks.size() if none.
size_t PrefixRank(std::string_view field, char domain_letter) {
	size_t rank = 0;
	while (rank < kLineageRanks.size() && !HasRankPrefix(field, RankLetter(rank, domain_letter))) {
		++rank;
	}
	return rank;
}

// Why a lineage's first field has the wrong prefix, with a hint when a parameter would
// read it.
std::string FirstFieldProblem(std::string_view field, char domain_letter, bool allow_subset) {
	std::string why;
	if (allow_subset) {
		why = "must start with a rank prefix, one of";
		for (size_t rank = 0; rank < kLineageRanks.size(); ++rank) {
			why += std::string(" ") + RankLetter(rank, domain_letter) + "__";
		}
	} else {
		why = "must start with '" + std::string(1, domain_letter) + "__'";
	}
	const char other = domain_letter == 'd' ? 'k' : 'd';
	why += " (domain_letter is '" + std::string(1, domain_letter) + "'";
	if (HasRankPrefix(field, other)) {
		const std::string o(1, other);
		why += "; pass domain_letter := '" + o + "' to read " + o + "__ lineages";
	} else if (!allow_subset && PrefixRank(field, domain_letter) < kLineageRanks.size()) {
		why += "; pass allow_subset := true to read lineages that start below the domain";
	}
	return why + ")";
}

// Every rank is written, empty or not, except the last (strain), which is written only
// when present.
std::string FormatLineage(const Lineage &l) {
	std::string s;
	for (size_t i = 0; i < kLineageRanks.size(); ++i) {
		const auto &rank = kLineageRanks[i];
		const std::string &name = l.*rank.field;
		if (i + 1 == kLineageRanks.size() && name.empty()) {
			break;
		}
		if (i > 0) {
			s += ';';
		}
		s += rank.letter;
		s += "__";
		s += name;
	}
	return s;
}

} // namespace

const LineageRank *FindRankColumn(std::string_view ncbi_rank) {
	if (ncbi_rank == "superkingdom") {
		return &kLineageRanks[0]; // domain
	}
	for (const auto &rank : kLineageRanks) {
		if (ncbi_rank == rank.column) {
			return &rank;
		}
	}
	return nullptr;
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

LineageRanks ParseLineageString(std::string_view lineage, char domain_letter, bool allow_subset) {
	if (domain_letter != 'd' && domain_letter != 'k') {
		throw std::invalid_argument("domain_letter must be 'd' or 'k', got '" + std::string(1, domain_letter) + "'");
	}
	auto fail = [&](size_t index, std::string_view field, const std::string &why) {
		return std::invalid_argument("field " + std::to_string(index + 1) + " ('" + std::string(field) + "') " + why +
		                             " in lineage '" + std::string(lineage) + "'");
	};

	LineageRanks ranks;
	size_t index = 0; // the field's place in the string
	size_t rank = 0;  // the rank it must carry
	size_t start = 0;
	while (true) {
		size_t stop = lineage.find(';', start);
		size_t end = stop == std::string_view::npos ? lineage.size() : stop;
		while (start < end && IsLineageSpace(lineage[start])) {
			++start;
		}
		while (end > start && IsLineageSpace(lineage[end - 1])) {
			--end;
		}
		const std::string_view field = lineage.substr(start, end - start);

		if (rank == ranks.size()) {
			// A subset reaches strain in fewer than eight fields, so its count would mislead.
			throw fail(index, field,
			           allow_subset
			               ? "is past strain, the deepest rank"
			               : "is past strain: a lineage has at most " + std::to_string(ranks.size()) + " fields");
		}
		if (field.empty()) {
			throw fail(index, field, "is empty");
		}
		if (index == 0) {
			// With allow_subset the first field's prefix picks the rank the lineage starts at.
			if (allow_subset) {
				rank = PrefixRank(field, domain_letter);
			}
			if (rank == ranks.size() || !HasRankPrefix(field, RankLetter(rank, domain_letter))) {
				throw fail(index, field, FirstFieldProblem(field, domain_letter, allow_subset));
			}
		} else if (!HasRankPrefix(field, RankLetter(rank, domain_letter))) {
			throw fail(index, field, "must start with '" + std::string(1, RankLetter(rank, domain_letter)) + "__'");
		}
		const std::string_view name = field.substr(kRankPrefixWidth);
		if (!name.empty()) {
			if (IsLineageSpace(name[0])) {
				throw fail(index, field, "has whitespace after its rank prefix");
			}
			if (name[0] >= 'a' && name[0] <= 'z' && HasRankPrefix(name, name[0])) {
				throw fail(index, field, "has a second rank prefix");
			}
			ranks[rank] = std::string(name);
		}
		++index;
		++rank;
		if (stop == std::string_view::npos) {
			break;
		}
		start = stop + 1;
	}
	return ranks;
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
