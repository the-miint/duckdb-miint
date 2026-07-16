#include "taxdump_parser.hpp"

#include <cstdlib>
#include <unordered_map>
#include <unordered_set>

namespace miint {

namespace {

// The name class in names.dmp that carries the canonical scientific name.
constexpr const char *kScientificName = "scientific name";

// Invoke `fn(line_view)` for each non-empty line in `buf` (split on '\n', trailing
// '\r' stripped). Manual iteration avoids stream overhead at taxdump scale (~3.5M
// lines in names.dmp).
template <typename Fn>
void ForEachLine(const std::string &buf, Fn &&fn) {
	const size_t n = buf.size();
	size_t start = 0;
	while (start < n) {
		size_t nl = buf.find('\n', start);
		size_t end = (nl == std::string::npos) ? n : nl;
		size_t line_end = end;
		if (line_end > start && buf[line_end - 1] == '\r') {
			--line_end;
		}
		if (line_end > start) {
			fn(std::string_view(buf.data() + start, line_end - start));
		}
		start = end + 1;
	}
}

int64_t ParseInt(const std::string &s) {
	return static_cast<int64_t>(std::strtoll(s.c_str(), nullptr, 10));
}

} // namespace

std::vector<std::string> TaxdumpParser::SplitFields(std::string_view sv) {
	// A dmp line is  FIELD ("\t|\t" FIELD)* "\t|"  optionally followed by "\n".
	if (!sv.empty() && sv.back() == '\n') {
		sv.remove_suffix(1);
	}
	if (!sv.empty() && sv.back() == '\r') {
		sv.remove_suffix(1);
	}
	// Strip the trailing "\t|" line terminator (leaving the final "\t|\t" separator
	// intact so a trailing empty field is still recovered by the split below).
	if (sv.size() >= 2 && sv[sv.size() - 2] == '\t' && sv[sv.size() - 1] == '|') {
		sv.remove_suffix(2);
	}

	std::vector<std::string> fields;
	constexpr std::string_view sep = "\t|\t";
	size_t pos = 0;
	while (true) {
		size_t next = sv.find(sep, pos);
		if (next == std::string_view::npos) {
			fields.emplace_back(sv.substr(pos));
			break;
		}
		fields.emplace_back(sv.substr(pos, next - pos));
		pos = next + sep.size();
	}
	return fields;
}

std::vector<TaxdumpNode> TaxdumpParser::ParseNodes(const std::string &nodes_dmp, const std::string &names_dmp) {
	// Pass 1: scientific name per taxid.
	std::unordered_map<int64_t, std::string> scientific_names;
	ForEachLine(names_dmp, [&](std::string_view line) {
		auto f = SplitFields(line);
		if (f.size() < 4 || f[3] != kScientificName) {
			return;
		}
		scientific_names.emplace(ParseInt(f[0]), f[1]);
	});

	// Pass 2: raw (taxid, parent, rank) plus the set of referenced parents (for is_tip).
	struct RawNode {
		int64_t taxid;
		int64_t parent;
		std::string rank;
	};
	std::vector<RawNode> raw;
	std::unordered_set<int64_t> parents;
	ForEachLine(nodes_dmp, [&](std::string_view line) {
		auto f = SplitFields(line);
		if (f.size() < 3) {
			return;
		}
		RawNode rn {ParseInt(f[0]), ParseInt(f[1]), f[2]};
		parents.insert(rn.parent);
		raw.push_back(std::move(rn));
	});

	// Assemble live nodes.
	std::vector<TaxdumpNode> nodes;
	nodes.reserve(raw.size());
	for (auto &rn : raw) {
		TaxdumpNode node;
		node.taxid = rn.taxid;
		node.parent_taxid = (rn.parent == rn.taxid) ? std::nullopt : std::optional<int64_t>(rn.parent);
		node.rank = std::move(rn.rank);
		auto it = scientific_names.find(rn.taxid);
		if (it != scientific_names.end()) {
			node.name = it->second;
		}
		node.is_tip = parents.find(rn.taxid) == parents.end();
		nodes.push_back(std::move(node));
	}
	return nodes;
}

std::vector<TaxdumpMerge> TaxdumpParser::ParseMerged(const std::string &merged_dmp) {
	std::vector<TaxdumpMerge> merged;
	ForEachLine(merged_dmp, [&](std::string_view line) {
		auto f = SplitFields(line);
		if (f.size() < 2) {
			return;
		}
		merged.push_back(TaxdumpMerge {ParseInt(f[0]), ParseInt(f[1])});
	});
	return merged;
}

std::vector<int64_t> TaxdumpParser::ParseDeleted(const std::string &delnodes_dmp) {
	std::vector<int64_t> deleted;
	ForEachLine(delnodes_dmp, [&](std::string_view line) {
		auto f = SplitFields(line);
		if (f.empty() || f[0].empty()) {
			return;
		}
		deleted.push_back(ParseInt(f[0]));
	});
	return deleted;
}

} // namespace miint
