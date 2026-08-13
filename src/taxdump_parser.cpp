#include "taxdump_parser.hpp"

#include <charconv>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <unordered_set>

namespace miint {

namespace {

// The name class in names.dmp that carries the canonical scientific name.
constexpr const char *kScientificName = "scientific name";

// Invoke `fn(line_view, line_number)` for each non-empty line in `buf` (split on '\n',
// trailing '\r' stripped). Manual iteration avoids stream overhead at taxdump scale
// (names.dmp is ~5M lines). `line_number` is 1-based and counts blank lines too, so it
// matches what a text editor shows for the offending row.
template <typename Fn>
void ForEachLine(const std::string &buf, Fn &&fn) {
	const size_t n = buf.size();
	size_t start = 0;
	size_t line_number = 0;
	while (start < n) {
		size_t nl = buf.find('\n', start);
		size_t end = (nl == std::string::npos) ? n : nl;
		size_t line_end = end;
		if (line_end > start && buf[line_end - 1] == '\r') {
			--line_end;
		}
		++line_number;
		if (line_end > start) {
			fn(std::string_view(buf.data() + start, line_end - start), line_number);
		}
		start = end + 1;
	}
}

// Parse a taxid field, rejecting anything that is not the whole field as an integer.
//
// from_chars (rather than strtoll) because it reports BOTH failure channels that
// strtoll discarded: `ec` distinguishes a non-numeric field and an out-of-range one,
// and `ptr` catches a partial parse. Those last two mattered most -- "123abc" used to
// yield 123 and an overflowing value used to clamp to INT64_MAX, both landing in the
// valid taxid range where no consumer could tell them from real data.
int64_t ParseTaxid(const std::string &field, const char *member, size_t line) {
	int64_t value = 0;
	const char *begin = field.data();
	const char *end = begin + field.size();
	auto [ptr, ec] = std::from_chars(begin, end, value);
	if (ec != std::errc {} || ptr != end) {
		throw std::runtime_error("taxdump: " + std::string(member) + " line " + std::to_string(line) + ": '" + field +
		                         "' is not a valid taxid");
	}
	return value;
}

// Enforce the field count NCBI's taxdump_readme.txt documents for a member whose every
// field we consume (names.dmp 4, merged.dmp 2, delnodes.dmp 1).
//
// This is an exact check, not a lower bound. A lower bound only catches truncation; it
// lets a row with an inserted column through to be read positionally, so a `name_class`
// that is really a `unique name` reaches the consumer looking perfectly well-formed. An
// exact check also rejects a merely-appended column, which is the conservative trade: a
// loud error naming the member and line is recoverable, a silently mislabelled name is
// not.
//
// nodes.dmp is excluded on purpose -- the readme documents 13 fields for taxdump and 18
// for new_taxdump, so its width legitimately varies. It gets a lower bound instead.
void RequireFieldCount(size_t actual, size_t expected, const char *member, size_t line) {
	if (actual != expected) {
		throw std::runtime_error("taxdump: " + std::string(member) + " line " + std::to_string(line) + " has " +
		                         std::to_string(actual) + " field(s); expected exactly " + std::to_string(expected) +
		                         ". The dump layout may have changed; miint's parser needs updating.");
	}
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
	ForEachLine(names_dmp, [&](std::string_view line, size_t line_number) {
		auto f = SplitFields(line);
		RequireFieldCount(f.size(), 4, "names.dmp", line_number);
		if (f[3] != kScientificName) {
			return;
		}
		scientific_names.emplace(ParseTaxid(f[0], "names.dmp", line_number), f[1]);
	});

	// Pass 2: raw (taxid, parent, rank) plus the set of referenced parents (for is_tip).
	struct RawNode {
		int64_t taxid;
		int64_t parent;
		std::string rank;
	};
	std::vector<RawNode> raw;
	std::unordered_set<int64_t> parents;
	ForEachLine(nodes_dmp, [&](std::string_view line, size_t line_number) {
		auto f = SplitFields(line);
		// Lower bound, not an exact count: nodes.dmp is 13 fields in taxdump and 18 in
		// new_taxdump, and only the first three are consumed, so appended columns must
		// keep parsing. Fewer than three is a truncated row and now fails loudly rather
		// than silently shortening the tree.
		if (f.size() < 3) {
			throw std::runtime_error("taxdump: nodes.dmp line " + std::to_string(line_number) + " has " +
			                         std::to_string(f.size()) +
			                         " field(s); expected at least 3 (tax_id, parent tax_id, rank)");
		}
		RawNode rn {ParseTaxid(f[0], "nodes.dmp", line_number), ParseTaxid(f[1], "nodes.dmp", line_number), f[2]};
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
	ForEachLine(merged_dmp, [&](std::string_view line, size_t line_number) {
		auto f = SplitFields(line);
		RequireFieldCount(f.size(), 2, "merged.dmp", line_number);
		// Both ids are checked: a bad new_tax_id would silently remap a retired taxid
		// onto the wrong live node, which is worse than failing the read.
		merged.push_back(
		    TaxdumpMerge {ParseTaxid(f[0], "merged.dmp", line_number), ParseTaxid(f[1], "merged.dmp", line_number)});
	});
	return merged;
}

std::vector<TaxdumpName> TaxdumpParser::ParseNames(const std::string &names_dmp) {
	std::vector<TaxdumpName> names;
	ForEachLine(names_dmp, [&](std::string_view line, size_t line_number) {
		auto f = SplitFields(line);
		RequireFieldCount(f.size(), 4, "names.dmp", line_number);
		int64_t taxid = ParseTaxid(f[0], "names.dmp", line_number);
		names.push_back(TaxdumpName {taxid, std::move(f[1]), std::move(f[2]), std::move(f[3])});
	});
	return names;
}

std::vector<int64_t> TaxdumpParser::ParseDeleted(const std::string &delnodes_dmp) {
	std::vector<int64_t> deleted;
	ForEachLine(delnodes_dmp, [&](std::string_view line, size_t line_number) {
		auto f = SplitFields(line);
		// The previous guard here (`f.empty() || f[0].empty()`) could not reject a wrong
		// width, because SplitFields always emplaces at least one element -- so any
		// non-blank line became a taxid, and a file that was not delnodes.dmp yielded
		// one bogus row per line instead of an error.
		RequireFieldCount(f.size(), 1, "delnodes.dmp", line_number);
		deleted.push_back(ParseTaxid(f[0], "delnodes.dmp", line_number));
	});
	return deleted;
}

} // namespace miint
