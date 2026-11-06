#include "TaxonomyStore.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <utility>

namespace miint {

namespace {

constexpr size_t kRanks = kLineageRanks.size();
constexpr size_t kMaxListedProblems = 10;
constexpr size_t kMaxListedParents = 3;

// Sets `out` to `s` with ASCII letters lowercased, reusing out's buffer. Spelled out
// because std::tolower depends on the locale.
void AsciiLower(std::string_view s, std::string &out) {
	out.assign(s);
	for (auto &c : out) {
		if (c >= 'A' && c <= 'Z') {
			c = static_cast<char>(c - 'A' + 'a');
		}
	}
}

// Why `value` cannot be a rank name, or nullptr if it can.
const char *ValueProblem(const std::string &value) {
	if (value.empty()) {
		return "is empty";
	}
	if (IsLineageSpace(value.front()) || IsLineageSpace(value.back())) {
		return "starts or ends with whitespace";
	}
	if (StartsWithRankPrefix(value)) {
		return "starts with a rank prefix such as 'p__'";
	}
	if (value.find(';') != std::string::npos) {
		return "contains ';'";
	}
	return nullptr;
}

// Counts every problem and keeps the first few, one per line in input order, so the error
// for a broken taxonomy stays a bounded number of lines. FormatIdList (absquant.cpp) and
// join_names (NewickTree.cpp) list sorted, deduplicated ids; each entry here is a whole
// problem with its reason, kept in input order so the caller can find the row.
class Problems {
public:
	void Add(std::string problem) {
		if (listed_.size() < kMaxListedProblems) {
			listed_.push_back(std::move(problem));
		}
		count_++;
	}

	void ThrowIfAny() const {
		if (count_ == 0) {
			return;
		}
		std::string message = "taxonomy: " + std::to_string(count_) + (count_ == 1 ? " problem" : " problems");
		message += count_ > listed_.size() ? "; the first " + std::to_string(listed_.size()) + ":" : ":";
		for (const auto &problem : listed_) {
			message += "\n  " + problem;
		}
		throw std::invalid_argument(message);
	}

private:
	std::vector<std::string> listed_;
	size_t count_ = 0;
};

// How a problem names its row: by feature_id, or by position when the id is empty.
std::string RowText(const std::string &feature_id, size_t row) {
	if (feature_id.empty()) {
		return "row " + std::to_string(row + 1) + " (empty feature_id)";
	}
	return "feature_id '" + feature_id + "'";
}

} // namespace

const std::vector<std::string> &DefaultBadNameTokens() {
	// BAD_NAMES in biocore/tax2tree's t2t/nlevel.py, in its order. The first token that
	// matches is the one recorded, so "environmental samples" is recorded as
	// "environmental sample". Then incertae sedis, spelled with a space and with an
	// underscore, which tax2tree reads as a name: it marks members of uncertain placement,
	// not a clade.
	static const std::vector<std::string> tokens = {
	    "environmental sample", "uncultured",    "UNNAMEABLE", "unclassified",
	    "unidentified",         "cluster",       "isolate",    "environmental samples",
	    "incertae sedis",       "incertae_sedis"};
	return tokens;
}

TaxonomyStore::TaxonomyStore(const std::vector<TaxonomyRow> &rows, const std::vector<std::string> &bad_name_tokens) {
	std::vector<std::string> lower_tokens;
	for (const auto &token : bad_name_tokens) {
		if (token.empty()) {
			throw std::invalid_argument("taxonomy: a bad-name token is empty, which would match every name");
		}
		lower_tokens.emplace_back();
		AsciiLower(token, lower_tokens.back());
	}

	// Row checks: ids, values, gaps. Each row is checked in full before any is cleaned.
	feature_ids_.reserve(rows.size());
	for (const auto &row : rows) {
		feature_ids_.push_back(row.feature_id);
	}
	// How many times each repeated id occurs; ids that occur once are not listed.
	std::unordered_map<std::string_view, size_t> repeats;
	for (size_t i = 0; i < feature_ids_.size(); i++) {
		if (!feature_ids_[i].empty() && !row_of_id_.emplace(feature_ids_[i], i).second) {
			repeats.emplace(feature_ids_[i], 1).first->second++;
		}
	}

	Problems problems;
	for (size_t i = 0; i < rows.size(); i++) {
		const auto &id = feature_ids_[i];
		const auto &ranks = rows[i].ranks;
		if (id.empty()) {
			problems.Add("row " + std::to_string(i + 1) + " has an empty feature_id");
		} else if (auto repeat = repeats.find(id); repeat != repeats.end() && row_of_id_.at(id) == i) {
			problems.Add(RowText(id, i) + " appears " + std::to_string(repeat->second) + " times");
		}
		for (size_t r = 0; r < kRanks; r++) {
			if (!ranks[r]) {
				continue;
			}
			if (const char *why = ValueProblem(*ranks[r])) {
				problems.Add(RowText(id, i) + ", column '" + kLineageRanks[r].column + "': value '" + *ranks[r] + "' " +
				             why);
			}
		}
		for (size_t r = 0; r + 1 < kRanks; r++) {
			if (ranks[r]) {
				continue;
			}
			for (size_t deeper = r + 1; deeper < kRanks; deeper++) {
				if (ranks[deeper]) {
					problems.Add(RowText(id, i) + ": '" + kLineageRanks[r].column + "' is NULL but deeper rank '" +
					             kLineageRanks[deeper].column + "' is not");
					break;
				}
			}
			break;
		}
	}
	problems.ThrowIfAny();

	// Clean, intern and link each lineage to its parent names.
	names_.assign(rows.size() * kRanks, kMissingName);
	std::array<std::unordered_map<std::string, NameId>, kRanks> name_ids;
	// The row that first interned each name, to show the lineage a parent came from.
	std::array<std::vector<size_t>, kRanks> first_row;
	// Names seen under more than one parent, in the order the second parent appeared. Each
	// distinct parent is kept with the row it first came from.
	struct Conflict {
		size_t rank;
		NameId name;
		std::vector<std::pair<NameId, size_t>> parents;
	};
	std::vector<Conflict> conflicts;
	std::map<std::pair<size_t, NameId>, size_t> conflict_index;

	std::string lower;
	for (size_t i = 0; i < rows.size(); i++) {
		const auto &ranks = rows[i].ranks;
		// With no gaps, the present ranks are a prefix of the lineage.
		size_t depth = 0;
		while (depth < kRanks && ranks[depth]) {
			depth++;
		}
		const size_t given_depth = depth;

		for (size_t r = 0; r < depth && !lower_tokens.empty(); r++) {
			AsciiLower(*ranks[r], lower);
			for (size_t t = 0; t < lower_tokens.size(); t++) {
				if (lower.find(lower_tokens[t]) != std::string::npos) {
					truncations_.push_back({i, r, bad_name_tokens[t]});
					depth = r;
					break;
				}
			}
		}
		if (depth < 2) {
			if (given_depth > 0) {
				wiped_rows_.push_back(i);
			}
			depth = 0;
		}

		for (size_t r = 0; r < depth; r++) {
			const NameId parent = r == 0 ? kMissingName : names_[i * kRanks + r - 1];
			auto [it, inserted] = name_ids[r].emplace(*ranks[r], static_cast<NameId>(name_text_[r].size()));
			const NameId name = it->second;
			if (inserted) {
				if (name == kMissingName) {
					throw std::invalid_argument("taxonomy: more than " + std::to_string(kMissingName) +
					                            " distinct names at rank '" + kLineageRanks[r].column + "'");
				}
				name_text_[r].push_back(*ranks[r]);
				parent_[r].push_back(parent);
				first_row[r].push_back(i);
			} else if (parent_[r][name] != parent) {
				auto [entry, added] = conflict_index.emplace(std::make_pair(r, name), conflicts.size());
				if (added) {
					conflicts.push_back({r, name, {{parent_[r][name], first_row[r][name]}}});
				}
				auto &seen = conflicts[entry->second].parents;
				auto same_parent = [&](const std::pair<NameId, size_t> &p) {
					return p.first == parent;
				};
				if (std::none_of(seen.begin(), seen.end(), same_parent)) {
					seen.emplace_back(parent, i);
				}
			}
			names_[i * kRanks + r] = name;
		}
	}

	for (const auto &conflict : conflicts) {
		const size_t listed = std::min(conflict.parents.size(), kMaxListedParents);
		std::string problem = "'" + name_text_[conflict.rank][conflict.name] + "' (" +
		                      kLineageRanks[conflict.rank].column + ") is under " +
		                      std::to_string(conflict.parents.size()) + " different " +
		                      kLineageRanks[conflict.rank - 1].column + " names:";
		for (size_t p = 0; p < listed; p++) {
			// The lineage above the name in the row the parent came from, joined with "; " as
			// tax2tree writes labels, for reading. (FormatLineage writes lineages for parsing.)
			const size_t row = conflict.parents[p].second;
			std::string lineage;
			for (size_t r = 0; r < conflict.rank; r++) {
				lineage += (r == 0 ? "" : "; ") + std::string(1, kLineageRanks[r].letter) + "__" +
				           name_text_[r][names_[row * kRanks + r]];
			}
			problem += (p == 0 ? " [" : ", [") + lineage + "]";
		}
		if (conflict.parents.size() > listed) {
			problem += ", and " + std::to_string(conflict.parents.size() - listed) + " more";
		}
		problems.Add(std::move(problem));
	}
	problems.ThrowIfAny();
}

size_t TaxonomyStore::RowCount() const {
	return feature_ids_.size();
}

const std::string &TaxonomyStore::FeatureId(size_t row) const {
	return feature_ids_[row];
}

std::optional<size_t> TaxonomyStore::FindRow(std::string_view feature_id) const {
	auto it = row_of_id_.find(feature_id);
	if (it == row_of_id_.end()) {
		return std::nullopt;
	}
	return it->second;
}

NameId TaxonomyStore::Name(size_t row, size_t rank) const {
	return names_[row * kRanks + rank];
}

bool TaxonomyStore::IsInformative(size_t row) const {
	return Name(row, 0) != kMissingName;
}

size_t TaxonomyStore::NameCount(size_t rank) const {
	return name_text_[rank].size();
}

const std::string &TaxonomyStore::NameText(size_t rank, NameId name) const {
	return name_text_[rank][name];
}

NameId TaxonomyStore::Parent(size_t rank, NameId name) const {
	return parent_[rank][name];
}

const std::vector<LineageTruncation> &TaxonomyStore::Truncations() const {
	return truncations_;
}

const std::vector<size_t> &TaxonomyStore::WipedRows() const {
	return wiped_rows_;
}

} // namespace miint
