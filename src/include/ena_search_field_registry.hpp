#pragma once

#include <set>
#include <string>

namespace miint {

// Curated allowlist of ENA Portal API sample-search field names that we know
// behave sensibly as a structured `/search?fields=...&query=...` equality
// predicate. Used by the filter-pushdown path in `read_ena_attributes` to
// decide whether a `WHERE tag='X'` filter can be translated into a structured
// ENA search or must fall back to the XML path.
//
// Policy: if ANY tag referenced in a user's filter is not in this set, we do
// NOT push anything down — we fall back entirely to the XML path. This keeps
// the structured path equivalent to the XML output (per the plan's decision
// #4 in localdocs/PLAN-ena-predicate-maxseqs.md).
//
// Lookup is case-insensitive and whitespace-tolerant. ENA's own field names
// are lowercase, but user `WHERE` clauses often use mixed case.
class ENASearchFieldRegistry {
public:
	// Returns true iff `tag`, after lowercasing + trimming, matches a known
	// searchable sample field.
	static bool IsSearchableSampleField(const std::string &tag);

	// Canonical lowercase set. Exposed for tests and for diagnostics that want
	// to enumerate the allowlist (e.g., to emit a helpful error message).
	static const std::set<std::string> &SearchableSampleFields();
};

} // namespace miint
