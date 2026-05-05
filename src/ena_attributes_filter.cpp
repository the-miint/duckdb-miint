#include "ena_attributes_filter.hpp"

#include "ena_parser.hpp"
#include "ena_search_field_registry.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace miint {

// ---- ENAFilterNode factories ----

std::unique_ptr<ENAFilterNode> ENAFilterNode::MakeEqual(const std::string &column, const std::string &value) {
	auto node = std::make_unique<ENAFilterNode>();
	node->kind = Kind::EQUAL;
	node->column = column;
	node->values = {value};
	return node;
}

std::unique_ptr<ENAFilterNode> ENAFilterNode::MakeIn(const std::string &column,
                                                     const std::vector<std::string> &values) {
	auto node = std::make_unique<ENAFilterNode>();
	node->kind = Kind::IN;
	node->column = column;
	node->values = values;
	return node;
}

std::unique_ptr<ENAFilterNode> ENAFilterNode::MakeAnd(std::vector<std::unique_ptr<ENAFilterNode>> children) {
	auto node = std::make_unique<ENAFilterNode>();
	node->kind = Kind::AND_CONJUNCTION;
	node->children = std::move(children);
	return node;
}

std::unique_ptr<ENAFilterNode> ENAFilterNode::MakeOr(std::vector<std::unique_ptr<ENAFilterNode>> children) {
	auto node = std::make_unique<ENAFilterNode>();
	node->kind = Kind::OR_CONJUNCTION;
	node->children = std::move(children);
	return node;
}

std::unique_ptr<ENAFilterNode> ENAFilterNode::MakeUnsupported() {
	auto node = std::make_unique<ENAFilterNode>();
	node->kind = Kind::UNSUPPORTED;
	return node;
}

// ---- Extractor ----

namespace {

std::string ToLower(std::string s) {
	for (auto &c : s) {
		c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
	}
	return s;
}

// Per-conjunct analysis result. Short-circuit flag (`reject`) lets the caller
// abort the whole pushdown on the first unsupported atom.
struct AtomResult {
	bool reject = false;
	bool is_tag_eq = false;
	bool is_tag_in = false;
	bool is_value_eq = false;
	std::string tag_value;               // for is_tag_eq
	std::vector<std::string> tag_values; // for is_tag_in (canonicalized)
	std::string pinned_value;            // for is_value_eq
};

AtomResult AnalyzeAtom(const ENAFilterNode &node) {
	AtomResult r;
	switch (node.kind) {
	case ENAFilterNode::Kind::EQUAL: {
		if (node.values.size() != 1) {
			r.reject = true;
			return r;
		}
		auto col = ToLower(node.column);
		if (col == "tag") {
			auto canonical = ToLower(node.values[0]);
			if (!ENASearchFieldRegistry::IsSearchableSampleField(canonical)) {
				r.reject = true;
				return r;
			}
			r.is_tag_eq = true;
			r.tag_value = canonical;
			return r;
		}
		if (col == "value") {
			r.is_value_eq = true;
			r.pinned_value = node.values[0];
			return r;
		}
		// Predicate on any other column (e.g., sample_accession) — can't
		// cleanly translate; reject per plan's safety policy.
		r.reject = true;
		return r;
	}
	case ENAFilterNode::Kind::IN: {
		auto col = ToLower(node.column);
		if (col != "tag") {
			// Only `tag IN (...)` is supported. `value IN (...)` without a
			// pinned tag is ambiguous; any other column is out of scope.
			r.reject = true;
			return r;
		}
		if (node.values.empty()) {
			r.reject = true;
			return r;
		}
		for (const auto &v : node.values) {
			auto canonical = ToLower(v);
			if (!ENASearchFieldRegistry::IsSearchableSampleField(canonical)) {
				r.reject = true;
				return r;
			}
			r.tag_values.push_back(canonical);
		}
		r.is_tag_in = true;
		return r;
	}
	case ENAFilterNode::Kind::AND_CONJUNCTION:
	case ENAFilterNode::Kind::OR_CONJUNCTION:
	case ENAFilterNode::Kind::UNSUPPORTED:
		// OR / nested AND / unrecognized — out of scope for the strict
		// extractor. DuckDB pre-splits top-level AND into conjuncts, so a
		// remaining AND node here would be inside an OR or similar construct
		// we've already chosen not to decompose.
		r.reject = true;
		return r;
	}
	r.reject = true;
	return r;
}

} // namespace

ENAAttributePushdown ExtractPushdownPredicates(const std::vector<std::unique_ptr<ENAFilterNode>> &conjuncts) {
	ENAAttributePushdown empty; // signal for fallback

	if (conjuncts.empty()) {
		return empty;
	}

	std::vector<std::string> tag_atoms_eq;              // from EQUAL on tag
	std::vector<std::vector<std::string>> tag_atoms_in; // from IN on tag
	std::vector<std::string> value_atoms_eq;            // from EQUAL on value

	for (const auto &node : conjuncts) {
		if (!node) {
			return empty;
		}
		auto atom = AnalyzeAtom(*node);
		if (atom.reject) {
			return empty;
		}
		if (atom.is_tag_eq) {
			tag_atoms_eq.push_back(std::move(atom.tag_value));
		} else if (atom.is_tag_in) {
			tag_atoms_in.push_back(std::move(atom.tag_values));
		} else if (atom.is_value_eq) {
			value_atoms_eq.push_back(std::move(atom.pinned_value));
		}
	}

	// A pinned `value` predicate only makes sense alongside a single-tag pin.
	// If no tag pin exists, the value is ambiguous — fallback.
	if (!value_atoms_eq.empty() && tag_atoms_eq.empty()) {
		return empty;
	}
	// Conjoined `tag=X AND tag=Y` would require both — impossible since one
	// row has one tag. Treat as fallback rather than silently returning zero.
	if (tag_atoms_eq.size() > 1) {
		return empty;
	}
	// `tag=X AND tag IN (...)` is expressible but semantically odd; fall back
	// for safety rather than tangling the query shape.
	if (!tag_atoms_eq.empty() && !tag_atoms_in.empty()) {
		return empty;
	}
	// Multiple `tag IN (...)` conjuncts: the SQL semantics is an intersection
	// (a row must satisfy EVERY IN list simultaneously), but collecting the
	// union of tags here would broaden the wire query sent to ENA. Rather
	// than emit a semantically wider URL and rely on DuckDB's re-applied
	// filter to trim, fall back. Matches the extractor's "strict or nothing"
	// personality: the structured path runs only when we can prove
	// equivalence.
	if (tag_atoms_in.size() > 1) {
		return empty;
	}

	// Collect deduplicated canonical tags in first-seen order.
	ENAAttributePushdown out;
	std::unordered_set<std::string> seen;

	auto push_tag = [&](const std::string &t) {
		if (seen.insert(t).second) {
			out.tags.push_back(t);
		}
	};

	for (const auto &t : tag_atoms_eq) {
		push_tag(t);
	}
	for (const auto &tag_list : tag_atoms_in) {
		for (const auto &t : tag_list) {
			push_tag(t);
		}
	}

	if (out.tags.empty()) {
		return empty;
	}

	// Attach pinned values. Only tag_atoms_eq contributes a (tag, value) pair
	// (we already rejected value-without-tag-eq above).
	//
	// When there are multiple `value=` conjuncts on the same pinned tag
	// (e.g. `tag='X' AND value='Y1' AND value='Y2'`), ENA's search cannot
	// conjoin equality on the same field in one query. We INTENTIONALLY
	// emit only `value_atoms_eq[0]` to the wire and DROP value_atoms_eq[1..n].
	// This is safe only because plan decision #5 is honored: we never prune
	// entries from the original `filters` vector, so DuckDB re-applies ALL
	// value= predicates as a LogicalFilter above the scan. Any future refactor
	// that starts pruning filters would need to revisit this branch — at that
	// point the correct fix is to fall back to XML when value_atoms_eq.size()
	// > 1, not to silently drop constraints.
	if (!tag_atoms_eq.empty() && !value_atoms_eq.empty()) {
		out.tag_value_pairs.emplace_back(tag_atoms_eq[0], value_atoms_eq[0]);
	}

	return out;
}

// ---- URL builder ----

namespace {

// Shared URL-encoder for the sample /search endpoint. Emits
// `<in_column> IN ("acc1","acc2",...) [AND tag="v"]*` and returns columns
// `sample_accession,<tags>`. Factored so the sample-batch path and the
// study-direct path share a single implementation; the only moving part
// is `in_column` ("sample_accession" vs "study_accession").
std::string BuildAttributeSearchURL(const char *in_column, const std::vector<std::string> &accessions,
                                    const std::vector<std::string> &tags,
                                    const std::vector<std::pair<std::string, std::string>> &tag_value_pairs,
                                    const char *error_label) {
	if (accessions.empty()) {
		throw std::invalid_argument(std::string(error_label) + ": accessions must not be empty");
	}
	if (tags.empty()) {
		throw std::invalid_argument(std::string(error_label) + ": tags must not be empty");
	}
	for (const auto &acc : accessions) {
		ENAParser::ValidateAccession(acc);
	}
	for (const auto &t : tags) {
		ENAParser::ValidateFields(t);
	}

	// `sample_accession` is always first in the column list so the row emitter
	// can locate it unambiguously; remaining tags are deduped in caller order.
	std::ostringstream fields;
	fields << "sample_accession";
	std::unordered_set<std::string> seen_fields = {"sample_accession"};
	for (const auto &t : tags) {
		if (seen_fields.insert(t).second) {
			fields << "," << t;
		}
	}

	std::ostringstream query;
	query << in_column << "%20IN%20%28";
	for (size_t i = 0; i < accessions.size(); i++) {
		if (i > 0) {
			query << "%2C";
		}
		query << "%22" << accessions[i] << "%22";
	}
	query << "%29";
	for (const auto &kv : tag_value_pairs) {
		// Precondition: every pinned tag must appear in `tags` so the response
		// contains a column to un-pivot. Callers that hand-build a pushdown
		// (rather than going through ExtractPushdownPredicates) can violate
		// this; reject early so we never emit a query clause for a field we
		// didn't ask ENA to return.
		if (seen_fields.count(kv.first) == 0) {
			throw std::invalid_argument(std::string(error_label) + ": pinned tag '" + kv.first +
			                            "' is not present in `tags`");
		}
		ENAParser::ValidateFields(kv.first);
		query << "%20AND%20" << kv.first << "%3D%22" << PercentEncodeQueryValue(kv.second) << "%22";
	}

	std::ostringstream url;
	url << ENAParser::PORTAL_BASE << "/search?"
	    << "result=sample&query=" << query.str() << "&fields=" << fields.str() << "&limit=0&format=tsv";
	return url.str();
}

} // namespace

std::string BuildStructuredSearchURL(const std::vector<std::string> &sample_accessions,
                                     const std::vector<std::string> &tags,
                                     const std::vector<std::pair<std::string, std::string>> &tag_value_pairs) {
	return BuildAttributeSearchURL("sample_accession", sample_accessions, tags, tag_value_pairs,
	                               "BuildStructuredSearchURL");
}

std::string BuildStudyDirectSearchURL(const std::vector<std::string> &study_accessions,
                                      const std::vector<std::string> &tags,
                                      const std::vector<std::pair<std::string, std::string>> &tag_value_pairs) {
	return BuildAttributeSearchURL("study_accession", study_accessions, tags, tag_value_pairs,
	                               "BuildStudyDirectSearchURL");
}

// ---- TSV unpivot ----

std::vector<SampleAttribute> UnpivotStructuredTSV(const ENATSVResult &parsed,
                                                  const std::vector<std::string> &requested_tags) {
	std::vector<SampleAttribute> out;
	if (parsed.rows.empty() || parsed.column_names.empty() || requested_tags.empty()) {
		return out;
	}

	// Locate sample_accession + each requested tag's column index
	// (case-insensitive match; ENA returns lowercase but be defensive).
	auto find_col = [&](const std::string &target) -> int {
		auto canonical = ToLower(target);
		for (size_t i = 0; i < parsed.column_names.size(); i++) {
			if (ToLower(parsed.column_names[i]) == canonical) {
				return static_cast<int>(i);
			}
		}
		return -1;
	};

	int acc_idx = find_col("sample_accession");
	if (acc_idx < 0) {
		// No sample_accession column — we can't attribute any row. Return
		// empty and let the caller decide (in practice this is a malformed
		// response and the caller will log/fail).
		return out;
	}

	std::vector<std::pair<std::string, int>> tag_cols;
	tag_cols.reserve(requested_tags.size());
	for (const auto &t : requested_tags) {
		tag_cols.emplace_back(t, find_col(t));
	}

	for (const auto &row : parsed.rows) {
		if (static_cast<size_t>(acc_idx) >= row.size()) {
			continue;
		}
		const auto &acc = row[acc_idx];
		if (acc.empty()) {
			continue;
		}
		for (const auto &tc : tag_cols) {
			if (tc.second < 0) {
				continue; // column absent from response
			}
			if (static_cast<size_t>(tc.second) >= row.size()) {
				continue;
			}
			const auto &val = row[tc.second];
			if (val.empty()) {
				continue; // matches XML path: absent attributes produce no row
			}
			out.push_back(SampleAttribute {acc, tc.first, val});
		}
	}

	return out;
}

} // namespace miint
