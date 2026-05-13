#include "unifrac_metadata.hpp"

#include <algorithm>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace miint::unifrac {

namespace {

NamedGrouping factorize_one(const std::string &variable, const std::vector<std::string> &ordered_sample_ids,
                            const std::unordered_map<std::string, const MetadataRow *> &sample_to_row) {
	NamedGrouping g;
	g.variable = variable;
	g.labels.reserve(ordered_sample_ids.size());

	// Label by first appearance in ordered_sample_ids — independent of input
	// row order. Same partition across variables ⇒ same labels.
	std::unordered_map<std::string, uint32_t> value_to_label;
	uint32_t next_label = 0;
	for (const auto &sid : ordered_sample_ids) {
		auto it = sample_to_row.find(sid);
		if (it == sample_to_row.end()) {
			throw std::invalid_argument("Metadata is missing variable '" + variable + "' for sample '" + sid + "'");
		}
		const std::string &value = it->second->value;
		auto label_it = value_to_label.find(value);
		uint32_t label;
		if (label_it == value_to_label.end()) {
			label = next_label++;
			value_to_label.emplace(value, label);
		} else {
			label = label_it->second;
		}
		g.labels.push_back(label);
	}
	g.n_groups = next_label;
	return g;
}

} // namespace

std::vector<NamedGrouping> BuildGroupings(const std::vector<MetadataRow> &rows,
                                          const std::vector<std::string> &ordered_sample_ids,
                                          const std::vector<std::string> &variables) {
	if (ordered_sample_ids.empty()) {
		throw std::invalid_argument("BuildGroupings: ordered_sample_ids is empty");
	}

	// Restrict to samples in the analysis set up front — metadata may
	// oversample, and "extra" rows should not error out.
	std::unordered_set<std::string> sample_set(ordered_sample_ids.begin(), ordered_sample_ids.end());

	// Index relevant rows by (variable, sample_id). Duplicate (variable,
	// sample) pairs are almost always user error (a botched join, accidental
	// cartesian, two metadata files concatenated). Silent last-write-wins
	// would corrupt PERMANOVA's F-stat; fail loud naming both keys.
	std::unordered_map<std::string, std::unordered_map<std::string, const MetadataRow *>> by_variable;
	for (const auto &r : rows) {
		if (sample_set.count(r.sample_id) == 0) {
			continue;
		}
		auto &sample_map = by_variable[r.variable];
		auto [it, inserted] = sample_map.emplace(r.sample_id, &r);
		if (!inserted && it->second->value != r.value) {
			throw std::invalid_argument("Metadata has conflicting values for sample '" + r.sample_id +
			                            "' on variable '" + r.variable + "': '" + it->second->value + "' vs '" +
			                            r.value + "'");
		}
	}

	// Determine the list of variables to factorize.
	std::vector<std::string> target_variables;
	if (variables.empty()) {
		// Lexicographic order for stable output cardinality and column order.
		std::set<std::string> distinct;
		for (const auto &r : rows) {
			distinct.insert(r.variable);
		}
		target_variables.assign(distinct.begin(), distinct.end());
	} else {
		target_variables = variables;
	}

	std::vector<NamedGrouping> result;
	result.reserve(target_variables.size());
	for (const auto &var : target_variables) {
		auto it = by_variable.find(var);
		if (it == by_variable.end()) {
			throw std::invalid_argument("Metadata does not contain variable '" + var + "'");
		}
		result.push_back(factorize_one(var, ordered_sample_ids, it->second));
	}
	return result;
}

} // namespace miint::unifrac
