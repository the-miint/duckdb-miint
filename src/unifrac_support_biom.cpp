#include "unifrac_support_biom.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace miint::unifrac {

namespace {

// Stable ID dictionary: sort the DISTINCT ids lexicographically, return both the
// ordered list and a name->index map.
//
// Takes the already-deduplicated set rather than the raw per-row ids on purpose.
// The obvious shape — copy every row's id into a vector, sort, unique — sorts nnz
// strings to discover a few thousand distinct ones: on a 3200-sample block that is
// ~2.4M strings sorted to find ~60k, twice (once for samples, once for features),
// which measured as the bulk of FromCoo's cost. Deduplicating during the single
// pass the caller already makes over the rows leaves only the survivors to sort.
std::pair<std::vector<std::string>, std::unordered_map<std::string, uint32_t>>
build_dictionary(const std::unordered_set<std::string> &distinct_ids) {
	std::vector<std::string> unique(distinct_ids.begin(), distinct_ids.end());
	std::sort(unique.begin(), unique.end());
	std::unordered_map<std::string, uint32_t> index;
	index.reserve(unique.size());
	for (uint32_t i = 0; i < unique.size(); ++i) {
		index.emplace(unique[i], i);
	}
	return {std::move(unique), std::move(index)};
}

} // namespace

UnifracSupportBiomView UnifracSupportBiomView::FromCoo(std::vector<CooRow> rows) {
	if (rows.empty()) {
		throw std::invalid_argument("UnifracSupportBiomView: input feature-table is empty");
	}

	// One pass to collect the distinct ids of both dimensions. An id already seen
	// costs a hash and a compare, not a copy, so this does not grow with how
	// duplicated the input is — unlike materializing one id per row per dimension.
	std::unordered_set<std::string> distinct_samples;
	std::unordered_set<std::string> distinct_features;
	for (const auto &r : rows) {
		distinct_samples.insert(r.sample_id);
		distinct_features.insert(r.feature_id);
	}
	auto [sample_ids, sample_index] = build_dictionary(distinct_samples);
	auto [feature_ids, feature_index] = build_dictionary(distinct_features);
	// The index maps carry everything needed from here; release the sets before
	// the entries vector below allocates.
	distinct_samples = {};
	distinct_features = {};

	// CSR layout (obs-major) per libssu's documented contract: outer loop is
	// over n_obs, inner indices are sample column indices. Sort by
	// (feature_idx, sample_idx) and aggregate duplicates.
	struct Entry {
		uint32_t feature;
		uint32_t sample;
		double value;
	};
	std::vector<Entry> entries;
	entries.reserve(rows.size());
	for (const auto &r : rows) {
		entries.push_back({feature_index.at(r.feature_id), sample_index.at(r.sample_id), r.count});
	}
	std::sort(entries.begin(), entries.end(), [](const Entry &a, const Entry &b) {
		if (a.feature != b.feature) {
			return a.feature < b.feature;
		}
		return a.sample < b.sample;
	});

	std::vector<uint32_t> indices;
	std::vector<uint32_t> indptr(feature_ids.size() + 1, 0);
	std::vector<double> data;
	indices.reserve(entries.size());
	data.reserve(entries.size());

	uint32_t cursor_feature = 0;
	for (size_t i = 0; i < entries.size();) {
		size_t j = i;
		double sum = 0.0;
		while (j < entries.size() && entries[j].feature == entries[i].feature &&
		       entries[j].sample == entries[i].sample) {
			sum += entries[j].value;
			++j;
		}
		if (sum != 0.0) {
			while (cursor_feature < entries[i].feature) {
				indptr[cursor_feature + 1] = static_cast<uint32_t>(data.size());
				++cursor_feature;
			}
			indices.push_back(entries[i].sample);
			data.push_back(sum);
		}
		i = j;
	}
	while (cursor_feature < feature_ids.size()) {
		indptr[cursor_feature + 1] = static_cast<uint32_t>(data.size());
		++cursor_feature;
	}

	UnifracSupportBiomView view;
	view.sample_ids_ = std::move(sample_ids);
	view.feature_ids_ = std::move(feature_ids);
	view.indices_ = std::move(indices);
	view.indptr_ = std::move(indptr);
	view.data_ = std::move(data);
	view.wire_biom_ptrs();
	return view;
}

void UnifracSupportBiomView::wire_biom_ptrs() {
	sample_id_ptrs_.clear();
	feature_id_ptrs_.clear();
	sample_id_ptrs_.reserve(sample_ids_.size());
	feature_id_ptrs_.reserve(feature_ids_.size());
	for (auto &s : sample_ids_) {
		sample_id_ptrs_.push_back(s.data());
	}
	for (auto &f : feature_ids_) {
		feature_id_ptrs_.push_back(f.data());
	}
	biom_.sample_ids = sample_id_ptrs_.data();
	biom_.obs_ids = feature_id_ptrs_.data();
	biom_.indices = indices_.data();
	biom_.indptr = indptr_.data();
	biom_.data = data_.data();
	biom_.n_samples = static_cast<int>(sample_ids_.size());
	biom_.n_obs = static_cast<int>(feature_ids_.size());
	biom_.nnz = static_cast<int>(indices_.size());
}

UnifracSupportBiomView::UnifracSupportBiomView(UnifracSupportBiomView &&other) noexcept
    : sample_ids_(std::move(other.sample_ids_)), feature_ids_(std::move(other.feature_ids_)),
      indices_(std::move(other.indices_)), indptr_(std::move(other.indptr_)), data_(std::move(other.data_)) {
	wire_biom_ptrs();
}

UnifracSupportBiomView &UnifracSupportBiomView::operator=(UnifracSupportBiomView &&other) noexcept {
	if (this != &other) {
		sample_ids_ = std::move(other.sample_ids_);
		feature_ids_ = std::move(other.feature_ids_);
		indices_ = std::move(other.indices_);
		indptr_ = std::move(other.indptr_);
		data_ = std::move(other.data_);
		wire_biom_ptrs();
	}
	return *this;
}

} // namespace miint::unifrac
