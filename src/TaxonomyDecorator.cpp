#include "TaxonomyDecorator.hpp"

#include <utility>

namespace miint {

// Phase 1 (red) stubs: every function returns a neutral value so the tests build and fail on
// their assertions. Phase 2 replaces them.

uint32_t F1Score::Threshold(uint32_t, uint32_t) {
	return 0;
}

std::strong_ordering F1Score::Compare(uint32_t, uint32_t, uint32_t, uint32_t, uint32_t) {
	return std::strong_ordering::equal;
}

double F1Score::ScoreForReport(uint32_t, uint32_t, uint32_t) {
	return 0.0;
}

CladeCountIndex::CladeCountIndex(TreeIndex tree, const TaxonomyStore &, const std::vector<std::optional<size_t>> &)
    : tree_(std::move(tree)) {
}

const TreeIndex &CladeCountIndex::Tree() const {
	return tree_;
}

uint32_t CladeCountIndex::CladeCount(uint32_t, size_t, NameId) const {
	return 0;
}

uint32_t CladeCountIndex::NumTips(uint32_t) const {
	return 0;
}

uint32_t CladeCountIndex::LeafCount(uint32_t) const {
	return 0;
}

uint32_t CladeCountIndex::Total(size_t, NameId) const {
	return 0;
}

} // namespace miint
