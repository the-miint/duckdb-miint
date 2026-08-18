#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "NewickTree.hpp"
#include "api.hpp" // libssu (support_bptree_t)

namespace miint::unifrac {

// Wraps libssu's support_bptree_t with backing storage:
// preorder-open / postorder-close BP encoding of a NewickTree, with names
// populated at open-paren positions and 0-length at close positions.
// Move-only; the support_bptree_t pointer's fields are valid for the view's
// lifetime.
class UnifracBptreeView {
public:
	static UnifracBptreeView FromNewickTree(const miint::NewickTree &tree);

	const support_bptree_t *support_bptree() const {
		return &bptree_;
	}

	UnifracBptreeView(UnifracBptreeView &&) noexcept;
	UnifracBptreeView &operator=(UnifracBptreeView &&) noexcept;
	UnifracBptreeView(const UnifracBptreeView &) = delete;
	UnifracBptreeView &operator=(const UnifracBptreeView &) = delete;
	~UnifracBptreeView() = default;

private:
	UnifracBptreeView() = default;
	void wire_bptree_ptrs();

	// sizeof(bool) == 1 universally; we back the BP structure as bytes and
	// reinterpret_cast to bool* for libssu. std::vector<bool> is unusable
	// here (packed, no contiguous bool* access).
	std::vector<unsigned char> structure_bytes_;
	std::vector<double> lengths_;
	std::vector<std::string> names_storage_;
	std::vector<char *> name_ptrs_;
	support_bptree_t bptree_ {};
};

// Throws std::invalid_argument naming the first missing feature_id when any
// feature is not a tip of `tree`. Internal node names are not tips and never
// satisfy this check. Extra tips in the tree (not in feature_ids) are fine.
void ValidateTreeCoversFeatures(const miint::NewickTree &tree, const std::vector<std::string> &feature_ids);

} // namespace miint::unifrac
