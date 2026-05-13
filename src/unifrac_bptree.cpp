#include "unifrac_bptree.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace miint::unifrac {

namespace {

// Iterative pre-order DFS emitting open-paren on entry and close-paren on exit.
// Recursion would overflow the stack on deep phylogenies.
void encode_bp(const miint::NewickTree &tree, uint32_t node, std::vector<unsigned char> &structure,
               std::vector<double> &lengths, std::vector<std::string> &names) {
	struct Frame {
		uint32_t node;
		size_t child_idx;
	};
	std::vector<Frame> stack;
	stack.push_back({node, 0});

	// Emit open for root.
	structure.push_back(1);
	double bl = tree.branch_length(node);
	lengths.push_back(std::isnan(bl) ? 0.0 : bl);
	names.push_back(tree.name(node));

	while (!stack.empty()) {
		auto &frame = stack.back();
		const auto &children = tree.children(frame.node);
		if (frame.child_idx < children.size()) {
			uint32_t child = children[frame.child_idx++];
			structure.push_back(1);
			double cbl = tree.branch_length(child);
			lengths.push_back(std::isnan(cbl) ? 0.0 : cbl);
			names.push_back(tree.name(child));
			stack.push_back({child, 0});
		} else {
			structure.push_back(0);
			lengths.push_back(0.0);
			names.push_back("");
			stack.pop_back();
		}
	}
}

} // namespace

UnifracBptreeView UnifracBptreeView::FromNewickTree(const miint::NewickTree &tree) {
	UnifracBptreeView view;
	const size_t n_nodes = tree.num_nodes();
	view.structure_bytes_.reserve(2 * n_nodes);
	view.lengths_.reserve(2 * n_nodes);
	view.names_storage_.reserve(2 * n_nodes);

	encode_bp(tree, tree.root(), view.structure_bytes_, view.lengths_, view.names_storage_);

	view.wire_bptree_ptrs();
	return view;
}

void UnifracBptreeView::wire_bptree_ptrs() {
	name_ptrs_.clear();
	name_ptrs_.reserve(names_storage_.size());
	for (auto &n : names_storage_) {
		name_ptrs_.push_back(n.data());
	}
	bptree_.structure = reinterpret_cast<bool *>(structure_bytes_.data());
	bptree_.lengths = lengths_.data();
	bptree_.names = name_ptrs_.data();
	bptree_.n_parens = static_cast<int>(structure_bytes_.size());
}

UnifracBptreeView::UnifracBptreeView(UnifracBptreeView &&other) noexcept
    : structure_bytes_(std::move(other.structure_bytes_)), lengths_(std::move(other.lengths_)),
      names_storage_(std::move(other.names_storage_)) {
	wire_bptree_ptrs();
}

UnifracBptreeView &UnifracBptreeView::operator=(UnifracBptreeView &&other) noexcept {
	if (this != &other) {
		structure_bytes_ = std::move(other.structure_bytes_);
		lengths_ = std::move(other.lengths_);
		names_storage_ = std::move(other.names_storage_);
		wire_bptree_ptrs();
	}
	return *this;
}

void ValidateTreeCoversFeatures(const miint::NewickTree &tree, const std::vector<std::string> &feature_ids) {
	auto tip_names = tree.tip_names();
	std::unordered_set<std::string> tips(tip_names.begin(), tip_names.end());
	for (const auto &fid : feature_ids) {
		if (tips.find(fid) == tips.end()) {
			throw std::invalid_argument("Tree does not contain feature_id '" + fid + "' (tree has " +
			                            std::to_string(tips.size()) + " tips)");
		}
	}
}

} // namespace miint::unifrac
