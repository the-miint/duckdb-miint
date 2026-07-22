#include "NewickTree.hpp"
#include <algorithm>
#include <array>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <stack>

// Eigen (MPL-2.0, header-only; EIGEN_MPL2_ONLY set in CMakeLists) backs the SYM/ARD
// Mk-ML P(t) = exp(Q*t): SYM via a symmetric eigendecomposition (SelfAdjointEigenSolver),
// ARD via a general matrix exponential (MatrixFunctions, scaling-and-squaring Pade, Higham
// 2005). See ancestral_ml_sym / ancestral_ml_ard below.
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace miint {

// Maximum number of nodes supported (uint32_t limit)
static constexpr size_t MAX_NODES = static_cast<size_t>(UINT32_MAX);

// ============================================================================
// Parser implementation
// ============================================================================

class NewickParser {
public:
	explicit NewickParser(std::string_view input) : input_(input), pos_(0) {
	}

	NewickTree parse() {
		skip_whitespace_and_comments();

		if (pos_ >= input_.size() || input_[pos_] == ';') {
			// Empty tree (just ";" or whitespace/comments followed by ";")
			// This is valid - represents an unnamed root with no children
			NewickTree tree;
			tree.root_ = 0;
			tree.nodes_.emplace_back("", std::numeric_limits<double>::quiet_NaN(), std::nullopt);
			if (pos_ < input_.size() && input_[pos_] == ';') {
				return tree;
			}
			throw std::runtime_error("Cannot parse empty Newick string");
		}

		NewickTree tree;
		tree.root_ = parse_node(tree);

		skip_whitespace_and_comments();

		if (pos_ >= input_.size() || input_[pos_] != ';') {
			throw std::runtime_error("Missing semicolon at end of Newick string");
		}

		return tree;
	}

private:
	std::string_view input_;
	size_t pos_;

	char peek() const {
		return pos_ < input_.size() ? input_[pos_] : '\0';
	}

	char consume() {
		return pos_ < input_.size() ? input_[pos_++] : '\0';
	}

	void skip_whitespace_and_comments() {
		while (pos_ < input_.size()) {
			char c = input_[pos_];
			if (std::isspace(static_cast<unsigned char>(c))) {
				++pos_;
			} else if (c == '[') {
				skip_comment();
			} else {
				break;
			}
		}
	}

	void skip_comment() {
		// Skip '[' and everything until matching ']', handling nesting
		size_t comment_start = pos_;
		consume(); // consume '['
		int depth = 1;
		while (pos_ < input_.size() && depth > 0) {
			char c = input_[pos_++];
			if (c == '[') {
				++depth;
			} else if (c == ']') {
				--depth;
			}
		}
		if (depth > 0) {
			throw std::runtime_error("Unclosed comment starting at position " + std::to_string(comment_start));
		}
	}

	// Keep simple whitespace skip for places where comments aren't expected
	void skip_whitespace() {
		while (pos_ < input_.size() && std::isspace(static_cast<unsigned char>(input_[pos_]))) {
			++pos_;
		}
	}

	uint32_t parse_node(NewickTree &tree) {
		skip_whitespace_and_comments();

		std::vector<uint32_t> children;

		// Check for children (subtree)
		if (peek() == '(') {
			consume(); // consume '('
			children = parse_children(tree);

			skip_whitespace_and_comments();
			if (peek() != ')') {
				throw std::runtime_error("Unmatched opening parenthesis in Newick string");
			}
			consume(); // consume ')'
		}

		skip_whitespace_and_comments();

		// Parse node label
		std::string name = parse_label();

		// Parse branch length (optional)
		double branch_length = std::numeric_limits<double>::quiet_NaN();
		skip_whitespace_and_comments();
		if (peek() == ':') {
			consume(); // consume ':'
			branch_length = parse_branch_length();
		}

		// Parse edge identifier (optional, jplace format)
		std::optional<int64_t> edge_id;
		skip_whitespace_and_comments();
		if (peek() == '{') {
			edge_id = parse_edge_id();
		}

		// Check for overflow before creating node
		if (tree.nodes_.size() >= MAX_NODES) {
			throw std::runtime_error("Tree too large: exceeds maximum of " + std::to_string(MAX_NODES) + " nodes");
		}

		// Create node
		uint32_t node_idx = static_cast<uint32_t>(tree.nodes_.size());
		tree.nodes_.emplace_back(std::move(name), branch_length, edge_id);

		// Set up parent-child relationships
		for (auto child : children) {
			tree.nodes_[child].parent = node_idx;
			tree.nodes_[node_idx].children.push_back(child);
		}

		return node_idx;
	}

	std::vector<uint32_t> parse_children(NewickTree &tree) {
		std::vector<uint32_t> children;

		while (true) {
			skip_whitespace_and_comments();
			children.push_back(parse_node(tree));

			skip_whitespace_and_comments();
			if (peek() == ',') {
				consume(); // consume ','
			} else {
				break;
			}
		}

		return children;
	}

	std::string parse_label() {
		skip_whitespace_and_comments();

		char c = peek();

		// Quoted label
		if (c == '\'' || c == '"') {
			return parse_quoted_label(c);
		}

		// Unquoted label - ends at special characters
		// Use substr for O(1) extraction instead of character-by-character O(n²)
		size_t start = pos_;
		while (pos_ < input_.size()) {
			c = input_[pos_];
			if (c == '(' || c == ')' || c == ',' || c == ':' || c == ';' || c == '{' || c == '[' ||
			    std::isspace(static_cast<unsigned char>(c))) {
				break;
			}
			++pos_;
		}

		return std::string(input_.substr(start, pos_ - start));
	}

	std::string parse_quoted_label(char quote_char) {
		consume(); // consume opening quote

		// First pass: find end and check if escaping needed
		size_t start = pos_;
		bool has_escapes = false;
		size_t label_len = 0;

		while (pos_ < input_.size()) {
			char c = input_[pos_];
			if (c == quote_char) {
				if (pos_ + 1 < input_.size() && input_[pos_ + 1] == quote_char) {
					// Escaped quote
					has_escapes = true;
					pos_ += 2;
					label_len += 1; // Only one quote in output
				} else {
					// End of quoted string
					break;
				}
			} else {
				++pos_;
				++label_len;
			}
		}

		if (pos_ >= input_.size()) {
			throw std::runtime_error("Unclosed quote in Newick label");
		}

		size_t end = pos_;
		consume(); // consume closing quote

		// If no escapes, just return substring
		if (!has_escapes) {
			return std::string(input_.substr(start, end - start));
		}

		// Otherwise, build unescaped string
		std::string label;
		label.reserve(label_len);
		for (size_t i = start; i < end; ++i) {
			char c = input_[i];
			if (c == quote_char && i + 1 < end && input_[i + 1] == quote_char) {
				label += quote_char;
				++i; // skip the second quote
			} else {
				label += c;
			}
		}

		return label;
	}

	double parse_branch_length() {
		skip_whitespace();

		size_t start = pos_;
		// Find end of token (stop at Newick structural characters)
		while (pos_ < input_.size()) {
			char c = input_[pos_];
			if (c == '(' || c == ')' || c == ',' || c == ';' || c == '{' || c == '[' ||
			    std::isspace(static_cast<unsigned char>(c))) {
				break;
			}
			++pos_;
		}

		if (pos_ == start) {
			throw std::runtime_error("Invalid branch length: expected number after ':'");
		}

		std::string_view num_str = input_.substr(start, pos_ - start);

		// Use strtod for platform compatibility (std::from_chars for float not available on macOS yet)
		// Need null-terminated string for strtod
		std::string num_cstr(num_str);
		char *endptr;
		errno = 0;
		double value = std::strtod(num_cstr.c_str(), &endptr);

		if (errno == ERANGE) {
			throw std::runtime_error("Invalid branch length: '" + std::string(num_str) + "' (out of range)");
		}
		if (endptr == num_cstr.c_str()) {
			throw std::runtime_error("Invalid branch length: '" + std::string(num_str) + "'");
		}
		if (endptr != num_cstr.c_str() + num_cstr.size()) {
			throw std::runtime_error("Invalid branch length: unexpected characters in '" + std::string(num_str) + "'");
		}

		return value;
	}

	int64_t parse_edge_id() {
		if (consume() != '{') {
			throw std::runtime_error("Expected '{' for edge identifier");
		}

		skip_whitespace();

		size_t start = pos_;
		// Find end of token (stop at '}' or other structural chars)
		while (pos_ < input_.size()) {
			char c = input_[pos_];
			if (c == '}' || c == '(' || c == ')' || c == ',' || c == ';' || c == '[' ||
			    std::isspace(static_cast<unsigned char>(c))) {
				break;
			}
			++pos_;
		}

		if (pos_ == start) {
			throw std::runtime_error("Invalid edge identifier: expected integer");
		}

		std::string_view num_str = input_.substr(start, pos_ - start);
		int64_t value;

		auto [ptr, ec] = std::from_chars(num_str.data(), num_str.data() + num_str.size(), value);

		if (ec != std::errc()) {
			throw std::runtime_error("Invalid edge identifier: '" + std::string(num_str) + "'");
		}
		if (ptr != num_str.data() + num_str.size()) {
			throw std::runtime_error("Invalid edge identifier: unexpected characters in '" + std::string(num_str) +
			                         "'");
		}

		skip_whitespace();

		if (pos_ >= input_.size() || input_[pos_] != '}') {
			throw std::runtime_error("Unclosed brace in edge identifier");
		}
		consume(); // consume '}'

		return value;
	}
};

NewickTree NewickTree::parse(std::string_view newick) {
	NewickParser parser(newick);
	return parser.parse();
}

// ============================================================================
// Serialization (iterative to avoid stack overflow on deep trees)
// ============================================================================

std::string NewickTree::to_newick() const {
	if (nodes_.empty()) {
		return ";";
	}

	std::string result;
	result.reserve(nodes_.size() * 40); // Rough estimate (name + branch length + edge_id + delimiters)

	// Iterative serialization using explicit stack
	// Each entry: (node_index, child_index, state)
	// state: 0 = start, 1 = processing children, 2 = done with children
	struct StackEntry {
		uint32_t node;
		size_t child_idx;
		int state;
	};

	std::stack<StackEntry> stack;
	stack.push({root_, 0, 0});

	while (!stack.empty()) {
		auto &entry = stack.top();
		const auto &n = nodes_[entry.node];

		if (entry.state == 0) {
			// Starting this node
			if (!n.children.empty()) {
				result += '(';
				entry.state = 1;
				entry.child_idx = 0;
			} else {
				entry.state = 2; // No children, go straight to serializing node info
			}
		} else if (entry.state == 1) {
			// Processing children
			if (entry.child_idx < n.children.size()) {
				if (entry.child_idx > 0) {
					result += ',';
				}
				// Push child onto stack
				stack.push({n.children[entry.child_idx], 0, 0});
				entry.child_idx++;
			} else {
				// Done with all children
				result += ')';
				entry.state = 2;
			}
		} else {
			// state == 2: Serialize node info and pop
			// Name - quote if contains special characters
			if (!n.name.empty()) {
				bool needs_quote = false;
				for (char c : n.name) {
					if (c == '(' || c == ')' || c == ',' || c == ':' || c == ';' || c == '{' || c == '}' || c == '\'' ||
					    c == '"' || c == '[' || c == ']' || std::isspace(static_cast<unsigned char>(c))) {
						needs_quote = true;
						break;
					}
				}

				if (needs_quote) {
					result += '\'';
					for (char c : n.name) {
						if (c == '\'') {
							result += "''"; // escape single quote
						} else {
							result += c;
						}
					}
					result += '\'';
				} else {
					result += n.name;
				}
			}

			// Branch length - use snprintf for platform compatibility
			if (!std::isnan(n.branch_length)) {
				result += ':';
				// Use snprintf for portability (std::to_chars for float not available on macOS yet)
				std::array<char, 32> buf;
				int len = std::snprintf(buf.data(), buf.size(), "%.15g", n.branch_length);
				if (len > 0 && len < static_cast<int>(buf.size())) {
					result.append(buf.data(), len);
				}
			}

			// Edge identifier
			if (n.edge_id.has_value()) {
				result += '{';
				result += std::to_string(n.edge_id.value());
				result += '}';
			}

			stack.pop();
		}
	}

	result += ';';
	return result;
}

// ============================================================================
// Statistics
// ============================================================================

size_t NewickTree::num_tips() const {
	size_t count = 0;
	for (const auto &node : nodes_) {
		if (node.children.empty()) {
			++count;
		}
	}
	return count;
}

// ============================================================================
// Traversal
// ============================================================================

std::vector<uint32_t> NewickTree::postorder() const {
	std::vector<uint32_t> result;
	result.reserve(nodes_.size());

	if (nodes_.empty()) {
		return result;
	}

	// Iterative post-order using two stacks
	std::stack<uint32_t> s1, s2;
	s1.push(root_);

	while (!s1.empty()) {
		uint32_t node = s1.top();
		s1.pop();
		s2.push(node);

		for (auto child : nodes_[node].children) {
			s1.push(child);
		}
	}

	while (!s2.empty()) {
		result.push_back(s2.top());
		s2.pop();
	}

	return result;
}

std::vector<uint32_t> NewickTree::preorder() const {
	std::vector<uint32_t> result;
	result.reserve(nodes_.size());

	if (nodes_.empty()) {
		return result;
	}

	std::stack<uint32_t> s;
	s.push(root_);

	while (!s.empty()) {
		uint32_t node = s.top();
		s.pop();
		result.push_back(node);

		// Push children in reverse order so they come out in correct order
		const auto &children = nodes_[node].children;
		for (auto it = children.rbegin(); it != children.rend(); ++it) {
			s.push(*it);
		}
	}

	return result;
}

// ============================================================================
// Queries
// ============================================================================

std::optional<uint32_t> NewickTree::find_node_by_name(std::string_view name) const {
	for (uint32_t i = 0; i < nodes_.size(); ++i) {
		if (nodes_[i].name == name) {
			return i;
		}
	}
	return std::nullopt;
}

std::optional<uint32_t> NewickTree::find_node_by_edge_id(int64_t edge_id) const {
	for (uint32_t i = 0; i < nodes_.size(); ++i) {
		if (nodes_[i].edge_id.has_value() && nodes_[i].edge_id.value() == edge_id) {
			return i;
		}
	}
	return std::nullopt;
}

// ============================================================================
// Build from node data
// ============================================================================

NewickTree NewickTree::build(const std::vector<NodeInput> &nodes) {
	if (nodes.empty()) {
		throw std::runtime_error("Cannot build tree from empty node list");
	}

	if (nodes.size() > MAX_NODES) {
		throw std::runtime_error("Too many nodes: " + std::to_string(nodes.size()) + " exceeds maximum of " +
		                         std::to_string(MAX_NODES));
	}

	// Build mapping from input node_id to index in input vector
	std::unordered_map<int64_t, size_t> id_to_input_idx;
	for (size_t i = 0; i < nodes.size(); ++i) {
		auto [it, inserted] = id_to_input_idx.emplace(nodes[i].node_id, i);
		if (!inserted) {
			throw std::runtime_error("Duplicate node_id: " + std::to_string(nodes[i].node_id));
		}
	}

	// Find root(s) and validate parent references
	std::vector<size_t> roots;
	for (size_t i = 0; i < nodes.size(); ++i) {
		const auto &node = nodes[i];
		if (!node.parent_id.has_value()) {
			roots.push_back(i);
		} else {
			if (id_to_input_idx.find(node.parent_id.value()) == id_to_input_idx.end()) {
				throw std::runtime_error("Node " + std::to_string(node.node_id) + " references non-existent parent " +
				                         std::to_string(node.parent_id.value()));
			}
		}
	}

	if (roots.empty()) {
		throw std::runtime_error("No root found (no node with null parent_id)");
	}
	if (roots.size() > 1) {
		throw std::runtime_error("Multiple roots found (" + std::to_string(roots.size()) +
		                         " nodes with null parent_id)");
	}

	// Cycle detection using DFS - verify all nodes are reachable from root
	// and no back edges exist
	std::vector<bool> visited(nodes.size(), false);
	std::vector<bool> in_stack(nodes.size(), false);

	// Build children map for traversal
	std::unordered_map<int64_t, std::vector<size_t>> children_map;
	for (size_t i = 0; i < nodes.size(); ++i) {
		if (nodes[i].parent_id.has_value()) {
			children_map[nodes[i].parent_id.value()].push_back(i);
		}
	}

	// Iterative DFS for cycle detection
	std::vector<std::pair<size_t, size_t>> stack; // (node_idx, child_iterator)
	stack.emplace_back(roots[0], 0);
	visited[roots[0]] = true;
	in_stack[roots[0]] = true;

	while (!stack.empty()) {
		auto &[node_idx, child_iter] = stack.back();
		const auto &node = nodes[node_idx];
		auto &children = children_map[node.node_id];

		if (child_iter < children.size()) {
			size_t child_idx = children[child_iter];
			child_iter++;

			if (in_stack[child_idx]) {
				throw std::runtime_error("Cycle detected involving node " + std::to_string(nodes[child_idx].node_id));
			}

			if (!visited[child_idx]) {
				visited[child_idx] = true;
				in_stack[child_idx] = true;
				stack.emplace_back(child_idx, 0);
			}
		} else {
			in_stack[node_idx] = false;
			stack.pop_back();
		}
	}

	// Check all nodes are reachable from root
	for (size_t i = 0; i < nodes.size(); ++i) {
		if (!visited[i]) {
			throw std::runtime_error("Node " + std::to_string(nodes[i].node_id) +
			                         " is not reachable from root (disconnected tree)");
		}
	}

	// Build the tree
	NewickTree tree;
	tree.nodes_.reserve(nodes.size());

	// Add all nodes (input order preserved for node indices)
	std::unordered_map<int64_t, uint32_t> id_to_tree_idx;
	for (size_t i = 0; i < nodes.size(); ++i) {
		const auto &node = nodes[i];
		uint32_t tree_idx = tree.add_node(node.name, node.branch_length, node.edge_id);
		id_to_tree_idx[node.node_id] = tree_idx;
	}

	// Set parent relationships and identify root
	for (size_t i = 0; i < nodes.size(); ++i) {
		const auto &node = nodes[i];
		uint32_t tree_idx = id_to_tree_idx[node.node_id];

		if (node.parent_id.has_value()) {
			uint32_t parent_tree_idx = id_to_tree_idx[node.parent_id.value()];
			tree.set_parent(tree_idx, parent_tree_idx);
		} else {
			// This is the root
			tree.root_ = tree_idx;
		}
	}

	return tree;
}

std::unordered_map<int64_t, uint32_t> NewickTree::build_edge_index() const {
	std::unordered_map<int64_t, uint32_t> index;
	for (uint32_t i = 0; i < nodes_.size(); ++i) {
		if (nodes_[i].edge_id.has_value()) {
			index[nodes_[i].edge_id.value()] = i;
		}
	}
	return index;
}

// ============================================================================
// Distance calculations
// ============================================================================

std::vector<uint32_t> NewickTree::tips() const {
	std::vector<uint32_t> result;
	for (uint32_t i = 0; i < nodes_.size(); ++i) {
		if (nodes_[i].children.empty()) {
			result.push_back(i);
		}
	}
	return result;
}

std::vector<std::string> NewickTree::tip_names() const {
	std::vector<std::string> result;
	for (const auto &node : nodes_) {
		if (node.children.empty()) {
			result.push_back(node.name);
		}
	}
	return result;
}

double NewickTree::distance_to_root(uint32_t node) const {
	if (node >= nodes_.size()) {
		throw std::out_of_range("Invalid node index: " + std::to_string(node));
	}

	double dist = 0.0;
	uint32_t current = node;
	while (current != root_) {
		double bl = nodes_[current].branch_length;
		if (!std::isnan(bl)) {
			dist += bl;
		}
		current = nodes_[current].parent;
		if (current == NO_PARENT) {
			break; // Reached root or disconnected node
		}
	}
	return dist;
}

uint32_t NewickTree::find_lca(uint32_t a, uint32_t b) const {
	if (a >= nodes_.size()) {
		throw std::out_of_range("Invalid node index: " + std::to_string(a));
	}
	if (b >= nodes_.size()) {
		throw std::out_of_range("Invalid node index: " + std::to_string(b));
	}

	// Collect ancestors of a (including a itself)
	std::vector<bool> is_ancestor_of_a(nodes_.size(), false);
	uint32_t current = a;
	while (current != NO_PARENT) {
		is_ancestor_of_a[current] = true;
		current = nodes_[current].parent;
	}

	// Find first ancestor of b that is also ancestor of a
	current = b;
	while (current != NO_PARENT) {
		if (is_ancestor_of_a[current]) {
			return current;
		}
		current = nodes_[current].parent;
	}

	// Should not happen in a valid tree
	throw std::runtime_error("Nodes have no common ancestor - tree may be disconnected");
}

double NewickTree::pairwise_distance(uint32_t a, uint32_t b) const {
	if (a == b) {
		return 0.0;
	}

	uint32_t lca = find_lca(a, b);

	// Sum distances from a to LCA
	double dist = 0.0;
	uint32_t current = a;
	while (current != lca) {
		double bl = nodes_[current].branch_length;
		if (!std::isnan(bl)) {
			dist += bl;
		}
		current = nodes_[current].parent;
	}

	// Sum distances from b to LCA
	current = b;
	while (current != lca) {
		double bl = nodes_[current].branch_length;
		if (!std::isnan(bl)) {
			dist += bl;
		}
		current = nodes_[current].parent;
	}

	return dist;
}

// ============================================================================
// Phylogenetic placement
// ============================================================================

void NewickTree::insert_fully_resolved(const std::vector<Placement> &placements) {
	if (placements.empty()) {
		return;
	}

	// Step 1: Build edge_id -> node_id index
	auto edge_index = build_edge_index();

	// Step 2: Validate ALL placements upfront before any processing
	// This ensures we catch and report all invalid data, not just the "winning" placements
	for (const auto &p : placements) {
		if (edge_index.find(p.edge_id) == edge_index.end()) {
			throw std::runtime_error("Unknown edge_id " + std::to_string(p.edge_id) + " for fragment '" +
			                         p.fragment_id + "'");
		}
		if (p.distal_length < 0) {
			throw std::runtime_error("Negative distal_length " + std::to_string(p.distal_length) + " for fragment '" +
			                         p.fragment_id + "'");
		}
		if (p.pendant_length < 0) {
			throw std::runtime_error("Negative pendant_length " + std::to_string(p.pendant_length) + " for fragment '" +
			                         p.fragment_id + "'");
		}
		// Validate distal_length against edge length
		uint32_t edge_node = edge_index.at(p.edge_id);
		double edge_length = nodes_[edge_node].branch_length;
		if (!std::isnan(edge_length) && p.distal_length > edge_length) {
			throw std::runtime_error("distal_length " + std::to_string(p.distal_length) + " exceeds edge length " +
			                         std::to_string(edge_length) + " for fragment '" + p.fragment_id + "'");
		}
	}

	// Step 3: Deduplicate placements by fragment_id
	// Keep the one with highest like_weight_ratio, then lowest pendant_length
	// Note: Using epsilon comparison for floating-point like_weight_ratio
	constexpr double EPSILON = 1e-9;
	std::unordered_map<std::string, const Placement *> best_placement;
	for (const auto &p : placements) {
		auto it = best_placement.find(p.fragment_id);
		if (it == best_placement.end()) {
			best_placement[p.fragment_id] = &p;
		} else {
			const Placement *existing = it->second;
			double diff = p.like_weight_ratio - existing->like_weight_ratio;
			// Prefer higher like_weight_ratio (with epsilon tolerance)
			if (diff > EPSILON) {
				it->second = &p;
			} else if (std::abs(diff) <= EPSILON) {
				// Tiebreaker: prefer lower pendant_length
				if (p.pendant_length < existing->pendant_length) {
					it->second = &p;
				}
			}
		}
	}

	// Step 4: Check node limit before any allocation
	// Each unique placement creates 2 nodes: internal + fragment
	size_t new_nodes_needed = best_placement.size() * 2;
	if (nodes_.size() + new_nodes_needed > MAX_NODES) {
		throw std::runtime_error("Too many placements: would create " + std::to_string(new_nodes_needed) +
		                         " new nodes, exceeding maximum of " + std::to_string(MAX_NODES) + " total nodes");
	}

	// Step 5: Group placements by edge_id
	std::unordered_map<int64_t, std::vector<const Placement *>> by_edge;
	for (const auto &[frag_id, p] : best_placement) {
		by_edge[p->edge_id].push_back(p);
	}

	// Step 6: Sort each edge's placements by distal_length descending
	for (auto &[edge_id, edge_placements] : by_edge) {
		std::sort(edge_placements.begin(), edge_placements.end(),
		          [](const Placement *a, const Placement *b) { return a->distal_length > b->distal_length; });
	}

	// Step 7: Reserve space for new nodes (2 per placement: internal + fragment)
	nodes_.reserve(nodes_.size() + best_placement.size() * 2);

	// Step 8: Process each edge
	for (const auto &[edge_id, edge_placements] : by_edge) {
		uint32_t edge_node = edge_index.at(edge_id);
		double original_length = nodes_[edge_node].branch_length;

		// Get original parent
		uint32_t original_parent = nodes_[edge_node].parent;
		// Note: If original_parent == NO_PARENT, we're inserting on the root's edge.
		// In this case, current_parent starts as NO_PARENT, and the first new internal
		// node will become the new root (handled at line ~760 below).

		// Remove edge_node from its parent
		if (original_parent != NO_PARENT) {
			remove_child(original_parent, edge_node);
			nodes_[edge_node].parent = NO_PARENT;
		}

		// Insert placements from highest distal_length to lowest
		// This creates a chain: original_parent -> new_internal_1 -> new_internal_2 -> ... -> edge_node
		double remaining_length = original_length;
		uint32_t current_parent = original_parent;

		for (size_t i = 0; i < edge_placements.size(); ++i) {
			const Placement *p = edge_placements[i];

			// Calculate the branch length from current_parent to new internal node
			// This is: remaining_length - distal_length
			double internal_branch_length = remaining_length - p->distal_length;
			if (std::isnan(remaining_length)) {
				internal_branch_length = std::numeric_limits<double>::quiet_NaN();
			}

			// Create new internal node
			uint32_t new_internal = add_node("", internal_branch_length, std::nullopt);

			// Create fragment node
			uint32_t fragment_node = add_node(p->fragment_id, p->pendant_length, std::nullopt);

			// Connect internal node to current parent
			if (current_parent != NO_PARENT) {
				set_parent(new_internal, current_parent);
			} else {
				// This internal node becomes the new root
				root_ = new_internal;
			}

			// Connect fragment to internal
			set_parent(fragment_node, new_internal);

			// Update for next iteration
			current_parent = new_internal;
			remaining_length = p->distal_length;
		}

		// Finally, connect edge_node to the last internal node
		nodes_[edge_node].branch_length = remaining_length;
		set_parent(edge_node, current_parent);
	}
}

// ============================================================================
// Modification
// ============================================================================

NewickTree NewickTree::shear(const std::unordered_set<std::string> &keep_names, bool collapse,
                             bool ignore_missing) const {
	const uint32_t n = static_cast<uint32_t>(nodes_.size());

	// Map tip name -> tip node indices (only leaves; internal labels do not
	// count as tips). Duplicate tip labels map to multiple indices.
	std::unordered_map<std::string_view, std::vector<uint32_t>> tip_name_to_indices;
	for (uint32_t i = 0; i < n; ++i) {
		if (nodes_[i].children.empty()) {
			tip_name_to_indices[nodes_[i].name].push_back(i);
		}
	}

	// Upward closure: mark every kept tip and all of its ancestors. Early-stop
	// when we reach an already-marked node keeps this O(num_nodes) overall.
	std::vector<bool> marked(n, false);
	std::vector<std::string> missing;
	size_t kept_tip_count = 0;
	for (const auto &want : keep_names) {
		auto it = tip_name_to_indices.find(want);
		if (it == tip_name_to_indices.end()) {
			missing.push_back(want);
			continue;
		}
		for (uint32_t tip : it->second) {
			++kept_tip_count;
			uint32_t cur = tip;
			while (cur != NO_PARENT && !marked[cur]) {
				marked[cur] = true;
				cur = nodes_[cur].parent;
			}
		}
	}

	if (!missing.empty() && !ignore_missing) {
		std::sort(missing.begin(), missing.end());
		std::string msg =
		    "shear: " + std::to_string(missing.size()) + " requested tip name(s) not found as tips in the tree: ";
		constexpr size_t cap = 10;
		for (size_t i = 0; i < missing.size() && i < cap; ++i) {
			if (i > 0) {
				msg += ", ";
			}
			msg += "'" + missing[i] + "'";
		}
		if (missing.size() > cap) {
			msg += ", ... (" + std::to_string(missing.size() - cap) + " more)";
		}
		throw std::runtime_error(msg);
	}

	if (kept_tip_count == 0) {
		throw std::runtime_error("shear: no tips from the requested set matched a tip in the tree");
	}

	// Build the result tree in place. The kept nodes already have dense indices
	// (0..n-1) that this function controls, so there is no need to round-trip
	// through build()'s arbitrary-id remapping (hash maps) or a second copy of
	// every node via a NodeInput vector. Two passes are required because a
	// node's parent has a larger original index and is therefore created later
	// in the ascending scan; the parent must exist before it can be wired.
	//
	// Ordering matches build(): nodes are created in ascending original index
	// (so new node k = the k-th kept node by original index), and children are
	// appended in child-creation order. Freshly add_node'd nodes have no parent
	// yet, so set_parent never takes its remove-old-parent branch and children
	// append cleanly -- giving byte-identical output to the old build() path.
	NewickTree result;
	std::vector<uint32_t> remap(n, NO_PARENT); // original index -> new index
	std::vector<uint32_t> parents_old;         // per new node: parent's original index (NO_PARENT for root)

	auto emit = [&](uint32_t orig, double branch_length, uint32_t parent_orig) {
		// add_node copies the name once (shear is const, so the source cannot be
		// moved from); the old path copied it twice.
		remap[orig] = result.add_node(nodes_[orig].name, branch_length, nodes_[orig].edge_id);
		parents_old.push_back(parent_orig);
	};

	if (!collapse) {
		// Preserve every marked node with its original edge and parent (the root
		// keeps its original branch length; parent == NO_PARENT marks it).
		for (uint32_t i = 0; i < n; ++i) {
			if (marked[i]) {
				emit(i, nodes_[i].branch_length, nodes_[i].parent);
			}
		}
	} else {
		// collapse=true: a marked node is "retained" iff it is a tip or has >= 2
		// marked children; single-child internal nodes are dropped and their
		// edges merged. Count marked children (every marked non-root node's
		// parent is marked, by construction of the upward closure).
		std::vector<uint32_t> marked_child_count(n, 0);
		for (uint32_t i = 0; i < n; ++i) {
			if (marked[i] && nodes_[i].parent != NO_PARENT) {
				++marked_child_count[nodes_[i].parent];
			}
		}
		auto retained = [&](uint32_t node) -> bool {
			return marked[node] && (nodes_[node].children.empty() || marked_child_count[node] >= 2);
		};

		for (uint32_t i = 0; i < n; ++i) {
			if (!retained(i)) {
				continue;
			}

			// Walk up through the collapsed (single-child) intermediates to the
			// nearest retained ancestor, summing branch lengths onto the
			// surviving edge. NaN (unspecified) contributes 0, but an all-NaN
			// chain stays NaN so topology-only trees are preserved.
			double acc = 0.0;
			bool any_finite = false;
			auto add_bl = [&](double v) {
				if (!std::isnan(v)) {
					acc += v;
					any_finite = true;
				}
			};
			add_bl(nodes_[i].branch_length); // edge i -> parent(i)
			uint32_t cur = nodes_[i].parent;
			while (cur != NO_PARENT && !retained(cur)) {
				add_bl(nodes_[cur].branch_length);
				cur = nodes_[cur].parent;
			}

			// cur == NO_PARENT -> i is the LCA of the kept tips, i.e. the new
			// root; it has no incoming edge, so its branch length is dropped.
			double branch_length = (cur == NO_PARENT || !any_finite) ? std::numeric_limits<double>::quiet_NaN() : acc;
			emit(i, branch_length, cur);
		}
	}

	// Wire pass (shared): every parent target was created above, so set_parent
	// resolves cleanly and never hits its remove-old-parent branch.
	for (uint32_t j = 0; j < parents_old.size(); ++j) {
		uint32_t parent_orig = parents_old[j];
		if (parent_orig == NO_PARENT) {
			result.root_ = j;
		} else {
			result.set_parent(j, remap[parent_orig]);
		}
	}

	// Provably unreachable (there is always exactly one no-parent node: the
	// original root when collapse=false, the kept-tips' LCA when collapse=true),
	// but guard loudly rather than leave root_ as an out-of-bounds sentinel.
	if (result.root_ == NO_PARENT) {
		throw std::runtime_error("shear: internal error, no root produced");
	}

	return result;
}

NewickTree NewickTree::resolve_multifurcations() const {
	const uint32_t n = static_cast<uint32_t>(nodes_.size());

	NewickTree result;

	// Pass 1: copy every original node in ascending index order, preserving
	// name / branch length / edge id. Because nodes are added in order, original
	// index i maps to result index i; new connector nodes are appended after, so
	// original edges and node identities are untouched.
	for (uint32_t i = 0; i < n; ++i) {
		result.add_node(nodes_[i].name, nodes_[i].branch_length, nodes_[i].edge_id);
	}
	// build()/parse() always yield >= 1 node with a valid root, so root_ is a real
	// index here; an empty source would copy NO_PARENT through, matching an empty
	// result, so no separate guard is needed (cf. shear()'s explicit check).
	result.root_ = root_;

	// Pass 2: wire children. A node with <= 2 children is attached directly in
	// its original order (bifurcations and single-child unifurcations pass
	// through unchanged). A node with m >= 3 children is resolved into a
	// deterministic left-comb: c0 stays on the node and each further child hangs
	// off a fresh zero-length connector, so c0,c1,...,c(m-1) become successive
	// "left" children and the final connector holds the last two.
	for (uint32_t i = 0; i < n; ++i) {
		const auto &ch = nodes_[i].children;
		const size_t m = ch.size();
		if (m <= 2) {
			for (uint32_t c : ch) {
				result.set_parent(c, i);
			}
			continue;
		}
		uint32_t current = i;
		for (size_t k = 0; k < m; ++k) {
			result.set_parent(ch[k], current);
			if (k < m - 2) {
				uint32_t connector = result.add_node("", 0.0, std::nullopt);
				result.set_parent(connector, current);
				current = connector;
			}
		}
	}

	return result;
}

namespace {

// Trait-independent scaffold for Felsenstein independent contrasts: the validated
// tree plus per-node variance-extended branch lengths and per-internal-node
// contrast variances. Built once, reused across every trait over the same tree.
struct ContrastScaffold {
	std::vector<uint32_t> postorder;                            // children-before-parents
	std::vector<double> extended_length;                        // v' per node (raw branch length at tips)
	std::vector<double> contrast_variance;                      // v_i + v_j per internal node (unused at tips)
	std::unordered_map<std::string_view, uint32_t> tip_by_name; // for trait-completeness checks
};

// Human-readable node identifier for error messages. Prefers the node's name;
// for unnamed nodes, reports the caller-facing id from `node_ids` when available
// (so messages reference the source table's node_index, not the internal dense
// index), else the dense index.
std::string node_label(const NewickTree &t, uint32_t i, const std::vector<int64_t> *node_ids) {
	if (!t.name(i).empty()) {
		return "node '" + t.name(i) + "'";
	}
	if (node_ids && i < node_ids->size()) {
		return "node " + std::to_string((*node_ids)[i]);
	}
	return "node " + std::to_string(i);
}

std::string join_names(std::vector<std::string> &v) {
	std::sort(v.begin(), v.end());
	std::string s;
	constexpr size_t cap = 10;
	for (size_t i = 0; i < v.size() && i < cap; ++i) {
		if (i) {
			s += ", ";
		}
		s += "'" + v[i] + "'";
	}
	if (v.size() > cap) {
		s += ", ... (" + std::to_string(v.size() - cap) + " more)";
	}
	return s;
}

// Validate structure + branch lengths (all trait-independent) and precompute the
// variance-extended branch lengths and per-node contrast variances in one
// post-order pass. Throws on any structural / branch-length violation.
ContrastScaffold build_contrast_scaffold(const NewickTree &t, const std::vector<int64_t> *node_ids) {
	const uint32_t n = static_cast<uint32_t>(t.num_nodes());
	ContrastScaffold s;
	s.extended_length.assign(n, 0.0);
	s.contrast_variance.assign(n, 0.0);

	// Bifurcation + unique tip names.
	for (uint32_t i = 0; i < n; ++i) {
		const size_t nc = t.children(i).size();
		if (nc == 0) {
			bool inserted = s.tip_by_name.emplace(t.name(i), i).second;
			if (!inserted) {
				throw std::runtime_error("independent_contrasts: duplicate tip name '" + t.name(i) +
				                         "' (tip names must be unique to key traits)");
			}
		} else if (nc > 2) {
			throw std::runtime_error(
			    "independent_contrasts: " + node_label(t, i, node_ids) + " has " + std::to_string(nc) +
			    " children; requires a strictly bifurcating tree — use tree_resolve_multifurcations to resolve "
			    "polytomies");
		} else if (nc == 1) {
			throw std::runtime_error(
			    "independent_contrasts: " + node_label(t, i, node_ids) +
			    " has 1 child; requires a strictly bifurcating tree — use shear_tree (collapse) to remove "
			    "unifurcations");
		}
	}

	// Finite, non-negative branch length on every non-root edge. (The root's own
	// branch length is never used.) A non-finite length (NaN from an unspecified
	// edge, or ±Inf from a literal 'inf'/'infinity') would poison sqrt()/the
	// weighting and silently emit 0 / NaN contrasts, so reject both; a negative
	// length is biologically meaningless even if a sibling masks it in the sum.
	// Zero is allowed here — only a zero contrast *variance* is fatal, checked
	// below, so resolver-inserted zero-length internal edges pass.
	for (uint32_t i = 0; i < n; ++i) {
		if (i == t.root()) {
			continue;
		}
		double bl = t.branch_length(i);
		if (!std::isfinite(bl)) {
			throw std::runtime_error("independent_contrasts: " + node_label(t, i, node_ids) +
			                         " has a non-finite (NaN or infinite) branch length; PIC requires a finite branch "
			                         "length on every non-root edge");
		}
		if (bl < 0.0) {
			throw std::runtime_error("independent_contrasts: " + node_label(t, i, node_ids) +
			                         " has a negative branch length (" + std::to_string(bl) + ")");
		}
	}

	// Variance-extended branch lengths + per-node contrast variances in one
	// post-order pass (none of this depends on the trait values).
	s.postorder = t.postorder();
	for (uint32_t node : s.postorder) {
		if (t.is_tip(node)) {
			s.extended_length[node] = t.branch_length(node);
			continue;
		}
		uint32_t i = t.children(node)[0];
		uint32_t j = t.children(node)[1];
		double vi = s.extended_length[i];
		double vj = s.extended_length[j];
		double denom = vi + vj;
		if (!(denom > 0.0)) {
			throw std::runtime_error("independent_contrasts: " + node_label(t, node, node_ids) +
			                         " has two children with zero total branch length (contrast variance is zero)");
		}
		s.contrast_variance[node] = denom;
		// The root has no incoming edge, so only the variance-adjustment term
		// applies; every other node adds its own branch length.
		double own_bl = (node == t.root()) ? 0.0 : t.branch_length(node);
		s.extended_length[node] = own_bl + (vi * vj) / denom;
	}

	return s;
}

// Per-trait pass: validate completeness against the scaffold's tip set, then a
// single post-order computing ancestral values + standardized contrasts using the
// precomputed variances (Felsenstein 1985).
// `x_workspace` is a caller-owned scratch buffer sized to num_nodes; every entry
// is overwritten during the post-order (each node is written before any ancestor
// reads it), so it needs no re-initialization and can be reused across traits.
std::vector<IndependentContrast> contrasts_for_trait(const NewickTree &t, const ContrastScaffold &s,
                                                     const std::unordered_map<std::string, double> &trait_values,
                                                     std::vector<double> &x_workspace) {
	// Trait completeness: exactly the tip set (no missing tips, no extras).
	std::vector<std::string> missing;
	for (const auto &entry : s.tip_by_name) {
		if (trait_values.find(std::string(entry.first)) == trait_values.end()) {
			missing.emplace_back(entry.first);
		}
	}
	if (!missing.empty()) {
		throw std::runtime_error("independent_contrasts: no trait value for tip(s): " + join_names(missing));
	}
	if (trait_values.size() != s.tip_by_name.size()) {
		std::vector<std::string> extra;
		for (const auto &[name, val] : trait_values) {
			if (s.tip_by_name.find(name) == s.tip_by_name.end()) {
				extra.push_back(name);
			}
		}
		throw std::runtime_error("independent_contrasts: trait value(s) for name(s) that are not a tip in the tree: " +
		                         join_names(extra));
	}
	// Finiteness: a NaN/infinite trait value would silently propagate into every contrast
	// and ancestral estimate on the path to the root (and any downstream aggregate), so
	// reject it — fail loud, matching ancestral_states_bm.
	for (const auto &[nm, val] : trait_values) {
		if (!std::isfinite(val)) {
			throw std::runtime_error("independent_contrasts: non-finite (NaN or infinite) trait value for tip '" + nm +
			                         "'");
		}
	}

	std::vector<double> &X = x_workspace; // reconstructed ancestral value (tip trait value at leaves)
	std::vector<IndependentContrast> result;
	result.reserve(s.tip_by_name.empty() ? 0 : s.tip_by_name.size() - 1);

	for (uint32_t node : s.postorder) {
		if (t.is_tip(node)) {
			X[node] = trait_values.at(t.name(node)); // present (validated above)
			continue;
		}
		uint32_t i = t.children(node)[0];
		uint32_t j = t.children(node)[1];
		double vi = s.extended_length[i];
		double vj = s.extended_length[j];
		double denom = s.contrast_variance[node]; // precomputed, guaranteed > 0
		double contrast = (X[i] - X[j]) / std::sqrt(denom);
		double anc = (X[i] * vj + X[j] * vi) / denom;
		X[node] = anc;
		result.push_back(IndependentContrast {node, contrast, anc, denom});
	}

	return result;
}

} // namespace

std::vector<IndependentContrast>
NewickTree::independent_contrasts(const std::unordered_map<std::string, double> &trait_values,
                                  const std::vector<int64_t> *node_ids) const {
	ContrastScaffold scaffold = build_contrast_scaffold(*this, node_ids);
	std::vector<double> x_workspace(num_nodes());
	return contrasts_for_trait(*this, scaffold, trait_values, x_workspace);
}

std::vector<std::vector<IndependentContrast>>
NewickTree::independent_contrasts(const std::vector<std::unordered_map<std::string, double>> &trait_values_list,
                                  const std::vector<int64_t> *node_ids) const {
	// Trait-independent validation + variances computed once, reused per trait; the
	// tree-sized ancestral-value scratch buffer is allocated once and reused too.
	ContrastScaffold scaffold = build_contrast_scaffold(*this, node_ids);
	std::vector<double> x_workspace(num_nodes());
	std::vector<std::vector<IndependentContrast>> results;
	results.reserve(trait_values_list.size());
	for (const auto &trait_values : trait_values_list) {
		results.push_back(contrasts_for_trait(*this, scaffold, trait_values, x_workspace));
	}
	return results;
}

namespace {

// Trait-independent scaffold for Brownian-motion ancestral reconstruction: the
// validated tree plus per-node structural (rate-independent) variances V_down
// (subtree-conditional) and V_up (rest-of-tree). Built once, reused across every
// trait over the same tree (mirrors ContrastScaffold for PIC). Requires every
// internal node to have >= 2 children.
struct BmScaffold {
	std::vector<uint32_t> postorder;                            // children before parents
	std::vector<uint32_t> preorder;                             // parents before children
	std::vector<double> V_down;                                 // subtree-conditional variance (0 at tips)
	std::vector<double> V_up;                                   // rest-of-tree variance (inf at root)
	std::unordered_map<std::string_view, uint32_t> tip_by_name; // for trait-completeness checks
	uint32_t n_tips = 0;
	uint32_t max_arity = 0; // largest child count (for scratch sizing)
};

// Validate structure + branch lengths (all trait-independent) and precompute the
// structural variances V_down (post-order) and V_up (pre-order). The up-pass
// leave-one-out uses prefix/suffix precision sums (never a total-minus-self
// subtraction), so it stays numerically robust even when one child's precision
// dwarfs the others (e.g. a zero-length resolver edge). Throws on any violation.
BmScaffold build_bm_scaffold(const NewickTree &t, const std::vector<int64_t> *node_ids) {
	const uint32_t n = static_cast<uint32_t>(t.num_nodes());
	BmScaffold s;
	s.V_down.assign(n, 0.0);
	s.V_up.assign(n, std::numeric_limits<double>::infinity());

	// Unique tip names; reject unifurcations. Multifurcations (>= 2 children) are
	// supported, but a single-child internal node carries no information and leaves
	// the leave-one-out cavity empty, so it is rejected with a shear_tree hint
	// (matching PIC's unifurcation handling).
	for (uint32_t i = 0; i < n; ++i) {
		const size_t nc = t.children(i).size();
		if (nc == 0) {
			s.n_tips++;
			if (!s.tip_by_name.emplace(t.name(i), i).second) {
				throw std::runtime_error("ancestral_states_bm: duplicate tip name '" + t.name(i) +
				                         "' (tip names must be unique to key traits)");
			}
		} else if (nc == 1) {
			throw std::runtime_error("ancestral_states_bm: " + node_label(t, i, node_ids) +
			                         " has 1 child; continuous ASR does not reconstruct unifurcations — use "
			                         "shear_tree (collapse) to remove them");
		} else if (nc > s.max_arity) {
			s.max_arity = static_cast<uint32_t>(nc);
		}
	}

	// Branch-length policy: finite, non-negative on every non-root edge; strictly
	// positive on every tip edge. A zero-length tip edge would pin the node with
	// infinite precision and break the leave-one-out, so it is rejected here (a
	// departure from PIC). Zero on internal edges is allowed (e.g. connectors from
	// tree_resolve_multifurcations); with positive tip edges every message variance
	// stays > 0.
	for (uint32_t i = 0; i < n; ++i) {
		if (i == t.root()) {
			continue;
		}
		double bl = t.branch_length(i);
		if (!std::isfinite(bl)) {
			throw std::runtime_error("ancestral_states_bm: " + node_label(t, i, node_ids) +
			                         " has a non-finite (NaN or infinite) branch length; continuous ASR requires a "
			                         "finite branch length on every non-root edge");
		}
		if (bl < 0.0) {
			throw std::runtime_error("ancestral_states_bm: " + node_label(t, i, node_ids) +
			                         " has a negative branch length (" + std::to_string(bl) + ")");
		}
		if (t.is_tip(i) && !(bl > 0.0)) {
			throw std::runtime_error("ancestral_states_bm: tip '" + t.name(i) +
			                         "' has a zero-length branch; continuous ASR requires a strictly positive branch "
			                         "length on every tip edge");
		}
	}

	s.postorder = t.postorder();
	s.preorder = t.preorder();

	// Down-pass: V_down[u] = 1 / sum_c 1/(V_down[c] + bl_c). Trait-independent.
	for (uint32_t u : s.postorder) {
		if (t.is_tip(u)) {
			continue; // V_down = 0
		}
		double prec = 0.0;
		for (uint32_t c : t.children(u)) {
			prec += 1.0 / (s.V_down[c] + t.branch_length(c));
		}
		s.V_down[u] = 1.0 / prec;
	}

	// Up-pass: V_up[c] = 1/cavity_precision + bl_c, where the cavity precision at
	// parent p for child c is p's up-precision plus every OTHER child's precision,
	// summed via prefix/suffix (no total-minus-self subtraction). Every internal node
	// has >= 2 children, so the cavity always contains at least one positive term.
	std::vector<double> pref(s.max_arity + 1, 0.0); // reused scratch: prefix precision sums
	for (uint32_t p : s.preorder) {
		if (t.is_tip(p)) {
			continue;
		}
		const auto &ch = t.children(p);
		const size_t m = ch.size();
		double up_prec = std::isfinite(s.V_up[p]) ? 1.0 / s.V_up[p] : 0.0;
		pref[0] = up_prec;
		for (size_t i = 0; i < m; ++i) {
			pref[i + 1] = pref[i] + 1.0 / (s.V_down[ch[i]] + t.branch_length(ch[i]));
		}
		double suf = 0.0;
		for (size_t i = m; i-- > 0;) {
			double cav_prec = pref[i] + suf; // up_prec + every child except i
			s.V_up[ch[i]] = 1.0 / cav_prec + t.branch_length(ch[i]);
			suf += 1.0 / (s.V_down[ch[i]] + t.branch_length(ch[i]));
		}
	}

	return s;
}

// Per-trait pass: validate completeness + finiteness against the scaffold's tip
// set, then compute means (down + up via the same prefix/suffix leave-one-out),
// the REML rate sigma^2, and the marginal estimate + variance + 95% CI at each
// internal node.
//
// REML BM rate: sum of squared standardized contrasts over internal nodes /
// (n_tips - 1). Justification for REML over full ML (divisor n_tips): the
// n_tips - 1 standardized contrasts are iid N(0, sigma^2), so Sum(contrast^2)/
// (n_tips - 1) is the UNBIASED estimator, whereas full ML (as in ape::ace
// method="ML") divides by n_tips and is biased low by a factor (n_tips - 1)/n_tips.
// REML is also PIC-consistent: PIC (which we already ship) is built on exactly
// these contrasts. A degree-m node yields m - 1 contrasts; the weighted residual
// sum Sum_c (m_down[c] - m_down[u])^2 / vp_c reduces to the squared standardized
// PIC contrast at a bifurcation, and summed over internal nodes has n_tips - 1 df
// on any resolution. `variance` and the CI are this REML rate times the structural
// (rate-independent) variance.
std::vector<AncestralStateBM> bm_states_for_trait(const NewickTree &t, const BmScaffold &s,
                                                  const std::unordered_map<std::string, double> &trait_values) {
	// Completeness: exactly the tip set (no missing tips, no extras).
	std::vector<std::string> missing;
	for (const auto &entry : s.tip_by_name) {
		if (trait_values.find(std::string(entry.first)) == trait_values.end()) {
			missing.emplace_back(entry.first);
		}
	}
	if (!missing.empty()) {
		throw std::runtime_error("ancestral_states_bm: no trait value for tip(s): " + join_names(missing));
	}
	if (trait_values.size() != s.tip_by_name.size()) {
		std::vector<std::string> extra;
		for (const auto &[nm, val] : trait_values) {
			if (s.tip_by_name.find(nm) == s.tip_by_name.end()) {
				extra.push_back(nm);
			}
		}
		throw std::runtime_error("ancestral_states_bm: trait value(s) for name(s) that are not a tip in the tree: " +
		                         join_names(extra));
	}
	// Finiteness: a NaN/infinite trait value would silently poison every estimate and
	// the shared REML rate (hence the whole output table), so reject it — fail loud.
	for (const auto &[nm, val] : trait_values) {
		if (!std::isfinite(val)) {
			throw std::runtime_error("ancestral_states_bm: non-finite (NaN or infinite) trait value for tip '" + nm +
			                         "'");
		}
	}

	const uint32_t n = static_cast<uint32_t>(t.num_nodes());
	std::vector<double> m_down(n, 0.0);
	std::vector<double> m_up(n, 0.0);

	// Down-pass means: m_down[u] = V_down[u] * sum_c m_down[c]/(V_down[c] + bl_c).
	for (uint32_t u : s.postorder) {
		if (t.is_tip(u)) {
			m_down[u] = trait_values.at(t.name(u));
			continue;
		}
		double wmean = 0.0;
		for (uint32_t c : t.children(u)) {
			wmean += m_down[c] / (s.V_down[c] + t.branch_length(c));
		}
		m_down[u] = s.V_down[u] * wmean;
	}

	// Up-pass means: cavity weighted-mean via prefix/suffix (same leave-one-out as
	// the scaffold's precision pass; the cavity precision is recomputed here so the
	// division uses matching operands).
	std::vector<double> pref_p(s.max_arity + 1, 0.0);
	std::vector<double> pref_wm(s.max_arity + 1, 0.0);
	for (uint32_t p : s.preorder) {
		if (t.is_tip(p)) {
			continue;
		}
		const auto &ch = t.children(p);
		const size_t m = ch.size();
		double up_prec = std::isfinite(s.V_up[p]) ? 1.0 / s.V_up[p] : 0.0;
		pref_p[0] = up_prec;
		pref_wm[0] = up_prec * m_up[p]; // 0 at the root
		for (size_t i = 0; i < m; ++i) {
			double pc = 1.0 / (s.V_down[ch[i]] + t.branch_length(ch[i]));
			pref_p[i + 1] = pref_p[i] + pc;
			pref_wm[i + 1] = pref_wm[i] + pc * m_down[ch[i]];
		}
		double suf_p = 0.0, suf_wm = 0.0;
		for (size_t i = m; i-- > 0;) {
			double pc = 1.0 / (s.V_down[ch[i]] + t.branch_length(ch[i]));
			m_up[ch[i]] = (pref_wm[i] + suf_wm) / (pref_p[i] + suf_p);
			suf_p += pc;
			suf_wm += pc * m_down[ch[i]];
		}
	}

	double sum_wrss = 0.0;
	for (uint32_t u : s.postorder) {
		if (t.is_tip(u)) {
			continue;
		}
		for (uint32_t c : t.children(u)) {
			double d = m_down[c] - m_down[u];
			sum_wrss += d * d / (s.V_down[c] + t.branch_length(c));
		}
	}
	double sigma2 = sum_wrss / static_cast<double>(s.n_tips - 1);

	// Marginal at each internal node: combine the subtree (down) and rest-of-tree
	// (up) messages by precision weighting; scale the structural variance by sigma^2.
	constexpr double z975 = 1.959963984540054; // qnorm(0.975)
	std::vector<AncestralStateBM> result;
	result.reserve(n - s.n_tips);
	for (uint32_t u = 0; u < n; ++u) {
		if (t.is_tip(u)) {
			continue;
		}
		double up_prec = std::isfinite(s.V_up[u]) ? 1.0 / s.V_up[u] : 0.0;
		double full_prec = 1.0 / s.V_down[u] + up_prec;
		double estimate = (m_down[u] / s.V_down[u] + up_prec * m_up[u]) / full_prec;
		double variance = sigma2 * (1.0 / full_prec);
		double half = z975 * std::sqrt(variance);
		result.push_back(AncestralStateBM {u, estimate, variance, estimate - half, estimate + half});
	}
	return result;
}

} // namespace

std::vector<AncestralStateBM>
NewickTree::ancestral_states_bm(const std::unordered_map<std::string, double> &trait_values,
                                const std::vector<int64_t> *node_ids) const {
	BmScaffold scaffold = build_bm_scaffold(*this, node_ids);
	if (scaffold.n_tips <= 1) {
		// A degenerate tree (<= 1 tip) has no internal nodes to reconstruct, so it returns
		// empty here -- before per-trait completeness validation runs. (independent_contrasts
		// always validates; for these methods a degenerate tree has nothing to validate
		// against, so the early-out is intentional and harmless.)
		return {};
	}
	return bm_states_for_trait(*this, scaffold, trait_values);
}

std::vector<std::vector<AncestralStateBM>>
NewickTree::ancestral_states_bm(const std::vector<std::unordered_map<std::string, double>> &trait_values_list,
                                const std::vector<int64_t> *node_ids) const {
	// Trait-independent validation + structural variances computed once, reused per
	// trait (mirrors the PIC batch overload).
	BmScaffold scaffold = build_bm_scaffold(*this, node_ids);
	std::vector<std::vector<AncestralStateBM>> results;
	results.reserve(trait_values_list.size());
	for (const auto &trait_values : trait_values_list) {
		if (scaffold.n_tips <= 1) {
			results.emplace_back();
		} else {
			results.push_back(bm_states_for_trait(*this, scaffold, trait_values));
		}
	}
	return results;
}

namespace {

// Trait-independent scaffold for Sankoff parsimony: the validated tree plus the
// post-/pre-order traversals and the tip-name index (mirrors ContrastScaffold /
// BmScaffold). Parsimony ignores branch lengths and supports any arity (including
// unifurcations), so the only structural requirement is unique tip names. Built
// once, reused across every trait over the same tree.
struct ParsimonyScaffold {
	std::vector<uint32_t> postorder; // children before parents
	std::vector<uint32_t> preorder;  // parents before children
	std::unordered_map<std::string_view, uint32_t> tip_by_name;
	uint32_t n_tips = 0;
};

ParsimonyScaffold build_parsimony_scaffold(const NewickTree &t) {
	ParsimonyScaffold s;
	const uint32_t n = static_cast<uint32_t>(t.num_nodes());
	for (uint32_t i = 0; i < n; ++i) {
		if (t.is_tip(i)) {
			s.n_tips++;
			if (!s.tip_by_name.emplace(t.name(i), i).second) {
				throw std::runtime_error("ancestral_parsimony: duplicate tip name '" + t.name(i) +
				                         "' (tip names must be unique to key traits)");
			}
		}
	}
	s.postorder = t.postorder();
	s.preorder = t.preorder();
	return s;
}

// Per-trait Sankoff (1975) DP; the unit cost matrix (0 on the diagonal, 1 elsewhere)
// gives Fitch (1971) parsimony. `tip_state` maps tip name -> state index (0..k-1);
// `cost` is the k*k row-major substitution cost matrix, cost[from*k + to] being the
// cost of an edge whose parent end is `from` and child end is `to`.
//
// Down-pass (post-order): g_u[s] = min cost of u's subtree given u=s
//   = sum over children c of min_{s'} (cost[s*k + s'] + g_c[s']); tips are 0 at the
//   observed state and +inf elsewhere. score = min_s g_root[s].
// Up-pass (pre-order): h_u[s] = min cost of the rest of the tree given u=s; h_root=0.
//   For a child c of p with msg_c[s'] = min_t (cost[s'*k + t] + g_c[t]) (the same
//   term g_p sums over its children), the sibling contribution is g_p[s'] - msg_c[s']
//   and above_c[s'] = h_p[s'] + that; then h_c[t] = min_{s'} (cost[s'*k + t] +
//   above_c[s']). Subtraction is safe here (unlike BM's precision-form up-pass): every
//   internal g is finite (each child sends a finite message), so g_p - msg_c never
//   touches an infinite/zero-precision regime.
// F_u[s] = g_u[s] + h_u[s] is the min total tree cost with u fixed to s; a state is
// in the MPR set iff F_u[s] == score within a hybrid absolute/relative tolerance
// (max(1e-9, 1e-6*|score|)): the absolute floor decides exactly for integer/unit
// costs (differences are whole numbers), while the relative term absorbs rounding at
// large real-valued cost scales.
std::vector<AncestralStateParsimony> parsimony_for_trait(const NewickTree &t, const ParsimonyScaffold &s,
                                                         const std::unordered_map<std::string, uint32_t> &tip_state,
                                                         uint32_t k, const std::vector<double> &cost) {
	// Cost-matrix validity (carried per trait because the unit-cost default builds one
	// per trait's alphabet): k >= 1, k*k entries, all finite and non-negative.
	if (k == 0) {
		throw std::runtime_error("ancestral_parsimony: state alphabet is empty");
	}
	if (cost.size() != static_cast<size_t>(k) * k) {
		throw std::runtime_error("ancestral_parsimony: cost matrix must have k*k entries");
	}
	for (double c : cost) {
		if (!std::isfinite(c) || c < 0.0) {
			throw std::runtime_error("ancestral_parsimony: cost matrix entries must be finite and non-negative");
		}
	}

	// Completeness: exactly the tip set (no missing tips, no extras); every state < k.
	std::vector<std::string> missing;
	for (const auto &entry : s.tip_by_name) {
		if (tip_state.find(std::string(entry.first)) == tip_state.end()) {
			missing.emplace_back(entry.first);
		}
	}
	if (!missing.empty()) {
		throw std::runtime_error("ancestral_parsimony: no state for tip(s): " + join_names(missing));
	}
	if (tip_state.size() != s.tip_by_name.size()) {
		std::vector<std::string> extra;
		for (const auto &[nm, st] : tip_state) {
			if (s.tip_by_name.find(nm) == s.tip_by_name.end()) {
				extra.push_back(nm);
			}
		}
		throw std::runtime_error("ancestral_parsimony: state(s) for name(s) that are not a tip in the tree: " +
		                         join_names(extra));
	}
	for (const auto &[nm, st] : tip_state) {
		if (st >= k) {
			throw std::runtime_error("ancestral_parsimony: tip '" + nm + "' has state index " + std::to_string(st) +
			                         " out of range for k=" + std::to_string(k));
		}
	}

	const uint32_t n = static_cast<uint32_t>(t.num_nodes());
	const double INF = std::numeric_limits<double>::infinity();
	std::vector<double> g(static_cast<size_t>(n) * k, 0.0);

	// Down-pass.
	for (uint32_t u : s.postorder) {
		double *gu = &g[static_cast<size_t>(u) * k];
		if (t.is_tip(u)) {
			uint32_t obs = tip_state.at(t.name(u));
			for (uint32_t st = 0; st < k; ++st) {
				gu[st] = (st == obs) ? 0.0 : INF;
			}
			continue;
		}
		for (uint32_t st = 0; st < k; ++st) {
			gu[st] = 0.0;
		}
		for (uint32_t c : t.children(u)) {
			const double *gc = &g[static_cast<size_t>(c) * k];
			for (uint32_t st = 0; st < k; ++st) {
				const double *row = &cost[static_cast<size_t>(st) * k];
				double best = INF;
				for (uint32_t sp = 0; sp < k; ++sp) {
					best = std::min(best, row[sp] + gc[sp]);
				}
				gu[st] += best;
			}
		}
	}

	double score = INF;
	{
		const double *gr = &g[static_cast<size_t>(t.root()) * k];
		for (uint32_t st = 0; st < k; ++st) {
			score = std::min(score, gr[st]);
		}
	}

	// Up-pass. h's root row stays 0.
	std::vector<double> h(static_cast<size_t>(n) * k, 0.0);
	std::vector<double> above(k);
	for (uint32_t p : s.preorder) {
		if (t.is_tip(p)) {
			continue;
		}
		const double *gp = &g[static_cast<size_t>(p) * k];
		const double *hp = &h[static_cast<size_t>(p) * k];
		for (uint32_t c : t.children(p)) {
			const double *gc = &g[static_cast<size_t>(c) * k];
			for (uint32_t sp = 0; sp < k; ++sp) {
				const double *row = &cost[static_cast<size_t>(sp) * k];
				double msg = INF; // msg_c[sp] = min_t (cost[sp*k + t] + g_c[t])
				for (uint32_t tt = 0; tt < k; ++tt) {
					msg = std::min(msg, row[tt] + gc[tt]);
				}
				above[sp] = hp[sp] + (gp[sp] - msg); // h_p + sibling contribution
			}
			double *hc = &h[static_cast<size_t>(c) * k];
			for (uint32_t tt = 0; tt < k; ++tt) {
				double best = INF;
				for (uint32_t sp = 0; sp < k; ++sp) {
					best = std::min(best, cost[static_cast<size_t>(sp) * k + tt] + above[sp]);
				}
				hc[tt] = best;
			}
		}
	}

	// Marginal per (internal node, state).
	const double tol = std::max(1e-9, 1e-6 * std::abs(score));
	std::vector<AncestralStateParsimony> result;
	result.reserve(static_cast<size_t>(n - s.n_tips) * k);
	for (uint32_t u = 0; u < n; ++u) {
		if (t.is_tip(u)) {
			continue;
		}
		const double *gu = &g[static_cast<size_t>(u) * k];
		const double *hu = &h[static_cast<size_t>(u) * k];
		for (uint32_t st = 0; st < k; ++st) {
			double f = gu[st] + hu[st];
			result.push_back(AncestralStateParsimony {u, st, f <= score + tol, f});
		}
	}
	return result;
}

} // namespace

std::vector<AncestralStateParsimony>
NewickTree::ancestral_parsimony(const std::unordered_map<std::string, uint32_t> &tip_states, uint32_t k,
                                const std::vector<double> &cost, const std::vector<int64_t> *node_ids) const {
	// Parsimony's validation errors reference tip NAMES (or the cost matrix), never an
	// unnamed internal node, so node_ids is accepted for interface symmetry but unused.
	(void)node_ids;
	ParsimonyScaffold scaffold = build_parsimony_scaffold(*this);
	if (scaffold.n_tips <= 1) {
		// A degenerate tree (<= 1 tip) has no internal nodes to reconstruct, so it returns
		// empty here -- before per-trait completeness validation runs. (independent_contrasts
		// always validates; for these methods a degenerate tree has nothing to validate
		// against, so the early-out is intentional and harmless.)
		return {};
	}
	return parsimony_for_trait(*this, scaffold, tip_states, k, cost);
}

std::vector<std::vector<AncestralStateParsimony>>
NewickTree::ancestral_parsimony(const std::vector<std::unordered_map<std::string, uint32_t>> &tip_states_list,
                                const std::vector<uint32_t> &k_list, const std::vector<std::vector<double>> &cost_list,
                                const std::vector<int64_t> *node_ids) const {
	(void)node_ids;
	if (tip_states_list.size() != k_list.size() || tip_states_list.size() != cost_list.size()) {
		throw std::runtime_error("ancestral_parsimony: tip_states_list, k_list, and cost_list must be parallel");
	}
	// Trait-independent validation + traversals computed once, reused per trait.
	ParsimonyScaffold scaffold = build_parsimony_scaffold(*this);
	std::vector<std::vector<AncestralStateParsimony>> results;
	results.reserve(tip_states_list.size());
	for (size_t i = 0; i < tip_states_list.size(); ++i) {
		if (scaffold.n_tips <= 1) {
			results.emplace_back();
		} else {
			results.push_back(parsimony_for_trait(*this, scaffold, tip_states_list[i], k_list[i], cost_list[i]));
		}
	}
	return results;
}

namespace {

// Trait-independent scaffold for Mk maximum-likelihood: the validated tree plus the
// post-/pre-order traversals and the tip-name index (mirrors ContrastScaffold /
// BmScaffold / ParsimonyScaffold). ML supports any arity; the branch-length policy
// matches BM (finite non-negative non-root edges; strictly positive tip edges — a
// zero-length tip edge puts a hard zero in the leave-one-out; zero internal edges are
// allowed). Built once, reused across every trait over the same tree.
struct MLScaffold {
	std::vector<uint32_t> postorder;
	std::vector<uint32_t> preorder;
	std::unordered_map<std::string_view, uint32_t> tip_by_name;
	uint32_t n_tips = 0;
};

MLScaffold build_ml_scaffold(const NewickTree &t, const std::vector<int64_t> *node_ids) {
	const uint32_t n = static_cast<uint32_t>(t.num_nodes());
	MLScaffold s;
	for (uint32_t i = 0; i < n; ++i) {
		if (t.is_tip(i)) {
			s.n_tips++;
			if (!s.tip_by_name.emplace(t.name(i), i).second) {
				throw std::runtime_error("ancestral_ml: duplicate tip name '" + t.name(i) +
				                         "' (tip names must be unique to key traits)");
			}
		}
	}
	for (uint32_t i = 0; i < n; ++i) {
		if (i == t.root()) {
			continue;
		}
		double bl = t.branch_length(i);
		if (!std::isfinite(bl)) {
			throw std::runtime_error("ancestral_ml: " + node_label(t, i, node_ids) +
			                         " has a non-finite (NaN or infinite) branch length; ML requires a finite branch "
			                         "length on every non-root edge");
		}
		if (bl < 0.0) {
			throw std::runtime_error("ancestral_ml: " + node_label(t, i, node_ids) + " has a negative branch length (" +
			                         std::to_string(bl) + ")");
		}
		if (t.is_tip(i) && !(bl > 0.0)) {
			throw std::runtime_error("ancestral_ml: tip '" + t.name(i) +
			                         "' has a zero-length branch; ML requires a strictly positive branch length on "
			                         "every tip edge");
		}
	}
	s.postorder = t.postorder();
	s.preorder = t.preorder();
	return s;
}

// Felsenstein (1981) pruning down-pass at ER rate r, with per-node rescaling. Fills L
// (n*k, each internal row rescaled so its max entry is 1) and returns the log-
// likelihood under a uniform root prior 1/k (the retained scalers are added back). For
// ER, P(t)[i][j] = 1/k + (delta_ij - 1/k) exp(-k r t), so the message from child c to
// its parent in parent-state s is M_c(s) = E*L_c[s] + ((1-E)/k)*sum_s' L_c[s'] with
// E = exp(-k r t_c) -- O(k) per child, no k*k matrix needed.
double ml_down_pass(const NewickTree &t, const MLScaffold &s, const std::vector<uint32_t> &obs_state, uint32_t k,
                    double r, std::vector<double> &L) {
	double logscale = 0.0;
	for (uint32_t u : s.postorder) {
		double *Lu = &L[static_cast<size_t>(u) * k];
		if (t.is_tip(u)) {
			uint32_t o = obs_state[u];
			for (uint32_t st = 0; st < k; ++st) {
				Lu[st] = (st == o) ? 1.0 : 0.0;
			}
			continue;
		}
		for (uint32_t st = 0; st < k; ++st) {
			Lu[st] = 1.0;
		}
		for (uint32_t c : t.children(u)) {
			const double *Lc = &L[static_cast<size_t>(c) * k];
			double SL = 0.0;
			for (uint32_t st = 0; st < k; ++st) {
				SL += Lc[st];
			}
			double E = std::exp(-static_cast<double>(k) * r * t.branch_length(c));
			double f = (1.0 - E) / k * SL;
			for (uint32_t st = 0; st < k; ++st) {
				Lu[st] *= (E * Lc[st] + f);
			}
		}
		double m = 0.0;
		for (uint32_t st = 0; st < k; ++st) {
			m = std::max(m, Lu[st]);
		}
		if (m > 0.0) {
			for (uint32_t st = 0; st < k; ++st) {
				Lu[st] /= m;
			}
			logscale += std::log(m);
		}
	}
	const double *Lr = &L[static_cast<size_t>(t.root()) * k];
	double liksum = 0.0;
	for (uint32_t st = 0; st < k; ++st) {
		liksum += Lr[st];
	}
	return std::log(liksum / k) + logscale; // uniform prior pi = 1/k
}

// Marginal posterior up-pass (Yang, Kumar & Nei 1995) from the down-pass L. G (n*k) is
// scratch for the "outside" likelihoods. The sibling contribution is obtained without
// a product accumulation as L_p[s_p] / M_c(s_p) (the full-children product over the
// per-child message) -- O(1) per (child, state), no degree-dependent underflow. For any
// branch with exp(-k r t) < 1 the message M_c(s) = E*L_c[s] + ((1-E)/k)*sum L_c is
// strictly positive, so the division is well defined. It can fail ONLY on a
// pathologically tiny (but validly > 0) branch whose exp(-k r t) rounds to EXACTLY 1.0
// (k*r*t <~ 1e-16): such an edge acts like a zero-length edge, collapsing a message to a
// hard delta, and conflicting states can then drive M_c(s) (and L_p[s]) to exactly 0 ->
// 0/0. We guard that at the division site and throw (fail loud) rather than emit a NaN
// probability. Per-node rescaling of G is harmless (it cancels in the normalization).
// Posteriors are written for internal nodes only.
void ml_posteriors(const NewickTree &t, const MLScaffold &s, uint32_t k, double r, const std::vector<double> &L,
                   std::vector<double> &G, std::vector<AncestralStateML> &out) {
	const uint32_t root = t.root();
	{
		double *Gr = &G[static_cast<size_t>(root) * k];
		for (uint32_t st = 0; st < k; ++st) {
			Gr[st] = 1.0 / k; // pi
		}
	}
	std::vector<double> W(k);
	for (uint32_t p : s.preorder) {
		if (t.is_tip(p)) {
			continue;
		}
		const double *Lp = &L[static_cast<size_t>(p) * k];
		const double *Gp = &G[static_cast<size_t>(p) * k];
		for (uint32_t c : t.children(p)) {
			const double *Lc = &L[static_cast<size_t>(c) * k];
			double SLc = 0.0;
			for (uint32_t st = 0; st < k; ++st) {
				SLc += Lc[st];
			}
			double E = std::exp(-static_cast<double>(k) * r * t.branch_length(c));
			double f = (1.0 - E) / k * SLc;
			double SW = 0.0;
			for (uint32_t sp = 0; sp < k; ++sp) {
				double Mc = E * Lc[sp] + f; // message from c to p at parent-state sp
				if (!(Mc > 0.0)) {
					// Only reachable when exp(-k r t_c) rounded to exactly 1.0 (an
					// effectively-zero branch) and the subtree excludes this state.
					throw std::runtime_error(
					    "ancestral_ml: a subtree message underflowed to zero at the fitted/fixed rate; a branch "
					    "length is too small relative to the rate (exp(-k*rate*t) rounds to 1) to reconstruct "
					    "posteriors -- remove or lengthen effectively-zero-length edges");
				}
				W[sp] = Gp[sp] * Lp[sp] / Mc; // outside * (siblings' product = L_p / M_c)
				SW += W[sp];
			}
			double *Gc = &G[static_cast<size_t>(c) * k];
			double fc = (1.0 - E) / k * SW;
			double m = 0.0;
			for (uint32_t st = 0; st < k; ++st) {
				Gc[st] = E * W[st] + fc; // sum_sp W[sp] * P(t_c)[sp][st], ER closed form
				m = std::max(m, Gc[st]);
			}
			if (m > 0.0) {
				for (uint32_t st = 0; st < k; ++st) {
					Gc[st] /= m;
				}
			}
		}
	}
	const uint32_t n = static_cast<uint32_t>(t.num_nodes());
	for (uint32_t u = 0; u < n; ++u) {
		if (t.is_tip(u)) {
			continue;
		}
		const double *Lu = &L[static_cast<size_t>(u) * k];
		const double *Gu = &G[static_cast<size_t>(u) * k];
		double norm = 0.0;
		for (uint32_t st = 0; st < k; ++st) {
			norm += Lu[st] * Gu[st];
		}
		for (uint32_t st = 0; st < k; ++st) {
			out.push_back(AncestralStateML {u, st, Lu[st] * Gu[st] / norm});
		}
	}
}

// Golden-section search for the argmax of a unimodal f on [a, b] (Kiefer 1953). Used
// over x = log(rate); no derivatives, robust, cleanroom (not from Numerical Recipes).
template <typename F>
double golden_section_argmax(F &&f, double a, double b) {
	const double phi = 0.6180339887498949; // (sqrt(5) - 1) / 2
	double c = b - phi * (b - a), d = a + phi * (b - a);
	double fc = f(c), fd = f(d);
	for (int it = 0; it < 200 && (b - a) > 1e-8; ++it) {
		if (fc > fd) {
			b = d;
			d = c;
			fd = fc;
			c = b - phi * (b - a);
			fc = f(c);
		} else {
			a = c;
			c = d;
			fc = fd;
			d = a + phi * (b - a);
			fd = f(d);
		}
	}
	return (fc > fd) ? c : d;
}

// Nelder-Mead downhill simplex minimizer (Nelder & Mead, Computer Journal 7:308-313,
// 1965). Derivative-free; cleanroom (standard reflection/expansion/contraction/shrink
// coefficients, NOT from Numerical Recipes). Minimizes f over R^p from x0 with initial
// simplex edge `step`; returns the best vertex found. Used to fit the SYM rate matrix.
template <typename F>
std::vector<double> nelder_mead_minimize(F &&f, const std::vector<double> &x0, double step, int max_iter) {
	const size_t p = x0.size();
	if (p == 0) {
		return x0;
	}
	const double alpha = 1.0, gamma = 2.0, rho = 0.5, sigma = 0.5;
	std::vector<std::vector<double>> X(p + 1, x0);
	for (size_t i = 0; i < p; ++i) {
		X[i + 1][i] += step;
	}
	std::vector<double> fv(p + 1);
	for (size_t i = 0; i <= p; ++i) {
		fv[i] = f(X[i]);
	}
	std::vector<size_t> ord(p + 1);
	std::vector<double> c(p), xr(p), xe(p), xc(p);
	for (int it = 0; it < max_iter; ++it) {
		for (size_t i = 0; i <= p; ++i) {
			ord[i] = i;
		}
		std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return fv[a] < fv[b]; });
		size_t best = ord.front(), worst = ord.back(), second_worst = ord[p - 1];
		if (std::abs(fv[worst] - fv[best]) <= 1e-12 * (std::abs(fv[best]) + std::abs(fv[worst]) + 1e-300)) {
			break; // simplex collapsed in function value
		}
		// Centroid of every vertex except the worst.
		for (size_t d = 0; d < p; ++d) {
			c[d] = 0.0;
		}
		for (size_t i = 0; i <= p; ++i) {
			if (i == worst) {
				continue;
			}
			for (size_t d = 0; d < p; ++d) {
				c[d] += X[i][d];
			}
		}
		for (size_t d = 0; d < p; ++d) {
			c[d] /= static_cast<double>(p);
		}
		for (size_t d = 0; d < p; ++d) {
			xr[d] = c[d] + alpha * (c[d] - X[worst][d]);
		}
		double fr = f(xr);
		if (fr < fv[best]) {
			for (size_t d = 0; d < p; ++d) {
				xe[d] = c[d] + gamma * (xr[d] - c[d]);
			}
			double fe = f(xe);
			if (fe < fr) {
				X[worst] = xe;
				fv[worst] = fe;
			} else {
				X[worst] = xr;
				fv[worst] = fr;
			}
		} else if (fr < fv[second_worst]) {
			X[worst] = xr;
			fv[worst] = fr;
		} else {
			bool outside = fr < fv[worst];
			for (size_t d = 0; d < p; ++d) {
				xc[d] = outside ? c[d] + rho * (xr[d] - c[d]) : c[d] + rho * (X[worst][d] - c[d]);
			}
			double fc = f(xc);
			if (fc < (outside ? fr : fv[worst])) {
				X[worst] = xc;
				fv[worst] = fc;
			} else {
				for (size_t i = 0; i <= p; ++i) {
					if (i == best) {
						continue;
					}
					for (size_t d = 0; d < p; ++d) {
						X[i][d] = X[best][d] + sigma * (X[i][d] - X[best][d]);
					}
					fv[i] = f(X[i]);
				}
			}
		}
	}
	size_t best = 0;
	for (size_t i = 1; i <= p; ++i) {
		if (fv[i] < fv[best]) {
			best = i;
		}
	}
	return X[best];
}

// Completeness/range validation of a discrete trait against the ML scaffold: exactly the
// tip set (no missing tips, no extras) and every state index < k. Shared by ER and SYM.
void validate_ml_trait(const MLScaffold &s, const std::unordered_map<std::string, uint32_t> &tip_state, uint32_t k) {
	std::vector<std::string> missing;
	for (const auto &entry : s.tip_by_name) {
		if (tip_state.find(std::string(entry.first)) == tip_state.end()) {
			missing.emplace_back(entry.first);
		}
	}
	if (!missing.empty()) {
		throw std::runtime_error("ancestral_ml: no state for tip(s): " + join_names(missing));
	}
	if (tip_state.size() != s.tip_by_name.size()) {
		std::vector<std::string> extra;
		for (const auto &[nm, st] : tip_state) {
			if (s.tip_by_name.find(nm) == s.tip_by_name.end()) {
				extra.push_back(nm);
			}
		}
		throw std::runtime_error("ancestral_ml: state(s) for name(s) that are not a tip in the tree: " +
		                         join_names(extra));
	}
	for (const auto &[nm, st] : tip_state) {
		if (st >= k) {
			throw std::runtime_error("ancestral_ml: tip '" + nm + "' has state index " + std::to_string(st) +
			                         " out of range for k=" + std::to_string(k));
		}
	}
}

// Fit the single ER rate by maximum likelihood: bracket on a coarse grid over log(rate)
// (scaled by the mean branch length so it covers the meaningful r*t range regardless of
// the tree's length units), then refine with golden-section in the bracketing cell (a
// plain golden-section over the whole range can walk onto the r->infinity plateau). `L`
// is used as pruning scratch. Shared by the ER fit and the SYM warm start.
double fit_er_rate(const NewickTree &t, const MLScaffold &s, const std::vector<uint32_t> &obs_state, uint32_t k,
                   std::vector<double> &L) {
	const uint32_t n = static_cast<uint32_t>(t.num_nodes());
	double bsum = 0.0;
	uint32_t bn = 0;
	for (uint32_t i = 0; i < n; ++i) {
		if (i == t.root()) {
			continue;
		}
		bsum += t.branch_length(i);
		bn++;
	}
	double bbar = (bn > 0 && bsum > 0.0) ? bsum / bn : 1.0;
	double xa = std::log(1e-6 / bbar);
	double xb = std::log(1e6 / bbar);
	auto eval = [&](double x) {
		return ml_down_pass(t, s, obs_state, k, std::exp(x), L);
	};
	const int GRID = 60;
	int best_i = 0;
	double best_f = eval(xa);
	for (int i = 1; i <= GRID; ++i) {
		double fx = eval(xa + (xb - xa) * i / GRID);
		if (fx > best_f) {
			best_f = fx;
			best_i = i;
		}
	}
	double lo = xa + (xb - xa) * std::max(0, best_i - 1) / GRID;
	double hi = xa + (xb - xa) * std::min(GRID, best_i + 1) / GRID;
	return std::exp(golden_section_argmax(eval, lo, hi));
}

// Per-trait ER Mk-ML: validate completeness/range, fit (or fix) the rate, then a final
// pruning pass + marginal up-pass.
AncestralMLResult ml_for_trait(const NewickTree &t, const MLScaffold &s,
                               const std::unordered_map<std::string, uint32_t> &tip_state, uint32_t k,
                               std::optional<double> fixed_rate) {
	validate_ml_trait(s, tip_state, k);

	const uint32_t n = static_cast<uint32_t>(t.num_nodes());

	// A single-state trait is monomorphic: posterior 1 for that state at every internal
	// node, logL 0, and the rate is irrelevant (report 0).
	if (k <= 1) {
		AncestralMLResult res;
		res.rate = 0.0;
		res.log_likelihood = 0.0;
		for (uint32_t u = 0; u < n; ++u) {
			if (!t.is_tip(u)) {
				res.states.push_back(AncestralStateML {u, 0, 1.0});
			}
		}
		return res;
	}

	std::vector<uint32_t> obs_state(n, 0);
	for (uint32_t i = 0; i < n; ++i) {
		if (t.is_tip(i)) {
			obs_state[i] = tip_state.at(t.name(i));
		}
	}

	std::vector<double> L(static_cast<size_t>(n) * k);

	double r;
	if (fixed_rate.has_value()) {
		if (!(fixed_rate.value() > 0.0)) {
			throw std::runtime_error("ancestral_ml: fixed rate must be strictly positive");
		}
		r = fixed_rate.value();
	} else {
		r = fit_er_rate(t, s, obs_state, k, L);
	}

	double logL = ml_down_pass(t, s, obs_state, k, r, L); // final pass fills L at r
	if (!std::isfinite(logL)) {
		// The whole-tree likelihood is zero at the chosen rate: some state is impossible
		// everywhere, which for the ER model only happens when an effectively-zero branch
		// (exp(-k*rate*t) rounded to exactly 1.0) hard-pins conflicting tips. Fail loud
		// rather than return -inf/NaN.
		throw std::runtime_error("ancestral_ml: the tree likelihood is zero at the fitted/fixed rate (log-"
		                         "likelihood is not finite); a branch length is too small relative to the rate "
		                         "(exp(-k*rate*t) rounds to 1) -- remove or lengthen effectively-zero-length edges");
	}
	AncestralMLResult res;
	res.rate = r;
	res.log_likelihood = logL;
	std::vector<double> G(static_cast<size_t>(n) * k);
	ml_posteriors(t, s, k, r, L, G, res.states);
	return res;
}

// ---------------------------------------------------------------------------
// SYM (symmetric-rates) Mk-ML: P(t) = exp(Q*t) via a symmetric eigendecomposition, and a
// Nelder-Mead fit of the k(k-1)/2 off-diagonal rates. Shares build_ml_scaffold /
// validate_ml_trait / fit_er_rate and mirrors the ER pruning/rescaling/posterior logic
// (and its two fail-loud guards), but with full k*k matrix messages instead of the ER
// closed form.
// ---------------------------------------------------------------------------

// Symmetric-Q eigendecomposition: Q = U diag(lambda) U^T (U orthonormal, lambda real).
struct SymEigen {
	uint32_t k = 0;
	std::vector<double> lambda; // k eigenvalues
	std::vector<double> U;      // k*k row-major: U[i*k+m] = component i of eigenvector m
};

// Build Q (symmetric, q_ij=q_ji=rates[pair], diagonal = -row sum) from k(k-1)/2 rates in
// pair order (i,j), i<j, and eigendecompose it with Eigen's SelfAdjointEigenSolver.
SymEigen build_sym_eigen(const std::vector<double> &rates, uint32_t k) {
	Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(k, k);
	size_t p = 0;
	for (uint32_t i = 0; i < k; ++i) {
		for (uint32_t j = i + 1; j < k; ++j) {
			double q = rates[p++];
			Q(i, j) = q;
			Q(j, i) = q;
		}
	}
	for (uint32_t i = 0; i < k; ++i) {
		double sum = 0.0;
		for (uint32_t j = 0; j < k; ++j) {
			if (j != i) {
				sum += Q(i, j);
			}
		}
		Q(i, i) = -sum;
	}
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Q);
	SymEigen S;
	S.k = k;
	S.lambda.resize(k);
	S.U.resize(static_cast<size_t>(k) * k);
	for (uint32_t m = 0; m < k; ++m) {
		S.lambda[m] = es.eigenvalues()(m);
	}
	for (uint32_t i = 0; i < k; ++i) {
		for (uint32_t m = 0; m < k; ++m) {
			S.U[static_cast<size_t>(i) * k + m] = es.eigenvectors()(i, m);
		}
	}
	return S;
}

// out[x] = sum_c P(t)[x][c] * v[c] with P(t) = U diag(e^{lambda t}) U^T. Compute
// a[m] = e^{lambda_m t} * (sum_c U[c][m] v[c]), then out[x] = sum_m U[x][m] a[m]: O(k^2),
// no k*k matrix materialized. Because P(t) is symmetric this single routine serves both
// the child->parent down-pass message (Mc[sp] = sum_sc P[sp][sc] Lc[sc]) and the up-pass
// (Gc[sc] = sum_sp P[sp][sc] W[sp] = sum_sp P[sc][sp] W[sp]). `a` is a k-length scratch.
// Every entry of P(t)*v is mathematically non-negative (P(t) >= 0 for a valid rate matrix
// and v >= 0), but the eigenbasis matvec is not a sum of non-negative terms like the ER
// closed form: with U only orthonormal to ~1e-15, an entry whose true value is ~0 (e.g. at
// a zero-length internal edge, where P(0)=U U^T ~ I) can round to a tiny NEGATIVE. We clamp
// at this single shared site so no negative can propagate into a likelihood or an emitted
// probability. A clamped-to-zero message still trips the ml_posteriors_sym Mc>0 guard.
void sym_apply(const SymEigen &S, double t, const double *v, double *out, double *a) {
	const uint32_t k = S.k;
	for (uint32_t m = 0; m < k; ++m) {
		double acc = 0.0;
		for (uint32_t c = 0; c < k; ++c) {
			acc += S.U[static_cast<size_t>(c) * k + m] * v[c];
		}
		a[m] = std::exp(S.lambda[m] * t) * acc;
	}
	for (uint32_t x = 0; x < k; ++x) {
		double acc = 0.0;
		for (uint32_t m = 0; m < k; ++m) {
			acc += S.U[static_cast<size_t>(x) * k + m] * a[m];
		}
		out[x] = std::max(0.0, acc);
	}
}

// A transition-probability model over one tree, shared by the SYM and ARD matrix pruning.
// forward(t, v) = P(t) * v (a child->parent message: M[sp] = sum_sc P(t)[sp][sc] v[sc]);
// backward(t, v) = P(t)^T * v (the up-pass "outside" step: G[sc] = sum_sp P(t)[sp][sc] v[sp]).
// Both outputs are mathematically non-negative (P(t) >= 0, v >= 0); each model clamps its
// output to >= 0 to absorb the sign of floating-point round-off near zero (see sym_apply /
// the ARD note), so no negative can propagate into a likelihood or posterior. A
// clamped-to-zero message still trips the ml_posteriors_matrix Mc>0 guard.

// SYM: P(t) symmetric, so forward == backward, both via the O(k^2) eigenbasis sym_apply
// (no k*k matrix is ever materialized).
struct SymTransition {
	SymEigen S;
	mutable std::vector<double> scratch;
	explicit SymTransition(SymEigen s) : S(std::move(s)), scratch(S.k) {
	}
	void forward(double t, const double *v, double *out) const {
		sym_apply(S, t, v, out, scratch.data());
	}
	void backward(double t, const double *v, double *out) const {
		sym_apply(S, t, v, out, scratch.data());
	}
};

// ARD: Q is a general (non-symmetric) generator, so P(t) = exp(Q*t) is built by a general
// matrix exponential (Eigen MatrixFunctions: scaling-and-squaring Pade, Higham 2005) and
// materialized as a k*k row-major matrix, memoized per distinct branch length WITHIN one
// ArdTransition instance. (A fresh ArdTransition is built for each likelihood evaluation
// during the rate fit, so the cache only amortizes sibling edges that share a length and
// the two passes of the final reconstruction, which share one instance -- there is no
// cross-evaluation reuse.) P is NOT symmetric, so forward uses P and backward uses P^T.
struct ArdTransition {
	uint32_t k;
	Eigen::MatrixXd Q;
	mutable std::map<double, std::vector<double>> pcache; // t -> row-major k*k P(t)

	// rates: the k*(k-1) off-diagonal entries in row order (i,j), i!=j; diagonal = -rowsum.
	ArdTransition(const std::vector<double> &rates, uint32_t k_) : k(k_), Q(Eigen::MatrixXd::Zero(k_, k_)) {
		size_t p = 0;
		for (uint32_t i = 0; i < k; ++i) {
			for (uint32_t j = 0; j < k; ++j) {
				if (j != i) {
					Q(i, j) = rates[p++];
				}
			}
		}
		for (uint32_t i = 0; i < k; ++i) {
			double sum = 0.0;
			for (uint32_t j = 0; j < k; ++j) {
				if (j != i) {
					sum += Q(i, j);
				}
			}
			Q(i, i) = -sum;
		}
	}

	const std::vector<double> &Pt(double t) const {
		auto it = pcache.find(t);
		if (it != pcache.end()) {
			return it->second;
		}
		Eigen::MatrixXd P = (Q * t).exp();
		std::vector<double> Pr(static_cast<size_t>(k) * k);
		for (uint32_t i = 0; i < k; ++i) {
			for (uint32_t j = 0; j < k; ++j) {
				Pr[static_cast<size_t>(i) * k + j] = P(i, j);
			}
		}
		return pcache.emplace(t, std::move(Pr)).first->second;
	}

	void forward(double t, const double *v, double *out) const {
		const std::vector<double> &P = Pt(t);
		for (uint32_t sp = 0; sp < k; ++sp) {
			double acc = 0.0;
			for (uint32_t sc = 0; sc < k; ++sc) {
				acc += P[static_cast<size_t>(sp) * k + sc] * v[sc];
			}
			out[sp] = std::max(0.0, acc);
		}
	}
	void backward(double t, const double *v, double *out) const {
		const std::vector<double> &P = Pt(t);
		for (uint32_t sc = 0; sc < k; ++sc) {
			double acc = 0.0;
			for (uint32_t sp = 0; sp < k; ++sp) {
				acc += P[static_cast<size_t>(sp) * k + sc] * v[sp];
			}
			out[sc] = std::max(0.0, acc);
		}
	}
};

// Felsenstein (1981) pruning down-pass with a general per-edge transition matrix (SYM or
// ARD), with per-node rescaling. Each child message is tr.forward(t_c, L_c) = P(t_c)*L_c.
// Returns logL under the uniform root prior 1/k. Generic over the Transition model.
template <typename Transition>
double ml_down_pass_matrix(const NewickTree &t, const MLScaffold &s, const std::vector<uint32_t> &obs_state, uint32_t k,
                           const Transition &tr, std::vector<double> &L) {
	std::vector<double> Mc(k);
	double logscale = 0.0;
	for (uint32_t u : s.postorder) {
		double *Lu = &L[static_cast<size_t>(u) * k];
		if (t.is_tip(u)) {
			uint32_t o = obs_state[u];
			for (uint32_t st = 0; st < k; ++st) {
				Lu[st] = (st == o) ? 1.0 : 0.0;
			}
			continue;
		}
		for (uint32_t st = 0; st < k; ++st) {
			Lu[st] = 1.0;
		}
		for (uint32_t c : t.children(u)) {
			const double *Lc = &L[static_cast<size_t>(c) * k];
			tr.forward(t.branch_length(c), Lc, Mc.data());
			for (uint32_t st = 0; st < k; ++st) {
				Lu[st] *= Mc[st];
			}
		}
		double m = 0.0;
		for (uint32_t st = 0; st < k; ++st) {
			m = std::max(m, Lu[st]);
		}
		if (m > 0.0) {
			for (uint32_t st = 0; st < k; ++st) {
				Lu[st] /= m;
			}
			logscale += std::log(m);
		}
	}
	const double *Lr = &L[static_cast<size_t>(t.root()) * k];
	double liksum = 0.0;
	for (uint32_t st = 0; st < k; ++st) {
		liksum += Lr[st];
	}
	return std::log(liksum / k) + logscale;
}

// Marginal posterior up-pass (Yang, Kumar & Nei 1995) with a general per-edge transition
// matrix. The sibling product is L_p[sp]/Mc(sp) (Mc = tr.forward, guarded > 0, fail loud),
// and the child "outside" likelihood Gc = P(t_c)^T * W = tr.backward(t_c, W). Generic over
// the Transition model. Posteriors are written for internal nodes only.
template <typename Transition>
void ml_posteriors_matrix(const NewickTree &t, const MLScaffold &s, uint32_t k, const Transition &tr,
                          const std::vector<double> &L, std::vector<double> &G, std::vector<AncestralStateML> &out) {
	const uint32_t root = t.root();
	{
		double *Gr = &G[static_cast<size_t>(root) * k];
		for (uint32_t st = 0; st < k; ++st) {
			Gr[st] = 1.0 / k; // pi
		}
	}
	std::vector<double> W(k), Mc(k);
	for (uint32_t p : s.preorder) {
		if (t.is_tip(p)) {
			continue;
		}
		const double *Lp = &L[static_cast<size_t>(p) * k];
		const double *Gp = &G[static_cast<size_t>(p) * k];
		for (uint32_t c : t.children(p)) {
			const double *Lc = &L[static_cast<size_t>(c) * k];
			double tc = t.branch_length(c);
			tr.forward(tc, Lc, Mc.data()); // message from c to p
			for (uint32_t sp = 0; sp < k; ++sp) {
				if (!(Mc[sp] > 0.0)) {
					throw std::runtime_error(
					    "ancestral_ml: a subtree message underflowed to zero at the fitted rates; a branch "
					    "length is too small relative to the rate to reconstruct posteriors -- remove or "
					    "lengthen effectively-zero-length edges");
				}
				W[sp] = Gp[sp] * Lp[sp] / Mc[sp]; // outside * (siblings' product = L_p / M_c)
			}
			double *Gc = &G[static_cast<size_t>(c) * k];
			tr.backward(tc, W.data(), Gc);
			double m = 0.0;
			for (uint32_t st = 0; st < k; ++st) {
				m = std::max(m, Gc[st]);
			}
			if (m > 0.0) {
				for (uint32_t st = 0; st < k; ++st) {
					Gc[st] /= m;
				}
			}
		}
	}
	const uint32_t n = static_cast<uint32_t>(t.num_nodes());
	for (uint32_t u = 0; u < n; ++u) {
		if (t.is_tip(u)) {
			continue;
		}
		const double *Lu = &L[static_cast<size_t>(u) * k];
		const double *Gu = &G[static_cast<size_t>(u) * k];
		double norm = 0.0;
		for (uint32_t st = 0; st < k; ++st) {
			norm += Lu[st] * Gu[st];
		}
		for (uint32_t st = 0; st < k; ++st) {
			out.push_back(AncestralStateML {u, st, Lu[st] * Gu[st] / norm});
		}
	}
}

// Per-trait SYM Mk-ML: validate, fit the k(k-1)/2 rates by Nelder-Mead (warm-started from
// the ER fit, so the SYM optimum is at least as good as ER), then a final pruning pass +
// marginal up-pass. Reports rate = NaN and the fitted rates vector.
AncestralMLResult sym_for_trait(const NewickTree &t, const MLScaffold &s,
                                const std::unordered_map<std::string, uint32_t> &tip_state, uint32_t k) {
	validate_ml_trait(s, tip_state, k);

	const uint32_t n = static_cast<uint32_t>(t.num_nodes());

	// A single-state trait is monomorphic: posterior 1 for that state at every internal
	// node, logL 0, no rates.
	if (k <= 1) {
		AncestralMLResult res;
		res.rate = std::numeric_limits<double>::quiet_NaN();
		res.log_likelihood = 0.0;
		for (uint32_t u = 0; u < n; ++u) {
			if (!t.is_tip(u)) {
				res.states.push_back(AncestralStateML {u, 0, 1.0});
			}
		}
		return res;
	}

	// SYM fits k*(k-1)/2 rates by simplex search, which becomes unreliable well before a
	// few dozen dimensions. Cap the state count and point high-k traits at ER (which fits a
	// single rate at any k) rather than silently returning a poorly-optimized rate matrix.
	constexpr uint32_t SYM_MAX_STATES = 8; // 8 states -> 28 free rates
	if (k > SYM_MAX_STATES) {
		throw std::runtime_error("ancestral_ml: model 'SYM' fits k*(k-1)/2 rates by simplex search, which is "
		                         "unreliable for many states; this trait has k=" +
		                         std::to_string(k) + " states (> " + std::to_string(SYM_MAX_STATES) +
		                         ") -- use model 'ER' for high-state-count traits");
	}

	std::vector<uint32_t> obs_state(n, 0);
	for (uint32_t i = 0; i < n; ++i) {
		if (t.is_tip(i)) {
			obs_state[i] = tip_state.at(t.name(i));
		}
	}

	std::vector<double> L(static_cast<size_t>(n) * k);
	const size_t np = static_cast<size_t>(k) * (k - 1) / 2;

	// Warm start: the ER rate (all off-diagonals equal) is a valid SYM point.
	double r0 = fit_er_rate(t, s, obs_state, k, L);
	if (!(r0 > 0.0) || !std::isfinite(r0)) {
		r0 = 1.0;
	}

	auto neg_logL = [&](const std::vector<double> &x) -> double {
		std::vector<double> rates(np);
		for (size_t i = 0; i < np; ++i) {
			rates[i] = std::exp(x[i]);
		}
		double ll = ml_down_pass_matrix(t, s, obs_state, k, SymTransition(build_sym_eigen(rates, k)), L);
		return std::isfinite(ll) ? -ll : 1e300; // penalize degenerate regions
	};

	std::vector<double> x0(np, std::log(r0));
	std::vector<double> xbest = nelder_mead_minimize(neg_logL, x0, 0.5, 400);
	xbest = nelder_mead_minimize(neg_logL, xbest, 0.2, 400); // restart to escape early contraction

	std::vector<double> rates(np);
	for (size_t i = 0; i < np; ++i) {
		rates[i] = std::exp(xbest[i]);
	}
	SymTransition tr(build_sym_eigen(rates, k));
	double logL = ml_down_pass_matrix(t, s, obs_state, k, tr, L);
	if (!std::isfinite(logL)) {
		throw std::runtime_error("ancestral_ml: the tree likelihood is zero at the fitted rates (log-likelihood is "
		                         "not finite); a branch length is too small relative to the rate -- remove or "
		                         "lengthen effectively-zero-length edges");
	}
	AncestralMLResult res;
	res.rate = std::numeric_limits<double>::quiet_NaN();
	res.rates = std::move(rates);
	res.log_likelihood = logL;
	std::vector<double> G(static_cast<size_t>(n) * k);
	ml_posteriors_matrix(t, s, k, tr, L, G, res.states);
	return res;
}

// Per-trait ARD Mk-ML: validate, fit the k*(k-1) off-diagonal rates by Nelder-Mead
// (warm-started from the ER fit, so the ARD optimum is at least as good as ER), then a
// final pruning pass + marginal up-pass. Reports rate = NaN and the fitted rates vector.
// ARD is non-reversible: the reconstruction depends on root placement (no pulley
// principle), and P(t) is a general matrix exponential rather than a symmetric one.
AncestralMLResult ard_for_trait(const NewickTree &t, const MLScaffold &s,
                                const std::unordered_map<std::string, uint32_t> &tip_state, uint32_t k) {
	validate_ml_trait(s, tip_state, k);

	const uint32_t n = static_cast<uint32_t>(t.num_nodes());

	if (k <= 1) {
		AncestralMLResult res;
		res.rate = std::numeric_limits<double>::quiet_NaN();
		res.log_likelihood = 0.0;
		for (uint32_t u = 0; u < n; ++u) {
			if (!t.is_tip(u)) {
				res.states.push_back(AncestralStateML {u, 0, 1.0});
			}
		}
		return res;
	}

	// ARD fits k*(k-1) rates -- twice SYM's -- by simplex search, so its state ceiling is
	// tighter (k=6 -> 30 free rates). High-k traits error and are pointed at ER.
	constexpr uint32_t ARD_MAX_STATES = 6; // 6 states -> 30 free rates
	if (k > ARD_MAX_STATES) {
		throw std::runtime_error("ancestral_ml: model 'ARD' fits k*(k-1) rates by simplex search, which is "
		                         "unreliable for many states; this trait has k=" +
		                         std::to_string(k) + " states (> " + std::to_string(ARD_MAX_STATES) +
		                         ") -- use model 'ER' for high-state-count traits");
	}

	std::vector<uint32_t> obs_state(n, 0);
	for (uint32_t i = 0; i < n; ++i) {
		if (t.is_tip(i)) {
			obs_state[i] = tip_state.at(t.name(i));
		}
	}

	std::vector<double> L(static_cast<size_t>(n) * k);
	const size_t np = static_cast<size_t>(k) * (k - 1); // all off-diagonal rates

	double r0 = fit_er_rate(t, s, obs_state, k, L);
	if (!(r0 > 0.0) || !std::isfinite(r0)) {
		r0 = 1.0;
	}

	auto neg_logL = [&](const std::vector<double> &x) -> double {
		std::vector<double> rates(np);
		for (size_t i = 0; i < np; ++i) {
			rates[i] = std::exp(x[i]);
		}
		double ll = ml_down_pass_matrix(t, s, obs_state, k, ArdTransition(rates, k), L);
		return std::isfinite(ll) ? -ll : 1e300;
	};

	// ARD's k*(k-1)-dimensional likelihood surface is multimodal, so a single warm-started
	// simplex is unreliable (risk noted in the plan). Multi-start: the ER warm start plus
	// several deterministic fixed-seed perturbations of it (so the fit is reproducible), each
	// refined with a restart; keep the best. This does not guarantee the global optimum --
	// gradient-based tools may reach boundary (rate->0) solutions this log-space search
	// cannot -- but it materially reduces the chance of parking at a poor local optimum.
	std::mt19937 rng(0x5EED5EEDu);
	std::uniform_real_distribution<double> jitter(-1.5, 1.5); // log-rate units
	const int N_START = 10;
	double best_neg = 1e300;
	std::vector<double> xbest(np, std::log(r0));
	for (int start = 0; start < N_START; ++start) {
		std::vector<double> x0(np, std::log(r0));
		if (start > 0) {
			for (double &xi : x0) {
				xi += jitter(rng); // start 0 is the plain ER warm start
			}
		}
		std::vector<double> xs = nelder_mead_minimize(neg_logL, x0, 0.5, 600);
		xs = nelder_mead_minimize(neg_logL, xs, 0.2, 600); // restart to escape early contraction
		double f = neg_logL(xs);
		if (f < best_neg) {
			best_neg = f;
			xbest = xs;
		}
	}

	std::vector<double> rates(np);
	for (size_t i = 0; i < np; ++i) {
		rates[i] = std::exp(xbest[i]);
	}
	ArdTransition tr(rates, k);
	double logL = ml_down_pass_matrix(t, s, obs_state, k, tr, L);
	if (!std::isfinite(logL)) {
		throw std::runtime_error("ancestral_ml: the tree likelihood is zero at the fitted rates (log-likelihood is "
		                         "not finite); a branch length is too small relative to the rate -- remove or "
		                         "lengthen effectively-zero-length edges");
	}
	AncestralMLResult res;
	res.rate = std::numeric_limits<double>::quiet_NaN();
	res.rates = std::move(rates);
	res.log_likelihood = logL;
	std::vector<double> G(static_cast<size_t>(n) * k);
	ml_posteriors_matrix(t, s, k, tr, L, G, res.states);
	return res;
}

} // namespace

AncestralMLResult NewickTree::ancestral_ml(const std::unordered_map<std::string, uint32_t> &tip_states, uint32_t k,
                                           std::optional<double> fixed_rate,
                                           const std::vector<int64_t> *node_ids) const {
	MLScaffold scaffold = build_ml_scaffold(*this, node_ids);
	if (scaffold.n_tips <= 1) {
		// A degenerate tree (<= 1 tip) has no internal nodes to reconstruct, so it returns
		// empty here -- before per-trait completeness validation runs. (independent_contrasts
		// always validates; for these methods a degenerate tree has nothing to validate
		// against, so the early-out is intentional and harmless.)
		return {};
	}
	return ml_for_trait(*this, scaffold, tip_states, k, fixed_rate);
}

std::vector<AncestralMLResult>
NewickTree::ancestral_ml(const std::vector<std::unordered_map<std::string, uint32_t>> &tip_states_list,
                         const std::vector<uint32_t> &k_list, std::optional<double> fixed_rate,
                         const std::vector<int64_t> *node_ids) const {
	if (tip_states_list.size() != k_list.size()) {
		throw std::runtime_error("ancestral_ml: tip_states_list and k_list must be parallel");
	}
	// Trait-independent validation + traversals computed once, reused per trait.
	MLScaffold scaffold = build_ml_scaffold(*this, node_ids);
	std::vector<AncestralMLResult> results;
	results.reserve(tip_states_list.size());
	for (size_t i = 0; i < tip_states_list.size(); ++i) {
		if (scaffold.n_tips <= 1) {
			results.emplace_back();
		} else {
			results.push_back(ml_for_trait(*this, scaffold, tip_states_list[i], k_list[i], fixed_rate));
		}
	}
	return results;
}

// An empty SYM/ARD result (single tip / empty tree): no internal nodes to reconstruct, and
// rate = NaN to hold the invariant (a rate matrix has no single scalar rate) for any direct
// C++ caller checking std::isnan(res.rate).
static AncestralMLResult empty_matrix_ml_result() {
	AncestralMLResult res;
	res.rate = std::numeric_limits<double>::quiet_NaN();
	return res;
}

AncestralMLResult NewickTree::ancestral_ml_sym(const std::unordered_map<std::string, uint32_t> &tip_states, uint32_t k,
                                               const std::vector<int64_t> *node_ids) const {
	MLScaffold scaffold = build_ml_scaffold(*this, node_ids);
	if (scaffold.n_tips <= 1) {
		return empty_matrix_ml_result();
	}
	return sym_for_trait(*this, scaffold, tip_states, k);
}

std::vector<AncestralMLResult>
NewickTree::ancestral_ml_sym(const std::vector<std::unordered_map<std::string, uint32_t>> &tip_states_list,
                             const std::vector<uint32_t> &k_list, const std::vector<int64_t> *node_ids) const {
	if (tip_states_list.size() != k_list.size()) {
		throw std::runtime_error("ancestral_ml: tip_states_list and k_list must be parallel");
	}
	MLScaffold scaffold = build_ml_scaffold(*this, node_ids);
	std::vector<AncestralMLResult> results;
	results.reserve(tip_states_list.size());
	for (size_t i = 0; i < tip_states_list.size(); ++i) {
		if (scaffold.n_tips <= 1) {
			results.push_back(empty_matrix_ml_result());
		} else {
			results.push_back(sym_for_trait(*this, scaffold, tip_states_list[i], k_list[i]));
		}
	}
	return results;
}

AncestralMLResult NewickTree::ancestral_ml_ard(const std::unordered_map<std::string, uint32_t> &tip_states, uint32_t k,
                                               const std::vector<int64_t> *node_ids) const {
	MLScaffold scaffold = build_ml_scaffold(*this, node_ids);
	if (scaffold.n_tips <= 1) {
		return empty_matrix_ml_result();
	}
	return ard_for_trait(*this, scaffold, tip_states, k);
}

std::vector<AncestralMLResult>
NewickTree::ancestral_ml_ard(const std::vector<std::unordered_map<std::string, uint32_t>> &tip_states_list,
                             const std::vector<uint32_t> &k_list, const std::vector<int64_t> *node_ids) const {
	if (tip_states_list.size() != k_list.size()) {
		throw std::runtime_error("ancestral_ml: tip_states_list and k_list must be parallel");
	}
	MLScaffold scaffold = build_ml_scaffold(*this, node_ids);
	std::vector<AncestralMLResult> results;
	results.reserve(tip_states_list.size());
	for (size_t i = 0; i < tip_states_list.size(); ++i) {
		if (scaffold.n_tips <= 1) {
			results.push_back(empty_matrix_ml_result());
		} else {
			results.push_back(ard_for_trait(*this, scaffold, tip_states_list[i], k_list[i]));
		}
	}
	return results;
}

uint32_t NewickTree::add_node(const std::string &name, double branch_length, std::optional<int64_t> edge_id) {
	if (nodes_.size() >= MAX_NODES) {
		throw std::runtime_error("Tree too large: exceeds maximum of " + std::to_string(MAX_NODES) + " nodes");
	}
	uint32_t idx = static_cast<uint32_t>(nodes_.size());
	nodes_.emplace_back(name, branch_length, edge_id);
	return idx;
}

void NewickTree::set_parent(uint32_t child, uint32_t parent) {
	if (child >= nodes_.size()) {
		throw std::out_of_range("Invalid child node index: " + std::to_string(child));
	}
	if (parent >= nodes_.size()) {
		throw std::out_of_range("Invalid parent node index: " + std::to_string(parent));
	}
	if (child == parent) {
		throw std::invalid_argument("Cannot make node its own parent");
	}

	// If child already has a parent, remove from that parent's children list
	if (nodes_[child].parent != NO_PARENT) {
		remove_child(nodes_[child].parent, child);
	}

	// Set new parent
	nodes_[child].parent = parent;
	nodes_[parent].children.push_back(child);
}

void NewickTree::remove_child(uint32_t parent, uint32_t child) {
	if (parent >= nodes_.size()) {
		throw std::out_of_range("Invalid parent node index: " + std::to_string(parent));
	}
	auto &children = nodes_[parent].children;
	auto it = std::find(children.begin(), children.end(), child);
	if (it != children.end()) {
		children.erase(it);
	}
}

} // namespace miint
