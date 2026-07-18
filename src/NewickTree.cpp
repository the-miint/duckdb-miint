#include "NewickTree.hpp"
#include <algorithm>
#include <array>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <stack>

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
