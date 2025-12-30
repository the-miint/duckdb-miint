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
					if (c == '(' || c == ')' || c == ',' || c == ':' || c == ';' || c == '{' || c == '}' ||
					    c == '\'' || c == '"' || c == '[' || c == ']' ||
					    std::isspace(static_cast<unsigned char>(c))) {
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
		throw std::runtime_error("Too many nodes: " + std::to_string(nodes.size()) +
		                         " exceeds maximum of " + std::to_string(MAX_NODES));
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
				throw std::runtime_error("Node " + std::to_string(node.node_id) +
				                         " references non-existent parent " + std::to_string(node.parent_id.value()));
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
				throw std::runtime_error("Cycle detected involving node " +
				                         std::to_string(nodes[child_idx].node_id));
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
			throw std::runtime_error("Unknown edge_id " + std::to_string(p.edge_id) +
			                         " for fragment '" + p.fragment_id + "'");
		}
		if (p.distal_length < 0) {
			throw std::runtime_error("Negative distal_length " + std::to_string(p.distal_length) +
			                         " for fragment '" + p.fragment_id + "'");
		}
		if (p.pendant_length < 0) {
			throw std::runtime_error("Negative pendant_length " + std::to_string(p.pendant_length) +
			                         " for fragment '" + p.fragment_id + "'");
		}
		// Validate distal_length against edge length
		uint32_t edge_node = edge_index.at(p.edge_id);
		double edge_length = nodes_[edge_node].branch_length;
		if (!std::isnan(edge_length) && p.distal_length > edge_length) {
			throw std::runtime_error("distal_length " + std::to_string(p.distal_length) +
			                         " exceeds edge length " + std::to_string(edge_length) +
			                         " for fragment '" + p.fragment_id + "'");
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
