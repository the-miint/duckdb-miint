#include "NewickTree.hpp"
#include <algorithm>
#include <charconv>
#include <cmath>
#include <stdexcept>
#include <stack>
#include <sstream>
#include <iomanip>

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
		double value;

		auto [ptr, ec] = std::from_chars(num_str.data(), num_str.data() + num_str.size(), value);

		if (ec != std::errc()) {
			throw std::runtime_error("Invalid branch length: '" + std::string(num_str) + "'");
		}
		if (ptr != num_str.data() + num_str.size()) {
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
	result.reserve(nodes_.size() * 20); // Rough estimate

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

			// Branch length - use max_digits10 for proper round-trip
			if (!std::isnan(n.branch_length)) {
				result += ':';
				std::ostringstream oss;
				oss << std::setprecision(std::numeric_limits<double>::max_digits10) << n.branch_length;
				result += oss.str();
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
