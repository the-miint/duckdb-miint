#include "massql_parser.hpp"

#include "formula_parser.hpp"
#include "mass_tables.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <unordered_map>

namespace miint {

// ── Token types ──────────────────────────────────────────────────────────────

enum class TokenType {
	IDENTIFIER,
	NUMBER,
	LPAREN,
	RPAREN,
	EQUALS,
	PLUS,
	MINUS,
	STAR,
	SLASH,
	COLON,
	GREATER,
	LESS,
	COMMA,
	END_OF_INPUT
};

struct Token {
	TokenType type;
	std::string text;
};

// ── Tokenizer ────────────────────────────────────────────────────────────────

class Tokenizer {
public:
	explicit Tokenizer(const std::string &input) : input_(input), pos_(0) {
	}

	Token next() {
		skip_whitespace();
		if (pos_ >= input_.size()) {
			return {TokenType::END_OF_INPUT, ""};
		}

		char c = input_[pos_];

		if (c == '(') {
			pos_++;
			return {TokenType::LPAREN, "("};
		}
		if (c == ')') {
			pos_++;
			return {TokenType::RPAREN, ")"};
		}
		if (c == '=') {
			pos_++;
			return {TokenType::EQUALS, "="};
		}
		if (c == '+') {
			pos_++;
			return {TokenType::PLUS, "+"};
		}
		if (c == '-') {
			pos_++;
			return {TokenType::MINUS, "-"};
		}
		if (c == '*') {
			pos_++;
			return {TokenType::STAR, "*"};
		}
		if (c == ':') {
			pos_++;
			return {TokenType::COLON, ":"};
		}
		if (c == '>') {
			pos_++;
			return {TokenType::GREATER, ">"};
		}
		if (c == '<') {
			pos_++;
			return {TokenType::LESS, "<"};
		}
		if (c == '/') {
			pos_++;
			return {TokenType::SLASH, "/"};
		}
		if (c == ',') {
			pos_++;
			return {TokenType::COMMA, ","};
		}

		if (std::isdigit(c) || c == '.') {
			return read_number();
		}

		if (std::isalpha(c) || c == '_') {
			return read_identifier();
		}

		throw std::runtime_error("MassQL parse error at position " + std::to_string(pos_) + ": unexpected character '" +
		                         std::string(1, c) + "'");
	}

	Token peek() {
		auto saved = pos_;
		auto tok = next();
		pos_ = saved;
		return tok;
	}

	size_t pos() const {
		return pos_;
	}

private:
	void skip_whitespace() {
		while (pos_ < input_.size() && std::isspace(input_[pos_])) {
			pos_++;
		}
	}

	Token read_number() {
		size_t start = pos_;
		bool has_digit = false;
		while (pos_ < input_.size() && (std::isdigit(input_[pos_]) || input_[pos_] == '.')) {
			if (std::isdigit(input_[pos_])) {
				has_digit = true;
			}
			pos_++;
		}
		if (!has_digit) {
			throw std::runtime_error("MassQL parse error: invalid number (no digits)");
		}
		return {TokenType::NUMBER, input_.substr(start, pos_ - start)};
	}

	Token read_identifier() {
		size_t start = pos_;
		while (pos_ < input_.size() && (std::isalnum(input_[pos_]) || input_[pos_] == '_')) {
			pos_++;
		}
		return {TokenType::IDENTIFIER, input_.substr(start, pos_ - start)};
	}

	std::string input_;
	size_t pos_;
};

// ── Helpers ──────────────────────────────────────────────────────────────────

static void expect_equals(Tokenizer &tokenizer, const std::string &context) {
	auto tok = tokenizer.next();
	if (tok.type != TokenType::EQUALS) {
		throw std::runtime_error("MassQL parse error: expected '=' after '" + context + "'");
	}
}

static std::string to_upper(const std::string &s) {
	std::string result = s;
	std::transform(result.begin(), result.end(), result.begin(), ::toupper);
	return result;
}

static std::string to_lower(const std::string &s) {
	std::string result = s;
	std::transform(result.begin(), result.end(), result.begin(), ::tolower);
	return result;
}

// ── "Did you mean?" suggestions ──────────────────────────────────────────────

static size_t edit_distance(const std::string &a, const std::string &b) {
	size_t m = a.size(), n = b.size();
	std::vector<size_t> prev(n + 1), curr(n + 1);
	for (size_t j = 0; j <= n; j++) {
		prev[j] = j;
	}
	for (size_t i = 1; i <= m; i++) {
		curr[0] = i;
		for (size_t j = 1; j <= n; j++) {
			size_t cost = (std::toupper(a[i - 1]) == std::toupper(b[j - 1])) ? 0 : 1;
			curr[j] = std::min({prev[j] + 1, curr[j - 1] + 1, prev[j - 1] + cost});
		}
		std::swap(prev, curr);
	}
	return prev[n];
}

template <typename MapType>
static std::string suggest_from_map(const std::string &input, const MapType &map, size_t max_dist = 3) {
	std::string best;
	size_t best_dist = max_dist + 1;
	for (const auto &entry : map) {
		size_t d = edit_distance(input, entry.first);
		if (d < best_dist) {
			best_dist = d;
			best = entry.first;
		}
	}
	if (!best.empty()) {
		return ". Did you mean '" + best + "'?";
	}
	return "";
}

static std::string suggest_agg_function(const std::string &input) {
	static const std::vector<std::string> agg_names = {"SCANINFO", "SCANSUM",    "SCANNUM",
	                                                   "SCANMZ",   "SCANMAXINT", "SCANRANGESUM"};
	std::string best;
	size_t best_dist = 4; // max_dist + 1 (max_dist=3, matching suggest_from_map)
	for (const auto &name : agg_names) {
		size_t d = edit_distance(input, name);
		if (d < best_dist) {
			best_dist = d;
			best = name;
		}
	}
	if (!best.empty()) {
		return ". Did you mean '" + to_lower(best) + "()'?";
	}
	return "";
}

// ── Field name lookup ────────────────────────────────────────────────────────

static const std::unordered_map<std::string, ConditionField> FIELD_MAP = {
    {"MS1MZ", ConditionField::MS1MZ},     {"MS2PROD", ConditionField::MS2PROD},
    {"MS2MZ", ConditionField::MS2PROD}, // alias
    {"MS2PREC", ConditionField::MS2PREC}, {"MS2NL", ConditionField::MS2NL},
    {"RTMIN", ConditionField::RTMIN},     {"RTMAX", ConditionField::RTMAX},
    {"SCANMIN", ConditionField::SCANMIN}, {"SCANMAX", ConditionField::SCANMAX},
    {"CHARGE", ConditionField::CHARGE},   {"POLARITY", ConditionField::POLARITY},
};

// ── Qualifier name set ───────────────────────────────────────────────────────

static const std::unordered_map<std::string, bool> QUALIFIER_MAP = {
    {"TOLERANCEMZ", true},
    {"TOLERANCEPPM", true},
    {"INTENSITYPERCENT", true},
    {"INTENSITYVALUE", true},
    {"INTENSITYTICPERCENT", true},
    {"MASSDEFECT", true},
    {"INTENSITYMATCH", true},
    {"INTENSITYMATCHREFERENCE", false},
    {"INTENSITYMATCHPERCENT", true},
    {"OTHERSCAN", true},
    {"EXCLUDED", false},
    {"CARDINALITY", true},
    {"MATCHCOUNT", true}, // alias for CARDINALITY
};

// ── Parser ───────────────────────────────────────────────────────────────────

static DataType parse_data_type(const std::string &text) {
	auto upper = to_upper(text);
	if (upper == "MS1DATA") {
		return DataType::MS1DATA;
	}
	if (upper == "MS2DATA") {
		return DataType::MS2DATA;
	}
	throw std::runtime_error("MassQL parse error: unknown data type '" + text + "'");
}

static AggFunction parse_agg_function(const std::string &text) {
	auto upper = to_upper(text);
	if (upper == "SCANINFO") {
		return AggFunction::SCANINFO;
	}
	if (upper == "SCANSUM") {
		return AggFunction::SCANSUM;
	}
	if (upper == "SCANNUM") {
		return AggFunction::SCANNUM;
	}
	if (upper == "SCANMZ") {
		return AggFunction::SCANMZ;
	}
	if (upper == "SCANMAXINT") {
		return AggFunction::SCANMAXINT;
	}
	if (upper == "SCANRANGESUM") {
		return AggFunction::SCANRANGESUM;
	}
	throw std::runtime_error("MassQL parse error: unknown aggregation function '" + text + "'" +
	                         suggest_agg_function(upper));
}

static double safe_stod(const std::string &text) {
	try {
		return std::stod(text);
	} catch (const std::exception &) {
		throw std::runtime_error("MassQL parse error: invalid number '" + text + "'");
	}
}

static int safe_stoi(const std::string &text) {
	try {
		return std::stoi(text);
	} catch (const std::exception &) {
		throw std::runtime_error("MassQL parse error: invalid integer '" + text + "'");
	}
}

// Parse formula(FORMULA_STRING) call — "formula" keyword already consumed
static double parse_formula_call(Tokenizer &tokenizer) {
	auto lp = tokenizer.next();
	if (lp.type != TokenType::LPAREN) {
		throw std::runtime_error("MassQL parse error: expected '(' after 'formula'");
	}
	auto formula_tok = tokenizer.next();
	if (formula_tok.type != TokenType::IDENTIFIER) {
		throw std::runtime_error("MassQL parse error: expected formula string inside formula()");
	}
	auto rp = tokenizer.next();
	if (rp.type != TokenType::RPAREN) {
		throw std::runtime_error("MassQL parse error: expected ')' after formula string");
	}
	try {
		return ParseFormula(formula_tok.text);
	} catch (const std::exception &e) {
		throw std::runtime_error("MassQL parse error: " + std::string(e.what()));
	}
}

// Parse aminoaciddelta(SEQUENCE) call — "aminoaciddelta" keyword already consumed
static double parse_aminoaciddelta_call(Tokenizer &tokenizer) {
	auto lp = tokenizer.next();
	if (lp.type != TokenType::LPAREN) {
		throw std::runtime_error("MassQL parse error: expected '(' after 'aminoaciddelta'");
	}
	auto seq_tok = tokenizer.next();
	if (seq_tok.type != TokenType::IDENTIFIER) {
		throw std::runtime_error("MassQL parse error: expected amino acid sequence inside aminoaciddelta()");
	}
	auto rp = tokenizer.next();
	if (rp.type != TokenType::RPAREN) {
		throw std::runtime_error("MassQL parse error: expected ')' after amino acid sequence");
	}
	try {
		return AminoacidDeltaMass(seq_tok.text);
	} catch (const std::exception &e) {
		throw std::runtime_error("MassQL parse error: " + std::string(e.what()));
	}
}

// Parse peptide(SEQUENCE,charge=N,ion=X) call — "peptide" keyword already consumed
static double parse_peptide_call(Tokenizer &tokenizer) {
	auto lp = tokenizer.next();
	if (lp.type != TokenType::LPAREN) {
		throw std::runtime_error("MassQL parse error: expected '(' after 'peptide'");
	}
	auto seq_tok = tokenizer.next();
	if (seq_tok.type != TokenType::IDENTIFIER) {
		throw std::runtime_error("MassQL parse error: expected amino acid sequence inside peptide()");
	}
	std::string sequence = seq_tok.text;

	// Parse ,charge=N
	auto comma1 = tokenizer.next();
	if (comma1.type != TokenType::COMMA) {
		throw std::runtime_error("MassQL parse error: expected ',' after sequence in peptide()");
	}
	auto charge_kw = tokenizer.next();
	if (to_lower(charge_kw.text) != "charge") {
		throw std::runtime_error("MassQL parse error: expected 'charge' parameter in peptide()");
	}
	expect_equals(tokenizer, "charge");
	auto charge_val = tokenizer.next();
	int charge = safe_stoi(charge_val.text);

	// Parse ,ion=X
	auto comma2 = tokenizer.next();
	if (comma2.type != TokenType::COMMA) {
		throw std::runtime_error("MassQL parse error: expected ',' after charge in peptide()");
	}
	auto ion_kw = tokenizer.next();
	if (to_lower(ion_kw.text) != "ion") {
		throw std::runtime_error("MassQL parse error: expected 'ion' parameter in peptide()");
	}
	expect_equals(tokenizer, "ion");
	auto ion_val = tokenizer.next();
	if (ion_val.text.size() != 1) {
		throw std::runtime_error("MassQL parse error: ion type must be a single letter (a,b,c,x,y,z)");
	}
	char ion_type = static_cast<char>(std::tolower(ion_val.text[0]));

	auto rp = tokenizer.next();
	if (rp.type != TokenType::RPAREN) {
		throw std::runtime_error("MassQL parse error: expected ')' after peptide() parameters");
	}
	try {
		return PeptideFragmentMass(sequence, charge, ion_type);
	} catch (const std::exception &e) {
		throw std::runtime_error("MassQL parse error: " + std::string(e.what()));
	}
}

// Parse a number, formula(), aminoaciddelta(), or peptide() call as an atomic numeric value
static double parse_atomic_value(Tokenizer &tokenizer) {
	auto tok = tokenizer.next();
	if (tok.type == TokenType::NUMBER) {
		return safe_stod(tok.text);
	}
	if (tok.type == TokenType::IDENTIFIER && to_upper(tok.text) == "FORMULA") {
		return parse_formula_call(tokenizer);
	}
	if (tok.type == TokenType::IDENTIFIER && to_upper(tok.text) == "AMINOACIDDELTA") {
		return parse_aminoaciddelta_call(tokenizer);
	}
	if (tok.type == TokenType::IDENTIFIER && to_upper(tok.text) == "PEPTIDE") {
		return parse_peptide_call(tokenizer);
	}
	throw std::runtime_error("MassQL parse error: expected number or formula(), got '" + tok.text + "'");
}

// Parse a single condition value (number, X-variable expression, or formula)
// Handles: 157.0857+10, X+164.9, X-formula(Fe), 2*(X-formula(Fe)), formula(Fe)
static ConditionValue parse_condition_value(Tokenizer &tokenizer) {
	ConditionValue val;
	auto tok = tokenizer.next();

	if (tok.type == TokenType::NUMBER) {
		double num_val = safe_stod(tok.text);

		// Handle division: NUMBER/NUMBER folds to constant
		if (tokenizer.peek().type == TokenType::SLASH) {
			tokenizer.next(); // consume /
			double divisor = parse_atomic_value(tokenizer);
			if (divisor == 0.0) {
				throw std::runtime_error("MassQL parse error: division by zero");
			}
			num_val /= divisor;
		}

		auto peek = tokenizer.peek();

		if (peek.type == TokenType::STAR) {
			// Coefficient: NUMBER * (X-expression)
			tokenizer.next(); // consume *
			auto lp = tokenizer.next();
			if (lp.type != TokenType::LPAREN) {
				throw std::runtime_error("MassQL parse error: expected '(' after coefficient '*'");
			}
			auto inner_tok = tokenizer.next();
			if (inner_tok.type != TokenType::IDENTIFIER || to_upper(inner_tok.text) != "X") {
				throw std::runtime_error("MassQL parse error: expected 'X' inside coefficient expression");
			}
			val.has_x_variable = true;
			val.x_coefficient = num_val;
			val.constant_value = 0.0;

			peek = tokenizer.peek();
			if (peek.type == TokenType::PLUS) {
				tokenizer.next(); // consume +
				val.constant_value = num_val * parse_atomic_value(tokenizer);
			} else if (peek.type == TokenType::MINUS) {
				tokenizer.next(); // consume -
				val.constant_value = -num_val * parse_atomic_value(tokenizer);
			}

			auto rp = tokenizer.next();
			if (rp.type != TokenType::RPAREN) {
				throw std::runtime_error("MassQL parse error: expected ')' after coefficient expression");
			}
		} else if (peek.type == TokenType::PLUS) {
			// Constant folding: number + (number or formula)
			tokenizer.next(); // consume +
			val.constant_value = num_val + parse_atomic_value(tokenizer);
			val.has_x_variable = false;
		} else if (peek.type == TokenType::MINUS) {
			// Constant folding: number - (number or formula)
			tokenizer.next(); // consume -
			val.constant_value = num_val - parse_atomic_value(tokenizer);
			val.has_x_variable = false;
		} else {
			val.constant_value = num_val;
			val.has_x_variable = false;
		}
	} else if (tok.type == TokenType::IDENTIFIER && to_upper(tok.text) == "ANY") {
		val.has_any_wildcard = true;
		val.has_x_variable = false;
		val.constant_value = 0.0;
	} else if (tok.type == TokenType::IDENTIFIER && to_upper(tok.text) == "X") {
		val.has_x_variable = true;
		val.x_coefficient = 1.0;
		val.constant_value = 0.0;

		// Check for X+offset or X-offset (number or formula)
		auto peek = tokenizer.peek();
		if (peek.type == TokenType::PLUS) {
			tokenizer.next(); // consume +
			val.constant_value = parse_atomic_value(tokenizer);
		} else if (peek.type == TokenType::MINUS) {
			tokenizer.next(); // consume -
			val.constant_value = -parse_atomic_value(tokenizer);
		}
	} else if (tok.type == TokenType::IDENTIFIER && to_upper(tok.text) == "FORMULA") {
		// formula(FORMULA_STRING) as constant value
		val.has_x_variable = false;
		val.constant_value = parse_formula_call(tokenizer);

		auto peek = tokenizer.peek();
		if (peek.type == TokenType::PLUS) {
			tokenizer.next(); // consume +
			val.constant_value += parse_atomic_value(tokenizer);
		} else if (peek.type == TokenType::MINUS) {
			tokenizer.next(); // consume -
			val.constant_value -= parse_atomic_value(tokenizer);
		}
	} else if (tok.type == TokenType::IDENTIFIER && to_upper(tok.text) == "AMINOACIDDELTA") {
		// aminoaciddelta(SEQUENCE) as constant value
		val.has_x_variable = false;
		val.constant_value = parse_aminoaciddelta_call(tokenizer);

		auto peek = tokenizer.peek();
		if (peek.type == TokenType::PLUS) {
			tokenizer.next(); // consume +
			val.constant_value += parse_atomic_value(tokenizer);
		} else if (peek.type == TokenType::MINUS) {
			tokenizer.next(); // consume -
			val.constant_value -= parse_atomic_value(tokenizer);
		}
	} else if (tok.type == TokenType::IDENTIFIER && to_upper(tok.text) == "PEPTIDE") {
		// peptide(SEQUENCE,charge=N,ion=X) as constant value
		val.has_x_variable = false;
		val.constant_value = parse_peptide_call(tokenizer);

		auto peek = tokenizer.peek();
		if (peek.type == TokenType::PLUS) {
			tokenizer.next(); // consume +
			val.constant_value += parse_atomic_value(tokenizer);
		} else if (peek.type == TokenType::MINUS) {
			tokenizer.next(); // consume -
			val.constant_value -= parse_atomic_value(tokenizer);
		}
	} else {
		throw std::runtime_error("MassQL parse error: expected number or X variable, got '" + tok.text + "'");
	}

	return val;
}

// Parse qualifiers after a condition value: :QUALIFIER=value or :QUALIFIER>value or :EXCLUDED
static void parse_qualifiers(Tokenizer &tokenizer, Condition &cond) {
	while (tokenizer.peek().type == TokenType::COLON) {
		tokenizer.next(); // consume ':'

		auto name_tok = tokenizer.next();
		if (name_tok.type != TokenType::IDENTIFIER) {
			throw std::runtime_error("MassQL parse error: expected qualifier name after ':'");
		}

		auto qname = to_upper(name_tok.text);
		if (qname == "MOBILITY") {
			throw std::runtime_error(
			    "MassQL parse error: MOBILITY is not supported (read_mzml does not provide ion mobility data)");
		}
		auto it = QUALIFIER_MAP.find(qname);
		if (it == QUALIFIER_MAP.end()) {
			throw std::runtime_error("MassQL parse error: unknown qualifier '" + name_tok.text + "'" +
			                         suggest_from_map(qname, QUALIFIER_MAP));
		}

		Qualifier qual;
		// Normalize MATCHCOUNT → CARDINALITY
		qual.name = (qname == "MATCHCOUNT") ? "CARDINALITY" : qname;

		if (it->second) {
			// Qualifier takes a value: =, >, or <
			auto op = tokenizer.next();
			if (op.type != TokenType::EQUALS && op.type != TokenType::GREATER && op.type != TokenType::LESS) {
				throw std::runtime_error("MassQL parse error: expected '=', '>', or '<' after qualifier name");
			}
			if (op.type == TokenType::GREATER) {
				qual.op = QualifierOp::GREATER_THAN;
			} else if (op.type == TokenType::LESS) {
				qual.op = QualifierOp::LESS_THAN;
			} else {
				qual.op = QualifierOp::EQUALS;
			}
			auto val_tok = tokenizer.next();

			// Handle INTENSITYMATCH=Y-expression
			if (qname == "INTENSITYMATCH" && val_tok.type == TokenType::IDENTIFIER && to_upper(val_tok.text) == "Y") {
				// Parse Y-expression: Y, Y*number, Y*(number), Y*(number+number*X)
				auto peek = tokenizer.peek();
				if (peek.type == TokenType::STAR) {
					tokenizer.next(); // consume *
					peek = tokenizer.peek();
					if (peek.type == TokenType::LPAREN) {
						tokenizer.next(); // consume (
						int open_parens = 1;
						auto num_tok = tokenizer.next();
						qual.y_expr_constant = safe_stod(num_tok.text);
						peek = tokenizer.peek();
						if (peek.type == TokenType::PLUS) {
							// Y*(constant+coeff*X)
							tokenizer.next(); // consume +
							peek = tokenizer.peek();
							if (peek.type == TokenType::LPAREN) {
								tokenizer.next(); // consume inner (
								open_parens++;
							}
							auto coeff_tok = tokenizer.next();
							double coeff = safe_stod(coeff_tok.text);
							auto star = tokenizer.next();  // consume *
							auto x_tok = tokenizer.next(); // consume X
							if (to_upper(x_tok.text) != "X") {
								throw std::runtime_error("MassQL parse error: expected X in Y expression");
							}
							qual.y_expr_x_coeff = coeff;
							qual.y_expr_has_x = true;
							// consume exactly the matching closing parens
							for (int i = 0; i < open_parens; i++) {
								auto rp = tokenizer.next();
								if (rp.type != TokenType::RPAREN) {
									throw std::runtime_error("MassQL parse error: expected ')' in Y expression");
								}
							}
						} else {
							auto rp = tokenizer.next(); // consume )
							if (rp.type != TokenType::RPAREN) {
								throw std::runtime_error("MassQL parse error: expected ')' in Y expression");
							}
						}
					} else {
						auto num_tok = tokenizer.next();
						qual.y_expr_constant = safe_stod(num_tok.text);
					}
				}
				// else: bare Y, constant stays 1.0
			} else if (val_tok.type == TokenType::IDENTIFIER &&
			           (to_lower(val_tok.text) == "range" || to_lower(val_tok.text) == "massdefect")) {
				// Handle range(min=a,max=b) for CARDINALITY and massdefect(min=a,max=b) for MASSDEFECT
				// Expect: (min=N,max=N) — consume character by character through the tokenizer
				auto lp = tokenizer.next(); // (
				if (lp.type != TokenType::LPAREN) {
					throw std::runtime_error("MassQL parse error: expected '(' after " + val_tok.text);
				}
				// Parse min=N
				auto min_kw = tokenizer.next();
				if (to_lower(min_kw.text) != "min") {
					throw std::runtime_error("MassQL parse error: expected 'min' in range()");
				}
				expect_equals(tokenizer, "min");
				auto min_val = tokenizer.next();
				qual.value = safe_stod(min_val.text);

				// Consume comma between min and max
				auto comma = tokenizer.next();
				if (comma.type != TokenType::COMMA) {
					throw std::runtime_error("MassQL parse error: expected ',' in range()");
				}

				auto max_kw = tokenizer.next();
				if (to_lower(max_kw.text) != "max") {
					throw std::runtime_error("MassQL parse error: expected 'max' in range()");
				}
				expect_equals(tokenizer, "max");
				auto max_val = tokenizer.next();
				qual.max_value = safe_stod(max_val.text);

				auto rp = tokenizer.next(); // )
				if (rp.type != TokenType::RPAREN) {
					throw std::runtime_error("MassQL parse error: expected ')' after range values");
				}
			} else if (val_tok.type == TokenType::IDENTIFIER && to_lower(val_tok.text) == "rtrange") {
				// Handle rtrange(left=N,right=N) for OTHERSCAN
				auto lp = tokenizer.next(); // (
				if (lp.type != TokenType::LPAREN) {
					throw std::runtime_error("MassQL parse error: expected '(' after rtrange");
				}
				auto left_kw = tokenizer.next();
				if (to_lower(left_kw.text) != "left") {
					throw std::runtime_error("MassQL parse error: expected 'left' in rtrange()");
				}
				expect_equals(tokenizer, "left");
				auto left_val = tokenizer.next();
				qual.value = safe_stod(left_val.text);

				auto comma = tokenizer.next();
				if (comma.type != TokenType::COMMA) {
					throw std::runtime_error("MassQL parse error: expected ',' in rtrange()");
				}

				auto right_kw = tokenizer.next();
				if (to_lower(right_kw.text) != "right") {
					throw std::runtime_error("MassQL parse error: expected 'right' in rtrange()");
				}
				expect_equals(tokenizer, "right");
				auto right_val = tokenizer.next();
				qual.max_value = safe_stod(right_val.text);

				auto rp = tokenizer.next(); // )
				if (rp.type != TokenType::RPAREN) {
					throw std::runtime_error("MassQL parse error: expected ')' after rtrange values");
				}
			} else if (val_tok.type != TokenType::NUMBER) {
				throw std::runtime_error("MassQL parse error: expected number for qualifier value");
			} else {
				qual.value = safe_stod(val_tok.text);
			}
		}

		cond.qualifiers.push_back(qual);
	}
}

// Parse a single condition: FIELD=value[:qualifiers]
static Condition parse_condition(Tokenizer &tokenizer) {
	Condition cond;

	auto field_tok = tokenizer.next();
	if (field_tok.type != TokenType::IDENTIFIER) {
		throw std::runtime_error("MassQL parse error: expected condition field name");
	}

	auto field_upper = to_upper(field_tok.text);
	auto it = FIELD_MAP.find(field_upper);
	if (it == FIELD_MAP.end()) {
		throw std::runtime_error("MassQL parse error: unknown condition field '" + field_tok.text + "'" +
		                         suggest_from_map(field_upper, FIELD_MAP));
	}
	cond.field = it->second;

	// Consume '='
	auto eq = tokenizer.next();
	if (eq.type != TokenType::EQUALS) {
		throw std::runtime_error("MassQL parse error: expected '=' after field name");
	}

	// POLARITY takes a string value
	if (cond.field == ConditionField::POLARITY) {
		auto val_tok = tokenizer.next();
		if (val_tok.type != TokenType::IDENTIFIER) {
			throw std::runtime_error("MassQL parse error: expected polarity value (Positive/Negative)");
		}
		cond.string_value = to_lower(val_tok.text);
	} else {
		// Check for parenthesized OR list: (val1 OR val2 OR ...)
		auto peek = tokenizer.peek();
		if (peek.type == TokenType::LPAREN) {
			tokenizer.next(); // consume '('
			cond.values.push_back(parse_condition_value(tokenizer));
			while (true) {
				auto next = tokenizer.peek();
				if (next.type == TokenType::IDENTIFIER && to_upper(next.text) == "OR") {
					tokenizer.next(); // consume OR
					cond.values.push_back(parse_condition_value(tokenizer));
				} else {
					break;
				}
			}
			auto rp = tokenizer.next();
			if (rp.type != TokenType::RPAREN) {
				throw std::runtime_error("MassQL parse error: expected ')' after OR list");
			}
		} else {
			cond.values.push_back(parse_condition_value(tokenizer));
		}
	}

	parse_qualifiers(tokenizer, cond);

	return cond;
}

// Parse condition list: condition [AND condition]*
static std::vector<Condition> parse_condition_list(Tokenizer &tokenizer) {
	std::vector<Condition> conditions;
	conditions.push_back(parse_condition(tokenizer));

	while (true) {
		auto peek = tokenizer.peek();
		if (peek.type == TokenType::IDENTIFIER && to_upper(peek.text) == "AND") {
			tokenizer.next(); // consume AND
			conditions.push_back(parse_condition(tokenizer));
		} else {
			break;
		}
	}

	return conditions;
}

MassQLQuery MassQLParser::parse(const std::string &query) {
	if (query.empty()) {
		throw std::runtime_error("MassQL parse error: empty query string");
	}

	Tokenizer tokenizer(query);
	MassQLQuery result;

	// Expect QUERY keyword
	auto tok = tokenizer.next();
	if (tok.type != TokenType::IDENTIFIER || to_upper(tok.text) != "QUERY") {
		throw std::runtime_error("MassQL parse error: expected QUERY keyword, got '" + tok.text + "'");
	}

	// Next token: either aggregation function or data type
	tok = tokenizer.next();
	if (tok.type != TokenType::IDENTIFIER) {
		throw std::runtime_error("MassQL parse error: expected data type or aggregation function");
	}

	// Check if this is an aggregation function: identifier followed by '('
	auto peek = tokenizer.peek();
	if (peek.type == TokenType::LPAREN) {
		result.agg_function = parse_agg_function(tok.text);
		tokenizer.next(); // consume '('
		tok = tokenizer.next();
		if (tok.type != TokenType::IDENTIFIER) {
			throw std::runtime_error("MassQL parse error: expected data type inside aggregation function");
		}
		result.data_type = parse_data_type(tok.text);
		// Check for optional parameters: scanrangesum(MS2DATA, TOLERANCE=0.5)
		auto next_tok = tokenizer.next();
		if (next_tok.type == TokenType::COMMA && result.agg_function == AggFunction::SCANRANGESUM) {
			auto param_name = tokenizer.next();
			if (param_name.type != TokenType::IDENTIFIER || to_upper(param_name.text) != "TOLERANCE") {
				throw std::runtime_error("MassQL parse error: expected TOLERANCE parameter in scanrangesum()");
			}
			auto eq = tokenizer.next();
			if (eq.type != TokenType::EQUALS) {
				throw std::runtime_error("MassQL parse error: expected '=' after TOLERANCE");
			}
			auto val = tokenizer.next();
			if (val.type != TokenType::NUMBER) {
				throw std::runtime_error("MassQL parse error: expected number for TOLERANCE value");
			}
			result.scanrangesum_tolerance = safe_stod(val.text);
			if (result.scanrangesum_tolerance <= 0.0) {
				throw std::runtime_error("MassQL parse error: TOLERANCE must be > 0");
			}
			next_tok = tokenizer.next(); // should be ')'
		}
		if (next_tok.type != TokenType::RPAREN) {
			throw std::runtime_error("MassQL parse error: expected ')' after data type");
		}
	} else {
		result.data_type = parse_data_type(tok.text);
	}

	// Optional WHERE clause
	peek = tokenizer.peek();
	if (peek.type == TokenType::IDENTIFIER && to_upper(peek.text) == "WHERE") {
		tokenizer.next(); // consume WHERE
		result.where_conditions = parse_condition_list(tokenizer);
	}

	// Optional FILTER clause
	peek = tokenizer.peek();
	if (peek.type == TokenType::IDENTIFIER && to_upper(peek.text) == "FILTER") {
		tokenizer.next(); // consume FILTER
		result.filter_conditions = parse_condition_list(tokenizer);
	}

	// Reject unexpected trailing tokens
	auto trailing = tokenizer.peek();
	if (trailing.type != TokenType::END_OF_INPUT) {
		throw std::runtime_error("MassQL parse error: unexpected token '" + trailing.text + "' after query");
	}

	return result;
}

} // namespace miint
