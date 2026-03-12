#include "massql_parser.hpp"

#include "formula_parser.hpp"

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
	COLON,
	GREATER,
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

		throw std::runtime_error("MassQL parse error: unexpected character '" + std::string(1, c) + "'");
	}

	Token peek() {
		auto saved = pos_;
		auto tok = next();
		pos_ = saved;
		return tok;
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
    {"TOLERANCEMZ", true},    {"TOLERANCEPPM", true}, {"INTENSITYPERCENT", true},
    {"INTENSITYVALUE", true}, {"EXCLUDED", false},    {"CARDINALITY", true},
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
	throw std::runtime_error("MassQL parse error: unknown aggregation function '" + text + "'");
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

// Parse a number or formula() call as an atomic numeric value
static double parse_atomic_value(Tokenizer &tokenizer) {
	auto tok = tokenizer.next();
	if (tok.type == TokenType::NUMBER) {
		return std::stod(tok.text);
	}
	if (tok.type == TokenType::IDENTIFIER && to_upper(tok.text) == "FORMULA") {
		return parse_formula_call(tokenizer);
	}
	throw std::runtime_error("MassQL parse error: expected number or formula(), got '" + tok.text + "'");
}

// Parse a single condition value (number, X-variable expression, or formula)
// Handles: 157.0857+10, X+164.9, X-formula(Fe), 2*(X-formula(Fe)), formula(Fe)
static ConditionValue parse_condition_value(Tokenizer &tokenizer) {
	ConditionValue val;
	auto tok = tokenizer.next();

	if (tok.type == TokenType::NUMBER) {
		double num_val = std::stod(tok.text);
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
		auto it = QUALIFIER_MAP.find(qname);
		if (it == QUALIFIER_MAP.end()) {
			throw std::runtime_error("MassQL parse error: unknown qualifier '" + name_tok.text + "'");
		}

		Qualifier qual;
		qual.name = qname;

		if (it->second) {
			// Qualifier takes a value: = or >
			auto op = tokenizer.next();
			if (op.type != TokenType::EQUALS && op.type != TokenType::GREATER) {
				throw std::runtime_error("MassQL parse error: expected '=' or '>' after qualifier name");
			}
			qual.op = (op.type == TokenType::GREATER) ? QualifierOp::GREATER_THAN : QualifierOp::EQUALS;
			auto val_tok = tokenizer.next();

			// Handle range(min=a,max=b) for CARDINALITY
			if (val_tok.type == TokenType::IDENTIFIER && to_lower(val_tok.text) == "range") {
				// Expect: (min=N,max=N) — consume character by character through the tokenizer
				auto lp = tokenizer.next(); // (
				if (lp.type != TokenType::LPAREN) {
					throw std::runtime_error("MassQL parse error: expected '(' after range");
				}
				// Parse min=N
				auto min_kw = tokenizer.next();
				if (to_lower(min_kw.text) != "min") {
					throw std::runtime_error("MassQL parse error: expected 'min' in range()");
				}
				tokenizer.next(); // consume =
				auto min_val = tokenizer.next();
				qual.value = std::stod(min_val.text);

				// Consume comma between min and max
				auto comma = tokenizer.next();
				if (comma.type != TokenType::COMMA) {
					throw std::runtime_error("MassQL parse error: expected ',' in range()");
				}

				auto max_kw = tokenizer.next();
				if (to_lower(max_kw.text) != "max") {
					throw std::runtime_error("MassQL parse error: expected 'max' in range()");
				}
				tokenizer.next(); // consume =
				auto max_val = tokenizer.next();
				qual.max_value = std::stod(max_val.text);

				auto rp = tokenizer.next(); // )
				if (rp.type != TokenType::RPAREN) {
					throw std::runtime_error("MassQL parse error: expected ')' after range values");
				}
			} else if (val_tok.type != TokenType::NUMBER) {
				throw std::runtime_error("MassQL parse error: expected number for qualifier value");
			} else {
				qual.value = std::stod(val_tok.text);
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
		throw std::runtime_error("MassQL parse error: unknown condition field '" + field_tok.text + "'");
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
		auto rparen = tokenizer.next();
		if (rparen.type != TokenType::RPAREN) {
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
