#include "sam_style_tags.hpp"

#include <stdexcept>

namespace miint {

namespace {

// Character classes and value grammars from the SAM specification's optional-field
// section (section 1.5), which GFA 1 adopts for its tags.

bool IsDigit(char c) {
	return c >= '0' && c <= '9';
}

bool IsAlpha(char c) {
	return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
}

// [!-~]
bool IsPrintable(char c) {
	return c >= '!' && c <= '~';
}

// [ !-~]
bool IsPrintableOrSpace(char c) {
	return c >= ' ' && c <= '~';
}

} // namespace

bool AllPrintableOrSpace(std::string_view s) {
	for (char c : s) {
		if (!IsPrintableOrSpace(c)) {
			return false;
		}
	}
	return true;
}

namespace {

// [-+]?[0-9]+
bool IsInteger(std::string_view s) {
	size_t i = 0;
	if (i < s.size() && (s[i] == '-' || s[i] == '+')) {
		++i;
	}
	if (i == s.size()) {
		return false;
	}
	for (; i < s.size(); ++i) {
		if (!IsDigit(s[i])) {
			return false;
		}
	}
	return true;
}

// `s` equals the lowercase `word`, ignoring the case of ASCII letters.
bool EqualsLowercase(std::string_view s, std::string_view word) {
	if (s.size() != word.size()) {
		return false;
	}
	for (size_t i = 0; i < s.size(); ++i) {
		const char c = (s[i] >= 'A' && s[i] <= 'Z') ? static_cast<char>(s[i] - 'A' + 'a') : s[i];
		if (c != word[i]) {
			return false;
		}
	}
	return true;
}

// [-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?, or an optional sign then nan, inf or infinity in
// any case.
// The specification's grammar stops at the first form, but htslib parses floats with strtod,
// which accepts the named ones, and tools write them. strtod's nan(chars) and hex floats are
// not accepted.
bool IsFloat(std::string_view s) {
	size_t i = 0;
	const size_t n = s.size();
	if (i < n && (s[i] == '-' || s[i] == '+')) {
		++i;
	}
	const std::string_view named = s.substr(i);
	if (EqualsLowercase(named, "nan") || EqualsLowercase(named, "inf") || EqualsLowercase(named, "infinity")) {
		return true;
	}
	size_t int_digits = 0;
	while (i < n && IsDigit(s[i])) {
		++i;
		++int_digits;
	}
	if (i < n && s[i] == '.') {
		++i;
		size_t frac_digits = 0;
		while (i < n && IsDigit(s[i])) {
			++i;
			++frac_digits;
		}
		if (frac_digits == 0) {
			return false;
		}
	} else if (int_digits == 0) {
		return false;
	}
	if (i < n && (s[i] == 'e' || s[i] == 'E')) {
		++i;
		if (i < n && (s[i] == '-' || s[i] == '+')) {
			++i;
		}
		size_t exp_digits = 0;
		while (i < n && IsDigit(s[i])) {
			++i;
			++exp_digits;
		}
		if (exp_digits == 0) {
			return false;
		}
	}
	return i == n;
}

// ([0-9A-F][0-9A-F])*
bool IsHexBytes(std::string_view s) {
	if (s.size() % 2 != 0) {
		return false;
	}
	for (char c : s) {
		if (!IsDigit(c) && !(c >= 'A' && c <= 'F')) {
			return false;
		}
	}
	return true;
}

// [cCsSiIf](,<element>)*: integer elements for the integer subtypes, IsFloat ones for f.
bool IsNumericArray(std::string_view s) {
	if (s.empty()) {
		return false;
	}
	const char subtype = s[0];
	bool is_float;
	switch (subtype) {
	case 'c':
	case 'C':
	case 's':
	case 'S':
	case 'i':
	case 'I':
		is_float = false;
		break;
	case 'f':
		is_float = true;
		break;
	default:
		return false;
	}
	size_t pos = 1;
	while (pos < s.size()) {
		if (s[pos] != ',') {
			return false;
		}
		size_t next = s.find(',', pos + 1);
		auto element = s.substr(pos + 1, next == std::string_view::npos ? std::string_view::npos : next - pos - 1);
		if (is_float ? !IsFloat(element) : !IsInteger(element)) {
			return false;
		}
		pos = next == std::string_view::npos ? s.size() : next;
	}
	return true;
}

bool IsValidValue(char type, std::string_view value) {
	switch (type) {
	case 'A':
		return value.size() == 1 && IsPrintable(value[0]);
	case 'i':
		return IsInteger(value);
	case 'f':
		return IsFloat(value);
	case 'Z':
		return AllPrintableOrSpace(value);
	case 'J':
		return !value.empty() && AllPrintableOrSpace(value);
	case 'H':
		return IsHexBytes(value);
	case 'B':
		return IsNumericArray(value);
	default:
		return false;
	}
}

// [A-Za-z_][A-Za-z0-9_]*, for typed and untyped fields alike (why: TagSyntax::Typed).
bool IsTagName(std::string_view key) {
	if (key.empty() || !(IsAlpha(key[0]) || key[0] == '_')) {
		return false;
	}
	for (char c : key) {
		if (!(IsAlpha(c) || IsDigit(c) || c == '_')) {
			return false;
		}
	}
	return true;
}

std::string Quoted(std::string_view field) {
	return "'" + std::string(field) + "'";
}

} // namespace

void ParseSamStyleTags(const std::vector<std::string_view> &fields, size_t first, TagSyntax syntax,
                       std::vector<SamStyleTag> &out) {
	out.clear();
	for (size_t i = first; i < fields.size(); ++i) {
		std::string_view field = fields[i];
		const size_t colon = field.find(':');
		const bool named = colon != std::string_view::npos && IsTagName(field.substr(0, colon));
		// NAME:TYPE:VALUE: a known type letter between the first two colons.
		const bool typed = named && colon + 2 < field.size() &&
		                   std::string_view("AifZJHB").find(field[colon + 1]) != std::string_view::npos &&
		                   field[colon + 2] == ':';
		std::string_view value;
		if (typed) {
			value = field.substr(colon + 3);
			if (!IsValidValue(field[colon + 1], value)) {
				throw std::runtime_error(Quoted(field) + " has a value that is not valid for type '" +
				                         std::string(1, field[colon + 1]) + "'");
			}
		} else {
			if (syntax != TagSyntax::TypedOrUntyped || !named || !AllPrintableOrSpace(field.substr(colon + 1))) {
				throw std::runtime_error(Quoted(field) + " is not a tag (expected TAG:TYPE:VALUE" +
				                         std::string(syntax == TagSyntax::TypedOrUntyped ? " or key:value" : "") + ")");
			}
			value = field.substr(colon + 1);
		}
		// Every unnamed field has thrown above.
		const std::string_view key = field.substr(0, colon);
		for (const auto &existing : out) {
			if (existing.key == key) {
				throw std::runtime_error("tag " + Quoted(key) + " appears more than once");
			}
		}
		out.push_back({std::string(key), std::string(value)});
	}
}

} // namespace miint
