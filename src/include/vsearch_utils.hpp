#pragma once

#include <stdexcept>
#include <string>

namespace miint {

// vsearch result structs use fixed 1024-char label buffers. Validate at input
// boundaries to throw early rather than silently truncate.
inline void validate_label_length(const std::string &label, const char *context) {
	if (label.size() > 1023) {
		throw std::runtime_error(std::string(context) + ": ID '" + label.substr(0, 40) +
		                         "...' exceeds 1023-character limit (got " + std::to_string(label.size()) + ")");
	}
}

// Convert RNA alphabet (U/u) to DNA (T/t) for vsearch compatibility.
// vsearch operates on DNA internally; RNA sequences must be normalized.
inline std::string normalize_rna(const std::string &seq) {
	std::string out = seq;
	for (auto &c : out) {
		if (c == 'U') {
			c = 'T';
		} else if (c == 'u') {
			c = 't';
		}
	}
	return out;
}

} // namespace miint
