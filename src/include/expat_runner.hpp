// SPDX-License-Identifier: MIT
//
// Shared expat plumbing: parser-create + handler-dispatch + RAII free,
// plus an attribute-lookup helper. Three callers wrap their per-XML state
// machine around RunExpatParse: ena_parser (ParseSampleAttributesXML),
// ena_receipt_parser (ParseReceiptXML), ena_checklist (ParseChecklistXML).
// Header-only so each caller's user-data state struct stays its own type.

#pragma once

#include <expat.h>

#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace miint {

// Look up an attribute value in expat's null-terminated key,value array.
// Returns nullptr if the key is absent.
inline const char *AttrLookup(const char **attrs, const char *key) {
	for (int i = 0; attrs[i]; i += 2) {
		if (std::strcmp(attrs[i], key) == 0) {
			return attrs[i + 1];
		}
	}
	return nullptr;
}

// Run an expat parse over `xml` with the given handlers, passing `state` as
// user data. Throws std::runtime_error on documents larger than 2 GB (expat's
// int-sized limit) and on parse errors, prefixing the error with
// `error_label` (e.g. "ENA receipt").
template <class State>
void RunExpatParse(const std::string &xml, State &state, XML_StartElementHandler start_handler,
                   XML_EndElementHandler end_handler, XML_CharacterDataHandler char_handler, const char *error_label) {
	if (xml.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
		throw std::runtime_error(std::string(error_label) + ": XML document too large for expat (>2 GB)");
	}
	XML_Parser parser = XML_ParserCreate(nullptr);
	if (!parser) {
		throw std::runtime_error(std::string(error_label) + ": failed to create XML parser");
	}
	struct ParserGuard {
		XML_Parser p;
		~ParserGuard() {
			XML_ParserFree(p);
		}
	} guard {parser};

	XML_SetUserData(parser, &state);
	XML_SetElementHandler(parser, start_handler, end_handler);
	if (char_handler) {
		XML_SetCharacterDataHandler(parser, char_handler);
	}

	if (XML_Parse(parser, xml.data(), static_cast<int>(xml.size()), XML_TRUE) == XML_STATUS_ERROR) {
		auto err = std::string(XML_ErrorString(XML_GetErrorCode(parser)));
		throw std::runtime_error(std::string(error_label) + ": XML parse error: " + err);
	}
}

} // namespace miint
