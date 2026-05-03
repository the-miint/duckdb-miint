// Phase 3 GREEN: implemented inline. See ena_receipt_parser.hpp.

#include "ena_receipt_parser.hpp"

#include <expat.h>

#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace miint {

namespace {

// Object element names per SRA.receipt.xsd line 75–86. Membership test for
// "is this element a (potentially accessioned) object container?".
bool IsObjectElement(const char *name) {
	static const char *kObjectElements[] = {"ANALYSIS", "EXPERIMENT", "RUN",     "SAMPLE",  "SAMPLEGROUP",
	                                         "STUDY",    "DAC",        "POLICY",  "DATASET", "PROJECT",
	                                         "CHECKLIST"};
	for (const auto *n : kObjectElements) {
		if (std::strcmp(name, n) == 0) {
			return true;
		}
	}
	return false;
}

const char *AttrLookup(const char **attrs, const char *key) {
	for (int i = 0; attrs[i]; i += 2) {
		if (std::strcmp(attrs[i], key) == 0) {
			return attrs[i + 1];
		}
	}
	return nullptr;
}

struct ParserState {
	ENAReceipt receipt;
	// When inside <ERROR>/<INFO>/<ACTIONS>, accumulate text.
	// SRA.receipt.xsd does NOT permit those elements to nest under each other,
	// so a single boolean is sufficient — no element-stack needed. If a future
	// XSD relaxes this, switch to a depth counter.
	std::string text_buffer;
	bool capture_text = false;
};

void XMLCALL StartElement(void *user_data, const XML_Char *name, const XML_Char **attrs) {
	auto *st = static_cast<ParserState *>(user_data);

	if (std::strcmp(name, "RECEIPT") == 0) {
		const char *success = AttrLookup(attrs, "success");
		st->receipt.success = success && (std::strcmp(success, "true") == 0);
		if (auto *v = AttrLookup(attrs, "receiptDate")) {
			st->receipt.receipt_date = v;
		}
		if (auto *v = AttrLookup(attrs, "submissionFile")) {
			st->receipt.submission_file = v;
		}
	} else if (std::strcmp(name, "SUBMISSION") == 0) {
		// SUBMISSION is special — its accession lives on the receipt root.
		if (auto *acc = AttrLookup(attrs, "accession")) {
			st->receipt.submission_accession = acc;
		}
	} else if (IsObjectElement(name)) {
		ENAObjectReceipt obj;
		obj.object_type = name;
		if (auto *v = AttrLookup(attrs, "alias")) {
			obj.alias = v;
		}
		if (auto *v = AttrLookup(attrs, "accession")) {
			obj.accession = v;
		}
		if (auto *v = AttrLookup(attrs, "status")) {
			obj.status = v;
		}
		if (auto *v = AttrLookup(attrs, "holdUntilDate")) {
			obj.hold_until_date = v;
		}
		st->receipt.objects.push_back(std::move(obj));
	} else if (std::strcmp(name, "EXT_ID") == 0) {
		// Attach to the most recent object on the receipt (ext_id is always
		// nested under one of the object elements per the XSD).
		if (!st->receipt.objects.empty()) {
			ENAExtId ext;
			if (auto *v = AttrLookup(attrs, "accession")) {
				ext.accession = v;
			}
			if (auto *v = AttrLookup(attrs, "type")) {
				ext.type = v;
			}
			st->receipt.objects.back().ext_ids.push_back(std::move(ext));
		}
	} else if (std::strcmp(name, "ERROR") == 0 || std::strcmp(name, "INFO") == 0 ||
	           std::strcmp(name, "ACTIONS") == 0) {
		st->capture_text = true;
		st->text_buffer.clear();
	}
}

void XMLCALL CharData(void *user_data, const XML_Char *s, int len) {
	auto *st = static_cast<ParserState *>(user_data);
	if (st->capture_text) {
		st->text_buffer.append(s, static_cast<size_t>(len));
	}
}

void XMLCALL EndElement(void *user_data, const XML_Char *name) {
	auto *st = static_cast<ParserState *>(user_data);

	if (st->capture_text) {
		if (std::strcmp(name, "ERROR") == 0) {
			st->receipt.errors.push_back(st->text_buffer);
		} else if (std::strcmp(name, "INFO") == 0) {
			st->receipt.info_messages.push_back(st->text_buffer);
		} else if (std::strcmp(name, "ACTIONS") == 0) {
			st->receipt.actions.push_back(st->text_buffer);
		}
		st->text_buffer.clear();
		st->capture_text = false;
	}
}

} // namespace

ENAReceipt ParseReceiptXML(const std::string &xml) {
	if (xml.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
		throw std::runtime_error("ENA receipt: XML document too large for expat (>2 GB)");
	}
	XML_Parser parser = XML_ParserCreate(nullptr);
	if (!parser) {
		throw std::runtime_error("ENA receipt: failed to create XML parser");
	}
	struct ParserGuard {
		XML_Parser p;
		~ParserGuard() {
			XML_ParserFree(p);
		}
	} guard {parser};

	ParserState state;
	XML_SetUserData(parser, &state);
	XML_SetElementHandler(parser, StartElement, EndElement);
	XML_SetCharacterDataHandler(parser, CharData);

	if (XML_Parse(parser, xml.data(), static_cast<int>(xml.size()), XML_TRUE) == XML_STATUS_ERROR) {
		auto err = std::string(XML_ErrorString(XML_GetErrorCode(parser)));
		throw std::runtime_error("ENA receipt: XML parse error: " + err);
	}

	return std::move(state.receipt);
}

} // namespace miint
