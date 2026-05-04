// SPDX-License-Identifier: MIT
//
// ENA Webin V2 receipt parser (XML form).
//
// Parses RECEIPT documents per SRA.receipt.xsd into a structured ENAReceipt.
// XML is the canonical form (XSD-governed); JSON receipt mirrors lowercase
// plurals and is intentionally not handled here — Phase 4+ will POST with
// Accept: application/xml so receipts always arrive as XML.

#pragma once

#include <cctype>
#include <string>
#include <vector>

namespace miint {

// Case-insensitive ASCII equality, used for matching ENA EXT_ID/@type values
// (the receipt XSD declares it as a free-form string and observed casings
// include "biosample", "BioSample", "BIOSAMPLE"). Inline so the pure-data
// layer doesn't have to link against the rest of the receipt parser.
inline bool EqualsIgnoreCase(const std::string &a, const std::string &b) {
	if (a.size() != b.size()) {
		return false;
	}
	for (size_t i = 0; i < a.size(); ++i) {
		if (std::tolower(static_cast<unsigned char>(a[i])) != std::tolower(static_cast<unsigned char>(b[i]))) {
			return false;
		}
	}
	return true;
}

struct ENAExtId {
	std::string accession;
	std::string type; // "biosample", "study", "Project", ... per SRA.receipt.xsd EXT_ID enum
};

struct ENAObjectReceipt {
	std::string object_type; // "PROJECT", "SAMPLE", "EXPERIMENT", "RUN", "ANALYSIS", ...
	std::string alias;       // required
	std::string accession;   // empty if the object failed validation
	std::string status;      // "PRIVATE", "DRAFT", "CANCELLED", "PUBLIC", ...
	std::string hold_until_date;
	std::vector<ENAExtId> ext_ids;
};

struct ENAReceipt {
	bool success = false; // defaults to false; only set true on explicit success="true"
	std::string receipt_date;
	std::string submission_file;
	std::string submission_accession; // ERA... — from <SUBMISSION accession="...">
	std::vector<ENAObjectReceipt> objects;
	std::vector<std::string> errors;
	std::vector<std::string> info_messages;
	std::vector<std::string> actions; // "ADD", "HOLD", "VALIDATE", ...
};

// Parse an ENA Webin V2 RECEIPT XML document. Throws std::runtime_error on
// malformed XML; never throws on missing/unexpected fields (returns the
// best-effort parse, with success=false if the receipt is unusable).
ENAReceipt ParseReceiptXML(const std::string &xml);

} // namespace miint
