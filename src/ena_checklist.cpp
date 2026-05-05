// SPDX-License-Identifier: MIT
//
// ENA checklist parser + validator + registry (Phase 8 Step 8b).
// Pure-data: no DuckDB linkage, so test/cpp/test_ena_checklist.cpp links the
// .cpp directly into the unit-test binary.

#include "ena_checklist.hpp"

#include <expat.h>

#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace miint {

namespace {

const char *AttrLookup(const char **attrs, const char *key) {
	for (int i = 0; attrs[i]; i += 2) {
		if (std::strcmp(attrs[i], key) == 0) {
			return attrs[i + 1];
		}
	}
	return nullptr;
}

// Parser state machine. The checklist XML has a fixed shape:
//   CHECKLIST_SET > CHECKLIST > DESCRIPTOR > FIELD_GROUP > FIELD
//                                                          > LABEL
//                                                          > NAME
//                                                          > MANDATORY
//                                                          > UNITS > UNIT
//                                                          > FIELD_TYPE
//                                                            > TEXT_CHOICE_FIELD
//                                                              > TEXT_VALUE > VALUE
// Tracking the element stack lets us disambiguate <NAME> inside CHECKLIST
// (the descriptor name) vs inside FIELD (the snake-case identifier).
struct ParserState {
	ChecklistDef checklist;
	std::vector<std::string> element_stack;

	// Per-field accumulator. Active only between <FIELD> and </FIELD>.
	bool in_field = false;
	ChecklistFieldDef current;

	// Text capture for any element whose text content matters.
	bool capture_text = false;
	std::string text_buffer;

	// Did we ever see a CHECKLIST element?
	bool saw_checklist = false;
};

bool StackEndsWith(const std::vector<std::string> &stack, std::initializer_list<const char *> suffix) {
	if (stack.size() < suffix.size()) {
		return false;
	}
	auto it = stack.end() - static_cast<std::ptrdiff_t>(suffix.size());
	for (const auto *s : suffix) {
		if (*it != s) {
			return false;
		}
		++it;
	}
	return true;
}

ChecklistFieldMandatory ParseMandatory(const std::string &raw) {
	if (raw == "mandatory") {
		return ChecklistFieldMandatory::MANDATORY;
	}
	if (raw == "recommended") {
		return ChecklistFieldMandatory::RECOMMENDED;
	}
	// Default to OPTIONAL on unknown / "optional" — ENA's enum is fixed but
	// being permissive on the read side is harmless.
	return ChecklistFieldMandatory::OPTIONAL;
}

void XMLCALL StartElement(void *user_data, const XML_Char *name, const XML_Char **attrs) {
	auto *st = static_cast<ParserState *>(user_data);
	st->element_stack.emplace_back(name);

	if (std::strcmp(name, "CHECKLIST") == 0) {
		st->saw_checklist = true;
		if (const auto *acc = AttrLookup(attrs, "accession")) {
			st->checklist.accession = acc;
		}
	} else if (std::strcmp(name, "FIELD") == 0 && StackEndsWith(st->element_stack, {"FIELD_GROUP", "FIELD"})) {
		st->in_field = true;
		st->current = ChecklistFieldDef {};
	} else if (st->in_field) {
		// Inside a FIELD: capture text for the children we care about.
		if (StackEndsWith(st->element_stack, {"FIELD", "LABEL"}) ||
		    StackEndsWith(st->element_stack, {"FIELD", "NAME"}) ||
		    StackEndsWith(st->element_stack, {"FIELD", "MANDATORY"}) ||
		    StackEndsWith(st->element_stack, {"UNITS", "UNIT"}) ||
		    StackEndsWith(st->element_stack, {"TEXT_VALUE", "VALUE"})) {
			st->capture_text = true;
			st->text_buffer.clear();
		}
	} else if (StackEndsWith(st->element_stack, {"DESCRIPTOR", "LABEL"})) {
		// Top-level checklist label (e.g. "GSC MIxS human gut").
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

	// Drain any captured text for this element BEFORE popping the stack so
	// StackEndsWith continues to match the in-field path.
	if (st->capture_text) {
		if (st->in_field) {
			if (StackEndsWith(st->element_stack, {"FIELD", "LABEL"})) {
				st->current.label = st->text_buffer;
			} else if (StackEndsWith(st->element_stack, {"FIELD", "NAME"})) {
				st->current.name = st->text_buffer;
			} else if (StackEndsWith(st->element_stack, {"FIELD", "MANDATORY"})) {
				st->current.mandatory = ParseMandatory(st->text_buffer);
			} else if (StackEndsWith(st->element_stack, {"UNITS", "UNIT"})) {
				st->current.allowed_units.push_back(st->text_buffer);
			} else if (StackEndsWith(st->element_stack, {"TEXT_VALUE", "VALUE"})) {
				st->current.allowed_values.push_back(st->text_buffer);
			}
		} else if (StackEndsWith(st->element_stack, {"DESCRIPTOR", "LABEL"})) {
			st->checklist.label = st->text_buffer;
		}
		st->text_buffer.clear();
		st->capture_text = false;
	}

	if (std::strcmp(name, "FIELD") == 0 && st->in_field) {
		// Only push fields with at least a LABEL — defensive against
		// malformed XML (a FIELD without a NAME/LABEL is nonsense).
		if (!st->current.label.empty()) {
			st->checklist.fields.push_back(std::move(st->current));
		}
		st->in_field = false;
		st->current = ChecklistFieldDef {};
	}

	if (!st->element_stack.empty()) {
		st->element_stack.pop_back();
	}
}

std::string TrimTrailingSlashes(std::string s) {
	while (!s.empty() && s.back() == '/') {
		s.pop_back();
	}
	return s;
}

} // namespace

ChecklistDef ParseChecklistXML(const std::string &xml) {
	if (xml.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
		throw std::runtime_error("ENA checklist: XML document too large for expat (>2 GB)");
	}
	XML_Parser parser = XML_ParserCreate(nullptr);
	if (!parser) {
		throw std::runtime_error("ENA checklist: failed to create XML parser");
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
		throw std::runtime_error("ENA checklist: XML parse error: " + err);
	}

	if (!state.saw_checklist) {
		throw std::runtime_error("ENA checklist: no <CHECKLIST> element in document");
	}

	return std::move(state.checklist);
}

std::vector<ChecklistValidationIssue>
ValidateAttributesAgainstChecklist(const ChecklistDef &checklist,
                                   const std::vector<std::pair<std::string, std::string>> &attributes,
                                   const std::vector<std::pair<std::string, std::string>> &units) {
	std::vector<ChecklistValidationIssue> issues;

	// Build label → FieldDef pointer for O(1) lookup in the user-attribute pass.
	std::unordered_map<std::string, const ChecklistFieldDef *> by_label;
	by_label.reserve(checklist.fields.size());
	for (const auto &f : checklist.fields) {
		by_label.emplace(f.label, &f);
	}

	// DuckDB MAP rejects duplicate keys at parse time, so the duplicate
	// case shouldn't reach us from SQL. We still build a defensive presence
	// map (last-wins) so callers passing pre-built attribute lists in C++
	// don't get spurious "missing mandatory" issues for an attribute that
	// also has a stale empty entry.
	std::unordered_map<std::string, std::string> attr_values;
	for (const auto &kv : attributes) {
		attr_values[kv.first] = kv.second;
	}

	std::unordered_map<std::string, std::string> unit_values;
	for (const auto &kv : units) {
		unit_values[kv.first] = kv.second;
	}

	// Pass 1: every mandatory field must be present with non-empty value.
	for (const auto &f : checklist.fields) {
		if (f.mandatory != ChecklistFieldMandatory::MANDATORY) {
			continue;
		}
		auto it = attr_values.find(f.label);
		if (it == attr_values.end() || it->second.empty()) {
			issues.push_back({f.label, "checklist '" + checklist.accession + "' marks field '" + f.label +
			                               "' as mandatory but it is missing or empty"});
		}
	}

	// Pass 2: each user attribute is checked against the checklist.
	for (const auto &kv : attributes) {
		const auto &label = kv.first;
		const auto &value = kv.second;
		auto it = by_label.find(label);
		if (it == by_label.end()) {
			issues.push_back({label, "attribute '" + label + "' is not in checklist '" + checklist.accession + "'"});
			continue;
		}
		const auto *field = it->second;
		// 2a: when a field declares units, the user must supply a unit value
		// from the allowed set. Empty user value → no unit is required (the
		// mandatory-presence pass above already flagged it if needed).
		if (!field->allowed_units.empty() && !value.empty()) {
			auto u_it = unit_values.find(label);
			if (u_it == unit_values.end() || u_it->second.empty()) {
				std::string allowed;
				for (size_t i = 0; i < field->allowed_units.size(); i++) {
					if (i > 0) {
						allowed += ", ";
					}
					allowed += "'" + field->allowed_units[i] + "'";
				}
				issues.push_back({label, "field '" + label + "' requires a unit (one of " + allowed +
				                             ") to be set in attribute_units"});
			} else {
				bool ok = false;
				for (const auto &u : field->allowed_units) {
					if (u == u_it->second) {
						ok = true;
						break;
					}
				}
				if (!ok) {
					std::string allowed;
					for (size_t i = 0; i < field->allowed_units.size(); i++) {
						if (i > 0) {
							allowed += ", ";
						}
						allowed += "'" + field->allowed_units[i] + "'";
					}
					issues.push_back({label, "field '" + label + "' has unit '" + u_it->second +
					                             "' which is not in allowed set " + allowed});
				}
			}
		}
		// 2b: controlled vocabulary check.
		if (!field->allowed_values.empty() && !value.empty()) {
			bool ok = false;
			for (const auto &v : field->allowed_values) {
				if (v == value) {
					ok = true;
					break;
				}
			}
			if (!ok) {
				issues.push_back({label, "value '" + value + "' for '" + label + "' is not in checklist CV (e.g. '" +
				                             field->allowed_values.front() + "', ...)"});
			}
		}
	}

	return issues;
}

std::string BuildChecklistFetchURL(const std::string &base, const std::string &accession) {
	if (accession.empty()) {
		throw std::invalid_argument("BuildChecklistFetchURL: accession must be non-empty");
	}
	// ENA accessions are alphanumeric (digits + letters). Reject anything
	// else so a stray value can't inject path segments or query parameters.
	for (unsigned char c : accession) {
		const bool ok = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9');
		if (!ok) {
			throw std::invalid_argument("BuildChecklistFetchURL: accession must be alphanumeric, got '" + accession +
			                            "'");
		}
	}
	return TrimTrailingSlashes(base) + "/" + accession;
}

std::string ResolveChecklistBaseFromEnv() {
	const char *override_base = std::getenv("MIINT_ENA_CHECKLIST_URL_BASE");
	const std::string base = (override_base != nullptr && *override_base != '\0')
	                             ? std::string(override_base)
	                             : std::string(DEFAULT_CHECKLIST_URL_BASE);
	return TrimTrailingSlashes(base);
}

ChecklistRegistry &ChecklistRegistry::Instance() {
	static ChecklistRegistry singleton;
	return singleton;
}

const ChecklistDef &ChecklistRegistry::GetOrFetch(const std::string &accession, const Fetcher &fetcher) {
	{
		std::lock_guard<std::mutex> guard(mu);
		auto it = cache.find(accession);
		if (it != cache.end()) {
			return it->second;
		}
	}
	if (!fetcher) {
		throw std::invalid_argument("ChecklistRegistry::GetOrFetch: fetcher must be set");
	}
	// Fetch outside the lock so we don't block other accessions while one
	// network round-trip is in flight. Two concurrent fetches for the same
	// accession would each hit the network — accepted: the cost is one
	// extra fetch in the rare race window, and the alternative
	// (per-accession lock map) buys us nothing for typical workloads.
	const auto url = BuildChecklistFetchURL(ResolveChecklistBaseFromEnv(), accession);
	const auto body = fetcher(url);
	auto parsed = ParseChecklistXML(body);

	// Override the parsed checklist's accession with the user-supplied key
	// so error messages reference the accession the user asked for, not
	// whatever the upstream XML labelled itself with (the test mock's
	// permissive fallback returns accession="ERC000000-mock" regardless).
	parsed.accession = accession;
	std::lock_guard<std::mutex> guard(mu);
	auto inserted = cache.emplace(accession, std::move(parsed));
	return inserted.first->second;
}

void ChecklistRegistry::Clear() {
	std::lock_guard<std::mutex> guard(mu);
	cache.clear();
}

} // namespace miint
