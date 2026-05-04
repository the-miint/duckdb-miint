// Phase 3 GREEN: implemented inline. See ena_envelope_builder.hpp.

#include "ena_envelope_builder.hpp"

#include <cstdio>
#include <stdexcept>

namespace miint {

namespace {

// JSON string escaping per RFC 8259 §7. Treats input as opaque bytes; multibyte
// UTF-8 characters pass through unchanged. Control bytes (< 0x20) are escaped
// as \uXXXX; named escapes used where defined.
void AppendJsonString(std::string &out, const std::string &s) {
	out.push_back('"');
	for (auto c : s) {
		const auto u = static_cast<unsigned char>(c);
		switch (u) {
		case '"':
			out.append("\\\"");
			break;
		case '\\':
			out.append("\\\\");
			break;
		case '\b':
			out.append("\\b");
			break;
		case '\f':
			out.append("\\f");
			break;
		case '\n':
			out.append("\\n");
			break;
		case '\r':
			out.append("\\r");
			break;
		case '\t':
			out.append("\\t");
			break;
		default:
			if (u < 0x20) {
				char buf[8];
				std::snprintf(buf, sizeof(buf), "\\u%04x", u);
				out.append(buf);
			} else {
				out.push_back(c);
			}
		}
	}
	out.push_back('"');
}

const char *ActionName(ENAAction a) {
	switch (a) {
	case ENAAction::ADD:
		return "ADD";
	case ENAAction::MODIFY:
		return "MODIFY";
	case ENAAction::CANCEL:
		return "CANCEL";
	case ENAAction::HOLD:
		return "HOLD";
	case ENAAction::RELEASE:
		return "RELEASE";
	case ENAAction::VALIDATE:
		return "VALIDATE";
	}
	throw std::logic_error("ENA envelope: unhandled ENAAction value");
}

void AppendActions(std::string &out, const SubmissionSpec &env) {
	// Invariant: HOLD action requires a date; date is set by adding a separate
	// HOLD entry alongside the user-chosen action (typically ADD). Setting both
	// `action=HOLD` and `hold_until_date` would produce a double-HOLD, so
	// reject it here.
	if (env.action == ENAAction::HOLD && env.hold_until_date.empty()) {
		throw std::runtime_error("ENA envelope: HOLD action requires hold_until_date");
	}
	if (env.action == ENAAction::HOLD && !env.hold_until_date.empty()) {
		throw std::runtime_error(
		    "ENA envelope: with hold_until_date, use action=ADD; the HOLD entry is added automatically");
	}
	out.append("\"actions\":[");
	out.append("{\"type\":");
	AppendJsonString(out, ActionName(env.action));
	out.push_back('}');
	if (!env.hold_until_date.empty()) {
		out.append(",{\"type\":\"HOLD\",\"holdUntilDate\":");
		AppendJsonString(out, env.hold_until_date);
		out.push_back('}');
	}
	out.push_back(']');
}

void AppendProject(std::string &out, const ProjectSpec &p) {
	if (p.alias.empty()) {
		throw std::runtime_error("ENA envelope: project alias must be non-empty");
	}
	out.push_back('{');
	out.append("\"alias\":");
	AppendJsonString(out, p.alias);
	out.append(",\"title\":");
	AppendJsonString(out, p.title);
	if (!p.description.empty()) {
		out.append(",\"description\":");
		AppendJsonString(out, p.description);
	}
	out.append(p.is_umbrella ? ",\"umbrellaProject\":{}" : ",\"sequencingProject\":{}");
	out.push_back('}');
}

void AppendSample(std::string &out, const SampleSpec &s) {
	if (s.alias.empty()) {
		throw std::runtime_error("ENA envelope: sample alias must be non-empty");
	}
	if (s.taxon_id <= 0) {
		throw std::runtime_error("ENA envelope: sample.taxon_id must be > 0 (got " + std::to_string(s.taxon_id) +
		                         " for alias '" + s.alias + "')");
	}
	out.push_back('{');
	out.append("\"alias\":");
	AppendJsonString(out, s.alias);
	if (!s.title.empty()) {
		out.append(",\"title\":");
		AppendJsonString(out, s.title);
	}
	if (!s.description.empty()) {
		out.append(",\"description\":");
		AppendJsonString(out, s.description);
	}
	out.append(",\"organism\":{\"taxonId\":");
	AppendJsonString(out, std::to_string(s.taxon_id));
	if (!s.scientific_name.empty()) {
		out.append(",\"scientificName\":");
		AppendJsonString(out, s.scientific_name);
	}
	out.push_back('}');

	const bool any_attrs = !s.checklist.empty() || !s.attributes.empty();
	if (any_attrs) {
		out.append(",\"attributes\":[");
		bool first = true;
		if (!s.checklist.empty()) {
			out.append("{\"tag\":\"ENA-CHECKLIST\",\"value\":");
			AppendJsonString(out, s.checklist);
			out.push_back('}');
			first = false;
		}
		for (const auto &kv : s.attributes) {
			if (!first) {
				out.push_back(',');
			}
			out.append("{\"tag\":");
			AppendJsonString(out, kv.first);
			out.append(",\"value\":");
			AppendJsonString(out, kv.second);
			out.push_back('}');
			first = false;
		}
		out.push_back(']');
	}
	out.push_back('}');
}

template <typename T, typename Appender>
void AppendArray(std::string &out, bool &needs_comma, const char *key, const std::vector<T> &items, Appender appender) {
	if (items.empty()) {
		return;
	}
	if (needs_comma) {
		out.push_back(',');
	}
	out.push_back('"');
	out.append(key);
	out.append("\":[");
	bool first = true;
	for (const auto &item : items) {
		if (!first) {
			out.push_back(',');
		}
		appender(out, item);
		first = false;
	}
	out.push_back(']');
	needs_comma = true;
}

} // namespace

std::string BuildEnvelopeJSON(const SubmissionSpec &env) {
	std::string out;
	out.reserve(256);
	out.push_back('{');
	out.append("\"submission\":{");
	AppendActions(out, env);
	out.push_back('}');

	bool needs_comma = true;
	AppendArray(out, needs_comma, "projects", env.projects, AppendProject);
	AppendArray(out, needs_comma, "samples", env.samples, AppendSample);

	out.push_back('}');
	return out;
}

} // namespace miint
