#pragma once

#include "sam_style_tags.hpp"

#include <string>
#include <string_view>
#include <vector>

namespace miint {

// One GFA record, split the way the GFA 1 specification defines its record type.
struct GfaRecord {
	std::string record_type;
	// The record type's mandatory positional fields, verbatim. For a record type the
	// specification does not define, every field after the type.
	std::vector<std::string> fields;
	// False for a record type the specification does not define: its fields are not split
	// into positional fields and tags, because nothing says where one ends.
	bool has_tags = false;
	std::vector<SamStyleTag> tags;
};

// Parses GFA 1.x lines (https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md). Holds
// scratch space, so one parser is reused across a file's lines.
class GfaParser {
public:
	// Parse one line into `record`. Returns false for a comment line, which carries no
	// record. Throws std::runtime_error describing the defect; the caller adds the file and
	// line number.
	bool Parse(std::string_view line, GfaRecord &record);

private:
	std::vector<std::string_view> fields_;
};

} // namespace miint
