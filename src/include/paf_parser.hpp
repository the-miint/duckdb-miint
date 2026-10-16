#pragma once

#include "sam_style_tags.hpp"

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace miint {

// One PAF record, with coordinates in miint's convention (docs/reading.md, "Coordinate
// conventions"): the read axis stays 0-based half-open as PAF writes it, and the
// reference axis becomes 1-based half-open.
struct PafRecord {
	std::string read_id;
	int64_t query_length = 0;
	int64_t query_start = 0; // read axis, 0-based half-open, on the read as sequenced
	int64_t query_stop = 0;
	std::string strand;    // "+" or "-"; "*" for an unmapped record
	std::string reference; // "*" for an unmapped record
	int64_t reference_length = 0;
	int64_t position = 0;      // reference axis, 1-based half-open; 0 for an unmapped record
	int64_t stop_position = 0; // 0 for an unmapped record
	int64_t residue_matches = 0;
	int64_t block_length = 0;
	uint8_t mapq = 0;
	std::vector<SamStyleTag> tags;
};

// Parses PAF lines (https://github.com/lh3/miniasm/blob/master/PAF.md). Holds scratch
// space, so one parser is reused across a file's lines.
class PafParser {
public:
	// Parse one line into `record`. Throws std::runtime_error describing the defect; the
	// caller adds the file and line number.
	void Parse(std::string_view line, PafRecord &record);

private:
	std::vector<std::string_view> fields_;
};

} // namespace miint
