#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

struct SliceInput {
	std::string cigar;
	int64_t position;             // 1-based
	int64_t stop_position;        // 1-based, exclusive
	std::string sequence;         // may be empty (not available)
	std::vector<uint8_t> quality; // may be empty (not available)
};

struct SliceResult {
	bool overlaps; // false = exclude from output
	std::string cigar;
	int64_t position;      // adjusted, 1-based
	int64_t stop_position; // adjusted, 1-based exclusive
	std::string sequence;
	std::vector<uint8_t> quality;
};

class AlignmentSlicer {
public:
	// region_start and region_stop are 1-based, half-open [start, stop)
	AlignmentSlicer(int64_t region_start, int64_t region_stop, bool include_deletions);

	SliceResult Slice(const SliceInput &input) const;

private:
	int64_t region_start;
	int64_t region_stop;
	bool include_deletions;
};

} // namespace miint
