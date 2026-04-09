#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

class CoverageDepthCalculator {
public:
	CoverageDepthCalculator(int64_t reference_length, bool include_deletions);

	void AddRead(int64_t position, int64_t stop_position, const std::string &cigar);
	void Combine(const CoverageDepthCalculator &other);
	const std::vector<uint32_t> &GetDepths() const;
	bool Empty() const;

private:
	void IncrementRange(int64_t start, int64_t end);
	std::vector<uint32_t> depths;
	int64_t reference_length;
	bool include_deletions;
	bool has_reads;
};

} // namespace miint
