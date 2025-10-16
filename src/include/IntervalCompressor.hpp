#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

namespace miint {

class IntervalCompressor {
public:
	std::vector<int64_t> starts;
	std::vector<int64_t> stops;

	IntervalCompressor();

	void Add(int64_t start, int64_t stop);
	void Compress();
	bool Empty() const;
	size_t Size() const;
};

} // namespace miint
