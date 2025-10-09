#include <cstddef>
#include <cstdint>
#include <span>

namespace miint {

std::size_t find_low_quality_window(const std::span<uint8_t> &data, uint8_t min_quality, std::size_t window_length) {
	const std::size_t n = data.size();
	if (window_length == 0 || window_length > n) {
		return n + 1;
	}

	std::size_t run = 0;
	for (std::size_t i = 0; i < n; ++i) {
		if (data[i] < min_quality) {
			if (++run >= window_length) {
				return i - window_length + 1;
			}
		} else {
			run = 0;
		}
	}
	return n + 1;
}
} // namespace miint
