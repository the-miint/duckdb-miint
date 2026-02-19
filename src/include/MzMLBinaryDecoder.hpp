#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

class MzMLBinaryDecoder {
public:
	[[nodiscard]] static std::vector<uint8_t> base64_decode(const std::string &input);
	[[nodiscard]] static std::vector<uint8_t> zlib_inflate(const std::vector<uint8_t> &compressed);
	[[nodiscard]] static std::vector<double> to_doubles_64(const std::vector<uint8_t> &bytes);
	[[nodiscard]] static std::vector<double> to_doubles_32(const std::vector<uint8_t> &bytes);
	[[nodiscard]] static std::vector<double> decode(const std::string &base64_text, bool is_compressed, bool is_64bit);

	// Test helper: compress data with zlib for round-trip testing
	[[nodiscard]] static std::vector<uint8_t> zlib_compress_for_test(const std::vector<uint8_t> &data);
};

} // namespace miint
