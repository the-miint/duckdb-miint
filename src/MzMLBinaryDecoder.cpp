#include "MzMLBinaryDecoder.hpp"
#include <cstring>
#include <limits>
#include <stdexcept>
#include <zlib.h>

// mzML binary arrays are defined as little-endian IEEE 754 (mzML spec section 5.9).
// memcpy-based reinterpretation only produces correct results on LE platforms.
static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__, "mzML binary decoding requires a little-endian platform");

namespace miint {

// Base64 decoding lookup table: maps ASCII char value to 6-bit value
// 255 = invalid, 254 = padding '='
// NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays)
static const uint8_t BASE64_DECODE_TABLE[256] = {
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 0-15
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 16-31
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 62,  255, 255, 255, 63,  // 32-47 ('+' = 62, '/' = 63)
    52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  255, 255, 255, 254, 255, 255, // 48-63 ('0'-'9', '=' = 254)
    255, 0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  // 64-79 ('A'-'O')
    15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  255, 255, 255, 255, 255, // 80-95 ('P'-'Z')
    255, 26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  // 96-111 ('a'-'o')
    41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  255, 255, 255, 255, 255, // 112-127 ('p'-'z')
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 128-143
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 144-159
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 160-175
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 176-191
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 192-207
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 208-223
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 224-239
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 240-255
};
// NOLINTEND(cppcoreguidelines-avoid-c-arrays)

std::vector<uint8_t> MzMLBinaryDecoder::base64_decode(const std::string &input) {
	if (input.empty()) {
		return {};
	}

	// Strip whitespace only if present (common case: no whitespace in mzML binary data)
	bool has_whitespace = false;
	for (char c : input) {
		if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
			has_whitespace = true;
			break;
		}
	}

	const std::string *data_ptr = &input;
	std::string cleaned;
	if (has_whitespace) {
		cleaned.reserve(input.size());
		for (char c : input) {
			if (c != ' ' && c != '\t' && c != '\n' && c != '\r') {
				cleaned.push_back(c);
			}
		}
		if (cleaned.empty()) {
			return {};
		}
		data_ptr = &cleaned;
	}

	const std::string &data = *data_ptr;

	// Validate length is multiple of 4
	if (data.size() % 4 != 0) {
		throw std::invalid_argument("base64_decode: input length is not a multiple of 4");
	}

	// Validate: padding only at end, and proper structure
	bool seen_padding = false;
	for (size_t i = 0; i < data.size(); i++) {
		uint8_t v = BASE64_DECODE_TABLE[static_cast<uint8_t>(data[i])];
		if (v == 254) {
			// Padding '=' — only valid at positions 2 or 3 of the last group
			size_t pos_in_group = i % 4;
			bool is_last_group = (i / 4 == (data.size() / 4) - 1);
			if (!is_last_group || pos_in_group < 2) {
				throw std::invalid_argument("base64_decode: padding '=' at invalid position");
			}
			seen_padding = true;
		} else if (v == 255) {
			throw std::invalid_argument("base64_decode: invalid character in input");
		} else if (seen_padding) {
			// Data character after padding
			throw std::invalid_argument("base64_decode: data after padding");
		}
	}

	// Calculate output size accounting for padding
	size_t output_size = (data.size() / 4) * 3;
	if (data.size() >= 1 && data[data.size() - 1] == '=') {
		output_size--;
	}
	if (data.size() >= 2 && data[data.size() - 2] == '=') {
		output_size--;
	}

	std::vector<uint8_t> result;
	result.reserve(output_size);

	for (size_t i = 0; i < data.size(); i += 4) {
		uint8_t a = BASE64_DECODE_TABLE[static_cast<uint8_t>(data[i])];
		uint8_t b = BASE64_DECODE_TABLE[static_cast<uint8_t>(data[i + 1])];
		uint8_t c = BASE64_DECODE_TABLE[static_cast<uint8_t>(data[i + 2])];
		uint8_t d = BASE64_DECODE_TABLE[static_cast<uint8_t>(data[i + 3])];

		// Treat padding as 0 for the triple calculation
		if (c == 254) {
			c = 0;
		}
		if (d == 254) {
			d = 0;
		}

		uint32_t triple = (static_cast<uint32_t>(a) << 18) | (static_cast<uint32_t>(b) << 12) |
		                  (static_cast<uint32_t>(c) << 6) | static_cast<uint32_t>(d);

		result.push_back(static_cast<uint8_t>((triple >> 16) & 0xFF));
		if (result.size() < output_size) {
			result.push_back(static_cast<uint8_t>((triple >> 8) & 0xFF));
		}
		if (result.size() < output_size) {
			result.push_back(static_cast<uint8_t>(triple & 0xFF));
		}
	}

	return result;
}

std::vector<uint8_t> MzMLBinaryDecoder::zlib_inflate(const std::vector<uint8_t> &compressed) {
	if (compressed.empty()) {
		return {};
	}

	z_stream strm = {};
	strm.next_in = const_cast<Bytef *>(compressed.data());
	if (compressed.size() > static_cast<size_t>(std::numeric_limits<uInt>::max())) {
		throw std::runtime_error("zlib_inflate: input too large for zlib (> 4GB)");
	}
	strm.avail_in = static_cast<uInt>(compressed.size());

	int ret = inflateInit(&strm);
	if (ret != Z_OK) {
		throw std::runtime_error("zlib_inflate: inflateInit failed");
	}

	std::vector<uint8_t> result;
	uint8_t buffer[16384];

	do {
		strm.next_out = buffer;
		strm.avail_out = sizeof(buffer);
		ret = inflate(&strm, Z_NO_FLUSH);

		if (ret != Z_OK && ret != Z_STREAM_END) {
			inflateEnd(&strm);
			throw std::runtime_error("zlib_inflate: corrupt or invalid compressed data");
		}

		size_t have = sizeof(buffer) - strm.avail_out;
		result.insert(result.end(), buffer, buffer + have);
	} while (ret != Z_STREAM_END);

	inflateEnd(&strm);
	return result;
}

std::vector<double> MzMLBinaryDecoder::to_doubles_64(const std::vector<uint8_t> &bytes) {
	if (bytes.size() % 8 != 0) {
		throw std::invalid_argument("to_doubles_64: byte count is not a multiple of 8");
	}

	size_t count = bytes.size() / 8;
	std::vector<double> result(count);

	for (size_t i = 0; i < count; i++) {
		// memcpy avoids strict aliasing violations
		std::memcpy(&result[i], bytes.data() + i * 8, 8);
	}

	return result;
}

std::vector<double> MzMLBinaryDecoder::to_doubles_32(const std::vector<uint8_t> &bytes) {
	if (bytes.size() % 4 != 0) {
		throw std::invalid_argument("to_doubles_32: byte count is not a multiple of 4");
	}

	size_t count = bytes.size() / 4;
	std::vector<double> result(count);

	for (size_t i = 0; i < count; i++) {
		float f;
		std::memcpy(&f, bytes.data() + i * 4, 4);
		result[i] = static_cast<double>(f);
	}

	return result;
}

std::vector<double> MzMLBinaryDecoder::decode(const std::string &base64_text, bool is_compressed, bool is_64bit) {
	auto raw_bytes = base64_decode(base64_text);

	if (is_compressed) {
		raw_bytes = zlib_inflate(raw_bytes);
	}

	if (is_64bit) {
		return to_doubles_64(raw_bytes);
	} else {
		return to_doubles_32(raw_bytes);
	}
}

std::vector<uint8_t> MzMLBinaryDecoder::zlib_compress_for_test(const std::vector<uint8_t> &data) {
	if (data.empty()) {
		return {};
	}

	z_stream strm = {};
	int ret = deflateInit(&strm, Z_DEFAULT_COMPRESSION);
	if (ret != Z_OK) {
		throw std::runtime_error("zlib_compress_for_test: deflateInit failed");
	}

	strm.next_in = const_cast<Bytef *>(data.data());
	if (data.size() > static_cast<size_t>(std::numeric_limits<uInt>::max())) {
		throw std::runtime_error("zlib_compress_for_test: input too large for zlib (> 4GB)");
	}
	strm.avail_in = static_cast<uInt>(data.size());

	std::vector<uint8_t> result;
	uint8_t buffer[16384];

	do {
		strm.next_out = buffer;
		strm.avail_out = sizeof(buffer);
		ret = deflate(&strm, Z_FINISH);

		if (ret < 0) {
			deflateEnd(&strm);
			throw std::runtime_error("zlib_compress_for_test: deflate failed");
		}

		size_t have = sizeof(buffer) - strm.avail_out;
		result.insert(result.end(), buffer, buffer + have);
	} while (ret != Z_STREAM_END);

	deflateEnd(&strm);
	return result;
}

} // namespace miint
