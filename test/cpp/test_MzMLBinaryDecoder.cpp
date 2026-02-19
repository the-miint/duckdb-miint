#include <MzMLBinaryDecoder.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstring>

using miint::MzMLBinaryDecoder;

// ===== base64_decode =====

TEST_CASE("base64_decode: empty string", "[mzml][binary]") {
	auto result = MzMLBinaryDecoder::base64_decode("");
	CHECK(result.empty());
}

TEST_CASE("base64_decode: whitespace only", "[mzml][binary]") {
	auto result = MzMLBinaryDecoder::base64_decode("  \n\t\r  ");
	CHECK(result.empty());
}

TEST_CASE("base64_decode: known value Hello", "[mzml][binary]") {
	// "SGVsbG8=" -> "Hello"
	auto result = MzMLBinaryDecoder::base64_decode("SGVsbG8=");
	REQUIRE(result.size() == 5);
	CHECK(result[0] == 'H');
	CHECK(result[1] == 'e');
	CHECK(result[2] == 'l');
	CHECK(result[3] == 'l');
	CHECK(result[4] == 'o');
}

TEST_CASE("base64_decode: handles embedded whitespace", "[mzml][binary]") {
	// Same as "SGVsbG8=" but with whitespace
	auto result = MzMLBinaryDecoder::base64_decode("SGVs\nbG8=");
	REQUIRE(result.size() == 5);
	std::string str(result.begin(), result.end());
	CHECK(str == "Hello");
}

TEST_CASE("base64_decode: handles + and / characters", "[mzml][binary]") {
	// "++++", decodes to bytes 0xfb 0xef 0xbe
	auto result = MzMLBinaryDecoder::base64_decode("++++");
	REQUIRE(result.size() == 3);
	CHECK(result[0] == 0xfb);
	CHECK(result[1] == 0xef);
	CHECK(result[2] == 0xbe);
}

TEST_CASE("base64_decode: invalid characters", "[mzml][binary]") {
	CHECK_THROWS_AS(MzMLBinaryDecoder::base64_decode("!!!@"), std::invalid_argument);
}

TEST_CASE("base64_decode: non-multiple-of-4 length", "[mzml][binary]") {
	CHECK_THROWS_AS(MzMLBinaryDecoder::base64_decode("SGVSB"), std::invalid_argument);
	CHECK_THROWS_AS(MzMLBinaryDecoder::base64_decode("AB"), std::invalid_argument);
}

// ===== zlib_inflate =====

TEST_CASE("zlib_inflate: round-trip compress/decompress", "[mzml][binary]") {
	std::vector<uint8_t> original = {1, 2, 3, 4, 5, 6, 7, 8};
	auto compressed = MzMLBinaryDecoder::zlib_compress_for_test(original);
	auto decompressed = MzMLBinaryDecoder::zlib_inflate(compressed);
	CHECK(decompressed == original);
}

TEST_CASE("zlib_inflate: corrupt data rejection", "[mzml][binary]") {
	std::vector<uint8_t> garbage = {0xFF, 0xFE, 0xFD, 0xFC};
	CHECK_THROWS_WITH(MzMLBinaryDecoder::zlib_inflate(garbage), Catch::Matchers::ContainsSubstring("corrupt"));
}

TEST_CASE("zlib_inflate: empty input", "[mzml][binary]") {
	auto result = MzMLBinaryDecoder::zlib_inflate({});
	CHECK(result.empty());
}

// ===== to_doubles_64 =====

TEST_CASE("to_doubles_64: known IEEE 754 LE value 100.0", "[mzml][binary]") {
	// 100.0 = 0x4059000000000000 in LE: 00 00 00 00 00 00 59 40
	std::vector<uint8_t> bytes = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x59, 0x40};
	auto result = MzMLBinaryDecoder::to_doubles_64(bytes);
	REQUIRE(result.size() == 1);
	CHECK(result[0] == 100.0);
}

TEST_CASE("to_doubles_64: multiple values", "[mzml][binary]") {
	// 100.0, 200.0, 300.0
	double vals[] = {100.0, 200.0, 300.0};
	std::vector<uint8_t> bytes(24);
	std::memcpy(bytes.data(), vals, 24);

	auto result = MzMLBinaryDecoder::to_doubles_64(bytes);
	REQUIRE(result.size() == 3);
	CHECK(result[0] == 100.0);
	CHECK(result[1] == 200.0);
	CHECK(result[2] == 300.0);
}

TEST_CASE("to_doubles_64: non-aligned rejection", "[mzml][binary]") {
	std::vector<uint8_t> bytes = {0x00, 0x00, 0x00}; // 3 bytes, not multiple of 8
	CHECK_THROWS_AS(MzMLBinaryDecoder::to_doubles_64(bytes), std::invalid_argument);
}

// ===== to_doubles_32 =====

TEST_CASE("to_doubles_32: 32-bit float to double upcast", "[mzml][binary]") {
	// 100.0f = 0x42C80000 in LE: 00 00 C8 42
	float val = 100.0f;
	std::vector<uint8_t> bytes(4);
	std::memcpy(bytes.data(), &val, 4);

	auto result = MzMLBinaryDecoder::to_doubles_32(bytes);
	REQUIRE(result.size() == 1);
	CHECK_THAT(result[0], Catch::Matchers::WithinRel(100.0, 1e-6));
}

TEST_CASE("to_doubles_32: multiple values", "[mzml][binary]") {
	float vals[] = {100.0f, 200.0f, 300.0f};
	std::vector<uint8_t> bytes(12);
	std::memcpy(bytes.data(), vals, 12);

	auto result = MzMLBinaryDecoder::to_doubles_32(bytes);
	REQUIRE(result.size() == 3);
	CHECK_THAT(result[0], Catch::Matchers::WithinRel(100.0, 1e-6));
	CHECK_THAT(result[1], Catch::Matchers::WithinRel(200.0, 1e-6));
	CHECK_THAT(result[2], Catch::Matchers::WithinRel(300.0, 1e-6));
}

TEST_CASE("to_doubles_32: non-aligned rejection", "[mzml][binary]") {
	std::vector<uint8_t> bytes = {0x00, 0x00, 0x00}; // 3 bytes, not multiple of 4
	CHECK_THROWS_AS(MzMLBinaryDecoder::to_doubles_32(bytes), std::invalid_argument);
}

// ===== decode (full pipeline) =====

TEST_CASE("decode: uncompressed 64-bit", "[mzml][binary]") {
	// base64 of 3 doubles (100.0, 200.0, 300.0) LE 64-bit
	auto result = MzMLBinaryDecoder::decode("AAAAAAAAWUAAAAAAAABpQAAAAAAAwHJA", false, true);
	REQUIRE(result.size() == 3);
	CHECK(result[0] == 100.0);
	CHECK(result[1] == 200.0);
	CHECK(result[2] == 300.0);
}

TEST_CASE("decode: compressed 64-bit", "[mzml][binary]") {
	// zlib-compressed base64 of 3 doubles (100.0, 200.0, 300.0)
	auto result = MzMLBinaryDecoder::decode("eJxjYACBSAcwxZAJoQ8UOQAAFFgCtQ==", true, true);
	REQUIRE(result.size() == 3);
	CHECK(result[0] == 100.0);
	CHECK(result[1] == 200.0);
	CHECK(result[2] == 300.0);
}

TEST_CASE("decode: uncompressed 32-bit", "[mzml][binary]") {
	// base64 of 3 floats (100.0, 200.0, 300.0) LE 32-bit
	auto result = MzMLBinaryDecoder::decode("AADIQgAASEMAAJZD", false, false);
	REQUIRE(result.size() == 3);
	CHECK_THAT(result[0], Catch::Matchers::WithinRel(100.0, 1e-6));
	CHECK_THAT(result[1], Catch::Matchers::WithinRel(200.0, 1e-6));
	CHECK_THAT(result[2], Catch::Matchers::WithinRel(300.0, 1e-6));
}

TEST_CASE("decode: empty base64", "[mzml][binary]") {
	auto result = MzMLBinaryDecoder::decode("", false, true);
	CHECK(result.empty());
}

// ===== base64_decode validation edge cases =====

TEST_CASE("base64_decode: rejects mid-string padding", "[mzml][binary]") {
	// "AA==" is valid (one output byte), but "AA==AAAA" has padding mid-string
	CHECK_THROWS_AS(MzMLBinaryDecoder::base64_decode("AA==AAAA"), std::invalid_argument);
}

TEST_CASE("base64_decode: rejects padding at position 0 or 1 in group", "[mzml][binary]") {
	// Padding only valid at positions 2 and 3 of the last 4-char group
	CHECK_THROWS_AS(MzMLBinaryDecoder::base64_decode("===="), std::invalid_argument);
	CHECK_THROWS_AS(MzMLBinaryDecoder::base64_decode("=AAA"), std::invalid_argument);
}

TEST_CASE("base64_decode: accepts valid double-padding", "[mzml][binary]") {
	// "QQ==" encodes single byte 'A' (0x41)
	auto result = MzMLBinaryDecoder::base64_decode("QQ==");
	REQUIRE(result.size() == 1);
	CHECK(result[0] == 0x41);
}

TEST_CASE("base64_decode: accepts valid single-padding", "[mzml][binary]") {
	// "QUI=" encodes two bytes 'AB' (0x41, 0x42)
	auto result = MzMLBinaryDecoder::base64_decode("QUI=");
	REQUIRE(result.size() == 2);
	CHECK(result[0] == 0x41);
	CHECK(result[1] == 0x42);
}
