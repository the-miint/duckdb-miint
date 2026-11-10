#include <catch2/catch_test_macros.hpp>

#include "sra_codecs.hpp"
#include "sra_test_helpers.hpp"

#include <zlib.h>

#include <cstring>
#include <string>
#include <vector>

using namespace miint::sra;
using miint_test::PackedRows;
using miint_test::ThrowsTypeMentioning;

// Hand-built inputs for the two codecs public runs use. The expected values follow ncbi-vdb's
// unzip.c and irzip.impl.h (the active two-series branch), which the oracles were decoded with;
// the fixture parity test (test_SraSpots.cpp) checks the same code on real data.

namespace {
constexpr int64_t DELTA_POS = 0x7ffffffffffffff0;
constexpr int64_t DELTA_NEG = 0x7ffffffffffffff1;
constexpr int64_t DELTA_BOTH = 0x7ffffffffffffff2;

std::vector<uint8_t> RawDeflate(const std::vector<uint8_t> &in) {
	z_stream zs {};
	REQUIRE(deflateInit2(&zs, Z_BEST_SPEED, Z_DEFLATED, -15, 9, Z_DEFAULT_STRATEGY) == Z_OK);
	std::vector<uint8_t> out(in.size() + 64);
	zs.next_in = const_cast<Bytef *>(in.data());
	zs.avail_in = uInt(in.size());
	zs.next_out = out.data();
	zs.avail_out = uInt(out.size());
	REQUIRE(deflate(&zs, Z_FINISH) == Z_STREAM_END);
	out.resize(zs.total_out);
	deflateEnd(&zs);
	return out;
}

// An irzip payload: one deflate stream per byte plane, low plane first.
std::vector<uint8_t> Planes(const std::vector<std::vector<uint8_t>> &planes) {
	std::vector<uint8_t> out;
	for (const auto &p : planes) {
		auto z = RawDeflate(p);
		out.insert(out.end(), z.begin(), z.end());
	}
	return out;
}

HeaderFrame IzipFrame(uint8_t version, uint8_t planes, std::vector<int64_t> args) {
	HeaderFrame f;
	f.version = version;
	f.ops = {planes};
	f.args = std::move(args);
	return f;
}

template <typename T>
std::vector<T> As(const std::vector<uint8_t> &bytes) {
	std::vector<T> v(bytes.size() / sizeof(T));
	std::memcpy(v.data(), bytes.data(), v.size() * sizeof(T));
	return v;
}

const ElementType U32 {32, false};
const ElementType I32 {32, true};
} // namespace

TEST_CASE("Unzip version 1 inflates to the size the frame records", "[sra]") {
	std::vector<uint8_t> raw {'h', 'e', 'l', 'l', 'o'};
	auto z = RawDeflate(raw);
	HeaderFrame f;
	f.version = 1;
	f.osize = 5;
	uint64_t bits = 0;
	CHECK(Unzip(f, z.data(), z.size(), bits) == raw);
	CHECK(bits == 40);

	// A frame that disagrees with the stream is damage, not a shorter blob.
	f.osize = 6;
	CHECK_THROWS_AS(Unzip(f, z.data(), z.size(), bits), SraFormatError);
}

TEST_CASE("Unzip version 2 drops the unused bits of the last byte", "[sra]") {
	auto z = RawDeflate({0xff, 0xc0});
	HeaderFrame f;
	f.version = 2;
	f.osize = 2;
	f.args = {2}; // two meaningful bits in the last byte
	uint64_t bits = 0;
	Unzip(f, z.data(), z.size(), bits);
	CHECK(bits == 10);
	f.args = {0}; // 0 means the last byte is whole
	Unzip(f, z.data(), z.size(), bits);
	CHECK(bits == 16);
}

TEST_CASE("Unzip refuses a frame version it does not know", "[sra]") {
	auto z = RawDeflate({1});
	HeaderFrame f;
	f.version = 3;
	f.osize = 1;
	uint64_t bits;
	CHECK_THROWS_AS(Unzip(f, z.data(), z.size(), bits), SraUnsupportedError);
}

TEST_CASE("Iunzip reassembles byte planes", "[sra]") {
	// 0x0102 = 258: plane 0 holds 0x02, plane 1 holds 0x01; planes absent from the mask are 0.
	auto src = Planes({{0x02, 0x05}, {0x01, 0x00}});
	auto out = As<uint32_t>(Iunzip(IzipFrame(2, 0b11, {0, 0}), src.data(), src.size(), U32, 2));
	CHECK(out == std::vector<uint32_t> {258, 5});
	// Only plane 2 (bits 16..23) present.
	src = Planes({{0x01}});
	out = As<uint32_t>(Iunzip(IzipFrame(2, 0b100, {0, 0}), src.data(), src.size(), U32, 1));
	CHECK(out == std::vector<uint32_t> {0x10000});
}

TEST_CASE("Iunzip applies an offset and a slope", "[sra]") {
	auto src = Planes({{0, 1, 2}});
	CHECK(As<uint32_t>(Iunzip(IzipFrame(2, 1, {10, 0}), src.data(), src.size(), U32, 3)) ==
	      std::vector<uint32_t> {10, 11, 12});
	// Linear: Y[i] = stored + min + i * slope.
	CHECK(As<uint32_t>(Iunzip(IzipFrame(2, 1, {5, 3}), src.data(), src.size(), U32, 3)) ==
	      std::vector<uint32_t> {5, 9, 13});
	// Version 1 carries no slope: an offset only.
	CHECK(As<uint32_t>(Iunzip(IzipFrame(1, 1, {10}), src.data(), src.size(), U32, 3)) ==
	      std::vector<uint32_t> {10, 11, 12});
	// A negative offset on a signed column.
	CHECK(As<int32_t>(Iunzip(IzipFrame(2, 1, {-10, 0}), src.data(), src.size(), I32, 3)) ==
	      std::vector<int32_t> {-10, -9, -8});
}

TEST_CASE("Iunzip undoes the three delta codings", "[sra]") {
	SECTION("increasing") {
		auto src = Planes({{0, 5, 1}});
		CHECK(As<uint32_t>(Iunzip(IzipFrame(2, 1, {100, DELTA_POS}), src.data(), src.size(), U32, 3)) ==
		      std::vector<uint32_t> {100, 105, 106});
	}
	SECTION("decreasing") {
		auto src = Planes({{0, 5, 1}});
		CHECK(As<uint32_t>(Iunzip(IzipFrame(2, 1, {100, DELTA_NEG}), src.data(), src.size(), U32, 3)) ==
		      std::vector<uint32_t> {100, 95, 94});
	}
	SECTION("either way: the low bit is the sign") {
		auto src = Planes({{0, (5 << 1) | 0, (3 << 1) | 1}});
		CHECK(As<uint32_t>(Iunzip(IzipFrame(2, 1, {100, DELTA_BOTH}), src.data(), src.size(), U32, 3)) ==
		      std::vector<uint32_t> {100, 105, 102});
	}
}

TEST_CASE("Iunzip separates two interleaved series", "[sra]") {
	// Coordinates of two tile rows interleave; the low bit says which series a value continues.
	// Main series from 1000, secondary from 50, both increasing.
	auto src = Planes({{0, (2 << 1) | 1, 3 << 1, (1 << 1) | 1}});
	auto out = As<int32_t>(Iunzip(IzipFrame(3, 1, {1000, DELTA_POS, 50, DELTA_POS}), src.data(), src.size(), I32, 4));
	CHECK(out == std::vector<int32_t> {1000, 52, 1003, 53});
}

TEST_CASE("Iunzip keeps the reference decoder's seeding of a leading secondary value", "[sra]") {
	// When the first value belongs to the secondary series, ncbi-vdb seeds the main series from
	// it and reseeds the secondary series at its next value. fasterq-dump prints what that code
	// produces, so the reader must too.
	auto src = Planes({{(4 << 1) | 1, 2 << 1, (1 << 1) | 1}});
	auto out = As<int32_t>(Iunzip(IzipFrame(3, 1, {1000, DELTA_POS, 50, DELTA_POS}), src.data(), src.size(), I32, 3));
	CHECK(out == std::vector<int32_t> {54, 56, 51});
}

TEST_CASE("Iunzip refuses legacy izip and reports damaged planes", "[sra]") {
	auto src = Planes({{1, 2}});
	CHECK(ThrowsTypeMentioning<SraUnsupportedError>(
	    [&] { Iunzip(IzipFrame(0, 1, {}), src.data(), src.size(), U32, 2); }, "izip"));
	std::vector<uint8_t> junk(8, 0xff);
	CHECK_THROWS_AS(Iunzip(IzipFrame(2, 1, {0, 0}), junk.data(), junk.size(), U32, 2), SraFormatError);
	// A plane that inflates to fewer values than the blob holds.
	CHECK_THROWS_AS(Iunzip(IzipFrame(2, 1, {0, 0}), src.data(), src.size(), U32, 3), SraFormatError);
}

TEST_CASE("DecodeBlob passes through a stage the writer skipped", "[sra]") {
	// Flag bit 0: compression did not pay, so the payload is the elements themselves. The survey
	// found such frames (with version 0) on zip-encoded columns.
	VdbBlob blob;
	blob.bytes = {7, 8, 9};
	blob.payload_offset = 0;
	blob.payload_bits = 24;
	HeaderFrame f;
	f.flags = 1;
	f.version = 0;
	f.osize = 3;
	blob.frames = {f};
	blob.page_map = PackedRows({3});
	DecodedBlob d = DecodeBlob(blob, "<INSDC:4na:bin>zip_encoding#1", ElementType {8, false});
	CHECK(d.data == std::vector<uint8_t> {7, 8, 9});
	CHECK(d.elem_count == 3);
}

TEST_CASE("DecodeBlob byte-swaps an unencoded big-endian payload of wide elements", "[sra]") {
	VdbBlob blob;
	blob.bytes = {0, 0, 1, 2};
	blob.payload_bits = 32;
	blob.big_endian = true;
	blob.page_map = PackedRows({1});
	DecodedBlob d = DecodeBlob(blob, "", U32);
	CHECK(As<uint32_t>(d.data) == std::vector<uint32_t> {0x0102});
}

TEST_CASE("DecodeBlob refuses a page map addressing past the decoded elements", "[sra]") {
	// Offsets come from the file. One near 2^64 plus its row's length wraps to a small number, so
	// each row is checked on its own: a wrapped sum must not pass the bound.
	VdbBlob blob;
	blob.bytes = {7, 8, 9};
	blob.payload_bits = 24;
	blob.page_map = PackedRows({3, 2});
	blob.page_map.runs[1].offset = UINT64_MAX - 1;
	CHECK(ThrowsTypeMentioning<SraFormatError>([&] { DecodeBlob(blob, "", ElementType {8, false}); }, "page map"));
	blob.page_map.runs[1].offset = 2;
	CHECK(ThrowsTypeMentioning<SraFormatError>([&] { DecodeBlob(blob, "", ElementType {8, false}); }, "page map"));
}

TEST_CASE("DecodeBlob refuses izip on elements narrower than a byte", "[sra]") {
	// izip encodes whole integers; a 2na column claiming it is refused, not divided by zero.
	VdbBlob blob;
	blob.bytes = {0};
	blob.payload_bits = 8;
	blob.frames = {IzipFrame(2, 1, {0, 0})};
	blob.frames[0].osize = 1;
	blob.page_map = PackedRows({4});
	CHECK(ThrowsTypeMentioning<SraUnsupportedError>(
	    [&] { DecodeBlob(blob, "izip_encoding#1", ElementTypeOf("INSDC:2na:packed")); }, "2-bit"));
}

TEST_CASE("A frame's decoded size is bounded by its compressed bytes", "[sra]") {
	// Deflate expands at most about 1032 to 1: a frame claiming a terabyte from a few bytes is
	// damage, refused before a buffer is sized from it.
	std::vector<uint8_t> src(16, 0);
	HeaderFrame zip;
	zip.version = 1;
	zip.osize = int64_t(1) << 40;
	uint64_t bits = 0;
	CHECK(ThrowsTypeMentioning<SraFormatError>([&] { Unzip(zip, src.data(), src.size(), bits); }, "inflate"));
	HeaderFrame iz = IzipFrame(2, 1, {0, 0});
	CHECK(ThrowsTypeMentioning<SraFormatError>(
	    [&] {
		    Iunzip(iz, src.data(), src.size(), ElementType {8, false}, uint64_t(1) << 40);
	    },
	    "inflate"));
	// With no byte planes the values come from the header alone (a constant column), so no byte
	// count bounds them: a size cap does.
	HeaderFrame constant = IzipFrame(2, 0, {5, 0});
	CHECK(ThrowsTypeMentioning<SraUnsupportedError>(
	    [&] {
		    Iunzip(constant, nullptr, 0, ElementType {8, false}, uint64_t(1) << 40);
	    },
	    "GiB"));
}

TEST_CASE("A column stored without encoding refuses a blob that carries one", "[sra]") {
	// The column's metadata records no encoding, yet the blob has an encoding header that the
	// writer did not mark as skipped: reading its compressed bytes as elements would return
	// garbage, so it is damage.
	VdbBlob blob;
	blob.bytes = {7, 8, 9};
	blob.payload_bits = 24;
	HeaderFrame f;
	f.version = 1;
	f.osize = 3;
	blob.frames = {f};
	blob.page_map = PackedRows({3});
	CHECK(ThrowsTypeMentioning<SraFormatError>(
	    [&] {
		    DecodeBlob(blob, "", ElementType {8, false});
	    },
	    "without encoding"));
}

TEST_CASE("The deflate bound admits the most compressible stream", "[sra]") {
	// MAX_DEFLATE_RATIO refuses a size no deflate stream can reach. Zeros at level 9 are deflate's
	// best case; a bound below their ratio would refuse real, highly compressed blobs as damage.
	const size_t size = size_t(64) << 20;
	std::vector<uint8_t> zeros(size, 0);
	z_stream zs {};
	REQUIRE(deflateInit2(&zs, 9, Z_DEFLATED, -15, 9, Z_DEFAULT_STRATEGY) == Z_OK);
	std::vector<uint8_t> z(deflateBound(&zs, uLong(size)));
	zs.next_in = zeros.data();
	zs.avail_in = uInt(size);
	zs.next_out = z.data();
	zs.avail_out = uInt(z.size());
	REQUIRE(deflate(&zs, Z_FINISH) == Z_STREAM_END);
	z.resize(zs.total_out);
	deflateEnd(&zs);
	INFO("ratio " << double(size) / double(z.size()));
	CHECK(uint64_t(z.size()) * MAX_DEFLATE_RATIO >= size);
	HeaderFrame f;
	f.version = 1;
	f.osize = int64_t(size);
	uint64_t bits = 0;
	CHECK(Unzip(f, z.data(), z.size(), bits).size() == size);
}

TEST_CASE("Iunzip refuses widths it cannot decode before dividing by them", "[sra]") {
	// DecodeBlob refuses izip on sub-byte columns first, but Iunzip is public: a 2-bit type must be
	// refused, not divided by zero.
	std::vector<uint8_t> src(16, 0);
	CHECK(ThrowsTypeMentioning<SraUnsupportedError>(
	    [&] {
		    Iunzip(IzipFrame(2, 1, {0, 0}), src.data(), src.size(), ElementType {2, false}, 4);
	    },
	    "2-bit"));
}
