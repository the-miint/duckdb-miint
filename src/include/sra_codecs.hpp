#pragma once

#include "sra_column.hpp"

#include <cstdint>
#include <vector>

namespace miint {
namespace sra {

// Decoding of a physical column's blobs: the encodings public runs use on the columns the reader
// reads (raw, zip_encoding, izip_encoding). Output is the column's stored elements in native
// (little-endian) byte order.

// zip_encoding: raw deflate. Version 1 inflates to the frame's osize bytes; version 2 also pops
// the count of meaningful bits in the last byte. Returns the inflated bytes; out_bits receives
// their length in bits.
std::vector<uint8_t> Unzip(const HeaderFrame &frame, const uint8_t *src, size_t n, uint64_t &out_bits);

// izip_encoding (VDB's irzip, header versions 1-3): byte planes, each a raw deflate stream, then
// an offset/slope, delta, or two-series delta reconstruction. elem_count values of the given
// width are produced. Version 0 (legacy izip) raises SraUnsupportedError: no public run surveyed
// used it.
std::vector<uint8_t> Iunzip(const HeaderFrame &frame, const uint8_t *src, size_t n, ElementType type,
                            uint64_t elem_count);

// A blob decoded to its column's elements. `data` holds `elem_count` elements of the column's
// width (bit-packed for widths under 8); the page map says where each row's elements are.
struct DecodedBlob {
	std::vector<uint8_t> data;
	uint64_t elem_count = 0;
	PageMap page_map;
};

// Decodes a blob of a column stored with `encoding` (empty when unencoded). Throws
// SraUnsupportedError for an encoding outside the supported set.
DecodedBlob DecodeBlob(VdbBlob blob, const std::string &encoding, ElementType type);

} // namespace sra
} // namespace miint
