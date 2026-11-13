#include "sra_codecs.hpp"

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <type_traits>
#include <zlib.h>

namespace miint {
namespace sra {

namespace {
// zlib's windowBits for a raw deflate stream (no zlib or gzip wrapper) with a 32 KiB window.
constexpr int RAW_DEFLATE_WINDOW_BITS = -15;
// The first buffer InflateRawGrowing tries; it doubles from there.
constexpr uint64_t INFLATE_FIRST_BUFFER = 4096;
} // namespace

std::vector<uint8_t> InflateRawGrowing(const uint8_t *src, size_t n, uint64_t limit, const char *what) {
	if (n > UINT32_MAX) {
		throw SraFormatError(std::string(what) + " is larger than 4 GiB");
	}
	z_stream zs {};
	if (inflateInit2(&zs, RAW_DEFLATE_WINDOW_BITS) != Z_OK) {
		throw SraFormatError(std::string("cannot start inflating ") + what);
	}
	// Ends the stream on every way out, a failed allocation included.
	struct EndStream {
		z_stream &zs;
		~EndStream() {
			inflateEnd(&zs);
		}
	} end_stream {zs};
	zs.next_in = const_cast<Bytef *>(src);
	zs.avail_in = uInt(n);
	std::vector<uint8_t> out;
	int rc = Z_OK;
	while (rc != Z_STREAM_END) {
		if (out.size() >= limit) {
			throw SraFormatError(std::string(what) + " inflates past the " + std::to_string(limit) +
			                     " bytes it can need");
		}
		size_t have = size_t(zs.total_out);
		out.resize(
		    size_t(std::min<uint64_t>(limit, std::max<uint64_t>(INFLATE_FIRST_BUFFER, uint64_t(out.size()) * 2))));
		zs.next_out = out.data() + have;
		zs.avail_out = uInt(std::min<size_t>(out.size() - have, UINT32_MAX));
		rc = inflate(&zs, Z_NO_FLUSH);
		if (rc != Z_OK && rc != Z_STREAM_END && !(rc == Z_BUF_ERROR && zs.avail_out == 0)) {
			throw SraFormatError(std::string(what) + " does not inflate");
		}
		if (rc == Z_OK && zs.avail_in == 0 && zs.avail_out != 0) {
			// The input ran out before the stream ended: truncated.
			throw SraFormatError(std::string(what) + " does not inflate");
		}
	}
	out.resize(size_t(zs.total_out));
	return out;
}

Inflated InflateRaw(const uint8_t *src, size_t n, uint8_t *dst, size_t dsize, const char *what) {
	if (n > UINT32_MAX || dsize > UINT32_MAX) {
		throw SraFormatError(std::string(what) + " is larger than 4 GiB");
	}
	z_stream zs {};
	if (inflateInit2(&zs, RAW_DEFLATE_WINDOW_BITS) != Z_OK) {
		throw SraFormatError(std::string("cannot start inflating ") + what);
	}
	zs.next_in = const_cast<Bytef *>(src);
	zs.avail_in = uInt(n);
	zs.next_out = dst;
	zs.avail_out = uInt(dsize);
	int rc = inflate(&zs, Z_FINISH);
	Inflated r {size_t(zs.total_in), size_t(zs.total_out)};
	inflateEnd(&zs);
	if (rc != Z_STREAM_END) {
		throw SraFormatError(std::string(what) + " does not inflate");
	}
	return r;
}

namespace {

// Inflates raw deflate into exactly dsize bytes and reports how many input bytes the stream used.
// The stream must end exactly there: a short or overlong stream is damage.
size_t InflateExact(const uint8_t *src, size_t n, uint8_t *dst, size_t dsize, const char *what) {
	Inflated r = InflateRaw(src, n, dst, dsize, what);
	if (r.produced != dsize) {
		throw SraFormatError(std::string(what) + " does not inflate to the " + std::to_string(dsize) +
		                     " bytes its header records");
	}
	return r.used;
}

// irzip's markers in the slope argument for delta-coded series (irzip.impl.h).
constexpr int64_t DELTA_POS = 0x7ffffffffffffff0;
constexpr int64_t DELTA_NEG = 0x7ffffffffffffff1;
constexpr int64_t DELTA_BOTH = 0x7ffffffffffffff2;

// Reconstruction after the byte planes are reassembled: a port of ncbi-vdb's DECODE in
// libs/vxf/irzip.impl.h (public domain; THIRD_PARTY_LICENSES.md), its active two-series branch
// included. Arithmetic is on the unsigned type so it wraps as the C code does.
template <typename U>
void IrzipReconstruct(U *Y, size_t N, const int64_t *min, const int64_t *slope, int series_count) {
	if (series_count == 2) {
		// The low bit of each value says which series it continues (0 main, 1 secondary). The
		// main series is seeded at index 0; the secondary series reseeds while its last index is 0,
		// which matters when the first value is secondary. Both reproduce the reference decoder.
		uint32_t last[2] = {0, 0};
		for (size_t i = 0; i < N; i++) {
			int k = (Y[i] & 1) ? 1 : 0;
			U val = U(Y[i] >> 1);
			if (k == 0) {
				if (i == 0) {
					Y[i] = U(min[0]);
					last[0] = uint32_t(i);
				}
			} else if (last[1] == 0) {
				Y[i] = U(min[1]);
				last[1] = uint32_t(i);
			}
			U base = Y[last[k]];
			if (slope[k] == DELTA_POS) {
				Y[i] = U(base + val);
			} else if (slope[k] == DELTA_NEG) {
				Y[i] = U(base - val);
			} else if ((val & 1) == 0) {
				Y[i] = U(base + U(val >> 1));
			} else {
				Y[i] = U(base - U(val >> 1));
			}
			last[k] = uint32_t(i);
		}
		return;
	}
	if (min[0] == 0 && slope[0] == 0) {
		return;
	}
	if (slope[0] == DELTA_POS || slope[0] == DELTA_NEG || slope[0] == DELTA_BOTH) {
		if (N == 0) {
			return;
		}
		Y[0] = U(min[0]);
		for (size_t i = 1; i < N; i++) {
			if (slope[0] == DELTA_POS) {
				Y[i] = U(Y[i - 1] + Y[i]);
			} else if (slope[0] == DELTA_NEG) {
				Y[i] = U(Y[i - 1] - Y[i]);
			} else { // DELTA_BOTH: the low bit is the sign
				U val = U(Y[i] >> 1);
				Y[i] = (Y[i] & 1) ? U(Y[i - 1] - val) : U(Y[i - 1] + val);
			}
		}
		return;
	}
	if (slope[0] == 0) {
		for (size_t i = 0; i < N; i++) {
			Y[i] = U(Y[i] + U(min[0]));
		}
		return;
	}
	uint64_t m = uint64_t(min[0]);
	for (size_t i = 0; i < N; i++) {
		Y[i] = U(Y[i] + U(m));
		m += uint64_t(slope[0]);
	}
}

// An izip header's arguments: the main series' minimum and slope, then, for two series, the
// secondary's. Its first op byte has a bit per byte plane stored, least significant plane first.
enum IzipArg : size_t { MIN_0 = 0, SLOPE_0 = 1, MIN_1 = 2, SLOPE_1 = 3, TWO_SERIES_ARGS = 4 };

// The plane loop and argument unpacking follow the first half of ncbi-vdb's DECODE in
// libs/vxf/irzip.impl.h and its caller iunzip in libs/vxf/irzip.c (public domain;
// THIRD_PARTY_LICENSES.md).
template <typename U>
std::vector<uint8_t> IrzipDecode(const HeaderFrame &frame, const uint8_t *src, size_t n, uint64_t count) {
	if (frame.ops.empty() || frame.args.empty()) {
		throw SraFormatError("an izip blob header lacks its parameters");
	}
	uint8_t planes = frame.ops[0];
	int64_t min[2] = {frame.args[MIN_0], 0};
	int64_t slope[2] = {frame.args.size() > SLOPE_0 ? frame.args[SLOPE_0] : 0, 0};
	int series_count = 1;
	if (frame.args.size() >= TWO_SERIES_ARGS) {
		min[1] = frame.args[MIN_1];
		slope[1] = frame.args[SLOPE_1];
		series_count = 2;
	}
	std::vector<U> Y(count, 0);
	std::vector<uint8_t> plane(planes != 0 ? count : 0); // only a blob with planes reads any
	size_t used = 0;
	for (unsigned p = 0; p < sizeof(planes) * BITS_PER_BYTE; p++) {
		if (!(planes & (1u << p))) {
			continue;
		}
		if (p >= sizeof(U)) {
			throw SraFormatError("an izip blob has a byte plane wider than its type");
		}
		used += InflateExact(src + used, n - used, plane.data(), count, "an izip byte plane");
		for (size_t i = 0; i < count; i++) {
			Y[i] = U(Y[i] | (U(plane[i]) << (BITS_PER_BYTE * p)));
		}
	}
	IrzipReconstruct(Y.data(), count, min, slope, series_count);
	std::vector<uint8_t> out(count * sizeof(U));
	std::memcpy(out.data(), Y.data(), out.size());
	return out;
}

void ByteSwap(std::vector<uint8_t> &data, uint32_t width) {
	for (size_t i = 0; i + width <= data.size(); i += width) {
		for (uint32_t a = 0, b = width - 1; a < b; a++, b--) {
			std::swap(data[i + a], data[i + b]);
		}
	}
}

enum class Encoding { Raw, Zip, Izip };
// The physical encodings read, as a column's expression names them. zip_encoding may carry
// template arguments after its name.
constexpr char IZIP_ENCODING[] = "izip_encoding#1";
constexpr char ZIP_ENCODING[] = "zip_encoding#1";

Encoding EncodingOf(const std::string &expr) {
	if (expr.empty()) {
		return Encoding::Raw;
	}
	std::string e;
	for (char c : expr) {
		if (c != ' ') {
			e += c;
		}
	}
	auto close = e.find('>');
	std::string name = e[0] == '<' && close != std::string::npos ? e.substr(close + 1) : e;
	if (name == IZIP_ENCODING) {
		return Encoding::Izip;
	}
	if (name.compare(0, sizeof(ZIP_ENCODING) - 1, ZIP_ENCODING) == 0) {
		return Encoding::Zip;
	}
	throw SraUnsupportedError("encoding '" + expr + "'");
}

} // namespace

namespace {
// A zip header's version: 1 decodes to whole bytes, 2 records the significant bits of the last.
constexpr uint8_t ZIP_WHOLE_BYTES = 1;
constexpr uint8_t ZIP_TRAILING_BITS = 2;
// izip header versions: 0 is a legacy layout, 1 to 3 the current one.
constexpr uint8_t IZIP_LEGACY = 0;
constexpr uint8_t IZIP_MAX_VERSION = 3;
// A header frame whose flags have this bit passes its data through unchanged.
constexpr uint8_t HEADER_FLAG_SKIPPED = 1;
constexpr int GIB_SHIFT = 30;
} // namespace

std::vector<uint8_t> Unzip(const HeaderFrame &frame, const uint8_t *src, size_t n, uint64_t &out_bits) {
	if (frame.version != ZIP_WHOLE_BYTES && frame.version != ZIP_TRAILING_BITS) {
		throw SraUnsupportedError("zip blob header version " + std::to_string(frame.version));
	}
	if (frame.osize < 0 || uint64_t(frame.osize) > uint64_t(n) * MAX_DEFLATE_RATIO) {
		throw SraFormatError("a zip blob claims " + std::to_string(frame.osize) + " bytes from " + std::to_string(n) +
		                     ", more than deflate can inflate");
	}
	std::vector<uint8_t> out(size_t(frame.osize));
	InflateExact(src, n, out.data(), out.size(), "a zip blob");
	out_bits = uint64_t(out.size()) * BITS_PER_BYTE;
	if (frame.version == ZIP_TRAILING_BITS) {
		if (frame.args.empty()) {
			throw SraFormatError("a zip blob header lacks its trailing-bit count");
		}
		// The last byte's significant bits, 0 meaning all of them.
		int64_t trailing = frame.args[0];
		if (trailing < 0 || trailing >= BITS_PER_BYTE || (trailing != 0 && out.empty())) {
			throw SraFormatError("a zip blob header records " + std::to_string(trailing) + " trailing bits");
		}
		if (trailing != 0) {
			out_bits -= uint64_t(BITS_PER_BYTE - trailing);
		}
	}
	return out;
}

std::vector<uint8_t> Iunzip(const HeaderFrame &frame, const uint8_t *src, size_t n, ElementType type,
                            uint64_t elem_count) {
	if (frame.version == IZIP_LEGACY) {
		throw SraUnsupportedError("legacy izip (blob header version 0), which no public run surveyed used");
	}
	if (frame.version > IZIP_MAX_VERSION) {
		throw SraUnsupportedError("izip blob header version " + std::to_string(frame.version));
	}
	constexpr uint32_t U8_BITS = sizeof(uint8_t) * BITS_PER_BYTE, U16_BITS = sizeof(uint16_t) * BITS_PER_BYTE,
	                   U32_BITS = sizeof(uint32_t) * BITS_PER_BYTE, U64_BITS = sizeof(uint64_t) * BITS_PER_BYTE;
	if (type.bits != U8_BITS && type.bits != U16_BITS && type.bits != U32_BITS && type.bits != U64_BITS) {
		throw SraUnsupportedError("izip on " + std::to_string(type.bits) + "-bit elements");
	}
	// Each byte plane inflates to elem_count bytes from at most n. A blob with no planes is all
	// offset and slope (a constant column), so only the size cap bounds it.
	const bool has_planes = !frame.ops.empty() && frame.ops[0] != 0;
	if (has_planes && elem_count > uint64_t(n) * MAX_DEFLATE_RATIO) {
		throw SraFormatError("an izip blob claims " + std::to_string(elem_count) + " elements from " +
		                     std::to_string(n) + " bytes, more than deflate can inflate");
	}
	if (elem_count > MAX_DECODED_BLOB / (type.bits / BITS_PER_BYTE)) {
		throw SraUnsupportedError("an izip blob that decodes to " + std::to_string(elem_count) +
		                          " elements, more than " + std::to_string(MAX_DECODED_BLOB >> GIB_SHIFT) + " GiB");
	}
	switch (type.bits) {
	case U8_BITS:
		return IrzipDecode<uint8_t>(frame, src, n, elem_count);
	case U16_BITS:
		return IrzipDecode<uint16_t>(frame, src, n, elem_count);
	case U32_BITS:
		return IrzipDecode<uint32_t>(frame, src, n, elem_count);
	case U64_BITS:
		return IrzipDecode<uint64_t>(frame, src, n, elem_count);
	default:
		throw std::logic_error("Iunzip: width checked above"); // unreachable
	}
}

// The decoded elements are in host order for a little-endian host, which every target miint ships
// is (x86_64, arm64, wasm32): a big-endian blob is swapped, a little-endian one is not.
DecodedBlob DecodeBlob(VdbBlob blob, const std::string &encoding, ElementType type) {
	DecodedBlob out;
	Encoding enc = EncodingOf(encoding);
	const bool skipped = !blob.frames.empty() && (blob.frames[0].flags & HEADER_FLAG_SKIPPED);
	if (enc == Encoding::Raw && !blob.frames.empty() && !skipped) {
		// The column records no encoding, so an encoded blob cannot be decoded as what it says it
		// is; reading its bytes as elements would return garbage.
		throw SraFormatError("a blob of a column stored without encoding carries an encoding header");
	}
	bool passthrough = enc == Encoding::Raw || skipped;
	uint64_t bits;
	if (passthrough) {
		if (blob.frames.size() > 1) {
			throw SraUnsupportedError("a blob with " + std::to_string(blob.frames.size()) + " encoding stages");
		}
		out.data.assign(blob.Payload(), blob.Payload() + blob.PayloadSize());
		bits = blob.payload_bits;
		if (blob.big_endian && type.bits > BITS_PER_BYTE) {
			ByteSwap(out.data, type.bits / BITS_PER_BYTE);
		}
	} else {
		if (blob.frames.size() != 1) {
			throw SraFormatError("an encoded blob has " + std::to_string(blob.frames.size()) + " encoding headers");
		}
		const HeaderFrame &frame = blob.frames[0];
		if (enc == Encoding::Zip) {
			out.data = Unzip(frame, blob.Payload(), blob.PayloadSize(), bits);
			if (blob.big_endian && type.bits > BITS_PER_BYTE) {
				ByteSwap(out.data, type.bits / BITS_PER_BYTE);
			}
		} else {
			if (type.bits < BITS_PER_BYTE) {
				throw SraUnsupportedError("izip_encoding on a column of " + std::to_string(type.bits) +
				                          "-bit elements");
			}
			if (frame.osize % (type.bits / BITS_PER_BYTE) != 0) {
				throw SraFormatError("an izip blob's size is not a whole number of elements");
			}
			uint64_t count = uint64_t(frame.osize) / (type.bits / BITS_PER_BYTE);
			out.data = Iunzip(frame, blob.Payload(), blob.PayloadSize(), type, count);
			bits = count * type.bits;
		}
	}
	if (bits % type.bits != 0) {
		throw SraFormatError("a blob does not decode to a whole number of elements");
	}
	out.elem_count = bits / type.bits;
	// ElementsSpanned saturates, so an offset near 2^64 cannot wrap past the bound.
	if (blob.page_map.ElementsSpanned() > out.elem_count) {
		throw SraFormatError("a blob's page map addresses elements past the " + std::to_string(out.elem_count) +
		                     " it decodes to");
	}
	out.page_map = std::move(blob.page_map);
	return out;
}

} // namespace sra
} // namespace miint
