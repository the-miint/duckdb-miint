#include "sra_codecs.hpp"

#include <cstdint>

#include <zlib.h>

#include <cstring>
#include <type_traits>

namespace miint {
namespace sra {

Inflated InflateRaw(const uint8_t *src, size_t n, uint8_t *dst, size_t dsize, const char *what) {
	if (n > UINT32_MAX || dsize > UINT32_MAX) {
		throw SraFormatError(std::string(what) + " is larger than 4 GiB");
	}
	z_stream zs {};
	if (inflateInit2(&zs, -15) != Z_OK) {
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

template <typename U>
std::vector<uint8_t> IrzipDecode(const HeaderFrame &frame, const uint8_t *src, size_t n, uint64_t count) {
	if (frame.ops.empty() || frame.args.empty()) {
		throw SraFormatError("an izip blob header lacks its parameters");
	}
	uint8_t planes = frame.ops[0];
	int64_t min[2] = {frame.args[0], 0};
	int64_t slope[2] = {frame.args.size() > 1 ? frame.args[1] : 0, 0};
	int series_count = 1;
	if (frame.args.size() >= 4) {
		min[1] = frame.args[2];
		slope[1] = frame.args[3];
		series_count = 2;
	}
	std::vector<U> Y(count, 0);
	std::vector<uint8_t> plane(count);
	size_t used = 0;
	for (unsigned p = 0; p < 8; p++) {
		if (!(planes & (1u << p))) {
			continue;
		}
		if (p >= sizeof(U)) {
			throw SraFormatError("an izip blob has a byte plane wider than its type");
		}
		used += InflateExact(src + used, n - used, plane.data(), count, "an izip byte plane");
		for (size_t i = 0; i < count; i++) {
			Y[i] = U(Y[i] | (U(plane[i]) << (8 * p)));
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
	if (name == "izip_encoding#1") {
		return Encoding::Izip;
	}
	if (name.compare(0, 14, "zip_encoding#1") == 0) {
		return Encoding::Zip;
	}
	throw SraUnsupportedError("encoding '" + expr + "'");
}

} // namespace

std::vector<uint8_t> Unzip(const HeaderFrame &frame, const uint8_t *src, size_t n, uint64_t &out_bits) {
	if (frame.version != 1 && frame.version != 2) {
		throw SraUnsupportedError("zip blob header version " + std::to_string(frame.version));
	}
	if (frame.osize < 0 || uint64_t(frame.osize) > uint64_t(n) * MAX_DEFLATE_RATIO) {
		throw SraFormatError("a zip blob claims " + std::to_string(frame.osize) + " bytes from " + std::to_string(n) +
		                     ", more than deflate can inflate");
	}
	std::vector<uint8_t> out(size_t(frame.osize));
	InflateExact(src, n, out.data(), out.size(), "a zip blob");
	out_bits = uint64_t(out.size()) * 8;
	if (frame.version == 2) {
		if (frame.args.empty()) {
			throw SraFormatError("a zip blob header lacks its trailing-bit count");
		}
		int64_t trailing = frame.args[0];
		if (trailing < 0 || trailing > 7 || (trailing != 0 && out.empty())) {
			throw SraFormatError("a zip blob header records " + std::to_string(trailing) + " trailing bits");
		}
		if (trailing != 0) {
			out_bits -= uint64_t(8 - trailing);
		}
	}
	return out;
}

std::vector<uint8_t> Iunzip(const HeaderFrame &frame, const uint8_t *src, size_t n, ElementType type,
                            uint64_t elem_count) {
	if (frame.version == 0) {
		throw SraUnsupportedError("legacy izip (blob header version 0), which no public run surveyed used");
	}
	if (frame.version > 3) {
		throw SraUnsupportedError("izip blob header version " + std::to_string(frame.version));
	}
	// Each byte plane inflates to elem_count bytes from at most n. A blob with no planes is all
	// offset and slope (a constant column), so only the size cap bounds it.
	const bool has_planes = !frame.ops.empty() && frame.ops[0] != 0;
	if (has_planes && elem_count > uint64_t(n) * MAX_DEFLATE_RATIO) {
		throw SraFormatError("an izip blob claims " + std::to_string(elem_count) + " elements from " +
		                     std::to_string(n) + " bytes, more than deflate can inflate");
	}
	if (elem_count > MAX_DECODED_BLOB / (type.bits / 8)) {
		throw SraUnsupportedError("an izip blob that decodes to " + std::to_string(elem_count) +
		                          " elements, more than " + std::to_string(MAX_DECODED_BLOB >> 30) + " GiB");
	}
	switch (type.bits) {
	case 8:
		return IrzipDecode<uint8_t>(frame, src, n, elem_count);
	case 16:
		return IrzipDecode<uint16_t>(frame, src, n, elem_count);
	case 32:
		return IrzipDecode<uint32_t>(frame, src, n, elem_count);
	case 64:
		return IrzipDecode<uint64_t>(frame, src, n, elem_count);
	default:
		throw SraUnsupportedError("izip on " + std::to_string(type.bits) + "-bit elements");
	}
}

// The decoded elements are in host order for a little-endian host, which every target miint ships
// is (x86_64, arm64, wasm32): a big-endian blob is swapped, a little-endian one is not.
DecodedBlob DecodeBlob(VdbBlob blob, const std::string &encoding, ElementType type) {
	DecodedBlob out;
	Encoding enc = EncodingOf(encoding);
	const bool skipped = !blob.frames.empty() && (blob.frames[0].flags & 1);
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
		if (blob.big_endian && type.bits >= 16) {
			ByteSwap(out.data, type.bits / 8);
		}
	} else {
		if (blob.frames.size() != 1) {
			throw SraFormatError("an encoded blob has " + std::to_string(blob.frames.size()) + " encoding headers");
		}
		const HeaderFrame &frame = blob.frames[0];
		if (enc == Encoding::Zip) {
			out.data = Unzip(frame, blob.Payload(), blob.PayloadSize(), bits);
			if (blob.big_endian && type.bits >= 16) {
				ByteSwap(out.data, type.bits / 8);
			}
		} else {
			if (type.bits < 8) {
				throw SraUnsupportedError("izip_encoding on a column of " + std::to_string(type.bits) +
				                          "-bit elements");
			}
			if (frame.osize % (type.bits / 8) != 0) {
				throw SraFormatError("an izip blob's size is not a whole number of elements");
			}
			uint64_t count = uint64_t(frame.osize) / (type.bits / 8);
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
