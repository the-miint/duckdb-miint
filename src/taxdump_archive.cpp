#include "taxdump_archive.hpp"

#include "microtar.h"
#include <zlib.h>

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace miint {

namespace {

// In-memory stream backing for microtar. microtar's read callback must behave
// like fread: read `size` bytes from the current position and advance it; seek
// repositions. The invariant pos <= len is preserved by both callbacks, so the
// `len - pos` subtractions never underflow.
struct MemStream {
	const char *buf;
	size_t len;
	size_t pos;
};

} // namespace

extern "C" {
static int MtarMemRead(mtar_t *tar, void *data, unsigned size) {
	auto *ms = static_cast<MemStream *>(tar->stream);
	if (static_cast<size_t>(size) > ms->len - ms->pos) {
		return MTAR_EREADFAIL;
	}
	std::memcpy(data, ms->buf + ms->pos, size);
	ms->pos += size;
	return MTAR_ESUCCESS;
}

static int MtarMemSeek(mtar_t *tar, unsigned pos) {
	auto *ms = static_cast<MemStream *>(tar->stream);
	if (static_cast<size_t>(pos) > ms->len) {
		return MTAR_ESEEKFAIL;
	}
	ms->pos = pos;
	return MTAR_ESUCCESS;
}

static int MtarMemClose(mtar_t *) {
	return MTAR_ESUCCESS;
}
} // extern "C"

std::string TaxdumpArchive::Gunzip(const std::string &gz_bytes) {
	if (gz_bytes.empty()) {
		return {};
	}

	z_stream zs;
	std::memset(&zs, 0, sizeof(zs));
	// 16 + MAX_WBITS selects gzip framing (as opposed to raw zlib).
	if (inflateInit2(&zs, 16 + MAX_WBITS) != Z_OK) {
		throw std::runtime_error("taxdump: failed to initialize gzip decompressor");
	}

	// Feed the input in <= 1 GiB windows so gz_bytes.size() > UINT_MAX never
	// overflows the uInt avail_in field. taxdump.tar.gz is a single gzip member,
	// so decoding is complete only when inflate reports Z_STREAM_END; if the input
	// runs out first the stream was truncated (e.g. a half-finished download) and
	// we fail loud rather than return a partial buffer.
	const char *in_ptr = gz_bytes.data();
	size_t in_remaining = gz_bytes.size();
	constexpr size_t kFeedCap = static_cast<size_t>(1) << 30; // 1 GiB, well under UINT_MAX

	std::string out;
	std::vector<char> chunk(static_cast<size_t>(1) << 20); // 1 MiB output window
	int ret = Z_OK;
	while (ret != Z_STREAM_END) {
		if (zs.avail_in == 0) {
			if (in_remaining == 0) {
				inflateEnd(&zs);
				throw std::runtime_error("taxdump: gzip stream is truncated or incomplete");
			}
			uInt take = static_cast<uInt>(std::min(in_remaining, kFeedCap));
			zs.next_in = reinterpret_cast<Bytef *>(const_cast<char *>(in_ptr));
			zs.avail_in = take;
			in_ptr += take;
			in_remaining -= take;
		}

		zs.next_out = reinterpret_cast<Bytef *>(chunk.data());
		zs.avail_out = static_cast<uInt>(chunk.size());
		ret = inflate(&zs, Z_NO_FLUSH);
		if (ret != Z_OK && ret != Z_STREAM_END) {
			std::string msg = zs.msg ? zs.msg : "corrupt stream";
			inflateEnd(&zs);
			throw std::runtime_error("taxdump: gzip inflate failed: " + msg);
		}
		out.append(chunk.data(), chunk.size() - zs.avail_out);
	}

	inflateEnd(&zs);
	return out;
}

std::unordered_map<std::string, std::string> TaxdumpArchive::ExtractTarMembers(const std::string &tar_bytes,
                                                                               const std::vector<std::string> &names) {
	MemStream ms {tar_bytes.data(), tar_bytes.size(), 0};
	mtar_t tar;
	std::memset(&tar, 0, sizeof(tar));
	tar.read = MtarMemRead;
	tar.seek = MtarMemSeek;
	tar.close = MtarMemClose;
	tar.stream = &ms;

	std::unordered_map<std::string, std::string> out;
	for (const auto &name : names) {
		mtar_header_t h;
		int err = mtar_find(&tar, name.c_str(), &h);
		if (err == MTAR_ENOTFOUND) {
			continue;
		}
		if (err != MTAR_ESUCCESS) {
			throw std::runtime_error("taxdump: tar read error for '" + name + "': " + mtar_strerror(err));
		}
		std::string data;
		data.resize(h.size);
		if (h.size > 0) {
			err = mtar_read_data(&tar, data.data(), h.size);
			if (err != MTAR_ESUCCESS) {
				throw std::runtime_error("taxdump: tar data read error for '" + name + "': " + mtar_strerror(err));
			}
		}
		out.emplace(name, std::move(data));
	}
	return out;
}

TaxdumpFiles TaxdumpArchive::ExtractTaxdump(const std::string &targz_bytes) {
	std::string tar_bytes = Gunzip(targz_bytes);
	auto members = ExtractTarMembers(tar_bytes, {"nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"});

	if (members.find("nodes.dmp") == members.end() || members.find("names.dmp") == members.end()) {
		throw std::runtime_error("taxdump: archive is missing required members (nodes.dmp and/or names.dmp)");
	}

	auto take = [&](const char *n) -> std::string {
		auto it = members.find(n);
		return it == members.end() ? std::string() : std::move(it->second);
	};
	TaxdumpFiles files;
	files.nodes = take("nodes.dmp");
	files.names = take("names.dmp");
	files.merged = take("merged.dmp");
	files.delnodes = take("delnodes.dmp");
	return files;
}

} // namespace miint
