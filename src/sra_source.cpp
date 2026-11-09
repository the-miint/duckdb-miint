#include "sra_source.hpp"

#include <algorithm>
#include <cstring>
#include <stdexcept>

namespace miint {
namespace sra {

std::vector<uint8_t> SraByteSource::Read(uint64_t pos, size_t n) const {
	std::vector<uint8_t> out(n);
	if (n > 0) {
		ReadAt(pos, out.data(), n);
	}
	return out;
}

namespace {
void CheckRange(uint64_t pos, size_t n, uint64_t size) {
	if (pos > size || n > size - pos) {
		throw SraFormatError("read of " + std::to_string(n) + " bytes at offset " + std::to_string(pos) +
		                     " runs past the end of the file (" + std::to_string(size) + " bytes)");
	}
}
} // namespace

SraFileSource::SraFileSource(const std::string &path) : in_(path, std::ios::binary) {
	if (!in_) {
		throw SraFormatError("cannot open " + path);
	}
	in_.seekg(0, std::ios::end);
	size_ = uint64_t(in_.tellg());
}

uint64_t SraFileSource::Size() const {
	return size_;
}

void SraFileSource::ReadAt(uint64_t pos, uint8_t *buf, size_t n) const {
	CheckRange(pos, n, size_);
	std::lock_guard<std::mutex> lock(mutex_);
	in_.seekg(std::streamoff(pos));
	in_.read(reinterpret_cast<char *>(buf), std::streamsize(n));
	if (size_t(in_.gcount()) != n) {
		in_.clear();
		throw SraFormatError("short read at offset " + std::to_string(pos));
	}
}

SraMemorySource::SraMemorySource(std::vector<uint8_t> bytes) : bytes_(std::move(bytes)) {
}

uint64_t SraMemorySource::Size() const {
	return bytes_.size();
}

void SraMemorySource::ReadAt(uint64_t pos, uint8_t *buf, size_t n) const {
	CheckRange(pos, n, bytes_.size());
	if (n > 0) {
		std::memcpy(buf, bytes_.data() + pos, n);
	}
}

uint16_t LoadU16(const uint8_t *p) {
	return uint16_t(p[0] | (p[1] << 8));
}

uint32_t LoadU32(const uint8_t *p) {
	return uint32_t(p[0]) | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16) | (uint32_t(p[3]) << 24);
}

uint64_t LoadU64(const uint8_t *p) {
	return uint64_t(LoadU32(p)) | (uint64_t(LoadU32(p + 4)) << 32);
}

} // namespace sra
} // namespace miint
