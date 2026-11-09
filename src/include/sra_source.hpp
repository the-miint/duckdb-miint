#pragma once

#include <cstdint>
#include <fstream>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace miint {
namespace sra {

// The bytes of a .sra file are damaged or truncated: a structure the format defines could not
// be read. Distinct from SraUnsupportedError so a caller can tell a bad download from a run
// miint deliberately does not read.
class SraFormatError : public std::runtime_error {
public:
	using std::runtime_error::runtime_error;
};

// The file is well formed but holds something miint does not read: an archive variant never
// seen in public runs, or a kind of run (cSRA, SRA Lite, ...) that is rejected by design. The
// message names what was found, without a prefix; the table function adds the run.
class SraUnsupportedError : public std::runtime_error {
public:
	using std::runtime_error::runtime_error;
};

// Random access to the bytes of one .sra file. ReadAt must be safe to call from several threads
// at once: the table functions decode different row ranges of one run in parallel.
class SraByteSource {
public:
	virtual ~SraByteSource() = default;
	virtual uint64_t Size() const = 0;
	// Reads exactly n bytes at pos, or throws SraFormatError if the file is shorter.
	virtual void ReadAt(uint64_t pos, uint8_t *buf, size_t n) const = 0;
	// A hint that these bytes will be read soon. A source that fetches ahead starts on them; others
	// ignore it.
	virtual void Prefetch(uint64_t pos, uint64_t n) const {
	}

	std::vector<uint8_t> Read(uint64_t pos, size_t n) const;
};

// A local file. Used by the C++ tests; the table functions read through DuckDB's file system.
class SraFileSource : public SraByteSource {
public:
	explicit SraFileSource(const std::string &path);
	uint64_t Size() const override;
	void ReadAt(uint64_t pos, uint8_t *buf, size_t n) const override;

private:
	mutable std::mutex mutex_;
	mutable std::ifstream in_;
	uint64_t size_;
};

// Bytes held in memory: lets a test damage a real file, or build a structure by hand.
class SraMemorySource : public SraByteSource {
public:
	explicit SraMemorySource(std::vector<uint8_t> bytes);
	uint64_t Size() const override;
	void ReadAt(uint64_t pos, uint8_t *buf, size_t n) const override;

private:
	std::vector<uint8_t> bytes_;
};

// Little-endian loads. Every archive variant seen in public runs is little-endian, and a
// big-endian one is rejected when the archive is opened, so the format code reads LE only.
uint16_t LoadU16(const uint8_t *p);
uint32_t LoadU32(const uint8_t *p);
uint64_t LoadU64(const uint8_t *p);

} // namespace sra
} // namespace miint
