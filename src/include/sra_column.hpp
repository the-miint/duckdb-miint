#pragma once

#include "sra_archive.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace miint {
namespace sra {

// Where one blob of a physical column lives: rows [start_id, start_id + id_count) are stored in
// `size` bytes at data offset pg * page_size (checksum included).
struct BlobLoc {
	int64_t start_id = 0;
	uint64_t id_count = 0;
	uint64_t pg = 0;
	uint32_t size = 0;
};

enum class ChecksumType { None, Crc32, Md5 };

// A column's index files (idx, idx0, idx1, idx2) resolved into the list of its blobs. Only
// column header version 3 (every column in the survey of public runs) is read; others raise
// SraUnsupportedError.
struct ColumnIndex {
	uint64_t data_eof = 0;
	uint32_t page_size = 0;
	ChecksumType checksum = ChecksumType::None;
	// Sorted by start_id, non-overlapping.
	std::vector<BlobLoc> blobs;

	// The bytes of each file, as stored (idx0 may be empty). Trailing bytes past the logical
	// ends the header records (an MD5 append trailer) are ignored.
	static ColumnIndex Parse(const std::vector<uint8_t> &idx, const std::vector<uint8_t> &idx0,
	                         const std::vector<uint8_t> &idx1, const std::vector<uint8_t> &idx2);
	// The blob holding row id; throws SraFormatError if no blob does.
	const BlobLoc &Locate(int64_t row_id) const;
};

// VDB's CRC32: polynomial 0x04C11DB7, MSB first, initial value 0, no reflection or final XOR.
// Not the zlib CRC.
uint32_t Crc32Ncbi(const uint8_t *data, size_t n, uint32_t crc = 0);

// One frame of a blob's encoding-header stack: the parameters of one decoding stage.
struct HeaderFrame {
	uint8_t flags = 0; // bit 0: the stage was skipped, the data passes through unchanged
	uint8_t version = 0;
	int64_t fmt = 0;
	int64_t osize = 0; // bytes the stage decodes to
	std::vector<uint8_t> ops;
	std::vector<int64_t> args;
};

// How a blob's rows map onto its decoded elements. Kept as the runs the page map stores, not per
// row: a blob of millions of empty rows (VDB cuts blobs by size) is a few bytes on disk and must
// stay small here.
struct PageMap {
	// Rows [first_row, first_row + rows) hold `length` elements each. The first starts at element
	// `offset` and each next one `step` elements after it: `length` when every row has its own
	// record, 0 when the rows share one.
	struct Run {
		uint64_t first_row = 0;
		uint64_t rows = 0;
		uint32_t length = 0;
		uint64_t offset = 0;
		uint64_t step = 0;
	};
	std::vector<Run> runs; // in row order, covering rows [0, row_count)
	uint64_t row_count = 0;

	// Row i's length and first element; i < row_count. `hint` is the run a caller found last:
	// rows are read in order, so it usually holds the next one too and saves the search.
	void Row(uint64_t i, uint32_t &length, uint64_t &offset, size_t &hint) const;
	void Row(uint64_t i, uint32_t &length, uint64_t &offset) const {
		size_t hint = 0;
		Row(i, length, offset, hint);
	}
	// The end of the furthest row, or UINT64_MAX when that sum would overflow: the bound a blob's
	// decoded elements must meet.
	uint64_t ElementsSpanned() const;
};

// A blob as stored in a column's data file, framing parsed but payload still encoded.
struct VdbBlob {
	bool big_endian = false;
	std::vector<HeaderFrame> frames; // empty for a v1 blob, or a v2 blob with no encoding
	PageMap page_map;
	std::vector<uint8_t> bytes; // the whole blob, checksum removed
	size_t payload_offset = 0;
	uint64_t payload_bits = 0; // payload length in bits, padding excluded

	const uint8_t *Payload() const {
		return bytes.data() + payload_offset;
	}
	size_t PayloadSize() const {
		return bytes.size() - payload_offset;
	}
};

// Parses a blob's framing: the v1 or v2 header, the encoding headers, and the page map for its
// row_count rows. elem_bits is the bit width of the column's stored element type; a v1 blob
// derives its page map from it.
VdbBlob ParseVdbBlob(std::vector<uint8_t> bytes, uint64_t row_count, uint32_t elem_bits);

// A stored element type as the decoders need it: width and signedness.
struct ElementType {
	uint32_t bits = 8;
	bool is_signed = false;
};
// Throws SraUnsupportedError for a type it does not know.
ElementType ElementTypeOf(const std::string &type);
// Bit width of a stored element type, e.g. 2 for INSDC:2na:packed.
uint32_t ElementBits(const std::string &type);

// A column directory of a table: its declared type and encoding, its index, and its blobs.
class PhysicalColumn {
public:
	PhysicalColumn(const SraTable &table, const std::string &name);

	const std::string &Name() const {
		return name_;
	}
	// The stored element type (the `type` attribute of the column's schema node).
	const std::string &Type() const {
		return type_;
	}
	// The physical encoding expression, e.g. "<INSDC:coord:len>izip_encoding#1"; empty when the
	// column is stored unencoded.
	const std::string &Encoding() const {
		return encoding_;
	}
	const ColumnIndex &Index() const {
		return index_;
	}
	// The blob's bytes with its checksum verified and removed. Throws SraFormatError on a
	// checksum mismatch, SraUnsupportedError for MD5 checksums.
	std::vector<uint8_t> ReadBlob(const BlobLoc &loc) const;
	VdbBlob ReadVdbBlob(const BlobLoc &loc) const;
	// Hints to the byte source that the blobs holding rows [first, end) will be read soon.
	void PrefetchRows(int64_t first, int64_t end) const;

private:
	const KarArchive *archive_;
	std::string name_;
	std::string dir_;
	std::string type_;
	std::string encoding_;
	ColumnIndex index_;
};

} // namespace sra
} // namespace miint
