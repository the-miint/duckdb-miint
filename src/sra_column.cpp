#include "sra_column.hpp"

#include "sra_codecs.hpp"

#include <algorithm>
#include <array>
#include <map>

namespace miint {
namespace sra {

namespace {

// The column header (idx): the KDB file header, then the data and idx2 ends, the idx0 entry and
// idx1 block counts, the page size and the checksum type.
constexpr uint32_t COLUMN_HEADER_VERSION = 3;
constexpr size_t COLUMN_HEADER_SIZE = 40;
constexpr size_t COLUMN_DATA_EOF_AT = 8;
constexpr size_t COLUMN_IDX2_EOF_AT = 16;
constexpr size_t COLUMN_IDX0_COUNT_AT = 24;
constexpr size_t COLUMN_BLOCK_COUNT_AT = 28;
constexpr size_t COLUMN_PAGE_SIZE_AT = 32;
constexpr size_t COLUMN_CHECKSUM_AT = 36;
enum StoredChecksum : uint8_t { STORED_NO_CHECKSUM = 0, STORED_CRC32 = 1, STORED_MD5 = 2 };

// A location record, 24 bytes: an idx0 entry locates one blob, an idx1 entry one idx2 block. Both
// hold a page, a 32-bit field of packed bits, a row count, and the first row id.
constexpr size_t LOC_SIZE = 24;
constexpr size_t LOC_GEN_AT = 8;
constexpr size_t LOC_ROWS_AT = 12;
constexpr size_t LOC_START_ID_AT = 16;
// An idx0 entry's packed field: 31 bits of blob size, then a flag that the entry removes a blob.
constexpr uint32_t BLOB_LOC_SIZE_MASK = (uint32_t(1) << 31) - 1;
constexpr int BLOB_LOC_REMOVE_SHIFT = 31;
// An idx1 entry's packed field: 27 bits of block size, the 2-bit row-id and offset block types,
// then a flag that the block is compressed.
constexpr uint32_t BLOCK_LOC_SIZE_MASK = (uint32_t(1) << 27) - 1;
constexpr int BLOCK_LOC_ID_TYPE_SHIFT = 27;
constexpr int BLOCK_LOC_PG_TYPE_SHIFT = 29;
constexpr uint32_t BLOCK_LOC_TYPE_MASK = 3;
constexpr int BLOCK_LOC_COMPRESSED_SHIFT = 31;
// An idx2 block's layout: a Magnitude series starts with its first value (8 bytes) and lists a
// 4-byte span per blob; a Predictable offset series is its first offset and the one blob size.
constexpr size_t BLOCK_MAGNITUDE_HEADER_SIZE = sizeof(uint64_t);
constexpr size_t BLOCK_MAGNITUDE_ENTRY_SIZE = sizeof(uint32_t);
constexpr size_t BLOCK_PREDICTABLE_PG_SIZE = sizeof(uint64_t) + sizeof(uint32_t);

void CheckKdbHeader(const std::vector<uint8_t> &file, const std::string &what) {
	if (file.size() < KDB_FILE_HEADER_SIZE) {
		throw SraFormatError(what + " is shorter than its header");
	}
	CheckByteOrder(LoadU32(file.data()), what);
}

std::string RowsText(int64_t start, uint64_t count) {
	return std::to_string(start) + ".." + std::to_string(start + int64_t(count) - 1);
}

// idx2 block types: how a block stores its blobs' row ids and data offsets. Only Magnitude (2)
// and Predictable (3) occur in public runs.
enum BlockType : uint32_t { RANDOM = 0, UNIFORM = 1, MAGNITUDE = 2, PREDICTABLE = 3 };

void DecodeBlock(const std::vector<uint8_t> &idx1, size_t at, const std::vector<uint8_t> &idx2, uint64_t idx2_eof,
                 std::vector<BlobLoc> &out) {
	uint64_t pg = LoadU64(idx1.data() + at);
	uint32_t gen = LoadU32(idx1.data() + at + LOC_GEN_AT);
	uint32_t id_range = LoadU32(idx1.data() + at + LOC_ROWS_AT);
	int64_t start_id = int64_t(LoadU64(idx1.data() + at + LOC_START_ID_AT));
	uint32_t size = gen & BLOCK_LOC_SIZE_MASK;
	uint32_t id_type = (gen >> BLOCK_LOC_ID_TYPE_SHIFT) & BLOCK_LOC_TYPE_MASK;
	uint32_t pg_type = (gen >> BLOCK_LOC_PG_TYPE_SHIFT) & BLOCK_LOC_TYPE_MASK;
	if (gen >> BLOCK_LOC_COMPRESSED_SHIFT) {
		throw SraUnsupportedError("a compressed idx2 block, which no public run surveyed used");
	}
	if (id_type < MAGNITUDE || pg_type < MAGNITUDE) {
		throw SraUnsupportedError("an idx2 block of type " + std::to_string(id_type) + "/" + std::to_string(pg_type) +
		                          ", which no public run surveyed used");
	}
	if (id_range == 0) {
		throw SraFormatError("an idx2 block covers no rows");
	}

	uint64_t count;
	const uint8_t *blk;
	size_t blk_size;
	if (id_type == PREDICTABLE && pg_type == PREDICTABLE) {
		// Both predictable: the locator's size field is the blob count, and the block is just
		// {first data offset, blob size}.
		count = size;
		blk_size = BLOCK_PREDICTABLE_PG_SIZE;
	} else {
		blk_size = size;
	}
	if (pg > idx2_eof || blk_size > idx2_eof - pg || idx2_eof > idx2.size()) {
		throw SraFormatError("an idx2 block lies outside idx2");
	}
	blk = idx2.data() + pg;

	size_t id_hdr = id_type == MAGNITUDE ? BLOCK_MAGNITUDE_HEADER_SIZE : 0;
	size_t pg_hdr = pg_type == MAGNITUDE ? BLOCK_MAGNITUDE_HEADER_SIZE : BLOCK_PREDICTABLE_PG_SIZE;
	size_t id_each = id_type == MAGNITUDE ? BLOCK_MAGNITUDE_ENTRY_SIZE : 0;
	size_t pg_each = pg_type == MAGNITUDE ? BLOCK_MAGNITUDE_ENTRY_SIZE : 0;
	if (!(id_type == PREDICTABLE && pg_type == PREDICTABLE)) {
		size_t per = id_each + pg_each;
		if (blk_size < id_hdr + pg_hdr || (blk_size - id_hdr - pg_hdr) % per != 0) {
			throw SraFormatError("an idx2 block's size does not match its layout");
		}
		count = (blk_size - id_hdr - pg_hdr) / per;
	}
	if (count == 0) {
		throw SraFormatError("an idx2 block holds no blobs");
	}
	// Every blob holds a row at least, so the block's rows bound its blob count before anything is
	// sized from it.
	if (count > id_range) {
		throw SraFormatError("an idx2 block lists more blobs (" + std::to_string(count) + ") than rows (" +
		                     std::to_string(id_range) + ")");
	}
	const uint8_t *id_s = blk + id_hdr + pg_hdr;
	const uint8_t *pg_s = id_s + id_each * count;

	std::vector<BlobLoc> blobs(count);
	// Row ids: Magnitude stores the first id and each blob's span; Predictable spreads the
	// block's id range evenly.
	if (id_type == MAGNITUDE) {
		int64_t next = int64_t(LoadU64(blk));
		for (uint64_t i = 0; i < count; i++) {
			blobs[i].start_id = next;
			blobs[i].id_count = LoadU32(id_s + BLOCK_MAGNITUDE_ENTRY_SIZE * i);
			next += int64_t(blobs[i].id_count);
		}
	} else {
		if (id_range % count != 0) {
			throw SraFormatError("an idx2 block's rows do not divide evenly among its blobs");
		}
		uint64_t per = id_range / count;
		for (uint64_t i = 0; i < count; i++) {
			blobs[i].start_id = start_id + int64_t(i * per);
			blobs[i].id_count = per;
		}
	}
	// Data offsets: Magnitude stores the first offset and each blob's size; Predictable has one
	// size for every blob.
	if (pg_type == MAGNITUDE) {
		uint64_t next = LoadU64(blk + id_hdr);
		for (uint64_t i = 0; i < count; i++) {
			blobs[i].pg = next;
			blobs[i].size = LoadU32(pg_s + BLOCK_MAGNITUDE_ENTRY_SIZE * i);
			next += blobs[i].size;
		}
	} else {
		uint64_t first = LoadU64(blk + id_hdr);
		uint32_t each = LoadU32(blk + id_hdr + sizeof(uint64_t));
		for (uint64_t i = 0; i < count; i++) {
			blobs[i].pg = first + i * each;
			blobs[i].size = each;
		}
	}
	if (blobs.front().start_id != start_id ||
	    blobs.back().start_id + int64_t(blobs.back().id_count) != start_id + int64_t(id_range)) {
		throw SraFormatError("an idx2 block's blobs do not cover rows " + RowsText(start_id, id_range));
	}
	out.insert(out.end(), blobs.begin(), blobs.end());
}

} // namespace

ColumnIndex ColumnIndex::Parse(const std::vector<uint8_t> &idx, const std::vector<uint8_t> &idx0,
                               const std::vector<uint8_t> &idx1, const std::vector<uint8_t> &idx2) {
	CheckKdbHeader(idx, "a column header");
	uint32_t version = LoadU32(idx.data() + sizeof(uint32_t));
	if (version != COLUMN_HEADER_VERSION) {
		throw SraUnsupportedError("column header version " + std::to_string(version) +
		                          " is not supported (public runs use version 3)");
	}
	if (idx.size() < COLUMN_HEADER_SIZE) {
		throw SraFormatError("a column header is truncated");
	}
	ColumnIndex ix;
	ix.data_eof = LoadU64(idx.data() + COLUMN_DATA_EOF_AT);
	uint64_t idx2_eof = LoadU64(idx.data() + COLUMN_IDX2_EOF_AT);
	uint32_t idx0_count = LoadU32(idx.data() + COLUMN_IDX0_COUNT_AT);
	uint32_t num_blocks = LoadU32(idx.data() + COLUMN_BLOCK_COUNT_AT);
	ix.page_size = LoadU32(idx.data() + COLUMN_PAGE_SIZE_AT);
	const uint8_t checksum = idx[COLUMN_CHECKSUM_AT];
	switch (checksum) {
	case STORED_NO_CHECKSUM:
		ix.checksum = ChecksumType::None;
		break;
	case STORED_CRC32:
		ix.checksum = ChecksumType::Crc32;
		break;
	case STORED_MD5:
		ix.checksum = ChecksumType::Md5;
		break;
	default:
		throw SraFormatError("a column header names checksum type " + std::to_string(checksum));
	}
	if (ix.page_size == 0) {
		throw SraFormatError("a column header has page size 0");
	}

	// idx0 is a journal: an entry for a row range already listed replaces it.
	std::map<int64_t, BlobLoc> journal;
	if (idx0.size() / LOC_SIZE < idx0_count) {
		throw SraFormatError("idx0 holds fewer entries than the column header records");
	}
	for (size_t i = 0; i < idx0_count; i++) {
		const uint8_t *p = idx0.data() + i * LOC_SIZE;
		BlobLoc loc;
		loc.pg = LoadU64(p);
		uint32_t gen = LoadU32(p + LOC_GEN_AT);
		loc.size = gen & BLOB_LOC_SIZE_MASK;
		loc.id_count = LoadU32(p + LOC_ROWS_AT);
		loc.start_id = int64_t(LoadU64(p + LOC_START_ID_AT));
		if (gen >> BLOB_LOC_REMOVE_SHIFT) {
			throw SraUnsupportedError("an idx0 entry removing a blob");
		}
		if (loc.id_count == 0) {
			throw SraFormatError("an idx0 entry covers no rows");
		}
		// Entries never overlap one another unless they cover exactly the same rows.
		int64_t end = loc.start_id + int64_t(loc.id_count);
		auto next = journal.lower_bound(loc.start_id);
		const BlobLoc *neighbours[] = {next != journal.end() ? &next->second : nullptr,
		                               next != journal.begin() ? &std::prev(next)->second : nullptr};
		for (const BlobLoc *e : neighbours) {
			if (e == nullptr) {
				continue;
			}
			bool overlaps = e->start_id < end && loc.start_id < e->start_id + int64_t(e->id_count);
			if (overlaps && (e->start_id != loc.start_id || e->id_count != loc.id_count)) {
				throw SraFormatError("idx0 entries for rows " + RowsText(e->start_id, e->id_count) + " and " +
				                     RowsText(loc.start_id, loc.id_count) + " overlap");
			}
		}
		journal[loc.start_id] = loc;
	}

	if (num_blocks > 0) {
		CheckKdbHeader(idx1, "idx1");
		if ((idx1.size() - KDB_FILE_HEADER_SIZE) / LOC_SIZE < num_blocks) {
			throw SraFormatError("idx1 holds fewer blocks than the column header records");
		}
	}
	std::vector<BlobLoc> from_blocks;
	for (size_t b = 0; b < num_blocks; b++) {
		DecodeBlock(idx1, KDB_FILE_HEADER_SIZE + b * LOC_SIZE, idx2, idx2_eof, from_blocks);
	}

	ix.blobs.reserve(journal.size() + from_blocks.size());
	for (const auto &kv : journal) {
		ix.blobs.push_back(kv.second);
	}
	ix.blobs.insert(ix.blobs.end(), from_blocks.begin(), from_blocks.end());
	std::sort(ix.blobs.begin(), ix.blobs.end(),
	          [](const BlobLoc &a, const BlobLoc &b) { return a.start_id < b.start_id; });
	for (size_t i = 0; i < ix.blobs.size(); i++) {
		const BlobLoc &b = ix.blobs[i];
		if (i > 0) {
			const BlobLoc &prev = ix.blobs[i - 1];
			if (prev.start_id + int64_t(prev.id_count) > b.start_id) {
				// Lookup would prefer idx0 for these rows; no public run surveyed needed that.
				throw SraUnsupportedError("blobs for rows " + RowsText(prev.start_id, prev.id_count) + " and " +
				                          RowsText(b.start_id, b.id_count) + " overlap");
			}
		}
		uint64_t at = b.pg * ix.page_size;
		if (at > ix.data_eof || b.size > ix.data_eof - at) {
			throw SraFormatError("the blob for rows " + RowsText(b.start_id, b.id_count) +
			                     " lies past the column's data");
		}
	}
	return ix;
}

const BlobLoc &ColumnIndex::Locate(int64_t row_id) const {
	auto it = std::upper_bound(blobs.begin(), blobs.end(), row_id,
	                           [](int64_t id, const BlobLoc &b) { return id < b.start_id; });
	if (it != blobs.begin()) {
		const BlobLoc &b = *std::prev(it);
		if (row_id < b.start_id + int64_t(b.id_count)) {
			return b;
		}
	}
	throw SraFormatError("no blob holds row " + std::to_string(row_id));
}

uint32_t Crc32Ncbi(const uint8_t *data, size_t n, uint32_t crc) {
	constexpr uint32_t POLYNOMIAL = 0x04C11DB7;
	constexpr uint32_t TOP_BIT = uint32_t(1) << 31;
	constexpr int TOP_BYTE_SHIFT = 32 - BITS_PER_BYTE;
	static const auto table = [] {
		std::array<uint32_t, UINT8_MAX + 1> t {};
		for (uint32_t i = 0; i < t.size(); i++) {
			uint32_t c = i << TOP_BYTE_SHIFT;
			for (int k = 0; k < BITS_PER_BYTE; k++) {
				c = (c & TOP_BIT) ? (c << 1) ^ POLYNOMIAL : (c << 1);
			}
			t[i] = c;
		}
		return t;
	}();
	for (size_t i = 0; i < n; i++) {
		crc = (crc << BITS_PER_BYTE) ^ table[((crc >> TOP_BYTE_SHIFT) ^ data[i]) & UINT8_MAX];
	}
	return crc;
}

// ---------------------------------------------------------------------------------------------
// Blob framing

namespace {

// VDB's unsigned variable-length integer: 7-bit groups, most significant first, the top bit set
// on every byte but the last. Ten bytes hold any 64-bit value.
constexpr int VLEN_GROUP_BITS = 7;
constexpr uint8_t VLEN_GROUP_MASK = 0x7f;
constexpr uint8_t VLEN_CONTINUES = 0x80;
constexpr size_t VLEN_MAX_BYTES = 10;
// The signed form's first byte holds a sign bit and 6 value bits.
constexpr uint8_t VLEN_SIGN = 0x40;
constexpr uint8_t VLEN_FIRST_MASK = 0x3f;

uint64_t ReadVlenU(const uint8_t *p, size_t len, size_t &pos, const char *what) {
	uint64_t v = 0;
	for (size_t i = 0; i < VLEN_MAX_BYTES; i++) {
		if (pos >= len) {
			throw SraFormatError(std::string(what) + " is truncated");
		}
		uint8_t b = p[pos++];
		v = (v << VLEN_GROUP_BITS) | (b & VLEN_GROUP_MASK);
		if (!(b & VLEN_CONTINUES)) {
			return v;
		}
	}
	throw SraFormatError(std::string(what) + " holds an over-long integer");
}

int64_t ReadVlenS(const uint8_t *p, size_t len, size_t &pos, const char *what) {
	if (pos >= len) {
		throw SraFormatError(std::string(what) + " is truncated");
	}
	uint8_t b = p[pos++];
	bool negative = b & VLEN_SIGN;
	uint64_t v = b & VLEN_FIRST_MASK;
	for (size_t i = 1; (b & VLEN_CONTINUES); i++) {
		if (i == VLEN_MAX_BYTES || pos >= len) {
			throw SraFormatError(std::string(what) + " holds a bad integer");
		}
		b = p[pos++];
		v = (v << VLEN_GROUP_BITS) | (b & VLEN_GROUP_MASK);
	}
	return negative ? -int64_t(v) : int64_t(v);
}

// A blob's encoding headers: a stack version byte, then frames of a flags byte, a version byte and
// vlen fields.
constexpr uint8_t HEADER_STACK_VERSION = 0;
constexpr size_t HEADER_FRAME_FIXED_SIZE = 2;

std::vector<HeaderFrame> ParseHeaders(const uint8_t *p, size_t len) {
	const char *what = "a blob's encoding headers";
	if (p[0] != HEADER_STACK_VERSION) {
		throw SraUnsupportedError("blob header stack version " + std::to_string(p[0]));
	}
	std::vector<HeaderFrame> frames;
	size_t pos = sizeof(HEADER_STACK_VERSION);
	while (pos < len) {
		HeaderFrame f;
		if (len - pos < HEADER_FRAME_FIXED_SIZE) {
			throw SraFormatError(std::string(what) + " are truncated");
		}
		f.flags = p[pos];
		f.version = p[pos + 1];
		pos += HEADER_FRAME_FIXED_SIZE;
		f.fmt = ReadVlenS(p, len, pos, what);
		f.osize = ReadVlenS(p, len, pos, what);
		int64_t op_count = ReadVlenS(p, len, pos, what);
		int64_t arg_count = ReadVlenS(p, len, pos, what);
		if (op_count < 0 || arg_count < 0 || uint64_t(op_count) > len - pos || f.osize < 0) {
			throw SraFormatError(std::string(what) + " are malformed");
		}
		f.ops.assign(p + pos, p + pos + op_count);
		pos += size_t(op_count);
		if (uint64_t(arg_count) > len - pos) {
			throw SraFormatError(std::string(what) + " are malformed");
		}
		f.args.reserve(size_t(arg_count));
		for (int64_t i = 0; i < arg_count; i++) {
			f.args.push_back(ReadVlenS(p, len, pos, what));
		}
		frames.push_back(std::move(f));
	}
	return frames;
}

// A page map tells, per row of a blob, how many elements the row has and where they start in the
// decoded data. Rows of equal length come as (length, run) pairs; consecutive rows sharing one
// stored record as data runs; version 2 instead records each row's offset ("random access").
// A page map's first byte: the variant in its low two bits, the version above them. Version 0
// stores the values plainly, 1 deflated, and 2 deflated with an offset per row.
constexpr uint8_t PAGE_MAP_VARIANT_MASK = 3;
constexpr int PAGE_MAP_VERSION_SHIFT = 2;
constexpr uint8_t PAGE_MAP_DEFLATED = 1;
constexpr uint8_t PAGE_MAP_RANDOM_ACCESS = 2;
// Which runs the map stores: one row length for every row or length runs, with or without data
// runs (records shared by consecutive rows).
enum PageMapVariant : uint8_t { ONE_LENGTH = 0, ONE_LENGTH_DATA_RUNS = 1, LENGTH_RUNS = 2, LENGTH_AND_DATA_RUNS = 3 };
// Far past any blob's: VDB cuts blobs by size, so even empty rows run to millions, not billions.
constexpr uint64_t MAX_PAGE_MAP_VALUES = uint64_t(1) << 32;
// Inflating stops at the space the claimed values could need, plus this much for the stream's end.
constexpr uint64_t PAGE_MAP_INFLATE_SLACK = 16;

PageMap ParsePageMap(const uint8_t *src, size_t len, uint64_t row_count) {
	const char *what = "a page map";
	uint8_t variant = src[0] & PAGE_MAP_VARIANT_MASK;
	uint8_t version = src[0] >> PAGE_MAP_VERSION_SHIFT;
	if (version > PAGE_MAP_RANDOM_ACCESS) {
		throw SraFormatError("page map version " + std::to_string(version));
	}
	const bool random_access = version == PAGE_MAP_RANDOM_ACCESS;
	const bool deflated = version >= PAGE_MAP_DEFLATED;
	const bool has_length_runs = variant == LENGTH_RUNS || variant == LENGTH_AND_DATA_RUNS;
	const bool has_data_runs = variant == ONE_LENGTH_DATA_RUNS || variant == LENGTH_AND_DATA_RUNS;
	size_t pos = 1;
	uint64_t row_len = 0, data_recs = 0, leng_recs = 0;
	if (has_length_runs) {
		leng_recs = ReadVlenU(src, len, pos, what);
	} else {
		row_len = ReadVlenU(src, len, pos, what);
	}
	if (has_data_runs) {
		data_recs = ReadVlenU(src, len, pos, what);
	}
	// Each count is bounded before they are summed, so the sum cannot wrap past the bound.
	// The row count adds to them only for a random-access map, which stores an offset per row.
	if (leng_recs > MAX_PAGE_MAP_VALUES || data_recs > MAX_PAGE_MAP_VALUES ||
	    (random_access && row_count > MAX_PAGE_MAP_VALUES)) {
		throw SraFormatError("a page map declares too many values");
	}
	// A length run is two values (the length and how many rows have it).
	uint64_t values = (has_length_runs ? 2 * leng_recs : 0) + (has_data_runs ? data_recs : 0) +
	                  (random_access && !has_data_runs ? row_count : 0);
	if (values > MAX_PAGE_MAP_VALUES) {
		throw SraFormatError("a page map declares too many values");
	}
	const uint8_t *body = src + pos;
	size_t body_len = len - pos;
	// Every value is a vlen of at least one byte, and a deflated body expands at most
	// MAX_DEFLATE_RATIO times: a count the bytes cannot hold is refused before anything is sized
	// from it.
	const uint64_t max_values = deflated ? uint64_t(body_len) * MAX_DEFLATE_RATIO : uint64_t(body_len);
	if (values > max_values) {
		throw SraFormatError("a page map claims " + std::to_string(values) + " values in " + std::to_string(body_len) +
		                     " bytes");
	}
	std::vector<uint8_t> inflated;
	if (deflated && values > 0) {
		// The stream never needs more than the longest vlen per value; the buffer grows with what
		// it actually produces rather than being sized from the claim.
		inflated = InflateRawGrowing(
		    body, body_len, std::min<uint64_t>(values * VLEN_MAX_BYTES + PAGE_MAP_INFLATE_SLACK, MAX_DECODED_BLOB),
		    "a page map");
		body = inflated.data();
		body_len = inflated.size();
	}
	size_t bpos = 0;
	auto take = [&](uint64_t n) {
		if (n > body_len - bpos) {
			throw SraFormatError("a page map's values run past its end");
		}
		std::vector<uint64_t> v(n);
		for (auto &x : v) {
			x = ReadVlenU(body, body_len, bpos, what);
		}
		return v;
	};
	auto length_of = [](uint64_t length) {
		if (length > UINT32_MAX) {
			throw SraFormatError("a page map row length overflows");
		}
		return uint32_t(length);
	};
	const std::string not_covered = "a page map's runs do not cover the blob's " + std::to_string(row_count) + " rows";

	PageMap pm;
	pm.row_count = row_count;
	uint64_t covered = 0;
	auto add = [&](uint64_t rows, uint32_t length, uint64_t offset, uint64_t step) {
		if (rows > row_count - covered) {
			throw SraFormatError(not_covered);
		}
		if (rows > 0) {
			pm.runs.push_back({covered, rows, length, offset, step});
			covered += rows;
		}
	};
	// The next record's first element, refused rather than wrapped.
	auto advance = [](uint64_t &offset, uint64_t length, uint64_t rows) {
		if (length != 0 && rows > (UINT64_MAX - offset) / length) {
			throw SraFormatError("a page map's records run past 2^64 elements");
		}
		offset += length * rows;
	};

	std::vector<uint64_t> length, leng_run;
	if (has_length_runs) {
		length = take(leng_recs);
		leng_run = take(leng_recs);
	}
	std::vector<uint64_t> data_run;
	if (has_data_runs) {
		data_run = take(data_recs);
	}

	if (random_access && !has_data_runs) {
		// One stored offset per row, kept as stored; the runs give only the lengths.
		pm.row_offset = take(row_count);
		if (variant == ONE_LENGTH) {
			add(row_count, length_of(row_len), 0, 0);
		} else {
			for (uint64_t i = 0; i < leng_recs; i++) {
				add(leng_run[i], length_of(length[i]), 0, 0);
			}
		}
	} else if (variant == ONE_LENGTH) {
		add(row_count, length_of(row_len), 0, row_len);
	} else if (variant == LENGTH_RUNS) {
		uint64_t offset = 0;
		for (uint64_t i = 0; i < leng_recs; i++) {
			add(leng_run[i], length_of(length[i]), offset, length[i]);
			advance(offset, length[i], leng_run[i]);
		}
	} else {
		// Each data run is one stored record shared by consecutive rows of equal length.
		uint64_t offset = 0;
		size_t li = 0;
		uint64_t left = variant == LENGTH_AND_DATA_RUNS && leng_recs > 0 ? leng_run[0] : 0;
		for (uint64_t run : data_run) {
			if (run == 0) {
				throw SraFormatError(not_covered);
			}
			uint32_t len = length_of(row_len);
			if (variant == LENGTH_AND_DATA_RUNS) {
				// The rows the record covers must all have one length, across length runs.
				uint64_t need = run;
				bool first = true;
				while (need > 0) {
					while (left == 0) {
						if (++li >= leng_recs) {
							throw SraFormatError(not_covered);
						}
						left = leng_run[li];
					}
					uint32_t here = length_of(length[li]);
					if (first) {
						len = here;
						first = false;
					} else if (here != len) {
						throw SraFormatError("a page map shares one record among rows of different lengths");
					}
					uint64_t used = std::min(need, left);
					need -= used;
					left -= used;
				}
			}
			add(run, len, offset, 0);
			advance(offset, len, 1);
		}
	}
	if (covered != row_count) {
		throw SraFormatError(not_covered);
	}
	if (has_length_runs) {
		uint64_t total = 0;
		for (uint64_t r : leng_run) {
			if (r > row_count - total) {
				throw SraFormatError(not_covered);
			}
			total += r;
		}
		if (total != row_count) {
			throw SraFormatError(not_covered);
		}
	}
	if (bpos != body_len) {
		throw SraFormatError("a page map has " + std::to_string(body_len - bpos) + " unexpected trailing bytes");
	}
	return pm;
}

} // namespace

uint64_t PageMap::ElementsSpanned() const {
	uint64_t end = 0;
	if (!row_offset.empty()) {
		for (const Run &r : runs) {
			for (uint64_t i = r.first_row; i < r.first_row + r.rows; i++) {
				if (r.length > UINT64_MAX - row_offset[i]) {
					return UINT64_MAX;
				}
				end = std::max(end, row_offset[i] + r.length);
			}
		}
		return end;
	}
	for (const Run &r : runs) {
		// offset + step * (rows - 1) + length, saturating: a wrapped sum must not pass a bound.
		uint64_t last = r.rows - 1;
		if (r.step != 0 && last > (UINT64_MAX - r.offset) / r.step) {
			return UINT64_MAX;
		}
		uint64_t start = r.offset + r.step * last;
		if (r.length > UINT64_MAX - start) {
			return UINT64_MAX;
		}
		end = std::max(end, start + r.length);
	}
	return end;
}

void PageMap::Row(uint64_t i, uint32_t &length, uint64_t &offset, size_t &hint) const {
	if (i >= row_count) {
		throw SraFormatError("a blob has no row " + std::to_string(i));
	}
	auto holds = [&](size_t k) {
		return k < runs.size() && runs[k].first_row <= i && i - runs[k].first_row < runs[k].rows;
	};
	if (!holds(hint)) {
		if (holds(hint + 1)) {
			hint++;
		} else {
			auto it = std::upper_bound(runs.begin(), runs.end(), i,
			                           [](uint64_t row, const Run &r) { return row < r.first_row; });
			hint = size_t(std::prev(it) - runs.begin());
		}
	}
	const Run &r = runs[hint];
	length = r.length;
	if (!row_offset.empty()) {
		offset = row_offset[i];
		return;
	}
	offset = r.offset + r.step * (i - r.first_row);
}

namespace {
// A blob's first byte. Its top bit tells a v2 blob from a v1 one. A v1 byte holds the byte order,
// the padding bits at the payload's end, and which of four widths its row length is stored in.
// A v2 byte holds the top two bits 10, which of four layouts its header sizes take, the byte order,
// and the payload's significant bits in its last byte.
constexpr uint8_t BLOB_V2 = 0x80;
constexpr uint8_t BLOB_V1_ORDER_MASK = 3;
constexpr uint8_t BLOB_V1_BIG_ENDIAN = 2;
constexpr int BLOB_V1_ADJUST_SHIFT = 2;
constexpr uint8_t BLOB_V1_ADJUST_MASK = 7;
constexpr int BLOB_V1_RLS_SHIFT = 5;
constexpr uint8_t BLOB_V1_RLS_MASK = 3;
// The bytes a v1 row length takes, by its 2-bit code; 0 means rows of one element.
constexpr uint8_t BLOB_V1_ROW_LEN_BYTES[] = {1, 2, 4, 0};
constexpr int BLOB_V2_TAG_SHIFT = 6;
constexpr uint8_t BLOB_V2_TAG = 2;
constexpr uint8_t BLOB_V2_LAST_BITS_MASK = 7;
constexpr int BLOB_V2_ORDER_SHIFT = 3;
constexpr int BLOB_V2_LAYOUT_SHIFT = 4;
constexpr uint8_t BLOB_V2_LAYOUT_MASK = 3;
// The first byte plus the header and page-map sizes: one byte each, one and two, one and four, or
// four each.
constexpr size_t BLOB_V2_PREFIX_SIZE[] = {3, 4, 6, 9};
} // namespace

VdbBlob ParseVdbBlob(std::vector<uint8_t> bytes, uint64_t row_count, uint32_t elem_bits) {
	if (bytes.empty()) {
		throw SraFormatError("an empty blob");
	}
	VdbBlob blob;
	const uint8_t b0 = bytes[0];
	if ((b0 & BLOB_V2) == 0) {
		// v1: fixed row length, no encoding headers and no page map. This branch follows ncbi-vdb's
		// VBlobCreateFromData_v1 (libs/vdb/blob.c; public domain, THIRD_PARTY_LICENSES.md), its
		// row-length repair included.
		uint8_t byte_order = b0 & BLOB_V1_ORDER_MASK;
		uint8_t adjust = (b0 >> BLOB_V1_ADJUST_SHIFT) & BLOB_V1_ADJUST_MASK;
		uint8_t rls = BLOB_V1_ROW_LEN_BYTES[(b0 >> BLOB_V1_RLS_SHIFT) & BLOB_V1_RLS_MASK];
		size_t offset = size_t(rls) + 1;
		if (bytes.size() < offset) {
			throw SraFormatError("a v1 blob is truncated");
		}
		uint64_t row_len = 1;
		if (rls != 0) {
			row_len = 0;
			for (size_t i = 0; i < rls; i++) {
				row_len |= uint64_t(bytes[1 + i]) << (BITS_PER_BYTE * i);
			}
		}
		uint64_t payload_bits = uint64_t(bytes.size() - offset) * BITS_PER_BYTE;
		if (payload_bits < adjust) {
			throw SraFormatError("a v1 blob is truncated");
		}
		payload_bits -= adjust;
		if (payload_bits % elem_bits != 0) {
			throw SraFormatError("a v1 blob's payload is not a whole number of elements");
		}
		// A 4-byte row length that disagrees with the payload is corrected from the row count,
		// as VDB does.
		if (rls == sizeof(uint32_t) && row_count != 0 && row_len * row_count != payload_bits / elem_bits) {
			row_len = (payload_bits / elem_bits) / row_count;
		}
		if (row_len == 0 || row_len > UINT32_MAX) {
			throw SraUnsupportedError("a legacy v1 blob without a row length");
		}
		blob.big_endian = byte_order == BLOB_V1_BIG_ENDIAN;
		blob.page_map.row_count = row_count;
		if (row_count > 0) {
			blob.page_map.runs.push_back({0, row_count, uint32_t(row_len), 0, row_len});
		}
		blob.payload_offset = offset;
		blob.payload_bits = payload_bits;
		blob.bytes = std::move(bytes);
		return blob;
	}

	// v2: header sizes, encoding headers, page map, payload.
	if ((b0 >> BLOB_V2_TAG_SHIFT) != BLOB_V2_TAG) {
		throw SraFormatError("a blob has an unknown header byte");
	}
	// The last byte's significant bits, 0 meaning all eight: the rest are padding.
	uint8_t adjust = (BITS_PER_BYTE - (b0 & BLOB_V2_LAST_BITS_MASK)) & BLOB_V2_LAST_BITS_MASK;
	blob.big_endian = (b0 >> BLOB_V2_ORDER_SHIFT) & 1;
	uint8_t layout = (b0 >> BLOB_V2_LAYOUT_SHIFT) & BLOB_V2_LAYOUT_MASK;
	size_t offset = BLOB_V2_PREFIX_SIZE[layout];
	if (bytes.size() < offset) {
		throw SraFormatError("a v2 blob is truncated");
	}
	// Both sizes follow the first byte, the header size first.
	const uint8_t *sizes = &bytes[1];
	uint64_t hdr_size, map_size;
	switch (layout) {
	case 0:
		hdr_size = sizes[0];
		map_size = sizes[1];
		break;
	case 1:
		hdr_size = sizes[0];
		map_size = LoadU16(sizes + 1);
		break;
	case 2:
		hdr_size = sizes[0];
		map_size = LoadU32(sizes + 1);
		break;
	default:
		hdr_size = LoadU32(sizes);
		map_size = LoadU32(sizes + sizeof(uint32_t));
		break;
	}
	if (hdr_size > bytes.size() - offset || map_size > bytes.size() - offset - hdr_size) {
		throw SraFormatError("a v2 blob's headers run past its end");
	}
	if (hdr_size > 0) {
		blob.frames = ParseHeaders(bytes.data() + offset, size_t(hdr_size));
	}
	if (map_size == 0) {
		throw SraUnsupportedError("a v2 blob without a page map, which no public run surveyed had");
	}
	blob.page_map = ParsePageMap(bytes.data() + offset + hdr_size, size_t(map_size), row_count);
	blob.payload_offset = size_t(offset + hdr_size + map_size);
	uint64_t payload_bits = uint64_t(bytes.size() - blob.payload_offset) * BITS_PER_BYTE;
	if (payload_bits < adjust) {
		throw SraFormatError("a v2 blob's payload is truncated");
	}
	blob.payload_bits = payload_bits - adjust;
	blob.bytes = std::move(bytes);
	return blob;
}

ElementType ElementTypeOf(const std::string &type) {
	static const std::map<std::string, ElementType> TYPES = {{"INSDC:2na:packed", {2, false}},
	                                                         {"INSDC:4na:packed", {4, false}},
	                                                         {"INSDC:4na:bin", {8, false}},
	                                                         {"INSDC:dna:text", {8, false}},
	                                                         {"INSDC:quality:phred", {8, false}},
	                                                         {"INSDC:coord:len", {32, false}},
	                                                         {"INSDC:coord:zero", {32, true}},
	                                                         {"INSDC:coord:one", {32, true}},
	                                                         {"INSDC:coord:val", {32, true}},
	                                                         {"INSDC:SRA:xread_type", {8, false}},
	                                                         {"INSDC:SRA:read_type", {8, false}},
	                                                         {"INSDC:SRA:spot_filter", {8, false}},
	                                                         {"INSDC:SRA:read_filter", {8, false}},
	                                                         {"INSDC:SRA:platform_id", {8, false}},
	                                                         {"ascii", {8, false}},
	                                                         {"U8", {8, false}},
	                                                         {"U16", {16, false}},
	                                                         {"U32", {32, false}},
	                                                         {"U64", {64, false}},
	                                                         {"I8", {8, true}},
	                                                         {"I16", {16, true}},
	                                                         {"I32", {32, true}},
	                                                         {"I64", {64, true}}};
	auto it = TYPES.find(type);
	if (it == TYPES.end()) {
		throw SraUnsupportedError("stored type '" + type + "'");
	}
	return it->second;
}

uint32_t ElementBits(const std::string &type) {
	return ElementTypeOf(type).bits;
}

// ---------------------------------------------------------------------------------------------
// Physical column

PhysicalColumn::PhysicalColumn(const SraTable &table, const std::string &name)
    : archive_(&table.Archive()), name_(name), dir_(table.ColumnDir(name)) {
	const KarArchive &ar = *archive_;
	if (!ar.IsDirectory(dir_)) {
		throw SraFormatError("the table has no column '" + name + "'");
	}
	const std::string md_path = dir_ + "/" + KDB_METADATA;
	if (ar.Find(md_path) == nullptr) {
		throw SraUnsupportedError("column '" + name + "' has no metadata (a legacy column)");
	}
	MetadataNode md = ParseMetadata(ar.ReadFile(md_path));
	const MetadataNode *schema = md.Find("schema");
	if (schema == nullptr) {
		throw SraUnsupportedError("column '" + name + "' has no schema node (a legacy v0/v1 column)");
	}
	auto type = schema->attrs.find("type");
	if (type == schema->attrs.end()) {
		throw SraFormatError("column '" + name + "' does not record its type");
	}
	type_ = type->second;
	auto expr = schema->attrs.find("expr");
	if (expr != schema->attrs.end()) {
		encoding_ = expr->second;
	}
	auto file = [&](const char *f) {
		std::string path = dir_ + "/" + f;
		return ar.Find(path) == nullptr ? std::vector<uint8_t>() : ar.ReadFile(path);
	};
	index_ = ColumnIndex::Parse(file("idx"), file("idx0"), file("idx1"), file("idx2"));
}

std::vector<uint8_t> PhysicalColumn::ReadBlob(const BlobLoc &loc) const {
	size_t cs_bytes = 0;
	switch (index_.checksum) {
	case ChecksumType::None:
		break;
	case ChecksumType::Crc32:
		cs_bytes = sizeof(uint32_t);
		break;
	case ChecksumType::Md5:
		throw SraUnsupportedError("column '" + name_ + "' uses MD5 blob checksums, which no public run surveyed did");
	}
	if (loc.size < cs_bytes) {
		throw SraFormatError("a blob of column '" + name_ + "' is shorter than its checksum");
	}
	auto bytes = archive_->ReadFile(dir_ + "/data", loc.pg * index_.page_size, loc.size);
	if (index_.checksum == ChecksumType::Crc32) {
		uint32_t stored = LoadU32(bytes.data() + bytes.size() - cs_bytes);
		bytes.resize(bytes.size() - cs_bytes);
		if (Crc32Ncbi(bytes.data(), bytes.size()) != stored) {
			throw SraFormatError("checksum mismatch in column '" + name_ + "', rows " +
			                     RowsText(loc.start_id, loc.id_count));
		}
	}
	return bytes;
}

void PhysicalColumn::PrefetchRows(int64_t first, int64_t end) const {
	const KarEntry *data = archive_->Find(dir_ + "/data");
	if (data == nullptr || first >= end) {
		return;
	}
	const auto &blobs = index_.blobs;
	auto it = std::upper_bound(blobs.begin(), blobs.end(), first,
	                           [](int64_t id, const BlobLoc &b) { return id < b.start_id; });
	if (it != blobs.begin() && std::prev(it)->start_id + int64_t(std::prev(it)->id_count) > first) {
		--it;
	}
	// Blobs of consecutive rows usually lie end to end: hint them as one range.
	uint64_t run_start = 0, run_end = 0;
	for (; it != blobs.end() && it->start_id < end; ++it) {
		uint64_t at = it->pg * index_.page_size;
		if (run_end != run_start && at != run_end) {
			archive_->Source().Prefetch(data->offset + run_start, run_end - run_start);
			run_start = run_end = at;
		}
		if (run_end == run_start) {
			run_start = run_end = at;
		}
		run_end = at + it->size;
	}
	if (run_end > run_start) {
		archive_->Source().Prefetch(data->offset + run_start, run_end - run_start);
	}
}

VdbBlob PhysicalColumn::ReadVdbBlob(const BlobLoc &loc) const {
	return ParseVdbBlob(ReadBlob(loc), loc.id_count, ElementBits(type_));
}

} // namespace sra
} // namespace miint
