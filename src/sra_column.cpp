#include "sra_column.hpp"

#include "sra_codecs.hpp"

#include <algorithm>
#include <array>
#include <map>

namespace miint {
namespace sra {

namespace {

constexpr uint32_t ORDER_NATIVE = 0x05031988;
constexpr uint32_t ORDER_SWAPPED = 0x88190305;
constexpr size_t LOC_SIZE = 24; // a KColBlobLoc / KColBlockLoc

void CheckKdbHeader(const std::vector<uint8_t> &file, const std::string &what) {
	if (file.size() < 8) {
		throw SraFormatError(what + " is shorter than its header");
	}
	uint32_t order = LoadU32(file.data());
	if (order == ORDER_SWAPPED) {
		throw SraUnsupportedError(what + " is big-endian, which no public run surveyed was");
	}
	if (order != ORDER_NATIVE) {
		throw SraFormatError(what + " has an invalid byte-order tag");
	}
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
	uint32_t gen = LoadU32(idx1.data() + at + 8);
	uint32_t id_range = LoadU32(idx1.data() + at + 12);
	int64_t start_id = int64_t(LoadU64(idx1.data() + at + 16));
	uint32_t size = gen & 0x7ffffff;
	uint32_t id_type = (gen >> 27) & 3;
	uint32_t pg_type = (gen >> 29) & 3;
	if (gen >> 31) {
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
		blk_size = 12;
	} else {
		blk_size = size;
	}
	if (pg > idx2_eof || blk_size > idx2_eof - pg || idx2_eof > idx2.size()) {
		throw SraFormatError("an idx2 block lies outside idx2");
	}
	blk = idx2.data() + pg;

	size_t id_hdr = id_type == MAGNITUDE ? 8 : 0;
	size_t pg_hdr = pg_type == MAGNITUDE ? 8 : 12;
	size_t id_each = id_type == MAGNITUDE ? 4 : 0;
	size_t pg_each = pg_type == MAGNITUDE ? 4 : 0;
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
			blobs[i].id_count = LoadU32(id_s + 4 * i);
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
			blobs[i].size = LoadU32(pg_s + 4 * i);
			next += blobs[i].size;
		}
	} else {
		uint64_t first = LoadU64(blk + id_hdr);
		uint32_t each = LoadU32(blk + id_hdr + 8);
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
	uint32_t version = LoadU32(idx.data() + 4);
	if (version != 3) {
		throw SraUnsupportedError("column header version " + std::to_string(version) +
		                          " is not supported (public runs use version 3)");
	}
	if (idx.size() < 40) {
		throw SraFormatError("a column header is truncated");
	}
	ColumnIndex ix;
	ix.data_eof = LoadU64(idx.data() + 8);
	uint64_t idx2_eof = LoadU64(idx.data() + 16);
	uint32_t idx0_count = LoadU32(idx.data() + 24);
	uint32_t num_blocks = LoadU32(idx.data() + 28);
	ix.page_size = LoadU32(idx.data() + 32);
	switch (idx[36]) {
	case 0:
		ix.checksum = ChecksumType::None;
		break;
	case 1:
		ix.checksum = ChecksumType::Crc32;
		break;
	case 2:
		ix.checksum = ChecksumType::Md5;
		break;
	default:
		throw SraFormatError("a column header names checksum type " + std::to_string(idx[36]));
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
		uint32_t gen = LoadU32(p + 8);
		loc.size = gen & 0x7fffffff;
		loc.id_count = LoadU32(p + 12);
		loc.start_id = int64_t(LoadU64(p + 16));
		if (gen >> 31) {
			throw SraUnsupportedError("an idx0 entry removing a blob");
		}
		if (loc.id_count == 0) {
			throw SraFormatError("an idx0 entry covers no rows");
		}
		// Entries never overlap one another unless they cover exactly the same rows.
		int64_t end = loc.start_id + int64_t(loc.id_count);
		auto next = journal.lower_bound(loc.start_id);
		std::vector<const BlobLoc *> neighbours;
		if (next != journal.end()) {
			neighbours.push_back(&next->second);
		}
		if (next != journal.begin()) {
			neighbours.push_back(&std::prev(next)->second);
		}
		for (const BlobLoc *e : neighbours) {
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
		if ((idx1.size() - 8) / LOC_SIZE < num_blocks) {
			throw SraFormatError("idx1 holds fewer blocks than the column header records");
		}
	}
	std::vector<BlobLoc> from_blocks;
	for (size_t b = 0; b < num_blocks; b++) {
		DecodeBlock(idx1, 8 + b * LOC_SIZE, idx2, idx2_eof, from_blocks);
	}

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
	static const auto table = [] {
		std::array<uint32_t, 256> t {};
		for (uint32_t i = 0; i < 256; i++) {
			uint32_t c = i << 24;
			for (int k = 0; k < 8; k++) {
				c = (c & 0x80000000u) ? (c << 1) ^ 0x04C11DB7u : (c << 1);
			}
			t[i] = c;
		}
		return t;
	}();
	for (size_t i = 0; i < n; i++) {
		crc = (crc << 8) ^ table[((crc >> 24) ^ data[i]) & 0xff];
	}
	return crc;
}

// ---------------------------------------------------------------------------------------------
// Blob framing

namespace {

// VDB's unsigned variable-length integer: 7-bit groups, most significant first, the top bit set
// on every byte but the last.
uint64_t ReadVlenU(const uint8_t *p, size_t len, size_t &pos, const char *what) {
	uint64_t v = 0;
	for (int i = 0; i < 10; i++) {
		if (pos >= len) {
			throw SraFormatError(std::string(what) + " is truncated");
		}
		uint8_t b = p[pos++];
		v = (v << 7) | (b & 0x7f);
		if (!(b & 0x80)) {
			return v;
		}
	}
	throw SraFormatError(std::string(what) + " holds an over-long integer");
}

// The signed form: the first byte carries a sign bit (0x40) and 6 value bits.
int64_t ReadVlenS(const uint8_t *p, size_t len, size_t &pos, const char *what) {
	if (pos >= len) {
		throw SraFormatError(std::string(what) + " is truncated");
	}
	uint8_t b = p[pos++];
	bool negative = b & 0x40;
	uint64_t v = b & 0x3f;
	for (int i = 0; (b & 0x80); i++) {
		if (i == 9 || pos >= len) {
			throw SraFormatError(std::string(what) + " holds a bad integer");
		}
		b = p[pos++];
		v = (v << 7) | (b & 0x7f);
	}
	return negative ? -int64_t(v) : int64_t(v);
}

std::vector<HeaderFrame> ParseHeaders(const uint8_t *p, size_t len) {
	const char *what = "a blob's encoding headers";
	if (p[0] != 0) {
		throw SraUnsupportedError("blob header stack version " + std::to_string(p[0]));
	}
	std::vector<HeaderFrame> frames;
	size_t pos = 1;
	while (pos < len) {
		HeaderFrame f;
		if (len - pos < 2) {
			throw SraFormatError(std::string(what) + " are truncated");
		}
		f.flags = p[pos];
		f.version = p[pos + 1];
		pos += 2;
		f.fmt = ReadVlenS(p, len, pos, what);
		f.osize = ReadVlenS(p, len, pos, what);
		int64_t op_count = ReadVlenS(p, len, pos, what);
		int64_t arg_count = ReadVlenS(p, len, pos, what);
		if (op_count < 0 || arg_count < 0 || uint64_t(op_count) > len - pos || f.osize < 0) {
			throw SraFormatError(std::string(what) + " are malformed");
		}
		f.ops.assign(p + pos, p + pos + op_count);
		pos += size_t(op_count);
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
PageMap ParsePageMap(const uint8_t *src, size_t len, uint64_t row_count) {
	const char *what = "a page map";
	uint8_t variant = src[0] & 3;
	uint8_t version = src[0] >> 2;
	if (version > 2) {
		throw SraFormatError("page map version " + std::to_string(version));
	}
	bool random_access = version == 2;
	size_t pos = 1;
	uint64_t row_len = 0, data_recs = 0, leng_recs = 0;
	switch (variant) {
	case 0:
		row_len = ReadVlenU(src, len, pos, what);
		break;
	case 1:
		row_len = ReadVlenU(src, len, pos, what);
		data_recs = ReadVlenU(src, len, pos, what);
		break;
	case 2:
		leng_recs = ReadVlenU(src, len, pos, what);
		break;
	default:
		leng_recs = ReadVlenU(src, len, pos, what);
		data_recs = ReadVlenU(src, len, pos, what);
		break;
	}
	uint64_t values = (variant >= 2 ? 2 * leng_recs : 0) + (variant == 1 || variant == 3 ? data_recs : 0) +
	                  (random_access && variant != 1 && variant != 3 ? row_count : 0);
	if (values > (uint64_t(1) << 32)) {
		throw SraFormatError("a page map declares too many values");
	}
	const uint8_t *body = src + pos;
	size_t body_len = len - pos;
	// Every value is a vlen of at least one byte, and a deflated body expands at most
	// MAX_DEFLATE_RATIO times: a count the bytes cannot hold is refused before anything is sized
	// from it.
	const uint64_t max_values = version >= 1 ? uint64_t(body_len) * MAX_DEFLATE_RATIO : uint64_t(body_len);
	if (values > max_values) {
		throw SraFormatError("a page map claims " + std::to_string(values) + " values in " + std::to_string(body_len) +
		                     " bytes");
	}
	std::vector<uint8_t> inflated;
	if (version >= 1 && values > 0) {
		inflated.resize(size_t(values) * 10 + 16);
		inflated.resize(InflateRaw(body, body_len, inflated.data(), inflated.size(), "a page map").produced);
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
	if (variant >= 2) {
		length = take(leng_recs);
		leng_run = take(leng_recs);
	}
	std::vector<uint64_t> data_run;
	if (variant == 1 || variant == 3) {
		data_run = take(data_recs);
	}

	if (random_access && (variant == 0 || variant == 2)) {
		// One stored offset per row: as many runs as offsets, so bounded by the map's bytes.
		auto offsets = take(row_count);
		uint64_t row = 0;
		auto per_row = [&](uint64_t rows, uint32_t len) {
			if (rows > row_count - row) {
				throw SraFormatError(not_covered);
			}
			for (uint64_t r = 0; r < rows; r++, row++) {
				add(1, len, offsets[row], len);
			}
		};
		if (variant == 0) {
			per_row(row_count, length_of(row_len));
		} else {
			for (uint64_t i = 0; i < leng_recs; i++) {
				per_row(leng_run[i], length_of(length[i]));
			}
		}
	} else if (variant == 0) {
		add(row_count, length_of(row_len), 0, row_len);
	} else if (variant == 2) {
		uint64_t offset = 0;
		for (uint64_t i = 0; i < leng_recs; i++) {
			add(leng_run[i], length_of(length[i]), offset, length[i]);
			advance(offset, length[i], leng_run[i]);
		}
	} else {
		// Each data run is one stored record shared by consecutive rows of equal length.
		uint64_t offset = 0;
		size_t li = 0;
		uint64_t left = variant == 3 && leng_recs > 0 ? leng_run[0] : 0;
		for (uint64_t run : data_run) {
			if (run == 0) {
				throw SraFormatError(not_covered);
			}
			uint32_t len = length_of(row_len);
			if (variant == 3) {
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
	if (variant >= 2) {
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
	offset = r.offset + r.step * (i - r.first_row);
}

VdbBlob ParseVdbBlob(std::vector<uint8_t> bytes, uint64_t row_count, uint32_t elem_bits) {
	if (bytes.empty()) {
		throw SraFormatError("an empty blob");
	}
	VdbBlob blob;
	const uint8_t b0 = bytes[0];
	if ((b0 & 0x80) == 0) {
		// v1: fixed row length, no encoding headers and no page map.
		uint8_t byte_order = b0 & 3;
		uint8_t adjust = (b0 >> 2) & 7;
		static const uint8_t RLS[4] = {1, 2, 4, 0};
		uint8_t rls = RLS[(b0 >> 5) & 3];
		size_t offset = size_t(rls) + 1;
		if (bytes.size() < offset) {
			throw SraFormatError("a v1 blob is truncated");
		}
		uint64_t row_len = 1;
		if (rls != 0) {
			row_len = 0;
			for (size_t i = 0; i < rls; i++) {
				row_len |= uint64_t(bytes[1 + i]) << (8 * i);
			}
		}
		uint64_t payload_bits = uint64_t(bytes.size() - offset) * 8;
		if (payload_bits < adjust) {
			throw SraFormatError("a v1 blob is truncated");
		}
		payload_bits -= adjust;
		if (payload_bits % elem_bits != 0) {
			throw SraFormatError("a v1 blob's payload is not a whole number of elements");
		}
		// A 4-byte row length that disagrees with the payload is corrected from the row count,
		// as VDB does.
		if (rls == 4 && row_count != 0 && row_len * row_count != payload_bits / elem_bits) {
			row_len = (payload_bits / elem_bits) / row_count;
		}
		if (row_len == 0 || row_len > UINT32_MAX) {
			throw SraUnsupportedError("a legacy v1 blob without a row length");
		}
		blob.big_endian = byte_order == 2;
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
	if ((b0 >> 6) != 2) {
		throw SraFormatError("a blob has an unknown header byte");
	}
	uint8_t adjust = (8 - (b0 & 7)) & 7;
	blob.big_endian = (b0 >> 3) & 1;
	uint8_t variant = (b0 >> 4) & 3;
	static const size_t OFFSET[4] = {3, 4, 6, 9};
	size_t offset = OFFSET[variant];
	if (bytes.size() < offset) {
		throw SraFormatError("a v2 blob is truncated");
	}
	uint64_t hdr_size, map_size;
	switch (variant) {
	case 0:
		hdr_size = bytes[1];
		map_size = bytes[2];
		break;
	case 1:
		hdr_size = bytes[1];
		map_size = LoadU16(&bytes[2]);
		break;
	case 2:
		hdr_size = bytes[1];
		map_size = LoadU32(&bytes[2]);
		break;
	default:
		hdr_size = LoadU32(&bytes[1]);
		map_size = LoadU32(&bytes[5]);
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
	uint64_t payload_bits = uint64_t(bytes.size() - blob.payload_offset) * 8;
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
    : archive_(&table.Archive()), name_(name), dir_(table.Prefix() + "col/" + name) {
	const KarArchive &ar = *archive_;
	if (!ar.IsDirectory(dir_)) {
		throw SraFormatError("the table has no column '" + name + "'");
	}
	if (ar.Find(dir_ + "/md/cur") == nullptr) {
		throw SraUnsupportedError("column '" + name + "' has no metadata (a legacy column)");
	}
	MetadataNode md = ParseMetadata(ar.ReadFile(dir_ + "/md/cur"));
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
		cs_bytes = 4;
		break;
	case ChecksumType::Md5:
		throw SraUnsupportedError("column '" + name_ + "' uses MD5 blob checksums, which no public run surveyed did");
	}
	if (loc.size < cs_bytes) {
		throw SraFormatError("a blob of column '" + name_ + "' is shorter than its checksum");
	}
	auto bytes = archive_->ReadFile(dir_ + "/data", loc.pg * index_.page_size, loc.size);
	if (cs_bytes == 4) {
		uint32_t stored = LoadU32(bytes.data() + bytes.size() - 4);
		bytes.resize(bytes.size() - 4);
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
