#include "sra_name.hpp"

#include "sra_archive.hpp"
#include "sra_source.hpp"

#include <algorithm>
#include <cctype>

namespace miint {
namespace sra {

namespace {

// The index file's header: the KDB file header, then the index type (0 for text), the first and
// last row ids, and the bit width of the packed row spans. The trie follows it.
constexpr uint32_t SKEY_VERSION = 4;
constexpr uint32_t SKEY_TEXT_INDEX = 0;
constexpr size_t SKEY_TYPE_AT = 8;
constexpr size_t SKEY_FIRST_AT = 16;
constexpr size_t SKEY_LAST_AT = 24;
constexpr size_t SKEY_SPAN_BITS_AT = 34;
constexpr size_t SKEY_HEADER_SIZE = 40;

// The persisted trie's header: the transition and value-node counts, the node data size (32 low
// bits, then 8 high bits after the key flags), the key flags, and the character map's width.
constexpr size_t TRIE_NUM_TRANS_AT = 0;
constexpr size_t TRIE_NUM_NODES_AT = 4;
constexpr size_t TRIE_DATA_SIZE_AT = 8;
constexpr size_t TRIE_KEY_FLAGS_AT = 12;
constexpr size_t TRIE_DATA_SIZE_HIGH_AT = 13;
constexpr size_t TRIE_WIDTH_AT = 14;
constexpr size_t TRIE_HEADER_SIZE = 16;
constexpr int TRIE_DATA_SIZE_HIGH_SHIFT = 32;
// The key flags: keys stored outside the trie, backtrace links to each node's parent, and how
// node ids split into a transition and a value (3 bits, 7 unused).
constexpr uint8_t TRIE_EXTERNAL_KEYS = 1;
constexpr uint8_t TRIE_BACKTRACE = 2;
constexpr int TRIE_ID_CODING_SHIFT = 2;
constexpr uint8_t TRIE_ID_CODING_MASK = 7;
constexpr uint32_t TRIE_ID_CODING_UNUSED = 7;
// A node id's low bits number a value within its transition: 8 bits, plus 2 per id-coding step.
constexpr uint32_t TRIE_ID_SHIFT_BASE = 8;
constexpr uint32_t TRIE_ID_SHIFT_STEP = 2;
// The trie addresses node data in 4-byte units. Offsets, character-map indexes and transition ids
// each take the narrowest of 1, 2 or 4 bytes that holds their range.
constexpr size_t TRIE_UNIT = 4;
constexpr uint32_t U8_RANGE = uint32_t(UINT8_MAX) + 1;
constexpr uint32_t U16_RANGE = uint32_t(UINT16_MAX) + 1;
// A dense projection, one node per row, is stored when the rows are fewer than twice the
// ordinals; a sparse one adds each ordinal's row span.
constexpr int64_t DENSE_PROJECTION_FACTOR = 2;
constexpr uint32_t MAX_SPAN_BITS = 64;

// The narrowest width, in bytes, whose values cover `range`.
uint32_t WidthFor(uint64_t range) {
	return range <= U8_RANGE ? sizeof(uint8_t) : range <= U16_RANGE ? sizeof(uint16_t) : sizeof(uint32_t);
}

// The trie's transition-node header, in index units (ptrie: RECORD_HDR_IDX and _DEPTH set).
enum TransHeader : uint32_t {
	HDR_IDX = 0,
	HDR_DEPTH = 1,
	HDR_TRANS_CNT = 2,
	HDR_IDX_CNT = 3,
	HDR_SEQ_LEN = 4,
	HDR_CHILD_CNT = 5,
	FIRST_IDX = 6
};

// UTF-8 (RFC 3629): a lead byte, then continuation bytes of 6 bits each under the 10 tag.
constexpr uint8_t UTF8_CONT = 0x80;
constexpr uint8_t UTF8_CONT_MASK = 0x3f;
constexpr uint8_t UTF8_CONT_MAX = 0xbf;
constexpr int UTF8_CONT_BITS = 6;
constexpr uint8_t UTF8_LEAD_2 = 0xc0, UTF8_LEAD_3 = 0xe0, UTF8_LEAD_4 = 0xf0;
// The first code point that needs 2, 3 and 4 bytes.
constexpr uint32_t UTF8_2_BYTE_MIN = 0x80, UTF8_3_BYTE_MIN = 0x800, UTF8_4_BYTE_MIN = 0x10000;

void AppendUtf8(std::string &out, uint32_t c) {
	auto cont = [&](int shift) {
		out += char(UTF8_CONT | ((c >> shift) & UTF8_CONT_MASK));
	};
	if (c < UTF8_2_BYTE_MIN) {
		out += char(c);
	} else if (c < UTF8_3_BYTE_MIN) {
		out += char(UTF8_LEAD_2 | (c >> UTF8_CONT_BITS));
		cont(0);
	} else if (c < UTF8_4_BYTE_MIN) {
		out += char(UTF8_LEAD_3 | (c >> (2 * UTF8_CONT_BITS)));
		cont(UTF8_CONT_BITS);
		cont(0);
	} else {
		out += char(UTF8_LEAD_4 | (c >> (3 * UTF8_CONT_BITS)));
		cont(2 * UTF8_CONT_BITS);
		cont(UTF8_CONT_BITS);
		cont(0);
	}
}

[[noreturn]] void Damaged(const std::string &what) {
	throw SraFormatError("the skey index is damaged: " + what);
}

} // namespace

// The trie header follows ncbi-vdb's PTrieMakeInt (libs/klib/ptrie.c), and the projection its
// KRPTrieIndexInit_v3_v4 and KPTrieIndexInitID2Ord (libs/kdb/rtrieidx-v2.c): the same checks in
// the same order (public domain; THIRD_PARTY_LICENSES.md).
SkeyIndex::SkeyIndex(std::vector<uint8_t> bytes) : bytes_(std::move(bytes)) {
	const auto &b = bytes_;
	if (b.size() < KDB_FILE_HEADER_SIZE) {
		Damaged("shorter than its header");
	}
	CheckByteOrder(LoadU32(b.data()), "the skey index");
	uint32_t version = LoadU32(b.data() + sizeof(uint32_t));
	if (version != SKEY_VERSION) {
		throw SraUnsupportedError("skey index version " + std::to_string(version) +
		                          " is not supported (public runs use version 4)");
	}
	if (b.size() < SKEY_HEADER_SIZE) {
		Damaged("shorter than its header");
	}
	if (LoadU32(b.data() + SKEY_TYPE_AT) != SKEY_TEXT_INDEX) {
		throw SraUnsupportedError("an skey index that is not a text index");
	}
	first_ = int64_t(LoadU64(b.data() + SKEY_FIRST_AT));
	last_ = int64_t(LoadU64(b.data() + SKEY_LAST_AT));
	uint32_t span_bits = LoadU16(b.data() + SKEY_SPAN_BITS_AT);

	// The persisted trie (P_Trie): counts, flags, the character map, node offsets, node data.
	trie_ = SKEY_HEADER_SIZE;
	if (b.size() < trie_ + TRIE_HEADER_SIZE) {
		Damaged("truncated trie header");
	}
	num_trans_ = LoadU32(&b[trie_ + TRIE_NUM_TRANS_AT]);
	uint32_t num_nodes = LoadU32(&b[trie_ + TRIE_NUM_NODES_AT]);
	data_size_ = LoadU32(&b[trie_ + TRIE_DATA_SIZE_AT]) +
	             (uint64_t(b[trie_ + TRIE_DATA_SIZE_HIGH_AT]) << TRIE_DATA_SIZE_HIGH_SHIFT);
	uint8_t keys = b[trie_ + TRIE_KEY_FLAGS_AT];
	width_ = LoadU16(&b[trie_ + TRIE_WIDTH_AT]);
	if (keys & TRIE_EXTERNAL_KEYS) {
		throw SraUnsupportedError("an skey index whose keys are stored outside it");
	}
	backtrace_ = (keys & TRIE_BACKTRACE) != 0;
	if (!backtrace_) {
		throw SraUnsupportedError("an skey index without backtrace links, whose keys cannot be rebuilt");
	}
	uint32_t id_coding = (keys >> TRIE_ID_CODING_SHIFT) & TRIE_ID_CODING_MASK;
	if (id_coding == TRIE_ID_CODING_UNUSED) {
		throw SraUnsupportedError("skey trie node id coding 7, which no public run surveyed used");
	}
	id_shift_ = TRIE_ID_SHIFT_BASE + TRIE_ID_SHIFT_STEP * id_coding;
	// The character map: the code point of each character index.
	rmap_ = trie_ + TRIE_HEADER_SIZE;
	if (b.size() < rmap_ + size_t(width_) * sizeof(uint32_t)) {
		Damaged("truncated character map");
	}
	for (uint32_t i = 0; i < width_; i++) {
		if (LoadU32(&b[rmap_ + sizeof(uint32_t) * i]) == 0) {
			unmapped_idx_ = i + 1;
		}
	}
	// One offset per transition, into node data counted in units, packed into whole units.
	trans_off_width_ = WidthFor((data_size_ + TRIE_UNIT - 1) / TRIE_UNIT);
	size_t trans_off_words = (size_t(num_trans_) * trans_off_width_ + TRIE_UNIT - 1) / TRIE_UNIT;
	trans_off_ = rmap_ + size_t(width_) * sizeof(uint32_t);
	data_ = trans_off_ + trans_off_words * TRIE_UNIT;
	size_t trie_end = data_ + size_t(data_size_);
	if (trie_end > b.size()) {
		Damaged("the trie runs past the file");
	}

	// The projection: a count, then the trie node of each ordinal, then (for a sparse index) the
	// bit-packed number of rows each ordinal covers.
	if (trie_end == b.size()) {
		return; // no projection: no row has a key
	}
	// Node ids, the count and the spans' bits are all 4-byte words or packed into them.
	constexpr size_t NODE_ID_SIZE = sizeof(uint32_t);
	size_t rest = b.size() - trie_end;
	if (rest < size_t(num_nodes) * NODE_ID_SIZE || rest < sizeof(uint32_t)) {
		Damaged("the projection is truncated");
	}
	count_ = LoadU32(&b[trie_end]);
	ord2node_ = trie_end + sizeof(uint32_t);
	rest -= sizeof(uint32_t);
	if (last_ - first_ < int64_t(count_) * DENSE_PROJECTION_FACTOR) {
		count_ = uint32_t(last_ - first_ + 1);
		if (rest != size_t(count_) * NODE_ID_SIZE) {
			Damaged("the projection does not have one node per row");
		}
		return;
	}
	if (!((rest == NODE_ID_SIZE && count_ == 1) || rest > size_t(count_) * NODE_ID_SIZE) || count_ == 0) {
		Damaged("the projection is malformed");
	}
	rest -= size_t(count_) * NODE_ID_SIZE;
	uint64_t n = count_ - 1;
	if (span_bits == 0 || span_bits > MAX_SPAN_BITS || span_bits * n > rest * BITS_PER_BYTE) {
		Damaged("the row spans are truncated");
	}
	// Spans are packed most significant bit first; the first ordinal starts at row 0.
	const uint8_t *packed = &b[ord2node_ + size_t(count_) * NODE_ID_SIZE];
	id2ord_.reserve(count_);
	id2ord_.assign(1, 0);
	uint64_t bit = 0;
	for (uint64_t i = 0; i < n; i++) {
		uint64_t v = 0;
		for (uint32_t k = 0; k < span_bits; k++, bit++) {
			const int from_top = BITS_PER_BYTE - 1 - int(bit % BITS_PER_BYTE);
			v = (v << 1) | ((packed[bit / BITS_PER_BYTE] >> from_top) & 1);
		}
		id2ord_.push_back(id2ord_.back() + v);
	}
}

uint32_t SkeyIndex::Ordinal(int64_t row_id) const {
	if (count_ == 0 || row_id < first_ || row_id > last_) {
		return 0;
	}
	uint64_t id = uint64_t(row_id - first_);
	if (id2ord_.empty()) {
		return uint32_t(id + 1);
	}
	// The ordinal whose first row is the last at or before id.
	return uint32_t(std::upper_bound(id2ord_.begin(), id2ord_.end(), id) - id2ord_.begin());
}

std::string SkeyIndex::Project(int64_t row_id) const {
	return KeyOfOrdinal(Ordinal(row_id));
}

std::string SkeyIndex::KeyOfOrdinal(uint32_t ord) const {
	if (ord == 0 || ord > count_) {
		return std::string();
	}
	uint32_t node = LoadU32(&bytes_[ord2node_ + sizeof(uint32_t) * size_t(ord - 1)]);
	std::string key = node == 0 ? std::string() : KeyOfNode(node);
	while (!key.empty() && key.back() == '\0') {
		key.pop_back();
	}
	return key;
}

// Follows ncbi-vdb's PTrieInitNode and PTNodeMakeKey (libs/klib/ptrie.c; public domain,
// THIRD_PARTY_LICENSES.md): the node layout is read, and the key rebuilt, step for step as there.
std::string SkeyIndex::KeyOfNode(uint32_t node_id) const {
	const auto &b = bytes_;
	struct Trans {
		uint32_t code = 0;
		uint32_t depth = 0;
		uint32_t dad = 0;
		size_t vals = 0; // 0: no values
	};
	const uint32_t idx_width = WidthFor(width_);
	const uint32_t tid_width = WidthFor(num_trans_);
	auto read = [&](size_t at, uint32_t w) -> uint32_t {
		if (at + w > b.size()) {
			Damaged("a trie node runs past the file");
		}
		return w == sizeof(uint8_t) ? b[at] : w == sizeof(uint16_t) ? LoadU16(&b[at]) : LoadU32(&b[at]);
	};
	auto trans = [&](uint32_t tid) {
		if (tid == 0 || tid > num_trans_) {
			Damaged("trie node " + std::to_string(tid) + " out of range");
		}
		size_t base =
		    data_ + TRIE_UNIT * size_t(read(trans_off_ + size_t(tid - 1) * trans_off_width_, trans_off_width_));
		auto idx = [&](uint32_t k) {
			return read(base + size_t(k) * idx_width, idx_width);
		};
		Trans t;
		t.code = idx(HDR_IDX);
		t.depth = idx(HDR_DEPTH);
		uint32_t icnt = idx(HDR_IDX_CNT);
		size_t dad;
		uint32_t child_count = 0;
		if (icnt == 0) {
			dad = base + size_t(HDR_SEQ_LEN) * idx_width;
		} else {
			uint32_t slen = idx(HDR_SEQ_LEN);
			child_count = idx(HDR_CHILD_CNT);
			// The sequence's bit map, rounded up to whole bytes.
			dad = base + size_t(FIRST_IDX + icnt) * idx_width + (slen + BITS_PER_BYTE - 1) / BITS_PER_BYTE;
		}
		dad = (dad + tid_width - 1) & ~size_t(tid_width - 1);
		t.dad = read(dad, tid_width);
		// The value tree starts at the next unit boundary; a nonzero byte in the padding says it is
		// there, and on a boundary it always is.
		size_t align2 = dad + tid_width + size_t(tid_width) * child_count;
		if (align2 % TRIE_UNIT != 0) {
			if (read(align2, sizeof(uint8_t)) != 0) {
				t.vals = (align2 | (TRIE_UNIT - 1)) + 1;
			}
		} else {
			t.vals = align2;
		}
		return t;
	};

	uint32_t tid = ((node_id - 1) >> id_shift_) + 1;
	uint32_t btid = ((node_id - 1) & ((1u << id_shift_) - 1)) + 1;
	Trans t = trans(tid);
	if (t.vals == 0) {
		Damaged("trie node " + std::to_string(tid) + " holds no values");
	}
	// The value node: a PBSTree entry holding the key's remaining text, NUL-terminated.
	uint32_t n = read(t.vals, sizeof(uint32_t));
	if (btid > n) {
		Damaged("value " + std::to_string(btid) + " out of range");
	}
	uint32_t data_size = read(t.vals + sizeof(uint32_t), sizeof(uint32_t));
	uint32_t w = WidthFor(data_size);
	const size_t offsets = t.vals + PBST_HEADER_SIZE;
	size_t payload = offsets + size_t(w) * n;
	uint32_t begin = read(offsets + size_t(w) * (btid - 1), w);
	uint32_t end = btid < n ? read(offsets + size_t(w) * btid, w) : data_size;
	if (begin > end || payload + end > b.size()) {
		Damaged("a value runs past the file");
	}
	std::string suffix;
	for (size_t i = payload + begin; i < payload + end && b[i] != 0; i++) {
		suffix += char(b[i]);
	}
	// The key's leading characters are the transitions from the root: walk back up.
	std::vector<uint32_t> chars(t.depth);
	for (uint32_t i = t.depth; i > 0; i--) {
		uint32_t z = t.code;
		chars[i - 1] =
		    z + 1 == unmapped_idx_ ? uint32_t('?') : read(rmap_ + sizeof(uint32_t) * size_t(z), sizeof(uint32_t));
		if (i > 1) {
			if (t.dad == 0) {
				Damaged("a trie path ends before the root");
			}
			t = trans(t.dad);
		}
	}
	std::string key;
	for (uint32_t c : chars) {
		AppendUtf8(key, c);
	}
	return key + suffix;
}

// ---------------------------------------------------------------------------------------------

namespace {
constexpr size_t NAME_BUFFER = 1024; // format-spot-name.c's sname[1024]
// The most characters %d prints for a 32-bit value: ten digits and a sign.
constexpr size_t MAX_INT32_DIGITS = 11;
// $X%n and $Y%n add n tiles of this many units to the coordinate.
constexpr uint32_t TILE_OFFSET = 24 * 1024;
// $Q packs X (20 bits) above Y (12 bits) and prints the result in base 36 as five characters.
constexpr char Q_ALPHABET[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
constexpr uint32_t Q_BASE = sizeof(Q_ALPHABET) - 1;
constexpr int Q_LENGTH = 5;
constexpr int Q_Y_BITS = 12;
constexpr uint32_t Q_MAX_Y = (uint32_t(1) << Q_Y_BITS) - 1;
constexpr uint32_t Q_MAX_X = (uint32_t(1) << 20) - 1;
constexpr uint32_t Q_LIMIT = Q_BASE * Q_BASE * Q_BASE * Q_BASE * Q_BASE; // Q_BASE^Q_LENGTH
// "$X" and the like: the '$' and the letter.
constexpr size_t TOKEN_SIZE = 2;

[[noreturn]] void TooLong(int64_t row_id) {
	throw SraFormatError("spot " + std::to_string(row_id) + ": its name would exceed " +
	                     std::to_string(NAME_BUFFER - 1) + " bytes");
}
} // namespace

std::string FormatSpotName(const std::string &fmt, std::optional<int32_t> x_cell, std::optional<int32_t> y_cell,
                           int64_t row_id, const std::string &spot_name) {
	if (fmt.empty()) {
		if (!spot_name.empty()) {
			return spot_name;
		}
		throw SraFormatError("spot " + std::to_string(row_id) + " has no name");
	}
	// A port of ncbi-vdb's libs/sraxf/format-spot-name.c (public domain; THIRD_PARTY_LICENSES.md),
	// including its bounds: the coordinates are uint32_t printed with %d, and a '$' before anything
	// but X, Y, Q or R is copied.
	uint32_t x = x_cell ? uint32_t(*x_cell) : 0;
	uint32_t y = y_cell ? uint32_t(*y_cell) : 0;
	const size_t size = fmt.size();
	std::string out;
	size_t i = 0;
	auto coordinate = [&](bool have, uint32_t value, const char *axis) {
		if (out.size() > NAME_BUFFER - MAX_INT32_DIGITS) {
			TooLong(row_id);
		}
		if (!have) {
			throw SraFormatError("spot " + std::to_string(row_id) + ": its name format needs " + axis +
			                     ", which the spot does not have");
		}
		i += TOKEN_SIZE;
		if (i < size - 1 && fmt[i] == '%' && std::isdigit(static_cast<unsigned char>(fmt[i + 1]))) {
			value += TILE_OFFSET * uint32_t(fmt[i + 1] - '0');
			i += TOKEN_SIZE;
		}
		out += std::to_string(int32_t(value));
	};
	while (i < size - 1) {
		if (fmt[i] == '$') {
			switch (fmt[i + 1]) {
			case 'x':
			case 'X':
				coordinate(x_cell.has_value(), x, "X");
				break;
			case 'y':
			case 'Y':
				coordinate(y_cell.has_value(), y, "Y");
				break;
			case 'q':
			case 'Q': {
				if (out.size() > NAME_BUFFER - Q_LENGTH) {
					TooLong(row_id);
				}
				if (y > Q_MAX_Y) {
					throw SraFormatError("spot " + std::to_string(row_id) + ": Y is too large for a $Q name");
				}
				uint32_t q = x << Q_Y_BITS | y;
				if (x > Q_MAX_X || q >= Q_LIMIT) {
					throw SraFormatError("spot " + std::to_string(row_id) + ": X is too large for a $Q name");
				}
				char code[Q_LENGTH];
				for (int k = Q_LENGTH - 1; k > 0; k--) {
					code[k] = Q_ALPHABET[q % Q_BASE];
					q /= Q_BASE;
				}
				code[0] = Q_ALPHABET[q];
				out.append(code, Q_LENGTH);
				i += TOKEN_SIZE;
				break;
			}
			case 'r':
			case 'R':
				out += std::to_string(row_id);
				i += TOKEN_SIZE;
				break;
			default:
				out += fmt[i++];
				break;
			}
		} else {
			if (out.size() >= NAME_BUFFER) {
				TooLong(row_id);
			}
			out += fmt[i++];
		}
		if (out.size() > NAME_BUFFER - 1) {
			TooLong(row_id);
		}
	}
	if (i == size - 1) {
		if (out.size() >= NAME_BUFFER) {
			TooLong(row_id);
		}
		out += fmt[i];
	}
	return out;
}

std::string FormatSpotNameNoCoord(const std::string &fmt, int64_t row_id, const std::string &spot_name) {
	// format_spot_name_no_coord from the same file: only $R is substituted.
	if (fmt.empty()) {
		if (!spot_name.empty()) {
			return spot_name;
		}
		throw SraFormatError("spot " + std::to_string(row_id) + " has no name");
	}
	const size_t size = fmt.size();
	std::string out;
	size_t i = 0;
	while (i < size - 1) {
		if (fmt[i] == '$' && (fmt[i + 1] == 'r' || fmt[i + 1] == 'R')) {
			out += std::to_string(row_id);
			i += TOKEN_SIZE;
		} else {
			if (out.size() >= NAME_BUFFER) {
				TooLong(row_id);
			}
			out += fmt[i++];
		}
		if (out.size() > NAME_BUFFER - 1) {
			TooLong(row_id);
		}
	}
	if (i == size - 1) {
		if (out.size() >= NAME_BUFFER) {
			TooLong(row_id);
		}
		out += fmt[i];
	}
	return out;
}

void CheckSpotName(const std::string &name, int64_t row_id) {
	// RFC 3629: no overlong forms, no surrogates, nothing past U+10FFFF. Written here rather than
	// with DuckDB's Utf8Proc because src/sra_*.cpp stays free of DuckDB (docs/internals/sra-reader.md).
	const auto *p = reinterpret_cast<const unsigned char *>(name.data());
	const size_t n = name.size();
	for (size_t i = 0; i < n;) {
		unsigned char c = p[i];
		size_t len;
		// RFC 3629's table of well-formed sequences: the lead bytes of each length, and the narrower
		// range a few leads allow for the second byte (no overlong forms, surrogates, or code points
		// past U+10FFFF).
		constexpr unsigned char LEAD_2_MIN = 0xc2, LEAD_2_MAX = 0xdf, LEAD_3_MAX = 0xef, LEAD_4_MAX = 0xf4;
		constexpr unsigned char LEAD_SURROGATES = 0xed; // its second byte must stay below the surrogates
		constexpr unsigned char E0_SECOND_MIN = 0xa0, ED_SECOND_MAX = 0x9f, F0_SECOND_MIN = 0x90, F4_SECOND_MAX = 0x8f;
		unsigned char lo = UTF8_CONT, hi = UTF8_CONT_MAX; // bounds of the second byte
		if (c < UTF8_2_BYTE_MIN) {
			i++;
			continue;
		} else if (c >= LEAD_2_MIN && c <= LEAD_2_MAX) {
			len = 2;
		} else if (c >= UTF8_LEAD_3 && c <= LEAD_3_MAX) {
			len = 3;
			lo = c == UTF8_LEAD_3 ? E0_SECOND_MIN : UTF8_CONT;
			hi = c == LEAD_SURROGATES ? ED_SECOND_MAX : UTF8_CONT_MAX;
		} else if (c >= UTF8_LEAD_4 && c <= LEAD_4_MAX) {
			len = 4;
			lo = c == UTF8_LEAD_4 ? F0_SECOND_MIN : UTF8_CONT;
			hi = c == LEAD_4_MAX ? F4_SECOND_MAX : UTF8_CONT_MAX;
		} else {
			len = 0;
		}
		bool valid = len != 0 && n - i >= len && p[i + 1] >= lo && p[i + 1] <= hi;
		for (size_t k = 2; valid && k < len; k++) {
			valid = p[i + k] >= UTF8_CONT && p[i + k] <= UTF8_CONT_MAX;
		}
		if (!valid) {
			// A run that stores Latin-1 names is well formed, but its names are not text DuckDB holds.
			throw SraUnsupportedError("spot " + std::to_string(row_id) + ": its name is not valid UTF-8 (byte " +
			                          std::to_string(i + 1) + ")");
		}
		i += len;
	}
}

} // namespace sra
} // namespace miint
