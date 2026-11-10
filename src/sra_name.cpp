#include "sra_name.hpp"

#include "sra_source.hpp"

#include <algorithm>
#include <cctype>

namespace miint {
namespace sra {

namespace {

constexpr uint32_t ORDER_NATIVE = 0x05031988;
constexpr uint32_t ORDER_SWAPPED = 0x88190305;

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

void AppendUtf8(std::string &out, uint32_t c) {
	if (c < 0x80) {
		out += char(c);
	} else if (c < 0x800) {
		out += char(0xc0 | (c >> 6));
		out += char(0x80 | (c & 0x3f));
	} else if (c < 0x10000) {
		out += char(0xe0 | (c >> 12));
		out += char(0x80 | ((c >> 6) & 0x3f));
		out += char(0x80 | (c & 0x3f));
	} else {
		out += char(0xf0 | (c >> 18));
		out += char(0x80 | ((c >> 12) & 0x3f));
		out += char(0x80 | ((c >> 6) & 0x3f));
		out += char(0x80 | (c & 0x3f));
	}
}

[[noreturn]] void Damaged(const std::string &what) {
	throw SraFormatError("the skey index is damaged: " + what);
}

} // namespace

SkeyIndex::SkeyIndex(std::vector<uint8_t> bytes) : bytes_(std::move(bytes)) {
	const auto &b = bytes_;
	if (b.size() < 8) {
		Damaged("shorter than its header");
	}
	uint32_t order = LoadU32(b.data());
	if (order == ORDER_SWAPPED) {
		throw SraUnsupportedError("a big-endian skey index, which no public run surveyed had");
	}
	if (order != ORDER_NATIVE) {
		Damaged("invalid byte-order tag");
	}
	uint32_t version = LoadU32(b.data() + 4);
	if (version != 4) {
		throw SraUnsupportedError("skey index version " + std::to_string(version) +
		                          " is not supported (public runs use version 4)");
	}
	if (b.size() < 40) {
		Damaged("shorter than its header");
	}
	if (LoadU32(b.data() + 8) != 0) {
		throw SraUnsupportedError("an skey index that is not a text index");
	}
	first_ = int64_t(LoadU64(b.data() + 16));
	last_ = int64_t(LoadU64(b.data() + 24));
	uint32_t span_bits = LoadU16(b.data() + 34);

	// The persisted trie (P_Trie): counts, flags, the character map, node offsets, node data.
	trie_ = 40;
	if (b.size() < trie_ + 16) {
		Damaged("truncated trie header");
	}
	num_trans_ = LoadU32(&b[trie_]);
	uint32_t num_nodes = LoadU32(&b[trie_ + 4]);
	data_size_ = LoadU32(&b[trie_ + 8]) + (uint64_t(b[trie_ + 13]) << 32);
	uint8_t keys = b[trie_ + 12];
	width_ = LoadU16(&b[trie_ + 14]);
	if (keys & 1) {
		throw SraUnsupportedError("an skey index whose keys are stored outside it");
	}
	backtrace_ = (keys & 2) != 0;
	if (!backtrace_) {
		throw SraUnsupportedError("an skey index without backtrace links, whose keys cannot be rebuilt");
	}
	uint32_t id_coding = (keys >> 2) & 7;
	if (id_coding == 7) {
		throw SraUnsupportedError("skey trie node id coding 7, which no public run surveyed used");
	}
	id_shift_ = 8 + 2 * id_coding;
	rmap_ = trie_ + 16;
	if (b.size() < rmap_ + size_t(width_) * 4) {
		Damaged("truncated character map");
	}
	for (uint32_t i = 0; i < width_; i++) {
		if (LoadU32(&b[rmap_ + 4 * i]) == 0) {
			unmapped_idx_ = i + 1;
		}
	}
	size_t trans_off_words;
	if (data_size_ <= 256 * 4) {
		trans_off_width_ = 1;
		trans_off_words = (num_trans_ + 3) >> 2;
	} else if (data_size_ <= 65536 * 4) {
		trans_off_width_ = 2;
		trans_off_words = (num_trans_ + 1) >> 1;
	} else {
		trans_off_width_ = 4;
		trans_off_words = num_trans_;
	}
	trans_off_ = rmap_ + size_t(width_) * 4;
	data_ = trans_off_ + trans_off_words * 4;
	size_t trie_end = data_ + size_t(data_size_);
	if (trie_end > b.size()) {
		Damaged("the trie runs past the file");
	}

	// The projection: a count, then the trie node of each ordinal, then (for a sparse index) the
	// bit-packed number of rows each ordinal covers.
	if (trie_end == b.size()) {
		return; // no projection: no row has a key
	}
	size_t rest = b.size() - trie_end;
	if (rest < size_t(num_nodes) * 4 || rest < 4) {
		Damaged("the projection is truncated");
	}
	count_ = LoadU32(&b[trie_end]);
	ord2node_ = trie_end + 4;
	rest -= 4;
	if (last_ - first_ < int64_t(count_) * 2) {
		count_ = uint32_t(last_ - first_ + 1);
		if (rest != size_t(count_) * 4) {
			Damaged("the projection does not have one node per row");
		}
		return;
	}
	if (!((rest == 4 && count_ == 1) || rest > size_t(count_) * 4) || count_ == 0) {
		Damaged("the projection is malformed");
	}
	rest -= size_t(count_) * 4;
	uint64_t n = count_ - 1;
	if (span_bits == 0 || span_bits > 64 || span_bits * n > rest * 8) {
		Damaged("the row spans are truncated");
	}
	// Spans are packed most significant bit first; the first ordinal starts at row 0.
	const uint8_t *packed = &b[ord2node_ + size_t(count_) * 4];
	id2ord_.assign(1, 0);
	uint64_t bit = 0;
	for (uint64_t i = 0; i < n; i++) {
		uint64_t v = 0;
		for (uint32_t k = 0; k < span_bits; k++, bit++) {
			v = (v << 1) | ((packed[bit >> 3] >> (7 - (bit & 7))) & 1);
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
	uint32_t node = LoadU32(&bytes_[ord2node_ + 4 * size_t(ord - 1)]);
	std::string key = node == 0 ? std::string() : KeyOfNode(node);
	while (!key.empty() && key.back() == '\0') {
		key.pop_back();
	}
	return key;
}

std::string SkeyIndex::KeyOfNode(uint32_t node_id) const {
	const auto &b = bytes_;
	struct Trans {
		uint32_t code = 0;
		uint32_t depth = 0;
		uint32_t dad = 0;
		size_t vals = 0; // 0: no values
	};
	const uint32_t idx_width = width_ <= 256 ? 1 : 2;
	const uint32_t tid_width = num_trans_ <= 256 ? 1 : num_trans_ <= 65536 ? 2 : 4;
	auto read = [&](size_t at, uint32_t w) -> uint32_t {
		if (at + w > b.size()) {
			Damaged("a trie node runs past the file");
		}
		return w == 1 ? b[at] : w == 2 ? LoadU16(&b[at]) : LoadU32(&b[at]);
	};
	auto trans = [&](uint32_t tid) {
		if (tid == 0 || tid > num_trans_) {
			Damaged("trie node " + std::to_string(tid) + " out of range");
		}
		size_t base = data_ + 4 * size_t(read(trans_off_ + size_t(tid - 1) * trans_off_width_, trans_off_width_));
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
			dad = base + size_t(FIRST_IDX + icnt) * idx_width + ((slen + 7) >> 3);
		}
		dad = (dad + tid_width - 1) & ~size_t(tid_width - 1);
		t.dad = read(dad, tid_width);
		size_t align2 = dad + tid_width + size_t(tid_width) * child_count;
		if (align2 & 3) {
			if (read(align2, 1) != 0) {
				t.vals = (align2 | 3) + 1;
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
	uint32_t n = read(t.vals, 4);
	if (btid > n) {
		Damaged("value " + std::to_string(btid) + " out of range");
	}
	uint32_t data_size = read(t.vals + 4, 4);
	uint32_t w = data_size <= 256 ? 1 : data_size <= 65536 ? 2 : 4;
	size_t payload = t.vals + 8 + size_t(w) * n;
	uint32_t begin = read(t.vals + 8 + size_t(w) * (btid - 1), w);
	uint32_t end = btid < n ? read(t.vals + 8 + size_t(w) * btid, w) : data_size;
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
		chars[i - 1] = z + 1 == unmapped_idx_ ? uint32_t('?') : read(rmap_ + 4 * size_t(z), 4);
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
const char *Q_ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

[[noreturn]] void TooLong(int64_t row_id) {
	throw SraFormatError("spot " + std::to_string(row_id) + ": its name would exceed 1023 bytes");
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
		if (out.size() > NAME_BUFFER - 11) {
			TooLong(row_id);
		}
		if (!have) {
			throw SraFormatError("spot " + std::to_string(row_id) + ": its name format needs " + axis +
			                     ", which the spot does not have");
		}
		i += 2;
		if (i < size - 1 && fmt[i] == '%' && std::isdigit(static_cast<unsigned char>(fmt[i + 1]))) {
			value += uint32_t(24 * 1024 * (fmt[i + 1] - '0'));
			i += 2;
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
				if (out.size() > NAME_BUFFER - 5) {
					TooLong(row_id);
				}
				if (y > 0xFFF) {
					throw SraFormatError("spot " + std::to_string(row_id) + ": Y is too large for a $Q name");
				}
				uint32_t q = x << 12 | y;
				if (x > 0xFFFFF || q >= 36u * 36 * 36 * 36 * 36) {
					throw SraFormatError("spot " + std::to_string(row_id) + ": X is too large for a $Q name");
				}
				char code[5];
				for (int k = 4; k > 0; k--) {
					code[k] = Q_ALPHABET[q % 36];
					q /= 36;
				}
				code[0] = Q_ALPHABET[q];
				out.append(code, 5);
				i += 2;
				break;
			}
			case 'r':
			case 'R':
				out += std::to_string(row_id);
				i += 2;
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
			i += 2;
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
		unsigned char lo = 0x80, hi = 0xbf; // bounds of the second byte
		if (c < 0x80) {
			i++;
			continue;
		} else if (c >= 0xc2 && c <= 0xdf) {
			len = 2;
		} else if (c >= 0xe0 && c <= 0xef) {
			len = 3;
			lo = c == 0xe0 ? 0xa0 : 0x80;
			hi = c == 0xed ? 0x9f : 0xbf;
		} else if (c >= 0xf0 && c <= 0xf4) {
			len = 4;
			lo = c == 0xf0 ? 0x90 : 0x80;
			hi = c == 0xf4 ? 0x8f : 0xbf;
		} else {
			len = 0;
		}
		bool valid = len != 0 && n - i >= len && p[i + 1] >= lo && p[i + 1] <= hi;
		for (size_t k = 2; valid && k < len; k++) {
			valid = p[i + k] >= 0x80 && p[i + k] <= 0xbf;
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
