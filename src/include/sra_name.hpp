#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace miint {
namespace sra {

// Spot names, as VDB reconstructs them for fasterq-dump's $sn.
//
// Most loaders do not store a name per spot. They store a format per run of similar names in the
// table's `idx/skey` index (e.g. "M00763:72:000000000-KK57D:1:1101:$X:$Y" for one tile) and the
// coordinates per spot, and the name is the format with the coordinates substituted.

// The `skey` index: a persisted trie of name formats with a projection from row id to the format
// covering it (KDB text index version 4, as every surveyed public run has). Other versions, and the
// id coding no surveyed index used, raise SraUnsupportedError.
class SkeyIndex {
public:
	explicit SkeyIndex(std::vector<uint8_t> bytes);

	// The key (name format) of the row, without trailing NULs; empty if the index has none. Safe to
	// call from several threads.
	std::string Project(int64_t row_id) const;
	// Project in two steps, for a caller that keeps the last key: consecutive rows usually share
	// one (a tile), so Ordinal changes far less often than the row. Ordinal 0 has no key.
	uint32_t Ordinal(int64_t row_id) const;
	std::string KeyOfOrdinal(uint32_t ordinal) const;

private:
	std::vector<uint8_t> bytes_;
	int64_t first_ = 0;
	int64_t last_ = 0;
	// The trie.
	size_t trie_ = 0;
	uint32_t num_trans_ = 0;
	uint16_t width_ = 0;
	size_t rmap_ = 0;
	size_t trans_off_ = 0;
	uint32_t trans_off_width_ = 0;
	size_t data_ = 0;
	uint64_t data_size_ = 0;
	bool backtrace_ = false;
	uint32_t id_shift_ = 0;
	uint32_t unmapped_idx_ = 0;
	// The projection: ordinal (one-based) -> trie node id, and for a sparse index the first row
	// (relative to first_) of each ordinal.
	size_t ord2node_ = 0;
	uint32_t count_ = 0;
	std::vector<uint64_t> id2ord_;

	std::string KeyOfNode(uint32_t node_id) const;
};

// INSDC:SRA:format_spot_name: the format with $X / $Y (optionally followed by %d, adding
// 24576 * d), $Q (454 base-36 of X and Y) and $R (the row id) substituted. An empty format yields
// spot_name verbatim. Throws SraFormatError where VDB fails the row: a coordinate the format needs
// is missing, the spot has neither format nor name, or the name would exceed 1023 bytes.
std::string FormatSpotName(const std::string &fmt, std::optional<int32_t> x, std::optional<int32_t> y, int64_t row_id,
                           const std::string &spot_name);

// Throws SraUnsupportedError, naming the spot, unless the name is valid UTF-8: names become VARCHAR
// values, which must be.
void CheckSpotName(const std::string &name, int64_t row_id);

// INSDC:SRA:format_spot_name_no_coord: only $R is substituted.
std::string FormatSpotNameNoCoord(const std::string &fmt, int64_t row_id, const std::string &spot_name);

} // namespace sra
} // namespace miint
