#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

// One per-base row emitted by PileupWalker::Walk. `ref_pos` is 1-based
// (matches SAM POS). `query_is_null` flags D/N positions where the read
// has no base; `qual_is_null` is true on D/N positions OR when the input
// qual list for the alignment is NULL.
struct PileupRow {
	std::string ref_id;
	std::int64_t ref_pos;
	std::string read_id;
	char ref_base;
	char query_base;
	std::uint8_t query_qual;
	bool query_is_null;
	bool qual_is_null;
};

// Walks a single alignment's CIGAR against a reference sequence and emits
// per-base PileupRow entries by appending to `out`.
//
// V1 op handling (Karst-protocol UMI use case — SNP positions only):
//   M/=/X : emit one row per ref pos; advance both ref and query cursors
//   D/N   : emit one row per ref pos with query_base/query_qual NULL;
//           advance ref only
//   I     : advance query only — insertion rows NOT emitted (deferred to v2;
//           insertions matter for indel-aware callers but not for the
//           Karst variant pipeline that only needs SNP positions)
//   S     : advance query only (clip drops from output)
//   H/P   : no-op
//
// Note on N: the SAM spec defines N as a reference skip (intron). V1 treats
// N identically to D for the UMI use case, which produces NULL query_base
// rows across intron spans. Downstream RNA-seq callers may want to filter
// these out or request a future N-aware mode.
class PileupWalker {
public:
	static void Walk(const std::string &cigar, const std::string &ref_id, const std::string &ref_seq,
	                 std::int64_t align_start_1based, const std::string &read_id, const std::string &seq,
	                 const std::uint8_t *qual_data, std::size_t qual_len, bool qual_is_null,
	                 std::vector<PileupRow> &out);
};

} // namespace miint
