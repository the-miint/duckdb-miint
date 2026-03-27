#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace miint {

// 8-mer inverted index for fast candidate parent search in UCHIME chimera detection.
//
// DNA sequences are indexed by their unique 8-mers (presence/absence per sequence,
// not multiplicity). Queries are partitioned into 4 equal chunks; for each chunk,
// the top 4 references by shared unique 8-mer count are retained. After deduplication
// across chunks, the top 16 candidates are returned.
//
// Sequence IDs are assigned internally by add_sequence() in sequential order
// starting from 0. The returned ID is the index into the caller's sequence vector.
//
// The flat array of 65536 entries (4^8 possible 8-mers) uses ~1.5MB of vector
// metadata. This is allocated once per index lifetime. For lookup performance during
// find_candidates(), a flat array beats an unordered_map because k-mer codes map
// directly to array indices with zero hashing overhead.
//
// Thread safety: add_sequence() is NOT thread-safe. find_candidates() is
// safe for concurrent calls after indexing is complete (read-only on the index;
// each call uses thread_local scratch buffers for hit counting).
class KmerIndex {
public:
	static constexpr int K = 8;
	static constexpr int TOP_PER_CHUNK = 4;
	static constexpr int NUM_CHUNKS = 4;
	static constexpr int MAX_CANDIDATES = TOP_PER_CHUNK * NUM_CHUNKS; // 16 after dedup

	KmerIndex() = default;

	// Add a reference sequence to the index. Returns the assigned seq_id
	// (sequential starting from 0). 8-mers containing ambiguous bases (non-ACGT)
	// are skipped. Posting lists store each seq_id at most once per k-mer
	// (presence/absence model, matching UCHIME's candidate search semantics).
	uint32_t add_sequence(const std::string &sequence);

	// Find candidate parent sequences for a query.
	// Partitions query into NUM_CHUNKS equal chunks, finds top TOP_PER_CHUNK
	// references per chunk by shared unique 8-mer count, deduplicates across chunks,
	// returns at most MAX_CANDIDATES seq_ids sorted by selected hit count (descending).
	//
	// "Selected hit count" means: for each chunk, only the top-4 candidates'
	// hit counts are accumulated. A sequence that is rank 5+ in a given chunk
	// does not get that chunk's contribution. This matches UCHIME's per-chunk
	// top-N selection strategy.
	[[nodiscard]] std::vector<uint32_t> find_candidates(const std::string &query) const;

	// Partition a query of given length into NUM_CHUNKS approximately equal chunks.
	// Returns vector of (start_position, length) pairs.
	[[nodiscard]] static std::vector<std::pair<size_t, size_t>> partition_query(size_t query_len);

	// Number of sequences indexed.
	[[nodiscard]] uint32_t size() const {
		return seq_count_;
	}

private:
	// Flat array for all possible 8-mers (4^8 = 65536).
	// Each entry is a sorted list of unique seq_ids containing that 8-mer.
	static constexpr size_t NUM_KMERS = 65536; // 4^8
	std::array<std::vector<uint32_t>, NUM_KMERS> index_;
	uint32_t seq_count_ = 0;

	// Encode an 8-mer starting at position pos in seq.
	// Returns false if any base is ambiguous (non-ACGT).
	// On success, sets out to the 16-bit encoded k-mer.
	// Caller must ensure pos + K <= seq.size().
	static bool encode_kmer(const std::string &seq, size_t pos, uint16_t &out);
};

} // namespace miint
