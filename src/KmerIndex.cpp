#include "KmerIndex.hpp"

#include <algorithm>
#include <unordered_set>

namespace miint {

bool KmerIndex::encode_kmer(const std::string &seq, size_t pos, uint16_t &out) {
	uint16_t code = 0;
	for (int i = 0; i < K; i++) {
		char c = seq[pos + i];
		uint16_t base;
		switch (c) {
		// Uppercase = unmasked DNA/RNA bases
		case 'A':
			base = 0;
			break;
		case 'C':
			base = 1;
			break;
		case 'G':
			base = 2;
			break;
		case 'T':
		case 'U':
			base = 3;
			break;
		// Lowercase = DUST-masked (low-complexity) bases — skip this k-mer.
		// Also reject ambiguous bases (N, R, Y, etc.).
		default:
			return false;
		}
		code = (code << 2) | base;
	}
	out = code;
	return true;
}

uint32_t KmerIndex::add_sequence(const std::string &sequence) {
	uint32_t seq_id = seq_count_++;

	if (sequence.size() < static_cast<size_t>(K)) {
		return seq_id; // Too short for any k-mers, but ID is still assigned
	}

	// Collect unique k-mer codes for this sequence, then add seq_id once per code.
	// This ensures posting lists contain each seq_id at most once (presence/absence).
	std::unordered_set<uint16_t> seen_codes;
	for (size_t i = 0; i <= sequence.size() - K; i++) {
		uint16_t code;
		if (encode_kmer(sequence, i, code)) {
			seen_codes.insert(code);
		}
	}
	for (uint16_t code : seen_codes) {
		index_[code].push_back(seq_id);
	}

	return seq_id;
}

std::vector<std::pair<size_t, size_t>> KmerIndex::partition_query(size_t query_len) {
	std::vector<std::pair<size_t, size_t>> chunks;
	chunks.reserve(NUM_CHUNKS);
	size_t remaining = query_len;
	size_t cursor = 0;
	for (int i = 0; i < NUM_CHUNKS; i++) {
		size_t parts_left = NUM_CHUNKS - i;
		size_t chunk_len = remaining / parts_left;
		chunks.emplace_back(cursor, chunk_len);
		cursor += chunk_len;
		remaining -= chunk_len;
	}
	return chunks;
}

std::vector<uint32_t> KmerIndex::find_candidates(const std::string &query) const {
	if (query.size() < static_cast<size_t>(K) || seq_count_ == 0) {
		return {};
	}

	auto chunks = partition_query(query.size());

	// Per-chunk top-N selection: only selected candidates get their hits accumulated.
	// A sequence that is rank 5+ in a chunk does not get that chunk's contribution.
	thread_local std::vector<uint32_t> selected_hits;
	selected_hits.assign(seq_count_, 0);

	// Scratch buffer for per-chunk hit counting. Reused across chunks by
	// zeroing only touched entries (not the full array).
	thread_local std::vector<uint32_t> chunk_hits;
	chunk_hits.assign(seq_count_, 0);

	// All seq_ids touched in current chunk (for sparse reset of chunk_hits).
	thread_local std::vector<uint32_t> chunk_touched;

	std::vector<uint32_t> all_candidates;
	// Track all selected seq_ids across chunks (for sparse reset of selected_hits).
	std::vector<uint32_t> all_selected;

	for (auto &[chunk_start, chunk_len] : chunks) {
		// Pre-compute loop bound to avoid size_t underflow when chunk_len < K
		if (chunk_len < static_cast<size_t>(K)) {
			continue;
		}
		size_t loop_end = chunk_start + chunk_len - K;

		chunk_touched.clear();

		// Count unique shared k-mers for this chunk.
		// Since posting lists store each seq_id at most once per k-mer,
		// chunk_hits[sid] counts the number of distinct k-mers shared.
		for (size_t i = chunk_start; i <= loop_end; i++) {
			uint16_t code;
			if (!encode_kmer(query, i, code)) {
				continue;
			}
			for (uint32_t sid : index_[code]) {
				if (chunk_hits[sid] == 0) {
					chunk_touched.push_back(sid);
				}
				chunk_hits[sid]++;
			}
		}

		// Select top TOP_PER_CHUNK from touched entries (no full-array scan)
		size_t select_count = chunk_touched.size();
		if (select_count > static_cast<size_t>(TOP_PER_CHUNK)) {
			std::partial_sort(chunk_touched.begin(), chunk_touched.begin() + TOP_PER_CHUNK, chunk_touched.end(),
			                  [&](uint32_t a, uint32_t b) { return chunk_hits[a] > chunk_hits[b]; });
			select_count = TOP_PER_CHUNK;
		}

		// Accumulate selected candidates' hits into selected_hits
		for (size_t j = 0; j < select_count; j++) {
			uint32_t sid = chunk_touched[j];
			selected_hits[sid] += chunk_hits[sid];
			all_selected.push_back(sid);
		}
		all_candidates.insert(all_candidates.end(), chunk_touched.begin(), chunk_touched.begin() + select_count);

		// Sparse reset: zero ALL touched entries in chunk_hits
		for (uint32_t sid : chunk_touched) {
			chunk_hits[sid] = 0;
		}
	}

	// Deduplicate candidates
	std::sort(all_candidates.begin(), all_candidates.end());
	all_candidates.erase(std::unique(all_candidates.begin(), all_candidates.end()), all_candidates.end());

	// Sort by selected hit count descending, seq_id ascending for deterministic tie-breaking
	std::sort(all_candidates.begin(), all_candidates.end(), [&](uint32_t a, uint32_t b) {
		if (selected_hits[a] != selected_hits[b]) {
			return selected_hits[a] > selected_hits[b];
		}
		return a < b;
	});

	// Cap at MAX_CANDIDATES
	if (all_candidates.size() > static_cast<size_t>(MAX_CANDIDATES)) {
		all_candidates.resize(MAX_CANDIDATES);
	}

	// Sparse reset of selected_hits (only entries we touched)
	for (uint32_t sid : all_selected) {
		selected_hits[sid] = 0;
	}

	return all_candidates;
}

} // namespace miint
