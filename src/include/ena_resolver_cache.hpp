#pragma once

#include "ena_run_info_extractor.hpp"
#include <cstddef>
#include <functional>
#include <list>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace miint {

// Callable that performs a single HTTP GET and returns the response body.
// Injectable so ResolveRunsBatch can be unit-tested without hitting the network.
using ENAFetcher = std::function<std::string(const std::string &url)>;

// Thread-safe LRU cache keyed on (accession, prefer_format).
// Value is the list of RunInfo records returned by ENA for that accession.
// Holds negative results (empty vectors) so repeated lookups of known-missing
// accessions don't re-hit the network.
class ENAResolverCache {
public:
	struct Key {
		std::string accession;
		std::string prefer_format;

		bool operator==(const Key &other) const {
			return accession == other.accession && prefer_format == other.prefer_format;
		}
	};

	struct KeyHash {
		size_t operator()(const Key &k) const noexcept {
			auto h1 = std::hash<std::string> {}(k.accession);
			auto h2 = std::hash<std::string> {}(k.prefer_format);
			return h1 ^ (h2 * 0x9e3779b97f4a7c15ULL);
		}
	};

	// capacity = 0 disables caching (Put is a no-op, Get always misses).
	explicit ENAResolverCache(size_t capacity = 256);

	// Returns true on hit; copies the stored vector into 'out'.
	bool Get(const Key &key, std::vector<ENARunInfo> &out);

	// Inserts or replaces. Evicts the least-recently-used entry if at capacity.
	void Put(const Key &key, std::vector<ENARunInfo> value);

	size_t Size();

private:
	size_t capacity;
	std::mutex m;
	std::list<Key> order; // MRU at front, LRU at back
	std::unordered_map<Key, std::pair<std::vector<ENARunInfo>, std::list<Key>::iterator>, KeyHash> map;
};

// Resolve one-or-more accessions to RunInfo records using batched ENA Portal queries.
//
// Behavior:
//   - Deduplicates input accessions.
//   - Returns cache hits without hitting the fetcher.
//   - Groups cache-miss accessions by detected accession type (RUN, STUDY, etc.)
//     and issues one compound "<col> IN (...)" query per group, splitting into
//     batches of at most max_batch_size to respect URL-length limits.
//   - UNKNOWN-type accessions fall back to single-accession queries (one fetch each).
//   - Negative results (zero rows returned for an accession) are cached to avoid
//     re-fetching.
//   - Every input accession appears as a key in the output map. The value is empty
//     when no runs were found.
//
// fetcher: called with a fully-formed URL; returns the TSV body.
// cache:   mutated in-place.
//
// Thread-safety / concurrency note:
//   There is a TOCTOU window between the per-accession Get(miss) and the later
//   Put(result). Two threads that both miss on the same accession will issue two
//   fetches and the second Put overwrites the first. The final result is idempotent,
//   so correctness is preserved, but the second fetch wastes an HTTP round-trip.
//   This is acceptable for Phase A because ResolveRuns is only invoked from Bind()
//   on a single thread. Phase B (lateral in_out_function) will call this from worker
//   threads — if workers can share a cache and query overlapping accessions, consider
//   adding a per-key in-flight map guarded by the same mutex so the first miss blocks
//   the duplicates until the fetch completes.
std::unordered_map<std::string, std::vector<ENARunInfo>>
ResolveRunsBatch(const ENAFetcher &fetcher, ENAResolverCache &cache, const std::vector<std::string> &accessions,
                 const std::string &prefer_format, size_t max_batch_size = 50);

} // namespace miint
