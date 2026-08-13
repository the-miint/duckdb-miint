#include "ena_resolver_cache.hpp"
#include "ena_parser.hpp"
#include <algorithm>
#include <stdexcept>
#include <unordered_set>

namespace miint {

ENAResolverCache::ENAResolverCache(size_t capacity_p) : capacity(capacity_p) {
}

bool ENAResolverCache::Get(const Key &key, std::vector<ENARunInfo> &out) {
	if (capacity == 0) {
		return false;
	}
	std::lock_guard<std::mutex> lock(m);
	auto it = map.find(key);
	if (it == map.end()) {
		return false;
	}
	// Bump to MRU
	order.erase(it->second.second);
	order.push_front(key);
	it->second.second = order.begin();
	out = it->second.first;
	return true;
}

void ENAResolverCache::Put(const Key &key, std::vector<ENARunInfo> value) {
	if (capacity == 0) {
		return;
	}
	std::lock_guard<std::mutex> lock(m);
	auto it = map.find(key);
	if (it != map.end()) {
		// Replace existing value; bump to MRU
		order.erase(it->second.second);
		order.push_front(key);
		it->second.first = std::move(value);
		it->second.second = order.begin();
		return;
	}
	// Evict LRU if at capacity
	if (map.size() >= capacity) {
		const auto &evict_key = order.back();
		map.erase(evict_key);
		order.pop_back();
	}
	order.push_front(key);
	map.emplace(key, std::make_pair(std::move(value), order.begin()));
}

size_t ENAResolverCache::Size() {
	std::lock_guard<std::mutex> lock(m);
	return map.size();
}

// ---- ResolveRunsBatch ----

// Field list we request from ENA for every batched read_run query.
// Includes run_accession plus all other columns ENARunInfoExtractor consumes, plus
// study_accession so callers can tell which input accession a returned row belongs to
// when querying by study.
static constexpr const char *BATCH_FIELDS = "run_accession,sample_accession,experiment_accession,study_accession,"
                                            "fastq_ftp,fastq_aspera,fastq_bytes,fastq_md5,library_layout,"
                                            "submitted_ftp,submitted_aspera,submitted_bytes,submitted_format";

// Column name on a returned TSV row that matches the input accession for a given type.
// STUDY-typed inputs are matched against the row's study_accession, etc.
static int GroupingColumn(ENAAccessionType type, const ENAColumnIndexMap &cols) {
	switch (type) {
	case ENAAccessionType::STUDY:
		return cols.study_accession;
	case ENAAccessionType::SAMPLE:
		return cols.sample_accession;
	case ENAAccessionType::RUN:
		return cols.run_accession;
	case ENAAccessionType::EXPERIMENT:
		return cols.experiment_accession;
	default:
		return -1;
	}
}

// Fetch + parse + group one batch of accessions of the same type. Returns a map keyed
// by the input accession (matched against the TSV row's grouping column).
static std::unordered_map<std::string, std::vector<ENARunInfo>> FetchBatch(const ENAFetcher &fetcher,
                                                                           const std::vector<std::string> &batch,
                                                                           ENAAccessionType type,
                                                                           const std::string &prefer_format) {
	std::unordered_map<std::string, std::vector<ENARunInfo>> out;
	for (const auto &acc : batch) {
		out[acc]; // ensure every input accession appears, even if no rows match
	}

	auto url = ENAParser::BuildSearchURLBatch(batch, type, "read_run", BATCH_FIELDS);
	auto tsv = fetcher(url);
	auto parsed = ENAParser::ParseTSV(tsv);
	if (parsed.column_names.empty()) {
		// Header-only (or truly empty) TSV — nothing to group. Leave the pre-populated empty
		// entries in place so callers negative-cache every requested accession.
		return out;
	}
	auto cols = ENAColumnIndexMap::FromHeader(parsed.column_names);
	int group_col = GroupingColumn(type, cols);
	if (group_col < 0) {
		throw std::runtime_error("ENA batched response is missing the expected grouping column for the requested "
		                         "accession type (unexpected API response or server error)");
	}

	for (const auto &row : parsed.rows) {
		auto row_runs = ENARunInfoExtractor::FromTSVRow(row, cols, prefer_format);
		if (row_runs.empty()) {
			continue;
		}
		std::string grouping_key = ENAColumnIndexMap::Get(row, group_col);
		if (grouping_key.empty()) {
			continue;
		}
		auto it = out.find(grouping_key);
		if (it == out.end()) {
			// Row's grouping key doesn't match any requested accession (unexpected ENA
			// response). Skip silently — do NOT cache under an unrequested key, since
			// that would consume a cache slot and could create a false-positive hit in
			// a future query for that accession.
			continue;
		}
		for (auto &info : row_runs) {
			it->second.push_back(std::move(info));
		}
	}
	return out;
}

std::unordered_map<std::string, std::vector<ENARunInfo>>
ResolveRunsBatch(const ENAFetcher &fetcher, ENAResolverCache &cache, const std::vector<std::string> &accessions,
                 const std::string &prefer_format, size_t max_batch_size) {
	std::unordered_map<std::string, std::vector<ENARunInfo>> result;
	if (max_batch_size == 0) {
		max_batch_size = 1;
	}

	// Preserve input order for callers that care, while deduping.
	std::unordered_set<std::string> seen;
	std::vector<std::string> unique_inputs;
	unique_inputs.reserve(accessions.size());
	for (const auto &acc : accessions) {
		if (seen.insert(acc).second) {
			unique_inputs.push_back(acc);
		}
	}

	// Step 1: serve cache hits; bucket misses by accession type.
	std::unordered_map<ENAAccessionType, std::vector<std::string>> by_type;
	for (const auto &acc : unique_inputs) {
		ENAResolverCache::Key key {acc, prefer_format};
		std::vector<ENARunInfo> cached;
		if (cache.Get(key, cached)) {
			result[acc] = std::move(cached);
			continue;
		}
		by_type[ENAParser::DetectAccessionType(acc)].push_back(acc);
	}

	// Step 2: fetch typed groups in batches.
	for (auto &kv : by_type) {
		if (kv.first == ENAAccessionType::UNKNOWN) {
			// UNKNOWN type: we cannot build a compound query (no known column). Fall back to
			// per-accession single-URL fetches.
			for (const auto &acc : kv.second) {
				auto url = ENAParser::BuildSearchURL(acc, "read_run", BATCH_FIELDS);
				auto tsv = fetcher(url);
				auto parsed = ENAParser::ParseTSV(tsv);
				auto cols = ENAColumnIndexMap::FromHeader(parsed.column_names);
				std::vector<ENARunInfo> acc_runs;
				for (const auto &row : parsed.rows) {
					auto row_runs = ENARunInfoExtractor::FromTSVRow(row, cols, prefer_format);
					for (auto &info : row_runs) {
						acc_runs.push_back(std::move(info));
					}
				}
				cache.Put({acc, prefer_format}, acc_runs);
				result[acc] = std::move(acc_runs);
			}
			continue;
		}

		const auto &all = kv.second;
		for (size_t start = 0; start < all.size(); start += max_batch_size) {
			size_t end = std::min(start + max_batch_size, all.size());
			std::vector<std::string> batch(all.begin() + start, all.begin() + end);
			auto batch_result = FetchBatch(fetcher, batch, kv.first, prefer_format);
			for (const auto &acc : batch) {
				auto it = batch_result.find(acc);
				auto runs = (it == batch_result.end()) ? std::vector<ENARunInfo> {} : std::move(it->second);
				cache.Put({acc, prefer_format}, runs); // negative cache when runs is empty
				result[acc] = std::move(runs);
			}
		}
	}

	return result;
}

} // namespace miint
