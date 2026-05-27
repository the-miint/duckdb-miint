#pragma once

#include "blast_parser.hpp"
#include "duckdb/common/http_util.hpp"
#include "duckdb/main/database.hpp"
#include <chrono>
#include <mutex>
#include <string>

namespace miint {

// HTTP client for the NCBI BLAST web API.
//
// The BLAST API is two-phase:
//   1. Submit (POST CMD=Put) → returns RID + estimated wait time
//   2. Poll (GET CMD=Get&FORMAT_OBJECT=SearchInfo) → Status=WAITING/READY/UNKNOWN
//   3. Retrieve (GET CMD=Get&FORMAT_TYPE=Text&ALIGNMENT_VIEW=Tabular) → outfmt 6 results
//
// Rate limiting: NCBI requests ≤3 req/s. The poll loop uses 15s intervals
// which naturally satisfies this. Individual requests also respect the
// per-request rate limit.
//
// Independent from NCBIClient — uses DuckDB's HTTPUtil directly.
class BlastClient {
public:
	static constexpr const char *BLAST_BASE = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi";
	static constexpr int POLL_INTERVAL_SECONDS = 15;
	static constexpr int MAX_POLL_ATTEMPTS = 240; // 60 min timeout
	static constexpr double RATE_LIMIT_NO_KEY = 3.0;
	static constexpr double RATE_LIMIT_WITH_KEY = 10.0;
	static constexpr int MAX_RETRIES = 6;
	static constexpr int INITIAL_RETRY_DELAY_MS = 1000;
	static constexpr int MAX_UNKNOWN_STATUS_RETRIES = 3;

	BlastClient(duckdb::DatabaseInstance &db, const std::string &api_key = "");
	BlastClient(const BlastClient &) = delete;
	BlastClient &operator=(const BlastClient &) = delete;

	BlastSubmitResult Submit(const std::string &program, const std::string &database, const std::string &fasta_query,
	                         double evalue, int max_targets, bool megablast);

	void WaitForCompletion(const std::string &rid, int rtoe_hint);

	std::string RetrieveResults(const std::string &rid, int max_targets);

private:
	duckdb::DatabaseInstance &db_;
	std::string api_key_;
	std::chrono::steady_clock::time_point last_request_time_;
	mutable std::mutex rate_limit_mutex_;

	std::string MakeGetRequest(const std::string &url);
	std::string MakePostRequest(const std::string &url, const std::string &body, const std::string &content_type);
	void RespectRateLimit();
	double GetRateLimit() const;
	static bool IsRetryableStatus(int status);
};

} // namespace miint
