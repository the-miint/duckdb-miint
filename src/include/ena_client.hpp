#pragma once

#include "ena_parser.hpp"
#include "duckdb/common/http_util.hpp"
#include "duckdb/main/database.hpp"
#include <chrono>
#include <mutex>
#include <string>
#include <vector>

namespace miint {

// HTTP client for EBI/ENA APIs.
// For parsing utilities, use ENAParser directly (ena_parser.hpp).
//
// Rate Limiting: ~3 requests/second (no API key required)
// Retry: HTTP 429, 500, 502, 503 with exponential backoff
class ENAClient {
public:
	static constexpr double RATE_LIMIT = 3.0;
	static constexpr int MAX_RETRIES = 3;
	static constexpr int INITIAL_RETRY_DELAY_MS = 1000;

	ENAClient(duckdb::DatabaseInstance &db);

	// HTTP methods — build URL via ENAParser, then fetch
	std::string Search(const std::string &accession, const std::string &result_type, const std::string &fields);
	std::string FetchXML(const std::vector<std::string> &accessions);

private:
	duckdb::DatabaseInstance &db;
	std::chrono::steady_clock::time_point last_request_time;
	mutable std::mutex rate_limit_mutex;

	std::string MakeRequest(const std::string &url);
	static bool IsRetryableStatus(int status);
	void RespectRateLimit();
};

} // namespace miint
