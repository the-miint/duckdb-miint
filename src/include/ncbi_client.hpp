#pragma once

#include "ncbi_parser.hpp"
#include "duckdb/common/http_util.hpp"
#include "duckdb/main/database.hpp"
#include <chrono>
#include <mutex>
#include <string>

namespace miint {

// HTTP client wrapper for NCBI APIs
// Handles both E-utilities and Datasets API v2
// For parsing utilities, see NCBIParser in ncbi_parser.hpp
//
// Rate Limiting Behavior:
// - Without API key: Limited to 3 requests/second
// - With API key: Limited to 10 requests/second
// - Rate limiting is enforced per NCBIClient instance
// - Thread-safe: Multiple threads can share a single client instance
//
// Retry Behavior:
// - Retries on HTTP 429 (Too Many Requests) and 503 (Service Unavailable)
// - Uses exponential backoff: 1s, 2s, 4s (max 3 retries)
// - Non-retryable errors (4xx except 429) fail immediately
class NCBIClient {
public:
	// Base URLs
	static constexpr const char *EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";
	static constexpr const char *DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2";

	// Rate limits (requests per second)
	static constexpr double RATE_LIMIT_NO_KEY = 3.0;
	static constexpr double RATE_LIMIT_WITH_KEY = 10.0;

	// Retry configuration
	static constexpr int MAX_RETRIES = 3;
	static constexpr int INITIAL_RETRY_DELAY_MS = 1000;

	NCBIClient(duckdb::DatabaseInstance &db, const std::string &api_key = "");

	// Accession type detection - delegates to NCBIParser
	static AccessionType DetectAccessionType(const std::string &accession) {
		return NCBIParser::DetectAccessionType(accession);
	}
	static bool IsAssemblyAccession(const std::string &accession) {
		return NCBIParser::IsAssemblyAccession(accession);
	}

	// E-utilities methods (primary) - return raw response strings
	std::string FetchGenBankXML(const std::string &accession);
	std::string FetchFasta(const std::string &accession);
	std::string FetchFeatureTable(const std::string &accession);

	// Datasets API (assembly metadata only, JSON response)
	std::string FetchAssemblyReport(const std::string &accession);

	// Parsing methods - delegate to NCBIParser
	static GenBankMetadata ParseGenBankXML(const std::string &xml) {
		return NCBIParser::ParseGenBankXML(xml);
	}
	static SequenceRecordBatch ParseFasta(const std::string &fasta_text) {
		return NCBIParser::ParseFasta(fasta_text);
	}

private:
	duckdb::DatabaseInstance &db;
	std::string api_key;
	std::chrono::steady_clock::time_point last_request_time;
	mutable std::mutex rate_limit_mutex; // Protects last_request_time

	// HTTP request helper with retry logic
	std::string MakeRequest(const std::string &url, bool use_api_key_header = false);

	// Check if HTTP status code is retryable
	static bool IsRetryableStatus(int status);

	// Rate limiting (thread-safe)
	void RespectRateLimit();
	double GetRateLimit() const;
};

} // namespace miint
