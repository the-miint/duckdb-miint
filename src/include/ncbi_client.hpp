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
// - Retries on HTTP 400, 408, 429, 500, 502, 503, 504. NCBI E-utilities occasionally
//   emits transient 400/408 under load — same URL succeeds on retry — so treat these
//   as retryable rather than aborting a long scan.
// - Uses exponential backoff: 1s, 2s, 4s, 8s, 16s, 32s (max 6 retries, ~63s total)
// - Other 4xx (e.g. 404 for an invalid accession) fail immediately
class NCBIClient {
public:
	// Base URLs
	static constexpr const char *EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";
	static constexpr const char *DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2";

	// Rate limits (requests per second)
	static constexpr double RATE_LIMIT_NO_KEY = 3.0;
	static constexpr double RATE_LIMIT_WITH_KEY = 10.0;

	// Retry configuration
	// 6 retries with exponential backoff (1s, 2s, 4s, 8s, 16s, 32s) gives ~63s of
	// total backoff before giving up — enough headroom for transient NCBI hiccups
	// during long scans (e.g. 10k accessions) without making a genuinely broken
	// request hang forever.
	static constexpr int MAX_RETRIES = 6;
	static constexpr int INITIAL_RETRY_DELAY_MS = 1000;

	NCBIClient(duckdb::DatabaseInstance &db, const std::string &api_key = "");

	// Accession type detection - delegates to NCBIParser
	static AccessionType DetectAccessionType(const std::string &accession) {
		return NCBIParser::DetectAccessionType(accession);
	}
	static bool IsAssemblyAccession(const std::string &accession) {
		return NCBIParser::IsAssemblyAccession(accession);
	}

	// E-utilities single-accession methods. read_ncbi_fasta and read_ncbi route
	// through the batched epost+efetch path below; only read_ncbi_annotation
	// still uses the single-fetch endpoint because feature-table batching is
	// a separate code path (different rettype, different parser).
	std::string FetchFeatureTable(const std::string &accession);

	// E-utilities batched methods. epost stages the ID list on NCBI's server and
	// returns a WebEnv+QueryKey handle; efetch then reads all records in one call.
	// This collapses an N-call sequence into 2 calls per batch and is the NCBI-
	// recommended path for >1 accession.
	// `db` selects the Entrez database (default "nuccore"; "taxonomy" for taxid lineages).
	EPostResult EPostIds(const std::vector<std::string> &accessions, const std::string &db = "nuccore");

	// Batched fetch of FASTA / GenBank XML. Empty input is a no-op (returns "") so
	// callers can pass already-partitioned lists (e.g. sequences vs. assemblies)
	// without guarding the empty case.
	std::string FetchFastaBatch(const std::vector<std::string> &accessions);
	std::string FetchGenBankXMLBatch(const std::vector<std::string> &accessions);

	// Batched fetch of taxonomy XML (db=taxonomy, retmode=xml) for a list of taxids.
	// Empty input is a no-op (returns "").
	std::string FetchTaxonomyXMLBatch(const std::vector<std::string> &taxids);

	// Datasets API methods
	std::string FetchAssemblyReport(const std::string &accession); // JSON metadata
	std::string FetchAssemblyFasta(const std::string &accession);  // ZIP download, extract FASTA

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

	// POST counterpart of MakeRequest. PostRequestInfo carries response state in
	// its buffer_out field (DuckDB issue #19062), so the request struct must be
	// rebuilt inside the retry loop rather than reused — see ena_client.cpp for
	// the precedent this mirrors.
	std::string MakeRequestPOST(const std::string &url, const std::string &body, const std::string &content_type);

	// Check if HTTP status code is retryable
	static bool IsRetryableStatus(int status);

	// Rate limiting (thread-safe)
	void RespectRateLimit();
	double GetRateLimit() const;
};

} // namespace miint
