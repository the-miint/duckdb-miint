#include "ncbi_client.hpp"
#include "zip_utils.hpp"
#include "duckdb/common/exception/http_exception.hpp"
#include <sstream>
#include <thread>

namespace miint {

NCBIClient::NCBIClient(duckdb::DatabaseInstance &db, const std::string &api_key)
    : db(db), api_key(api_key), last_request_time(std::chrono::steady_clock::now() - std::chrono::seconds(1)) {
}

double NCBIClient::GetRateLimit() const {
	return api_key.empty() ? RATE_LIMIT_NO_KEY : RATE_LIMIT_WITH_KEY;
}

void NCBIClient::RespectRateLimit() {
	std::lock_guard<std::mutex> lock(rate_limit_mutex);

	auto now = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_request_time);

	// Calculate minimum interval between requests
	double rate_limit = GetRateLimit();
	auto min_interval = std::chrono::milliseconds(static_cast<int>(1000.0 / rate_limit));

	if (elapsed < min_interval) {
		std::this_thread::sleep_for(min_interval - elapsed);
	}

	last_request_time = std::chrono::steady_clock::now();
}

bool NCBIClient::IsRetryableStatus(int status) {
	// Retry on rate limiting (429) and service unavailable (503).
	// Also retry on 500 (internal server error), 502 (bad gateway), 504 (gateway timeout).
	// NCBI E-utilities occasionally emits transient 400/408 under load — the same URL
	// succeeds on retry, so treat these as retryable too rather than aborting a long scan.
	return status == 400 || status == 408 || status == 429 || status == 500 || status == 502 || status == 503 ||
	       status == 504;
}

std::string NCBIClient::MakeRequest(const std::string &url, bool use_api_key_header) {
	duckdb::HTTPHeaders headers(db);

	// For Datasets API, use header for API key
	if (use_api_key_header && !api_key.empty()) {
		headers.Insert("api-key", api_key);
	}

	auto &http_util = duckdb::HTTPUtil::Get(db);
	auto params = http_util.InitializeParameters(db, url);

	duckdb::GetRequestInfo get_request(url, headers, *params, nullptr, nullptr);
	get_request.try_request = true;

	int retry_delay_ms = INITIAL_RETRY_DELAY_MS;

	for (int attempt = 0; attempt <= MAX_RETRIES; attempt++) {
		// Rate-limit each attempt, not just the first: after backing off for a
		// 429/5xx retry we still need to honour the per-request cadence so the
		// next retry doesn't fire faster than NCBI's quota allows.
		RespectRateLimit();
		auto response = http_util.Request(get_request);

		if (response->Success()) {
			return response->body;
		}

		// Check if we should retry
		if (attempt < MAX_RETRIES && !response->HasRequestError() && IsRetryableStatus(int(response->status))) {
			std::this_thread::sleep_for(std::chrono::milliseconds(retry_delay_ms));
			retry_delay_ms *= 2; // Exponential backoff
			continue;
		}

		// Non-retryable error or max retries exceeded
		if (response->HasRequestError()) {
			throw duckdb::IOException("NCBI request failed: %s (URL: %s)", response->GetRequestError(), url);
		}
		throw duckdb::HTTPException(*response, "NCBI request failed with HTTP %d (URL: %s)", int(response->status),
		                            url);
	}

	// Should not reach here, but just in case
	throw duckdb::IOException("NCBI request failed after %d retries (URL: %s)", MAX_RETRIES, url);
}

std::string NCBIClient::FetchFeatureTable(const std::string &accession) {
	std::ostringstream url;
	url << EUTILS_BASE << "/efetch.fcgi?db=nuccore&id=" << accession << "&rettype=ft&retmode=text";
	if (!api_key.empty()) {
		url << "&api_key=" << api_key;
	}
	return MakeRequest(url.str(), false);
}

std::string NCBIClient::MakeRequestPOST(const std::string &url, const std::string &body,
                                        const std::string &content_type) {
	duckdb::HTTPHeaders headers(db);
	headers.Insert("Content-Type", content_type);

	auto &http_util = duckdb::HTTPUtil::Get(db);
	auto params = http_util.InitializeParameters(db, url);

	const auto *buf_ptr = reinterpret_cast<duckdb::const_data_ptr_t>(body.data());

	int retry_delay_ms = INITIAL_RETRY_DELAY_MS;

	// PostRequestInfo carries response state in `buffer_out` (DuckDB issue
	// duckdb/duckdb#19062), so rebuild it each iteration rather than reusing a
	// single instance — a retried POST must start with a clean output buffer.
	// This mirrors ena_client.cpp:PostBody, which hit the same problem first.
	// Unlike ENA (which intentionally skips rate-limiting on its write API),
	// NCBI's epost serves the read endpoint and must honour the per-request
	// cadence on every attempt, not just the first.
	for (int attempt = 0; attempt <= MAX_RETRIES; attempt++) {
		RespectRateLimit();
		duckdb::PostRequestInfo post_request(url, headers, *params, buf_ptr, body.size());
		post_request.try_request = true;
		auto response = http_util.Request(post_request);

		if (response->Success()) {
			return post_request.buffer_out;
		}

		if (attempt < MAX_RETRIES && !response->HasRequestError() && IsRetryableStatus(int(response->status))) {
			std::this_thread::sleep_for(std::chrono::milliseconds(retry_delay_ms));
			retry_delay_ms *= 2;
			continue;
		}

		if (response->HasRequestError()) {
			throw duckdb::IOException("NCBI POST failed: %s (URL: %s)", response->GetRequestError(), url);
		}
		throw duckdb::HTTPException(*response, "NCBI POST failed with HTTP %d (URL: %s)", int(response->status), url);
	}

	throw duckdb::IOException("NCBI POST failed after %d retries (URL: %s)", MAX_RETRIES, url);
}

EPostResult NCBIClient::EPostIds(const std::vector<std::string> &accessions, const std::string &db) {
	if (accessions.empty()) {
		throw duckdb::IOException("NCBIClient::EPostIds called with empty accession list");
	}

	std::ostringstream body;
	body << "db=" << db << "&id=";
	for (size_t i = 0; i < accessions.size(); i++) {
		if (i > 0) {
			body << ",";
		}
		body << accessions[i];
	}
	if (!api_key.empty()) {
		body << "&api_key=" << api_key;
	}

	std::string url = std::string(EUTILS_BASE) + "/epost.fcgi";
	std::string response = MakeRequestPOST(url, body.str(), "application/x-www-form-urlencoded");
	auto result = NCBIParser::ParseEPostResponse(response);

	if (result.webenv.empty() || result.query_key.empty()) {
		// Surface the server's <ERROR> verbatim if present; otherwise embed the
		// raw response so the failure is debuggable from a single log line
		// (Rule 10: fail loud with the actual server message, not a stock one).
		auto err = NCBIParser::ExtractXMLTagValue(response, "ERROR");
		auto snippet = response.substr(0, std::min<size_t>(response.size(), 400));
		if (!err.empty()) {
			throw duckdb::IOException("NCBI epost returned ERROR: %s (response: %s)", err.c_str(), snippet.c_str());
		}
		throw duckdb::IOException("NCBI epost response missing <WebEnv> or <QueryKey> (response: %s)", snippet.c_str());
	}

	return result;
}

std::string NCBIClient::FetchFastaBatch(const std::vector<std::string> &accessions) {
	if (accessions.empty()) {
		return ""; // No-op so callers needn't guard the empty partition case.
	}
	auto post = EPostIds(accessions);
	std::ostringstream url;
	url << EUTILS_BASE << "/efetch.fcgi?db=nuccore&query_key=" << post.query_key << "&WebEnv=" << post.webenv
	    << "&rettype=fasta&retmode=text";
	if (!api_key.empty()) {
		url << "&api_key=" << api_key;
	}
	return MakeRequest(url.str(), false);
}

std::string NCBIClient::FetchGenBankXMLBatch(const std::vector<std::string> &accessions) {
	if (accessions.empty()) {
		return "";
	}
	auto post = EPostIds(accessions);
	std::ostringstream url;
	url << EUTILS_BASE << "/efetch.fcgi?db=nuccore&query_key=" << post.query_key << "&WebEnv=" << post.webenv
	    << "&rettype=gb&retmode=xml";
	if (!api_key.empty()) {
		url << "&api_key=" << api_key;
	}
	return MakeRequest(url.str(), false);
}

std::string NCBIClient::FetchTaxonomyXMLBatch(const std::vector<std::string> &taxids) {
	if (taxids.empty()) {
		return "";
	}
	auto post = EPostIds(taxids, "taxonomy");
	std::ostringstream url;
	url << EUTILS_BASE << "/efetch.fcgi?db=taxonomy&query_key=" << post.query_key << "&WebEnv=" << post.webenv
	    << "&retmode=xml";
	if (!api_key.empty()) {
		url << "&api_key=" << api_key;
	}
	return MakeRequest(url.str(), false);
}

std::string NCBIClient::FetchAssemblyReport(const std::string &accession) {
	std::ostringstream url;
	url << DATASETS_BASE << "/genome/accession/" << accession << "/dataset_report";
	return MakeRequest(url.str(), true); // Use API key header for Datasets API
}

// NOTE: This loads the entire ZIP file and extracted FASTA into memory simultaneously.
// For large assemblies (e.g., human genome ~3GB), this will use ~6GB RAM.
// TODO: Implement streaming extraction for memory efficiency if needed.
std::string NCBIClient::FetchAssemblyFasta(const std::string &accession) {
	// Build download URL for assembly genome FASTA
	std::ostringstream url;
	url << DATASETS_BASE << "/genome/accession/" << accession << "/download?include_annotation_type=GENOME_FASTA";

	// Download ZIP file as binary data
	std::string zip_data = MakeRequest(url.str(), true); // Use API key header for Datasets API

	// Validate ZIP magic bytes (PK\x03\x04)
	if (zip_data.size() < 4 || zip_data[0] != 'P' || zip_data[1] != 'K' || zip_data[2] != 0x03 || zip_data[3] != 0x04) {
		throw duckdb::IOException("Datasets API returned invalid ZIP file for assembly '%s'. "
		                          "Response starts with: %s",
		                          accession.c_str(),
		                          zip_data.substr(0, std::min<size_t>(100, zip_data.size())).c_str());
	}

	// Extract *_genomic.fna from ZIP
	try {
		return ExtractFromZip(zip_data, "_genomic.fna");
	} catch (const std::runtime_error &e) {
		throw duckdb::IOException("Failed to extract assembly FASTA from ZIP: %s", e.what());
	}
}

} // namespace miint
