#include "ncbi_client.hpp"
#include "duckdb/common/exception/http_exception.hpp"
#include "miniz.hpp"
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
	// Retry on rate limiting (429) and service unavailable (503)
	// Also retry on 500 (internal server error) and 502 (bad gateway)
	return status == 429 || status == 500 || status == 502 || status == 503;
}

std::string NCBIClient::MakeRequest(const std::string &url, bool use_api_key_header) {
	RespectRateLimit();

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

std::string NCBIClient::FetchGenBankXML(const std::string &accession) {
	std::ostringstream url;
	url << EUTILS_BASE << "/efetch.fcgi?db=nuccore&id=" << accession << "&rettype=gb&retmode=xml";
	if (!api_key.empty()) {
		url << "&api_key=" << api_key;
	}
	return MakeRequest(url.str(), false);
}

std::string NCBIClient::FetchFasta(const std::string &accession) {
	std::ostringstream url;
	url << EUTILS_BASE << "/efetch.fcgi?db=nuccore&id=" << accession << "&rettype=fasta&retmode=text";
	if (!api_key.empty()) {
		url << "&api_key=" << api_key;
	}
	return MakeRequest(url.str(), false);
}

std::string NCBIClient::FetchFeatureTable(const std::string &accession) {
	std::ostringstream url;
	url << EUTILS_BASE << "/efetch.fcgi?db=nuccore&id=" << accession << "&rettype=ft&retmode=text";
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

std::string NCBIClient::ExtractFromZip(const std::string &zip_data, const std::string &pattern) {
	duckdb_miniz::mz_zip_archive zip;
	memset(&zip, 0, sizeof(zip));

	// Initialize ZIP reader from memory
	if (!duckdb_miniz::mz_zip_reader_init_mem(&zip, zip_data.data(), zip_data.size(), 0)) {
		throw duckdb::IOException("Failed to read ZIP archive from Datasets API");
	}

	// RAII guard to ensure ZIP reader cleanup on exception or return
	auto cleanup = [&zip](void *) {
		duckdb_miniz::mz_zip_reader_end(&zip);
	};
	std::unique_ptr<void, decltype(cleanup)> guard(reinterpret_cast<void *>(1), cleanup);

	std::string result;

	// Iterate through files in ZIP
	duckdb_miniz::mz_uint num_files = duckdb_miniz::mz_zip_reader_get_num_files(&zip);
	for (duckdb_miniz::mz_uint i = 0; i < num_files; i++) {
		duckdb_miniz::mz_zip_archive_file_stat file_stat;
		if (!duckdb_miniz::mz_zip_reader_file_stat(&zip, i, &file_stat)) {
			continue;
		}

		std::string filename = file_stat.m_filename;

		// Check if filename matches pattern (e.g., "*_genomic.fna")
		if (filename.find(pattern) != std::string::npos) {
			// Extract file to memory
			size_t uncompressed_size;
			void *data = duckdb_miniz::mz_zip_reader_extract_to_heap(&zip, i, &uncompressed_size, 0);
			if (!data) {
				throw duckdb::IOException("Failed to extract '%s' from ZIP: %s", filename.c_str(),
				                          duckdb_miniz::mz_zip_get_error_string(zip.m_last_error));
			}
			result.assign(static_cast<char *>(data), uncompressed_size);
			free(data); // Use standard free(), not mz_free() - allocated by pZip->m_pAlloc (malloc)
			break;
		}
	}

	if (result.empty()) {
		throw duckdb::IOException("No FASTA file found in assembly package");
	}

	return result;
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
	if (zip_data.size() < 4 || zip_data[0] != 'P' || zip_data[1] != 'K' || zip_data[2] != 0x03 ||
	    zip_data[3] != 0x04) {
		throw duckdb::IOException("Datasets API returned invalid ZIP file for assembly '%s'. "
		                          "Response starts with: %s",
		                          accession.c_str(),
		                          zip_data.substr(0, std::min<size_t>(100, zip_data.size())).c_str());
	}

	// Extract *_genomic.fna from ZIP using miniz
	return ExtractFromZip(zip_data, "_genomic.fna");
}

} // namespace miint
