#include "blast_client.hpp"
#include "duckdb/common/exception/http_exception.hpp"
#include "duckdb/common/printer.hpp"
#include <sstream>
#include <thread>

namespace miint {

BlastClient::BlastClient(duckdb::DatabaseInstance &db, const std::string &api_key)
    : db_(db), api_key_(api_key), last_request_time_(std::chrono::steady_clock::now() - std::chrono::seconds(1)) {
}

double BlastClient::GetRateLimit() const {
	return api_key_.empty() ? RATE_LIMIT_NO_KEY : RATE_LIMIT_WITH_KEY;
}

void BlastClient::RespectRateLimit() {
	std::lock_guard<std::mutex> lock(rate_limit_mutex_);
	auto now = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_request_time_);
	auto min_interval = std::chrono::milliseconds(static_cast<int>(1000.0 / GetRateLimit()));
	if (elapsed < min_interval) {
		std::this_thread::sleep_for(min_interval - elapsed);
	}
	last_request_time_ = std::chrono::steady_clock::now();
}

bool BlastClient::IsRetryableStatus(int status) {
	return status == 408 || status == 429 || status == 500 || status == 502 || status == 503 || status == 504;
}

std::string BlastClient::MakeGetRequest(const std::string &url) {
	duckdb::HTTPHeaders headers(db_);
	auto &http_util = duckdb::HTTPUtil::Get(db_);
	auto params = http_util.InitializeParameters(db_, url);
	duckdb::GetRequestInfo get_request(url, headers, *params, nullptr, nullptr);
	get_request.try_request = true;

	int retry_delay_ms = INITIAL_RETRY_DELAY_MS;
	for (int attempt = 0; attempt <= MAX_RETRIES; attempt++) {
		RespectRateLimit();
		auto response = http_util.Request(get_request);
		if (response->Success()) {
			return response->body;
		}
		if (attempt < MAX_RETRIES && !response->HasRequestError() && IsRetryableStatus(int(response->status))) {
			std::this_thread::sleep_for(std::chrono::milliseconds(retry_delay_ms));
			retry_delay_ms *= 2;
			continue;
		}
		if (response->HasRequestError()) {
			throw duckdb::IOException("blast: request failed: %s (URL: %s)", response->GetRequestError(), url);
		}
		throw duckdb::HTTPException(*response, "blast: request failed with HTTP %d (URL: %s)", int(response->status),
		                            url);
	}
	throw duckdb::IOException("blast: request failed after %d retries (URL: %s)", MAX_RETRIES, url);
}

std::string BlastClient::MakePostRequest(const std::string &url, const std::string &body,
                                         const std::string &content_type) {
	duckdb::HTTPHeaders headers(db_);
	headers.Insert("Content-Type", content_type);
	auto &http_util = duckdb::HTTPUtil::Get(db_);
	auto params = http_util.InitializeParameters(db_, url);
	const auto *buf_ptr = reinterpret_cast<duckdb::const_data_ptr_t>(body.data());

	int retry_delay_ms = INITIAL_RETRY_DELAY_MS;
	for (int attempt = 0; attempt <= MAX_RETRIES; attempt++) {
		RespectRateLimit();
		// Rebuild each iteration — PostRequestInfo carries response state in
		// buffer_out (DuckDB issue #19062).
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
			throw duckdb::IOException("blast: POST failed: %s (URL: %s)", response->GetRequestError(), url);
		}
		throw duckdb::HTTPException(*response, "blast: POST failed with HTTP %d (URL: %s)", int(response->status), url);
	}
	throw duckdb::IOException("blast: POST failed after %d retries (URL: %s)", MAX_RETRIES, url);
}

BlastSubmitResult BlastClient::Submit(const std::string &program, const std::string &database,
                                      const std::string &fasta_query, double evalue, int max_targets, bool megablast) {
	std::string body = BlastParser::BuildSubmitBody(program, database, fasta_query, evalue, max_targets, megablast);
	if (!api_key_.empty()) {
		body += "&api_key=" + BlastParser::UrlEncode(api_key_);
	}
	std::string response = MakePostRequest(BLAST_BASE, body, "application/x-www-form-urlencoded");
	auto result = BlastParser::ParseSubmitResponse(response);
	if (!result.error.empty()) {
		throw duckdb::IOException("blast: NCBI returned error: %s", result.error);
	}
	if (result.rid.empty()) {
		auto snippet = response.substr(0, std::min<size_t>(response.size(), 400));
		throw duckdb::IOException("blast: failed to obtain RID from NCBI (response: %s)", snippet);
	}
	return result;
}

void BlastClient::WaitForCompletion(const std::string &rid, int rtoe_hint) {
	if (rid.empty()) {
		throw duckdb::IOException("blast: cannot poll with empty RID");
	}

	std::string encoded_rid = BlastParser::UrlEncode(rid);

	const int max_initial_wait = MAX_POLL_ATTEMPTS * POLL_INTERVAL_SECONDS;
	if (rtoe_hint > max_initial_wait) {
		rtoe_hint = max_initial_wait;
	}
	if (rtoe_hint > 0) {
		duckdb::Printer::Print("blast: waiting " + std::to_string(rtoe_hint) + "s (NCBI estimate) for RID " + rid);
		std::this_thread::sleep_for(std::chrono::seconds(rtoe_hint));
	}

	int unknown_count = 0;
	for (int attempt = 0; attempt < MAX_POLL_ATTEMPTS; attempt++) {
		std::ostringstream url;
		url << BLAST_BASE << "?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" << encoded_rid;
		std::string response = MakeGetRequest(url.str());
		auto status = BlastParser::ParseStatusResponse(response);

		if (status == BlastStatus::READY) {
			return;
		}
		if (status == BlastStatus::UNKNOWN) {
			unknown_count++;
			if (unknown_count > MAX_UNKNOWN_STATUS_RETRIES) {
				auto snippet = response.substr(0, std::min<size_t>(response.size(), 400));
				throw duckdb::IOException(
				    "blast: job %s returned UNKNOWN status %d times — may have expired or failed (response: %s)", rid,
				    unknown_count, snippet);
			}
			duckdb::Printer::Print("blast: RID " + rid + " returned UNKNOWN status (attempt " +
			                       std::to_string(unknown_count) + "/" + std::to_string(MAX_UNKNOWN_STATUS_RETRIES) +
			                       "), retrying...");
		} else {
			unknown_count = 0;
		}

		int elapsed = rtoe_hint + (attempt + 1) * POLL_INTERVAL_SECONDS;
		if (status == BlastStatus::WAITING) {
			duckdb::Printer::Print("blast: RID " + rid + " still running (~" + std::to_string(elapsed) + "s elapsed)");
		}
		std::this_thread::sleep_for(std::chrono::seconds(POLL_INTERVAL_SECONDS));
	}

	throw duckdb::IOException("blast: timed out after %d poll attempts (~%d minutes) for RID %s. "
	                          "Check status at: %s?CMD=Get&RID=%s",
	                          MAX_POLL_ATTEMPTS, (MAX_POLL_ATTEMPTS * POLL_INTERVAL_SECONDS) / 60, rid, BLAST_BASE,
	                          rid);
}

std::string BlastClient::RetrieveResults(const std::string &rid, int max_targets) {
	std::string encoded_rid = BlastParser::UrlEncode(rid);
	std::ostringstream url;
	url << BLAST_BASE << "?CMD=Get&FORMAT_TYPE=Text&ALIGNMENT_VIEW=Tabular&HITLIST_SIZE=" << max_targets
	    << "&RID=" << encoded_rid;
	return MakeGetRequest(url.str());
}

} // namespace miint
