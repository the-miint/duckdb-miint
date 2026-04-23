#include "ena_client.hpp"
#include "duckdb/common/exception/http_exception.hpp"
#include <thread>

namespace miint {

ENAClient::ENAClient(duckdb::DatabaseInstance &db)
    : db(db), last_request_time(std::chrono::steady_clock::now() - std::chrono::seconds(1)) {
}

void ENAClient::RespectRateLimit() {
	// Reserve the next request slot under the lock, then sleep outside it.
	// Sleeping under the lock would serialize every caller behind one thread's
	// sleep and risks UB on recursive acquisition (std::mutex is non-recursive).
	const auto min_interval = std::chrono::milliseconds(static_cast<int>(1000.0 / RATE_LIMIT));
	std::chrono::steady_clock::time_point slot;
	{
		std::lock_guard<std::mutex> lock(rate_limit_mutex);
		auto earliest_next = last_request_time + min_interval;
		auto now = std::chrono::steady_clock::now();
		slot = (earliest_next > now) ? earliest_next : now;
		last_request_time = slot;
	}
	if (slot > std::chrono::steady_clock::now()) {
		std::this_thread::sleep_until(slot);
	}
}

bool ENAClient::IsRetryableStatus(int status) {
	return status == 429 || status == 500 || status == 502 || status == 503;
}

std::string ENAClient::MakeRequest(const std::string &url) {
	RespectRateLimit();

	duckdb::HTTPHeaders headers(db);
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

		if (attempt < MAX_RETRIES && !response->HasRequestError() && IsRetryableStatus(int(response->status))) {
			std::this_thread::sleep_for(std::chrono::milliseconds(retry_delay_ms));
			retry_delay_ms *= 2;
			continue;
		}

		if (response->HasRequestError()) {
			throw duckdb::IOException("ENA request failed: %s (URL: %s)", response->GetRequestError(), url);
		}
		throw duckdb::HTTPException(*response, "ENA request failed with HTTP %d (URL: %s)", int(response->status), url);
	}

	throw duckdb::IOException("ENA request failed after %d retries (URL: %s)", MAX_RETRIES, url);
}

std::string ENAClient::Search(const std::string &accession, const std::string &result_type, const std::string &fields) {
	auto url = ENAParser::BuildSearchURL(accession, result_type, fields);
	return MakeRequest(url);
}

std::string ENAClient::FetchXML(const std::vector<std::string> &accessions) {
	auto url = ENAParser::BuildXMLURL(accessions);
	return MakeRequest(url);
}

std::string ENAClient::FetchURL(const std::string &url) {
	return MakeRequest(url);
}

} // namespace miint
