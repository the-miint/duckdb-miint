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

// Shared GET retry loop. Caller supplies a pre-built HTTPHeaders (empty for
// the unauthenticated path, with an Authorization: Basic header for the
// authenticated path). Caller also handles rate-limiting — the public read
// path needs the global 3 req/s cap, but the authenticated submission path
// is bursty and runs against the user's own account, so it bypasses.
std::string ENAClient::ExecuteGet(const std::string &url, duckdb::HTTPHeaders &headers, const char *failure_label) {
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
			throw duckdb::IOException("%s failed: %s (URL: %s)", failure_label, response->GetRequestError(), url);
		}
		throw duckdb::HTTPException(*response, "%s failed with HTTP %d (URL: %s)", failure_label, int(response->status),
		                            url);
	}

	throw duckdb::IOException("%s failed after %d retries (URL: %s)", failure_label, MAX_RETRIES, url);
}

std::string ENAClient::MakeRequest(const std::string &url) {
	RespectRateLimit();
	duckdb::HTTPHeaders headers(db);
	return ExecuteGet(url, headers, "ENA request");
}

// PostBody intentionally does NOT call RespectRateLimit(): the submission
// endpoints are authenticated and run against the user's own Webin account,
// where the public 3 req/s cap is overly conservative.
std::string ENAClient::PostBody(const std::string &url, const std::string &body, const std::string &content_type,
                                const std::string &accept_type, const std::string &user, const std::string &password) {
	duckdb::HTTPHeaders headers(db);
	headers.Insert("Authorization", BuildBasicAuthHeader(user, password));
	headers.Insert("Content-Type", content_type);
	headers.Insert("Accept", accept_type);

	auto &http_util = duckdb::HTTPUtil::Get(db);
	auto params = http_util.InitializeParameters(db, url);

	const auto *buf_ptr = reinterpret_cast<duckdb::const_data_ptr_t>(body.data());
	int retry_delay_ms = INITIAL_RETRY_DELAY_MS;

	// PostRequestInfo carries response state in `buffer_out` (per DuckDB issue
	// duckdb/duckdb#19062), so we rebuild it each iteration rather than
	// reusing a single instance — a retried POST should start with a clean
	// output buffer. Loop body always exits via return or throw; the post-
	// loop throw below is unreachable but kept to make the contract explicit.
	for (int attempt = 0; attempt <= MAX_RETRIES; attempt++) {
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
			throw duckdb::IOException("ENA POST failed: %s (URL: %s)", response->GetRequestError(), url);
		}
		throw duckdb::HTTPException(*response, "ENA POST failed with HTTP %d (URL: %s)", int(response->status), url);
	}

	throw duckdb::IOException("ENA POST failed after %d retries (URL: %s)", MAX_RETRIES, url); // unreachable
}

std::string ENAClient::PostJSON(const std::string &url, const std::string &body, const std::string &user,
                                const std::string &password) {
	return PostBody(url, body, "application/json", "application/json", user, password);
}

std::string ENAClient::PostXML(const std::string &url, const std::string &body, const std::string &user,
                               const std::string &password) {
	return PostBody(url, body, "application/xml", "application/xml", user, password);
}

std::string ENAClient::PostJSONReceiveXML(const std::string &url, const std::string &body, const std::string &user,
                                          const std::string &password) {
	return PostBody(url, body, "application/json", "application/xml", user, password);
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

std::string ENAClient::AuthenticatedGet(const std::string &url, const std::string &user, const std::string &password) {
	duckdb::HTTPHeaders headers(db);
	headers.Insert("Authorization", BuildBasicAuthHeader(user, password));
	return ExecuteGet(url, headers, "ENA authenticated GET");
}

} // namespace miint
