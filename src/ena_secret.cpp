// SPDX-License-Identifier: MIT
//
// ENA Webin secret type. See src/include/ena_secret.hpp.

#include "ena_secret.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/main/secret/secret.hpp"

namespace miint {

namespace {

// Read a single named option, lowercased, from CreateSecretInput.
// Returns empty string when the option is not provided.
std::string GetOption(const duckdb::CreateSecretInput &input, const std::string &key) {
	auto it = input.options.find(key);
	if (it == input.options.end()) {
		return "";
	}
	if (it->second.IsNull()) {
		return "";
	}
	return it->second.ToString();
}

duckdb::unique_ptr<duckdb::BaseSecret> CreateENASecret(duckdb::ClientContext &, duckdb::CreateSecretInput &input) {
	const auto user = GetOption(input, "user");
	if (user.empty()) {
		throw duckdb::InvalidInputException("ENA secret requires 'user' (your Webin-XXXXX account id)");
	}

	const auto password_literal = GetOption(input, "password");
	const auto password_env = GetOption(input, "password_env");
	const auto password_file = GetOption(input, "password_file");
	std::string password_resolved;
	try {
		password_resolved = ResolvePasswordIndirection(password_literal, password_env, password_file);
	} catch (const std::exception &e) {
		throw duckdb::InvalidInputException(e.what());
	}

	auto endpoint = GetOption(input, "endpoint");
	if (endpoint.empty()) {
		endpoint = "test";
	} else if (endpoint != "test" && endpoint != "production") {
		throw duckdb::InvalidInputException("ENA secret: endpoint must be 'test' or 'production' (got '%s')", endpoint);
	}

	// Optional override of the resolved base URL — primarily for the local
	// mock server in tests, but also useful if EBI moves the V2 endpoint.
	const auto endpoint_url = GetOption(input, "endpoint_url");

	auto secret = duckdb::make_uniq<duckdb::KeyValueSecret>(input.scope, input.type, input.provider, input.name);
	secret->secret_map["user"] = duckdb::Value(user);
	secret->secret_map["password"] = duckdb::Value(password_resolved);
	secret->secret_map["endpoint"] = duckdb::Value(endpoint);
	if (!endpoint_url.empty()) {
		secret->secret_map["endpoint_url"] = duckdb::Value(endpoint_url);
	}
	secret->redact_keys = {"password"};
	return secret;
}

} // namespace

void RegisterENASecretType(duckdb::ExtensionLoader &loader) {
	duckdb::SecretType secret_type;
	secret_type.name = "ena";
	secret_type.deserializer = duckdb::KeyValueSecret::Deserialize<duckdb::KeyValueSecret>;
	secret_type.default_provider = "config";
	loader.RegisterSecretType(secret_type);

	duckdb::CreateSecretFunction fn;
	fn.secret_type = "ena";
	fn.provider = "config";
	fn.function = CreateENASecret;
	fn.named_parameters["user"] = duckdb::LogicalType::VARCHAR;
	fn.named_parameters["password"] = duckdb::LogicalType::VARCHAR;
	fn.named_parameters["password_env"] = duckdb::LogicalType::VARCHAR;
	fn.named_parameters["password_file"] = duckdb::LogicalType::VARCHAR;
	fn.named_parameters["endpoint"] = duckdb::LogicalType::VARCHAR;
	fn.named_parameters["endpoint_url"] = duckdb::LogicalType::VARCHAR;
	// Bearer-token auth is not supported by Webin V2 (HTTP Basic only); see
	// localdocs/ena-research-webin-v2-deep.md §3. Add `bearer_token` here
	// only if/when that changes.
	loader.RegisterFunction(fn);
}

} // namespace miint
