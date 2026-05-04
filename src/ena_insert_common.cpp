// SPDX-License-Identifier: MIT
//
// Implementation of shared helpers used by the ENA INSERT operators.
// See ena_insert_common.hpp for the contract.

#include "ena_insert_common.hpp"

#include "ena_storage.hpp"

#include "duckdb/catalog/catalog_transaction.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/date.hpp"
#include "duckdb/common/types/timestamp.hpp"
#include "duckdb/common/types/uuid.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/secret/secret.hpp"
#include "duckdb/main/secret/secret_manager.hpp"

namespace duckdb {

namespace {

unique_ptr<SecretEntry> LookupENASecret(ClientContext &context, const string &secret_name) {
	auto &secret_manager = SecretManager::Get(context);
	auto transaction = CatalogTransaction::GetSystemCatalogTransaction(context);
	// Empty storage string asks the manager to search every registered backend
	// (memory, local_file, plus any extension-provided backends). Narrowing to
	// memory+local_file would miss secrets stored elsewhere — confusing for
	// users whose `duckdb_secrets()` shows the secret but `INSERT` complains
	// it's missing.
	return secret_manager.GetSecretByName(transaction, secret_name);
}

} // namespace

ResolvedENACredentials ResolveENACredentials(ClientContext &context, ENACatalog &catalog) {
	ResolvedENACredentials creds;
	creds.endpoint = catalog.GetEndpoint();
	creds.endpoint_url = catalog.GetEndpointURL();
	creds.secret_name = catalog.GetSecretName();
	if (creds.secret_name.empty()) {
		throw BinderException("ENA INSERT requires a SECRET — re-attach with (TYPE ENA, SECRET 'name')");
	}
	auto entry = LookupENASecret(context, creds.secret_name);
	if (!entry) {
		throw BinderException("ENA secret '%s' not found", creds.secret_name);
	}
	auto kv = dynamic_cast<const KeyValueSecret *>(entry->secret.get());
	if (!kv) {
		throw InvalidInputException("ENA secret '%s' is not a KeyValueSecret — was it created with TYPE ENA?",
		                            creds.secret_name);
	}
	auto user_val = kv->TryGetValue("user");
	auto password_val = kv->TryGetValue("password");
	if (user_val.IsNull() || password_val.IsNull()) {
		throw BinderException("ENA secret '%s' is missing user or password", creds.secret_name);
	}
	creds.user = user_val.ToString();
	creds.password = password_val.ToString();
	auto endpoint_url_val = kv->TryGetValue("endpoint_url");
	if (!endpoint_url_val.IsNull()) {
		creds.endpoint_url = endpoint_url_val.ToString();
	}
	return creds;
}

idx_t ResolveInputColumn(const physical_index_vector_t<idx_t> &column_index_map, idx_t input_chunk_columns,
                         idx_t table_column_index) {
	if (column_index_map.empty()) {
		return table_column_index < input_chunk_columns ? table_column_index : DConstants::INVALID_INDEX;
	}
	// `physical_index_vector_t` is sparse; out-of-range reads silently return
	// `INVALID_INDEX` (the "not provided" sentinel). Treat an undersized map
	// the same way to keep behaviour deterministic if a future refactor of
	// LogicalInsert ever truncates the map.
	if (table_column_index >= column_index_map.size()) {
		return DConstants::INVALID_INDEX;
	}
	const auto mapped = column_index_map[PhysicalIndex(table_column_index)];
	if (mapped == DConstants::INVALID_INDEX) {
		return DConstants::INVALID_INDEX;
	}
	return mapped;
}

string ValueToVarchar(const Value &v) {
	if (v.IsNull()) {
		return string();
	}
	return v.ToString();
}

string ValueToDateString(const Value &v) {
	if (v.IsNull()) {
		return string();
	}
	if (v.type().id() != LogicalTypeId::DATE) {
		throw InvalidInputException("ENA insert: hold_until_date must be a DATE");
	}
	return Date::ToString(v.GetValue<date_t>());
}

string GenerateSubmissionId() {
	return UUID::ToString(UUID::GenerateRandomUUID());
}

void RecordSubmissionLog(ENACatalog &catalog, const ResolvedENACredentials &creds,
                         const SubmissionLogPayload &payload) {
	ENASubmissionLogRow row;
	row.submission_id = GenerateSubmissionId();
	row.submitted_at = timestamp_tz_t(Timestamp::GetCurrentTimestamp());
	row.endpoint = creds.endpoint;
	row.secret_name = creds.secret_name;
	row.action = payload.action;
	row.object_type = payload.object_type;
	row.n_objects = payload.n_objects;
	row.success = payload.success;
	row.era_accession = payload.era_accession;
	row.request_payload = payload.envelope_payload;
	row.receipt = payload.raw_receipt;
	row.error_messages.clear();
	row.error_messages.reserve(payload.error_messages.size());
	for (const auto &m : payload.error_messages) {
		row.error_messages.push_back(m);
	}
	row.duration_ms = payload.duration_ms;
	catalog.GetSubmissionLog().Append(row);
}

} // namespace duckdb
