// SPDX-License-Identifier: MIT
//
// Storage extension that exposes the ENA Webin V2 object graph as a virtual
// catalog. ATTACH '' AS ena (TYPE ENA[, SECRET 'name']) registers a single
// schema "main" with five hard-coded submission tables (projects, samples,
// experiments, runs, analyses) plus the bookkeeping table submission_log.
//
// INSERT INTO ena.projects (Phase 4) builds an envelope, POSTs it through
// ENAClient, parses the receipt, returns the server-assigned columns via
// RETURNING, and appends a row to submission_log.
//
// SELECT against the submission tables is not supported in Phase 4 (deferred
// to Phase 5 via the Reports API). SELECT against submission_log returns the
// in-memory append log.

#pragma once

#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/schema_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/common/mutex.hpp"
#include "duckdb/storage/storage_extension.hpp"
#include "duckdb/transaction/transaction.hpp"
#include "duckdb/transaction/transaction_manager.hpp"

#include <memory>
#include <string>

namespace duckdb {

class ENACatalog;

// Identifies one of the predefined ENA tables. Used to dispatch INSERT and
// scan logic from the otherwise-uniform ENATableEntry.
enum class ENATableKind : uint8_t {
	PROJECTS,
	SAMPLES,
	EXPERIMENTS,
	RUNS,
	ANALYSES,
	SUBMISSION_LOG,
};

struct ENATableSchema {
	ENATableKind kind;
	const char *table_name;
};

const vector<ENATableSchema> &ENATables();

// Build a CreateTableInfo with the predefined columns for a given ENA table.
unique_ptr<CreateTableInfo> BuildENATableInfo(SchemaCatalogEntry &schema, ENATableKind kind);

class ENATableEntry : public TableCatalogEntry {
public:
	ENATableEntry(Catalog &catalog, SchemaCatalogEntry &schema, CreateTableInfo &info, ENATableKind kind);

public:
	unique_ptr<BaseStatistics> GetStatistics(ClientContext &context, column_t column_id) override;
	TableFunction GetScanFunction(ClientContext &context, unique_ptr<FunctionData> &bind_data) override;
	TableStorageInfo GetStorageInfo(ClientContext &context) override;

	ENATableKind GetKind() const {
		return kind;
	}

private:
	ENATableKind kind;
};

class ENASchemaEntry : public SchemaCatalogEntry {
public:
	ENASchemaEntry(Catalog &catalog, CreateSchemaInfo &info);

public:
	optional_ptr<CatalogEntry> CreateTable(CatalogTransaction transaction, BoundCreateTableInfo &info) override;
	optional_ptr<CatalogEntry> CreateFunction(CatalogTransaction transaction, CreateFunctionInfo &info) override;
	optional_ptr<CatalogEntry> CreateIndex(CatalogTransaction transaction, CreateIndexInfo &info,
	                                       TableCatalogEntry &table) override;
	optional_ptr<CatalogEntry> CreateView(CatalogTransaction transaction, CreateViewInfo &info) override;
	optional_ptr<CatalogEntry> CreateSequence(CatalogTransaction transaction, CreateSequenceInfo &info) override;
	optional_ptr<CatalogEntry> CreateTableFunction(CatalogTransaction transaction,
	                                               CreateTableFunctionInfo &info) override;
	optional_ptr<CatalogEntry> CreateCopyFunction(CatalogTransaction transaction,
	                                              CreateCopyFunctionInfo &info) override;
	optional_ptr<CatalogEntry> CreatePragmaFunction(CatalogTransaction transaction,
	                                                CreatePragmaFunctionInfo &info) override;
	optional_ptr<CatalogEntry> CreateCollation(CatalogTransaction transaction, CreateCollationInfo &info) override;
	optional_ptr<CatalogEntry> CreateType(CatalogTransaction transaction, CreateTypeInfo &info) override;

	void Alter(CatalogTransaction transaction, AlterInfo &info) override;
	void Scan(ClientContext &context, CatalogType type, const std::function<void(CatalogEntry &)> &callback) override;
	void Scan(CatalogType type, const std::function<void(CatalogEntry &)> &callback) override;
	void DropEntry(ClientContext &context, DropInfo &info) override;
	optional_ptr<CatalogEntry> LookupEntry(CatalogTransaction transaction, const EntryLookupInfo &lookup_info) override;

private:
	// One ENATableEntry per predefined ENA table, indexed by ENATableKind.
	vector<unique_ptr<ENATableEntry>> tables;
};

class ENATransaction : public Transaction {
public:
	ENATransaction(TransactionManager &manager, ClientContext &context);
};

class ENATransactionManager : public TransactionManager {
public:
	explicit ENATransactionManager(AttachedDatabase &db);
	Transaction &StartTransaction(ClientContext &context) override;
	ErrorData CommitTransaction(ClientContext &context, Transaction &transaction) override;
	void RollbackTransaction(Transaction &transaction) override;
	void Checkpoint(ClientContext &context, bool force = false) override;

private:
	mutex transaction_lock;
	reference_map_t<Transaction, unique_ptr<ENATransaction>> transactions;
};

// Append-only in-memory backing for ena.submission_log. Rows are written by
// the INSERT physical operators and read back by the table scan function.
struct ENASubmissionLogRow {
	string submission_id;
	timestamp_tz_t submitted_at;
	string endpoint;
	string secret_name;
	string action;
	string object_type;
	int32_t n_objects;
	bool success;
	string era_accession;
	string request_payload;
	string receipt;
	vector<string> error_messages;
	int64_t duration_ms;
};

class ENASubmissionLog {
public:
	ENASubmissionLog();
	void Append(const ENASubmissionLogRow &row);
	// Snapshot copy; safe to scan without holding the lock.
	vector<ENASubmissionLogRow> Snapshot();

private:
	mutex log_lock;
	vector<ENASubmissionLogRow> rows;
};

class ENACatalog : public Catalog {
public:
	ENACatalog(AttachedDatabase &db, string secret_name, string endpoint, string endpoint_url);
	~ENACatalog() override;

	const string &GetSecretName() const {
		return secret_name;
	}
	const string &GetEndpoint() const {
		return endpoint;
	}
	const string &GetEndpointURL() const {
		return endpoint_url;
	}
	ENASubmissionLog &GetSubmissionLog() {
		return *submission_log;
	}

	void Initialize(bool load_builtin) override;
	string GetCatalogType() override {
		return "ena";
	}
	string GetDefaultSchema() const override {
		return "main";
	}

	optional_ptr<CatalogEntry> CreateSchema(CatalogTransaction transaction, CreateSchemaInfo &info) override;
	optional_ptr<SchemaCatalogEntry> LookupSchema(CatalogTransaction transaction, const EntryLookupInfo &schema_lookup,
	                                              OnEntryNotFound if_not_found) override;
	void ScanSchemas(ClientContext &context, std::function<void(SchemaCatalogEntry &)> callback) override;

	PhysicalOperator &PlanCreateTableAs(ClientContext &context, PhysicalPlanGenerator &planner, LogicalCreateTable &op,
	                                    PhysicalOperator &plan) override;
	PhysicalOperator &PlanInsert(ClientContext &context, PhysicalPlanGenerator &planner, LogicalInsert &op,
	                             optional_ptr<PhysicalOperator> plan) override;
	PhysicalOperator &PlanDelete(ClientContext &context, PhysicalPlanGenerator &planner, LogicalDelete &op,
	                             PhysicalOperator &plan) override;
	PhysicalOperator &PlanUpdate(ClientContext &context, PhysicalPlanGenerator &planner, LogicalUpdate &op,
	                             PhysicalOperator &plan) override;

	DatabaseSize GetDatabaseSize(ClientContext &context) override;
	bool InMemory() override {
		return true;
	}
	string GetDBPath() override {
		return endpoint;
	}

	CatalogLookupBehavior CatalogTypeLookupRule(CatalogType type) const override {
		switch (type) {
		case CatalogType::TABLE_ENTRY:
		case CatalogType::SCHEMA_ENTRY:
			return CatalogLookupBehavior::STANDARD;
		default:
			return CatalogLookupBehavior::NEVER_LOOKUP;
		}
	}

private:
	void DropSchema(ClientContext &context, DropInfo &info) override;

private:
	string secret_name;
	string endpoint;     // "test" or "production"
	string endpoint_url; // resolved base URL (override or derived)
	unique_ptr<ENASchemaEntry> main_schema;
	unique_ptr<ENASubmissionLog> submission_log;
};

class ENAStorageExtension : public StorageExtension {
public:
	ENAStorageExtension();
};

} // namespace duckdb
