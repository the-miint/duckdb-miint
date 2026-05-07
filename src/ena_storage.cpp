// SPDX-License-Identifier: MIT
//
// ENA Webin V2 storage extension. See src/include/ena_storage.hpp.

#include "ena_storage.hpp"

#include "ena_experiments_insert_op.hpp"
#include "ena_projects_insert_op.hpp"
#include "ena_runs_insert_op.hpp"
#include "ena_samples_insert_op.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/timestamp.hpp"
#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/attached_database.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/parser/parsed_data/attach_info.hpp"
#include "duckdb/parser/parsed_data/create_schema_info.hpp"
#include "duckdb/parser/parsed_data/create_table_info.hpp"
#include "duckdb/parser/parsed_data/drop_info.hpp"
#include "duckdb/planner/operator/logical_insert.hpp"
#include "duckdb/storage/database_size.hpp"
#include "duckdb/storage/table_storage_info.hpp"

namespace duckdb {

//===--------------------------------------------------------------------===//
// Predefined table schemas
//===--------------------------------------------------------------------===//

const vector<ENATableSchema> &ENATables() {
	static const vector<ENATableSchema> tables = {
	    {ENATableKind::PROJECTS, "projects"},       {ENATableKind::SAMPLES, "samples"},
	    {ENATableKind::EXPERIMENTS, "experiments"}, {ENATableKind::RUNS, "runs"},
	    {ENATableKind::ANALYSES, "analyses"},       {ENATableKind::SUBMISSION_LOG, "submission_log"},
	};
	return tables;
}

// Schema-free column emitter for `submission_log`. Defined as a free function
// so it can be called both by `BuildENATableInfo` (which builds a
// `CreateTableInfo` for catalog registration) and by the standalone scan
// function (which only needs a name+type list). Single source of truth.
void AddSubmissionLogColumns(ColumnList &columns) {
	auto add = [&](const char *name, LogicalType type) {
		columns.AddColumn(ColumnDefinition(name, std::move(type)));
	};
	// VARCHAR rather than UUID: the row already carries the string form, and
	// the scan path uses StringVector::AddString. Going through UUID would
	// require converting through `hugeint_t` on every emit for no user-visible
	// gain (the log isn't a primary-key/join target).
	add("submission_id", LogicalType::VARCHAR);
	add("submitted_at", LogicalType::TIMESTAMP_TZ);
	add("endpoint", LogicalType::VARCHAR);
	add("secret_name", LogicalType::VARCHAR);
	add("action", LogicalType::VARCHAR);
	add("object_type", LogicalType::VARCHAR);
	add("n_objects", LogicalType::INTEGER);
	add("success", LogicalType::BOOLEAN);
	add("era_accession", LogicalType::VARCHAR);
	add("request_payload", LogicalType::VARCHAR);
	add("receipt", LogicalType::VARCHAR);
	add("error_messages", LogicalType::LIST(LogicalType::VARCHAR));
	// Duration is int64 ms — slow ENA submissions over a saturated wwwdev or a
	// large multi-object batch can plausibly exceed the int32 ~35-minute cap.
	add("duration_ms", LogicalType::BIGINT);
	// Lifecycle target accession or refname. Empty for ADD / MODIFY / VALIDATE
	// (those identify their objects via the body); populated for CANCEL /
	// RELEASE / HOLD with the value sent on `target=`.
	add("target", LogicalType::VARCHAR);
}

unique_ptr<CreateTableInfo> BuildENATableInfo(SchemaCatalogEntry &schema, ENATableKind kind) {
	auto info = make_uniq<CreateTableInfo>(schema, "");

	auto add = [&](const char *name, LogicalType type) {
		info->columns.AddColumn(ColumnDefinition(name, std::move(type)));
	};

	switch (kind) {
	case ENATableKind::PROJECTS:
		info->table = "projects";
		add("alias", LogicalType::VARCHAR);
		add("title", LogicalType::VARCHAR);
		add("description", LogicalType::VARCHAR);
		add("project_type", LogicalType::VARCHAR);
		add("umbrella_children", LogicalType::LIST(LogicalType::VARCHAR));
		add("attributes", LogicalType::MAP(LogicalType::VARCHAR, LogicalType::VARCHAR));
		add("prjeb_accession", LogicalType::VARCHAR);
		add("erp_accession", LogicalType::VARCHAR);
		add("hold_until_date", LogicalType::DATE);
		break;
	case ENATableKind::SAMPLES:
		info->table = "samples";
		// Column order matters — BuildFromBuffer uses positional COL_*
		// constants in src/ena_samples_insert_op.cpp.
		add("alias", LogicalType::VARCHAR);
		add("title", LogicalType::VARCHAR);
		add("description", LogicalType::VARCHAR);
		add("taxon_id", LogicalType::INTEGER);
		add("scientific_name", LogicalType::VARCHAR);
		add("checklist", LogicalType::VARCHAR);
		add("attributes", LogicalType::MAP(LogicalType::VARCHAR, LogicalType::VARCHAR));
		// Sparse per-attribute units. Some checklist attributes (e.g.
		// ERC000015 lat/lon → 'DD') are rejected by the server without
		// units. Keys here are matched by tag against `attributes`; tags
		// not present in this map emit no `units` JSON field.
		add("attribute_units", LogicalType::MAP(LogicalType::VARCHAR, LogicalType::VARCHAR));
		add("ers_accession", LogicalType::VARCHAR);
		add("samea_accession", LogicalType::VARCHAR);
		break;
	case ENATableKind::EXPERIMENTS:
		info->table = "experiments";
		// Column order matters — BuildFromBuffer uses positional COL_*
		// constants in src/ena_experiments_insert_op.cpp.
		add("alias", LogicalType::VARCHAR);
		add("title", LogicalType::VARCHAR);
		add("study_ref", LogicalType::VARCHAR);
		add("sample_descriptor", LogicalType::VARCHAR);
		add("design_description", LogicalType::VARCHAR);
		add("library_name", LogicalType::VARCHAR);
		add("library_strategy", LogicalType::VARCHAR);
		add("library_source", LogicalType::VARCHAR);
		add("library_selection", LogicalType::VARCHAR);
		add("library_layout", LogicalType::VARCHAR); // "SINGLE" or "PAIRED"
		add("platform", LogicalType::VARCHAR);
		add("instrument_model", LogicalType::VARCHAR);
		add("erx_accession", LogicalType::VARCHAR);
		break;
	case ENATableKind::RUNS:
		info->table = "runs";
		// Column order matters — BuildFromBuffer uses positional COL_*
		// constants in src/ena_runs_insert_op.cpp. The `files` LIST(STRUCT)
		// shape matches the SRA.run.xsd <FILES> element: each entry carries
		// filename + filetype + md5. Server re-computes MD5 after upload and
		// compares against this value.
		add("alias", LogicalType::VARCHAR);
		add("experiment_ref", LogicalType::VARCHAR);
		add("title", LogicalType::VARCHAR);
		add("files", LogicalType::LIST(LogicalType::STRUCT({{"filename", LogicalType::VARCHAR},
		                                                    {"filetype", LogicalType::VARCHAR},
		                                                    {"md5", LogicalType::VARCHAR}})));
		add("err_accession", LogicalType::VARCHAR);
		break;
	case ENATableKind::ANALYSES:
		info->table = "analyses";
		add("alias", LogicalType::VARCHAR);
		add("study_ref", LogicalType::VARCHAR);
		add("analysis_type", LogicalType::VARCHAR);
		add("accession", LogicalType::VARCHAR);
		break;
	case ENATableKind::SUBMISSION_LOG:
		info->table = "submission_log";
		AddSubmissionLogColumns(info->columns);
		break;
	}
	return info;
}

//===--------------------------------------------------------------------===//
// Submission log scan function
//===--------------------------------------------------------------------===//
namespace {

struct ENASubmissionLogBindData : public TableFunctionData {
	explicit ENASubmissionLogBindData(ENASubmissionLog &log) : log(log) {
	}
	ENASubmissionLog &log;
};

struct ENASubmissionLogState : public GlobalTableFunctionState {
	vector<ENASubmissionLogRow> snapshot;
	idx_t cursor = 0;
};

void EmitSubmissionLogColumns(vector<LogicalType> &return_types, vector<string> &names) {
	// Single source of truth for the submission_log column list — defer to the
	// declarations in `AddSubmissionLogColumns` so the catalog and the scan
	// function cannot drift apart.
	ColumnList columns;
	AddSubmissionLogColumns(columns);
	for (const auto &col : columns.Logical()) {
		names.push_back(col.GetName());
		return_types.push_back(col.GetType());
	}
}

unique_ptr<FunctionData> ENASubmissionLogBind(ClientContext &, TableFunctionBindInput &,
                                              vector<LogicalType> &return_types, vector<string> &names) {
	EmitSubmissionLogColumns(return_types, names);
	// Bind data must be supplied via GetScanFunction's bind_data out-parameter,
	// not here; this bind callback exists only to satisfy DuckDB's table-scan
	// pipeline (it needs *a* bind, even when the "real" data was set up earlier).
	return nullptr;
}

unique_ptr<GlobalTableFunctionState> ENASubmissionLogInit(ClientContext &, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<ENASubmissionLogBindData>();
	auto state = make_uniq<ENASubmissionLogState>();
	state->snapshot = bind_data.log.Snapshot();
	return std::move(state);
}

void ENASubmissionLogScan(ClientContext &, TableFunctionInput &data, DataChunk &output) {
	auto &state = data.global_state->Cast<ENASubmissionLogState>();
	const idx_t available = state.snapshot.size() - state.cursor;
	const idx_t produce = MinValue<idx_t>(STANDARD_VECTOR_SIZE, available);
	if (produce == 0) {
		output.SetCardinality(0);
		return;
	}
	output.SetCardinality(produce);

	auto submission_id = FlatVector::GetData<string_t>(output.data[0]);
	auto submitted_at = FlatVector::GetData<timestamp_tz_t>(output.data[1]);
	auto endpoint = FlatVector::GetData<string_t>(output.data[2]);
	auto secret_name = FlatVector::GetData<string_t>(output.data[3]);
	auto action = FlatVector::GetData<string_t>(output.data[4]);
	auto object_type = FlatVector::GetData<string_t>(output.data[5]);
	auto n_objects = FlatVector::GetData<int32_t>(output.data[6]);
	auto success = FlatVector::GetData<bool>(output.data[7]);
	auto era_accession = FlatVector::GetData<string_t>(output.data[8]);
	auto request_payload = FlatVector::GetData<string_t>(output.data[9]);
	auto receipt = FlatVector::GetData<string_t>(output.data[10]);
	auto duration_ms = FlatVector::GetData<int64_t>(output.data[12]);
	auto target = FlatVector::GetData<string_t>(output.data[13]);

	auto &error_messages = output.data[11];
	ListVector::SetListSize(error_messages, 0);

	for (idx_t i = 0; i < produce; i++) {
		const auto &row = state.snapshot[state.cursor + i];
		submission_id[i] = StringVector::AddString(output.data[0], row.submission_id);
		submitted_at[i] = row.submitted_at;
		endpoint[i] = StringVector::AddString(output.data[2], row.endpoint);
		secret_name[i] = StringVector::AddString(output.data[3], row.secret_name);
		action[i] = StringVector::AddString(output.data[4], row.action);
		object_type[i] = StringVector::AddString(output.data[5], row.object_type);
		n_objects[i] = row.n_objects;
		success[i] = row.success;
		era_accession[i] = StringVector::AddString(output.data[8], row.era_accession);
		request_payload[i] = StringVector::AddString(output.data[9], row.request_payload);
		receipt[i] = StringVector::AddString(output.data[10], row.receipt);
		duration_ms[i] = row.duration_ms;
		target[i] = StringVector::AddString(output.data[13], row.target);
	}

	auto child_offset = ListVector::GetListSize(error_messages);
	for (idx_t i = 0; i < produce; i++) {
		const auto &row = state.snapshot[state.cursor + i];
		const auto count = static_cast<idx_t>(row.error_messages.size());
		ListVector::Reserve(error_messages, child_offset + count);
		auto &child_vec = ListVector::GetEntry(error_messages);
		auto child_data = FlatVector::GetData<string_t>(child_vec);
		for (idx_t j = 0; j < count; j++) {
			child_data[child_offset + j] = StringVector::AddString(child_vec, row.error_messages[j]);
		}
		auto entries = ListVector::GetData(error_messages);
		entries[i].offset = child_offset;
		entries[i].length = count;
		child_offset += count;
	}
	ListVector::SetListSize(error_messages, child_offset);

	state.cursor += produce;
}

// Bind throws before scan is ever called, so a paired scan-function is dead
// code; route both pointers at the same throwing helper to keep them in sync.
unique_ptr<FunctionData> NotImplementedBind(ClientContext &, TableFunctionBindInput &input,
                                            vector<LogicalType> &return_types, vector<string> &names) {
	throw NotImplementedException("Reading from ena.%s is not supported in this build "
	                              "(SELECT support is planned for a future phase).",
	                              input.table_function.name);
}

void NotImplementedScanFn(ClientContext &, TableFunctionInput &, DataChunk &) {
	// Unreachable: NotImplementedBind throws before this is invoked.
	throw NotImplementedException("ENA virtual table scan invoked unexpectedly");
}

} // namespace

//===--------------------------------------------------------------------===//
// ENATableEntry
//===--------------------------------------------------------------------===//

ENATableEntry::ENATableEntry(Catalog &catalog, SchemaCatalogEntry &schema, CreateTableInfo &info, ENATableKind kind_p)
    : TableCatalogEntry(catalog, schema, info), kind(kind_p) {
}

unique_ptr<BaseStatistics> ENATableEntry::GetStatistics(ClientContext &, column_t) {
	return nullptr;
}

TableFunction ENATableEntry::GetScanFunction(ClientContext &, unique_ptr<FunctionData> &bind_data) {
	if (kind == ENATableKind::SUBMISSION_LOG) {
		auto &ena_catalog = catalog.Cast<ENACatalog>();
		bind_data = make_uniq<ENASubmissionLogBindData>(ena_catalog.GetSubmissionLog());
		TableFunction fn("ena_submission_log", {}, ENASubmissionLogScan, ENASubmissionLogBind, ENASubmissionLogInit);
		fn.projection_pushdown = false;
		return fn;
	}
	bind_data = nullptr;
	TableFunction fn("ena_virtual_scan", {}, NotImplementedScanFn);
	fn.bind = NotImplementedBind;
	return fn;
}

TableStorageInfo ENATableEntry::GetStorageInfo(ClientContext &) {
	TableStorageInfo info;
	info.cardinality = 0;
	return info;
}

//===--------------------------------------------------------------------===//
// ENASchemaEntry
//===--------------------------------------------------------------------===//

ENASchemaEntry::ENASchemaEntry(Catalog &catalog, CreateSchemaInfo &info) : SchemaCatalogEntry(catalog, info) {
	for (auto &t : ENATables()) {
		auto table_info = BuildENATableInfo(*this, t.kind);
		tables.push_back(make_uniq<ENATableEntry>(catalog, *this, *table_info, t.kind));
	}
}

[[noreturn]] static void ThrowFixedSchema(const string &what) {
	throw BinderException("ENA catalog has a fixed schema; %s is not supported", what);
}

optional_ptr<CatalogEntry> ENASchemaEntry::CreateTable(CatalogTransaction, BoundCreateTableInfo &) {
	ThrowFixedSchema("CREATE TABLE");
}
optional_ptr<CatalogEntry> ENASchemaEntry::CreateFunction(CatalogTransaction, CreateFunctionInfo &) {
	ThrowFixedSchema("CREATE FUNCTION");
}
optional_ptr<CatalogEntry> ENASchemaEntry::CreateIndex(CatalogTransaction, CreateIndexInfo &, TableCatalogEntry &) {
	ThrowFixedSchema("CREATE INDEX");
}
optional_ptr<CatalogEntry> ENASchemaEntry::CreateView(CatalogTransaction, CreateViewInfo &) {
	ThrowFixedSchema("CREATE VIEW");
}
optional_ptr<CatalogEntry> ENASchemaEntry::CreateSequence(CatalogTransaction, CreateSequenceInfo &) {
	ThrowFixedSchema("CREATE SEQUENCE");
}
optional_ptr<CatalogEntry> ENASchemaEntry::CreateTableFunction(CatalogTransaction, CreateTableFunctionInfo &) {
	ThrowFixedSchema("CREATE TABLE FUNCTION");
}
optional_ptr<CatalogEntry> ENASchemaEntry::CreateCopyFunction(CatalogTransaction, CreateCopyFunctionInfo &) {
	ThrowFixedSchema("CREATE COPY FUNCTION");
}
optional_ptr<CatalogEntry> ENASchemaEntry::CreatePragmaFunction(CatalogTransaction, CreatePragmaFunctionInfo &) {
	ThrowFixedSchema("CREATE PRAGMA FUNCTION");
}
optional_ptr<CatalogEntry> ENASchemaEntry::CreateCollation(CatalogTransaction, CreateCollationInfo &) {
	ThrowFixedSchema("CREATE COLLATION");
}
optional_ptr<CatalogEntry> ENASchemaEntry::CreateType(CatalogTransaction, CreateTypeInfo &) {
	ThrowFixedSchema("CREATE TYPE");
}

void ENASchemaEntry::Alter(CatalogTransaction, AlterInfo &) {
	ThrowFixedSchema("ALTER");
}

void ENASchemaEntry::Scan(ClientContext &, CatalogType type, const std::function<void(CatalogEntry &)> &callback) {
	if (type != CatalogType::TABLE_ENTRY) {
		return;
	}
	for (auto &t : tables) {
		callback(*t);
	}
}

void ENASchemaEntry::Scan(CatalogType type, const std::function<void(CatalogEntry &)> &callback) {
	if (type != CatalogType::TABLE_ENTRY) {
		return;
	}
	for (auto &t : tables) {
		callback(*t);
	}
}

void ENASchemaEntry::DropEntry(ClientContext &, DropInfo &) {
	ThrowFixedSchema("DROP");
}

optional_ptr<CatalogEntry> ENASchemaEntry::LookupEntry(CatalogTransaction, const EntryLookupInfo &lookup_info) {
	if (lookup_info.GetCatalogType() != CatalogType::TABLE_ENTRY) {
		return nullptr;
	}
	const auto &name = lookup_info.GetEntryName();
	for (auto &t : tables) {
		if (StringUtil::CIEquals(t->name, name)) {
			return t.get();
		}
	}
	return nullptr;
}

//===--------------------------------------------------------------------===//
// ENATransaction & ENATransactionManager
//===--------------------------------------------------------------------===//

ENATransaction::ENATransaction(TransactionManager &manager, ClientContext &context) : Transaction(manager, context) {
}

ENATransactionManager::ENATransactionManager(AttachedDatabase &db) : TransactionManager(db) {
}

Transaction &ENATransactionManager::StartTransaction(ClientContext &context) {
	auto transaction = make_uniq<ENATransaction>(*this, context);
	auto &result = *transaction;
	lock_guard<mutex> l(transaction_lock);
	transactions[result] = std::move(transaction);
	return result;
}

ErrorData ENATransactionManager::CommitTransaction(ClientContext &, Transaction &transaction) {
	lock_guard<mutex> l(transaction_lock);
	transactions.erase(transaction);
	return ErrorData();
}

void ENATransactionManager::RollbackTransaction(Transaction &transaction) {
	lock_guard<mutex> l(transaction_lock);
	transactions.erase(transaction);
}

void ENATransactionManager::Checkpoint(ClientContext &, bool) {
}

//===--------------------------------------------------------------------===//
// ENASubmissionLog
//===--------------------------------------------------------------------===//

ENASubmissionLog::ENASubmissionLog() = default;

void ENASubmissionLog::Append(const ENASubmissionLogRow &row) {
	lock_guard<mutex> l(log_lock);
	rows.push_back(row);
}

vector<ENASubmissionLogRow> ENASubmissionLog::Snapshot() {
	lock_guard<mutex> l(log_lock);
	return rows;
}

//===--------------------------------------------------------------------===//
// ENACatalog
//===--------------------------------------------------------------------===//

ENACatalog::ENACatalog(AttachedDatabase &db, string secret_name_p, string endpoint_p, string endpoint_url_p)
    : Catalog(db), secret_name(std::move(secret_name_p)), endpoint(std::move(endpoint_p)),
      endpoint_url(std::move(endpoint_url_p)), submission_log(make_uniq<ENASubmissionLog>()) {
}

ENACatalog::~ENACatalog() = default;

void ENACatalog::Initialize(bool) {
	CreateSchemaInfo info;
	info.schema = "main";
	info.internal = true;
	main_schema = make_uniq<ENASchemaEntry>(*this, info);
}

optional_ptr<CatalogEntry> ENACatalog::CreateSchema(CatalogTransaction, CreateSchemaInfo &) {
	throw BinderException("ENA catalog: cannot CREATE SCHEMA (only the predefined 'main' schema exists)");
}

void ENACatalog::DropSchema(ClientContext &, DropInfo &) {
	throw BinderException("ENA catalog: cannot DROP SCHEMA");
}

void ENACatalog::ScanSchemas(ClientContext &, std::function<void(SchemaCatalogEntry &)> callback) {
	if (main_schema) {
		callback(*main_schema);
	}
}

optional_ptr<SchemaCatalogEntry> ENACatalog::LookupSchema(CatalogTransaction, const EntryLookupInfo &schema_lookup,
                                                          OnEntryNotFound if_not_found) {
	const auto &schema_name = schema_lookup.GetEntryName();
	if (StringUtil::CIEquals(schema_name, "main")) {
		return main_schema.get();
	}
	if (if_not_found == OnEntryNotFound::THROW_EXCEPTION) {
		throw BinderException("ENA catalog: schema '%s' not found (only 'main' is defined)", schema_name);
	}
	return nullptr;
}

PhysicalOperator &ENACatalog::PlanCreateTableAs(ClientContext &, PhysicalPlanGenerator &, LogicalCreateTable &,
                                                PhysicalOperator &) {
	throw BinderException("ENA catalog: CREATE TABLE AS is not supported");
}

PhysicalOperator &ENACatalog::PlanInsert(ClientContext &context, PhysicalPlanGenerator &planner, LogicalInsert &op,
                                         optional_ptr<PhysicalOperator> plan) {
	if (op.on_conflict_info.action_type != OnConflictAction::THROW) {
		throw BinderException("ENA catalog: ON CONFLICT clause is not supported");
	}
	auto &table_entry = op.table.Cast<ENATableEntry>();
	if (!plan) {
		throw BinderException("ENA catalog: INSERT requires a child plan");
	}
	switch (table_entry.GetKind()) {
	case ENATableKind::PROJECTS: {
		auto &op_ref = planner.Make<ENAProjectsInsert>(op, table_entry, op.column_index_map, op.return_chunk);
		op_ref.children.push_back(*plan);
		return op_ref;
	}
	case ENATableKind::SAMPLES: {
		auto &op_ref = planner.Make<ENASamplesInsert>(op, table_entry, op.column_index_map, op.return_chunk);
		op_ref.children.push_back(*plan);
		return op_ref;
	}
	case ENATableKind::EXPERIMENTS: {
		auto &op_ref = planner.Make<ENAExperimentsInsert>(op, table_entry, op.column_index_map, op.return_chunk);
		op_ref.children.push_back(*plan);
		return op_ref;
	}
	case ENATableKind::RUNS: {
		auto &op_ref = planner.Make<ENARunsInsert>(op, table_entry, op.column_index_map, op.return_chunk);
		op_ref.children.push_back(*plan);
		return op_ref;
	}
	default:
		throw BinderException("ENA catalog: INSERT INTO ena.%s is not implemented in this build", table_entry.name);
	}
}

PhysicalOperator &ENACatalog::PlanDelete(ClientContext &, PhysicalPlanGenerator &, LogicalDelete &,
                                         PhysicalOperator &) {
	throw BinderException("ENA catalog: DELETE is not supported");
}

PhysicalOperator &ENACatalog::PlanUpdate(ClientContext &, PhysicalPlanGenerator &, LogicalUpdate &,
                                         PhysicalOperator &) {
	throw BinderException("ENA catalog: UPDATE is not supported");
}

DatabaseSize ENACatalog::GetDatabaseSize(ClientContext &) {
	DatabaseSize size;
	size.bytes = 0;
	return size;
}

//===--------------------------------------------------------------------===//
// ENAStorageExtension
//===--------------------------------------------------------------------===//
string ResolveDefaultENAEndpointURL(const string &endpoint) {
	if (endpoint == "production") {
		return "https://www.ebi.ac.uk/ena/submit/webin-v2";
	}
	return "https://wwwdev.ebi.ac.uk/ena/submit/webin-v2";
}

namespace {

unique_ptr<Catalog> ENAAttach(optional_ptr<StorageExtensionInfo>, ClientContext &, AttachedDatabase &db, const string &,
                              AttachInfo &info, AttachOptions &options) {
	string secret_name;
	string endpoint;
	string endpoint_url_override;
	for (auto &entry : options.options) {
		auto lower_name = StringUtil::Lower(entry.first);
		if (lower_name == "secret") {
			secret_name = entry.second.ToString();
		} else if (lower_name == "endpoint") {
			endpoint = StringUtil::Lower(entry.second.ToString());
		} else if (lower_name == "endpoint_url") {
			endpoint_url_override = entry.second.ToString();
		} else {
			throw BinderException("Unrecognized option for ENA attach: %s", entry.first);
		}
	}
	// Path can be empty, "test", "production", or an http(s):// URL. Path
	// takes precedence over the ENDPOINT option only when it is non-empty.
	const string &path = info.path;
	if (!path.empty()) {
		const string lowered = StringUtil::Lower(path);
		if (lowered == "test" || lowered == "production") {
			endpoint = lowered;
		} else if (StringUtil::StartsWith(lowered, "http://") || StringUtil::StartsWith(lowered, "https://")) {
			endpoint_url_override = path;
		} else {
			throw BinderException("ENA attach: path must be empty, 'test', 'production', "
			                      "or an http(s):// URL (got '%s')",
			                      path);
		}
	}
	if (endpoint.empty()) {
		endpoint = "test";
	}
	if (endpoint != "test" && endpoint != "production") {
		throw BinderException("ENA attach: endpoint must be 'test' or 'production' (got '%s')", endpoint);
	}
	const string endpoint_url =
	    endpoint_url_override.empty() ? ResolveDefaultENAEndpointURL(endpoint) : endpoint_url_override;

	return make_uniq<ENACatalog>(db, std::move(secret_name), std::move(endpoint), endpoint_url);
}

unique_ptr<TransactionManager> ENACreateTransactionManager(optional_ptr<StorageExtensionInfo>, AttachedDatabase &db,
                                                           Catalog &) {
	return make_uniq<ENATransactionManager>(db);
}

} // namespace

ENAStorageExtension::ENAStorageExtension() {
	attach = ENAAttach;
	create_transaction_manager = ENACreateTransactionManager;
}

} // namespace duckdb
