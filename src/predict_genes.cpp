#include "predict_genes.hpp"

// Before the headers that declare duckdb::miint, which would shadow the global ::miint this one names unqualified.
#include "sequence_table_reader.hpp"

#include "align_bowtie2_daemon_common.hpp"
#include "catalog_utils.hpp"
#include "gpl_boundary/arrow_ipc.hpp"
#include "gpl_boundary/process.hpp"
#include "gpl_boundary/session.hpp"
#include "gpl_boundary_version.hpp"
#include "id_column_utils.hpp"
#include "miint_named_parameter.hpp"
#include "miint_streaming_query.hpp"

#include "duckdb/common/arrow/arrow_wrapper.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/column/column_data_collection.hpp"
#include "duckdb/common/vector/constant_vector.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/common/vector/string_vector.hpp"
#include "duckdb/function/table/arrow.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include "yyjson.hpp"

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace duckdb {

namespace {

namespace gb = ::duckdb::miint::gpl_boundary;
namespace yj = duckdb_yyjson;

constexpr const char *kFunction = "predict_genes";

// gpl-boundary 0.5.0 is the first with prodigal-train, proteins and the per-contig output (model GC included), at
// these schema versions. Both tools must report Prodigal kProdigalVersion (RequireProdigalVersion).
constexpr GplBoundaryRelease kMinRelease {0, 5, 0};
constexpr const char *kPredictTool = "prodigal";
constexpr const char *kTrainTool = "prodigal-train";
constexpr uint32_t kProdigalSchemaVersion = 2;
constexpr uint32_t kProdigalTrainSchemaVersion = 1;

// The status column's values.
constexpr const char *kStatusOk = "ok";
constexpr const char *kStatusNoGenes = "no_genes";
constexpr const char *kStatusTrainingFailed = "training_failed";

// Ends an error about gpl-boundary output this miint cannot read.
constexpr const char *kOutputMismatch =
    "This gpl-boundary is 0.5.0 or later, but returned output that this miint cannot read: upgrade miint, or set "
    "MIINT_GPL_BOUNDARY_PATH to a gpl-boundary of the release that this miint was built for";

// Bases per prodigal batch. Prodigal calls genes on each contig on its own, with the bin it picks (metagenomic mode)
// or the one model (single-genome mode), so this bounds memory, not results. (Training takes a genome whole.)
constexpr idx_t kBatchBases = idx_t(64) << 20;

// Output columns.
enum OutputColumn : idx_t {
	COL_GENOME_ID,
	COL_READ_ID,
	COL_GENE_INDEX,
	COL_POSITION,
	COL_STOP_POSITION,
	COL_STRAND,
	COL_PARTIAL_LEFT,
	COL_PARTIAL_RIGHT,
	COL_START_TYPE,
	COL_RBS_MOTIF,
	COL_RBS_SPACER,
	COL_SCORE,
	COL_CSCORE,
	COL_SSCORE,
	COL_RSCORE,
	COL_USCORE,
	COL_TSCORE,
	COL_CONFIDENCE,
	COL_GENE_GC_CONTENT,
	COL_PROTEIN,
	COL_CONTIG_LENGTH,
	COL_CONTIG_GC_CONTENT,
	COL_TRANSLATION_TABLE,
	COL_USES_SD,
	COL_META_BIN,
	COL_MODEL_DESC,
	COL_MODEL_GC_CONTENT,
	COL_STATUS,
	COL_ERROR,
	COL_COUNT
};

// The daemon's outputs, as gpl-boundary 0.5.0's --describe prodigal lists them. A batch's outputs are checked
// against these by name and type, so a reordered or retyped column fails loudly.
const vector<string> kGenesNames = {"seq_index",  "seq_name",     "gene_index",    "begin",      "end",
                                    "strand",     "partial_left", "partial_right", "start_type", "score",
                                    "cscore",     "sscore",       "rscore",        "uscore",     "tscore",
                                    "confidence", "gc_cont",      "rbs_motif",     "rbs_spacer", "protein"};
const vector<LogicalType> kGenesTypes = {
    LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER,
    LogicalType::INTEGER, LogicalType::BOOLEAN, LogicalType::BOOLEAN, LogicalType::INTEGER, LogicalType::DOUBLE,
    LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE,
    LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR};
const vector<string> kContigsNames = {"seq_index",         "seq_name", "length",   "gc_content", "model_gc_content",
                                      "translation_table", "uses_sd",  "meta_bin", "model_desc", "n_genes"};
const vector<LogicalType> kContigsTypes = {
    LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::DOUBLE,  LogicalType::DOUBLE,
    LogicalType::INTEGER, LogicalType::BOOLEAN, LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::INTEGER};

// The outputs' columns, by position in kGenesNames and kContigsNames.
enum GenesColumn : idx_t {
	G_SEQ_INDEX,
	G_SEQ_NAME,
	G_GENE_INDEX,
	G_BEGIN,
	G_END,
	G_STRAND,
	G_PARTIAL_LEFT,
	G_PARTIAL_RIGHT,
	G_START_TYPE,
	G_SCORE,
	G_CSCORE,
	G_SSCORE,
	G_RSCORE,
	G_USCORE,
	G_TSCORE,
	G_CONFIDENCE,
	G_GC_CONT,
	G_RBS_MOTIF,
	G_RBS_SPACER,
	G_PROTEIN
};
enum ContigsColumn : idx_t {
	C_SEQ_INDEX,
	C_SEQ_NAME,
	C_LENGTH,
	C_GC_CONTENT,
	C_MODEL_GC_CONTENT,
	C_TRANSLATION_TABLE,
	C_USES_SD,
	C_META_BIN,
	C_MODEL_DESC,
	C_N_GENES
};

// Prodigal's start codon types (the genes output's start_type), as its GFF names them.
const char *const kStartTypes[] = {"ATG", "GTG", "TTG", "Edge"};

// =============================================================================
// Bind
// =============================================================================

struct PredictGenesBindData : public TableFunctionData {
	std::string relation;
	std::string genome_column;
	std::string order_column = "sequence_index";
	bool order_set = false;
	// Chunked input: each contig is the rows of one key (read_id or feature_idx), its sequence their chunk_data in
	// chunk_index order; mapping_table maps a key to its genome.
	bool chunked = false;
	std::string key_column;
	std::string mapping_table;
	LogicalType chunk_data_type;
	bool meta_mode = false;
	int32_t trans_table = 11;
	bool trans_table_set = false;
	bool closed_ends = false;
	bool mask_regions = false;
	std::string training_path;
	bool has_sequence2 = false;
	LogicalType genome_type;
	LogicalType read_id_type;
	vector<LogicalType> return_types;

	unique_ptr<FunctionData> Copy() const override {
		return make_uniq<PredictGenesBindData>(*this);
	}

	bool Equals(const FunctionData &other_p) const override {
		const auto &o = other_p.Cast<PredictGenesBindData>();
		return relation == o.relation && genome_column == o.genome_column && order_column == o.order_column &&
		       order_set == o.order_set && chunked == o.chunked && key_column == o.key_column &&
		       mapping_table == o.mapping_table && chunk_data_type == o.chunk_data_type && meta_mode == o.meta_mode &&
		       trans_table == o.trans_table && trans_table_set == o.trans_table_set && closed_ends == o.closed_ends &&
		       mask_regions == o.mask_regions && training_path == o.training_path && has_sequence2 == o.has_sequence2 &&
		       genome_type == o.genome_type && read_id_type == o.read_id_type;
	}
};

void ReadNamedParameters(const named_argument_map_t &params, PredictGenesBindData &data) {
	for (const auto &kv : params) {
		const auto name = StringUtil::Lower(kv.first.GetIdentifierName());
		const Value &value = kv.second;
		if (value.IsNull()) {
			throw BinderException("%s: %s must not be NULL", kFunction, name);
		}
		if (name == "genome_id") {
			data.genome_column = StringValue::Get(value);
		} else if (name == "mapping_table") {
			data.mapping_table = StringValue::Get(value);
		} else if (name == "method") {
			const auto method = StringValue::Get(value);
			if (!StringUtil::CIEquals(method, kPredictTool)) {
				throw BinderException("%s: method '%s' is not supported; the only method is '%s'", kFunction, method,
				                      kPredictTool);
			}
		} else if (name == "order_by") {
			data.order_column = StringValue::Get(value);
			data.order_set = true;
		} else if (name == "meta_mode") {
			data.meta_mode = BooleanValue::Get(value);
		} else if (name == "trans_table") {
			data.trans_table = IntegerValue::Get(value);
			data.trans_table_set = true;
		} else if (name == "closed_ends") {
			data.closed_ends = BooleanValue::Get(value);
		} else if (name == "mask_regions") {
			data.mask_regions = BooleanValue::Get(value);
		} else if (name == "training_path") {
			data.training_path = StringValue::Get(value);
			if (data.training_path.empty()) {
				throw BinderException("%s: training_path must not be empty", kFunction);
			}
		}
	}
	if (data.meta_mode && data.trans_table != 11) {
		throw BinderException("%s: metagenomic mode (meta_mode := true) takes trans_table 11 only: its models choose "
		                      "table 11 or 4 per contig. Other tables need single-genome mode",
		                      kFunction);
	}
	if (data.meta_mode && !data.training_path.empty()) {
		throw BinderException("%s: training_path is a single-genome model; metagenomic mode (meta_mode := true) "
		                      "uses Prodigal's own models",
		                      kFunction);
	}
}

// The column names of `relation`, lowercased. A LIMIT 0 probe works for tables and views.
std::set<std::string> ColumnNames(ClientContext &context, const std::string &relation) {
	auto conn = MakeReadOnlyHelperConnection(context);
	auto probe = conn.Query("SELECT * FROM " + SQLIdentifier::ToString(relation) + " LIMIT 0");
	if (probe->HasError()) {
		throw BinderException("%s: cannot read '%s': %s", kFunction, relation, probe->GetError());
	}
	std::set<std::string> names;
	for (const auto &n : probe->GetNames()) {
		names.insert(StringUtil::Lower(n.GetIdentifierName()));
	}
	return names;
}

// Contigs, one per row: read_id and sequence1, grouped into genomes by the genome_id column.
void BindContigs(ClientContext &context, PredictGenesBindData &data) {
	if (data.genome_column.empty()) {
		throw BinderException("%s: genome_id (the name of the column that groups contigs into genomes) is required",
		                      kFunction);
	}
	if (!data.mapping_table.empty()) {
		throw BinderException("%s: mapping_table is for chunked input (a relation with chunk_data); '%s' has a "
		                      "contig per row, grouped into genomes by genome_id",
		                      kFunction, data.relation);
	}
	// read_id (VARCHAR, BIGINT or UUID) and sequence1 (VARCHAR).
	const auto schema = ValidateSequenceTableSchema(context, data.relation, /*allow_bigint=*/true);
	data.read_id_type = schema.id_type;
	data.has_sequence2 = schema.has_sequence2;

	auto conn = MakeReadOnlyHelperConnection(context);
	auto probe = conn.Query("SELECT " + SQLIdentifier::ToString(data.genome_column) + ", " +
	                        SQLIdentifier::ToString(data.order_column) + " FROM " +
	                        SQLIdentifier::ToString(data.relation) + " LIMIT 0");
	if (probe->HasError()) {
		throw BinderException("%s: cannot read the genome_id column '%s' or the order_by column '%s' of '%s' "
		                      "(order_by defaults to sequence_index, the column read_fastx adds; name another "
		                      "integer column for a relation without it): %s",
		                      kFunction, data.genome_column, data.order_column, data.relation, probe->GetError());
	}
	data.genome_type = probe->GetTypes()[0];
	if (!IsAllowedIdType(data.genome_type)) {
		throw BinderException("%s: genome_id column '%s' must be %s (got %s)", kFunction, data.genome_column,
		                      AllowedIdTypeList(), data.genome_type.ToString());
	}
	if (!probe->GetTypes()[1].IsIntegral()) {
		throw BinderException("%s: order_by column '%s' must be an integer type (got %s)", kFunction, data.order_column,
		                      probe->GetTypes()[1].ToString());
	}
}

// Chunked contigs, as genomes are stored in blocks: rows of (read_id or feature_idx, chunk_index, chunk_data), a
// contig's sequence its rows' chunk_data in chunk_index order. mapping_table maps each key to its genome_id.
void BindChunked(ClientContext &context, PredictGenesBindData &data, const std::set<std::string> &columns) {
	if (data.mapping_table.empty()) {
		throw BinderException("%s: '%s' is chunked (it has chunk_data), so mapping_table, a table or view of "
		                      "(read_id or feature_idx, genome_id) that maps each contig to its genome, is required",
		                      kFunction, data.relation);
	}
	if (!data.genome_column.empty()) {
		throw BinderException("%s: genome_id is for a relation of a contig per row; '%s' is chunked (it has "
		                      "chunk_data), so its genomes come from mapping_table",
		                      kFunction, data.relation);
	}
	if (data.order_set) {
		throw BinderException("%s: order_by is for a relation of a contig per row; a chunked genome's contigs are "
		                      "in the order of their keys",
		                      kFunction);
	}
	const bool has_read_id = columns.count("read_id") > 0;
	const bool has_feature_idx = columns.count("feature_idx") > 0;
	if (has_read_id == has_feature_idx) {
		throw BinderException("%s: chunked '%s' must have one key column, read_id or feature_idx, naming each "
		                      "chunk's contig (it has %s)",
		                      kFunction, data.relation, has_read_id ? "both" : "neither");
	}
	data.key_column = has_read_id ? "read_id" : "feature_idx";
	const auto key = SQLIdentifier::ToString(data.key_column);

	auto conn = MakeReadOnlyHelperConnection(context);
	auto chunks = conn.Query("SELECT " + key + ", chunk_index, chunk_data FROM " +
	                         SQLIdentifier::ToString(data.relation) + " LIMIT 0");
	if (chunks->HasError()) {
		throw BinderException("%s: chunked '%s' must have %s, chunk_index and chunk_data: %s", kFunction, data.relation,
		                      data.key_column, chunks->GetError());
	}
	data.read_id_type = chunks->GetTypes()[0];
	if (!IsAllowedIdType(data.read_id_type)) {
		throw BinderException("%s: %s column of '%s' must be %s (got %s)", kFunction, data.key_column, data.relation,
		                      AllowedIdTypeList(), data.read_id_type.ToString());
	}
	if (!chunks->GetTypes()[1].IsIntegral()) {
		throw BinderException("%s: chunk_index column of '%s' must be an integer type (got %s)", kFunction,
		                      data.relation, chunks->GetTypes()[1].ToString());
	}
	data.chunk_data_type = chunks->GetTypes()[2];
	if (data.chunk_data_type.id() != LogicalTypeId::VARCHAR && data.chunk_data_type.id() != LogicalTypeId::BLOB) {
		throw BinderException("%s: chunk_data column of '%s' must be VARCHAR or BLOB (got %s)", kFunction,
		                      data.relation, data.chunk_data_type.ToString());
	}

	auto mapping =
	    conn.Query("SELECT " + key + ", genome_id FROM " + SQLIdentifier::ToString(data.mapping_table) + " LIMIT 0");
	if (mapping->HasError()) {
		throw BinderException("%s: mapping_table '%s' must have %s and genome_id: %s", kFunction, data.mapping_table,
		                      data.key_column, mapping->GetError());
	}
	if (mapping->GetTypes()[0] != data.read_id_type) {
		throw BinderException("%s: %s is %s in '%s' but %s in mapping_table '%s'; the two must be the same type",
		                      kFunction, data.key_column, data.read_id_type.ToString(), data.relation,
		                      mapping->GetTypes()[0].ToString(), data.mapping_table);
	}
	data.genome_type = mapping->GetTypes()[1];
	if (!IsAllowedIdType(data.genome_type)) {
		throw BinderException("%s: genome_id column of mapping_table '%s' must be %s (got %s)", kFunction,
		                      data.mapping_table, AllowedIdTypeList(), data.genome_type.ToString());
	}
}

unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input, vector<LogicalType> &return_types,
                              vector<Identifier> &names) {
	if (input.inputs.empty() || input.inputs[0].IsNull()) {
		throw BinderException("%s: a non-NULL relation name is required", kFunction);
	}
	auto data = make_uniq<PredictGenesBindData>();
	data->relation = StringValue::Get(input.inputs[0]);
	RejectCTERelationName(input, data->relation);
	ReadNamedParameters(input.named_parameters, *data);

	const auto columns = ColumnNames(context, data->relation);
	data->chunked = columns.count("chunk_data") > 0;
	if (data->chunked) {
		if (!data->mapping_table.empty()) {
			RejectCTERelationName(input, data->mapping_table);
		}
		BindChunked(context, *data, columns);
	} else {
		BindContigs(context, *data);
	}

	names = {"genome_id",
	         "read_id",
	         "gene_index",
	         "position",
	         "stop_position",
	         "strand",
	         "partial_left",
	         "partial_right",
	         "start_type",
	         "rbs_motif",
	         "rbs_spacer",
	         "score",
	         "cscore",
	         "sscore",
	         "rscore",
	         "uscore",
	         "tscore",
	         "confidence",
	         "gene_gc_content",
	         "protein",
	         "contig_length",
	         "contig_gc_content",
	         "translation_table",
	         "uses_sd",
	         "meta_bin",
	         "model_desc",
	         "model_gc_content",
	         "status",
	         "error"};
	return_types = {
	    data->genome_type,    data->read_id_type,   LogicalType::INTEGER, LogicalType::BIGINT,  LogicalType::BIGINT,
	    LogicalType::VARCHAR, LogicalType::BOOLEAN, LogicalType::BOOLEAN, LogicalType::VARCHAR, LogicalType::VARCHAR,
	    LogicalType::VARCHAR, LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE,
	    LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::VARCHAR,
	    LogicalType::BIGINT,  LogicalType::DOUBLE,  LogicalType::INTEGER, LogicalType::BOOLEAN, LogicalType::INTEGER,
	    LogicalType::VARCHAR, LogicalType::DOUBLE,  LogicalType::VARCHAR, LogicalType::VARCHAR};
	D_ASSERT(names.size() == COL_COUNT && return_types.size() == COL_COUNT);
	data->return_types = return_types;
	return std::move(data);
}

// =============================================================================
// Global state
// =============================================================================

// One contig as read from the relation.
struct Contig {
	Value genome;
	Value read_id;
	std::string name; // read_id as text, Prodigal's sequence name
	std::string sequence;
};

// A model file and its SHA-256, for prodigal's training_path / training_sha256.
struct Model {
	std::string path;
	std::string sha256;
	bool closed_ends = false;
	bool mask_regions = false;
};

struct PredictGenesGlobalState : public GlobalTableFunctionState {
	std::unique_ptr<gb::Session> session;

	// The relation, streamed one genome at a time in (genome, order_by, read_id) order.
	unique_ptr<Connection> conn;
	unique_ptr<StreamingQuery> stream;
	unique_ptr<DataChunk> chunk;
	idx_t chunk_row = 0;
	bool has_pending = false;
	Contig pending;

	// Single-genome mode: a model the user gave (fixed), else where each genome's model is written.
	bool have_fixed_model = false;
	Model fixed_model;
	std::string train_path;
	bool write_user_model = false;
	std::string temp_dir;

	// Metagenomic mode streams a genome in batches; its gene count spans them.
	bool in_genome = false;
	Value current_genome;
	idx_t current_genome_genes = 0;

	// Rows ready to emit.
	unique_ptr<ColumnDataCollection> rows;
	ColumnDataScanState scan;
	bool scanning = false;

	idx_t MaxThreads() const override {
		return 1;
	}

	~PredictGenesGlobalState() override {
		if (session) {
			try {
				session->Shutdown();
			} catch (...) {
				// Destructor must not throw.
			}
			session.reset();
		}
		if (!temp_dir.empty()) {
			std::error_code ec;
			std::filesystem::remove_all(temp_dir, ec); // Nothing useful to do on failure.
		}
	}
};

// Refuse a gpl-boundary whose `--version` does not list prodigal and prodigal-train at kProdigalVersion, the version
// COPY FORMAT GFF's DIALECT 'prodigal' names. Probed on every call, so a binary replaced in place is seen by the next
// query.
void RequireProdigalVersion(const std::string &gpl_path) {
	std::string text;
	try {
		const auto result = gb::RunCapture({gpl_path, "--version"});
		if (result.spawned && result.exited_normally && result.exit_code == 0) {
			text = result.stdout_text;
		}
	} catch (...) {
		// Undetermined, so refused below.
	}
	std::map<std::string, std::string> versions;
	std::unique_ptr<yj::yyjson_doc, decltype(&yj::yyjson_doc_free)> doc(yj::yyjson_read(text.c_str(), text.size(), 0),
	                                                                    &yj::yyjson_doc_free);
	if (doc) {
		auto *tools = yj::yyjson_obj_get(yj::yyjson_doc_get_root(doc.get()), "tools");
		size_t idx, max;
		yj::yyjson_val *tool;
		yyjson_arr_foreach(tools, idx, max, tool) {
			const char *name = yj::yyjson_get_str(yj::yyjson_obj_get(tool, "name"));
			const char *version = yj::yyjson_get_str(yj::yyjson_obj_get(tool, "version"));
			if (name && version) {
				versions[name] = version;
			}
		}
	}
	for (const char *tool : {kPredictTool, kTrainTool}) {
		auto it = versions.find(tool);
		if (it != versions.end() && it->second == kProdigalVersion) {
			continue;
		}
		const std::string found = it == versions.end() ? std::string("no version for it") : "version " + it->second;
		throw IOException("%s requires a gpl-boundary whose %s is version %s (the version COPY FORMAT GFF's DIALECT "
		                  "'prodigal' names), "
		                  "but the binary at %s reports %s. Point "
		                  "MIINT_GPL_BOUNDARY_PATH at a gpl-boundary that ships %s %s, or use a miint release that "
		                  "supports this one.",
		                  kFunction, tool, kProdigalVersion, gpl_path, found, tool, kProdigalVersion);
	}
}

// Refuse a trans_table that gpl-boundary's prodigal-train does not take, as its `--describe` lists them
// (config_params' trans_table allowed_values), rather than keep a list of its own.
void RequireTransTable(const std::string &gpl_path, int32_t trans_table) {
	std::string text;
	bool ran = false;
	try {
		const auto result = gb::RunCapture({gpl_path, "--describe", kTrainTool});
		ran = result.spawned && result.exited_normally && result.exit_code == 0;
		text = result.stdout_text;
	} catch (...) {
		// Undetermined, so refused below.
	}
	if (!ran) {
		throw IOException("%s: could not run `%s --describe %s` to read the translation tables it takes", kFunction,
		                  gpl_path, kTrainTool);
	}
	std::unique_ptr<yj::yyjson_doc, decltype(&yj::yyjson_doc_free)> doc(yj::yyjson_read(text.c_str(), text.size(), 0),
	                                                                    &yj::yyjson_doc_free);
	yj::yyjson_val *allowed = nullptr;
	if (doc) {
		auto *params = yj::yyjson_obj_get(yj::yyjson_doc_get_root(doc.get()), "config_params");
		size_t idx, max;
		yj::yyjson_val *param;
		yyjson_arr_foreach(params, idx, max, param) {
			const char *name = yj::yyjson_get_str(yj::yyjson_obj_get(param, "name"));
			if (name && std::string(name) == "trans_table") {
				allowed = yj::yyjson_obj_get(param, "allowed_values");
			}
		}
	}
	if (!yj::yyjson_is_arr(allowed)) {
		throw IOException("%s: cannot read the translation tables gpl-boundary's %s takes from `%s --describe %s`. %s",
		                  kFunction, kTrainTool, gpl_path, kTrainTool, kOutputMismatch);
	}
	std::string list;
	size_t idx, max;
	yj::yyjson_val *value;
	yyjson_arr_foreach(allowed, idx, max, value) {
		if (yj::yyjson_is_int(value) && yj::yyjson_get_int(value) == trans_table) {
			return;
		}
		list += (list.empty() ? "" : ", ") + std::to_string(yj::yyjson_get_int(value));
	}
	throw InvalidInputException("%s: trans_table %d is not one gpl-boundary's %s takes (%s)", kFunction, trans_table,
	                            kTrainTool, list);
}

std::string MakeTempModelDir() {
	const char *tmp = std::getenv("TMPDIR");
	if (!tmp || !*tmp) {
		tmp = "/tmp";
	}
	std::string tmpl = std::string(tmp) + "/miint-prodigal-XXXXXX";
	std::vector<char> buf(tmpl.begin(), tmpl.end());
	buf.push_back('\0');
	if (::mkdtemp(buf.data()) == nullptr) {
		throw IOException("%s: failed to create a temporary directory for Prodigal's models under '%s' (errno=%d)",
		                  kFunction, std::string(tmp), errno);
	}
	return std::string(buf.data());
}

// The SHA-256 of a model file, in lowercase hexadecimal as prodigal-train reports it.
std::string HashFile(Connection &conn, const std::string &path) {
	auto result = conn.Query("SELECT sha256(content) FROM read_blob(" + SQLString::ToString(path) + ")");
	if (result->HasError()) {
		throw IOException("%s: cannot read the model at training_path '%s': %s", kFunction, path, result->GetError());
	}
	auto chunk = result->Fetch();
	if (!chunk || chunk->size() != 1) {
		throw IOException("%s: training_path '%s' does not name one model file", kFunction, path);
	}
	return chunk->GetValue(0, 0).ToString();
}

// =============================================================================
// Reading the relation
// =============================================================================

std::string GenomeText(const Value &genome) {
	return genome.ToString();
}

// Prodigal drops any character outside 'A'..'z', shifting every coordinate after it, so gpl-boundary refuses one;
// refusing it here names the contig and its genome before a daemon round trip.
void CheckContig(const Contig &c, const PredictGenesBindData &data) {
	if (c.sequence.empty()) {
		throw InvalidInputException("%s: contig '%s' of genome '%s' has an empty %s", kFunction, c.name,
		                            GenomeText(c.genome), data.chunked ? "chunk_data" : "sequence1");
	}
	for (idx_t i = 0; i < c.sequence.size(); i++) {
		const char ch = c.sequence[i];
		if (ch < 'A' || ch > 'z') {
			const auto code = static_cast<unsigned>(static_cast<unsigned char>(ch));
			throw InvalidInputException("%s: contig '%s' of genome '%s' has character 0x%02x at position %llu, which "
			                            "Prodigal would drop (shifting every coordinate after it); only characters "
			                            "'A'..'z' are allowed",
			                            kFunction, c.name, GenomeText(c.genome), code,
			                            static_cast<unsigned long long>(i + 1));
		}
	}
}

bool FetchContig(PredictGenesGlobalState &g, const PredictGenesBindData &data, Contig &out) {
	while (!g.chunk || g.chunk_row >= g.chunk->size()) {
		g.chunk = g.stream->Fetch();
		g.chunk_row = 0;
		if (!g.chunk) {
			if (g.stream->HasError()) {
				throw InvalidInputException("%s: cannot read '%s': %s", kFunction, data.relation, g.stream->GetError());
			}
			return false;
		}
	}
	auto &chunk = *g.chunk;
	const idx_t row = g.chunk_row++;
	out.genome = chunk.GetValue(0, row);
	out.read_id = chunk.GetValue(1, row);
	const Value sequence = chunk.GetValue(2, row);
	if (out.read_id.IsNull()) {
		throw InvalidInputException("%s: a contig of genome '%s' has a NULL read_id", kFunction,
		                            GenomeText(out.genome));
	}
	out.name = out.read_id.ToString();
	if (sequence.IsNull()) {
		throw InvalidInputException("%s: contig '%s' of genome '%s' has a NULL sequence1", kFunction, out.name,
		                            GenomeText(out.genome));
	}
	if (data.has_sequence2 && !chunk.GetValue(3, row).IsNull()) {
		throw InvalidInputException("%s: contig '%s' of genome '%s' has a sequence2; contigs are single sequences",
		                            kFunction, out.name, GenomeText(out.genome));
	}
	out.sequence = StringValue::Get(sequence);
	CheckContig(out, data);
	return true;
}

// The next contig, or false at the end of the relation.
bool NextContig(PredictGenesGlobalState &g, const PredictGenesBindData &data, Contig &out) {
	if (g.has_pending) {
		out = std::move(g.pending);
		g.has_pending = false;
		return true;
	}
	return FetchContig(g, data, out);
}

// Whether the next contig belongs to `genome` (false at the end of the relation).
bool NextIsSameGenome(PredictGenesGlobalState &g, const PredictGenesBindData &data, const Value &genome) {
	if (!g.has_pending) {
		g.has_pending = FetchContig(g, data, g.pending);
	}
	return g.has_pending && g.pending.genome == genome;
}

// =============================================================================
// Calling the daemon
// =============================================================================

std::vector<uint8_t> EncodeContigs(const std::vector<Contig> &contigs, idx_t begin, idx_t end) {
	bt2_daemon::LoadedSubjects input;
	input.names.reserve(end - begin);
	input.sequences.reserve(end - begin);
	for (idx_t i = begin; i < end; i++) {
		input.names.push_back(contigs[i].name);
		input.sequences.push_back(contigs[i].sequence);
	}
	return bt2_daemon::BuildSubjectsIpc(input, kFunction);
}

gb::SubmitResult Submit(PredictGenesGlobalState &g, const std::string &tool, const std::string &config,
                        const std::vector<uint8_t> &input) {
	try {
		return g.session->Submit(tool, config, input.data(), input.size());
	} catch (const gb::BatchRefused &) {
		throw; // The caller decides what a refusal means.
	} catch (const std::runtime_error &e) {
		throw IOException("%s: %s", kFunction, e.what());
	}
}

// Decode one output's Arrow IPC stream into chunks, after checking its columns are `names` / `types`.
std::vector<unique_ptr<DataChunk>> DecodeOutput(ClientContext &context, const gb::SubmitResult &result,
                                                const std::string &label, const vector<string> &names,
                                                const vector<LogicalType> &types) {
	const gb::SubmitOutput *output = nullptr;
	for (const auto &o : result.outputs) {
		if (o.label == label) {
			output = &o;
		}
	}
	if (!output) {
		throw IOException("%s: gpl-boundary's prodigal returned no '%s' output. %s", kFunction, label, kOutputMismatch);
	}
	std::vector<unique_ptr<DataChunk>> chunks;
	try {
		gb::IpcStreamDecoder decoder(output->bytes(), output->size_bytes());
		ArrowSchemaWrapper schema;
		decoder.GetSchema(&schema.arrow_schema);
		ArrowTableSchema table;
		ArrowTableFunction::PopulateArrowTableSchema(context, table, schema.arrow_schema);
		vector<string> got_names;
		for (const auto &n : table.GetNames()) {
			got_names.push_back(n);
		}
		if (got_names != names || table.GetTypes() != types) {
			throw IOException("%s: gpl-boundary's prodigal '%s' output has columns (%s), not (%s). %s", kFunction,
			                  label, StringUtil::Join(got_names, ", "), StringUtil::Join(names, ", "), kOutputMismatch);
		}
		const auto &columns = table.GetColumns();
		while (true) {
			auto batch = make_shared_ptr<ArrowArrayWrapper>();
			if (!decoder.NextBatch(&batch->arrow_array)) {
				break;
			}
			const auto length = static_cast<idx_t>(batch->arrow_array.length);
			for (idx_t offset = 0; offset < length; offset += STANDARD_VECTOR_SIZE) {
				const idx_t n = MinValue<idx_t>(STANDARD_VECTOR_SIZE, length - offset);
				auto chunk = make_uniq<DataChunk>();
				chunk->Initialize(context, types);
				chunk->SetChildCardinality(n);
				for (idx_t c = 0; c < types.size(); c++) {
					auto &array = *batch->arrow_array.children[c];
					ArrowArrayScanState state(context);
					state.owned_data = batch;
					ArrowToDuckDBConversion::SetValidityMask(chunk->data[c], array, offset, n,
					                                         batch->arrow_array.offset, -1);
					ArrowToDuckDBConversion::ColumnArrowToDuckDB(chunk->data[c], array, offset, state, n,
					                                             *columns.at(c));
					FlatVector::SetSize(chunk->data[c], count_t(n));
				}
				chunk->Flatten();
				chunks.push_back(std::move(chunk));
			}
		}
	} catch (const Exception &) {
		throw;
	} catch (const std::runtime_error &e) {
		throw IOException("%s: cannot decode gpl-boundary's prodigal '%s' output: %s", kFunction, label, e.what());
	}
	return chunks;
}

// The per-contig values a gene row repeats.
struct ContigInfo {
	int64_t length = 0;
	double gc_content = 0;
	int32_t translation_table = 0;
	bool uses_sd = false;
	bool has_meta_bin = false;
	int32_t meta_bin = 0;
	bool has_model_desc = false;
	std::string model_desc;
	double model_gc_content = 0;
	int32_t n_genes = 0;
};

// A row of one genome's status, with every gene and contig column NULL.
void AppendStatusRow(ClientContext &context, PredictGenesGlobalState &g, const Value &genome, const char *status,
                     const std::string &error) {
	DataChunk out;
	out.Initialize(context, g.rows->Types());
	for (idx_t c = 0; c < COL_COUNT; c++) {
		FlatVector::SetNull(out.data[c], 0, true);
	}
	out.data[COL_GENOME_ID].SetValue(0, genome);
	out.data[COL_STATUS].SetValue(0, Value(status));
	if (!error.empty()) {
		out.data[COL_ERROR].SetValue(0, Value(error));
	}
	out.SetChildCardinality(1);
	g.rows->Append(out);
}

// A row of a contig without genes: its contig columns, every gene column NULL. COPY FORMAT GFF writes its header
// lines from it, as the CLI does for such a contig.
void AppendBarrenContigRow(ClientContext &context, PredictGenesGlobalState &g, const Contig &c, const ContigInfo &ci) {
	DataChunk out;
	out.Initialize(context, g.rows->Types());
	for (idx_t col = 0; col < COL_COUNT; col++) {
		FlatVector::SetNull(out.data[col], 0, true);
	}
	out.data[COL_GENOME_ID].SetValue(0, c.genome);
	out.data[COL_READ_ID].SetValue(0, c.read_id);
	out.data[COL_CONTIG_LENGTH].SetValue(0, Value::BIGINT(ci.length));
	out.data[COL_CONTIG_GC_CONTENT].SetValue(0, Value::DOUBLE(ci.gc_content));
	out.data[COL_TRANSLATION_TABLE].SetValue(0, Value::INTEGER(ci.translation_table));
	out.data[COL_USES_SD].SetValue(0, Value::BOOLEAN(ci.uses_sd));
	if (ci.has_meta_bin) {
		out.data[COL_META_BIN].SetValue(0, Value::INTEGER(ci.meta_bin));
	}
	if (ci.has_model_desc) {
		out.data[COL_MODEL_DESC].SetValue(0, Value(ci.model_desc));
	}
	out.data[COL_MODEL_GC_CONTENT].SetValue(0, Value::DOUBLE(ci.model_gc_content));
	out.data[COL_STATUS].SetValue(0, Value(kStatusNoGenes));
	out.SetChildCardinality(1);
	g.rows->Append(out);
}

// Call genes on contigs [begin, end) with `config`, and append a row per gene, and one per contig without genes,
// in contig order. Returns the number of genes.
idx_t Predict(ClientContext &context, PredictGenesGlobalState &g, const std::string &config,
              const std::vector<Contig> &contigs, idx_t begin, idx_t end) {
	gb::SubmitResult result;
	try {
		result = Submit(g, kPredictTool, config, EncodeContigs(contigs, begin, end));
	} catch (const gb::BatchRefused &e) {
		// Input is checked before it is sent, so this is the call itself refused (a model file it cannot use, say).
		throw InvalidInputException("%s: prodigal refused the call: %s", kFunction, e.reason());
	}
	const auto gene_chunks = DecodeOutput(context, result, "genes", kGenesNames, kGenesTypes);
	const auto contig_chunks = DecodeOutput(context, result, "contigs", kContigsNames, kContigsTypes);

	// seq_index is the contig's 1-based position in the batch.
	std::vector<ContigInfo> info(end - begin);
	for (const auto &chunk : contig_chunks) {
		const auto *seq_index = FlatVector::GetData<int32_t>(chunk->data[C_SEQ_INDEX]);
		const auto *length = FlatVector::GetData<int32_t>(chunk->data[C_LENGTH]);
		const auto *gc = FlatVector::GetData<double>(chunk->data[C_GC_CONTENT]);
		// The GC of the model that called the contig: its bin's, or the single-genome model's (the CLI's Model Data
		// gc_cont / 100)
		const auto *model_gc = FlatVector::GetData<double>(chunk->data[C_MODEL_GC_CONTENT]);
		const auto *table = FlatVector::GetData<int32_t>(chunk->data[C_TRANSLATION_TABLE]);
		const auto *uses_sd = FlatVector::GetData<bool>(chunk->data[C_USES_SD]);
		const auto *meta_bin = FlatVector::GetData<int32_t>(chunk->data[C_META_BIN]);
		const auto *model_desc = FlatVector::GetData<string_t>(chunk->data[C_MODEL_DESC]);
		const auto *contig_genes = FlatVector::GetData<int32_t>(chunk->data[C_N_GENES]);
		for (idx_t r = 0; r < chunk->size(); r++) {
			if (seq_index[r] < 1 || static_cast<idx_t>(seq_index[r]) > info.size()) {
				throw IOException("%s: gpl-boundary's prodigal returned contig seq_index %d for a batch of %llu",
				                  kFunction, seq_index[r], static_cast<unsigned long long>(info.size()));
			}
			auto &ci = info[seq_index[r] - 1];
			ci.length = length[r];
			ci.gc_content = gc[r];
			ci.translation_table = table[r];
			ci.uses_sd = uses_sd[r];
			ci.has_meta_bin = !FlatVector::IsNull(chunk->data[C_META_BIN], r);
			ci.meta_bin = ci.has_meta_bin ? meta_bin[r] : 0;
			ci.has_model_desc = !FlatVector::IsNull(chunk->data[C_MODEL_DESC], r);
			ci.model_desc = ci.has_model_desc ? model_desc[r].GetString() : std::string();
			ci.model_gc_content = model_gc[r];
			ci.n_genes = contig_genes[r];
		}
	}

	// Contigs [next_contig, upto) have no genes, as the genes output (in contig order) skips them.
	idx_t next_contig = 0;
	const auto append_barren_contigs = [&](idx_t upto) {
		for (; next_contig < upto; next_contig++) {
			if (info[next_contig].n_genes != 0) {
				throw IOException(
				    "%s: gpl-boundary's prodigal reported %d genes on contig %llu of a batch but returned "
				    "none. %s",
				    kFunction, info[next_contig].n_genes, static_cast<unsigned long long>(next_contig + 1),
				    kOutputMismatch);
			}
			AppendBarrenContigRow(context, g, contigs[begin + next_contig], info[next_contig]);
		}
	};

	idx_t n_genes = 0;
	DataChunk out;
	out.Initialize(context, g.rows->Types());
	for (const auto &chunk : gene_chunks) {
		out.Reset();
		const idx_t n = chunk->size();
		const auto *seq_index = FlatVector::GetData<int32_t>(chunk->data[G_SEQ_INDEX]);
		const auto *gene_index = FlatVector::GetData<int32_t>(chunk->data[G_GENE_INDEX]);
		const auto *gene_begin = FlatVector::GetData<int32_t>(chunk->data[G_BEGIN]);
		const auto *gene_end = FlatVector::GetData<int32_t>(chunk->data[G_END]);
		const auto *strand = FlatVector::GetData<int32_t>(chunk->data[G_STRAND]);
		const auto *partial_left = FlatVector::GetData<bool>(chunk->data[G_PARTIAL_LEFT]);
		const auto *partial_right = FlatVector::GetData<bool>(chunk->data[G_PARTIAL_RIGHT]);
		const auto *start_type = FlatVector::GetData<int32_t>(chunk->data[G_START_TYPE]);

		auto *o_gene_index = FlatVector::GetDataMutable<int32_t>(out.data[COL_GENE_INDEX]);
		auto *o_position = FlatVector::GetDataMutable<int64_t>(out.data[COL_POSITION]);
		auto *o_stop = FlatVector::GetDataMutable<int64_t>(out.data[COL_STOP_POSITION]);
		auto *o_partial_left = FlatVector::GetDataMutable<bool>(out.data[COL_PARTIAL_LEFT]);
		auto *o_partial_right = FlatVector::GetDataMutable<bool>(out.data[COL_PARTIAL_RIGHT]);
		auto *o_length = FlatVector::GetDataMutable<int64_t>(out.data[COL_CONTIG_LENGTH]);
		auto *o_contig_gc = FlatVector::GetDataMutable<double>(out.data[COL_CONTIG_GC_CONTENT]);
		auto *o_table = FlatVector::GetDataMutable<int32_t>(out.data[COL_TRANSLATION_TABLE]);
		auto *o_uses_sd = FlatVector::GetDataMutable<bool>(out.data[COL_USES_SD]);
		auto *o_meta_bin = FlatVector::GetDataMutable<int32_t>(out.data[COL_META_BIN]);
		auto *o_model_gc = FlatVector::GetDataMutable<double>(out.data[COL_MODEL_GC_CONTENT]);

		for (idx_t r = 0; r < n; r++) {
			if (seq_index[r] < 1 || static_cast<idx_t>(seq_index[r]) > info.size()) {
				throw IOException("%s: gpl-boundary's prodigal returned gene seq_index %d for a batch of %llu",
				                  kFunction, seq_index[r], static_cast<unsigned long long>(info.size()));
			}
			const idx_t contig = static_cast<idx_t>(seq_index[r] - 1);
			const auto &c = contigs[begin + contig];
			const auto &ci = info[contig];
			out.data[COL_GENOME_ID].SetValue(r, c.genome);
			out.data[COL_READ_ID].SetValue(r, c.read_id);
			o_gene_index[r] = gene_index[r];
			// Prodigal's begin/end are 1-based and inclusive; miint's are 1-based half-open.
			o_position[r] = gene_begin[r];
			o_stop[r] = static_cast<int64_t>(gene_end[r]) + 1;
			out.data[COL_STRAND].SetValue(r, Value(strand[r] > 0 ? "+" : "-"));
			o_partial_left[r] = partial_left[r];
			o_partial_right[r] = partial_right[r];
			if (start_type[r] < 0 || start_type[r] > 3) {
				throw IOException("%s: gpl-boundary's prodigal returned start_type %d (expected 0..3)", kFunction,
				                  start_type[r]);
			}
			out.data[COL_START_TYPE].SetValue(r, Value(kStartTypes[start_type[r]]));
			o_length[r] = ci.length;
			o_contig_gc[r] = ci.gc_content;
			o_table[r] = ci.translation_table;
			o_uses_sd[r] = ci.uses_sd;
			if (ci.has_meta_bin) {
				o_meta_bin[r] = ci.meta_bin;
			} else {
				FlatVector::SetNull(out.data[COL_META_BIN], r, true);
			}
			if (ci.has_model_desc) {
				out.data[COL_MODEL_DESC].SetValue(r, Value(ci.model_desc));
			} else {
				FlatVector::SetNull(out.data[COL_MODEL_DESC], r, true);
			}
			o_model_gc[r] = ci.model_gc_content;
			FlatVector::SetNull(out.data[COL_ERROR], r, true);
		}
		// Columns that pass through as they are: the scores, the RBS fields and the protein.
		out.data[COL_SCORE].Reference(chunk->data[G_SCORE]);
		out.data[COL_CSCORE].Reference(chunk->data[G_CSCORE]);
		out.data[COL_SSCORE].Reference(chunk->data[G_SSCORE]);
		out.data[COL_RSCORE].Reference(chunk->data[G_RSCORE]);
		out.data[COL_USCORE].Reference(chunk->data[G_USCORE]);
		out.data[COL_TSCORE].Reference(chunk->data[G_TSCORE]);
		out.data[COL_CONFIDENCE].Reference(chunk->data[G_CONFIDENCE]);
		out.data[COL_GENE_GC_CONTENT].Reference(chunk->data[G_GC_CONT]);
		out.data[COL_RBS_MOTIF].Reference(chunk->data[G_RBS_MOTIF]);
		out.data[COL_RBS_SPACER].Reference(chunk->data[G_RBS_SPACER]);
		out.data[COL_PROTEIN].Reference(chunk->data[G_PROTEIN]);
		ConstantVector::Reference(out.data[COL_STATUS], Value(kStatusOk), count_t(n));
		out.SetChildCardinality(n);
		// The genes in order, with a row for each contig without genes between them.
		idx_t from = 0;
		for (idx_t r = 0; r < n; r++) {
			const idx_t contig = static_cast<idx_t>(seq_index[r] - 1);
			if (contig + 1 < next_contig) {
				throw IOException("%s: gpl-boundary's prodigal returned genes out of contig order. %s", kFunction,
				                  kOutputMismatch);
			}
			if (contig > next_contig) {
				if (r > from) {
					DataChunk part;
					part.InitializeEmpty(out.GetTypes());
					part.Slice(out, from, r);
					g.rows->Append(part);
					from = r;
				}
				append_barren_contigs(contig);
			}
			next_contig = contig + 1;
		}
		if (from == 0) {
			g.rows->Append(out);
		} else {
			DataChunk part;
			part.InitializeEmpty(out.GetTypes());
			part.Slice(out, from, n);
			g.rows->Append(part);
		}
		n_genes += n;
	}
	append_barren_contigs(info.size());
	return n_genes;
}

std::string PredictConfig(const PredictGenesBindData &data, const Model *model) {
	bt2_daemon::ConfigJsonBuilder cfg;
	cfg.append_bool("meta_mode", data.meta_mode);
	if (model) {
		cfg.append_str("training_path", model->path);
		cfg.append_str("training_sha256", model->sha256);
		cfg.append_bool("closed_ends", model->closed_ends);
		cfg.append_bool("mask_regions", model->mask_regions);
	} else {
		cfg.append_bool("closed_ends", data.closed_ends);
		cfg.append_bool("mask_regions", data.mask_regions);
	}
	return cfg.build();
}

// Call genes on a genome's contigs in batches of up to kBatchBases. Returns the number of genes.
idx_t PredictInBatches(ClientContext &context, PredictGenesGlobalState &g, const PredictGenesBindData &data,
                       const std::vector<Contig> &contigs, const Model *model) {
	const auto config = PredictConfig(data, model);
	idx_t n_genes = 0;
	idx_t begin = 0;
	while (begin < contigs.size()) {
		idx_t end = begin;
		idx_t bases = 0;
		while (end < contigs.size() && (end == begin || bases + contigs[end].sequence.size() <= kBatchBases)) {
			bases += contigs[end].sequence.size();
			end++;
		}
		n_genes += Predict(context, g, config, contigs, begin, end);
		begin = end;
	}
	return n_genes;
}

// The model prodigal-train reports in its result JSON.
Model ParseTrainResult(const std::string &json) {
	auto *doc = yj::yyjson_read(json.c_str(), json.size(), 0);
	if (!doc) {
		throw IOException("%s: prodigal-train's result is not JSON: %s", kFunction, json);
	}
	Model m;
	auto *root = yj::yyjson_doc_get_root(doc);
	auto *path = yj::yyjson_obj_get(root, "training_path");
	auto *sha = yj::yyjson_obj_get(root, "training_sha256");
	auto *closed = yj::yyjson_obj_get(root, "closed_ends");
	auto *mask = yj::yyjson_obj_get(root, "mask_regions");
	const bool ok =
	    yj::yyjson_is_str(path) && yj::yyjson_is_str(sha) && yj::yyjson_is_bool(closed) && yj::yyjson_is_bool(mask);
	if (ok) {
		m.path = yj::yyjson_get_str(path);
		m.sha256 = yj::yyjson_get_str(sha);
		m.closed_ends = yj::yyjson_get_bool(closed);
		m.mask_regions = yj::yyjson_get_bool(mask);
	}
	yj::yyjson_doc_free(doc);
	if (!ok) {
		throw IOException("%s: prodigal-train's result lacks training_path, training_sha256, closed_ends or "
		                  "mask_regions: %s. %s",
		                  kFunction, json, kOutputMismatch);
	}
	return m;
}

// Single-genome mode: the next genome, whole. Trains its model unless one is fixed, then calls its genes.
bool ProduceSingleGenome(ClientContext &context, PredictGenesGlobalState &g, const PredictGenesBindData &data) {
	Contig first;
	if (!NextContig(g, data, first)) {
		return false;
	}
	const Value genome = first.genome;
	std::vector<Contig> contigs;
	contigs.push_back(std::move(first));
	while (NextIsSameGenome(g, data, genome)) {
		Contig c;
		NextContig(g, data, c);
		contigs.push_back(std::move(c));
	}

	Model model;
	if (g.have_fixed_model) {
		model = g.fixed_model;
	} else {
		bt2_daemon::ConfigJsonBuilder cfg;
		cfg.append_str("training_path", g.train_path);
		if (data.trans_table_set) {
			cfg.append_int("trans_table", data.trans_table);
		}
		cfg.append_bool("closed_ends", data.closed_ends);
		cfg.append_bool("mask_regions", data.mask_regions);
		gb::SubmitResult trained;
		try {
			trained = Submit(g, kTrainTool, cfg.build(), EncodeContigs(contigs, 0, contigs.size()));
		} catch (const gb::BatchRefused &e) {
			// Only a refusal: the daemon dying or breaking the protocol fails the query. Input and training_path's
			// directory are checked before (see InitGlobal), so this is Prodigal declining to train on this genome
			// (too little sequence, for one). The caller can rerun it in metagenomic mode.
			AppendStatusRow(context, g, genome, kStatusTrainingFailed, e.reason());
			return true;
		}
		model = ParseTrainResult(trained.result_json);
		if (g.write_user_model) {
			// The user's model now exists; it is this call's fixed model from here on.
			g.fixed_model = model;
			g.have_fixed_model = true;
		}
	}
	if (PredictInBatches(context, g, data, contigs, &model) == 0) {
		AppendStatusRow(context, g, genome, kStatusNoGenes, std::string());
	}
	return true;
}

// Metagenomic mode: the next batch of up to kBatchBases of the current genome.
bool ProduceMetaBatch(ClientContext &context, PredictGenesGlobalState &g, const PredictGenesBindData &data) {
	Contig first;
	if (!NextContig(g, data, first)) {
		return false;
	}
	if (!g.in_genome || !(first.genome == g.current_genome)) {
		g.in_genome = true;
		g.current_genome = first.genome;
		g.current_genome_genes = 0;
	}
	const Value genome = first.genome;
	std::vector<Contig> batch;
	idx_t bases = first.sequence.size();
	batch.push_back(std::move(first));
	bool genome_continues = false;
	while ((genome_continues = NextIsSameGenome(g, data, genome)) && bases + g.pending.sequence.size() <= kBatchBases) {
		Contig c;
		NextContig(g, data, c);
		bases += c.sequence.size();
		batch.push_back(std::move(c));
	}
	g.current_genome_genes += Predict(context, g, PredictConfig(data, nullptr), batch, 0, batch.size());
	if (!genome_continues) {
		if (g.current_genome_genes == 0) {
			AppendStatusRow(context, g, genome, kStatusNoGenes, std::string());
		}
		g.in_genome = false;
	}
	return true;
}

// The contigs, (genome, read_id, sequence1[, sequence2]), one genome after another, each's contigs in order: order_by
// order (the file's, for read_fastx) as the CLI reads them, or for chunked input the keys' order, each contig's
// chunks joined in chunk_index order.
std::string ContigsQuery(const PredictGenesBindData &data) {
	const auto relation = SQLIdentifier::ToString(data.relation);
	if (!data.chunked) {
		const auto genome = SQLIdentifier::ToString(data.genome_column);
		std::string sql = "SELECT " + genome + ", read_id, sequence1";
		if (data.has_sequence2) {
			sql += ", sequence2";
		}
		return sql + " FROM " + relation + " WHERE " + genome + " IS NOT NULL ORDER BY " + genome + ", " +
		       SQLIdentifier::ToString(data.order_column) + ", read_id";
	}
	const auto key = SQLIdentifier::ToString(data.key_column);
	// A BLOB's bytes, as they are: a ::VARCHAR cast would write a byte DuckDB escapes as the four characters '\xHH'.
	const std::string chunk_text =
	    data.chunk_data_type.id() == LogicalTypeId::BLOB ? "decode(c.chunk_data)" : "c.chunk_data";
	return StringUtil::Format(
	    "SELECT m.genome_id, c.%s AS read_id, "
	    "CASE WHEN count(c.chunk_data) <> count(*) "
	    "     THEN error(printf('%s: contig ''%%s'' has a NULL chunk_data', c.%s::VARCHAR)) "
	    "     WHEN min(c.chunk_index) <> 0 OR max(c.chunk_index) + 1 <> count(*) "
	    "       OR count(DISTINCT c.chunk_index) <> count(*) "
	    "     THEN error(printf('%s: the chunk_index values of contig ''%%s'' are not 0, 1, 2, ... with none missing "
	    "or repeated', c.%s::VARCHAR)) "
	    "     ELSE string_agg(%s, '' ORDER BY c.chunk_index) END AS sequence1 "
	    "FROM %s c JOIN %s m ON c.%s = m.%s WHERE m.genome_id IS NOT NULL "
	    "GROUP BY m.genome_id, c.%s ORDER BY m.genome_id, c.%s",
	    key, kFunction, key, kFunction, key, chunk_text, relation, SQLIdentifier::ToString(data.mapping_table), key,
	    key, key, key);
}

// The number of genomes ContigsQuery returns.
std::string GenomeCountQuery(const PredictGenesBindData &data) {
	const auto relation = SQLIdentifier::ToString(data.relation);
	if (!data.chunked) {
		return "SELECT count(DISTINCT " + SQLIdentifier::ToString(data.genome_column) + ") FROM " + relation;
	}
	const auto key = SQLIdentifier::ToString(data.key_column);
	return "SELECT count(DISTINCT genome_id) FROM " + SQLIdentifier::ToString(data.mapping_table) + " WHERE " + key +
	       " IN (SELECT " + key + " FROM " + relation + ")";
}

// A key mapped twice would put its contig in two genomes, or twice in one. Only the relation's keys, mapped to a
// genome, are checked: ContigsQuery reads no other.
void RequireOneGenomePerKey(Connection &conn, const PredictGenesBindData &data) {
	const auto key = SQLIdentifier::ToString(data.key_column);
	auto result =
	    conn.Query("SELECT " + key + "::VARCHAR FROM " + SQLIdentifier::ToString(data.mapping_table) +
	               " WHERE genome_id IS NOT NULL AND " + key + " IN (SELECT " + key + " FROM " +
	               SQLIdentifier::ToString(data.relation) + ") GROUP BY " + key + " HAVING count(*) > 1 LIMIT 1");
	if (result->HasError()) {
		throw InvalidInputException("%s: cannot read mapping_table '%s': %s", kFunction, data.mapping_table,
		                            result->GetError());
	}
	auto chunk = result->Fetch();
	if (chunk && chunk->size() > 0) {
		throw InvalidInputException("%s: mapping_table '%s' maps %s '%s' more than once; each contig belongs to one "
		                            "genome",
		                            kFunction, data.mapping_table, data.key_column, chunk->GetValue(0, 0).ToString());
	}
}

unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	const auto &data = input.bind_data->Cast<PredictGenesBindData>();
	auto g = make_uniq<PredictGenesGlobalState>();

	const std::string gpl_path = FindGplBoundaryOrThrow(kFunction);
	RequireGplBoundaryAtLeast(gpl_path, kFunction, kMinRelease,
	                          "Prodigal's prodigal-train, proteins and per-contig model data");
	RequireProdigalVersion(gpl_path);

	auto &fs = FileSystem::GetFileSystem(context);
	auto conn = make_uniq<Connection>(MakeReadOnlyHelperConnection(context));
	if (data.chunked) {
		RequireOneGenomePerKey(*conn, data);
	}
	if (!data.meta_mode) {
		if (!data.training_path.empty() && fs.FileExists(data.training_path)) {
			if (data.trans_table_set) {
				throw InvalidInputException("%s: trans_table cannot be set with an existing model at training_path "
				                            "'%s': the model holds its translation table",
				                            kFunction, data.training_path);
			}
			g->fixed_model.path = data.training_path;
			g->fixed_model.sha256 = HashFile(*conn, data.training_path);
			// A model does not store how it was trained on edges and N runs; it is called as the user says.
			g->fixed_model.closed_ends = data.closed_ends;
			g->fixed_model.mask_regions = data.mask_regions;
			g->have_fixed_model = true;
		} else if (!data.training_path.empty()) {
			auto count = conn->Query(GenomeCountQuery(data));
			if (count->HasError()) {
				throw InvalidInputException("%s: cannot read '%s': %s", kFunction, data.relation, count->GetError());
			}
			auto chunk = count->Fetch();
			const auto n_genomes = chunk ? chunk->GetValue(0, 0).GetValue<int64_t>() : 0;
			if (n_genomes > 1) {
				throw InvalidInputException(
				    "%s: training_path '%s' does not exist, so a model would be trained and written there, which "
				    "takes exactly one genome; '%s' has %lld. Train on one genome first, or omit training_path to "
				    "train a model per genome",
				    kFunction, data.training_path, data.relation, static_cast<long long>(n_genomes));
			}
			// Checked here, so that prodigal-train refusing to write it is not taken for the genome failing to train.
			const auto parent = std::filesystem::path(data.training_path).parent_path();
			const auto dir = parent.empty() ? std::string(".") : parent.string();
			if (!fs.DirectoryExists(dir)) {
				throw InvalidInputException("%s: training_path '%s' does not exist, so a model would be written there, "
				                            "but its directory '%s' does not exist",
				                            kFunction, data.training_path, dir);
			}
			if (::access(dir.c_str(), W_OK) != 0) {
				throw InvalidInputException("%s: training_path '%s' does not exist, so a model would be written there, "
				                            "but its directory '%s' is not writable",
				                            kFunction, data.training_path, dir);
			}
			g->train_path = data.training_path;
			g->write_user_model = true;
		} else {
			g->temp_dir = MakeTempModelDir();
			g->train_path = g->temp_dir + "/genome.trn";
		}
		if (data.trans_table_set && !g->have_fixed_model) {
			RequireTransTable(gpl_path, data.trans_table);
		}
	}

	std::vector<std::string> argv = {gpl_path};
	g->session = std::make_unique<gb::Session>(gb::ChildProcess(argv));
	try {
		g->session->Initialize();
	} catch (const std::runtime_error &e) {
		throw IOException("%s: %s", kFunction, e.what());
	}
	const auto check_tool = [&](const char *tool, uint32_t expected) {
		const auto got = g->session->tool_schema_version(tool);
		if (got != expected) {
			throw IOException("%s: gpl-boundary's %s has schema_version %u; miint is built for %u. %s", kFunction, tool,
			                  got, expected, kOutputMismatch);
		}
	};
	check_tool(kPredictTool, kProdigalSchemaVersion);
	if (!data.meta_mode && !g->have_fixed_model) {
		check_tool(kTrainTool, kProdigalTrainSchemaVersion);
	}

	g->conn = std::move(conn);
	g->stream = SubmitStream(*g->conn, ContigsQuery(data));
	if (g->stream->HasError()) {
		throw InvalidInputException("%s: cannot read '%s': %s", kFunction, data.relation, g->stream->GetError());
	}

	g->rows = make_uniq<ColumnDataCollection>(context, data.return_types);
	return std::move(g);
}

void Execute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	const auto &data = input.bind_data->Cast<PredictGenesBindData>();
	auto &g = input.global_state->Cast<PredictGenesGlobalState>();
	while (true) {
		if (g.scanning) {
			g.rows->Scan(g.scan, output);
			if (output.size() > 0) {
				return;
			}
			g.scanning = false;
		}
		context.BlockingInterruptCheck();
		g.rows->Reset();
		const bool produced =
		    data.meta_mode ? ProduceMetaBatch(context, g, data) : ProduceSingleGenome(context, g, data);
		if (!produced) {
			output.SetChildCardinality(0);
			return;
		}
		g.rows->InitializeScan(g.scan);
		g.scanning = true;
	}
}

} // namespace

TableFunction PredictGenesTableFunction::GetFunction() {
	TableFunction fn("predict_genes", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal);
	AddNamedParameter(fn, "genome_id", LogicalType::VARCHAR);
	AddNamedParameter(fn, "mapping_table", LogicalType::VARCHAR);
	AddNamedParameter(fn, "method", LogicalType::VARCHAR);
	AddNamedParameter(fn, "order_by", LogicalType::VARCHAR);
	AddNamedParameter(fn, "meta_mode", LogicalType::BOOLEAN);
	AddNamedParameter(fn, "trans_table", LogicalType::INTEGER);
	AddNamedParameter(fn, "closed_ends", LogicalType::BOOLEAN);
	AddNamedParameter(fn, "mask_regions", LogicalType::BOOLEAN);
	AddNamedParameter(fn, "training_path", LogicalType::VARCHAR);
	return fn;
}

void PredictGenesTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
