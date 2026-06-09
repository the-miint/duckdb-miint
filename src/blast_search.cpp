#include "blast_search.hpp"
#include "id_column_utils.hpp"
#include "miint_log.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/database.hpp"

namespace duckdb {

static std::vector<std::string> GetBlastOutputNames() {
	return {"query_id",    "subject_id", "pct_identity",  "alignment_length", "mismatches", "gap_opens",
	        "query_start", "query_end",  "subject_start", "subject_end",      "evalue",     "bit_score"};
}

// query_id (col 0) mirrors the query table's read_id type. subject_id (col 1)
// stays VARCHAR — it is an NCBI accession, not an id drawn from a user table.
static std::vector<LogicalType> GetBlastOutputTypes(const LogicalType &query_id_type) {
	return {query_id_type,        LogicalType::VARCHAR, LogicalType::DOUBLE, LogicalType::INTEGER,
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::BIGINT, LogicalType::BIGINT,
	        LogicalType::BIGINT,  LogicalType::BIGINT,  LogicalType::DOUBLE, LogicalType::DOUBLE};
}

unique_ptr<FunctionData> BlastSearchTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                        vector<LogicalType> &return_types, vector<std::string> &names) {
	auto data = make_uniq<Data>();

	if (input.inputs[0].IsNull()) {
		throw BinderException("blast: table name cannot be NULL");
	}
	data->query_table = input.inputs[0].ToString();
	if (data->query_table.empty()) {
		throw BinderException("blast: table name cannot be empty");
	}

	auto get_string = [&](const std::string &name, std::string &out) {
		auto it = input.named_parameters.find(name);
		if (it != input.named_parameters.end() && !it->second.IsNull()) {
			out = it->second.ToString();
		}
	};
	auto get_double = [&](const std::string &name, double &out) {
		auto it = input.named_parameters.find(name);
		if (it != input.named_parameters.end() && !it->second.IsNull()) {
			out = it->second.GetValue<double>();
		}
	};
	auto get_int = [&](const std::string &name, int &out) {
		auto it = input.named_parameters.find(name);
		if (it != input.named_parameters.end() && !it->second.IsNull()) {
			out = it->second.GetValue<int>();
		}
	};
	auto get_bool = [&](const std::string &name, bool &out) {
		auto it = input.named_parameters.find(name);
		if (it != input.named_parameters.end() && !it->second.IsNull()) {
			out = it->second.GetValue<bool>();
		}
	};

	get_string("program", data->program);
	get_string("database", data->database);
	get_double("evalue", data->evalue);
	get_int("max_targets", data->max_targets);
	get_bool("megablast", data->megablast);
	get_string("api_key", data->api_key);

	if (!miint::BlastParser::ValidateProgram(data->program)) {
		throw BinderException(
		    "blast: Invalid BLAST program '%s'. Must be one of: blastn, blastp, blastx, tblastn, tblastx",
		    data->program);
	}
	if (data->evalue <= 0) {
		throw BinderException("blast: evalue must be positive (got %g)", data->evalue);
	}
	if (data->max_targets <= 0) {
		throw BinderException("blast: max_targets must be positive (got %d)", data->max_targets);
	}
	if (data->megablast && data->program != "blastn") {
		throw BinderException("blast: megablast is only valid with program='blastn' (got program='%s')", data->program);
	}

	// Capture the query read_id type (VARCHAR/BIGINT/UUID) so the output query_id
	// mirrors it; subject_id stays VARCHAR.
	auto query_schema = ValidateSequenceTableSchema(context, data->query_table, /*allow_bigint=*/true);
	data->query_id_type = query_schema.id_type;

	data->names = GetBlastOutputNames();
	data->types = GetBlastOutputTypes(data->query_id_type);
	for (auto &n : data->names) {
		names.emplace_back(n);
	}
	for (auto &t : data->types) {
		return_types.emplace_back(t);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> BlastSearchTableFunction::InitGlobal(ClientContext &context,
                                                                          TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto &db = DatabaseInstance::GetDatabase(context);
	auto gstate = make_uniq<GlobalState>();

	gstate->client = make_uniq<miint::BlastClient>(db, data.api_key);
	gstate->program = data.program;
	gstate->database = data.database;
	gstate->evalue = data.evalue;
	gstate->max_targets = data.max_targets;
	gstate->megablast = data.megablast;

	auto loaded = LoadSingleEndSequences(context, data.query_table, "blast", /*strict=*/true);
	gstate->batches = miint::BlastParser::SplitIntoBatches(loaded.labels, loaded.sequences, DEFAULT_MAX_BATCH_BYTES);

	if (gstate->batches.size() > 1) {
		miint::EmitWarning(context, "blast: splitting %d sequences into %d batches (~%zu bytes per batch limit)",
		                   static_cast<int>(loaded.labels.size()), static_cast<int>(gstate->batches.size()),
		                   DEFAULT_MAX_BATCH_BYTES);
	}

	return gstate;
}

unique_ptr<LocalTableFunctionState> BlastSearchTableFunction::InitLocal(ExecutionContext &context,
                                                                        TableFunctionInitInput &input,
                                                                        GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

bool BlastSearchTableFunction::GlobalState::FetchNextBatch(ClientContext &context) {
	while (batch_cursor < batches.size()) {
		auto &batch = batches[batch_cursor];
		batch_cursor++;

		auto fasta = miint::BlastParser::BuildFastaPayload(batch.ids, batch.sequences);
		auto submit_result = client->Submit(program, database, fasta, evalue, max_targets, megablast);
		client->WaitForCompletion(submit_result.rid, submit_result.rtoe);
		auto raw_results = client->RetrieveResults(submit_result.rid, max_targets);
		auto hits = miint::BlastParser::ParseTabularResults(raw_results);

		if (!hits.empty()) {
			result_buffer = std::move(hits);
			return true;
		}
	}
	return false;
}

static void OutputBlastHits(DataChunk &output, const std::vector<miint::BlastHit> &hits, idx_t offset, idx_t count,
                            const LogicalType &query_id_type) {
	for (idx_t i = 0; i < count; i++) {
		const auto &hit = hits[offset + i];
		// query_id mirrors the query table's id type; subject_id is always VARCHAR.
		EmitIdCell(output.data[0], i, hit.query_id, query_id_type);
		FlatVector::GetData<string_t>(output.data[1])[i] = StringVector::AddString(output.data[1], hit.subject_id);
		FlatVector::GetData<double>(output.data[2])[i] = hit.pct_identity;
		FlatVector::GetData<int32_t>(output.data[3])[i] = hit.alignment_length;
		FlatVector::GetData<int32_t>(output.data[4])[i] = hit.mismatches;
		FlatVector::GetData<int32_t>(output.data[5])[i] = hit.gap_opens;
		FlatVector::GetData<int64_t>(output.data[6])[i] = hit.query_start;
		FlatVector::GetData<int64_t>(output.data[7])[i] = hit.query_end;
		FlatVector::GetData<int64_t>(output.data[8])[i] = hit.subject_start;
		FlatVector::GetData<int64_t>(output.data[9])[i] = hit.subject_end;
		FlatVector::GetData<double>(output.data[10])[i] = hit.evalue;
		FlatVector::GetData<double>(output.data[11])[i] = hit.bit_score;
	}
	output.SetCardinality(count);
}

void BlastSearchTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();
	// Lock held during FetchNextBatch (which blocks on HTTP). Safe because
	// MaxThreads()=1 guarantees a single caller. Same pattern as read_ncbi.cpp.
	lock_guard<mutex> guard(gstate.lock);

	while (gstate.result_offset >= gstate.result_buffer.size()) {
		gstate.result_buffer.clear();
		gstate.result_offset = 0;
		if (!gstate.FetchNextBatch(context)) {
			output.SetCardinality(0);
			return;
		}
	}

	idx_t remaining = gstate.result_buffer.size() - gstate.result_offset;
	idx_t count = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);
	OutputBlastHits(output, gstate.result_buffer, gstate.result_offset, count, data.query_id_type);
	gstate.result_offset += count;
}

TableFunction BlastSearchTableFunction::GetFunction() {
	auto tf = TableFunction("blast", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["program"] = LogicalType::VARCHAR;
	tf.named_parameters["database"] = LogicalType::VARCHAR;
	tf.named_parameters["evalue"] = LogicalType::DOUBLE;
	tf.named_parameters["max_targets"] = LogicalType::INTEGER;
	tf.named_parameters["megablast"] = LogicalType::BOOLEAN;
	tf.named_parameters["api_key"] = LogicalType::VARCHAR;
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void BlastSearchTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
