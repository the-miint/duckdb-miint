#include "woltka_batched.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/parser/keyword_helper.hpp"
#include "duckdb/parser/qualified_name.hpp"
#include "duckdb/parser/parsed_data/create_table_function_info.hpp"

namespace duckdb {

unique_ptr<FunctionData> WoltkaOguPerSampleBatchedTableFunction::Bind(ClientContext &context,
                                                                       TableFunctionBindInput &input,
                                                                       vector<LogicalType> &return_types,
                                                                       vector<std::string> &names) {
	// Validate argument count
	if (input.inputs.size() != 3) {
		throw InvalidInputException(
		    "woltka_ogu_per_sample_batched requires 3 arguments: relation, sample_id_field, sequence_id_field");
	}

	// Parse relation name (first argument - can be table name or expression)
	if (input.inputs[0].type().id() != LogicalTypeId::VARCHAR) {
		throw InvalidInputException("woltka_ogu_per_sample_batched: first argument (relation) must be VARCHAR");
	}
	std::string relation_name = input.inputs[0].ToString();

	// Parse sample_id_field (second argument)
	if (input.inputs[1].type().id() != LogicalTypeId::VARCHAR) {
		throw InvalidInputException("woltka_ogu_per_sample_batched: second argument (sample_id_field) must be VARCHAR");
	}
	std::string sample_id_field = input.inputs[1].ToString();

	// Parse sequence_id_field (third argument)
	if (input.inputs[2].type().id() != LogicalTypeId::VARCHAR) {
		throw InvalidInputException(
		    "woltka_ogu_per_sample_batched: third argument (sequence_id_field) must be VARCHAR");
	}
	std::string sequence_id_field = input.inputs[2].ToString();

	// Parse batch_size parameter (required for now)
	auto batch_param = input.named_parameters.find("batch_size");
	if (batch_param == input.named_parameters.end() || batch_param->second.IsNull()) {
		throw InvalidInputException("woltka_ogu_per_sample_batched: batch_size parameter is required");
	}

	if (batch_param->second.type().id() != LogicalTypeId::BIGINT &&
	    batch_param->second.type().id() != LogicalTypeId::INTEGER) {
		throw InvalidInputException("woltka_ogu_per_sample_batched: batch_size must be an integer");
	}

	int64_t batch_size = batch_param->second.GetValue<int64_t>();
	if (batch_size <= 0) {
		throw InvalidInputException("woltka_ogu_per_sample_batched: batch_size must be a positive integer");
	}

	auto data = make_uniq<Data>(relation_name, sample_id_field, sequence_id_field, static_cast<idx_t>(batch_size));

	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState>
WoltkaOguPerSampleBatchedTableFunction::InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	// Create persistent connection for streaming batch results
	// LIFETIME SAFETY: Connection internally holds a reference to the DatabaseInstance,
	// ensuring the database remains valid for the lifetime of this Connection.
	// DuckDB guarantees that ClientContext (and thus context.db) outlives GlobalState,
	// so this Connection will remain valid throughout the table function's execution.
	auto &instance = *context.db;
	gstate->connection = make_uniq<Connection>(instance);

	// Safely quote identifiers to prevent SQL injection
	auto safe_relation = QualifiedName::Parse(data.relation_name).ToString();
	auto safe_sample_field = KeywordHelper::WriteOptionallyQuoted(data.sample_id_field);

	// Query to get all distinct samples
	std::string sample_query = StringUtil::Format("SELECT DISTINCT %s FROM %s", safe_sample_field, safe_relation);

	// Execute query to get samples
	auto result = gstate->connection->Query(sample_query);

	if (result->HasError()) {
		throw InvalidInputException("Failed to query samples from relation '%s': %s", data.relation_name,
		                            result->GetError());
	}

	// Extract sample IDs from result
	while (true) {
		auto chunk = result->Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}

		auto &sample_vector = chunk->data[0];
		auto &validity = FlatVector::Validity(sample_vector);
		auto sample_data = FlatVector::GetData<string_t>(sample_vector);

		for (idx_t i = 0; i < chunk->size(); i++) {
			if (!validity.RowIsValid(i)) {
				throw InvalidInputException("woltka_ogu_per_sample_batched: sample_id cannot be NULL");
			}
			gstate->all_samples.push_back(sample_data[i].GetString());
		}
	}

	if (gstate->all_samples.empty()) {
		throw InvalidInputException("No samples found in relation '%s'", data.relation_name);
	}

	return gstate;
}

void WoltkaOguPerSampleBatchedTableFunction::ProcessChunkToOutput(const unique_ptr<DataChunk> &chunk,
                                                                    DataChunk &output) {
	// Get input vectors
	auto &sample_vec = chunk->data[0];
	auto &feature_vec = chunk->data[1];
	auto &value_vec = chunk->data[2];

	// Get source data pointers
	auto sample_src = FlatVector::GetData<string_t>(sample_vec);
	auto feature_src = FlatVector::GetData<string_t>(feature_vec);
	auto value_src = FlatVector::GetData<double>(value_vec);

	// Get output vectors
	auto &sample_out = output.data[0];
	auto &feature_out = output.data[1];
	auto &value_out = output.data[2];

	// Get destination data pointers
	auto sample_dst = FlatVector::GetData<string_t>(sample_out);
	auto feature_dst = FlatVector::GetData<string_t>(feature_out);
	auto value_dst = FlatVector::GetData<double>(value_out);

	// Copy data - no NULL validation needed as query structure guarantees non-NULL:
	// - sample_id validated non-NULL in InitGlobal, filtered by IN clause
	// - feature_id (reference) always has value per SAM spec (may be '*' but not NULL)
	// - value is SUM(1.0/COUNT(*)) which can never be NULL
	for (idx_t i = 0; i < chunk->size(); i++) {
		sample_dst[i] = StringVector::AddString(sample_out, sample_src[i]);
		feature_dst[i] = StringVector::AddString(feature_out, feature_src[i]);
		value_dst[i] = value_src[i];
	}

	output.SetCardinality(chunk->size());
}

bool WoltkaOguPerSampleBatchedTableFunction::TryStartNextBatch(GlobalState &gstate, const Data &data,
                                                                 DataChunk &output) {
	// Check if we've exhausted all batches
	if (gstate.current_batch_idx >= gstate.all_samples.size()) {
		return false;
	}

	idx_t batch_start = gstate.current_batch_idx;
	idx_t batch_end = std::min(batch_start + data.batch_size, static_cast<idx_t>(gstate.all_samples.size()));

	// Build sample list for IN clause (values, not identifiers - use single-quote escaping)
	// SQL INJECTION SAFETY: Sample IDs are data values, not identifiers, so they're embedded
	// as string literals in the IN clause. DuckDB follows standard SQL string literal rules:
	// single quotes are escaped by doubling them ('foo''s' represents the value "foo's").
	// This escaping is complete - no other characters need escaping in DuckDB string literals.
	// Identifiers (table/column names) are protected separately via KeywordHelper and
	// QualifiedName::Parse below, which handle quoting per DuckDB identifier rules.
	std::string sample_in_clause = "(";
	for (idx_t i = batch_start; i < batch_end; i++) {
		if (i > batch_start) {
			sample_in_clause += ", ";
		}
		// Escape single quotes in sample ID if needed (these are string literals, not identifiers)
		const std::string &sample = gstate.all_samples[i];
		if (sample.find('\'') != std::string::npos) {
			// Sample contains single quote - needs escaping
			sample_in_clause += "'" + StringUtil::Replace(sample, "'", "''") + "'";
		} else {
			// No escaping needed - avoid allocation
			sample_in_clause += "'" + sample + "'";
		}
	}
	sample_in_clause += ")";

	// Safely quote identifiers to prevent SQL injection
	auto safe_relation = QualifiedName::Parse(data.relation_name).ToString();
	auto safe_sequence_field = KeywordHelper::WriteOptionallyQuoted(data.sequence_id_field);
	auto safe_sample_field = KeywordHelper::WriteOptionallyQuoted(data.sample_id_field);

	// Build batch query using woltka logic with sample filter
	std::string batch_query = StringUtil::Format(
	    R"(
        WITH base AS (
            SELECT DISTINCT
                %s AS query_local_id_field,
                %s AS query_local_sample_id,
                reference AS feature_id,
                alignment_is_read1(flags::USMALLINT) AS is_fwd
            FROM %s
            WHERE %s IN %s
        ),
        with_counts AS (
            SELECT
                query_local_sample_id,
                feature_id,
                1.0 / COUNT(*) OVER (PARTITION BY query_local_id_field, is_fwd) AS local_value
            FROM base
        )
        SELECT
            query_local_sample_id AS sample_id,
            feature_id,
            SUM(local_value) AS value
        FROM with_counts
        GROUP BY query_local_sample_id, feature_id
    )",
	    safe_sequence_field, safe_sample_field, safe_relation, safe_sample_field, sample_in_clause);

	// Execute batch query using persistent connection
	auto result = gstate.connection->Query(batch_query);

	if (result->HasError()) {
		throw InvalidInputException("Failed to execute batch query for samples %zu-%zu (%zu samples): %s", batch_start,
		                            batch_end - 1, batch_end - batch_start, result->GetError());
	}

	// Move to next batch index
	gstate.current_batch_idx = batch_end;

	// Try to fetch first chunk to check if batch has results
	auto chunk = result->Fetch();
	if (chunk && chunk->size() > 0) {
		// Store result for streaming and process this first chunk
		gstate.current_batch_result = std::move(result);
		ProcessChunkToOutput(chunk, output);
		return true;
	}

	// Batch was empty
	return false;
}

void WoltkaOguPerSampleBatchedTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p,
                                                      DataChunk &output) {
	auto &data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	// Try to stream from current batch result if available
	if (gstate.current_batch_result) {
		auto chunk = gstate.current_batch_result->Fetch();
		if (chunk && chunk->size() > 0) {
			ProcessChunkToOutput(chunk, output);
			return;
		}
		// Current batch exhausted
		gstate.current_batch_result.reset();
	}

	// Start new batches until we find one with results or exhaust all batches
	// Use a loop instead of recursion to avoid stack overflow with consecutive empty batches
	while (gstate.current_batch_idx < gstate.all_samples.size()) {
		if (TryStartNextBatch(gstate, data, output)) {
			// Found batch with data, first chunk already in output
			return;
		}
		// Batch was empty, continue to next batch
	}

	// All batches exhausted
	output.SetCardinality(0);
}

TableFunction WoltkaOguPerSampleBatchedTableFunction::GetFunction() {
	TableFunction func("woltka_ogu_per_sample_batched",
	                   {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal);

	func.named_parameters["batch_size"] = LogicalType::BIGINT;

	return func;
}

void WoltkaOguPerSampleBatchedTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
