#include "align_mafft.hpp"
#include "MafftAligner.hpp"
#include "per_sample_table_function.hpp"
#include "sequence_table_reader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/parallel/task_scheduler.hpp"
#include <unordered_set>

namespace duckdb {

struct AlignMafftData : public TableFunctionData {
	std::string table_name;

	bool has_sample_id = false;
	PerSampleBindInfo sample_info;
};

struct AlignMafftGlobalState : public PerSampleGlobalState {
	// Non-sample path only: the full MAFFT run is done once at InitGlobal and drained here.
	// MAFFT needs every sequence in memory simultaneously before producing any output.
	std::vector<std::string> names;
	std::vector<std::string> sequences;
	std::vector<int32_t> original_lengths;
	int32_t aligned_length = 0;

	// Forwarded to MAFFT's internal pthread pool (cfg.n_threads). Distinct
	// from PerSampleGlobalState::max_threads, which controls duckdb-side
	// concurrency over the process-wide MAFFT mutex.
	int mafft_threads = 1;
};

struct AlignMafftLocalState : public LocalTableFunctionState {
	unique_ptr<Connection> conn; // sample_id path only

	// Current buffer (non-sample path: drains gstate once via shared data; sample path:
	// refilled per claimed sample).
	std::vector<std::string> names;
	std::vector<std::string> sequences;
	std::vector<int32_t> original_lengths;
	int32_t aligned_length = 0;
	idx_t current_row = 0;
	Value sample_value;
};

// Validate align_mafft's preconditions and run MAFFT on the loaded set. Writes the
// aligned result into the supplied vectors. `sample_literal` is empty in the non-sample
// path and non-empty in the per-sample path; error messages include the sample suffix
// only when present so single-sample users get clean diagnostics.
static void ValidateAndAlignInto(const std::string &sample_literal, LoadedSingleEndSequences &loaded,
                                 std::vector<std::string> &out_names, std::vector<std::string> &out_sequences,
                                 std::vector<int32_t> &out_original_lengths, int32_t &out_aligned_length,
                                 int n_threads) {
	bool is_sample = !sample_literal.empty();
	auto sample_ctx = is_sample ? (" in sample " + sample_literal) : std::string();

	if (loaded.sequences.size() < 2) {
		throw InvalidInputException("align_mafft requires at least 2 sequences%s, got %d", sample_ctx,
		                            static_cast<int>(loaded.sequences.size()));
	}

	for (size_t i = 0; i < loaded.sequences.size(); i++) {
		if (loaded.sequences[i].size() < 6) {
			throw InvalidInputException("align_mafft: sequence for read_id '%s'%s is too short (%d chars, minimum 6)",
			                            loaded.labels[i], sample_ctx, static_cast<int>(loaded.sequences[i].size()));
		}
	}

	std::unordered_set<std::string> seen;
	seen.reserve(loaded.labels.size());
	for (const auto &id : loaded.labels) {
		if (!seen.insert(id).second) {
			throw InvalidInputException("align_mafft: duplicate read_id '%s'%s in input table", id, sample_ctx);
		}
	}

	std::vector<std::string> comments(loaded.labels.size(), "");
	miint::MafftAligner aligner;
	auto result = aligner.align(loaded.labels, comments, loaded.sequences, n_threads);

	out_names = std::move(result.names);
	out_sequences = std::move(result.sequences);
	out_aligned_length = static_cast<int32_t>(result.aligned_length);
	out_original_lengths.clear();
	out_original_lengths.reserve(result.original_lengths.size());
	for (auto len : result.original_lengths) {
		out_original_lengths.push_back(static_cast<int32_t>(len));
	}
}

static unique_ptr<FunctionData> AlignMafftBind(ClientContext &context, TableFunctionBindInput &input,
                                               vector<LogicalType> &return_types, vector<string> &names) {
	if (input.inputs.empty() || input.inputs[0].IsNull()) {
		throw BinderException("align_mafft requires a sequence table name argument");
	}
	auto table_name = input.inputs[0].ToString();

	auto schema = ValidateSequenceTableSchema(context, table_name);
	if (schema.has_sequence2) {
		throw BinderException("align_mafft does not support paired-end data");
	}

	auto data = make_uniq<AlignMafftData>();
	data->table_name = std::move(table_name);

	auto sample_it = input.named_parameters.find("sample_id");
	if (sample_it != input.named_parameters.end()) {
		data->has_sample_id = true;
		data->sample_info.sample_id_col = sample_it->second.GetValue<string>();
	}

	if (data->has_sample_id) {
		auto &db = DatabaseInstance::GetDatabase(context);
		Connection conn(db);
		DiscoverSamples(conn, data->table_name, data->sample_info.sample_id_col,
		                {"sequence_index", "read_id", "aligned_sequence", "original_length", "aligned_length"},
		                "align_mafft", data->sample_info);
		names.push_back(data->sample_info.sample_id_col);
		return_types.push_back(data->sample_info.sample_id_type);
	}

	names.emplace_back("sequence_index");
	names.emplace_back("read_id");
	names.emplace_back("aligned_sequence");
	names.emplace_back("original_length");
	names.emplace_back("aligned_length");
	return_types.emplace_back(LogicalType::BIGINT);
	return_types.emplace_back(LogicalType::VARCHAR);
	return_types.emplace_back(LogicalType::VARCHAR);
	return_types.emplace_back(LogicalType::INTEGER);
	return_types.emplace_back(LogicalType::INTEGER);

	return data;
}

static unique_ptr<GlobalTableFunctionState> AlignMafftInitGlobal(ClientContext &context,
                                                                 TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<AlignMafftData>();
	auto gstate = make_uniq<AlignMafftGlobalState>();

	// Hand the duckdb-configured thread budget to MAFFT's internal pool.
	// Independent of max_threads (which gates duckdb-side concurrency over
	// the MAFFT mutex); MAFFT's own pthreads run inside a single align call.
	auto db_threads = TaskScheduler::GetScheduler(context).NumberOfThreads();
	gstate->mafft_threads = db_threads > 0 ? NumericCast<int>(db_threads) : 1;

	if (data.has_sample_id) {
		// MAFFT holds a process-wide mutex — running multiple duckdb threads
		// would only queue them behind the same lock. Force serial to avoid
		// the illusion of parallelism. (MAFFT itself still threads internally
		// via gstate->mafft_threads.)
		InitPerSampleGlobal(context, *gstate, data.sample_info.sample_values.size(), /*max_threads_hint=*/1);
		return gstate;
	}

	gstate->max_threads = 1;

	auto loaded = LoadSingleEndSequences(context, data.table_name, "align_mafft", /*strict=*/true);
	ValidateAndAlignInto(/*sample_literal=*/"", loaded, gstate->names, gstate->sequences, gstate->original_lengths,
	                     gstate->aligned_length, gstate->mafft_threads);

	return gstate;
}

static unique_ptr<LocalTableFunctionState> AlignMafftInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                               GlobalTableFunctionState * /*global_state*/) {
	auto &data = input.bind_data->Cast<AlignMafftData>();
	auto lstate = make_uniq<AlignMafftLocalState>();
	if (data.has_sample_id) {
		auto &db = DatabaseInstance::GetDatabase(context.client);
		lstate->conn = make_uniq<Connection>(db);
	}
	return lstate;
}

// Emit up to STANDARD_VECTOR_SIZE rows from the supplied align-result slice. When
// has_sample_id, col 0 is the sample column, filled via ConstantVector reference.
static void EmitRows(const AlignMafftData &data, AlignMafftLocalState &lstate, const std::vector<std::string> &names,
                     const std::vector<std::string> &sequences, const std::vector<int32_t> &original_lengths,
                     int32_t aligned_length, DataChunk &output) {
	idx_t total = names.size();
	idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - lstate.current_row);

	idx_t col = 0;
	if (data.has_sample_id) {
		output.data[col++].Reference(lstate.sample_value);
	}
	auto &sequence_index_vec = output.data[col++];
	auto &read_id_vec = output.data[col++];
	auto &aligned_seq_vec = output.data[col++];
	auto &original_length_vec = output.data[col++];
	auto &aligned_length_vec = output.data[col++];
	auto sequence_index_data = FlatVector::GetData<int64_t>(sequence_index_vec);
	auto read_id_data = FlatVector::GetData<string_t>(read_id_vec);
	auto aligned_seq_data = FlatVector::GetData<string_t>(aligned_seq_vec);
	auto original_length_data = FlatVector::GetData<int32_t>(original_length_vec);
	auto aligned_length_data = FlatVector::GetData<int32_t>(aligned_length_vec);
	for (idx_t i = 0; i < count; i++) {
		idx_t row = lstate.current_row + i;
		sequence_index_data[i] = static_cast<int64_t>(row);
		read_id_data[i] = StringVector::AddString(read_id_vec, names[row]);
		aligned_seq_data[i] = StringVector::AddString(aligned_seq_vec, sequences[row]);
		original_length_data[i] = original_lengths[row];
		aligned_length_data[i] = aligned_length;
	}

	lstate.current_row += count;
	output.SetCardinality(count);
}

static void AlignMafftExecute(ClientContext & /*context*/, TableFunctionInput &data_p, DataChunk &output) {
	auto &data = data_p.bind_data->Cast<AlignMafftData>();
	auto &gstate = data_p.global_state->Cast<AlignMafftGlobalState>();
	auto &lstate = data_p.local_state->Cast<AlignMafftLocalState>();

	if (!data.has_sample_id) {
		if (lstate.current_row >= gstate.names.size()) {
			output.SetCardinality(0);
			return;
		}
		EmitRows(data, lstate, gstate.names, gstate.sequences, gstate.original_lengths, gstate.aligned_length, output);
		return;
	}

	while (true) {
		if (lstate.current_row < lstate.names.size()) {
			EmitRows(data, lstate, lstate.names, lstate.sequences, lstate.original_lengths, lstate.aligned_length,
			         output);
			return;
		}
		idx_t sample_idx;
		if (!ClaimNextSample(gstate, data.sample_info.sample_values.size(), sample_idx)) {
			output.SetCardinality(0);
			return;
		}
		lstate.sample_value = data.sample_info.sample_values[sample_idx];
		auto sample_literal = lstate.sample_value.ToSQLString();
		auto q_col = KeywordHelper::WriteOptionallyQuoted(data.sample_info.sample_id_col);
		// Same CAST-as-VARCHAR equality as the other per-sample call sites; see the note
		// in deblur_table_function.cpp for the DECIMAL caveat.
		auto where_sql = "CAST(" + q_col + " AS VARCHAR) = CAST(" + sample_literal + " AS VARCHAR)";
		auto loaded = LoadSingleEndSequences(*lstate.conn, data.table_name, "align_mafft", /*strict=*/true, where_sql);
		ValidateAndAlignInto(sample_literal, loaded, lstate.names, lstate.sequences, lstate.original_lengths,
		                     lstate.aligned_length, gstate.mafft_threads);
		lstate.current_row = 0;
	}
}

TableFunction AlignMafftTableFunction::GetFunction() {
	auto tf = TableFunction("align_mafft", {LogicalType::VARCHAR}, AlignMafftExecute, AlignMafftBind,
	                        AlignMafftInitGlobal, AlignMafftInitLocal);
	tf.named_parameters["sample_id"] = LogicalType::VARCHAR;
	// Match sibling aligners — allow downstream CTAS pipelines to parallelize.
	// Callers that want deterministic order should ORDER BY (sample_id,) sequence_index.
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void AlignMafftTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
