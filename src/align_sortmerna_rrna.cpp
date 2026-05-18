#include "align_sortmerna_rrna.hpp"

#include "align_sortmerna_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

namespace duckdb {

unique_ptr<FunctionData> AlignSortMeRNARRNATableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                               vector<LogicalType> &return_types,
                                                               vector<std::string> &names) {
	auto data = make_uniq<Data>();

	if (input.inputs.empty() || input.inputs[0].IsNull()) {
		throw BinderException("align_sortmerna_rrna: query_table is required");
	}
	data->query_table = input.inputs[0].ToString();
	data->ref_paths = ParseSortMeRNARefPaths(input.named_parameters, "align_sortmerna_rrna");
	ParseSortMeRNAConfigParams(input.named_parameters, data->config);

	data->query_schema = ValidateSequenceTableSchema(context, data->query_table, /*allow_bigint=*/true);
	if (data->query_schema.has_sequence2 != data->config.paired) {
		throw BinderException("align_sortmerna_rrna: query table paired-ness (sequence2 %s) does not match "
		                      "paired=%s; set the paired parameter to match or reshape the query",
		                      data->query_schema.has_sequence2 ? "present" : "absent",
		                      data->config.paired ? "true" : "false");
	}

	// Output schema: read_id mirrors the query column type; the other 12
	// columns have fixed types (ref_name is free-form FASTA header text).
	data->types = GetSortMeRNARRNAOutputTypes(data->query_schema.id_type);

	for (const auto &n : data->names)
		names.emplace_back(n);
	for (const auto &t : data->types)
		return_types.emplace_back(t);
	return data;
}

unique_ptr<GlobalTableFunctionState> AlignSortMeRNARRNATableFunction::InitGlobal(ClientContext &context,
                                                                                 TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	miint::SortMeRNAConfig cfg = data.config;
	if (cfg.num_threads <= 0) {
		cfg.num_threads = NumericCast<int32_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	}
	gstate->aligner = std::make_unique<miint::SortMeRNAAligner>(cfg, data.ref_paths);
	gstate->query_stream = std::make_unique<QuerySequenceStream>(context, data.query_table, data.query_schema);
	// Execute holds gstate.lock across align() on the assumption that
	// MaxThreads() is 1. If that invariant ever changes, the shared
	// result_buffer's deposit path needs a rethink before this body is run
	// concurrently.
	D_ASSERT(gstate->MaxThreads() == 1);
	return gstate;
}

unique_ptr<LocalTableFunctionState>
AlignSortMeRNARRNATableFunction::InitLocal(ExecutionContext &, TableFunctionInitInput &, GlobalTableFunctionState *) {
	return make_uniq<LocalState>();
}

void AlignSortMeRNARRNATableFunction::Execute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();
	auto &bind_data = data_p.bind_data->Cast<Data>();

	// MaxThreads() == 1 means only one DuckDB thread ever reaches this body;
	// the lock is defensive against future changes. Holding it across the
	// align() call is deliberate — the streaming invariant (one sub-batch in
	// flight at a time, shared buffer consumed in order) requires it, and
	// sortmerna's g_run_mutex would serialize concurrent align() calls anyway.
	std::lock_guard<std::mutex> lock(gstate.lock);

	while (true) {
		// 1. Emit buffered rows first.
		idx_t available = gstate.result_buffer.size() - gstate.buffer_offset;
		if (available > 0) {
			idx_t count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			OutputSortMeRNARRNABatch(output, gstate.result_buffer, gstate.buffer_offset, count,
			                         bind_data.query_schema.id_type);
			gstate.buffer_offset += count;
			return;
		}

		// 2. Buffer exhausted — pull next sub-batch from the query stream.
		gstate.result_buffer.clear();
		gstate.buffer_offset = 0;

		auto query_batch = gstate.query_stream->FetchSubBatch();
		if (query_batch.empty()) {
			output.SetCardinality(0);
			return;
		}

		// 3. Convert SequenceRecordBatch → SortMeRNAQueryBatch.
		miint::SortMeRNAQueryBatch queries;
		queries.read_ids = std::move(query_batch.read_ids);
		queries.sequences = std::move(query_batch.sequences1);
		if (bind_data.config.paired) {
			queries.sequences2 = std::move(query_batch.sequences2);
		}

		// 4. Align. Rethrow library errors as IOException so DuckDB's
		//    table-function error path handles them cleanly rather than
		//    letting std::runtime_error cross the ABI boundary.
		try {
			gstate.aligner->align(queries, gstate.result_buffer);
		} catch (const std::exception &e) {
			throw IOException("align_sortmerna_rrna: %s", e.what());
		}
	}
}

TableFunction AlignSortMeRNARRNATableFunction::GetFunction() {
	auto tf = TableFunction("align_sortmerna_rrna", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);
	RegisterSortMeRNANamedParameters(tf);
	// Alignment output order depends on sortmerna's internal batch traversal,
	// not on the input row order. Declare NO_ORDER so DuckDB doesn't try to
	// preserve insertion order for CTAS pipelines.
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void AlignSortMeRNARRNATableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
