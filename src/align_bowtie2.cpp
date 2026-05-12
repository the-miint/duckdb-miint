#include "align_bowtie2.hpp"
#include "align_result_utils.hpp"
#include "duckdb/common/vector_size.hpp"
#include "gpl_boundary/process.hpp"

#include <fcntl.h>
#include <mutex>
#include <sys/wait.h>
#include <unistd.h>

namespace duckdb {

// Batch size for reading queries
static constexpr idx_t QUERY_BATCH_SIZE = 1024;

unique_ptr<FunctionData> AlignBowtie2TableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                         vector<LogicalType> &return_types,
                                                         vector<std::string> &names) {
	auto data = make_uniq<Data>();

	// Required positional parameters: query_table, subject_table
	if (input.inputs.size() < 2) {
		throw BinderException("align_bowtie2 requires query_table and subject_table parameters");
	}

	data->query_table = input.inputs[0].ToString();
	data->subject_table = input.inputs[1].ToString();

	// Validate query and subject tables/views exist
	data->query_schema = ValidateSequenceTableSchema(context, data->query_table);
	ValidateSequenceTableSchema(context, data->subject_table);

	// Parse optional named parameters
	auto preset_param = input.named_parameters.find("preset");
	if (preset_param != input.named_parameters.end() && !preset_param->second.IsNull()) {
		data->config.preset = preset_param->second.ToString();
	}

	auto local_param = input.named_parameters.find("local");
	if (local_param != input.named_parameters.end() && !local_param->second.IsNull()) {
		data->config.local = local_param->second.GetValue<bool>();
	}

	auto threads_param = input.named_parameters.find("threads");
	if (threads_param != input.named_parameters.end() && !threads_param->second.IsNull()) {
		data->config.threads = threads_param->second.GetValue<int32_t>();
	}

	auto max_secondary_param = input.named_parameters.find("max_secondary");
	if (max_secondary_param != input.named_parameters.end() && !max_secondary_param->second.IsNull()) {
		data->config.max_secondary = max_secondary_param->second.GetValue<int32_t>();
	}

	auto extra_args_param = input.named_parameters.find("extra_args");
	if (extra_args_param != input.named_parameters.end() && !extra_args_param->second.IsNull()) {
		data->config.extra_args = extra_args_param->second.ToString();
	}

	auto quiet_param = input.named_parameters.find("quiet");
	if (quiet_param != input.named_parameters.end() && !quiet_param->second.IsNull()) {
		data->config.quiet = quiet_param->second.GetValue<bool>();
	}

	// Pre-load all subjects at bind time (required for indexing)
	data->subjects = ReadSubjectTable(context, data->subject_table);

	// Set output schema
	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> AlignBowtie2TableFunction::InitGlobal(ClientContext &context,
                                                                           TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	// Create aligner with config
	gstate->aligner = std::make_unique<::miint::Bowtie2Aligner>(data.config);

	// Build index from all subjects
	gstate->aligner->build_index(data.subjects);

	return gstate;
}

unique_ptr<LocalTableFunctionState> AlignBowtie2TableFunction::InitLocal(ExecutionContext &context,
                                                                         TableFunctionInitInput &input,
                                                                         GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

void AlignBowtie2TableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	std::lock_guard<std::mutex> lock(global_state.lock);

	// Check if we're done
	if (global_state.done) {
		output.SetCardinality(0);
		return;
	}

	// Determine how many results we can output this call
	idx_t available = global_state.result_buffer.size() - global_state.buffer_offset;

	// If buffer is empty or exhausted, fill it with more alignments
	while (available == 0) {
		// Clear buffer for new batch
		global_state.result_buffer.clear();
		global_state.buffer_offset = 0;

		// If we've already finished aligning, we're done
		if (global_state.finished_aligning) {
			global_state.done = true;
			output.SetCardinality(0);
			return;
		}

		// Read next batch of queries
		::miint::SequenceRecordBatch query_batch;
		bool has_more = ReadQueryBatch(context, bind_data.query_table, bind_data.query_schema, QUERY_BATCH_SIZE,
		                               global_state.current_query_offset, query_batch);

		if (query_batch.empty() && !has_more) {
			// No more queries - call finish() to get remaining results
			global_state.aligner->finish(global_state.result_buffer);
			global_state.finished_aligning = true;
		} else if (!query_batch.empty()) {
			// Align this batch
			global_state.aligner->align(query_batch, global_state.result_buffer);
		}

		available = global_state.result_buffer.size() - global_state.buffer_offset;

		// If still no results after finish(), we're done
		if (available == 0 && global_state.finished_aligning) {
			global_state.done = true;
			output.SetCardinality(0);
			return;
		}
	}

	// Output up to STANDARD_VECTOR_SIZE results
	idx_t output_count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
	idx_t offset = global_state.buffer_offset;
	auto &batch = global_state.result_buffer;

	// Set result vectors using shared utilities
	idx_t field_idx = 0;
	SetAlignResultString(output.data[field_idx++], batch.read_ids, offset, output_count);
	SetAlignResultUInt16(output.data[field_idx++], batch.flags, offset, output_count);
	SetAlignResultString(output.data[field_idx++], batch.references, offset, output_count);
	SetAlignResultInt64(output.data[field_idx++], batch.positions, offset, output_count);
	SetAlignResultInt64(output.data[field_idx++], batch.stop_positions, offset, output_count);
	SetAlignResultUInt8(output.data[field_idx++], batch.mapqs, offset, output_count);
	SetAlignResultString(output.data[field_idx++], batch.cigars, offset, output_count);
	SetAlignResultString(output.data[field_idx++], batch.mate_references, offset, output_count);
	SetAlignResultInt64(output.data[field_idx++], batch.mate_positions, offset, output_count);
	SetAlignResultInt64(output.data[field_idx++], batch.template_lengths, offset, output_count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_as_values, offset, output_count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_xs_values, offset, output_count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_ys_values, offset, output_count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_xn_values, offset, output_count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_xm_values, offset, output_count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_xo_values, offset, output_count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_xg_values, offset, output_count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_nm_values, offset, output_count);
	SetAlignResultStringNullable(output.data[field_idx++], batch.tag_yt_values, offset, output_count);
	SetAlignResultStringNullable(output.data[field_idx++], batch.tag_md_values, offset, output_count);
	SetAlignResultStringNullable(output.data[field_idx++], batch.tag_sa_values, offset, output_count);

	output.SetCardinality(output_count);
	global_state.buffer_offset += output_count;
}

TableFunction AlignBowtie2TableFunction::GetFunction() {
	auto tf = TableFunction("align_bowtie2", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal,
	                        InitLocal);

	// Named parameters matching Bowtie2Config options
	tf.named_parameters["preset"] = LogicalType::VARCHAR;        // --very-fast, --fast, --sensitive, --very-sensitive
	tf.named_parameters["local"] = LogicalType::BOOLEAN;         // --local mode
	tf.named_parameters["threads"] = LogicalType::INTEGER;       // -p parameter
	tf.named_parameters["max_secondary"] = LogicalType::INTEGER; // -k parameter
	tf.named_parameters["extra_args"] = LogicalType::VARCHAR;    // Additional arguments
	tf.named_parameters["quiet"] = LogicalType::BOOLEAN;         // Suppress stderr output (default: true)

	// Alignment output order is non-deterministic — NO_ORDER enables parallel CTAS.
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;

	return tf;
}

void AlignBowtie2TableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

// Scalar function: returns true when the gpl-boundary daemon advertises the
// `bowtie2-align` and `bowtie2-build` tools. Mirrors
// `PhylogenyFastTreeAvailableImpl` (src/phylogeny_fasttree.cpp:882-914) —
// cached once per process to keep the cost out of hot SQL paths. We probe via
// `gpl-boundary --list-tools` rather than opening a Session because the
// scalar must remain cheap (no shm setup) and stateless.
static void Bowtie2AvailableFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	(void)args;
	(void)state;
	static std::once_flag flag;
	static bool available = false;

	std::call_once(flag, []() {
		// Fully qualified `::duckdb::miint::gpl_boundary::` because the
		// existing direct-subprocess code in this file references
		// `::miint::Bowtie2Aligner` (top-level miint), and a `namespace
		// gb = ...` alias inside `namespace duckdb` would put
		// `duckdb::miint` in scope, hiding it. Phase 7 deletes the
		// direct-subprocess code and this nesting concern goes with it.
		const std::string path = ::duckdb::miint::gpl_boundary::FindGplBoundary();
		if (path.empty()) {
			available = false;
			return;
		}
		try {
			std::vector<std::string> argv = {path, "--list-tools"};
			::duckdb::miint::gpl_boundary::ChildProcess child(argv);
			// Drain stdout fully. `--list-tools` writes a small JSON array
			// (a few hundred bytes for the v0.2.0 registry) and exits — well
			// below the pipe-buffer limit, so no risk of the child blocking
			// here. Then drain stderr to keep the child from blocking on a
			// stderr write before exit. `gpl-boundary --list-tools` doesn't
			// write stderr in the success path (see GPL-boundary
			// `src/main.rs:88-94`), but error paths do; the second drain is
			// defensive and bounded.
			auto drain = [](int fd, std::string &out) {
				if (fd < 0) {
					return;
				}
				char buf[256];
				ssize_t n;
				while ((n = ::read(fd, buf, sizeof(buf))) > 0) {
					out.append(buf, static_cast<size_t>(n));
				}
			};
			std::string out;
			std::string err;
			drain(child.stdout_fd(), out);
			drain(child.stderr_fd(), err);
			const int status = child.Wait();
			// Both bowtie2-align (alignment) and bowtie2-build (index build)
			// are required; we use both downstream. Require both in the
			// registry rather than just `bowtie2-align` so a stripped-down
			// daemon doesn't pass this check then fail at bind time.
			available = WIFEXITED(status) && WEXITSTATUS(status) == 0 &&
			            out.find("bowtie2-align") != std::string::npos &&
			            out.find("bowtie2-build") != std::string::npos;
		} catch (...) {
			available = false;
		}
	});

	result.SetVectorType(VectorType::CONSTANT_VECTOR);
	auto &constant_validity = ConstantVector::Validity(result);
	constant_validity.SetAllValid(1);
	*ConstantVector::GetData<bool>(result) = available;
}

void RegisterBowtie2AvailableFunction(ExtensionLoader &loader) {
	ScalarFunction bowtie2_available("bowtie2_available", {}, LogicalType::BOOLEAN, Bowtie2AvailableFunction);
	loader.RegisterFunction(bowtie2_available);
}

} // namespace duckdb
