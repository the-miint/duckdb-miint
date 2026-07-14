// =============================================================================
// sylph_profile() — DuckDB table function for FracMinHash relative-abundance
// profiling of shotgun metagenomic reads.
//
// Wraps Shaw & Yu's sylph (Nature Biotechnology, 2024) via its Rust C/Arrow
// FFI. The function accepts a paired-end reads table/view and an on-disk
// `.syldb` reference database, sketches the reads in-memory through the
// FFI's streaming builder, runs sylph_profile, and zero-copies the resulting
// Arrow batch into DuckDB output Vectors via ArrowToDuckDBConversion (same
// pattern as the rype_* table functions; see internals/arrow-zero-copy.md).
//
// Two execution paths are dispatched in Execute() based on whether
// `sample_id` is set in Bind():
//
//   1. Single-sample path: GlobalState owns a single SylphSketch* + a
//      QuerySequenceStream draining the source table, and an ArrowOutputState
//      that holds the FFI's returned Arrow batch. RunProfile runs once on
//      the first Execute() call; subsequent calls page through the cached
//      batch via EmitFromArrow.
//
//   2. Per-sample path (sample_id := 'col'): GlobalState extends
//      PerSampleGlobalState to drive the standard claim-next-sample
//      iterator. Each LocalState owns a per-thread Connection (so it can
//      CREATE OR REPLACE TEMP VIEW __sylph_per_sample filtered to the
//      claimed sample), its own sketch, and its own ArrowOutputState.
//      Output prepends a sample_id column matching the source column's
//      type. Cross-sample parallelism is real throughput here — sylph's
//      profile path is mutex-free and the threading auto-balance in
//      BuildProfileParams keeps total core utilization aligned with the
//      DuckDB scheduler thread budget.
//
// Output schema is locked: adding columns at the end is non-breaking;
// reordering is not. See SetupBaseOutputSchema() and the RegisterDocumented*
// description string at the bottom for the user-facing contract.
// =============================================================================

#ifdef MIINT_HAS_SYLPH

#include "sylph_profile.hpp"

#include "SylphDatabase.hpp"
#include "sylph.h"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/parallel/task_scheduler.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <cstring>
#include <new>
#include <optional>
#include <stdexcept>

// DuckDB ships ArrowArray / ArrowSchema (the Arrow C Data Interface ABI
// structs) in duckdb/common/arrow/arrow.hpp — re-use them rather than
// declaring our own to avoid duplicate-definition errors.
#include "duckdb/common/arrow/arrow.hpp"

namespace duckdb {

// =============================================================================
// Output schema
// =============================================================================

// 9 sylph profile output columns; per-sample mode prepends a 10th sample_id
// column. Adding columns at the end is non-breaking; reordering is.
static void SetupBaseOutputSchema(vector<LogicalType> &types, vector<std::string> &names) {
	names = {
	    "genome_index", "genome_name", "contig_name", "sequence_abundance", "taxonomic_abundance",
	    "adjusted_ani", "eff_cov",     "naive_ani",   "kmers_reassigned",
	};
	types = {
	    LogicalType::UINTEGER, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::DOUBLE,  LogicalType::DOUBLE,
	    LogicalType::DOUBLE,   LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::UBIGINT,
	};
}

// =============================================================================
// State destructors
// =============================================================================
SylphProfileTableFunction::GlobalState::~GlobalState() {
	if (sketch != nullptr) {
		sylph_sketch_free(sketch);
		sketch = nullptr;
	}
}

SylphProfileTableFunction::LocalState::~LocalState() {
	if (sketch != nullptr) {
		sylph_sketch_free(sketch);
		sketch = nullptr;
	}
}

// =============================================================================
// Bind
// =============================================================================
unique_ptr<FunctionData> SylphProfileTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                         vector<LogicalType> &return_types,
                                                         vector<std::string> &names) {
	auto data = make_uniq<Data>();
	data->source_table = input.inputs[0].GetValue<std::string>();
	data->syldb_path = input.inputs[1].GetValue<std::string>();

	// Accept BIGINT/UUID read_id (allow_bigint=true), consistent with the other
	// tools and with sylph_index_create — so a Qiita BIGINT-keyed reads table
	// works index->profile end to end. QuerySequenceStream threads schema.id_type
	// and stringifies ids via the id-column codec; profiling itself is on the
	// sequences, so the id type only affects how read_id is read, not the result.
	data->schema = ValidateSequenceTableSchema(context, data->source_table, /*allow_bigint=*/true);

	// Seed FFI parameter structs from sylph defaults — single source of truth.
	if (sylph_profile_params_default(&data->profile_params) != 0) {
		throw IOException("sylph_profile: sylph_profile_params_default failed");
	}
	if (sylph_sketch_params_default(&data->sketch_params) != 0) {
		throw IOException("sylph_profile: sylph_sketch_params_default failed");
	}

	// Named-parameter overrides. Each helper returns the user value if the
	// param was passed, std::nullopt otherwise — so we only touch FFI fields
	// the caller explicitly asked to change.
	auto get_double = [&](const std::string &param, double lo, double hi) -> std::optional<double> {
		auto it = input.named_parameters.find(param);
		if (it == input.named_parameters.end())
			return std::nullopt;
		double v = it->second.GetValue<double>();
		if (v < lo || v > hi) {
			throw InvalidInputException("sylph_profile: %s must be in [%g, %g] (got %g)", param.c_str(), lo, hi, v);
		}
		return v;
	};
	auto get_uint = [&](const std::string &param) -> std::optional<uint32_t> {
		auto it = input.named_parameters.find(param);
		if (it == input.named_parameters.end())
			return std::nullopt;
		auto v = it->second.GetValue<int64_t>();
		if (v < 0) {
			throw InvalidInputException("sylph_profile: %s must be >= 0 (got %lld)", param.c_str(), (long long)v);
		}
		return static_cast<uint32_t>(v);
	};
	auto get_bool = [&](const std::string &param) -> std::optional<bool> {
		auto it = input.named_parameters.find(param);
		if (it == input.named_parameters.end())
			return std::nullopt;
		return it->second.GetValue<bool>();
	};

	// min_ani: user passes a fraction [0,1]; FFI takes percent [0,100].
	if (auto v = get_double("min_ani", 0.0, 1.0))
		data->profile_params.minimum_ani = *v * 100.0;
	if (auto v = get_uint("min_number_kmers"))
		data->profile_params.min_number_kmers = static_cast<double>(*v);
	if (auto v = get_double("min_count_correct", 0.0, 1e9))
		data->profile_params.min_count_correct = *v;
	if (auto v = get_bool("estimate_unknown"))
		data->profile_params.estimate_unknown = *v ? 1 : 0;
	if (auto v = get_bool("dedup_paired_reads"))
		data->sketch_params.dedup = *v ? 1 : 0;
	if (auto v = get_double("dedup_fpr", 0.0, 1.0))
		data->sketch_params.dedup_fpr = *v;
	if (auto v = get_uint("threads"))
		data->user_threads = *v;

	// Per-sample mode: explicit sample_id parameter wins; empty string disables.
	auto sample_it = input.named_parameters.find("sample_id");
	if (sample_it != input.named_parameters.end()) {
		auto val = sample_it->second.GetValue<std::string>();
		if (!val.empty()) {
			data->has_sample_id = true;
			data->sample_info.sample_id_col = val;
		}
	}

	SetupBaseOutputSchema(data->output_types, data->output_names);

	// Per-sample mode prepends a sample_id column to the output. DiscoverSamples
	// validates the column, captures its type, and collects the distinct values.
	if (data->has_sample_id) {
		Connection conn(*context.db);
		std::vector<std::string> reserved_lower;
		reserved_lower.reserve(data->output_names.size());
		for (auto &n : data->output_names) {
			reserved_lower.push_back(StringUtil::Lower(n));
		}
		DiscoverSamples(conn, data->source_table, data->sample_info.sample_id_col, reserved_lower, "sylph_profile",
		                data->sample_info);
		// Prepend sample_id column to outputs.
		vector<std::string> out_names = {"sample_id"};
		vector<LogicalType> out_types = {data->sample_info.sample_id_type};
		out_names.insert(out_names.end(), data->output_names.begin(), data->output_names.end());
		out_types.insert(out_types.end(), data->output_types.begin(), data->output_types.end());
		data->output_names = std::move(out_names);
		data->output_types = std::move(out_types);
	}

	names = data->output_names;
	return_types = data->output_types;
	return std::move(data);
}

// =============================================================================
// InitGlobal
// =============================================================================
unique_ptr<GlobalTableFunctionState> SylphProfileTableFunction::InitGlobal(ClientContext &context,
                                                                           TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	try {
		gstate->db = std::make_unique<miint::SylphDatabaseHandle>(data.syldb_path);
	} catch (const std::runtime_error &e) {
		throw IOException("sylph_profile: %s", e.what());
	}

	if (data.has_sample_id) {
		// Per-sample mode: max_threads_hint=0 — let DuckDB's scheduler
		// run as many samples in parallel as it has worker threads (clamped
		// to num_samples). Each sample's inner sylph_profile call gets a
		// rayon thread count auto-balanced in BuildProfileParams so that
		// (outer × inner) ≈ db_threads — same total core utilization as
		// `sylph profile -t db_threads`. NOTE: this differs from
		// sortmerna's max_threads_hint=1 because sortmerna has a
		// process-wide g_run_mutex that serializes everything; sylph's
		// profile path is mutex-free, so cross-sample parallelism is real
		// throughput.
		InitPerSampleGlobal(context, *gstate, data.sample_info.sample_values.size(), /*max_threads_hint=*/0);
	} else {
		// Single-sample mode: build the streaming reader and an empty sketch.
		gstate->stream =
		    std::make_unique<QuerySequenceStream>(context, data.source_table, data.schema, STANDARD_VECTOR_SIZE);

		gstate->sketch = sylph_sketch_builder_create(&data.sketch_params);
		if (gstate->sketch == nullptr) {
			const char *err = sylph_get_last_error();
			throw IOException("sylph_profile: sketch builder create failed: %s", err ? err : "<unknown>");
		}
	}

	return std::move(gstate);
}

// =============================================================================
// InitLocal — per-sample mode only allocates per-thread state.
// =============================================================================
unique_ptr<LocalTableFunctionState> SylphProfileTableFunction::InitLocal(ExecutionContext &context,
                                                                         TableFunctionInitInput &input,
                                                                         GlobalTableFunctionState *) {
	auto &data = input.bind_data->Cast<Data>();
	auto lstate = make_uniq<LocalState>();
	if (data.has_sample_id) {
		lstate->conn = make_uniq<Connection>(*context.client.db);
	}
	return std::move(lstate);
}

// =============================================================================
// ArrowOutputState methods — destructor + Reset + GetArrayState
// =============================================================================
SylphProfileTableFunction::ArrowOutputState::~ArrowOutputState() {
	Reset();
}

void SylphProfileTableFunction::ArrowOutputState::Reset() {
	// Tear-down order matters (mirrors RYpe's pattern):
	//   1. Drop the wrapper first — DuckDB Vectors that ref-counted into
	//      arrow_array via owned_data have already released by the time
	//      we get here at end-of-Execute or end-of-sample.
	//   2. Clear arrow_table next — it holds pointers into output_schema's
	//      data, so it must be reset before the schema is released.
	//   3. Release output_schema last.
	current_chunk.reset();
	batch_offset = 0;
	array_states.clear();
	arrow_table = ArrowTableSchema();
	if (output_schema.release) {
		output_schema.release(&output_schema);
		output_schema = ArrowSchema {};
	}
}

ArrowArrayScanState &SylphProfileTableFunction::ArrowOutputState::GetArrayState(idx_t col_idx, ClientContext &context) {
	auto it = array_states.find(col_idx);
	if (it == array_states.end()) {
		auto state = make_uniq<ArrowArrayScanState>(context);
		auto &ref = *state;
		array_states.emplace(col_idx, std::move(state));
		return ref;
	}
	return *it->second;
}

namespace {

// =============================================================================
// Helpers
// =============================================================================

// Resolve the sylph_profile FFI params for one Execute call. data.profile_params
// already has sylph's defaults plus any user named-param overrides from Bind();
// we just layer on miint's threading auto-balance: per-sample mode splits
// db_threads across outer DuckDB workers × inner rayon threads (mirrors sylph
// CLI's `step` formula), single-sample leaves num_threads at the FFI default
// (0 = sylph's global rayon pool).
static SylphProfileParams BuildProfileParams(const SylphProfileTableFunction::Data &data, idx_t outer_threads,
                                             idx_t db_threads) {
	SylphProfileParams pp = data.profile_params;
	if (data.user_threads != 0) {
		pp.num_threads = data.user_threads;
	} else if (data.has_sample_id) {
		idx_t safe_outer = std::max<idx_t>(outer_threads, 1);
		idx_t inner = std::max<idx_t>(1, db_threads / safe_outer);
		pp.num_threads = static_cast<uint32_t>(inner);
	}
	return pp;
}

// Drain `stream` into the FFI sketch builder, finalize, run sylph_profile,
// and stash the resulting Arrow batch + schema in `arrow`. Releases the
// sketch on success or throw.
static void RunProfile(QuerySequenceStream &stream, ::SylphSketch *&sketch, const ::SylphDatabase *db,
                       const SylphProfileParams &pp, SylphProfileTableFunction::ArrowOutputState &arrow,
                       ClientContext &context) {
	while (true) {
		auto batch = stream.FetchSubBatch();
		if (batch.empty()) {
			break;
		}
		for (size_t i = 0; i < batch.size(); i++) {
			const auto &r1 = batch.sequences1[i];
			if (r1.empty()) {
				continue;
			}
			const unsigned char *r1p = reinterpret_cast<const unsigned char *>(r1.data());
			size_t r1len = r1.size();
			const unsigned char *r2p = nullptr;
			size_t r2len = 0;
			if (batch.is_paired && i < batch.sequences2.size() && !batch.sequences2[i].empty()) {
				r2p = reinterpret_cast<const unsigned char *>(batch.sequences2[i].data());
				r2len = batch.sequences2[i].size();
			}
			int rc = sylph_sketch_builder_add_pair(sketch, r1p, r1len, r2p, r2len);
			if (rc != 0) {
				const char *err = sylph_get_last_error();
				throw IOException("sylph_profile: add_pair failed: %s", err ? err : "<unknown>");
			}
		}
	}
	if (sylph_sketch_builder_finalize(sketch) != 0) {
		const char *err = sylph_get_last_error();
		throw IOException("sylph_profile: finalize failed: %s", err ? err : "<unknown>");
	}

	// Tear down any prior batch (per-sample mode reuses the state across
	// samples) before installing the new one.
	arrow.Reset();

	auto wrapper = make_shared_ptr<ArrowArrayWrapper>();
	int rc = sylph_profile(db, sketch, &pp, &wrapper->arrow_array, &arrow.output_schema);
	sylph_sketch_free(sketch);
	sketch = nullptr;
	if (rc != 0) {
		const char *err = sylph_get_last_error();
		throw IOException("sylph_profile: profile failed: %s", err ? err : "<unknown>");
	}
	if (wrapper->arrow_array.n_children != 9) {
		throw IOException("sylph_profile: FFI returned %lld columns, expected 9",
		                  (long long)wrapper->arrow_array.n_children);
	}

	// Parse the schema for DuckDB's converter; ArrowTableSchema holds
	// pointers into output_schema's data, so it must be (re)populated each
	// time output_schema is reset above.
	ArrowTableFunction::PopulateArrowTableSchema(context, arrow.arrow_table, arrow.output_schema);
	arrow.current_chunk = std::move(wrapper);
	arrow.batch_offset = 0;
}

// Emit up to STANDARD_VECTOR_SIZE rows from `arrow.current_chunk` into
// `output`, starting at column index `start_col` (1 in per-sample mode where
// column 0 is the sample_id; 0 otherwise). Uses DuckDB's Arrow→DuckDB
// conversion: zero-copy for fixed-width columns, copy for the two strings.
// Returns 0 when the chunk is fully drained.
static idx_t EmitFromArrow(SylphProfileTableFunction::ArrowOutputState &arrow, DataChunk &output, idx_t start_col,
                           ClientContext &context) {
	if (!arrow.current_chunk) {
		return 0;
	}
	auto &batch = arrow.current_chunk->arrow_array;
	const idx_t total = static_cast<idx_t>(batch.length);
	if (arrow.batch_offset >= total) {
		return 0;
	}
	const idx_t to_emit = MinValue<idx_t>(total - arrow.batch_offset, STANDARD_VECTOR_SIZE);
	auto &arrow_columns = arrow.arrow_table.GetColumns();

	for (idx_t c = 0; c < 9; c++) {
		auto &array = *batch.children[c];
		auto &arrow_type = *arrow_columns.at(c);
		auto &array_state = arrow.GetArrayState(c, context);
		// owned_data ref-counts the wrapper so any zero-copy buffers in the
		// output Vector outlive this Execute() call.
		array_state.owned_data = arrow.current_chunk;
		ArrowToDuckDBConversion::SetValidityMask(output.data[start_col + c], array, arrow.batch_offset, to_emit,
		                                         batch.offset, -1);
		ArrowToDuckDBConversion::ColumnArrowToDuckDB(output.data[start_col + c], array, arrow.batch_offset, array_state,
		                                             to_emit, arrow_type);
	}
	arrow.batch_offset += to_emit;
	return to_emit;
}

} // namespace

// =============================================================================
// Execute
// =============================================================================
void SylphProfileTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();
	auto &data = data_p.bind_data->Cast<Data>();

	// Resolve thread split once per Execute call. db_threads is the underlying
	// scheduler thread budget; gstate.max_threads is the outer per-sample
	// concurrency (set by InitPerSampleGlobal, capped to num_samples).
	const idx_t db_threads = static_cast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	const SylphProfileParams pp = BuildProfileParams(data, gstate.max_threads, db_threads);

	if (!data.has_sample_id) {
		// ---- Single-sample path ----
		if (!gstate.profile_done) {
			RunProfile(*gstate.stream, gstate.sketch, gstate.db->raw(), pp, gstate.arrow, context);
			gstate.profile_done = true;
		}
		idx_t emitted = EmitFromArrow(gstate.arrow, output, /*start_col=*/0, context);
		output.SetCardinality(emitted);
		return;
	}

	// ---- Per-sample path ----
	auto &lstate = data_p.local_state->Cast<LocalState>();
	const idx_t num_samples = data.sample_info.sample_values.size();

	// Loop: drain current sample's buffered batch, then claim the next sample.
	// A sample that produces zero passing genomes (every genome fell below the
	// ANI cutoff) must NOT short-circuit the function — we have to claim the
	// next sample on the same Execute call, otherwise DuckDB takes the empty
	// chunk as "done" and skips remaining samples.
	while (true) {
		idx_t emitted = EmitFromArrow(lstate.arrow, output, /*start_col=*/1, context);
		if (emitted > 0) {
			output.data[0].Reference(lstate.sample_value);
			output.SetCardinality(emitted);
			return;
		}

		idx_t sample_idx = 0;
		if (!ClaimNextSample(gstate, num_samples, sample_idx)) {
			output.SetCardinality(0);
			return;
		}
		lstate.sample_value = data.sample_info.sample_values[sample_idx];

		// Per-sample TEMP VIEW filtered to this sample. CAST-as-VARCHAR equality
		// matches uchime_ref's per-sample path. The view name is fixed because
		// each LocalState owns its own connection — TEMP VIEWs are scoped to the
		// connection, so concurrent threads don't collide.
		auto src = KeywordHelper::WriteOptionallyQuoted(data.source_table);
		auto col = KeywordHelper::WriteOptionallyQuoted(data.sample_info.sample_id_col);
		auto sample_lit = lstate.sample_value.ToSQLString();
		auto view_sql = "CREATE OR REPLACE TEMP VIEW __sylph_per_sample AS SELECT * FROM " + src + " WHERE CAST(" +
		                col + " AS VARCHAR) = CAST(" + sample_lit + " AS VARCHAR)";
		auto view_result = lstate.conn->Query(view_sql);
		if (view_result->HasError()) {
			throw InvalidInputException("sylph_profile: failed to create per-sample view for %s: %s", sample_lit,
			                            view_result->GetError());
		}
		lstate.stream = std::make_unique<QuerySequenceStream>(*lstate.conn, "__sylph_per_sample", data.schema,
		                                                      STANDARD_VECTOR_SIZE);

		lstate.sketch = sylph_sketch_builder_create(&data.sketch_params);
		if (lstate.sketch == nullptr) {
			const char *err = sylph_get_last_error();
			throw IOException("sylph_profile: sketch builder create failed: %s", err ? err : "<unknown>");
		}

		RunProfile(*lstate.stream, lstate.sketch, gstate.db->raw(), pp, lstate.arrow, context);
		// Loop back: emit if non-empty, claim next sample if empty.
	}
}

// =============================================================================
// Registration
// =============================================================================
TableFunction SylphProfileTableFunction::GetFunction() {
	auto tf = TableFunction("sylph_profile", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal,
	                        InitLocal);

	tf.named_parameters["sample_id"] = LogicalType::VARCHAR;
	tf.named_parameters["min_ani"] = LogicalType::DOUBLE;
	tf.named_parameters["min_number_kmers"] = LogicalType::UINTEGER;
	tf.named_parameters["min_count_correct"] = LogicalType::DOUBLE;
	tf.named_parameters["estimate_unknown"] = LogicalType::BOOLEAN;
	tf.named_parameters["dedup_paired_reads"] = LogicalType::BOOLEAN;
	tf.named_parameters["dedup_fpr"] = LogicalType::DOUBLE;
	tf.named_parameters["threads"] = LogicalType::UINTEGER;

	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void SylphProfileTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb

#endif // MIINT_HAS_SYLPH
