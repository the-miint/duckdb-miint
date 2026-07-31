#include "rype_index_create.hpp"
#include "catalog_utils.hpp"
#include "rype_common.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <algorithm>
#include <utility>

namespace duckdb {

namespace {

// A hand-rolled multiplexing ArrowArrayStream that splices several windowed query
// results into one continuous stream for rype_index_build_from_arrow.
//
// Why windows: building from a single `SELECT ... ORDER BY feature_idx,
// chunk_index` sorts the entire corpus of 64 KB chunk BLOBs in one blocking
// operator, which DuckDB cannot spill — it OOMs on real data (many genomes,
// measured: OOM at the memory limit). Pinning the scan to one thread to stream
// pre-sorted data without a sort is not an option either: `threads` is global and
// resizing it from InitGlobal (which runs on a scheduler worker thread) joins the
// current worker — "Resource deadlock avoided".
//
// Instead we partition the features into windows (contiguous feature_idx ranges)
// and run one `... WHERE feature_idx >= lo AND feature_idx <= hi ORDER BY
// feature_idx, chunk_index` per window. Each window's sort is bounded (multi-
// threaded, spillable); the windows are emitted back-to-back as one stream. RYpe's
// ordering requirement — each feature contiguous, ascending, 0-based, gap-free —
// holds: a feature lives wholly inside one window, windows ascend by feature_idx,
// and the per-window ORDER BY groups each feature's chunks. The chunk_table may be
// in ANY physical order. Only one window's StreamQueryResult is open at a time,
// satisfying the one-active-stream-per-Connection rule.
//
// Ownership mirrors ResultArrowArrayStreamWrapper: rype_index_build_from_arrow
// takes the stream and invokes release() synchronously during the build, which
// deletes this object.
struct WindowedChunkStream {
	ArrowArrayStream stream;
	Connection &conn;
	std::string table_quoted;
	std::vector<std::pair<int64_t, int64_t>> windows; // inclusive [lo, hi] feature_idx ranges, ascending
	idx_t batch_size;
	idx_t next_window = 0; // index of the next window to open
	unique_ptr<ResultArrowArrayStreamWrapper> active;
	std::string last_error;

	WindowedChunkStream(Connection &conn_p, std::string table_quoted_p,
	                    std::vector<std::pair<int64_t, int64_t>> windows_p, idx_t batch_size_p)
	    : conn(conn_p), table_quoted(std::move(table_quoted_p)), windows(std::move(windows_p)),
	      batch_size(batch_size_p) {
		stream.private_data = this;
		stream.get_schema = GetSchema;
		stream.get_next = GetNext;
		stream.release = Release;
		stream.get_last_error = GetLastError;
	}

	bool OpenQuery(const std::string &sql) {
		auto result = conn.SendQuery(sql);
		if (result->HasError()) {
			last_error = result->GetError();
			return false;
		}
		active = make_uniq<ResultArrowArrayStreamWrapper>(std::move(result), batch_size);
		return true;
	}

	std::string WindowSQL(int64_t lo, int64_t hi) const {
		return "SELECT feature_idx::BIGINT AS feature_idx, chunk_index::INTEGER AS chunk_index, chunk_data FROM " +
		       table_quoted + " WHERE feature_idx >= " + std::to_string(lo) +
		       " AND feature_idx <= " + std::to_string(hi) + " ORDER BY feature_idx, chunk_index";
	}

	// Open windows[next_window] (advancing next_window). Caller guarantees the index
	// is in range.
	bool OpenNextWindow() {
		auto &w = windows[next_window++];
		return OpenQuery(WindowSQL(w.first, w.second));
	}

	// Zero-row probe used only when there are no windows (empty chunk_table) so that
	// get_schema can still return a valid schema.
	bool OpenSchemaProbe() {
		return OpenQuery("SELECT feature_idx::BIGINT AS feature_idx, chunk_index::INTEGER AS chunk_index, "
		                 "chunk_data FROM " +
		                 table_quoted + " WHERE false");
	}

	static int GetSchema(ArrowArrayStream *stream, ArrowSchema *out) {
		auto self = reinterpret_cast<WindowedChunkStream *>(stream->private_data);
		if (!self->active) {
			bool ok = self->windows.empty() ? self->OpenSchemaProbe() : self->OpenNextWindow();
			if (!ok) {
				return -1;
			}
		}
		return self->active->stream.get_schema(&self->active->stream, out);
	}

	static int GetNext(ArrowArrayStream *stream, ArrowArray *out) {
		auto self = reinterpret_cast<WindowedChunkStream *>(stream->private_data);
		while (true) {
			if (!self->active) {
				if (self->next_window >= self->windows.size()) {
					out->release = nullptr; // end of the multiplexed stream
					return 0;
				}
				if (!self->OpenNextWindow()) {
					return -1; // last_error set by OpenQuery
				}
			}
			ArrowArray tmp;
			tmp.release = nullptr;
			if (self->active->stream.get_next(&self->active->stream, &tmp) != 0) {
				self->last_error = self->active->stream.get_last_error(&self->active->stream);
				return -1;
			}
			if (tmp.release != nullptr) {
				*out = tmp; // a batch from the current window
				return 0;
			}
			self->active.reset(); // window exhausted; advance to the next
		}
	}

	static void Release(ArrowArrayStream *stream) {
		if (!stream || !stream->release) {
			return;
		}
		stream->release = nullptr;
		delete reinterpret_cast<WindowedChunkStream *>(stream->private_data);
	}

	static const char *GetLastError(ArrowArrayStream *stream) {
		auto self = reinterpret_cast<WindowedChunkStream *>(stream->private_data);
		return self->last_error.c_str();
	}
};

// Choose how many features go in each feed window. Targets ~256 MiB of chunk data
// per window (the sort's working set inflates several-fold over the raw bytes, so a
// modest data budget keeps the per-window sort comfortably bounded) by sampling the
// average feature size over the first few features (feature_idx and chunk_data only
// -> cheap). Returns at least 1 so a single oversized feature still makes progress.
// `features` is the sorted distinct feature_idx list.
idx_t AutoWindowFeatures(Connection &conn, const std::string &table_quoted, const std::vector<int64_t> &features) {
	constexpr idx_t kSampleFeatures = 128;
	constexpr int64_t kBudgetBytes = 256LL << 20; // ~256 MiB of chunk data per window
	constexpr idx_t kFallback = 256;
	if (features.empty()) {
		return 1;
	}
	idx_t k = std::min(kSampleFeatures, static_cast<idx_t>(features.size()));
	std::string sql = "SELECT AVG(b)::BIGINT FROM (SELECT SUM(LENGTH(chunk_data))::BIGINT AS b FROM " + table_quoted +
	                  " WHERE feature_idx >= " + std::to_string(features.front()) +
	                  " AND feature_idx <= " + std::to_string(features[k - 1]) + " GROUP BY feature_idx)";
	auto result = conn.Query(sql);
	if (result->HasError()) {
		return kFallback;
	}
	auto value = result->GetValue(0, 0);
	if (value.IsNull()) {
		return kFallback;
	}
	int64_t avg_bytes = value.GetValue<int64_t>();
	if (avg_bytes <= 0) {
		return kFallback;
	}
	return static_cast<idx_t>(std::max<int64_t>(1, kBudgetBytes / avg_bytes));
}

} // namespace

// ============================================================================
// Bind
// ============================================================================
unique_ptr<FunctionData> RypeIndexCreateTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                            vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<Data>();

	// Required positional parameters: chunk_table, output_path
	if (input.inputs.size() < 2) {
		throw BinderException("rype_index_create requires chunk_table and output_path parameters");
	}
	data->chunk_table = input.inputs[0].ToString();
	data->output_path = input.inputs[1].ToString();

	// Optional: mapping_table (feature_idx, bucket_name). Empty => single bucket.
	auto mapping_param = input.named_parameters.find("mapping_table");
	if (mapping_param != input.named_parameters.end() && !mapping_param->second.IsNull()) {
		data->mapping_table = mapping_param->second.ToString();
	}

	auto k_param = input.named_parameters.find("k");
	if (k_param != input.named_parameters.end() && !k_param->second.IsNull()) {
		data->k = k_param->second.GetValue<int32_t>();
	}

	auto w_param = input.named_parameters.find("w");
	if (w_param != input.named_parameters.end() && !w_param->second.IsNull()) {
		data->w = w_param->second.GetValue<int32_t>();
	}

	auto salt_param = input.named_parameters.find("salt");
	if (salt_param != input.named_parameters.end() && !salt_param->second.IsNull()) {
		data->salt = salt_param->second.GetValue<uint64_t>();
	}

	auto orient_param = input.named_parameters.find("orient");
	if (orient_param != input.named_parameters.end() && !orient_param->second.IsNull()) {
		data->orient = orient_param->second.GetValue<bool>();
	}

	auto mem_param = input.named_parameters.find("max_memory");
	if (mem_param != input.named_parameters.end() && !mem_param->second.IsNull()) {
		data->max_memory = mem_param->second.GetValue<int64_t>();
	}

	// Advanced: number of features per feed window (see WindowedChunkStream).
	// 0 (default) = auto-size from a per-window memory budget. Mainly useful for
	// tuning or for tests that need to force multiple windows over small inputs.
	auto window_param = input.named_parameters.find("feed_window_features");
	if (window_param != input.named_parameters.end() && !window_param->second.IsNull()) {
		data->feed_window_features = window_param->second.GetValue<int64_t>();
		if (data->feed_window_features < 0) {
			throw BinderException("rype_index_create: feed_window_features must be >= 0 (got %lld)",
			                      static_cast<long long>(data->feed_window_features));
		}
	}

	// Validate at bind time (fail fast, before any build side effects). The chunk
	// table is checked first so a bad chunk table is not misattributed to the
	// synthesized mapping query (which is derived from it).
	ValidateTableHasColumns(context, data->chunk_table, {"feature_idx", "chunk_index", "chunk_data"}, "chunk table");
	if (!data->mapping_table.empty()) {
		ValidateTableHasColumns(context, data->mapping_table, {"feature_idx", "bucket_name"}, "mapping table");
	}

	// k must be one of RYpe's supported RY-space k-mer sizes. RYpe also validates
	// this, but checking here fails before the (non-atomic) build writes anything.
	if (data->k != 16 && data->k != 32 && data->k != 64) {
		throw BinderException("rype_index_create: k must be 16, 32, or 64 (got %d)", data->k);
	}

	// w is the minimizer window size and must be >= 1. RYpe does NOT validate this:
	// w <= 0 silently builds a minimizer-free, useless index that still reports
	// success, and a negative w would wrap to a huge size_t via the FFI cast.
	if (data->w < 1) {
		throw BinderException("rype_index_create: w must be >= 1 (got %d)", data->w);
	}

	// Set output schema
	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

// ============================================================================
// InitGlobal — performs the build (a synchronous side effect)
// ============================================================================
unique_ptr<GlobalTableFunctionState> RypeIndexCreateTableFunction::InitGlobal(ClientContext &context,
                                                                              TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	// Sub-connection to read the input tables as Arrow (avoids the context.Query
	// re-entrancy deadlock). Held in GlobalState so it outlives RYpe's lazy
	// consumption of the streamed chunk cursor during the synchronous build.
	auto &db = DatabaseInstance::GetDatabase(context);
	gstate->input_connection = make_uniq<Connection>(db);
	InheritTempObjects(context, *gstate->input_connection);
	auto &conn = *gstate->input_connection;

	// Arrow v1.4 (Utf8View/BinaryView) for VARCHAR — no i32 offset 2 GiB cap. The
	// rype build path accepts Utf8/LargeUtf8/Binary/LargeBinary/Utf8View/BinaryView.
	conn.Query("SET arrow_output_version = '1.4'");

	std::string chunk_quoted = KeywordHelper::WriteOptionallyQuoted(bind_data.chunk_table);

	// Mapping stream: feature_idx Int64, bucket_name Utf8. Materialized — it is
	// small (one row per feature) and RYpe reads the whole mapping into memory
	// before processing chunks, so streaming it would buy nothing. Done BEFORE the
	// streaming cursor below so it does not occupy the connection's stream slot.
	// When no mapping_table is supplied, every feature goes into a single bucket
	// named 'unnamed-bucket', synthesized from the distinct feature_idx in the
	// chunk table.
	std::string mapping_sql;
	if (bind_data.mapping_table.empty()) {
		mapping_sql = "SELECT DISTINCT feature_idx::BIGINT AS feature_idx, "
		              "'unnamed-bucket'::VARCHAR AS bucket_name FROM " +
		              chunk_quoted;
	} else {
		std::string mapping_quoted = KeywordHelper::WriteOptionallyQuoted(bind_data.mapping_table);
		mapping_sql =
		    "SELECT feature_idx::BIGINT AS feature_idx, bucket_name::VARCHAR AS bucket_name FROM " + mapping_quoted;
	}
	auto mapping_result = conn.Query(mapping_sql);
	if (mapping_result->HasError()) {
		throw InvalidInputException("Failed to read mapping table '%s': %s",
		                            bind_data.mapping_table.empty() ? bind_data.chunk_table : bind_data.mapping_table,
		                            mapping_result->GetError());
	}
	auto mapping_wrapper = make_uniq<ResultArrowArrayStreamWrapper>(std::move(mapping_result), STANDARD_VECTOR_SIZE);

	// Chunk stream: feature_idx Int64, chunk_index Int32, chunk_data. RYpe requires
	// each feature's chunks to arrive contiguously in ascending, 0-based, gap-free
	// chunk_index order and reassembles them internally. We feed it as a sequence of
	// bounded, independently-sorted windows over feature_idx ranges (see
	// WindowedChunkStream): a single `ORDER BY` over the whole corpus OOMs because
	// DuckDB cannot spill the 64 KB BLOB sort payload. The chunk_table may be in any
	// physical order.

	// Sorted distinct feature ids (feature_idx only -> cheap projection pushdown).
	auto feat_result = conn.Query("SELECT DISTINCT feature_idx::BIGINT AS f FROM " + chunk_quoted + " ORDER BY f");
	if (feat_result->HasError()) {
		throw InvalidInputException("Failed to read feature ids from chunk table '%s': %s", bind_data.chunk_table,
		                            feat_result->GetError());
	}
	std::vector<int64_t> features;
	features.reserve(feat_result->RowCount());
	while (auto chunk = feat_result->Fetch()) {
		for (idx_t i = 0; i < chunk->size(); i++) {
			auto v = chunk->data[0].GetValue(i);
			if (!v.IsNull()) {
				features.push_back(v.GetValue<int64_t>());
			}
		}
	}

	// Window size in #features: explicit override (feed_window_features), else auto
	// from a per-window byte budget.
	idx_t window_features = bind_data.feed_window_features > 0 ? static_cast<idx_t>(bind_data.feed_window_features)
	                                                           : AutoWindowFeatures(conn, chunk_quoted, features);

	// Partition the sorted distinct ids into inclusive [lo, hi] feature_idx windows.
	std::vector<std::pair<int64_t, int64_t>> windows;
	for (idx_t i = 0; i < features.size(); i += window_features) {
		idx_t last = std::min(i + window_features, static_cast<idx_t>(features.size())) - 1;
		windows.emplace_back(features[i], features[last]);
	}

	auto chunk_mux = make_uniq<WindowedChunkStream>(conn, chunk_quoted, std::move(windows), STANDARD_VECTOR_SIZE);

	int rc =
	    rype_index_build_from_arrow(bind_data.output_path.c_str(), &chunk_mux->stream, &mapping_wrapper->stream,
	                                static_cast<size_t>(bind_data.k), static_cast<size_t>(bind_data.w), bind_data.salt,
	                                bind_data.orient ? 1 : 0, static_cast<size_t>(bind_data.max_memory));

	// RYpe takes ownership of both (non-NULL) streams and invokes their release
	// callbacks during the build — release the unique_ptrs to avoid a double-free.
	(void)chunk_mux.release();
	(void)mapping_wrapper.release();

	if (rc != 0) {
		const char *err = rype_get_last_error();
		throw IOException("rype_index_create failed to build '%s': %s", bind_data.output_path,
		                  err ? err : "unknown error");
	}

	return gstate;
}

// ============================================================================
// InitLocal
// ============================================================================
unique_ptr<LocalTableFunctionState> RypeIndexCreateTableFunction::InitLocal(ExecutionContext &context,
                                                                            TableFunctionInitInput &input,
                                                                            GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// ============================================================================
// Execute — emit the single status row
// ============================================================================
void RypeIndexCreateTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	if (gstate.done) {
		output.SetCardinality(0);
		return;
	}

	output.data[0].SetValue(0, Value(bind_data.output_path));
	output.data[1].SetValue(0, Value::INTEGER(bind_data.k));
	output.data[2].SetValue(0, Value::INTEGER(bind_data.w));
	output.data[3].SetValue(0, Value("ok"));

	output.SetCardinality(1);
	gstate.done = true;
}

// ============================================================================
// Function registration
// ============================================================================
TableFunction RypeIndexCreateTableFunction::GetFunction() {
	TableFunction tf("rype_index_create", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal,
	                 InitLocal);

	tf.named_parameters["mapping_table"] = LogicalType::VARCHAR;
	tf.named_parameters["k"] = LogicalType::INTEGER;
	tf.named_parameters["w"] = LogicalType::INTEGER;
	tf.named_parameters["salt"] = LogicalType::UBIGINT;
	tf.named_parameters["orient"] = LogicalType::BOOLEAN;
	tf.named_parameters["max_memory"] = LogicalType::BIGINT;
	tf.named_parameters["feed_window_features"] = LogicalType::BIGINT;

	return tf;
}

void RypeIndexCreateTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
