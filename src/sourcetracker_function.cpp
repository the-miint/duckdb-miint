#include "sourcetracker_function.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "sourcetracker_dataset.hpp"
#include "unifrac_function_common.hpp"

#include "duckdb/common/arrow/arrow.hpp"
#include "duckdb/common/arrow/arrow_appender.hpp"
#include "duckdb/common/arrow/arrow_converter.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table/arrow/arrow_duck_schema.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/client_properties.hpp"

// st3.h refers to the Arrow C Data Interface structs by pointer and does not
// declare them; duckdb/common/arrow/arrow.hpp above provides the definitions.
#include "st3.h"

namespace duckdb {
namespace {

using miint::sourcetracker::Dataset;
using unifrac_internal::ProbeFeatureTableIdType;
using unifrac_internal::ReadFeatureTable;
using unifrac_internal::ReadWideMetadata;
using unifrac_internal::ResolveThreadsParameter;

constexpr const char *kCaller = "sourcetracker";

// ---------------------------------------------------------------------------
// Bind: parameters and schema only. The relations are read, and the sampler
// run, in InitGlobal (docs/internals/architecture.md, "no work in Bind").
// ---------------------------------------------------------------------------

struct SourcetrackerBindData : public TableFunctionData {
	std::string table_name;
	std::string metadata_name;
	// Everything st3 needs except the seed, which is drawn per execution when the
	// user did not fix it.
	St3Config config {};
	int64_t seed = -1;
	bool loo = false;
	bool assignments = false;
	LogicalType sink_id_type = LogicalType::VARCHAR;
	LogicalType feature_id_type = LogicalType::VARCHAR;
};

// One (feature, mean count) assignment cell; the feature is an index into the
// dataset's feature dictionary.
using TallyCell = std::pair<int32_t, double>;

struct SourcetrackerGlobalState : public GlobalTableFunctionState {
	std::vector<std::string> sink_ids; // st3's row order
	std::vector<std::string> sources;  // st3's column order: environments, then Unknown
	std::vector<double> means;         // [sink][source], row-major
	std::vector<double> stds;          // same shape
	// Only when assignments were requested: one cell list per (sink, source),
	// same row-major order as the means, sorted by feature, nonzero cells only.
	bool assignments = false;
	std::vector<std::string> feature_ids;
	std::vector<std::vector<TallyCell>> tally;
	idx_t cursor = 0;
	LogicalType sink_id_type = LogicalType::VARCHAR;
	LogicalType feature_id_type = LogicalType::VARCHAR;
	idx_t MaxThreads() const override {
		return 1;
	}
};

// SourceTracker2's defaults (sourcetracker/_gibbs_defaults.py), which the SQL
// surface follows. `collapse` is what SourceTracker2 hard-codes in sink mode.
St3Config DefaultConfig() {
	St3Config c {};
	c.struct_version = ST3_CONFIG_V1;
	c.struct_size = static_cast<uint32_t>(sizeof(St3Config));
	c.seed = 0;
	c.jobs = 1;
	c.source_rarefaction_depth = 1000;
	c.sink_rarefaction_depth = 1000;
	c.with_replacement = 0;
	c.collapse = ST3_COLLAPSE_MEAN;
	c.loo = 0;
	c.contingency = 0;
	c.estimator = ST3_ESTIMATOR_KIND_GIBBS_COLLAPSED;
	c.alpha1 = 0.001;
	c.alpha2 = 0.1;
	c.beta = 10.0;
	c.restarts = 10;
	c.draws_per_restart = 1;
	c.burnin = 100;
	c.delay = 1;
	return c;
}

unique_ptr<FunctionData> SourcetrackerBind(ClientContext &context, TableFunctionBindInput &input,
                                           vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<SourcetrackerBindData>();
	data->table_name = input.inputs[0].GetValue<string>();
	data->metadata_name = input.inputs[1].GetValue<string>();
	if (data->table_name.empty()) {
		throw BinderException("sourcetracker: feature-table name must not be empty");
	}
	if (data->metadata_name.empty()) {
		throw BinderException("sourcetracker: sample-metadata name must not be empty");
	}

	St3Config &c = data->config;
	c = DefaultConfig();
	int32_t threads = 0;

	auto at_least_one = [](const std::string &key, const Value &v) -> uint32_t {
		const auto x = v.GetValue<int32_t>();
		if (x < 1) {
			throw BinderException("sourcetracker: %s must be >= 1 (got %s)", key, std::to_string(x));
		}
		return static_cast<uint32_t>(x);
	};
	auto depth = [](const std::string &key, const Value &v) -> int32_t {
		const auto x = v.GetValue<int32_t>();
		if (x < 0) {
			throw BinderException("sourcetracker: %s must be >= 0 (0 disables rarefaction; got %s)", key,
			                      std::to_string(x));
		}
		return x;
	};
	auto prior = [](const std::string &key, const Value &v) -> double {
		const auto x = v.GetValue<double>();
		if (!std::isfinite(x) || x < 0) {
			throw BinderException("sourcetracker: %s must be finite and >= 0 (got %s)", key, v.ToString());
		}
		return x;
	};

	for (const auto &kv : input.named_parameters) {
		const auto key = StringUtil::Lower(kv.first);
		if (key == "loo") {
			data->loo = kv.second.GetValue<bool>();
		} else if (key == "assignments") {
			data->assignments = kv.second.GetValue<bool>();
		} else if (key == "alpha1") {
			c.alpha1 = prior(key, kv.second);
		} else if (key == "alpha2") {
			c.alpha2 = prior(key, kv.second);
		} else if (key == "beta") {
			c.beta = prior(key, kv.second);
		} else if (key == "restarts") {
			c.restarts = at_least_one(key, kv.second);
		} else if (key == "draws_per_restart") {
			c.draws_per_restart = at_least_one(key, kv.second);
		} else if (key == "burnin") {
			c.burnin = at_least_one(key, kv.second);
		} else if (key == "delay") {
			c.delay = at_least_one(key, kv.second);
		} else if (key == "source_rarefaction_depth") {
			c.source_rarefaction_depth = depth(key, kv.second);
		} else if (key == "sink_rarefaction_depth") {
			c.sink_rarefaction_depth = depth(key, kv.second);
		} else if (key == "with_replacement") {
			c.with_replacement = kv.second.GetValue<bool>() ? 1 : 0;
		} else if (key == "collapse") {
			const auto mode = StringUtil::Lower(kv.second.GetValue<string>());
			if (mode == "mean") {
				c.collapse = ST3_COLLAPSE_MEAN;
			} else if (mode == "sum") {
				c.collapse = ST3_COLLAPSE_SUM;
			} else {
				throw BinderException("sourcetracker: collapse must be 'mean' or 'sum' (got '%s')",
				                      kv.second.GetValue<string>());
			}
		} else if (key == "seed") {
			data->seed = kv.second.GetValue<int64_t>();
			if (data->seed < -1) {
				throw BinderException("sourcetracker: seed must be >= 0, or -1 for a fresh draw (got %s)",
				                      std::to_string(data->seed));
			}
		} else if (key == "threads") {
			threads = kv.second.GetValue<int32_t>();
		}
	}
	if (data->loo && data->assignments) {
		throw BinderException("sourcetracker: assignments := true is not available with loo := true "
		                      "(leave-one-out produces no per-feature assignment tally)");
	}
	c.loo = data->loo ? 1 : 0;
	c.contingency = data->assignments ? 1 : 0;

	// st3's own `jobs = 0` means every core; DuckDB's thread count decides here.
	int jobs = ResolveThreadsParameter(context, threads, kCaller);
#ifdef __EMSCRIPTEN__
	// st3 builds a rayon pool per call whenever jobs != 1, and there are no
	// worker threads to build it from under wasm (unifrac_omp_scope.cpp does the
	// same for OpenMP).
	jobs = 1;
#endif
	c.jobs = jobs;

	// Schema, without reading a row: the LIMIT 0 probe proves the feature table's
	// columns exist and cast, and the catalog gives both id types.
	data->sink_id_type = ProbeFeatureTableIdType(context, data->table_name, kCaller, &data->feature_id_type);
	auto cols = GetTableOrViewColumns(context, data->metadata_name, "sample-metadata");
	for (const char *required : {"sample_id", "source_sink", "env"}) {
		if (!HasColumn(cols, required)) {
			throw BinderException("sourcetracker: sample-metadata '%s' must expose (sample_id, source_sink, env); "
			                      "column '%s' is missing",
			                      data->metadata_name, required);
		}
	}

	names.emplace_back("sink_id");
	return_types.emplace_back(data->sink_id_type);
	names.emplace_back("source");
	return_types.emplace_back(LogicalType::VARCHAR);
	names.emplace_back("proportion");
	return_types.emplace_back(LogicalType::DOUBLE);
	names.emplace_back("proportion_std");
	return_types.emplace_back(LogicalType::DOUBLE);
	names.emplace_back("assignments");
	return_types.emplace_back(LogicalType::MAP(data->feature_id_type, LogicalType::DOUBLE));
	return std::move(data);
}

// ---------------------------------------------------------------------------
// st3 error and handle plumbing
// ---------------------------------------------------------------------------

struct St3TableDeleter {
	void operator()(St3Table *t) const {
		st3_table_free(t);
	}
};
struct St3ResultDeleter {
	void operator()(St3Result *r) const {
		st3_result_free(r);
	}
};

// Every st3 entry point clears the thread-local error slot, so the message is
// copied out before anything else is called.
[[noreturn]] void ThrowSt3(const char *what, St3Status status) {
	const char *raw = st3_last_error();
	const std::string message = (raw && *raw) ? std::string(raw)
	                                          : std::string(what) + " failed with st3 status " +
	                                                std::to_string(static_cast<int>(status)) + " and no message";
	switch (status) {
	case ST3_STATUS_ERR_INVALID_INPUT:
	case ST3_STATUS_ERR_SHALLOW_SAMPLE:
	case ST3_STATUS_ERR_NO_SOURCES:
		throw InvalidInputException("sourcetracker: %s", message);
	default:
		// A panic or allocation failure inside st3 is a library bug, but it is not
		// DuckDB's: an InternalException would invalidate the whole database
		// (ClientContext marks it invalid on ExceptionType::INTERNAL), which is
		// out of proportion for one failed query.
		throw InvalidInputException("sourcetracker: st3 internal error during %s: %s", what, message);
	}
}

void CheckInterrupt(ClientContext &context) {
	if (context.interrupted) {
		throw InterruptException();
	}
}

// ---------------------------------------------------------------------------
// DuckDB -> Arrow marshaling
//
// An Arrow array exported by ArrowAppender and its schema, released on every
// exit path unless st3 has taken the array. st3_table_from_arrow consumes its
// three arrays on EVERY status once all six pointers are valid, leaving the
// caller's struct with a release pointer that must not be called again; Taken()
// records that.
// ---------------------------------------------------------------------------

struct ArrowExport {
	ArrowArray array;
	ArrowSchema schema;
	ArrowExport() {
		array.Init();
		schema.Init();
	}
	~ArrowExport() {
		if (array.release) {
			array.release(&array);
		}
		if (schema.release) {
			schema.release(&schema);
		}
	}
	ArrowExport(const ArrowExport &) = delete;
	ArrowExport &operator=(const ArrowExport &) = delete;
	void Taken() {
		array.release = nullptr;
	}
};

// st3-arrow reads Utf8 / Int32 / Float64 with 32-bit offsets. The session's own
// Arrow settings (string views, list views, large offsets, format version) must
// not leak into what st3 is handed, so they are pinned on a copy.
ClientProperties St3ExportProperties(ClientContext &context) {
	ClientProperties props = context.GetClientProperties();
	props.arrow_offset_size = ArrowOffsetSize::REGULAR;
	props.arrow_use_list_view = false;
	props.produce_arrow_string_view = false;
	props.arrow_output_version = ArrowFormatVersion::V1_0;
	return props;
}

// Export `n` rows of `types` as one Arrow struct array (one child per column)
// through ArrowAppender, fed in vector-sized slices; `fill(chunk, from, count)`
// writes rows [from, from + count) into the chunk's flat vectors.
template <typename Fill>
void ExportColumns(ClientContext &context, const vector<LogicalType> &types, const vector<string> &names, idx_t n,
                   Fill fill, ArrowExport &out) {
	ClientProperties props = St3ExportProperties(context);
	auto extension_types = ArrowTypeExtensionData::GetExtensionTypes(context, types);
	ArrowAppender appender(types, MaxValue<idx_t>(n, 1), props, extension_types);
	DataChunk chunk;
	chunk.Initialize(context, types);
	for (idx_t from = 0; from < n; from += STANDARD_VECTOR_SIZE) {
		const idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, n - from);
		chunk.Reset();
		fill(chunk, from, count);
		chunk.SetCardinality(count);
		appender.Append(chunk, 0, count, count);
	}
	out.array = appender.Finalize();
	ArrowConverter::ToArrowSchema(&out.schema, types, names, props);
}

void SetString(Vector &vec, idx_t i, const std::string &s) {
	FlatVector::GetData<string_t>(vec)[i] = StringVector::AddString(vec, s);
}

// Hand the dataset to st3 and return the table handle.
std::unique_ptr<St3Table, St3TableDeleter> ImportDataset(ClientContext &context, const Dataset &ds) {
	ArrowExport coo;
	ExportColumns(
	    context, {LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::DOUBLE}, {"row", "col", "val"},
	    ds.vals.size(),
	    [&](DataChunk &chunk, idx_t from, idx_t count) {
		    auto rows = FlatVector::GetData<int32_t>(chunk.data[0]);
		    auto cols = FlatVector::GetData<int32_t>(chunk.data[1]);
		    auto vals = FlatVector::GetData<double>(chunk.data[2]);
		    for (idx_t i = 0; i < count; ++i) {
			    rows[i] = ds.rows[from + i];
			    cols[i] = ds.cols[from + i];
			    vals[i] = ds.vals[from + i];
		    }
	    },
	    coo);

	ArrowExport features;
	ExportColumns(
	    context, {LogicalType::VARCHAR}, {"feature_id"}, ds.feature_ids.size(),
	    [&](DataChunk &chunk, idx_t from, idx_t count) {
		    for (idx_t i = 0; i < count; ++i) {
			    SetString(chunk.data[0], i, ds.feature_ids[from + i]);
		    }
	    },
	    features);

	ArrowExport meta;
	ExportColumns(
	    context, {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR}, {"sample_id", "role", "env"},
	    ds.sample_ids.size(),
	    [&](DataChunk &chunk, idx_t from, idx_t count) {
		    for (idx_t i = 0; i < count; ++i) {
			    const idx_t s = from + i;
			    SetString(chunk.data[0], i, ds.sample_ids[s]);
			    SetString(chunk.data[1], i, ds.is_source[s] ? "source" : "sink");
			    if (ds.envs[s].has_value()) {
				    SetString(chunk.data[2], i, *ds.envs[s]);
			    } else {
				    FlatVector::SetNull(chunk.data[2], i, true);
			    }
		    }
	    },
	    meta);

	// st3 wants the feature ids as a plain string array, not a one-column struct.
	// The child is moved out of the root: st3 consumes the copy (children carry
	// their own release callback), and nulling the original's release makes
	// ArrowAppender::ReleaseArray skip it when the root is released below. This
	// happens right before the call so no exception can leave the copy dangling.
	St3Table *raw = nullptr;
	ArrowArray feature_child = *features.array.children[0];
	features.array.children[0]->release = nullptr;
	const St3Status status = st3_table_from_arrow(&coo.array, &coo.schema, &feature_child, features.schema.children[0],
	                                              &meta.array, &meta.schema, &raw);
	// Consumed on every status: the three arrays are st3's now; the schemas and
	// the feature-id root (child already gone) are still ours to release.
	coo.Taken();
	meta.Taken();
	feature_child.release = nullptr;
	if (status != ST3_STATUS_OK) {
		ThrowSt3("st3_table_from_arrow", status);
	}
	return std::unique_ptr<St3Table, St3TableDeleter>(raw);
}

// ---------------------------------------------------------------------------
// Arrow -> DuckDB: st3's dense mixing batches
// ---------------------------------------------------------------------------

struct DenseBatch {
	std::vector<std::string> sink_ids;
	std::vector<std::string> columns;
	std::vector<double> values; // [row][column]
};

[[noreturn]] void BadExport(const char *what, const std::string &why) {
	throw InvalidInputException("sourcetracker: st3 returned an unexpected %s batch: %s", what, why);
}

DenseBatch ReadDenseBatch(const St3Result *result,
                          St3Status (*exporter)(const St3Result *, ArrowArray *, ArrowSchema *), const char *what) {
	ArrowExport out;
	const St3Status status = exporter(result, &out.array, &out.schema);
	if (status != ST3_STATUS_OK) {
		ThrowSt3(what, status);
	}
	const ArrowArray &arr = out.array;
	const ArrowSchema &schema = out.schema;
	if (!schema.format || std::strcmp(schema.format, "+s") != 0 || schema.n_children < 1 ||
	    arr.n_children != schema.n_children) {
		BadExport(what, "not a struct batch with a sink_id column");
	}
	const auto n_rows = static_cast<idx_t>(arr.length);
	const int64_t root_offset = arr.offset;

	DenseBatch batch;
	{
		const ArrowArray &ids = *arr.children[0];
		const ArrowSchema &ids_schema = *schema.children[0];
		const bool large = ids_schema.format && std::strcmp(ids_schema.format, "U") == 0;
		if (!ids_schema.format || (std::strcmp(ids_schema.format, "u") != 0 && !large) || ids.n_buffers != 3 ||
		    ids.null_count != 0) {
			BadExport(what, "sink_id is not a non-null string column");
		}
		const auto data = static_cast<const char *>(ids.buffers[2]);
		batch.sink_ids.reserve(n_rows);
		for (idx_t i = 0; i < n_rows; ++i) {
			const int64_t k = root_offset + ids.offset + static_cast<int64_t>(i);
			int64_t begin, end;
			if (large) {
				const auto offsets = static_cast<const int64_t *>(ids.buffers[1]);
				begin = offsets[k];
				end = offsets[k + 1];
			} else {
				const auto offsets = static_cast<const int32_t *>(ids.buffers[1]);
				begin = offsets[k];
				end = offsets[k + 1];
			}
			batch.sink_ids.emplace_back(data + begin, static_cast<size_t>(end - begin));
		}
	}
	const auto n_cols = static_cast<idx_t>(schema.n_children - 1);
	batch.columns.reserve(n_cols);
	batch.values.assign(n_rows * n_cols, 0.0);
	for (idx_t c = 0; c < n_cols; ++c) {
		const ArrowArray &col = *arr.children[c + 1];
		const ArrowSchema &col_schema = *schema.children[c + 1];
		if (!col_schema.format || std::strcmp(col_schema.format, "g") != 0 || col.n_buffers != 2 ||
		    col.null_count != 0) {
			BadExport(what, "environment columns are not non-null Float64");
		}
		batch.columns.emplace_back(col_schema.name ? col_schema.name : "");
		const auto values = static_cast<const double *>(col.buffers[1]);
		for (idx_t i = 0; i < n_rows; ++i) {
			batch.values[i * n_cols + c] = values[root_offset + col.offset + static_cast<int64_t>(i)];
		}
	}
	return batch;
}

// ---------------------------------------------------------------------------
// Arrow -> DuckDB: st3's per-sink assignment tallies
//
// The stream yields one batch per sink over [sink INT32, source INT32,
// feature INT32, value DOUBLE]: `sink` is the row in the means batch, `source`
// the column, and `feature` the position in the feature-id array st3 was given,
// which is the dataset dictionary. Only nonzero cells are streamed.
// ---------------------------------------------------------------------------

struct ArrowStreamHolder {
	ArrowArrayStream stream {};
	~ArrowStreamHolder() {
		if (stream.release) {
			stream.release(&stream);
		}
	}
	ArrowStreamHolder(const ArrowStreamHolder &) = delete;
	ArrowStreamHolder &operator=(const ArrowStreamHolder &) = delete;
	ArrowStreamHolder() = default;
	[[noreturn]] void Fail(const char *what) {
		// A failed stream call reports through the stream itself, not st3's
		// thread-local slot.
		const char *raw = stream.get_last_error ? stream.get_last_error(&stream) : nullptr;
		throw InvalidInputException("sourcetracker: st3 contingency stream: %s failed: %s", what,
		                            (raw && *raw) ? raw : "no message");
	}
};

template <typename T>
const T *ChildValues(const ArrowArray &batch, idx_t child, const char *format, const ArrowSchema &schema) {
	const ArrowArray &col = *batch.children[child];
	const ArrowSchema &col_schema = *schema.children[child];
	if (!col_schema.format || std::strcmp(col_schema.format, format) != 0 || col.n_buffers != 2 ||
	    col.null_count != 0 || col.length != batch.length) {
		BadExport("contingency", "columns are not non-null (INT32, INT32, INT32, DOUBLE)");
	}
	return static_cast<const T *>(col.buffers[1]) + batch.offset + col.offset;
}

std::vector<std::vector<TallyCell>> ReadContingency(const St3Result *result, idx_t n_sinks, idx_t n_sources,
                                                    idx_t n_features) {
	ArrowStreamHolder holder;
	ArrowArrayStream &stream = holder.stream;
	const St3Status status = st3_result_contingency_stream(result, &stream);
	if (status != ST3_STATUS_OK) {
		ThrowSt3("st3_result_contingency_stream", status);
	}
	ArrowExport schema_holder;
	if (stream.get_schema(&stream, &schema_holder.schema) != 0) {
		holder.Fail("get_schema");
	}
	const ArrowSchema &schema = schema_holder.schema;
	if (!schema.format || std::strcmp(schema.format, "+s") != 0 || schema.n_children != 4) {
		BadExport("contingency", "not a four-column struct batch");
	}

	std::vector<std::vector<TallyCell>> tally(n_sinks * n_sources);
	while (true) {
		ArrowExport batch;
		if (stream.get_next(&stream, &batch.array) != 0) {
			holder.Fail("get_next");
		}
		if (!batch.array.release) {
			break; // end of stream
		}
		const ArrowArray &arr = batch.array;
		if (arr.n_children != 4 || arr.null_count != 0) {
			BadExport("contingency", "not a four-column struct batch");
		}
		const auto sinks = ChildValues<int32_t>(arr, 0, "i", schema);
		const auto sources = ChildValues<int32_t>(arr, 1, "i", schema);
		const auto features = ChildValues<int32_t>(arr, 2, "i", schema);
		const auto values = ChildValues<double>(arr, 3, "g", schema);
		for (int64_t i = 0; i < arr.length; ++i) {
			const int32_t s = sinks[i];
			const int32_t e = sources[i];
			const int32_t f = features[i];
			if (s < 0 || static_cast<idx_t>(s) >= n_sinks || e < 0 || static_cast<idx_t>(e) >= n_sources || f < 0 ||
			    static_cast<idx_t>(f) >= n_features) {
				BadExport("contingency", "an index is outside the sink, source or feature range");
			}
			tally[static_cast<idx_t>(s) * n_sources + static_cast<idx_t>(e)].emplace_back(f, values[i]);
		}
	}
	// A MAP needs each key once; the tally is one cell per (source, feature) by
	// construction, and this makes that a checked fact rather than an assumption.
	for (auto &cells : tally) {
		std::sort(cells.begin(), cells.end(), [](const TallyCell &a, const TallyCell &b) { return a.first < b.first; });
		const auto dup = std::adjacent_find(cells.begin(), cells.end(),
		                                    [](const TallyCell &a, const TallyCell &b) { return a.first == b.first; });
		if (dup != cells.end()) {
			BadExport("contingency", "a feature is tallied twice for one (sink, source)");
		}
	}
	return tally;
}

// ---------------------------------------------------------------------------
// InitGlobal: read both relations once, validate, run st3, keep the result.
// ---------------------------------------------------------------------------

unique_ptr<GlobalTableFunctionState> SourcetrackerInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<SourcetrackerBindData>();
	auto gstate = make_uniq<SourcetrackerGlobalState>();
	gstate->sink_id_type = data.sink_id_type;
	gstate->feature_id_type = data.feature_id_type;
	gstate->assignments = data.assignments;

	auto cells = ReadFeatureTable(context, data.table_name, kCaller);
	auto metadata = ReadWideMetadata(context, data.metadata_name, {"source_sink", "env"}, kCaller);

	Dataset ds;
	try {
		ds = miint::sourcetracker::IngestDataset(cells, metadata.rows);
		if (!data.loo && std::find(ds.is_source.begin(), ds.is_source.end(), false) == ds.is_source.end()) {
			throw std::invalid_argument("no sink samples in the metadata; sink mode needs at least one sample with "
			                            "source_sink = 'sink' (use loo := true to predict source samples instead)");
		}
		miint::sourcetracker::CheckDepths(ds, data.config.source_rarefaction_depth, data.config.sink_rarefaction_depth,
		                                  data.loo);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("sourcetracker: %s", e.what());
	}
	CheckInterrupt(context);

	auto table = ImportDataset(context, ds);

	St3Config config = data.config;
	config.seed = data.seed >= 0 ? static_cast<uint64_t>(data.seed) : std::random_device {}();
	St3Result *result_raw = nullptr;
	// One blocking call; Ctrl-C is honoured only between phases (st3 has no
	// cancellation hook yet).
	const St3Status run_status = st3_run(table.get(), &config, &result_raw);
	if (run_status != ST3_STATUS_OK) {
		ThrowSt3("st3_run", run_status);
	}
	std::unique_ptr<St3Result, St3ResultDeleter> result(result_raw);
	table.reset();
	CheckInterrupt(context);

	DenseBatch means = ReadDenseBatch(result.get(), st3_result_means, "st3_result_means");
	DenseBatch stds = ReadDenseBatch(result.get(), st3_result_stds, "st3_result_stds");
	if (stds.sink_ids != means.sink_ids || stds.columns != means.columns) {
		BadExport("st3_result_stds", "shape differs from the means batch");
	}
	gstate->sink_ids = std::move(means.sink_ids);
	gstate->sources = std::move(means.columns);
	gstate->means = std::move(means.values);
	gstate->stds = std::move(stds.values);
	if (data.assignments) {
		gstate->tally =
		    ReadContingency(result.get(), gstate->sink_ids.size(), gstate->sources.size(), ds.feature_ids.size());
		gstate->feature_ids = std::move(ds.feature_ids);
	}
	return std::move(gstate);
}

// ---------------------------------------------------------------------------
// Execute: one row per (sink, source). The assignments MAP is NULL when the
// tally was not requested, otherwise one entry per nonzero feature (so a
// source with no mass in a sink has an empty map).
// ---------------------------------------------------------------------------

void SourcetrackerExecute(ClientContext &, TableFunctionInput &input, DataChunk &output) {
	auto &g = input.global_state->Cast<SourcetrackerGlobalState>();
	const idx_t n_sources = g.sources.size();
	const idx_t total = g.sink_ids.size() * n_sources;
	if (g.cursor >= total) {
		output.SetCardinality(0);
		return;
	}
	const idx_t n = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - g.cursor);

	auto &sink_vec = output.data[0];
	auto &source_vec = output.data[1];
	auto proportion = FlatVector::GetData<double>(output.data[2]);
	auto proportion_std = FlatVector::GetData<double>(output.data[3]);
	auto &map_vec = output.data[4];
	auto map_entries = FlatVector::GetData<list_entry_t>(map_vec);

	for (idx_t i = 0; i < n; ++i) {
		const idx_t k = g.cursor + i;
		EmitIdCell(sink_vec, i, g.sink_ids[k / n_sources], g.sink_id_type);
		SetString(source_vec, i, g.sources[k % n_sources]);
		proportion[i] = g.means[k];
		proportion_std[i] = g.stds[k];
	}

	if (!g.assignments) {
		for (idx_t i = 0; i < n; ++i) {
			map_entries[i] = list_entry_t {0, 0};
			FlatVector::SetNull(map_vec, i, true);
		}
		ListVector::SetListSize(map_vec, 0);
	} else {
		idx_t total_entries = 0;
		for (idx_t i = 0; i < n; ++i) {
			total_entries += g.tally[g.cursor + i].size();
		}
		auto &keys = MapVector::GetKeys(map_vec);
		auto &values = MapVector::GetValues(map_vec);
		idx_t offset = ListVector::GetListSize(map_vec);
		ListVector::Reserve(map_vec, offset + total_entries);
		for (idx_t i = 0; i < n; ++i) {
			const auto &cells = g.tally[g.cursor + i];
			map_entries[i] = list_entry_t {offset, cells.size()};
			for (const auto &cell : cells) {
				// Fetched per cell: Reserve above may have moved the child buffers.
				EmitIdCell(keys, offset, g.feature_ids[cell.first], g.feature_id_type);
				FlatVector::GetData<double>(values)[offset] = cell.second;
				++offset;
			}
		}
		ListVector::SetListSize(map_vec, offset);
	}
	g.cursor += n;
	output.SetCardinality(n);
}

} // namespace

void RegisterSourcetracker(ExtensionLoader &loader) {
	TableFunction fn("sourcetracker", {LogicalType::VARCHAR, LogicalType::VARCHAR}, SourcetrackerExecute,
	                 SourcetrackerBind, SourcetrackerInitGlobal);
	fn.named_parameters["loo"] = LogicalType::BOOLEAN;
	fn.named_parameters["assignments"] = LogicalType::BOOLEAN;
	fn.named_parameters["alpha1"] = LogicalType::DOUBLE;
	fn.named_parameters["alpha2"] = LogicalType::DOUBLE;
	fn.named_parameters["beta"] = LogicalType::DOUBLE;
	fn.named_parameters["restarts"] = LogicalType::INTEGER;
	fn.named_parameters["draws_per_restart"] = LogicalType::INTEGER;
	fn.named_parameters["burnin"] = LogicalType::INTEGER;
	fn.named_parameters["delay"] = LogicalType::INTEGER;
	fn.named_parameters["source_rarefaction_depth"] = LogicalType::INTEGER;
	fn.named_parameters["sink_rarefaction_depth"] = LogicalType::INTEGER;
	fn.named_parameters["with_replacement"] = LogicalType::BOOLEAN;
	fn.named_parameters["collapse"] = LogicalType::VARCHAR;
	fn.named_parameters["seed"] = LogicalType::BIGINT;
	fn.named_parameters["threads"] = LogicalType::INTEGER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
