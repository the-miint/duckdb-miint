#include "skb_arrow_functions.hpp"

// Before the headers that declare duckdb::miint, which would shadow the
// global ::miint this one names unqualified.
#include "unifrac_function_common.hpp"

#include "catalog_utils.hpp"
#include "gpl_boundary/arrow_ipc.hpp"
#include "miint_arrow_stream.hpp"
#include "miint_log.hpp"
#include "miint_named_parameter.hpp"
#include "skb_arrow/client.hpp"

#include "duckdb/common/arrow/arrow_wrapper.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/column/column_data_collection.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"
#include "duckdb/function/table/arrow.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

namespace duckdb {

namespace {

namespace sk = ::duckdb::miint::skb_arrow;
namespace gb = ::duckdb::miint::gpl_boundary;
using unifrac_internal::ProbeDistanceTableIdType;
using unifrac_internal::ProbeFeatureTableIdType;
using unifrac_internal::ResolveSampleIdOutputType;
using unifrac_internal::ResolveThreadsParameter;

// Rows per Arrow record batch sent to skb-arrow, and so between interrupt checks.
constexpr idx_t kRowsPerBatch = 65536;

// =============================================================================
// One call: bind data, and the bridge between DuckDB relations and skb-arrow
// =============================================================================

// An input table: skb-arrow's name for it, the user's relation, and the query
// that reads the relation in the shape skb-arrow takes.
struct SkbArrowInput {
	std::string table;
	std::string relation;
	std::string sql;

	bool operator==(const SkbArrowInput &other) const {
		return table == other.table && relation == other.relation && sql == other.sql;
	}
};

// Everything a call needs, resolved at bind; nothing is read or run until
// InitGlobal. A prepared statement re-runs the call on each EXECUTE, and
// EXPLAIN never starts a host.
struct SkbArrowCallData : public TableFunctionData {
	std::string function;
	std::string capability;
	std::vector<SkbArrowInput> inputs;
	// The named params the user set, forwarded verbatim; skb-arrow fills the
	// rest with its defaults and checks every value.
	std::vector<std::pair<std::string, Value>> params;
	int threads = 1;
	vector<string> output_names;
	// The output's types as skb-arrow sends them, and as the function returns
	// them; they differ where an ID returns to the input's type.
	vector<LogicalType> wire_types;
	vector<LogicalType> output_types;

	unique_ptr<FunctionData> Copy() const override {
		return make_uniq<SkbArrowCallData>(*this);
	}

	bool Equals(const FunctionData &other_p) const override {
		const auto &other = other_p.Cast<SkbArrowCallData>();
		if (function != other.function || capability != other.capability || inputs != other.inputs ||
		    threads != other.threads || output_names != other.output_names || wire_types != other.wire_types ||
		    output_types != other.output_types || params.size() != other.params.size()) {
			return false;
		}
		for (idx_t i = 0; i < params.size(); i++) {
			if (params[i].first != other.params[i].first ||
			    !Value::NotDistinctFrom(params[i].second, other.params[i].second)) {
				return false;
			}
		}
		return true;
	}
};

struct SkbArrowGlobalState : public GlobalTableFunctionState {
	unique_ptr<ColumnDataCollection> rows;
	ColumnDataScanState scan;

	idx_t MaxThreads() const override {
		return 1;
	}
};

// skb-arrow's failures as DuckDB exceptions, prefixed with the function's
// name. Never InternalException: on DuckDB v2 one invalidates the whole
// database, and the host's bug is not the database's.
[[noreturn]] void ThrowAsDuckDB(const std::string &function, const sk::Error &e) {
	const std::string &kind = e.kind();
	const std::string message = function + ": " + e.what();
	if (kind == "invalid_input" || kind == "invalid_param") {
		throw InvalidInputException(message);
	}
	if (kind == "unsupported") {
		throw NotImplementedException(message);
	}
	if (kind == "resource") {
		if (message.find("MIINT_SKB_ARROW_SEGMENT_DIR") != std::string::npos) {
			throw IOException(message);
		}
		throw IOException(message + ". skb-arrow's data passes through " + sk::SegmentDirParent() +
		                  "; set MIINT_SKB_ARROW_SEGMENT_DIR to a directory with more room (in Docker, /dev/shm "
		                  "is 64 MiB unless the container runs with --shm-size)");
	}
	if (kind == "host_incompatible") {
		throw IOException(message);
	}
	std::string internal = function + ": skb-arrow internal error (please report): " + e.what();
	if (!e.traceback().empty()) {
		internal += "\n" + e.traceback();
	}
	throw IOException(internal);
}

std::string LocateSkbArrow(const std::string &function) {
	const std::string binary = sk::FindSkbArrow();
	if (binary.empty() || ::access(binary.c_str(), X_OK) != 0) {
		throw IOException(function + ": skb-arrow not found" + (binary.empty() ? "" : " at " + binary) +
		                  ". Install it with: SELECT install_skb_arrow();  (it needs uv). If skb-arrow is "
		                  "installed elsewhere, set MIINT_SKB_ARROW_PATH=<absolute path>.");
	}
	return binary;
}

sk::Params ToParams(const std::vector<std::pair<std::string, Value>> &params) {
	sk::Params out;
	for (const auto &param : params) {
		const std::string &name = param.first;
		const Value &value = param.second;
		switch (value.type().id()) {
		case LogicalTypeId::BOOLEAN:
			out.SetBool(name, BooleanValue::Get(value));
			break;
		case LogicalTypeId::BIGINT:
			out.SetInt(name, BigIntValue::Get(value));
			break;
		case LogicalTypeId::DOUBLE:
			out.SetNumber(name, DoubleValue::Get(value));
			break;
		case LogicalTypeId::VARCHAR:
			out.SetString(name, StringValue::Get(value));
			break;
		case LogicalTypeId::LIST: {
			std::vector<std::string> items;
			for (const auto &item : ListValue::GetChildren(value)) {
				items.push_back(StringValue::Get(item));
			}
			out.SetStrings(name, std::move(items));
			break;
		}
		default:
			throw NotImplementedException("skb-arrow param '%s' of type %s", name, value.type().ToString());
		}
	}
	return out;
}

// Stream `input`'s query into segments in `dir`. One relation at a time: a
// connection streams one result.
std::vector<std::string> WriteInput(ClientContext &context, Connection &conn, const std::string &function,
                                    const SkbArrowInput &input, const std::string &dir) {
	auto query = SubmitStream(conn, input.sql);
	if (query->HasError()) {
		throw InvalidInputException(function + ": cannot read '" + input.relation + "': " + query->GetError());
	}
	// What nanoarrow's IPC writer and skb-arrow both take: Arrow 1.0 layouts
	// with 32-bit offsets. The helper connection inherits the user's settings,
	// so pin every one that changes the export.
	ClientProperties properties = query->GetClientProperties();
	properties.arrow_offset_size = ArrowOffsetSize::REGULAR;
	properties.arrow_output_version = ArrowFormatVersion::V1_0;
	properties.produce_arrow_string_view = false;
	properties.arrow_use_list_view = false;
	properties.arrow_lossless_conversion = false;
	// Releasing the stream deletes its ChunkSourceArrowStream, and the query with it.
	ArrowArrayStreamWrapper stream;
	stream.arrow_array_stream =
	    (new ChunkSourceArrowStream<StreamingQuery>(std::move(query), kRowsPerBatch, std::move(properties)))->stream;
	ArrowArrayStream &s = stream.arrow_array_stream;
	const auto read_error = [&]() {
		return InvalidInputException(function + ": cannot read '" + input.relation + "': " + s.get_last_error(&s));
	};

	ArrowSchemaWrapper schema;
	if (s.get_schema(&s, &schema.arrow_schema) != 0) {
		throw read_error();
	}
	sk::SegmentWriter writer(dir, input.table, schema.arrow_schema);
	while (true) {
		context.BlockingInterruptCheck();
		ArrowArrayWrapper batch;
		if (s.get_next(&s, &batch.arrow_array) != 0) {
			throw read_error();
		}
		if (!batch.arrow_array.release) {
			break;
		}
		writer.Append(batch.arrow_array);
	}
	return writer.Finish();
}

std::string DescribeColumns(const vector<string> &names, const vector<LogicalType> &types) {
	std::string out;
	for (idx_t i = 0; i < names.size(); i++) {
		out += (i ? ", " : "") + names[i] + " " + types[i].ToString();
	}
	return out;
}

// Decode one output segment's batches into `out`, after checking its columns
// are the ones `data` declares.
void ReadOutputSegment(ClientContext &context, const SkbArrowCallData &data, const sk::MappedSegment &segment,
                       ColumnDataCollection &out) {
	gb::IpcStreamDecoder decoder(segment.data(), segment.size());
	ArrowSchemaWrapper schema;
	decoder.GetSchema(&schema.arrow_schema);
	ArrowTableSchema table;
	ArrowTableFunction::PopulateArrowTableSchema(context, table, schema.arrow_schema);
	if (table.GetNames() != data.output_names || table.GetTypes() != data.wire_types) {
		throw sk::Error("host_incompatible", data.capability + " returned (" +
		                                         DescribeColumns(table.GetNames(), table.GetTypes()) + "), not (" +
		                                         DescribeColumns(data.output_names, data.wire_types) + ")");
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
			DataChunk chunk;
			chunk.Initialize(context, data.wire_types);
			chunk.SetChildCardinality(n);
			// As DuckDB's own Arrow scan (ArrowTableFunction::ArrowToDuckDB) converts a column.
			for (idx_t c = 0; c < data.output_types.size(); c++) {
				auto &array = *batch->arrow_array.children[c];
				const ArrowType &arrow_type = *columns.at(c);
				ArrowArrayScanState state(context);
				state.owned_data = batch;
				switch (arrow_type.GetPhysicalType()) {
				case ArrowArrayPhysicalType::DICTIONARY_ENCODED:
					ArrowToDuckDBConversion::ColumnArrowToDuckDBDictionary(chunk.data[c], array, offset, state, n,
					                                                       arrow_type);
					break;
				case ArrowArrayPhysicalType::RUN_END_ENCODED:
					ArrowToDuckDBConversion::ColumnArrowToDuckDBRunEndEncoded(chunk.data[c], array, offset, state, n,
					                                                          arrow_type);
					break;
				default:
					ArrowToDuckDBConversion::SetValidityMask(chunk.data[c], array, offset, n, batch->arrow_array.offset,
					                                         -1);
					ArrowToDuckDBConversion::ColumnArrowToDuckDB(chunk.data[c], array, offset, state, n, arrow_type);
					break;
				}
				FlatVector::SetSize(chunk.data[c], count_t(n));
			}
			if (data.wire_types == data.output_types) {
				out.Append(chunk);
				continue;
			}
			DataChunk cast;
			cast.Initialize(context, data.output_types);
			cast.SetChildCardinality(n);
			for (idx_t c = 0; c < data.output_types.size(); c++) {
				if (data.wire_types[c] == data.output_types[c]) {
					cast.data[c].Reference(chunk.data[c]);
				} else {
					VectorOperations::Cast(context, chunk.data[c], cast.data[c], n);
				}
			}
			out.Append(cast);
		}
	}
}

// Start a host, send it the inputs, make the call, relay its warnings, and
// decode its output into `out`.
void RunCall(ClientContext &context, const SkbArrowCallData &data, ColumnDataCollection &out) {
	sk::HostOptions options;
	options.binary = LocateSkbArrow(data.function);
	options.segment_parent = sk::SegmentDirParent();
	options.threads = data.threads;
	// The schema_version miint was built against, from the table Probe checks.
	for (const auto &capability : sk::MiintCapabilities()) {
		if (capability.first == data.capability) {
			options.required = {capability};
		}
	}
	if (options.required.empty()) {
		ThrowAsDuckDB(data.function,
		              sk::Error("internal", "miint has no schema_version for capability '" + data.capability + "'"));
	}
	const sk::InterruptCheck check = [&context] {
		context.BlockingInterruptCheck();
	};
	try {
		sk::Host host(options, check);
		auto conn = MakeReadOnlyHelperConnection(context);
		std::vector<std::pair<std::string, std::vector<std::string>>> tables;
		tables.reserve(data.inputs.size());
		for (const auto &input : data.inputs) {
			tables.emplace_back(input.table, WriteInput(context, conn, data.function, input, host.dir()));
		}
		const sk::CallResult result = host.Call(data.capability, ToParams(data.params), tables, check);
		for (const auto &warning : result.warnings) {
			std::string message = data.function + ": " + warning.category + ": " + warning.message;
			if (warning.count > 1) {
				message += " (x" + std::to_string(warning.count) + ")";
			}
			::miint::EmitWarning(context, message);
		}
		for (const auto &name : result.output) {
			const sk::MappedSegment segment = host.MapOutput(name);
			try {
				ReadOutputSegment(context, data, segment, out);
			} catch (const sk::Error &) {
				throw;
			} catch (const Exception &) {
				throw; // DuckDB's own (e.g. out of memory): not the host's doing
			} catch (const std::runtime_error &e) {
				// The IPC decoder's: the bytes aren't the Arrow stream the contract promises.
				throw sk::Error("host_incompatible", "cannot decode output segment " + name + ": " + e.what());
			}
		}
		// The answer is in, and the host is stateless: waiting for Python's
		// teardown would add tens of ms to every call for nothing.
		host.Kill();
	} catch (const sk::Error &e) {
		ThrowAsDuckDB(data.function, e);
	}
}

unique_ptr<GlobalTableFunctionState> SkbArrowInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	const auto &data = input.bind_data->Cast<SkbArrowCallData>();
	auto state = make_uniq<SkbArrowGlobalState>();
	state->rows = make_uniq<ColumnDataCollection>(context, data.output_types);
	RunCall(context, data, *state->rows);
	state->rows->InitializeScan(state->scan);
	return std::move(state);
}

void SkbArrowExecute(ClientContext &, TableFunctionInput &input, DataChunk &output) {
	auto &state = input.global_state->Cast<SkbArrowGlobalState>();
	state.rows->Scan(state.scan, output);
}

// Named params: `threads` is miint's; every other one is forwarded verbatim,
// in name order. A NULL one is unset, as skb-arrow reads null: its default.
void ReadNamedParams(ClientContext &context, TableFunctionBindInput &input, SkbArrowCallData &data) {
	int32_t threads = 0;
	for (const auto &kv : input.named_parameters) {
		const auto name = StringUtil::Lower(kv.first.GetIdentifierName());
		const Value &value = kv.second;
		if (value.IsNull()) {
			continue;
		}
		if (name == "threads") {
			threads = value.GetValue<int32_t>();
			continue;
		}
		// NaN and Infinity aren't JSON, and a null would mean "the default".
		if (value.type().id() == LogicalTypeId::DOUBLE && !std::isfinite(DoubleValue::Get(value))) {
			throw BinderException(data.function + ": " + name + " must be finite (got " + value.ToString() + ")");
		}
		if (value.type().id() == LogicalTypeId::LIST) {
			for (const auto &item : ListValue::GetChildren(value)) {
				if (item.IsNull()) {
					throw BinderException(data.function + ": " + name + " must not contain NULL");
				}
			}
		}
		data.params.emplace_back(name, value);
	}
	std::sort(data.params.begin(), data.params.end(),
	          [](const std::pair<std::string, Value> &a, const std::pair<std::string, Value> &b) {
		          return a.first < b.first;
	          });
	data.threads = ResolveThreadsParameter(context, threads, data.function);
}

// =============================================================================
// Inputs: miint's relations as skb-arrow's tables
// =============================================================================

// skb-arrow takes integer or string IDs. Integers of up to 64 bits go as
// BIGINT (a UBIGINT beyond BIGINT is refused: see IdExpression); anything
// else (UUID, HUGEINT, ...) as text.
LogicalType IdWireType(const LogicalType &type) {
	switch (type.id()) {
	case LogicalTypeId::TINYINT:
	case LogicalTypeId::SMALLINT:
	case LogicalTypeId::INTEGER:
	case LogicalTypeId::BIGINT:
	case LogicalTypeId::UTINYINT:
	case LogicalTypeId::USMALLINT:
	case LogicalTypeId::UINTEGER:
	case LogicalTypeId::UBIGINT:
		return LogicalType::BIGINT;
	default:
		return LogicalType::VARCHAR;
	}
}

// The SQL that reads ID column `column` (of `type`) as IdWireType sends it.
// A UBIGINT beyond BIGINT fails the read by name, never wrapping to a
// negative ID.
std::string IdExpression(const std::string &column, const LogicalType &type) {
	if (type.id() == LogicalTypeId::UBIGINT) {
		return "CASE WHEN " + column + " > 9223372036854775807::UBIGINT THEN error('" + column + " ' || " + column +
		       "::VARCHAR || ' is beyond BIGINT, the widest integer ID skb-arrow takes; cast the IDs to VARCHAR') "
		       "ELSE " +
		       column + "::BIGINT END";
	}
	return column + "::" + IdWireType(type).ToString();
}

// miint's condensed (sample_a, sample_b, distance) relation, as skb-arrow's
// distance table (id_a, id_b, distance DOUBLE). The cast absorbs DECIMAL and
// HUGEINT distances, which skb-arrow refuses.
SkbArrowInput DistanceInput(ClientContext &context, const std::string &function, const std::string &table,
                            const std::string &relation) {
	// Columns, castability, and one output id type across sample_a and sample_b.
	LogicalType sample_a_type;
	LogicalType sample_b_type;
	ProbeDistanceTableIdType(context, relation, function, nullptr, &sample_a_type, &sample_b_type);
	// Each column goes as its own kind of ID, so they must be one kind: an
	// INTEGER beside a VARCHAR would otherwise read '007' as 7.
	if (IdWireType(sample_a_type) != IdWireType(sample_b_type)) {
		throw BinderException(function + ": distance-table '" + relation + "' has sample_a " +
		                      sample_a_type.ToString() + " and sample_b " + sample_b_type.ToString() +
		                      ", which skb-arrow takes as different kinds of ID (integer and text); cast both to "
		                      "one type in a view");
	}
	return {table, relation,
	        "SELECT " + IdExpression("sample_a", sample_a_type) + " AS id_a, " +
	            IdExpression("sample_b", sample_b_type) + " AS id_b, distance::DOUBLE AS distance FROM " +
	            SQLIdentifier::ToString(relation)};
}

// miint's long (sample_id, feature_id, value) feature table, as skb-arrow's.
// `feature_id_type` receives the feature_id column's own type.
SkbArrowInput FeatureTableInput(ClientContext &context, const std::string &function, const std::string &relation,
                                LogicalType &feature_id_type) {
	// Columns and castability.
	LogicalType sample_id_type;
	ProbeFeatureTableIdType(context, relation, function, &sample_id_type, &feature_id_type);
	return {"table", relation,
	        "SELECT " + IdExpression("sample_id", sample_id_type) + " AS sample_id, " +
	            IdExpression("feature_id", feature_id_type) + " AS feature_id, value::DOUBLE AS value FROM " +
	            SQLIdentifier::ToString(relation)};
}

// The type a metadata covariate is sent as. skb-arrow takes booleans, integers,
// floating point and strings; DECIMAL and HUGEINT go as DOUBLE. ENUM is
// refused: nanoarrow's IPC writer can't send a dictionary, and text would
// quietly move the reference level from the ENUM's first to the alphabetical first.
std::string CovariateWireType(const std::string &function, const std::string &column, const LogicalType &type) {
	switch (type.id()) {
	case LogicalTypeId::BOOLEAN:
	case LogicalTypeId::TINYINT:
	case LogicalTypeId::SMALLINT:
	case LogicalTypeId::INTEGER:
	case LogicalTypeId::BIGINT:
	case LogicalTypeId::UTINYINT:
	case LogicalTypeId::USMALLINT:
	case LogicalTypeId::UINTEGER:
	case LogicalTypeId::UBIGINT:
	case LogicalTypeId::FLOAT:
	case LogicalTypeId::DOUBLE:
	case LogicalTypeId::VARCHAR:
		// Its own type, without any alias or collation.
		return LogicalType(type.id()).ToString();
	case LogicalTypeId::DECIMAL:
	case LogicalTypeId::HUGEINT:
	case LogicalTypeId::UHUGEINT:
		return "DOUBLE";
	case LogicalTypeId::ENUM:
		throw BinderException(function + ": metadata column '" + column +
		                      "' is an ENUM, which skb-arrow can't take; cast it to VARCHAR in a view (" + column +
		                      "::VARCHAR AS " + column +
		                      "). Its reference level then becomes the alphabetically first, not the ENUM's first");
	default:
		throw BinderException(function + ": metadata column '" + column + "' is " + type.ToString() +
		                      "; a covariate must be BOOLEAN, an integer, FLOAT, DOUBLE, DECIMAL or VARCHAR");
	}
}

// The names in `formula`: the metadata columns it can mean. skb-arrow parses
// the formula itself; this only picks which columns to send, so a column the
// formula doesn't name never matters, whatever its type.
std::set<std::string> FormulaNames(const std::string &formula) {
	std::set<std::string> names;
	std::string name;
	for (const char c : formula) {
		const auto u = static_cast<unsigned char>(c);
		// Python identifiers: ASCII letters, digits and '_', and any non-ASCII.
		if ((u >= 'a' && u <= 'z') || (u >= 'A' && u <= 'Z') || (u >= '0' && u <= '9') || u == '_' || u >= 0x80) {
			name += c;
		} else if (!name.empty()) {
			names.insert(name);
			name.clear();
		}
	}
	if (!name.empty()) {
		names.insert(name);
	}
	return names;
}

// miint's wide metadata (sample_id, then a column per covariate), as
// skb-arrow's: sample_id plus the columns `formula` or `grouping` name.
SkbArrowInput MetadataInput(ClientContext &context, const std::string &function, const std::string &relation,
                            const std::string &formula, const std::string &grouping) {
	auto wanted = FormulaNames(formula);
	if (!grouping.empty()) {
		wanted.insert(grouping);
	}
	const auto columns = GetTableOrViewColumns(context, relation, function + ": metadata table");
	const LogicalType *sample_id_type = nullptr;
	std::string covariates;
	for (idx_t i = 0; i < columns.names.size(); i++) {
		const std::string &column = columns.names[i];
		if (StringUtil::Lower(column) == "sample_id") {
			sample_id_type = &columns.types[i];
		} else if (wanted.count(column)) {
			const auto quoted = SQLIdentifier::ToString(column);
			covariates +=
			    ", " + quoted + "::" + CovariateWireType(function, column, columns.types[i]) + " AS " + quoted;
		}
	}
	if (!sample_id_type) {
		throw BinderException(function + ": metadata '" + relation + "' has no sample_id column");
	}
	return {"metadata", relation,
	        "SELECT " + IdExpression("sample_id", *sample_id_type) + " AS sample_id" + covariates + " FROM " +
	            SQLIdentifier::ToString(relation)};
}

// =============================================================================
// skb_mantel
// =============================================================================

unique_ptr<FunctionData> SkbMantelBind(ClientContext &context, TableFunctionBindInput &input,
                                       vector<LogicalType> &return_types, vector<Identifier> &names) {
	auto data = make_uniq<SkbArrowCallData>();
	data->function = "skb_mantel";
	data->capability = "mantel";
	const std::string x = input.inputs[0].GetValue<string>();
	const std::string y = input.inputs[1].GetValue<string>();
	// Every name before any is resolved: a CTE named by y must not hide behind a missing x.
	RejectCTERelationName(input, x);
	RejectCTERelationName(input, y);
	data->inputs.push_back(DistanceInput(context, data->function, "x", x));
	data->inputs.push_back(DistanceInput(context, data->function, "y", y));
	ReadNamedParams(context, input, *data);
	data->output_names = {"statistic", "pvalue", "n"};
	data->output_types = {LogicalType::DOUBLE, LogicalType::DOUBLE, LogicalType::BIGINT};
	data->wire_types = data->output_types;
	return_types = data->output_types;
	for (const auto &name : data->output_names) {
		names.emplace_back(name);
	}
	return std::move(data);
}

// =============================================================================
// skb_ancombc
// =============================================================================

unique_ptr<FunctionData> SkbAncombcBind(ClientContext &context, TableFunctionBindInput &input,
                                        vector<LogicalType> &return_types, vector<Identifier> &names) {
	auto data = make_uniq<SkbArrowCallData>();
	data->function = "skb_ancombc";
	data->capability = "ancombc";
	const std::string table = input.inputs[0].GetValue<string>();
	const std::string metadata = input.inputs[1].GetValue<string>();
	const std::string formula = input.inputs[2].GetValue<string>();
	RejectCTERelationName(input, table);
	RejectCTERelationName(input, metadata);
	data->params.emplace_back("formula", Value(formula));
	ReadNamedParams(context, input, *data);
	std::string grouping;
	for (const auto &param : data->params) {
		if (param.first == "grouping") {
			grouping = StringValue::Get(param.second);
		}
	}
	LogicalType feature_id_type = LogicalType::VARCHAR;
	data->inputs.push_back(FeatureTableInput(context, data->function, table, feature_id_type));
	data->inputs.push_back(MetadataInput(context, data->function, metadata, formula, grouping));

	data->output_names = {"feature_id", "test", "term", "lfc", "se", "w", "pvalue", "qvalue", "signif"};
	data->wire_types = {IdWireType(feature_id_type), LogicalType::VARCHAR, LogicalType::VARCHAR,
	                    LogicalType::DOUBLE,         LogicalType::DOUBLE,  LogicalType::DOUBLE,
	                    LogicalType::DOUBLE,         LogicalType::DOUBLE,  LogicalType::BOOLEAN};
	// feature_id returns as miint's id convention: BIGINT and UUID mirrored,
	// anything else VARCHAR.
	data->output_types = data->wire_types;
	data->output_types[0] = ResolveSampleIdOutputType(feature_id_type);
	return_types = data->output_types;
	for (const auto &name : data->output_names) {
		names.emplace_back(name);
	}
	return std::move(data);
}

} // namespace

void RegisterSkbArrowFunctions(ExtensionLoader &loader) {
	TableFunction mantel("skb_mantel", {LogicalType::VARCHAR, LogicalType::VARCHAR}, SkbArrowExecute, SkbMantelBind,
	                     SkbArrowInitGlobal);
	AddNamedParameter(mantel, "method", LogicalType::VARCHAR);
	AddNamedParameter(mantel, "permutations", LogicalType::BIGINT);
	AddNamedParameter(mantel, "alternative", LogicalType::VARCHAR);
	AddNamedParameter(mantel, "seed", LogicalType::BIGINT);
	AddNamedParameter(mantel, "threads", LogicalType::INTEGER);
	loader.RegisterFunction(mantel);

	TableFunction ancombc("skb_ancombc", {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR},
	                      SkbArrowExecute, SkbAncombcBind, SkbArrowInitGlobal);
	AddNamedParameter(ancombc, "grouping", LogicalType::VARCHAR);
	AddNamedParameter(ancombc, "posthoc", LogicalType::LIST(LogicalType::VARCHAR));
	AddNamedParameter(ancombc, "pseudocount", LogicalType::DOUBLE);
	AddNamedParameter(ancombc, "max_iter", LogicalType::BIGINT);
	AddNamedParameter(ancombc, "tol", LogicalType::DOUBLE);
	AddNamedParameter(ancombc, "alpha", LogicalType::DOUBLE);
	AddNamedParameter(ancombc, "p_adjust", LogicalType::VARCHAR);
	AddNamedParameter(ancombc, "bootstraps", LogicalType::BIGINT);
	AddNamedParameter(ancombc, "seed", LogicalType::BIGINT);
	AddNamedParameter(ancombc, "threads", LogicalType::INTEGER);
	loader.RegisterFunction(ancombc);
}

} // namespace duckdb
