#include "phylogeny_fasttree.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/arrow/arrow.hpp"
#include "duckdb/common/arrow/arrow_wrapper.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/materialized_query_result.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include "catalog_utils.hpp"
#include "gpl_boundary/arrow_ipc.hpp"
#include "gpl_boundary/process.hpp"
#include "gpl_boundary/session.hpp"

#include "nanoarrow/nanoarrow.h"

#include <cstring>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h> // read(2) on macOS isn't reachable via <cstdio>; need this for ::read
#include <unordered_set>
#include <vector>

namespace duckdb {

namespace {

namespace gb = ::duckdb::miint::gpl_boundary;

// =============================================================================
// Bind: parameter validation + JSON config builder
// =============================================================================
//
// gpl-boundary commit 19306f6's `apply_json_to_config` function in
// `src/tools/fasttree.rs` defines the wire-level config schema for FastTree.
// We mirror its 22 wired knobs verbatim and re-implement its 9 mutual-exclusion
// checks at miint's Bind time so users see config errors at SQL compile time
// rather than after a daemon round-trip.
//
// `fastest=true` is rejected unconditionally (Linus-flagged behavioral
// difference between the C library API and the FastTree CLI).

struct ConfigJsonBuilder {
	// Append `"key":value,` to the JSON object. Caller is responsible for
	// stripping the trailing comma and wrapping in `{...}`.
	std::ostringstream os;
	bool first = true;

	void append_raw(const std::string &key, const std::string &raw_value) {
		if (!first) {
			os << ",";
		}
		first = false;
		os << "\"" << key << "\":" << raw_value;
	}
	void append_str(const std::string &key, const std::string &value) {
		append_raw(key, "\"" + value + "\"");
	}
	void append_int(const std::string &key, int64_t value) {
		append_raw(key, std::to_string(value));
	}
	void append_double(const std::string &key, double value) {
		// max_digits10 (17 for double) round-trips exactly; std::to_string
		// would truncate to 6 fractional digits and silently lose precision
		// on parameters like topm or GTR rates.
		std::ostringstream s;
		s << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
		append_raw(key, s.str());
	}
	void append_bool(const std::string &key, bool value) {
		append_raw(key, value ? "true" : "false");
	}
	std::string build() const {
		return "{" + os.str() + "}";
	}
};

struct PhylogenyFastTreeBindData : public TableFunctionData {
	std::string table_name;
	std::string config_json;
};

// Helper: throw if a named parameter wasn't recognized by FastTree's wire
// schema. Catches typos at Bind time.
const std::unordered_set<std::string> kKnownParams = {
    "seq_type", "seed",  "verbose", "bootstrap", "nosupport", "pseudo", "pseudo_weight", "nni",     "spr",
    "mlnni",    "mlacc", "cat",     "noml",      "threads",   "model",  "gtrrates",      "gtrfreq", "slow",
    "bionj",    "nj",    "top",     "notop",     "topm",      "quote",  "fastest",       "gamma"};

// Convert a DuckDB Value containing an integer-like type to int64. Throws
// InvalidInputException if the value isn't integer-typed.
int64_t value_to_int(const std::string &param, const Value &v) {
	auto t = v.type().id();
	if (t == LogicalTypeId::INTEGER || t == LogicalTypeId::BIGINT || t == LogicalTypeId::HUGEINT ||
	    t == LogicalTypeId::SMALLINT || t == LogicalTypeId::TINYINT || t == LogicalTypeId::UINTEGER ||
	    t == LogicalTypeId::UBIGINT || t == LogicalTypeId::USMALLINT || t == LogicalTypeId::UTINYINT) {
		return v.GetValue<int64_t>();
	}
	throw InvalidInputException("phylogeny_fasttree: parameter '%s' expects an integer, got %s", param,
	                            v.type().ToString());
}

double value_to_double(const std::string &param, const Value &v) {
	auto t = v.type().id();
	if (t == LogicalTypeId::DOUBLE || t == LogicalTypeId::FLOAT || t == LogicalTypeId::DECIMAL) {
		return v.GetValue<double>();
	}
	if (t == LogicalTypeId::INTEGER || t == LogicalTypeId::BIGINT || t == LogicalTypeId::SMALLINT ||
	    t == LogicalTypeId::TINYINT) {
		return static_cast<double>(v.GetValue<int64_t>());
	}
	throw InvalidInputException("phylogeny_fasttree: parameter '%s' expects a number, got %s", param,
	                            v.type().ToString());
}

bool value_to_bool(const std::string &param, const Value &v) {
	if (v.type().id() == LogicalTypeId::BOOLEAN) {
		return v.GetValue<bool>();
	}
	throw InvalidInputException("phylogeny_fasttree: parameter '%s' expects a boolean, got %s", param,
	                            v.type().ToString());
}

std::string value_to_str(const std::string &param, const Value &v) {
	if (v.type().id() == LogicalTypeId::VARCHAR) {
		return v.GetValue<std::string>();
	}
	throw InvalidInputException("phylogeny_fasttree: parameter '%s' expects a string, got %s", param,
	                            v.type().ToString());
}

// Build the JSON config from named_parameters, applying every check that
// gpl-boundary's `apply_json_to_config` performs. Failures throw before the
// daemon spawns. Returns the serialized JSON object string ("{...}").
std::string BuildConfigJson(const named_parameter_map_t &named_params) {
	ConfigJsonBuilder cfg;

	// First pass: reject unknown params.
	for (const auto &kv : named_params) {
		if (kKnownParams.find(kv.first) == kKnownParams.end()) {
			throw InvalidInputException("phylogeny_fasttree: unknown named parameter '%s'. "
			                            "See `docs/phylogeny.md` for the supported list.",
			                            kv.first);
		}
	}

	auto get = [&](const std::string &k) -> const Value * {
		auto it = named_params.find(k);
		return it == named_params.end() ? nullptr : &it->second;
	};

	// === Pre-resolve fields that participate in cross-knob checks. ===
	std::string seq_type = "auto";
	if (auto *v = get("seq_type")) {
		seq_type = value_to_str("seq_type", *v);
		if (seq_type != "auto" && seq_type != "nucleotide" && seq_type != "protein") {
			throw InvalidInputException("phylogeny_fasttree: invalid seq_type '%s' (expected auto/nucleotide/protein)",
			                            seq_type);
		}
		cfg.append_str("seq_type", seq_type);
	}

	if (auto *v = get("seed")) {
		cfg.append_int("seed", value_to_int("seed", *v));
	}

	if (auto *v = get("verbose")) {
		cfg.append_bool("verbose", value_to_bool("verbose", *v));
	}

	// bootstrap + nosupport
	const Value *boot_v = get("bootstrap");
	const Value *nosup_v = get("nosupport");
	bool nosupport = nosup_v && value_to_bool("nosupport", *nosup_v);
	if (boot_v) {
		const int64_t b = value_to_int("bootstrap", *boot_v);
		if (b < 0) {
			throw InvalidInputException("phylogeny_fasttree: bootstrap must be >= 0 (got %lld)",
			                            static_cast<long long>(b));
		}
		if (nosupport && b != 0) {
			throw InvalidInputException("phylogeny_fasttree: bootstrap=%lld cannot be combined with nosupport=true "
			                            "(use bootstrap=0 or omit nosupport)",
			                            static_cast<long long>(b));
		}
		cfg.append_int("bootstrap", b);
	}
	if (nosup_v) {
		cfg.append_bool("nosupport", nosupport);
	}

	// pseudo + pseudo_weight
	const Value *pseudo_v = get("pseudo");
	const Value *pweight_v = get("pseudo_weight");
	if (pweight_v) {
		if (!pseudo_v || !value_to_bool("pseudo", *pseudo_v)) {
			throw InvalidInputException("phylogeny_fasttree: pseudo_weight requires pseudo=true "
			                            "(set pseudo=true to enable, or omit pseudo_weight)");
		}
		const double w = value_to_double("pseudo_weight", *pweight_v);
		if (w < 0.0) {
			throw InvalidInputException("phylogeny_fasttree: pseudo_weight must be >= 0 (got %f)", w);
		}
		cfg.append_double("pseudo_weight", w);
	}
	if (pseudo_v) {
		cfg.append_bool("pseudo", value_to_bool("pseudo", *pseudo_v));
	}

	// nni
	if (auto *v = get("nni")) {
		const int64_t n = value_to_int("nni", *v);
		if (n < 0) {
			throw InvalidInputException("phylogeny_fasttree: nni must be >= 0 (got %lld); "
			                            "omit the parameter to use FastTree's auto default",
			                            static_cast<long long>(n));
		}
		cfg.append_int("nni", n);
	}

	// spr
	if (auto *v = get("spr")) {
		const int64_t n = value_to_int("spr", *v);
		if (n < 0) {
			throw InvalidInputException("phylogeny_fasttree: spr must be >= 0 (got %lld)", static_cast<long long>(n));
		}
		cfg.append_int("spr", n);
	}

	// mlnni + noml
	const Value *mlnni_v = get("mlnni");
	const Value *noml_v = get("noml");
	const bool noml = noml_v && value_to_bool("noml", *noml_v);
	if (mlnni_v) {
		const int64_t n = value_to_int("mlnni", *mlnni_v);
		if (n < 0) {
			throw InvalidInputException("phylogeny_fasttree: mlnni must be >= 0 (got %lld); "
			                            "omit the parameter to use FastTree's auto default",
			                            static_cast<long long>(n));
		}
		if (noml && n > 0) {
			throw InvalidInputException("phylogeny_fasttree: mlnni=%lld cannot be combined with noml=true "
			                            "(use mlnni=0 or omit noml)",
			                            static_cast<long long>(n));
		}
		cfg.append_int("mlnni", n);
	}
	if (noml_v) {
		cfg.append_bool("noml", noml);
	}

	// mlacc
	if (auto *v = get("mlacc")) {
		const int64_t n = value_to_int("mlacc", *v);
		if (n < 1) {
			throw InvalidInputException("phylogeny_fasttree: mlacc must be >= 1 (got %lld)", static_cast<long long>(n));
		}
		cfg.append_int("mlacc", n);
	}

	// cat
	if (auto *v = get("cat")) {
		const int64_t n = value_to_int("cat", *v);
		if (n < 1) {
			throw InvalidInputException("phylogeny_fasttree: cat must be >= 1 (got %lld)", static_cast<long long>(n));
		}
		cfg.append_int("cat", n);
	}

	// threads
	if (auto *v = get("threads")) {
		const int64_t n = value_to_int("threads", *v);
		if (n < 1) {
			throw InvalidInputException("phylogeny_fasttree: threads must be >= 1 (got %lld)",
			                            static_cast<long long>(n));
		}
		cfg.append_int("threads", n);
	}

	// model + GTR cross-checks
	std::string model = "auto";
	const Value *model_v = get("model");
	if (model_v) {
		model = value_to_str("model", *model_v);
		const bool ok =
		    model == "auto" || model == "jtt" || model == "lg" || model == "wag" || model == "jc" || model == "gtr";
		if (!ok) {
			throw InvalidInputException("phylogeny_fasttree: invalid model '%s' "
			                            "(expected one of auto/jtt/lg/wag/jc/gtr)",
			                            model);
		}
		const bool is_protein_model = (model == "jtt" || model == "lg" || model == "wag");
		const bool is_nucleotide_model = (model == "jc" || model == "gtr");
		if (is_protein_model && seq_type == "nucleotide") {
			throw InvalidInputException("phylogeny_fasttree: model='%s' is a protein model and "
			                            "cannot be combined with seq_type='nucleotide'",
			                            model);
		}
		if (is_nucleotide_model && seq_type == "protein") {
			throw InvalidInputException("phylogeny_fasttree: model='%s' is a nucleotide model and "
			                            "cannot be combined with seq_type='protein'",
			                            model);
		}
		cfg.append_str("model", model);
	}

	// gtrrates / gtrfreq — only valid with model=gtr
	auto build_array_of_numbers = [&](const std::string &key, const Value &v, size_t expected_len) -> std::string {
		if (v.type().id() != LogicalTypeId::LIST) {
			throw InvalidInputException("phylogeny_fasttree: %s must be a list of numbers", key);
		}
		auto &children = ListValue::GetChildren(v);
		if (children.size() != expected_len) {
			throw InvalidInputException("phylogeny_fasttree: %s must have exactly %d entries (got %d)", key,
			                            static_cast<int>(expected_len), static_cast<int>(children.size()));
		}
		std::ostringstream arr;
		arr << "[";
		arr << std::setprecision(std::numeric_limits<double>::max_digits10);
		for (size_t i = 0; i < children.size(); ++i) {
			const double f = value_to_double(key, children[i]);
			if (f < 0.0) {
				throw InvalidInputException("phylogeny_fasttree: %s[%d] must be >= 0 (got %f)", key,
				                            static_cast<int>(i), f);
			}
			if (i > 0) {
				arr << ",";
			}
			arr << f;
		}
		arr << "]";
		return arr.str();
	};

	const Value *gtrrates_v = get("gtrrates");
	const Value *gtrfreq_v = get("gtrfreq");
	if ((gtrrates_v || gtrfreq_v) && model != "gtr") {
		throw InvalidInputException("phylogeny_fasttree: gtrrates/gtrfreq are only valid with model='gtr'");
	}
	if (gtrrates_v) {
		cfg.append_raw("gtrrates", build_array_of_numbers("gtrrates", *gtrrates_v, 6));
	}
	if (gtrfreq_v) {
		cfg.append_raw("gtrfreq", build_array_of_numbers("gtrfreq", *gtrfreq_v, 4));
	}

	// slow
	const Value *slow_v = get("slow");
	const bool slow_set = slow_v && value_to_bool("slow", *slow_v);
	if (slow_v) {
		cfg.append_bool("slow", slow_set);
	}

	// bionj + nj
	const Value *bionj_v = get("bionj");
	const Value *nj_v = get("nj");
	const bool bionj = bionj_v && value_to_bool("bionj", *bionj_v);
	const bool nj = nj_v && value_to_bool("nj", *nj_v);
	if (bionj && nj) {
		throw InvalidInputException("phylogeny_fasttree: bionj=true and nj=true are mutually exclusive "
		                            "(BIONJ and NJ are different join algorithms)");
	}
	if (bionj_v) {
		cfg.append_bool("bionj", bionj);
	}
	if (nj_v) {
		cfg.append_bool("nj", nj);
	}

	// top + notop + topm
	const Value *top_v = get("top");
	const Value *notop_v = get("notop");
	const Value *topm_v = get("topm");
	const bool top_set = top_v && value_to_bool("top", *top_v);
	const bool notop_set = notop_v && value_to_bool("notop", *notop_v);
	if (top_set && notop_set) {
		throw InvalidInputException("phylogeny_fasttree: top=true and notop=true are mutually exclusive");
	}
	if (slow_set && top_set) {
		throw InvalidInputException("phylogeny_fasttree: slow=true disables top hits; cannot also set top=true");
	}
	if ((notop_set || slow_set) && topm_v) {
		throw InvalidInputException("phylogeny_fasttree: topm is only meaningful when top hits are enabled "
		                            "(top hits are off when notop=true or slow=true)");
	}
	if (top_v) {
		cfg.append_bool("top", top_set);
	}
	if (notop_v) {
		cfg.append_bool("notop", notop_set);
	}
	if (topm_v) {
		const double m = value_to_double("topm", *topm_v);
		if (m <= 0.0) {
			throw InvalidInputException("phylogeny_fasttree: topm must be > 0 (got %f)", m);
		}
		cfg.append_double("topm", m);
	}

	// quote
	if (auto *v = get("quote")) {
		cfg.append_bool("quote", value_to_bool("quote", *v));
	}

	// fastest — Bind-side rejection (matches the daemon's rejection but
	// surfaces at SQL-compile time, before any process is spawned).
	// `fastest=false` is intentionally NOT emitted: the daemon's `fastest`
	// implementation is a strict subset of the CLI's, and explicitly emitting
	// `false` could activate future "I know about this knob" behavior we
	// don't want. Omitting the field preserves the daemon's default.
	if (auto *v = get("fastest")) {
		const bool fastest = value_to_bool("fastest", *v);
		if (fastest) {
			throw InvalidInputException("phylogeny_fasttree: fastest=true is not yet supported. "
			                            "The ext/fasttree library API exposes only a subset of the CLI's "
			                            "-fastest flag (missing tophitsRefresh and useTopHits2nd); wiring "
			                            "fastest=1 alone would silently produce a tree distinct from "
			                            "FastTree -fastest. Omit this parameter until the upstream API "
			                            "catches up.");
		}
		// fastest=false: accepted as a no-op. Do not emit anything.
	}

	// gamma
	if (auto *v = get("gamma")) {
		cfg.append_bool("gamma", value_to_bool("gamma", *v));
	}

	return cfg.build();
}

// =============================================================================
// Bind
// =============================================================================
unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input, vector<LogicalType> &return_types,
                              vector<string> &names) {
	if (input.inputs.empty() || input.inputs[0].IsNull()) {
		throw InvalidInputException("phylogeny_fasttree: a non-NULL input table name is required");
	}
	const std::string table_name = input.inputs[0].GetValue<std::string>();

	auto bind_data = make_uniq<PhylogenyFastTreeBindData>();
	bind_data->table_name = table_name;
	bind_data->config_json = BuildConfigJson(input.named_parameters);

	// Output schema — matches gpl-boundary's `--describe fasttree`
	// `output_schema` (commit 19306f6).
	names = {"node_index", "parent_index", "edge_id", "branch_length", "support", "is_tip", "name", "n_children"};
	return_types = {LogicalType::BIGINT, LogicalType::BIGINT,  LogicalType::BIGINT,  LogicalType::DOUBLE,
	                LogicalType::DOUBLE, LogicalType::BOOLEAN, LogicalType::VARCHAR, LogicalType::BIGINT};
	return std::move(bind_data);
}

// =============================================================================
// Global state
// =============================================================================

// Holds everything that needs to live across Execute calls. Lifetime order
// matters: the `ArrowArrayWrapper`s reference into the IPC bytes that
// `submit_result` owns via `OutputShmRegion`. Declaring `submit_result` first
// guarantees it's destroyed last (member destructors run in reverse).
struct PhylogenyFastTreeGlobalState : public GlobalTableFunctionState {
	// Always populated by InitGlobal; held as unique_ptr only because
	// gb::Session is move-only and we want `PhylogenyFastTreeGlobalState`
	// to be default-constructible before InitGlobal builds it up.
	std::unique_ptr<gb::Session> session;
	gb::SubmitResult submit_result;
	std::unique_ptr<gb::IpcStreamDecoder> decoder;

	// Decoded schema + arrays from the IPC stream. ArrowArrayWrapper releases
	// on destruction.
	ArrowSchemaWrapper schema;
	std::vector<ArrowArrayWrapper> batches;

	// Iteration cursor.
	idx_t batch_index = 0;
	idx_t row_in_batch = 0;

	idx_t MaxThreads() const override {
		return 1; // single-threaded drain
	}
};

struct PhylogenyFastTreeLocalState : public LocalTableFunctionState {};

// Read the input table into two parallel vectors (names, sequences). Throws
// on missing / wrong-shaped columns or empty input.
struct LoadedInput {
	std::vector<std::string> names;
	std::vector<std::string> sequences;
};

LoadedInput LoadInputTable(ClientContext &context, const std::string &table_name) {
	auto conn = MakeReadOnlyHelperConnection(context);
	const std::string sql = "SELECT name, sequence FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);
	auto result = conn.Query(sql);
	if (result->HasError()) {
		throw InvalidInputException("phylogeny_fasttree: failed to read input table '%s' "
		                            "(must have columns name VARCHAR, sequence VARCHAR): %s",
		                            table_name, result->GetError());
	}

	LoadedInput out;
	auto &materialized = result->Cast<MaterializedQueryResult>();
	out.names.reserve(materialized.RowCount());
	out.sequences.reserve(materialized.RowCount());
	while (auto chunk = materialized.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); i++) {
			auto name_val = chunk->GetValue(0, i);
			auto seq_val = chunk->GetValue(1, i);
			if (name_val.IsNull() || seq_val.IsNull()) {
				throw InvalidInputException("phylogeny_fasttree: NULL name or sequence in input table '%s'",
				                            table_name);
			}
			out.names.push_back(name_val.GetValue<std::string>());
			out.sequences.push_back(seq_val.GetValue<std::string>());
		}
	}
	if (out.names.empty()) {
		throw InvalidInputException("phylogeny_fasttree: input table '%s' is empty", table_name);
	}
	if (out.names.size() < 3) {
		throw InvalidInputException("phylogeny_fasttree: input table '%s' has only %d row(s); "
		                            "FastTree requires at least 3 sequences to produce a meaningful tree",
		                            table_name, static_cast<int>(out.names.size()));
	}
	return out;
}

// Build an Arrow IPC stream from a (names, sequences) pair matching the
// gpl-boundary fasttree input_schema:
//   { name: utf8 not null, sequence: utf8 not null }
std::vector<uint8_t> BuildInputIpcStream(const std::vector<std::string> &names,
                                         const std::vector<std::string> &sequences) {
	ArrowSchema schema {};
	auto rc = ArrowSchemaInitFromType(&schema, NANOARROW_TYPE_STRUCT);
	if (rc != NANOARROW_OK) {
		throw InternalException("phylogeny_fasttree: ArrowSchemaInit failed");
	}
	rc = ArrowSchemaAllocateChildren(&schema, 2);
	if (rc != NANOARROW_OK) {
		schema.release(&schema);
		throw InternalException("phylogeny_fasttree: ArrowSchemaAllocateChildren failed");
	}
	ArrowSchemaInitFromType(schema.children[0], NANOARROW_TYPE_STRING);
	ArrowSchemaSetName(schema.children[0], "name");
	ArrowSchemaInitFromType(schema.children[1], NANOARROW_TYPE_STRING);
	ArrowSchemaSetName(schema.children[1], "sequence");

	ArrowArray array {};
	ArrowError err {};
	if (ArrowArrayInitFromSchema(&array, &schema, &err) != NANOARROW_OK) {
		schema.release(&schema);
		throw InternalException("phylogeny_fasttree: ArrowArrayInit failed: %s", err.message);
	}
	if (ArrowArrayStartAppending(&array) != NANOARROW_OK) {
		array.release(&array);
		schema.release(&schema);
		throw InternalException("phylogeny_fasttree: ArrowArrayStartAppending failed");
	}
	for (size_t i = 0; i < names.size(); ++i) {
		ArrowStringView nv {names[i].data(), static_cast<int64_t>(names[i].size())};
		ArrowStringView sv {sequences[i].data(), static_cast<int64_t>(sequences[i].size())};
		if (ArrowArrayAppendString(array.children[0], nv) != NANOARROW_OK ||
		    ArrowArrayAppendString(array.children[1], sv) != NANOARROW_OK ||
		    ArrowArrayFinishElement(&array) != NANOARROW_OK) {
			array.release(&array);
			schema.release(&schema);
			throw InternalException("phylogeny_fasttree: array append failed at row %d", static_cast<int>(i));
		}
	}
	if (ArrowArrayFinishBuildingDefault(&array, &err) != NANOARROW_OK) {
		array.release(&array);
		schema.release(&schema);
		throw InternalException("phylogeny_fasttree: ArrowArrayFinishBuilding failed: %s", err.message);
	}

	std::vector<uint8_t> bytes;
	try {
		bytes = gb::EncodeIpcStream(&schema, &array, 1);
	} catch (...) {
		array.release(&array);
		schema.release(&schema);
		throw;
	}
	if (array.release) {
		array.release(&array);
	}
	if (schema.release) {
		schema.release(&schema);
	}
	return bytes;
}

unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<PhylogenyFastTreeBindData>();
	auto state = make_uniq<PhylogenyFastTreeGlobalState>();

	// 1. Load the input alignment.
	auto loaded = LoadInputTable(context, bind_data.table_name);

	// 2. Encode it as Arrow IPC.
	const auto ipc_bytes = BuildInputIpcStream(loaded.names, loaded.sequences);

	// 3. Spawn the gpl-boundary daemon. Fail loud (with a hint) when the
	//    binary isn't found in any of FindGplBoundary's lookup locations
	//    (MIINT_GPL_BOUNDARY_PATH, miint's install cache, PATH).
	const std::string gpl_path = gb::FindGplBoundary();
	if (gpl_path.empty()) {
		throw IOException("phylogeny_fasttree: gpl-boundary binary not found. To install:\n"
		                  "  Easiest:  SELECT install_gpl_boundary();   "
		                  "-- downloads a prebuilt binary into miint's cache dir\n"
		                  "  Manual:   curl -fsSL "
		                  "https://github.com/the-miint/GPL-boundary/releases/latest/download/install.sh | sh\n"
		                  "  Source:   https://github.com/the-miint/GPL-boundary#building-from-source\n"
		                  "Currently supported prebuilt platforms: Linux x86_64, macOS arm64. macOS Intel users "
		                  "must build from source.\n"
		                  "If gpl-boundary is installed at a non-standard location, set "
		                  "MIINT_GPL_BOUNDARY_PATH=<absolute path>.");
	}
	std::vector<std::string> argv = {gpl_path}; // no args → daemon mode
	gb::ChildProcess child(argv);
	state->session = std::make_unique<gb::Session>(std::move(child));
	state->session->Initialize();

	// 4. Submit the FastTree batch. Daemon returns one shm_outputs entry
	//    (the SOA tree as Arrow IPC bytes).
	state->submit_result =
	    state->session->Submit("fasttree", bind_data.config_json, ipc_bytes.data(), ipc_bytes.size());
	if (state->submit_result.outputs.empty()) {
		throw IOException("phylogeny_fasttree: daemon returned a successful response with zero "
		                  "shm_outputs (this is a daemon bug — please file an issue against "
		                  "gpl-boundary)");
	}

	// 5. Decode the IPC stream from the first output. The OutputShmRegion in
	//    submit_result.outputs[0].shm owns the mmap; we read through it.
	const auto &out0 = state->submit_result.outputs[0];
	state->decoder = std::make_unique<gb::IpcStreamDecoder>(out0.bytes(), out0.size_bytes());
	state->decoder->GetSchema(&state->schema.arrow_schema);

	// Drain all batches up front. fasttree sends a single batch in practice,
	// but the API allows multiple; pull them all so we can iterate strictly.
	for (;;) {
		ArrowArrayWrapper w;
		if (!state->decoder->NextBatch(&w.arrow_array)) {
			break;
		}
		state->batches.push_back(std::move(w));
	}

	if (state->batches.empty()) {
		throw IOException("phylogeny_fasttree: daemon's IPC stream contained zero record batches");
	}

	// 6. Sanity-check the output schema. If the daemon bumps schema_version
	//    in a way that adds/removes/reorders/renames columns, fail loud here
	//    rather than reading the wrong column with the wrong type in Execute.
	if (state->schema.arrow_schema.n_children != 8) {
		throw IOException("phylogeny_fasttree: daemon returned unexpected schema "
		                  "(%lld columns, expected 8). gpl-boundary version drift?",
		                  static_cast<long long>(state->schema.arrow_schema.n_children));
	}
	static const char *const kExpectedNames[8] = {"node_index", "parent_index", "edge_id", "branch_length",
	                                              "support",    "n_children",   "is_tip",  "name"};
	for (int c = 0; c < 8; ++c) {
		const auto *child = state->schema.arrow_schema.children[c];
		if (!child || !child->name || std::strcmp(child->name, kExpectedNames[c]) != 0) {
			throw IOException("phylogeny_fasttree: schema drift at column %d — expected '%s', "
			                  "got '%s'. gpl-boundary version mismatch?",
			                  c, kExpectedNames[c], (child && child->name) ? child->name : "(null)");
		}
	}

	return std::move(state);
}

unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                              GlobalTableFunctionState *gstate) {
	(void)context;
	(void)input;
	(void)gstate;
	return make_uniq<PhylogenyFastTreeLocalState>();
}

// =============================================================================
// Execute: drain decoded ArrowArrays into the output DataChunk
// =============================================================================
//
// Output columns (must match Bind's order):
//   0: node_index    BIGINT  ← Int64
//   1: parent_index  BIGINT  ← Int64 (nullable)
//   2: edge_id       BIGINT  ← Int64 (nullable)
//   3: branch_length DOUBLE  ← Float64 (nullable)
//   4: support       DOUBLE  ← Float64 (nullable)
//   5: is_tip        BOOLEAN ← Boolean
//   6: name          VARCHAR ← Utf8 (nullable)
//   7: n_children    BIGINT  ← Int64

// Read one row at logical index `i` from a fixed-width child column. Honors
// the child array's own .offset (Arrow C Data Interface allows children to
// carry independent offsets in addition to the parent's offset).
template <typename T>
T read_fixed(const ArrowArray &col, idx_t logical_index) {
	const auto *buf = reinterpret_cast<const T *>(col.buffers[1]);
	return buf[static_cast<idx_t>(col.offset) + logical_index];
}

bool read_bool_bit(const ArrowArray &col, idx_t logical_index) {
	const auto *bitmap = static_cast<const uint8_t *>(col.buffers[1]);
	const idx_t a = static_cast<idx_t>(col.offset) + logical_index;
	return (bitmap[a / 8] & (1u << (a % 8))) != 0;
}

bool is_null_at(const ArrowArray &arr, idx_t logical_index) {
	// Arrow C Data Interface (columnar-0.4): `null_count == -1` means
	// "unknown — consult the bitmap". The IPC decoder may leave it at -1.
	// Treating -1 as "no nulls" is a silent-wrong-NULLs bug.
	if (!arr.buffers[0]) {
		return false; // no validity buffer ⇒ no nulls
	}
	if (arr.null_count == 0) {
		return false;
	}
	// null_count > 0 OR null_count == -1: consult bitmap
	const auto *bitmap = static_cast<const uint8_t *>(arr.buffers[0]);
	const idx_t abs = static_cast<idx_t>(arr.offset) + logical_index;
	return (bitmap[abs / 8] & (1u << (abs % 8))) == 0;
}

void Execute(ClientContext &context, TableFunctionInput &data, DataChunk &output) {
	auto &gstate = data.global_state->Cast<PhylogenyFastTreeGlobalState>();
	(void)context;

	if (gstate.batch_index >= gstate.batches.size()) {
		output.SetCardinality(0);
		return;
	}

	auto &batch = gstate.batches[gstate.batch_index].arrow_array;

	// nanoarrow's IPC decoder produces fresh, unsliced batches with
	// batch.offset == 0. We rely on that — if a future decoder change starts
	// producing sliced parents, we'd need to add `batch.offset` to every
	// child index below. Fail loud if the assumption breaks.
	if (batch.offset != 0) {
		throw IOException("phylogeny_fasttree: decoded batch has non-zero parent offset (%lld); "
		                  "Execute does not yet handle sliced struct arrays",
		                  static_cast<long long>(batch.offset));
	}

	const idx_t total = static_cast<idx_t>(batch.length);
	const idx_t remaining = total - gstate.row_in_batch;
	const idx_t to_emit = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);
	output.SetCardinality(to_emit);

	// gpl-boundary `19306f6` `output_schema` (verified at
	// `ext/GPL-boundary/src/tools/fasttree.rs:145-152`):
	//   wire[0] node_index    Int64   not null
	//   wire[1] parent_index  Int64   nullable
	//   wire[2] edge_id       Int64   nullable
	//   wire[3] branch_length Float64 nullable
	//   wire[4] support       Float64 nullable
	//   wire[5] n_children    Int32   not null   ← narrow, widen to Int64 for SQL
	//   wire[6] is_tip        Boolean not null
	//   wire[7] name          Utf8    nullable
	//
	// miint exposes the user-facing columns in a different order (is_tip and
	// name before n_children) for readability. Map wire→output here.
	const auto &col_node_idx = *batch.children[0];
	const auto &col_parent_idx = *batch.children[1];
	const auto &col_edge_id = *batch.children[2];
	const auto &col_branch_len = *batch.children[3];
	const auto &col_support = *batch.children[4];
	const auto &col_n_children = *batch.children[5];
	const auto &col_is_tip = *batch.children[6];
	const auto &col_name = *batch.children[7];

	auto *out_node_idx = FlatVector::GetData<int64_t>(output.data[0]);
	auto *out_parent_idx = FlatVector::GetData<int64_t>(output.data[1]);
	auto *out_edge_id = FlatVector::GetData<int64_t>(output.data[2]);
	auto *out_branch_len = FlatVector::GetData<double>(output.data[3]);
	auto *out_support = FlatVector::GetData<double>(output.data[4]);
	auto *out_is_tip = FlatVector::GetData<bool>(output.data[5]);
	auto *out_n_children = FlatVector::GetData<int64_t>(output.data[7]);

	auto &v_name = output.data[6];

	auto &mask_node_idx = FlatVector::Validity(output.data[0]);
	auto &mask_parent_idx = FlatVector::Validity(output.data[1]);
	auto &mask_edge_id = FlatVector::Validity(output.data[2]);
	auto &mask_branch_len = FlatVector::Validity(output.data[3]);
	auto &mask_support = FlatVector::Validity(output.data[4]);
	auto &mask_name = FlatVector::Validity(v_name);
	auto &mask_n_children = FlatVector::Validity(output.data[7]);

	// For variable-length Utf8: offsets is `length + 1` int32s, indexed at
	// the column's own offset. data is a contiguous byte buffer addressed by
	// those offset values.
	const auto *p_name_offsets = static_cast<const int32_t *>(col_name.buffers[1]);
	const auto *p_name_data = static_cast<const char *>(col_name.buffers[2]);

	for (idx_t i = 0; i < to_emit; ++i) {
		const idx_t li = gstate.row_in_batch + i; // logical row index within child columns

		if (is_null_at(col_node_idx, li)) {
			mask_node_idx.SetInvalid(i);
		} else {
			out_node_idx[i] = read_fixed<int64_t>(col_node_idx, li);
		}
		if (is_null_at(col_parent_idx, li)) {
			mask_parent_idx.SetInvalid(i);
		} else {
			out_parent_idx[i] = read_fixed<int64_t>(col_parent_idx, li);
		}
		if (is_null_at(col_edge_id, li)) {
			mask_edge_id.SetInvalid(i);
		} else {
			out_edge_id[i] = read_fixed<int64_t>(col_edge_id, li);
		}
		if (is_null_at(col_branch_len, li)) {
			mask_branch_len.SetInvalid(i);
		} else {
			out_branch_len[i] = read_fixed<double>(col_branch_len, li);
		}
		if (is_null_at(col_support, li)) {
			mask_support.SetInvalid(i);
		} else {
			out_support[i] = read_fixed<double>(col_support, li);
		}

		out_is_tip[i] = read_bool_bit(col_is_tip, li);

		if (is_null_at(col_name, li) || !p_name_offsets || !p_name_data) {
			mask_name.SetInvalid(i);
		} else {
			const idx_t a = static_cast<idx_t>(col_name.offset) + li;
			const int32_t start = p_name_offsets[a];
			const int32_t end = p_name_offsets[a + 1];
			const int32_t len = end - start;
			if (len < 0) {
				throw IOException("phylogeny_fasttree: corrupt utf8 offsets at row %lld (start=%d end=%d)",
				                  static_cast<long long>(a), start, end);
			}
			FlatVector::GetData<string_t>(v_name)[i] =
			    StringVector::AddString(v_name, p_name_data + start, static_cast<idx_t>(len));
		}

		// n_children is Int32 on the wire; widen to Int64 for the user-facing
		// BIGINT column. n_children should never be null (gpl-boundary marks
		// it not-null in the schema), but defend anyway.
		if (is_null_at(col_n_children, li)) {
			mask_n_children.SetInvalid(i);
		} else {
			out_n_children[i] = static_cast<int64_t>(read_fixed<int32_t>(col_n_children, li));
		}
	}

	gstate.row_in_batch += to_emit;
	if (gstate.row_in_batch >= total) {
		gstate.batch_index += 1;
		gstate.row_in_batch = 0;
	}
}

// =============================================================================
// phylogeny_fasttree_available() scalar
// =============================================================================
//
// Mirrors `bowtie2_available()` (`src/align_bowtie2.cpp:205-260`):
// `std::call_once`-cached PATH lookup so the cost is paid at most once per
// process.

void PhylogenyFastTreeAvailableImpl(DataChunk &args, ExpressionState &state, Vector &result) {
	(void)args;
	(void)state;
	static std::once_flag flag;
	static bool available = false;
	std::call_once(flag, [&]() {
		const std::string path = gb::FindGplBoundary();
		if (path.empty()) {
			available = false;
			return;
		}
		// Confirm the binary actually responds to --list-tools and reports fasttree.
		try {
			std::vector<std::string> argv = {path, "--list-tools"};
			gb::ChildProcess child(argv);
			std::string out;
			char buf[256];
			ssize_t n;
			while ((n = ::read(child.stdout_fd(), buf, sizeof(buf))) > 0) {
				out.append(buf, static_cast<size_t>(n));
			}
			const int status = child.Wait();
			available = WIFEXITED(status) && WEXITSTATUS(status) == 0 && out.find("fasttree") != std::string::npos;
		} catch (...) {
			available = false;
		}
	});

	result.SetVectorType(VectorType::CONSTANT_VECTOR);
	auto &constant_validity = ConstantVector::Validity(result);
	constant_validity.SetAllValid(1);
	*ConstantVector::GetData<bool>(result) = available;
}

} // namespace

// =============================================================================
// Public registration
// =============================================================================

TableFunction PhylogenyFastTreeTableFunction::GetFunction() {
	TableFunction fn("phylogeny_fasttree", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);
	fn.named_parameters["seq_type"] = LogicalType::VARCHAR;
	fn.named_parameters["seed"] = LogicalType::BIGINT;
	fn.named_parameters["verbose"] = LogicalType::BOOLEAN;
	fn.named_parameters["bootstrap"] = LogicalType::BIGINT;
	fn.named_parameters["nosupport"] = LogicalType::BOOLEAN;
	fn.named_parameters["pseudo"] = LogicalType::BOOLEAN;
	fn.named_parameters["pseudo_weight"] = LogicalType::DOUBLE;
	fn.named_parameters["nni"] = LogicalType::BIGINT;
	fn.named_parameters["spr"] = LogicalType::BIGINT;
	fn.named_parameters["mlnni"] = LogicalType::BIGINT;
	fn.named_parameters["mlacc"] = LogicalType::BIGINT;
	fn.named_parameters["cat"] = LogicalType::BIGINT;
	fn.named_parameters["noml"] = LogicalType::BOOLEAN;
	fn.named_parameters["threads"] = LogicalType::BIGINT;
	fn.named_parameters["model"] = LogicalType::VARCHAR;
	fn.named_parameters["gtrrates"] = LogicalType::LIST(LogicalType::DOUBLE);
	fn.named_parameters["gtrfreq"] = LogicalType::LIST(LogicalType::DOUBLE);
	fn.named_parameters["slow"] = LogicalType::BOOLEAN;
	fn.named_parameters["bionj"] = LogicalType::BOOLEAN;
	fn.named_parameters["nj"] = LogicalType::BOOLEAN;
	fn.named_parameters["top"] = LogicalType::BOOLEAN;
	fn.named_parameters["notop"] = LogicalType::BOOLEAN;
	fn.named_parameters["topm"] = LogicalType::DOUBLE;
	fn.named_parameters["quote"] = LogicalType::BOOLEAN;
	fn.named_parameters["fastest"] = LogicalType::BOOLEAN;
	fn.named_parameters["gamma"] = LogicalType::BOOLEAN;
	fn.order_preservation_type = OrderPreservationType::NO_ORDER;
	return fn;
}

void PhylogenyFastTreeTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

ScalarFunction PhylogenyFastTreeAvailableScalar::GetFunction() {
	return ScalarFunction("phylogeny_fasttree_available", {}, LogicalType::BOOLEAN, PhylogenyFastTreeAvailableImpl);
}

void PhylogenyFastTreeAvailableScalar::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
