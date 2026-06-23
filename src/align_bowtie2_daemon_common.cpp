#include "align_bowtie2_daemon_common.hpp"

#include "align_common.hpp"
#include "fastq_encoder.hpp"
#include "gpl_boundary/arrow_ipc.hpp"
#include "gpl_boundary/session.hpp"
// This translation unit references both `::miint::` (codec, pulled in via
// align_common.hpp → id_column_utils.hpp → id_column_codec.hpp) and
// `duckdb::miint::gpl_boundary` (via the `gb` alias below). Direct codec
// calls inside the `duckdb` namespace must use the absolute `::miint::`
// qualifier to avoid resolving to duckdb::miint first.

#include "gpl_boundary/process.hpp"
#include "nanoarrow/nanoarrow.h"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/materialized_query_result.hpp"
#include "duckdb/parser/keyword_helper.hpp"
#include "yyjson.hpp"

#include <cstring>
#include <memory>
#include <mutex>
#include <tuple>
#include <unistd.h>

namespace duckdb {
namespace bt2_daemon {

namespace gb = ::duckdb::miint::gpl_boundary;
namespace yj = duckdb_yyjson;

// =============================================================================
// Config builder + param mapping (used by both align_bowtie2 callers)
// =============================================================================

void ConfigJsonBuilder::append_raw(const std::string &key, const std::string &raw_value) {
	if (!first) {
		os << ",";
	}
	first = false;
	os << "\"" << key << "\":" << raw_value;
}

void ConfigJsonBuilder::append_str(const std::string &key, const std::string &value) {
	// Route every string through gb::JsonEscape so the builder is safe for
	// user-controlled values (rg_id, score_min, ...) — no caller has to
	// remember whether their input is trusted.
	append_raw(key, "\"" + gb::JsonEscape(value) + "\"");
}

void ConfigJsonBuilder::append_int(const std::string &key, int64_t value) {
	append_raw(key, std::to_string(value));
}

void ConfigJsonBuilder::append_bool(const std::string &key, bool value) {
	append_raw(key, value ? "true" : "false");
}

std::string ConfigJsonBuilder::build() const {
	return "{" + os.str() + "}";
}

bool ValueAsBool(const char *caller, const std::string &name, const Value &v) {
	if (v.type().id() != LogicalTypeId::BOOLEAN) {
		throw InvalidInputException("%s: parameter '%s' expects a boolean, got %s", caller, name, v.type().ToString());
	}
	return v.GetValue<bool>();
}

int64_t ValueAsInt(const char *caller, const std::string &name, const Value &v) {
	auto t = v.type().id();
	if (t == LogicalTypeId::INTEGER || t == LogicalTypeId::BIGINT || t == LogicalTypeId::SMALLINT ||
	    t == LogicalTypeId::TINYINT || t == LogicalTypeId::HUGEINT || t == LogicalTypeId::UINTEGER ||
	    t == LogicalTypeId::UBIGINT || t == LogicalTypeId::USMALLINT || t == LogicalTypeId::UTINYINT) {
		return v.GetValue<int64_t>();
	}
	throw InvalidInputException("%s: parameter '%s' expects an integer, got %s", caller, name, v.type().ToString());
}

std::string ValueAsStr(const char *caller, const std::string &name, const Value &v) {
	if (v.type().id() != LogicalTypeId::VARCHAR) {
		throw InvalidInputException("%s: parameter '%s' expects a string, got %s", caller, name, v.type().ToString());
	}
	return v.GetValue<std::string>();
}

const std::unordered_set<std::string> kKnownPresets = {"very-fast", "fast", "sensitive", "very-sensitive"};
const std::unordered_set<std::string> kKnownMateOrientations = {"fr", "rf", "ff"};

// Every bowtie2-align knob the daemon's `--describe` v0.2.0 advertises EXCEPT
// the per-call internals (`index_path`, `nthreads`) and the caller-owned
// miint-side knobs (`threads`, `shard_directory`, `read_to_shard`,
// `max_threads_per_shard`, `include_shard_name`).
//
// `no_unal` is intentionally not exposed: the sharded path hard-codes it to
// true (matches the pre-migration FilterMappedOnly contract) and the
// non-sharded path leaves it at the daemon default (false) so users can see
// unaligned rows when they want them. If someone needs explicit no_unal
// control, expose it then.
const std::unordered_set<std::string> kCommonAlignParams = {
    "preset",
    "local",
    "max_secondary",
    "quiet",
    "seed",
    "trim5",
    "trim3",
    "match_bonus",
    "mismatch_penalty",
    "mismatch_penalty_min",
    "n_penalty",
    "read_gap_open",
    "read_gap_extend",
    "ref_gap_open",
    "ref_gap_extend",
    "score_min",
    "min_insert",
    "max_insert",
    "mate_orientation",
    "no_mixed",
    "no_discordant",
    "dovetail",
    "no_contain",
    "no_overlap",
    "nofw",
    "norc",
    "seed_mismatches",
    "seed_length",
    "max_dp_failures",
    "max_seed_rounds",
    "report_all",
    "xeq",
    "rg_id",
    "ignore_quals",
    "reorder",
    "no_exact_upfront",
    "no_1mm_upfront",
    "deterministic_seeds",
    "lowseeds",
    "memory_mapped",
};

void AppendBowtie2AlignParams(ConfigJsonBuilder &cfg, const named_parameter_map_t &named_params, const char *caller) {
	auto get = [&](const std::string &k) -> const Value * {
		auto it = named_params.find(k);
		return (it == named_params.end() || it->second.IsNull()) ? nullptr : &it->second;
	};

	// 1. Composed preset+local handling — must run before plain preset is
	//    appended because we may need to compose `<preset>-local`.
	bool local = false;
	if (auto *v = get("local")) {
		local = ValueAsBool(caller, "local", *v);
		cfg.append_bool("local_align", local);
	}
	if (auto *v = get("preset")) {
		const std::string p = ValueAsStr(caller, "preset", *v);
		if (kKnownPresets.find(p) == kKnownPresets.end()) {
			throw InvalidInputException("%s: invalid preset '%s' "
			                            "(expected one of very-fast/fast/sensitive/very-sensitive)",
			                            caller, p);
		}
		const std::string preset_value = local ? (p + "-local") : p;
		cfg.append_str("preset", preset_value);
	}

	// 2. miint-renamed knobs — `max_secondary` → daemon `k`.
	if (auto *v = get("max_secondary")) {
		const int64_t n = ValueAsInt(caller, "max_secondary", *v);
		if (n < 0) {
			throw InvalidInputException("%s: max_secondary must be >= 0 (got %lld)", caller, static_cast<long long>(n));
		}
		cfg.append_int("k", n);
	}

	// 3. quiet (miint default true) inverts to daemon verbose=false. Note that
	//    quiet=false has no user-visible effect: the daemon's verbose stderr is
	//    drained and discarded by Session::PumpStderr (never surfaced to SQL).
	//    It only adds overhead — see the `quiet` param docs in table-functions.md.
	bool quiet = true;
	if (auto *v = get("quiet")) {
		quiet = ValueAsBool(caller, "quiet", *v);
	}
	cfg.append_bool("verbose", !quiet);

	// 4. Pass-through integer knobs that are nullable on the daemon side
	//    (omitting them yields the daemon's mode-dependent defaults; only
	//    write when the user supplies a value).
	auto append_nullable_int = [&](const char *param_name, idx_t min_value = 0) {
		if (auto *v = get(param_name)) {
			const int64_t n = ValueAsInt(caller, param_name, *v);
			if (n < static_cast<int64_t>(min_value)) {
				throw InvalidInputException("%s: %s must be >= %lld (got %lld)", caller, param_name,
				                            static_cast<long long>(min_value), static_cast<long long>(n));
			}
			cfg.append_int(param_name, n);
		}
	};
	append_nullable_int("seed");
	append_nullable_int("trim5");
	append_nullable_int("trim3");
	append_nullable_int("match_bonus");
	// match_bonus / mismatch_penalty / n_penalty / gap penalties can in
	// principle be negative on the daemon side (the daemon documents them
	// as nullable integers, no range constraint). We don't reject negatives
	// — bowtie2 itself validates the score function makes sense.
	if (auto *v = get("mismatch_penalty")) {
		cfg.append_int("mismatch_penalty", ValueAsInt(caller, "mismatch_penalty", *v));
	}
	if (auto *v = get("mismatch_penalty_min")) {
		cfg.append_int("mismatch_penalty_min", ValueAsInt(caller, "mismatch_penalty_min", *v));
	}
	if (auto *v = get("n_penalty")) {
		cfg.append_int("n_penalty", ValueAsInt(caller, "n_penalty", *v));
	}
	if (auto *v = get("read_gap_open")) {
		cfg.append_int("read_gap_open", ValueAsInt(caller, "read_gap_open", *v));
	}
	if (auto *v = get("read_gap_extend")) {
		cfg.append_int("read_gap_extend", ValueAsInt(caller, "read_gap_extend", *v));
	}
	if (auto *v = get("ref_gap_open")) {
		cfg.append_int("ref_gap_open", ValueAsInt(caller, "ref_gap_open", *v));
	}
	if (auto *v = get("ref_gap_extend")) {
		cfg.append_int("ref_gap_extend", ValueAsInt(caller, "ref_gap_extend", *v));
	}
	if (auto *v = get("score_min")) {
		// String form like "L,-0.6,-0.6"; daemon validates the format.
		cfg.append_str("score_min", ValueAsStr(caller, "score_min", *v));
	}
	append_nullable_int("min_insert");
	append_nullable_int("max_insert");
	if (auto *v = get("mate_orientation")) {
		const std::string mo = ValueAsStr(caller, "mate_orientation", *v);
		if (kKnownMateOrientations.find(mo) == kKnownMateOrientations.end()) {
			throw InvalidInputException("%s: invalid mate_orientation '%s' (expected one of fr/rf/ff)", caller, mo);
		}
		cfg.append_str("mate_orientation", mo);
	}

	// 5. Pass-through bools, daemon defaults all false.
	auto pass_bool = [&](const char *param_name) {
		if (auto *v = get(param_name)) {
			cfg.append_bool(param_name, ValueAsBool(caller, param_name, *v));
		}
	};
	pass_bool("no_mixed");
	pass_bool("no_discordant");
	pass_bool("dovetail");
	pass_bool("no_contain");
	pass_bool("no_overlap");
	pass_bool("nofw");
	pass_bool("norc");

	// 6. Seed/DP tuning — nullable integers per the daemon schema.
	if (auto *v = get("seed_mismatches")) {
		const int64_t n = ValueAsInt(caller, "seed_mismatches", *v);
		if (n < 0 || n > 1) {
			throw InvalidInputException("%s: seed_mismatches must be 0 or 1 (got %lld)", caller,
			                            static_cast<long long>(n));
		}
		cfg.append_int("seed_mismatches", n);
	}
	if (auto *v = get("seed_length")) {
		const int64_t n = ValueAsInt(caller, "seed_length", *v);
		if (n < 1 || n > 32) {
			throw InvalidInputException("%s: seed_length must be 1-32 (got %lld)", caller, static_cast<long long>(n));
		}
		cfg.append_int("seed_length", n);
	}
	append_nullable_int("max_dp_failures", 1);
	append_nullable_int("max_seed_rounds", 1);

	// 7. Output formatting / misc.
	pass_bool("report_all");
	pass_bool("xeq");
	if (auto *v = get("rg_id")) {
		cfg.append_str("rg_id", ValueAsStr(caller, "rg_id", *v));
	}
	pass_bool("ignore_quals");
	pass_bool("reorder");

	// 8. v0.4 seed-behavior knobs (bowtie2 -d / --no-exact-upfront /
	//    --no-1mm-upfront / -l). `deterministic_seeds` gives reproducible seed
	//    selection but bowtie2 couples it: it requires report_all AND both
	//    upfront-disable flags set, and is incompatible with -k. We forward the
	//    values verbatim and let the daemon's config parser enforce the coupling
	//    (it returns a clear "deterministic_seeds requires ..." error) rather
	//    than duplicate that interdependency here.
	pass_bool("no_exact_upfront");
	pass_bool("no_1mm_upfront");
	pass_bool("deterministic_seeds");
	if (auto *v = get("lowseeds")) {
		cfg.append_str("lowseeds", ValueAsStr(caller, "lowseeds", *v));
	}

	// 9. Index loading (gpl-boundary v0.4.2+). `--mm` is the daemon's default, so
	//    UNLIKE the section-5 bools we only emit a value when the user sets one —
	//    omitting it preserves the daemon's mmap-on behavior. A diagnostic/perf
	//    lever: `memory_mapped:=false` makes the worker fread the whole index
	//    sequentially (vs lazy random page-faults), which can win on a cold
	//    network FS at the cost of anonymous RSS. Output-invariant. Silently
	//    ignored by pre-0.4.2 daemons (which always mmap).
	pass_bool("memory_mapped");
}

void RegisterBowtie2AlignNamedParameterTypes(TableFunction &tf) {
	tf.named_parameters["preset"] = LogicalType::VARCHAR;
	tf.named_parameters["local"] = LogicalType::BOOLEAN;
	tf.named_parameters["max_secondary"] = LogicalType::INTEGER;
	tf.named_parameters["quiet"] = LogicalType::BOOLEAN;

	tf.named_parameters["seed"] = LogicalType::INTEGER;
	tf.named_parameters["trim5"] = LogicalType::INTEGER;
	tf.named_parameters["trim3"] = LogicalType::INTEGER;
	tf.named_parameters["match_bonus"] = LogicalType::INTEGER;
	tf.named_parameters["mismatch_penalty"] = LogicalType::INTEGER;
	tf.named_parameters["mismatch_penalty_min"] = LogicalType::INTEGER;
	tf.named_parameters["n_penalty"] = LogicalType::INTEGER;
	tf.named_parameters["read_gap_open"] = LogicalType::INTEGER;
	tf.named_parameters["read_gap_extend"] = LogicalType::INTEGER;
	tf.named_parameters["ref_gap_open"] = LogicalType::INTEGER;
	tf.named_parameters["ref_gap_extend"] = LogicalType::INTEGER;
	tf.named_parameters["score_min"] = LogicalType::VARCHAR;
	tf.named_parameters["min_insert"] = LogicalType::INTEGER;
	tf.named_parameters["max_insert"] = LogicalType::INTEGER;
	tf.named_parameters["mate_orientation"] = LogicalType::VARCHAR;
	tf.named_parameters["no_mixed"] = LogicalType::BOOLEAN;
	tf.named_parameters["no_discordant"] = LogicalType::BOOLEAN;
	tf.named_parameters["dovetail"] = LogicalType::BOOLEAN;
	tf.named_parameters["no_contain"] = LogicalType::BOOLEAN;
	tf.named_parameters["no_overlap"] = LogicalType::BOOLEAN;
	tf.named_parameters["nofw"] = LogicalType::BOOLEAN;
	tf.named_parameters["norc"] = LogicalType::BOOLEAN;
	tf.named_parameters["seed_mismatches"] = LogicalType::INTEGER;
	tf.named_parameters["seed_length"] = LogicalType::INTEGER;
	tf.named_parameters["max_dp_failures"] = LogicalType::INTEGER;
	tf.named_parameters["max_seed_rounds"] = LogicalType::INTEGER;
	tf.named_parameters["report_all"] = LogicalType::BOOLEAN;
	tf.named_parameters["xeq"] = LogicalType::BOOLEAN;
	tf.named_parameters["rg_id"] = LogicalType::VARCHAR;
	tf.named_parameters["ignore_quals"] = LogicalType::BOOLEAN;
	tf.named_parameters["reorder"] = LogicalType::BOOLEAN;
	tf.named_parameters["no_exact_upfront"] = LogicalType::BOOLEAN;
	tf.named_parameters["no_1mm_upfront"] = LogicalType::BOOLEAN;
	tf.named_parameters["deterministic_seeds"] = LogicalType::BOOLEAN;
	tf.named_parameters["lowseeds"] = LogicalType::VARCHAR;
	tf.named_parameters["memory_mapped"] = LogicalType::BOOLEAN;
}

const char *const kOutputColumnNames[kNumOutputColumns] = {
    "read_id",        "flags",         "reference",       "position", "stop_position", "mapq",   "cigar",
    "mate_reference", "mate_position", "template_length", "tag_as",   "tag_xs",        "tag_ys", "tag_xn",
    "tag_xm",         "tag_xo",        "tag_xg",          "tag_nm",   "tag_yt",        "tag_md", "tag_sa"};

// Arrow C Data Interface format codes:
//   "u" = Utf8, "S" = uint16, "C" = uint8, "l" = int64, "i" = int32.
// EmitChunkRows reads fixed-width columns as the listed type — a daemon-side
// type drift without a schema_version bump would corrupt output if we trusted
// the name alone.
const char *const kOutputColumnFormats[kNumOutputColumns] = {
    "u", // read_id          — Utf8
    "S", // flags            — uint16
    "u", // reference        — Utf8
    "l", // position         — int64
    "l", // stop_position    — int64
    "C", // mapq             — uint8
    "u", // cigar            — Utf8
    "u", // mate_reference   — Utf8
    "l", // mate_position    — int64
    "l", // template_length  — int64
    "i", // tag_as           — int32 (widened to BIGINT in EmitChunkRows)
    "i", // tag_xs           — int32 (widened)
    "i", // tag_ys           — int32 (widened)
    "i", // tag_xn           — int32 (widened)
    "i", // tag_xm           — int32 (widened)
    "i", // tag_xo           — int32 (widened)
    "i", // tag_xg           — int32 (widened)
    "i", // tag_nm           — int32 (widened)
    "u", // tag_yt           — Utf8
    "u", // tag_md           — Utf8
    "u", // tag_sa           — Utf8
};

void PopulateOutputSchema(std::vector<std::string> &names, std::vector<LogicalType> &types,
                          const LogicalType &query_id_type, const LogicalType &subject_id_type) {
	// Delegate to the canonical helpers in align_common.hpp so the bowtie2
	// daemon path doesn't carry a second copy of the alignment schema. The
	// names match `kOutputColumnNames` (which stays as the wire-validation
	// truth source for ValidateOutputSchema); GetAlignmentOutputNames returns
	// the same 21 names in the same order.
	names = GetAlignmentOutputNames();
	types = GetAlignmentOutputTypes(query_id_type, subject_id_type);
}

void DecodeListQualToPhred33(const Value &v, const char *col_name, const std::string &query_table, std::string &out,
                             std::vector<uint8_t> &scratch) {
	const auto &children = ListValue::GetChildren(v);
	scratch.clear();
	scratch.reserve(children.size());
	for (auto &c : children) {
		if (c.IsNull()) {
			throw InvalidInputException("align_bowtie2: NULL element inside '%s' array in query table '%s' (per-base "
			                            "quality must be fully populated)",
			                            col_name, query_table);
		}
		scratch.push_back(c.GetValue<uint8_t>());
	}
	EncodePhred33Quality(scratch.data(), scratch.size(), 33, out);
}

void ValidateOutputSchema(const ArrowSchema &schema) {
	// Generic "bowtie2-align" rather than "align_bowtie2:" because both
	// align_bowtie2 and align_bowtie2_sharded reach this code; surfacing the
	// daemon tool name keeps the diagnostic accurate regardless of caller.
	if (schema.n_children != kNumOutputColumns) {
		throw IOException("bowtie2-align: daemon returned unexpected schema (%lld columns, expected %d)",
		                  static_cast<long long>(schema.n_children), kNumOutputColumns);
	}
	for (int c = 0; c < kNumOutputColumns; ++c) {
		const auto *child = schema.children[c];
		if (!child || !child->name || std::strcmp(child->name, kOutputColumnNames[c]) != 0) {
			throw IOException("bowtie2-align: schema drift at column %d — expected name '%s', got '%s'", c,
			                  kOutputColumnNames[c], (child && child->name) ? child->name : "(null)");
		}
		if (!child->format || std::strcmp(child->format, kOutputColumnFormats[c]) != 0) {
			// Type drift is worse than name drift: name-only checks miss a
			// silent int64→int32 swap that `read_fixed<int64_t>` would
			// reinterpret as 4 bytes of garbage from the next column.
			throw IOException("bowtie2-align: schema drift at column '%s' (idx %d) — expected Arrow format '%s', "
			                  "got '%s'",
			                  kOutputColumnNames[c], c, kOutputColumnFormats[c],
			                  child->format ? child->format : "(null)");
		}
	}
}

// =============================================================================
// BuildQueryIpc — nanoarrow encoder for {read_id, sequence1, sequence2?,
// qual1?, qual2?}. Matches gpl-boundary's bowtie2-align input schema.
// =============================================================================

std::vector<uint8_t> BuildQueryIpc(const QueryBatch &qb, const QueryArrowSchema &schema_flags) {
	const int n_cols = schema_flags.num_columns();
	ArrowSchema schema {};
	auto rc = ArrowSchemaInitFromType(&schema, NANOARROW_TYPE_STRUCT);
	if (rc != NANOARROW_OK) {
		throw InternalException("align_bowtie2: ArrowSchemaInit failed");
	}
	rc = ArrowSchemaAllocateChildren(&schema, n_cols);
	if (rc != NANOARROW_OK) {
		schema.release(&schema);
		throw InternalException("align_bowtie2: ArrowSchemaAllocateChildren failed");
	}
	int idx = 0;
	ArrowSchemaInitFromType(schema.children[idx], NANOARROW_TYPE_STRING);
	ArrowSchemaSetName(schema.children[idx], "read_id");
	++idx;
	ArrowSchemaInitFromType(schema.children[idx], NANOARROW_TYPE_STRING);
	ArrowSchemaSetName(schema.children[idx], "sequence1");
	++idx;
	if (schema_flags.has_sequence2) {
		ArrowSchemaInitFromType(schema.children[idx], NANOARROW_TYPE_STRING);
		ArrowSchemaSetName(schema.children[idx], "sequence2");
		++idx;
	}
	if (schema_flags.has_qual1) {
		ArrowSchemaInitFromType(schema.children[idx], NANOARROW_TYPE_STRING);
		ArrowSchemaSetName(schema.children[idx], "qual1");
		++idx;
	}
	if (schema_flags.has_qual2) {
		ArrowSchemaInitFromType(schema.children[idx], NANOARROW_TYPE_STRING);
		ArrowSchemaSetName(schema.children[idx], "qual2");
		++idx;
	}

	ArrowArray array {};
	ArrowError err {};
	if (ArrowArrayInitFromSchema(&array, &schema, &err) != NANOARROW_OK) {
		schema.release(&schema);
		throw InternalException("align_bowtie2: ArrowArrayInit failed: %s", err.message);
	}
	if (ArrowArrayStartAppending(&array) != NANOARROW_OK) {
		array.release(&array);
		schema.release(&schema);
		throw InternalException("align_bowtie2: ArrowArrayStartAppending failed");
	}

	const idx_t n_rows = qb.read_ids.size();
	for (idx_t row = 0; row < n_rows; ++row) {
		int c = 0;
		ArrowStringView rv {qb.read_ids[row].data(), static_cast<int64_t>(qb.read_ids[row].size())};
		if (ArrowArrayAppendString(array.children[c++], rv) != NANOARROW_OK) {
			array.release(&array);
			schema.release(&schema);
			throw InternalException("align_bowtie2: read_id append failed at row %lld", static_cast<long long>(row));
		}
		ArrowStringView s1v {qb.sequence1[row].data(), static_cast<int64_t>(qb.sequence1[row].size())};
		if (ArrowArrayAppendString(array.children[c++], s1v) != NANOARROW_OK) {
			array.release(&array);
			schema.release(&schema);
			throw InternalException("align_bowtie2: sequence1 append failed at row %lld", static_cast<long long>(row));
		}
		if (schema_flags.has_sequence2) {
			ArrowArray *child = array.children[c++];
			if (row < qb.sequence2_valid.size() && qb.sequence2_valid[row]) {
				ArrowStringView v {qb.sequence2[row].data(), static_cast<int64_t>(qb.sequence2[row].size())};
				if (ArrowArrayAppendString(child, v) != NANOARROW_OK) {
					array.release(&array);
					schema.release(&schema);
					throw InternalException("align_bowtie2: sequence2 append failed");
				}
			} else {
				if (ArrowArrayAppendNull(child, 1) != NANOARROW_OK) {
					array.release(&array);
					schema.release(&schema);
					throw InternalException("align_bowtie2: sequence2 null append failed");
				}
			}
		}
		if (schema_flags.has_qual1) {
			ArrowArray *child = array.children[c++];
			if (row < qb.qual1_valid.size() && qb.qual1_valid[row]) {
				ArrowStringView v {qb.qual1[row].data(), static_cast<int64_t>(qb.qual1[row].size())};
				if (ArrowArrayAppendString(child, v) != NANOARROW_OK) {
					array.release(&array);
					schema.release(&schema);
					throw InternalException("align_bowtie2: qual1 append failed");
				}
			} else {
				if (ArrowArrayAppendNull(child, 1) != NANOARROW_OK) {
					array.release(&array);
					schema.release(&schema);
					throw InternalException("align_bowtie2: qual1 null append failed");
				}
			}
		}
		if (schema_flags.has_qual2) {
			ArrowArray *child = array.children[c++];
			if (row < qb.qual2_valid.size() && qb.qual2_valid[row]) {
				ArrowStringView v {qb.qual2[row].data(), static_cast<int64_t>(qb.qual2[row].size())};
				if (ArrowArrayAppendString(child, v) != NANOARROW_OK) {
					array.release(&array);
					schema.release(&schema);
					throw InternalException("align_bowtie2: qual2 append failed");
				}
			} else {
				if (ArrowArrayAppendNull(child, 1) != NANOARROW_OK) {
					array.release(&array);
					schema.release(&schema);
					throw InternalException("align_bowtie2: qual2 null append failed");
				}
			}
		}
		if (ArrowArrayFinishElement(&array) != NANOARROW_OK) {
			array.release(&array);
			schema.release(&schema);
			throw InternalException("align_bowtie2: FinishElement failed at row %lld", static_cast<long long>(row));
		}
	}
	if (ArrowArrayFinishBuildingDefault(&array, &err) != NANOARROW_OK) {
		array.release(&array);
		schema.release(&schema);
		throw InternalException("align_bowtie2: ArrowArrayFinishBuilding failed: %s", err.message);
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

// =============================================================================
// EmitChunkRows — extract one DuckDB chunk's worth of rows from a decoded
// Arrow batch. Matches phylogeny_fasttree.cpp:711-737 decoder pattern.
// =============================================================================

namespace {

template <typename T>
T read_fixed(const ArrowArray &col, idx_t logical_index) {
	const auto *buf = reinterpret_cast<const T *>(col.buffers[1]);
	return buf[static_cast<idx_t>(col.offset) + logical_index];
}

bool is_null_at(const ArrowArray &arr, idx_t logical_index) {
	if (!arr.buffers[0]) {
		return false;
	}
	if (arr.null_count == 0) {
		return false;
	}
	const auto *bitmap = static_cast<const uint8_t *>(arr.buffers[0]);
	const idx_t abs = static_cast<idx_t>(arr.offset) + logical_index;
	return (bitmap[abs / 8] & (1u << (abs % 8))) == 0;
}

void emit_string(Vector &out, idx_t out_row, const ArrowArray &col, idx_t logical_index) {
	const auto *offsets = static_cast<const int32_t *>(col.buffers[1]);
	const auto *data = static_cast<const char *>(col.buffers[2]);
	const idx_t a = static_cast<idx_t>(col.offset) + logical_index;
	const int32_t start = offsets[a];
	const int32_t end = offsets[a + 1];
	const int32_t len = end - start;
	if (len < 0) {
		throw IOException("align_bowtie2: corrupt utf8 offsets at row %lld (start=%d end=%d)",
		                  static_cast<long long>(a), start, end);
	}
	FlatVector::GetData<string_t>(out)[out_row] = StringVector::AddString(out, data + start, static_cast<idx_t>(len));
}

// Read an Arrow utf8 cell into an owned std::string. Used by the BIGINT id-
// column path where the value also needs to be passed to the codec parser
// (which takes `const std::string&`) and, for mate_reference, possibly
// substituted in place during "=" resolution.
std::string read_arrow_string_owned(const ArrowArray &col, idx_t logical_index) {
	const auto *offsets = static_cast<const int32_t *>(col.buffers[1]);
	const auto *data = static_cast<const char *>(col.buffers[2]);
	const idx_t a = static_cast<idx_t>(col.offset) + logical_index;
	const int32_t start = offsets[a];
	const int32_t end = offsets[a + 1];
	const int32_t len = end - start;
	if (len < 0) {
		throw IOException("align_bowtie2: corrupt utf8 offsets at row %lld (start=%d end=%d)",
		                  static_cast<long long>(a), start, end);
	}
	return std::string(data + start, static_cast<size_t>(len));
}

} // namespace

void EmitChunkRows(DataChunk &output, idx_t to_emit, idx_t row_start, const ArrowArray &batch,
                   const LogicalType &query_id_type, const LogicalType &subject_id_type) {
	if (!IsAllowedIdType(query_id_type) || !IsAllowedIdType(subject_id_type)) {
		throw InternalException("align_bowtie2 EmitChunkRows: id types must be %s (got query=%s, subject=%s)",
		                        AllowedIdTypeList(), query_id_type.ToString(), subject_id_type.ToString());
	}
	// Non-VARCHAR ids (BIGINT, UUID) round-trip through the codec: read the
	// Arrow string and let EmitIdCell parse it back to the storage type. VARCHAR
	// ids copy straight from Arrow via emit_string.
	const bool query_needs_codec = query_id_type.id() != LogicalTypeId::VARCHAR;
	const bool subject_needs_codec = subject_id_type.id() != LogicalTypeId::VARCHAR;
	auto &v_read_id = output.data[0];
	auto *out_flags = FlatVector::GetData<uint16_t>(output.data[1]);
	auto &v_reference = output.data[2];
	auto *out_position = FlatVector::GetData<int64_t>(output.data[3]);
	auto *out_stop = FlatVector::GetData<int64_t>(output.data[4]);
	auto *out_mapq = FlatVector::GetData<uint8_t>(output.data[5]);
	auto &v_cigar = output.data[6];
	auto &v_mate_ref = output.data[7];
	auto *out_mate_pos = FlatVector::GetData<int64_t>(output.data[8]);
	auto *out_tlen = FlatVector::GetData<int64_t>(output.data[9]);
	auto *out_tag_as = FlatVector::GetData<int64_t>(output.data[10]);
	auto *out_tag_xs = FlatVector::GetData<int64_t>(output.data[11]);
	auto *out_tag_ys = FlatVector::GetData<int64_t>(output.data[12]);
	auto *out_tag_xn = FlatVector::GetData<int64_t>(output.data[13]);
	auto *out_tag_xm = FlatVector::GetData<int64_t>(output.data[14]);
	auto *out_tag_xo = FlatVector::GetData<int64_t>(output.data[15]);
	auto *out_tag_xg = FlatVector::GetData<int64_t>(output.data[16]);
	auto *out_tag_nm = FlatVector::GetData<int64_t>(output.data[17]);
	auto &v_tag_yt = output.data[18];
	auto &v_tag_md = output.data[19];
	auto &v_tag_sa = output.data[20];

	auto &mask_tag_as = FlatVector::Validity(output.data[10]);
	auto &mask_tag_xs = FlatVector::Validity(output.data[11]);
	auto &mask_tag_ys = FlatVector::Validity(output.data[12]);
	auto &mask_tag_xn = FlatVector::Validity(output.data[13]);
	auto &mask_tag_xm = FlatVector::Validity(output.data[14]);
	auto &mask_tag_xo = FlatVector::Validity(output.data[15]);
	auto &mask_tag_xg = FlatVector::Validity(output.data[16]);
	auto &mask_tag_nm = FlatVector::Validity(output.data[17]);
	auto &mask_tag_yt = FlatVector::Validity(v_tag_yt);
	auto &mask_tag_md = FlatVector::Validity(v_tag_md);
	auto &mask_tag_sa = FlatVector::Validity(v_tag_sa);

	const auto &col_read_id = *batch.children[0];
	const auto &col_flags = *batch.children[1];
	const auto &col_reference = *batch.children[2];
	const auto &col_position = *batch.children[3];
	const auto &col_stop = *batch.children[4];
	const auto &col_mapq = *batch.children[5];
	const auto &col_cigar = *batch.children[6];
	const auto &col_mate_ref = *batch.children[7];
	const auto &col_mate_pos = *batch.children[8];
	const auto &col_tlen = *batch.children[9];
	const auto &col_tag_as = *batch.children[10];
	const auto &col_tag_xs = *batch.children[11];
	const auto &col_tag_ys = *batch.children[12];
	const auto &col_tag_xn = *batch.children[13];
	const auto &col_tag_xm = *batch.children[14];
	const auto &col_tag_xo = *batch.children[15];
	const auto &col_tag_xg = *batch.children[16];
	const auto &col_tag_nm = *batch.children[17];
	const auto &col_tag_yt = *batch.children[18];
	const auto &col_tag_md = *batch.children[19];
	const auto &col_tag_sa = *batch.children[20];

	for (idx_t i = 0; i < to_emit; ++i) {
		const idx_t li = row_start + i;

		// read_id: dispatch on query id type.
		if (query_needs_codec) {
			EmitIdCell(v_read_id, i, read_arrow_string_owned(col_read_id, li), query_id_type);
		} else {
			emit_string(v_read_id, i, col_read_id, li);
		}
		out_flags[i] = read_fixed<uint16_t>(col_flags, li);

		// reference + mate_reference: dispatch on subject id type. For non-VARCHAR
		// subjects, cache the reference string so mate_reference "=" can resolve
		// into it without re-reading the Arrow cell.
		if (subject_needs_codec) {
			const std::string ref_str = read_arrow_string_owned(col_reference, li);
			EmitIdCell(v_reference, i, ref_str, subject_id_type);
			std::string mr_str = read_arrow_string_owned(col_mate_ref, li);
			if (mr_str == "=") {
				// SAM "=" sentinel: mate maps to the same reference as the
				// primary alignment. The literal "=" has no BIGINT/UUID encoding,
				// so substitute the row's reference value before parsing.
				mr_str = ref_str;
			}
			EmitIdCell(v_mate_ref, i, mr_str, subject_id_type);
		} else {
			emit_string(v_reference, i, col_reference, li);
			emit_string(v_mate_ref, i, col_mate_ref, li);
		}
		out_position[i] = read_fixed<int64_t>(col_position, li);
		out_stop[i] = read_fixed<int64_t>(col_stop, li);
		out_mapq[i] = read_fixed<uint8_t>(col_mapq, li);
		emit_string(v_cigar, i, col_cigar, li);
		out_mate_pos[i] = read_fixed<int64_t>(col_mate_pos, li);
		out_tlen[i] = read_fixed<int64_t>(col_tlen, li);

		auto widen = [&](const ArrowArray &col, int64_t *out, ValidityMask &mask) {
			if (is_null_at(col, li)) {
				mask.SetInvalid(i);
			} else {
				out[i] = static_cast<int64_t>(read_fixed<int32_t>(col, li));
			}
		};
		widen(col_tag_as, out_tag_as, mask_tag_as);
		widen(col_tag_xs, out_tag_xs, mask_tag_xs);
		widen(col_tag_ys, out_tag_ys, mask_tag_ys);
		widen(col_tag_xn, out_tag_xn, mask_tag_xn);
		widen(col_tag_xm, out_tag_xm, mask_tag_xm);
		widen(col_tag_xo, out_tag_xo, mask_tag_xo);
		widen(col_tag_xg, out_tag_xg, mask_tag_xg);
		widen(col_tag_nm, out_tag_nm, mask_tag_nm);

		auto emit_nullable_str = [&](Vector &v, ValidityMask &mask, const ArrowArray &col) {
			if (is_null_at(col, li)) {
				mask.SetInvalid(i);
			} else {
				emit_string(v, i, col, li);
			}
		};
		emit_nullable_str(v_tag_yt, mask_tag_yt, col_tag_yt);
		emit_nullable_str(v_tag_md, mask_tag_md, col_tag_md);
		emit_nullable_str(v_tag_sa, mask_tag_sa, col_tag_sa);
	}
}

void RequireGplBoundaryVersion(const std::string &binary_path, const char *caller) {
	// The daemon's release version is a process-level constant, but the sharded
	// SpawnAndCheckSession runs once PER worker thread (InitLocal) — so validate
	// once and cache the success. A failed check throws and leaves the once_flag
	// unset, so a later query (e.g. after the user upgrades the binary) re-checks;
	// a too-old daemon therefore always fails loud, never silently proceeds.
	static std::once_flag flag;
	std::call_once(flag, [&]() {
		// Minimum gpl-boundary release that supports the bowtie2 `memory_mapped`
		// (--mm toggle) config field. Below this, the field is silently ignored.
		constexpr int kReqMajor = 0, kReqMinor = 4, kReqPatch = 2;

		int major = 0, minor = 0, patch = 0;
		bool parsed = false;
		try {
			// Run `<binary> --version` and drain its small JSON stdout
			// ({"gpl_boundary":"X.Y.Z", ...}). Mirrors the drain pattern in
			// align_bowtie2.cpp's bowtie2_available() probe.
			std::vector<std::string> argv = {binary_path, "--version"};
			gb::ChildProcess child(argv);
			auto drain = [](int fd, std::string &dst) {
				if (fd < 0) {
					return;
				}
				char buf[256];
				ssize_t n;
				while ((n = ::read(fd, buf, sizeof(buf))) > 0) {
					dst.append(buf, static_cast<size_t>(n));
				}
			};
			std::string out;
			std::string err;
			drain(child.stdout_fd(), out);
			drain(child.stderr_fd(), err);
			child.Wait();
			parsed = gb::ParseGplBoundaryVersion(out, major, minor, patch);
		} catch (...) {
			// Any spawn/read failure → treat as "version undetermined" and fail
			// loud below, rather than proceeding against an unknown daemon.
			parsed = false;
		}

		const bool at_least =
		    parsed && std::tie(major, minor, patch) >= std::make_tuple(kReqMajor, kReqMinor, kReqPatch);
		if (!at_least) {
			const std::string found =
			    parsed ? (std::to_string(major) + "." + std::to_string(minor) + "." + std::to_string(patch))
			           : std::string("unknown (could not parse `gpl-boundary --version`)");
			throw IOException(
			    "%s requires gpl-boundary >= 0.4.2 (for the bowtie2 `memory_mapped` index-load control), but the "
			    "binary at %s reports version %s. Upgrade with `SELECT install_gpl_boundary();`, or point "
			    "MIINT_GPL_BOUNDARY_PATH at a >= 0.4.2 binary.",
			    caller, binary_path, found);
		}
	});
}

// =============================================================================
// bowtie2-build glue (shared by align_bowtie2 + save_bowtie2_index)
// =============================================================================

LoadedSubjects LoadSingleEndSubjects(ClientContext &context, const std::string &table_name, const char *caller) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	// Probe for an optional sequence2 column first (paired subjects are
	// rejected), then issue the actual SELECT. read_id may be VARCHAR or BIGINT;
	// the implicit Value::GetValue<std::string>() cast below handles either.
	const std::string columns_sql = "SELECT column_name FROM (DESCRIBE " +
	                                KeywordHelper::WriteOptionallyQuoted(table_name) +
	                                ") WHERE column_name IN ('sequence2')";
	auto columns_res = conn.Query(columns_sql);
	if (columns_res->HasError()) {
		throw InvalidInputException("%s: failed to introspect subject table '%s': %s", caller, table_name,
		                            columns_res->GetError());
	}
	const bool has_sequence2 = columns_res->RowCount() > 0;

	std::string select_sql = "SELECT read_id, sequence1";
	if (has_sequence2) {
		select_sql += ", sequence2";
	}
	select_sql += " FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);

	auto result = conn.Query(select_sql);
	if (result->HasError()) {
		throw InvalidInputException("%s: failed to read subject table '%s' "
		                            "(must have columns read_id VARCHAR or BIGINT, sequence1 VARCHAR): %s",
		                            caller, table_name, result->GetError());
	}

	LoadedSubjects out;
	auto &materialized = result->Cast<MaterializedQueryResult>();
	out.names.reserve(materialized.RowCount());
	out.sequences.reserve(materialized.RowCount());
	while (auto chunk = materialized.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); ++i) {
			auto name_val = chunk->GetValue(0, i);
			auto seq_val = chunk->GetValue(1, i);
			if (name_val.IsNull() || seq_val.IsNull()) {
				throw InvalidInputException("%s: NULL read_id or sequence1 in subject table '%s'", caller, table_name);
			}
			if (has_sequence2) {
				auto s2_val = chunk->GetValue(2, i);
				if (!s2_val.IsNull()) {
					throw InvalidInputException("%s: subject table '%s' has non-NULL sequence2 — "
					                            "subjects must be single-end (sequence2 must be NULL)",
					                            caller, table_name);
				}
			}
			out.names.push_back(name_val.GetValue<std::string>());
			out.sequences.push_back(seq_val.GetValue<std::string>());
		}
	}
	if (out.names.empty()) {
		throw InvalidInputException("%s: subject table '%s' is empty", caller, table_name);
	}
	return out;
}

std::vector<uint8_t> BuildSubjectsIpc(const LoadedSubjects &subjects, const char *caller) {
	ArrowSchema schema {};
	auto rc = ArrowSchemaInitFromType(&schema, NANOARROW_TYPE_STRUCT);
	if (rc != NANOARROW_OK) {
		throw InternalException("%s: ArrowSchemaInit failed", caller);
	}
	rc = ArrowSchemaAllocateChildren(&schema, 2);
	if (rc != NANOARROW_OK) {
		schema.release(&schema);
		throw InternalException("%s: ArrowSchemaAllocateChildren failed", caller);
	}
	ArrowSchemaInitFromType(schema.children[0], NANOARROW_TYPE_STRING);
	ArrowSchemaSetName(schema.children[0], "name");
	ArrowSchemaInitFromType(schema.children[1], NANOARROW_TYPE_STRING);
	ArrowSchemaSetName(schema.children[1], "sequence");

	ArrowArray array {};
	ArrowError err {};
	if (ArrowArrayInitFromSchema(&array, &schema, &err) != NANOARROW_OK) {
		schema.release(&schema);
		throw InternalException("%s: ArrowArrayInit failed: %s", caller, err.message);
	}
	if (ArrowArrayStartAppending(&array) != NANOARROW_OK) {
		array.release(&array);
		schema.release(&schema);
		throw InternalException("%s: ArrowArrayStartAppending failed", caller);
	}
	for (size_t i = 0; i < subjects.names.size(); ++i) {
		ArrowStringView nv {subjects.names[i].data(), static_cast<int64_t>(subjects.names[i].size())};
		ArrowStringView sv {subjects.sequences[i].data(), static_cast<int64_t>(subjects.sequences[i].size())};
		if (ArrowArrayAppendString(array.children[0], nv) != NANOARROW_OK ||
		    ArrowArrayAppendString(array.children[1], sv) != NANOARROW_OK ||
		    ArrowArrayFinishElement(&array) != NANOARROW_OK) {
			array.release(&array);
			schema.release(&schema);
			throw InternalException("%s: subjects append failed at row %d", caller, static_cast<int>(i));
		}
	}
	if (ArrowArrayFinishBuildingDefault(&array, &err) != NANOARROW_OK) {
		array.release(&array);
		schema.release(&schema);
		throw InternalException("%s: ArrowArrayFinishBuilding failed: %s", caller, err.message);
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

std::string BuildBowtie2BuildConfigJson(const std::string &index_basename, int64_t nthreads) {
	ConfigJsonBuilder cfg;
	cfg.append_str("index_path", index_basename);
	if (nthreads > 1) {
		cfg.append_int("nthreads", nthreads);
	}
	return cfg.build();
}

std::vector<std::string> ParseBowtie2BuildIndexFiles(const std::string &result_json, const char *caller) {
	using YyjsonDocPtr = std::unique_ptr<yj::yyjson_doc, decltype(&yj::yyjson_doc_free)>;
	YyjsonDocPtr doc(yj::yyjson_read(result_json.data(), result_json.size(), 0), &yj::yyjson_doc_free);
	if (!doc) {
		throw IOException("%s: bowtie2-build result_json was not valid JSON: %s", caller, result_json);
	}
	yj::yyjson_val *root = yj::yyjson_doc_get_root(doc.get());
	if (!yj::yyjson_is_obj(root)) {
		throw IOException("%s: bowtie2-build result was not a JSON object: %s", caller, result_json);
	}
	yj::yyjson_val *arr = yj::yyjson_obj_get(root, "index_files");
	if (!arr || !yj::yyjson_is_arr(arr)) {
		throw IOException("%s: bowtie2-build result missing 'index_files' array: %s", caller, result_json);
	}
	std::vector<std::string> out;
	const size_t n = yj::yyjson_arr_size(arr);
	out.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		yj::yyjson_val *item = yj::yyjson_arr_get(arr, i);
		if (!yj::yyjson_is_str(item)) {
			throw IOException("%s: bowtie2-build index_files contained non-string entry", caller);
		}
		out.emplace_back(yj::yyjson_get_str(item));
	}
	return out;
}

} // namespace bt2_daemon
} // namespace duckdb
