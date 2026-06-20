#pragma once

// Shared helpers between align_bowtie2 and align_bowtie2_sharded once both
// route through the gpl-boundary daemon's bowtie2-align tool. Lives in its
// own translation unit so neither caller has to drag in the other's
// table-function plumbing.

#include "duckdb/common/arrow/arrow.hpp"
#include "duckdb/common/arrow/arrow_wrapper.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"
#include "duckdb/function/table_function.hpp"

#include <cstdint>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

// nanoarrow C ABI is opaque here; the source file pulls in the full header.
struct ArrowSchema;
struct ArrowArray;

namespace duckdb {
namespace bt2_daemon {

// gpl-boundary bowtie2-align output schema_version we wire against
// (v0.2.0 commit a28cf56). Strict equality check at session bind.
constexpr uint32_t kAlignSchemaVersion = 2;

constexpr int kNumOutputColumns = 21;

// User-facing column order; MUST match the daemon's bowtie2-align v2 wire
// schema so wire column i → output column i.
extern const char *const kOutputColumnNames[kNumOutputColumns];

// Arrow C Data Interface format string per column; matches the daemon's
// bowtie2-align v2 wire schema. Validated structurally so a type drift on the
// daemon side fails loud rather than reinterpreting bytes through
// `read_fixed<T>` with the wrong width.
extern const char *const kOutputColumnFormats[kNumOutputColumns];

// Populate the 21-column output schema. `query_id_type` drives `read_id`;
// `subject_id_type` drives `reference` and `mate_reference`. Both must be
// VARCHAR or BIGINT — no defaults so every caller (align_bowtie2 +
// align_bowtie2_sharded) explicitly commits to a type, matching the
// no-default convention from PR 1 of the id-codec work.
void PopulateOutputSchema(std::vector<std::string> &names, std::vector<LogicalType> &types,
                          const LogicalType &query_id_type, const LogicalType &subject_id_type);

// Validate that an Arrow IPC schema returned by the daemon matches the
// expected 21 named columns. Throws IOException on drift.
void ValidateOutputSchema(const ArrowSchema &schema);

// Describe the column layout for a query Arrow batch. The query Arrow
// batch always has read_id+sequence1; sequence2/qual1/qual2 are optional
// based on the user's input table.
struct QueryArrowSchema {
	bool has_sequence2 = false;
	bool has_qual1 = false;
	bool has_qual2 = false;
	int num_columns() const {
		return 2 + (has_sequence2 ? 1 : 0) + (has_qual1 ? 1 : 0) + (has_qual2 ? 1 : 0);
	}
};

// Row-oriented query batch. Optional columns use parallel `_valid` vectors
// (1 = present, 0 = NULL). When the column flag is false in QueryArrowSchema,
// the corresponding vectors are empty / unused.
struct QueryBatch {
	std::vector<std::string> read_ids;
	std::vector<std::string> sequence1;
	std::vector<std::string> sequence2;
	std::vector<int8_t> sequence2_valid;
	std::vector<std::string> qual1;
	std::vector<int8_t> qual1_valid;
	std::vector<std::string> qual2;
	std::vector<int8_t> qual2_valid;
};

// Encode a QueryBatch as an Arrow IPC stream consumable by the daemon's
// bowtie2-align input schema. Throws InternalException on encoder failure.
std::vector<uint8_t> BuildQueryIpc(const QueryBatch &qb, const QueryArrowSchema &schema_flags);

// Decode one UTINYINT[] cell (DuckDB `Value` of type LIST(UTINYINT)) into a
// Phred+33 ASCII quality string suitable for the daemon's bowtie2-align
// `qual1` / `qual2` columns. `out` is overwritten; `scratch` is a caller-
// owned reusable buffer (reused across rows to amortize allocations in the
// hot loop).
//
// Throws InvalidInputException if any element inside the list is NULL —
// quality arrays must be fully populated (mirroring the FASTQ writer's
// invariant; bowtie2 requires len(qual) == len(seq) and silent gaps would
// misalign the parser).
//
// Centralized here so align_bowtie2 and align_bowtie2_sharded share one
// canonical decoder + error message — and so this stays in lockstep with
// `EncodePhred33Quality` (the byte-level conversion) from fastq_encoder.hpp.
void DecodeListQualToPhred33(const Value &v, const char *col_name, const std::string &query_table, std::string &out,
                             std::vector<uint8_t> &scratch);

// Drain `to_emit` rows starting at `row_start` from a decoded daemon Arrow
// batch into a DuckDB DataChunk. Assumes `output` has 21 columns matching
// PopulateOutputSchema's types. Caller is responsible for SetCardinality
// before calling (we just write into the vectors).
//
// Tag widening (Int32 → BIGINT) and nullable-Utf8 decoding match the
// daemon's wire schema.
// Emit one DuckDB chunk's worth of rows from a decoded Arrow batch into
// `output`. `query_id_type` drives `read_id`; `subject_id_type` drives
// `reference` and `mate_reference`. Both must be VARCHAR, BIGINT, or UUID. For
// non-VARCHAR (BIGINT/UUID) subjects, the SAM "=" mate-reference sentinel is
// resolved to the row's reference value before invoking the codec (the literal
// "=" has no BIGINT/UUID encoding); VARCHAR output preserves "=" verbatim,
// matching pre-existing user-observable behavior. No defaults — every caller
// must explicitly commit to id types.
void EmitChunkRows(DataChunk &output, idx_t to_emit, idx_t row_start, const ArrowArray &batch,
                   const LogicalType &query_id_type, const LogicalType &subject_id_type);

// =============================================================================
// Config-JSON builder + parameter validation/mapping for the daemon's
// bowtie2-align tool.
//
// Both align_bowtie2 and align_bowtie2_sharded build the same shape of
// `config_json` per Submit. The user-facing parameter set is identical
// modulo three caller-owned knobs:
//   - align_bowtie2 owns `threads` (maps to `nthreads`); the sharded variant
//     warns about it instead.
//   - align_bowtie2_sharded owns `shard_directory`, `read_to_shard`,
//     `max_threads_per_shard`, `include_shard_name`.
// Everything else (preset, local, max_secondary, plus the ~30 tunables we
// added from the daemon's `--describe`) flows through `AppendBowtie2AlignParams`
// so there's one place to validate / map them.
// =============================================================================

// Minimal manual JSON object builder. We hand-format rather than yyjson-build
// because the schema is fixed and the cost of pulling in a full JSON writer
// for one per-batch message isn't justified. All string values are routed
// through gb::JsonEscape so the builder is safe for user-controlled strings.
struct ConfigJsonBuilder {
	std::ostringstream os;
	bool first = true;

	void append_raw(const std::string &key, const std::string &raw_value);
	void append_str(const std::string &key, const std::string &value);
	void append_int(const std::string &key, int64_t value);
	void append_bool(const std::string &key, bool value);
	std::string build() const;
};

// Param-coercion helpers. Caller name is spliced into the error message so
// users see `align_bowtie2_sharded: parameter 'preset' expects a string`
// rather than a generic miint diagnostic.
bool ValueAsBool(const char *caller, const std::string &name, const Value &v);
int64_t ValueAsInt(const char *caller, const std::string &name, const Value &v);
std::string ValueAsStr(const char *caller, const std::string &name, const Value &v);

// Allowed values for the `preset` parameter — short-form names; the
// daemon's `*-local` variants are composed by `AppendBowtie2AlignParams`
// from `preset` + `local`.
extern const std::unordered_set<std::string> kKnownPresets;

// Allowed values for the `mate_orientation` parameter.
extern const std::unordered_set<std::string> kKnownMateOrientations;

// The universe of bowtie2-align user-facing parameter names that flow to the
// daemon, minus the per-caller miint-side knobs (threads, shard_directory,
// read_to_shard, max_threads_per_shard, include_shard_name). Each caller
// validates a param is in `kCommonAlignParams ∪ <caller-specific>`.
extern const std::unordered_set<std::string> kCommonAlignParams;

// Walk `named_params` and append every recognized bowtie2-align knob to
// `cfg` with the daemon's wire-name. Caller is responsible for setting
// `index_path` and `nthreads` (per-call internals derived from the table
// function's own state, not the user's params). `caller` is spliced into
// validation error messages.
//
// Composition rules preserved from the direct-subprocess era:
//   - `local=true` + `preset=X` → daemon `preset=X-local`.
//   - `quiet=true` (miint default) → daemon `verbose=false`.
//   - `max_secondary` → daemon `k`.
//   - `local` → daemon `local_align`.
// Everything else uses the daemon's wire-name verbatim.
void AppendBowtie2AlignParams(ConfigJsonBuilder &cfg, const named_parameter_map_t &named_params, const char *caller);

// Apply LogicalType declarations for the bowtie2-align typed parameters to a
// `TableFunction::named_parameters` map. Both align_bowtie2 and
// align_bowtie2_sharded call this so the SQL surface stays in lockstep.
void RegisterBowtie2AlignNamedParameterTypes(TableFunction &tf);

// Fail-loud version gate for the bowtie2 callers. The `memory_mapped` knob (and
// the sharded mm-off default) require gpl-boundary >= 0.4.2, but the IPC
// handshake cannot report the daemon's release version (it carries only
// protocol_version + per-tool schema_version, and bowtie2's schema_version was
// not bumped for 0.4.2). So we parse the `--version` CLI of the resolved binary
// and throw IOException (naming `caller`) when it is older than 0.4.2 or its
// version cannot be determined — rather than letting an old daemon silently
// ignore `memory_mapped` and reintroduce the cold-FS regression. `binary_path`
// is the path returned by gpl_boundary::FindGplBoundary().
void RequireGplBoundaryVersion(const std::string &binary_path, const char *caller);

} // namespace bt2_daemon
} // namespace duckdb
