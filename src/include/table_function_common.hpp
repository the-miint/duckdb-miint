#pragma once
#include <string>
#include <vector>
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/named_parameter_map.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/storage/statistics/node_statistics.hpp"
#include "SpectrumBatch.hpp"
#include "QualScore.hpp"

namespace duckdb {

// Check if a path refers to stdin (-, /dev/stdin, /dev/fd/0, /proc/self/fd/0)
bool IsStdinPath(const std::string &path);

// Check if a path has .gz extension (for gzip-compressed files)
bool IsGzipped(const std::string &path);

// Parse include_filepath named parameter (optional BOOLEAN, default false)
bool ParseIncludeFilepathParameter(const named_argument_map_t &named_parameters);

// Row-count estimate for a table function that streams files, so the query optimizer knows
// roughly how big the scan is. Without one, DuckDB assumes a tiny scan and can pick it as a
// join's build side -- hashing a multi-GB file into memory instead of streaming it.
//
// It errs high: an over-estimate makes the optimizer probe the scan, always the
// safe choice for a streaming source; an under-estimate makes it build a hash table on it.
// Heuristic: bytes on disk / `bytes_per_row`, with gzip counted at 4x its size. A remote or
// stdin path, or one whose size cannot be read, gives a large constant instead.
unique_ptr<NodeStatistics> EstimateFileScanCardinality(const std::vector<std::string> &paths, idx_t bytes_per_row);

// Expand a glob pattern into a sorted list of file paths
// - If pattern contains glob characters (*, ?, []), expands and sorts alphabetically
// - If pattern is a literal path, returns it as DuckDB's file system expands it, as its glob does for
//   the paths it returns: a file: URI becomes the local path it names, and a leading ~ the home
//   directory. Every other path, remote URIs and stdin included, is returned as given.
// - Throws IOException if glob matches zero files
// - Throws InvalidInputException if pattern could match stdin paths
std::vector<std::string> ExpandGlobPattern(FileSystem &fs, ClientContext &context, const std::string &pattern);

// The paths of a VARCHAR[] argument: literal paths (no glob expansion), each expanded as a literal
// path is above. Throws InvalidInputException, naming function_name (and list_name for a NULL), for a
// NULL element or an empty list.
std::vector<std::string> ExpandPathList(FileSystem &fs, const Value &list, const std::string &function_name,
                                        const std::string &list_name = "file path list");

// Result of glob expansion with metadata
struct GlobExpansionResult {
	std::vector<std::string> paths;
	bool is_glob; // Whether the original pattern was a glob
};

// Expand glob pattern and return whether it was a glob pattern
// Useful for paired-end read validation where both must be globs or both literals
GlobExpansionResult ExpandGlobPatternWithInfo(FileSystem &fs, ClientContext &context, const std::string &pattern);

// A reader's first argument: a VARCHAR path or glob (expanded as ExpandGlobPatternWithInfo does), or
// a VARCHAR[] of literal paths (as ExpandPathList does; is_glob is false). Throws
// InvalidInputException, naming function_name, for NULL or any other type.
GlobExpansionResult ExpandPathArgument(FileSystem &fs, ClientContext &context, const Value &arg,
                                       const std::string &function_name);

// Throws IOException("File not found: <path>") for the first local path that does not exist. Stdin
// and remote paths are not checked: they are opened when the scan starts.
void RequireLocalFilesExist(FileSystem &fs, const std::vector<std::string> &paths);

// The columns read_fastx returns, which every sequence reader returns first: sequence_index,
// read_id, comment, sequence1, sequence2, qual1 and qual2.
std::vector<std::string> SequenceColumnNames();
std::vector<LogicalType> SequenceColumnTypes();

// A BIGINT extension option, registered in LoadInternal, that must not be negative: `if_null` when
// it is NULL. `zero_means` says, in the error for a negative value, what 0 does.
idx_t GetNonNegativeOption(ClientContext &context, const char *name, idx_t if_null, const char *zero_means);

// --- Result vector helpers ---
// Shared utilities for populating DuckDB result vectors from batch data.
// Used by read_fastx, read_alignments, read_sequences_sam, and other table functions.

void SetResultVectorNull(Vector &result_vector);
void SetResultVectorString(Vector &result_vector, const std::vector<std::string> &values);
void SetResultVectorStringNullable(Vector &result_vector, const std::vector<std::string> &values);
void SetResultVectorFilepath(Vector &result_vector, const std::string &filepath);
void SetResultVectorUInt8(Vector &result_vector, const std::vector<uint8_t> &values);
void SetResultVectorUInt16(Vector &result_vector, const std::vector<uint16_t> &values);
void SetResultVectorInt64(Vector &result_vector, const std::vector<int64_t> &values);
void SetResultVectorInt64Nullable(Vector &result_vector, const std::vector<int64_t> &values,
                                  const std::vector<bool> &valid);
void SetResultVectorListUInt8(Vector &result_vector, const std::vector<miint::QualScore> &values, uint8_t qual_offset);

// Locate the child slice for one row of a LIST(UTINYINT) input: writes the
// pointer to the row's first byte and the row's length. The pointer is valid
// for the lifetime of the underlying ListVector child buffer; the caller must
// not free it.
void GetListUInt8Slice(Vector &list_vec, UnifiedVectorFormat &list_data, idx_t row_idx, const uint8_t *&out_data,
                       idx_t &out_length);
// Mark one row of a STRUCT-returning function's result NULL, children included.
//
// FlatVector::Validity(result).SetInvalid(row_idx) is NOT sufficient for a STRUCT. It
// marks only the struct's own validity, while struct_extract hands back a bare
// reference to the child without applying its parent's validity -- so `f(x) IS NULL`
// reads true while `(f(x)).some_field` returns whatever was left in the child buffer.
// DuckDB treats "a NULL struct entry implies NULL children" as an invariant and
// asserts it in debug builds, so the shortcut is a debug-build abort as well as a
// wrong release-build answer.
//
// This is a one-line forward to duckdb::FlatVector::SetNull(result, row_idx, true), which
// already recurses into STRUCT and ARRAY children. It exists to give the rule one name and
// somewhere to record the reasoning, NOT because core lacks the behaviour -- a caller that
// would need a heavyweight include to reach this header (see pairwise_align_shared.hpp)
// calls FlatVector::SetNull directly and is equally correct.
//
// It does NOT reset a LIST child's list_entry_t. That is not a correctness gap, because a
// NULL row's list_entry is not read through the null parent, but it is why several callers
// still zero their own list offsets alongside this call -- they are not compensating for a
// deficiency here.
//
// Only a BARE STRUCT return needs this. A LIST(STRUCT) return does not leak: measured on
// cigar_query_intervals, sequence_split, compress_intervals and cumulative_coverage, a NULL
// list row gives NULL from len(), from list_extract(), and from (list_extract(...)).field.
// The full inventory and the evidence are in test/sql/struct_null_children.test.
void SetStructRowNull(Vector &result, idx_t row_idx);

void SetResultVectorInt32(Vector &result_vector, const std::vector<int32_t> &values);
void SetResultVectorInt32Nullable(Vector &result_vector, const std::vector<int32_t> &values,
                                  const std::vector<bool> &valid);
void SetResultVectorDouble(Vector &result_vector, const std::vector<double> &values);
void SetResultVectorDoubleNullable(Vector &result_vector, const std::vector<double> &values,
                                   const std::vector<bool> &valid);
void SetResultVectorListDouble(Vector &result_vector, const std::vector<std::vector<double>> &values);

// Populate a DataChunk from a SpectrumBatch (shared by read_mzml and read_mzxml).
// Maps the 27-column schema from MzMLSpectrumBatch to output vectors.
void PopulateSpectrumBatchOutput(DataChunk &output, const miint::MzMLSpectrumBatch &batch, bool include_filepath,
                                 const std::string &filepath);

} // namespace duckdb
