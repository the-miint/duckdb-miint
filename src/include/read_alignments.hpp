#pragma once
#include "SAMReader.hpp"
#include "SAMRecord.hpp"
#include "remote_file_helper.hpp"
#include "hfile_duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <optional>
#include <thread>
#include <unordered_map>
#include <vector>

namespace duckdb {
class ReadAlignmentsTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> sam_paths;
		std::optional<std::string> reference_lengths_table;
		bool include_filepath;
		bool include_seq_qual;

		std::vector<std::string> names;
		std::vector<LogicalType> types;
		std::vector<miint::SAMRecordField> fields;

		explicit Data(const std::vector<std::string> &paths, const std::optional<std::string> &ref_table,
		              bool include_fp, bool include_sq)
		    : sam_paths(paths), reference_lengths_table(ref_table), include_filepath(include_fp),
		      include_seq_qual(include_sq),
		      names({"read_id", "flags",          "reference",     "position",        "stop_position", "mapq",
		             "cigar",   "mate_reference", "mate_position", "template_length", "tag_as",        "tag_xs",
		             "tag_ys",  "tag_xn",         "tag_xm",        "tag_xo",          "tag_xg",        "tag_nm",
		             "tag_yt",  "tag_md",         "tag_sa"}),
		      types({LogicalType::VARCHAR,   // read_id
		             LogicalType::USMALLINT, // flags
		             LogicalType::VARCHAR,   // reference
		             LogicalType::BIGINT,    // position
		             LogicalType::BIGINT,    // stop_position
		             LogicalType::UTINYINT,  // mapq
		             LogicalType::VARCHAR,   // cigar
		             LogicalType::VARCHAR,   // mate_reference
		             LogicalType::BIGINT,    // mate_position
		             LogicalType::BIGINT,    // template_length
		             LogicalType::BIGINT,    // tag_as
		             LogicalType::BIGINT,    // tag_xs
		             LogicalType::BIGINT,    // tag_ys
		             LogicalType::BIGINT,    // tag_xn
		             LogicalType::BIGINT,    // tag_xm
		             LogicalType::BIGINT,    // tag_xo
		             LogicalType::BIGINT,    // tag_xg
		             LogicalType::BIGINT,    // tag_nm
		             LogicalType::VARCHAR,   // tag_yt
		             LogicalType::VARCHAR,   // tag_md
		             LogicalType::VARCHAR}), // tag_sa
		      fields({miint::SAMRecordField::READ_ID,       miint::SAMRecordField::FLAGS,
		              miint::SAMRecordField::REFERENCE,     miint::SAMRecordField::POSITION,
		              miint::SAMRecordField::STOP_POSITION, miint::SAMRecordField::MAPQ,
		              miint::SAMRecordField::CIGAR,         miint::SAMRecordField::MATE_REFERENCE,
		              miint::SAMRecordField::MATE_POSITION, miint::SAMRecordField::TEMPLATE_LENGTH,
		              miint::SAMRecordField::TAG_AS,        miint::SAMRecordField::TAG_XS,
		              miint::SAMRecordField::TAG_YS,        miint::SAMRecordField::TAG_XN,
		              miint::SAMRecordField::TAG_XM,        miint::SAMRecordField::TAG_XO,
		              miint::SAMRecordField::TAG_XG,        miint::SAMRecordField::TAG_NM,
		              miint::SAMRecordField::TAG_YT,        miint::SAMRecordField::TAG_MD,
		              miint::SAMRecordField::TAG_SA}) {
			if (include_seq_qual) {
				names.emplace_back("sequence");
				types.emplace_back(LogicalType::VARCHAR);
				names.emplace_back("qual");
				types.emplace_back(LogicalType::LIST(LogicalType::UTINYINT));
			}
			if (include_filepath) {
				names.emplace_back("filepath");
				types.emplace_back(LogicalType::VARCHAR);
			}
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		// Readers are created lazily when a thread claims a file, not all upfront.
		// This avoids opening many HTTP connections simultaneously for remote file arrays.
		std::vector<std::unique_ptr<miint::SAMReader>> readers;
		std::vector<std::string> filepaths; // Original paths (for include_filepath)
		size_t next_file_idx;
		FileSystem &fs;
		std::optional<std::unordered_map<std::string, uint64_t>> ref_lengths;
		bool include_seq_qual;

		idx_t MaxThreads() const override {
			auto hw = std::thread::hardware_concurrency();
			if (hw == 0) {
				hw = 1;
			}
			return std::min<idx_t>(filepaths.size(), hw);
		}

		GlobalState(const std::vector<std::string> &paths, FileSystem &fs,
		            std::optional<std::unordered_map<std::string, uint64_t>> ref_lengths, bool include_seq_qual)
		    : filepaths(paths), next_file_idx(0), fs(fs), ref_lengths(std::move(ref_lengths)),
		      include_seq_qual(include_seq_qual) {
			// Pre-allocate slots but don't open connections yet
			readers.resize(paths.size());
		}

		// Open reader for a specific file index. Called under lock when a thread claims the file.
		void OpenReader(size_t file_idx) {
			if (readers[file_idx]) {
				return; // Already opened
			}
			const auto &path = filepaths[file_idx];
			try {
				if (miint::RemoteFileHelper::IsRemotePath(path)) {
					hFILE *hf = miint::hfile_duckdb_open(fs, path);
					if (!hf) {
						throw IOException("Failed to open remote file: " + path);
					}
					if (ref_lengths.has_value()) {
						readers[file_idx] =
						    std::make_unique<miint::SAMReader>(hf, path, ref_lengths.value(), include_seq_qual);
					} else {
						readers[file_idx] = std::make_unique<miint::SAMReader>(hf, path, include_seq_qual);
					}
				} else {
					if (ref_lengths.has_value()) {
						readers[file_idx] =
						    std::make_unique<miint::SAMReader>(path, ref_lengths.value(), include_seq_qual);
					} else {
						readers[file_idx] = std::make_unique<miint::SAMReader>(path, include_seq_qual);
					}
				}
			} catch (std::exception &e) {
				throw IOException("Error opening '%s': %s", path, e.what());
			}
		}
	};

	struct LocalState : public LocalTableFunctionState {
		size_t current_file_idx;
		bool has_file;

		LocalState() : current_file_idx(0), has_file(false) {
		}
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};
}; // namespace duckdb
