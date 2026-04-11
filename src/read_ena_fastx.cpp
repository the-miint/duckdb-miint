#include "read_ena_fastx.hpp"
#include "SequenceRecord.hpp"
#include "duckdb/common/vector_size.hpp"
#include <fstream>

namespace duckdb {

// ---- Data ----

ReadENAFastxTableFunction::Data::Data(std::vector<RunInfo> runs, bool include_fp, uint8_t offset)
    : runs(std::move(runs)), include_filepath(include_fp), qual_offset(offset), use_aspera(false),
      names({"sequence_index", "read_id", "comment", "sequence1", "sequence2", "qual1", "qual2", "run_accession",
             "sample_accession", "experiment_accession"}),
      types({LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
             LogicalType::VARCHAR, LogicalType::LIST(LogicalType::UTINYINT), LogicalType::LIST(LogicalType::UTINYINT),
             LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR}) {
	if (include_filepath) {
		names.emplace_back("filepath");
		types.emplace_back(LogicalType::VARCHAR);
	}
}

// ---- GlobalState ----

ReadENAFastxTableFunction::GlobalState::GlobalState(FileSystem &fs, const std::vector<RunInfo> &runs, bool use_aspera)
    : runs(runs), next_run_idx(0), fs(fs), use_aspera(use_aspera) {
	readers.resize(runs.size());
	run_sequence_counters.resize(runs.size(), 1);
}

ReadENAFastxTableFunction::GlobalState::~GlobalState() {
#if MIINT_ASPERA_SUPPORTED
	// Clean up temp file if one exists
	if (!temp_file_path.empty()) {
		std::remove(temp_file_path.c_str());
		temp_file_path.clear();
	}
	// AsperaProcess destructor handles killing child + closing pipe
#endif
}

void ReadENAFastxTableFunction::GlobalState::OpenReader(size_t run_idx) {
	if (readers[run_idx]) {
		return;
	}
	const auto &run = runs[run_idx];
	if (run.fastq_urls.empty()) {
		throw IOException("read_ena_fastx: no FASTQ URLs available for run '%s'", run.run_accession);
	}

	// Serialize file opens to avoid overwhelming remote servers with
	// simultaneous HTTP HEAD requests (e.g., EBI rejects concurrent connections).
	// Reading from opened files is still fully parallel — only the open is serialized.
	lock_guard<mutex> open_guard(open_mutex);

	auto *s1 = CreateDuckDBSeqStream(fs, run.fastq_urls[0]);
	miint::DuckDBSeqStream *s2 = nullptr;
	if (run.is_paired && run.fastq_urls.size() >= 2) {
		s2 = CreateDuckDBSeqStream(fs, run.fastq_urls[1]);
	}

	try {
		readers[run_idx] = std::make_unique<miint::SequenceReader>(s1, static_cast<miint::DuckDBSeqStream *>(s2), true);
	} catch (const std::runtime_error &e) {
		if (s2 && std::string(e.what()) == "Empty stream (sequence2)") {
			// ENA metadata may report PAIRED layout but the second FASTQ file
			// can be empty (e.g., single-end data deposited with PAIRED metadata).
			// The failed constructor consumed the original streams, so re-open
			// just the first file and fall back to single-end for this run.
			auto *s1_retry = CreateDuckDBSeqStream(fs, run.fastq_urls[0]);
			readers[run_idx] =
			    std::make_unique<miint::SequenceReader>(s1_retry, static_cast<miint::DuckDBSeqStream *>(nullptr), true);
		} else {
			throw;
		}
	}
}

#if MIINT_ASPERA_SUPPORTED

static std::string GetTempDir() {
	const char *tmpdir = getenv("TMPDIR");
	if (tmpdir && tmpdir[0] != '\0') {
		return std::string(tmpdir);
	}
	return "/tmp";
}

void ReadENAFastxTableFunction::GlobalState::OpenReaderAspera(size_t run_idx) {
	if (readers[run_idx]) {
		return;
	}
	const auto &run = runs[run_idx];

	std::string filename;
	size_t file_size;

	if (!aspera_process->NextFile(filename, file_size)) {
		throw IOException("read_ena_fastx: Aspera stream ended unexpectedly (expected files for run '%s')",
		                  run.run_accession);
	}

	// Determine gzip from tar filename, fall back to aspera path if available
	std::string gz_hint = filename;
	if (gz_hint.empty() && !run.aspera_paths.empty()) {
		gz_hint = run.aspera_paths[0].remote_path;
	}
	bool is_gz = IsGzipped(gz_hint);

	if (!run.is_paired) {
		auto *s1 = miint::CreateAsperaSeqStream(aspera_process.get(), is_gz);
		readers[run_idx] =
		    std::make_unique<miint::SequenceReader>(s1, static_cast<miint::AsperaSeqStream *>(nullptr), true);
	} else {
		// Paired-end: stream R1 from pipe to temp file in chunks, then stream R2 live
		static constexpr size_t COPY_BUF_SIZE = 1024 * 1024; // 1 MB

		std::string tmpl = GetTempDir() + "/miint_aspera_r1_XXXXXX";
		std::vector<char> tmpl_buf(tmpl.begin(), tmpl.end());
		tmpl_buf.push_back('\0');
		int fd = mkstemp(tmpl_buf.data());
		if (fd == -1) {
			throw IOException("read_ena_fastx: failed to create temp file for paired-end Aspera buffering");
		}
		temp_file_path = std::string(tmpl_buf.data());

		// Stream-copy R1 from pipe to temp file in 1MB chunks (no full-file buffering)
		std::vector<char> copy_buf(COPY_BUF_SIZE);
		while (true) {
			int n = aspera_process->ReadBounded(copy_buf.data(), COPY_BUF_SIZE);
			if (n <= 0) {
				break;
			}
			size_t written = 0;
			while (written < static_cast<size_t>(n)) {
				ssize_t w = write(fd, copy_buf.data() + written, static_cast<size_t>(n) - written);
				if (w <= 0) {
					close(fd);
					std::remove(temp_file_path.c_str());
					temp_file_path.clear();
					throw IOException("read_ena_fastx: failed to write temp file for paired-end Aspera buffering");
				}
				written += static_cast<size_t>(w);
			}
		}
		close(fd);

		// Advance tar stream to R2
		std::string filename2;
		size_t file_size2;
		if (!aspera_process->NextFile(filename2, file_size2)) {
			std::remove(temp_file_path.c_str());
			temp_file_path.clear();
			throw IOException("read_ena_fastx: Aspera stream ended unexpectedly (expected R2 for run '%s')",
			                  run.run_accession);
		}

		std::string gz_hint2 = filename2;
		if (gz_hint2.empty() && run.aspera_paths.size() >= 2) {
			gz_hint2 = run.aspera_paths[1].remote_path;
		}
		bool is_gz2 = IsGzipped(gz_hint2);

		auto *s1 = CreateDuckDBSeqStream(fs, temp_file_path, is_gz);
		auto *s2 = miint::CreateAsperaSeqStream(aspera_process.get(), is_gz2);
		readers[run_idx] = std::make_unique<miint::SequenceReader>(s1, s2, true);
	}
}
#endif // MIINT_ASPERA_SUPPORTED

// ---- Bind ----

// Resolve accessions to run info via ENA Portal API
static std::vector<ReadENAFastxTableFunction::RunInfo> ResolveRuns(miint::ENAClient &client,
                                                                   const std::vector<std::string> &accessions) {
	std::vector<ReadENAFastxTableFunction::RunInfo> runs;

	for (const auto &acc : accessions) {
		auto tsv =
		    client.Search(acc, "read_run",
		                  "run_accession,sample_accession,experiment_accession,fastq_ftp,fastq_aspera,library_layout");
		auto parsed = miint::ENAParser::ParseTSV(tsv);

		// Find column indices by name
		int run_col = -1, sample_col = -1, exp_col = -1, ftp_col = -1, aspera_col = -1, layout_col = -1;
		for (size_t i = 0; i < parsed.column_names.size(); i++) {
			if (parsed.column_names[i] == "run_accession")
				run_col = static_cast<int>(i);
			else if (parsed.column_names[i] == "sample_accession")
				sample_col = static_cast<int>(i);
			else if (parsed.column_names[i] == "experiment_accession")
				exp_col = static_cast<int>(i);
			else if (parsed.column_names[i] == "fastq_ftp")
				ftp_col = static_cast<int>(i);
			else if (parsed.column_names[i] == "fastq_aspera")
				aspera_col = static_cast<int>(i);
			else if (parsed.column_names[i] == "library_layout")
				layout_col = static_cast<int>(i);
		}

		for (const auto &row : parsed.rows) {
			ReadENAFastxTableFunction::RunInfo info;
			info.run_accession = (run_col >= 0 && run_col < (int)row.size()) ? row[run_col] : "";
			info.sample_accession = (sample_col >= 0 && sample_col < (int)row.size()) ? row[sample_col] : "";
			info.experiment_accession = (exp_col >= 0 && exp_col < (int)row.size()) ? row[exp_col] : "";

			std::string ftp_field = (ftp_col >= 0 && ftp_col < (int)row.size()) ? row[ftp_col] : "";
			info.fastq_urls = miint::ENAParser::FTPtoHTTPS(ftp_field);

			std::string aspera_field = (aspera_col >= 0 && aspera_col < (int)row.size()) ? row[aspera_col] : "";
			info.aspera_paths = miint::AsperaUtils::ParseAsperaPaths(aspera_field);

			std::string layout = (layout_col >= 0 && layout_col < (int)row.size()) ? row[layout_col] : "";
			info.is_paired = (layout == "PAIRED");

			if (!info.fastq_urls.empty()) {
				runs.push_back(std::move(info));
			}
		}
	}

	return runs;
}

unique_ptr<FunctionData> ReadENAFastxTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                         vector<LogicalType> &return_types,
                                                         vector<std::string> &names) {
	std::vector<std::string> accessions;
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		accessions.push_back(input.inputs[0].ToString());
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			accessions.push_back(child.ToString());
		}
	} else {
		throw InvalidInputException("read_ena_fastx: first argument must be VARCHAR or VARCHAR[]");
	}

	if (accessions.empty()) {
		throw InvalidInputException("read_ena_fastx: at least one accession must be provided");
	}
	for (auto &acc : accessions) {
		// Trim whitespace
		auto start = acc.find_first_not_of(" \t\n\r");
		if (start == std::string::npos) {
			acc.clear();
		} else {
			acc = acc.substr(start, acc.find_last_not_of(" \t\n\r") - start + 1);
		}
		if (acc.empty()) {
			throw InvalidInputException("read_ena_fastx: accession cannot be empty");
		}
	}

	bool include_filepath = false;
	auto fp_param = input.named_parameters.find("include_filepath");
	if (fp_param != input.named_parameters.end() && !fp_param->second.IsNull()) {
		include_filepath = fp_param->second.GetValue<bool>();
	}

	uint8_t qual_offset = 33;
	auto qo_param = input.named_parameters.find("qual_offset");
	if (qo_param != input.named_parameters.end() && !qo_param->second.IsNull()) {
		qual_offset = static_cast<uint8_t>(qo_param->second.GetValue<int64_t>());
	}

	std::string download_method = "auto";
	auto dm_param = input.named_parameters.find("download_method");
	if (dm_param != input.named_parameters.end() && !dm_param->second.IsNull()) {
		download_method = dm_param->second.ToString();
		if (download_method != "auto" && download_method != "aspera" && download_method != "http") {
			throw InvalidInputException(
			    "read_ena_fastx: download_method must be 'auto', 'aspera', or 'http' (got '%s')", download_method);
		}
	}

	// Resolve accessions to run info via ENA Portal API
	auto &db = DatabaseInstance::GetDatabase(context);
	miint::ENAClient client(db);
	auto runs = ResolveRuns(client, accessions);

	if (runs.empty()) {
		throw IOException("read_ena_fastx: no FASTQ data found for the provided accession(s)");
	}

	auto data = make_uniq<Data>(std::move(runs), include_filepath, qual_offset);

	// Determine whether to use Aspera
	data->use_aspera = false;
#if MIINT_ASPERA_SUPPORTED
	if (download_method == "aspera" || download_method == "auto") {
		bool aspera_available = miint::AsperaUtils::IsAvailable();
		if (download_method == "aspera" && !aspera_available) {
			throw IOException("read_ena_fastx: download_method='aspera' but ascp not found in PATH. "
			                  "Install IBM Aspera CLI or use download_method='http'");
		}
		if (aspera_available) {
			// Check that ALL runs have aspera paths (required for single-process stdio-tar)
			bool all_have_aspera = true;
			for (const auto &run : data->runs) {
				if (run.aspera_paths.empty()) {
					all_have_aspera = false;
					break;
				}
			}
			if (download_method == "aspera" && !all_have_aspera) {
				throw IOException("read_ena_fastx: download_method='aspera' but not all runs have Aspera paths");
			}
			if (all_have_aspera) {
				std::string ascp_path = miint::AsperaUtils::FindAscp();
				std::string key_path = miint::AsperaUtils::ResolveKey(db, ascp_path, download_method == "aspera");
				if (!key_path.empty()) {
					data->use_aspera = true;
					data->aspera_config = miint::AsperaUtils::BuildConfig(ascp_path, key_path);
				}
			}
		}
	}
#else
	if (download_method == "aspera") {
		throw IOException("read_ena_fastx: Aspera is not supported on this platform");
	}
#endif

	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

// ---- InitGlobal / InitLocal ----

unique_ptr<GlobalTableFunctionState> ReadENAFastxTableFunction::InitGlobal(ClientContext &context,
                                                                           TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	FileSystem &fs = FileSystem::GetFileSystem(context);
	auto state = make_uniq<GlobalState>(fs, data.runs, data.use_aspera);

#if MIINT_ASPERA_SUPPORTED
	if (data.use_aspera) {
		// Collect all remote paths across all runs (ordered: for each run, R1 then R2)
		std::vector<std::string> remote_paths;
		std::string host;
		for (const auto &run : data.runs) {
			for (const auto &ap : run.aspera_paths) {
				if (host.empty()) {
					host = ap.host;
				}
				remote_paths.push_back(ap.remote_path);
			}
		}

		auto config = data.aspera_config;
		if (!host.empty()) {
			config.host = host;
		}
		state->aspera_process = std::make_unique<miint::AsperaProcess>(config, remote_paths);
	}
#endif

	return state;
}

unique_ptr<LocalTableFunctionState> ReadENAFastxTableFunction::InitLocal(ExecutionContext &context,
                                                                         TableFunctionInitInput &input,
                                                                         GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// ---- Execute ----

void ReadENAFastxTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	miint::SequenceRecordBatch batch;
	size_t current_run_idx;

	// Claim-read-release loop (same pattern as read_fastx)
	while (true) {
		if (!local_state.has_run) {
			{
				lock_guard<mutex> read_lock(global_state.lock);
				if (global_state.next_run_idx >= global_state.runs.size()) {
#if MIINT_ASPERA_SUPPORTED
					// All runs consumed — check ascp exit code
					if (global_state.aspera_process) {
						int exit_code = global_state.aspera_process->WaitForExit();
						global_state.aspera_process.reset();
						if (exit_code != 0 && exit_code != -1) {
							throw IOException("read_ena_fastx: ascp exited with code %d", exit_code);
						}
					}
#endif
					output.SetCardinality(0);
					return;
				}
				local_state.current_run_idx = global_state.next_run_idx;
				global_state.next_run_idx++;
				local_state.has_run = true;
			} // Lock released before I/O

			// Open reader without holding the lock — safe because this thread
			// has exclusive ownership of the claimed run index
#if MIINT_ASPERA_SUPPORTED
			if (global_state.use_aspera) {
				global_state.OpenReaderAspera(local_state.current_run_idx);
			} else {
				global_state.OpenReader(local_state.current_run_idx);
			}
#else
			global_state.OpenReader(local_state.current_run_idx);
#endif
		}

		current_run_idx = local_state.current_run_idx;
		batch = global_state.readers[current_run_idx]->read(STANDARD_VECTOR_SIZE);

		if (batch.empty()) {
			// Release the reader to free the HTTP connection (or Aspera pipe slot)
			global_state.readers[current_run_idx].reset();
#if MIINT_ASPERA_SUPPORTED
			// Clean up temp file for paired-end Aspera
			if (!global_state.temp_file_path.empty()) {
				std::remove(global_state.temp_file_path.c_str());
				global_state.temp_file_path.clear();
			}
#endif
			local_state.has_run = false;
			continue;
		}
		break;
	}

	idx_t count = batch.size();
	idx_t field_idx = 0;
	const auto &run = global_state.runs[current_run_idx];
	auto &seq_counter = global_state.run_sequence_counters[current_run_idx];

	// sequence_index (column 0)
	auto seq_idx_data = FlatVector::GetData<int64_t>(output.data[field_idx++]);
	for (idx_t i = 0; i < count; i++) {
		seq_idx_data[i] = static_cast<int64_t>(seq_counter++);
	}

	// read_id (column 1)
	SetResultVectorString(output.data[field_idx++], batch.read_ids);

	// comment (column 2)
	SetResultVectorStringNullable(output.data[field_idx++], batch.comments);

	// sequence1 (column 3)
	SetResultVectorString(output.data[field_idx++], batch.sequences1);

	// sequence2 (column 4)
	if (batch.is_paired) {
		SetResultVectorString(output.data[field_idx++], batch.sequences2);
	} else {
		SetResultVectorNull(output.data[field_idx++]);
	}

	// qual1 (column 5)
	SetResultVectorListUInt8(output.data[field_idx++], batch.quals1, bind_data.qual_offset);

	// qual2 (column 6)
	if (batch.is_paired) {
		SetResultVectorListUInt8(output.data[field_idx++], batch.quals2, bind_data.qual_offset);
	} else {
		SetResultVectorNull(output.data[field_idx++]);
	}

	// run_accession (column 7) - constant for all rows in this chunk
	SetResultVectorFilepath(output.data[field_idx++], run.run_accession);

	// sample_accession (column 8) - constant
	SetResultVectorFilepath(output.data[field_idx++], run.sample_accession);

	// experiment_accession (column 9) - constant
	SetResultVectorFilepath(output.data[field_idx++], run.experiment_accession);

	// filepath (column 10) - optional, semicolon-separated for paired-end
	if (bind_data.include_filepath) {
		std::string filepath;
		if (bind_data.use_aspera && !run.aspera_paths.empty()) {
			filepath = run.aspera_paths[0].host + ":" + run.aspera_paths[0].remote_path;
			for (size_t u = 1; u < run.aspera_paths.size(); u++) {
				filepath += ";" + run.aspera_paths[u].host + ":" + run.aspera_paths[u].remote_path;
			}
		} else {
			filepath = run.fastq_urls[0];
			for (size_t u = 1; u < run.fastq_urls.size(); u++) {
				filepath += ";" + run.fastq_urls[u];
			}
		}
		SetResultVectorFilepath(output.data[field_idx++], filepath);
	}

	output.SetCardinality(count);
}

// ---- Registration ----

TableFunction ReadENAFastxTableFunction::GetFunction() {
	auto tf = TableFunction("read_ena_fastx", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	tf.named_parameters["qual_offset"] = LogicalType::BIGINT;
	tf.named_parameters["download_method"] = LogicalType::VARCHAR;
	return tf;
}

void ReadENAFastxTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
