#include "read_ena_sequences.hpp"
#include "SequenceRecord.hpp"
#include "duckdb/common/printer.hpp"
#include "duckdb/common/vector_size.hpp"
#include <cerrno>
#include <fstream>

namespace duckdb {

// ---- Data ----

ReadENASequencesTableFunction::Data::Data(std::vector<RunInfo> runs, bool include_fp, uint8_t offset, bool trim,
                                          const std::string &prefer_format)
    : runs(std::move(runs)), include_filepath(include_fp), qual_offset(offset), use_aspera(false), trim(trim),
      prefer_format(prefer_format), names({"sequence_index", "read_id", "comment", "sequence1", "sequence2", "qual1",
                                           "qual2", "run_accession", "sample_accession", "experiment_accession"}),
      types({LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
             LogicalType::VARCHAR, LogicalType::LIST(LogicalType::UTINYINT), LogicalType::LIST(LogicalType::UTINYINT),
             LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR}) {
	if (include_filepath) {
		names.emplace_back("filepath");
		types.emplace_back(LogicalType::VARCHAR);
	}
}

// ---- GlobalState ----

ReadENASequencesTableFunction::GlobalState::GlobalState(FileSystem &fs, const std::vector<RunInfo> &runs,
                                                        bool use_aspera, bool trim)
    : runs(runs), next_run_idx(0), fs(fs), use_aspera(use_aspera), trim(trim), runs_completed(0),
      total_runs(runs.size()), bytes_completed(0), skipped_warned(false) {
	uint64_t sum_bytes = 0;
	for (const auto &r : runs) {
		sum_bytes += r.total_bytes;
	}
	total_bytes = sum_bytes;
	readers.resize(runs.size());
	sff_readers.resize(runs.size());
	sff_temp_paths.resize(runs.size());
	run_sequence_counters.resize(runs.size(), 1);
#if MIINT_ASPERA_SUPPORTED
	aspera_processes.resize(runs.size());
	temp_file_paths.resize(runs.size());
#endif
}

ReadENASequencesTableFunction::GlobalState::~GlobalState() {
	// Clean up SFF temp files
	for (auto &path : sff_temp_paths) {
		if (!path.empty()) {
			std::remove(path.c_str());
			path.clear();
		}
	}
#if MIINT_ASPERA_SUPPORTED
	for (auto &path : temp_file_paths) {
		if (!path.empty()) {
			std::remove(path.c_str());
			path.clear();
		}
	}
	// AsperaProcess destructors handle killing children + closing pipes
#endif
}

void ReadENASequencesTableFunction::GlobalState::OpenReader(size_t run_idx) {
	if (readers[run_idx]) {
		return;
	}
	const auto &run = runs[run_idx];
	if (run.fastq_urls.empty()) {
		throw IOException("read_ena_sequences: no FASTQ URLs available for run '%s'", run.run_accession);
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

void ReadENASequencesTableFunction::GlobalState::OpenReaderSFF(size_t run_idx) {
	if (sff_readers[run_idx]) {
		return;
	}
	const auto &run = runs[run_idx];
	if (run.sff_url.empty()) {
		throw IOException("read_ena_sequences: no SFF URL available for run '%s'", run.run_accession);
	}

	// Create temp file for the SFF download
	auto temp_dir = miint::GetTempDir();
	auto temp_path = temp_dir + "/miint_ena_sff_XXXXXX";
	std::vector<char> tmpl_buf(temp_path.begin(), temp_path.end());
	tmpl_buf.push_back('\0');
	int fd = mkstemp(tmpl_buf.data());
	if (fd == -1) {
		throw IOException("read_ena_sequences: failed to create temp file for SFF download");
	}
	temp_path = std::string(tmpl_buf.data());
	sff_temp_paths[run_idx] = temp_path;

	// RAII guard: ensures fd is closed exactly once and temp file is removed on error
	bool fd_closed = false;
	auto cleanup_on_error = [&]() {
		if (!fd_closed) {
			close(fd);
			fd_closed = true;
		}
		std::remove(temp_path.c_str());
		sff_temp_paths[run_idx].clear();
	};

	try {
		// Serialize only the HTTP connection setup, not the full download.
		// This avoids blocking all threads for the duration of multi-MB downloads
		// while still rate-limiting simultaneous connection opens to remote servers.
		unique_ptr<FileHandle> source;
		{
			lock_guard<mutex> open_guard(open_mutex);
			source = fs.OpenFile(run.sff_url, FileOpenFlags(FileOpenFlags::FILE_FLAGS_READ));
		}

		static constexpr size_t DOWNLOAD_BUF_SIZE = 1048576; // 1MB
		std::vector<char> buf(DOWNLOAD_BUF_SIZE);
		while (true) {
			auto bytes_read = source->Read(buf.data(), DOWNLOAD_BUF_SIZE);
			if (bytes_read == 0) {
				break;
			}
			size_t written = 0;
			while (written < bytes_read) {
				ssize_t w = write(fd, buf.data() + written, bytes_read - written);
				if (w < 0) {
					if (errno == EINTR) {
						continue;
					}
					int saved_errno = errno;
					cleanup_on_error();
					throw IOException("read_ena_sequences: failed writing SFF temp file for run '%s' (errno=%d: %s)",
					                  run.run_accession, saved_errno, strerror(saved_errno));
				}
				written += static_cast<size_t>(w);
			}
		}
		close(fd);
		fd_closed = true;
	} catch (...) {
		cleanup_on_error();
		throw;
	}

	sff_readers[run_idx] = std::make_unique<miint::SFFReader>(temp_path, trim);
}

#if MIINT_ASPERA_SUPPORTED

void ReadENASequencesTableFunction::GlobalState::OpenReaderAspera(size_t run_idx) {
	if (readers[run_idx]) {
		return;
	}
	const auto &run = runs[run_idx];

	// Collect remote paths for this run only
	std::vector<std::string> remote_paths;
	for (const auto &ap : run.aspera_paths) {
		remote_paths.push_back(ap.remote_path);
	}

	auto config = aspera_config;
	if (!run.aspera_paths.empty()) {
		config.host = run.aspera_paths[0].host;
	}

	// Serialize ascp launches (fork/exec) — pipe reads happen on each thread's own pipe
	{
		lock_guard<mutex> open_guard(open_mutex);
		aspera_processes[run_idx] = std::make_unique<miint::AsperaProcess>(config, remote_paths);
	}
	auto *proc = aspera_processes[run_idx].get();

	std::string filename;
	size_t file_size;

	if (!proc->NextFile(filename, file_size)) {
		throw IOException("read_ena_sequences: Aspera stream ended unexpectedly (expected files for run '%s')",
		                  run.run_accession);
	}

	// Determine gzip from tar filename, fall back to aspera path if available
	std::string gz_hint = filename;
	if (gz_hint.empty() && !run.aspera_paths.empty()) {
		gz_hint = run.aspera_paths[0].remote_path;
	}
	bool is_gz = IsGzipped(gz_hint);

	if (!run.is_paired) {
		auto *s1 = miint::CreateAsperaSeqStream(proc, is_gz);
		try {
			readers[run_idx] =
			    std::make_unique<miint::SequenceReader>(s1, static_cast<miint::AsperaSeqStream *>(nullptr), true);
		} catch (...) {
			delete s1;
			throw;
		}
	} else {
		// Paired-end: stream R1 from pipe to temp file in chunks, then stream R2 live
		static constexpr size_t COPY_BUF_SIZE = 1024 * 1024; // 1 MB

		std::string tmpl = miint::GetTempDir() + "/miint_aspera_r1_XXXXXX";
		std::vector<char> tmpl_buf(tmpl.begin(), tmpl.end());
		tmpl_buf.push_back('\0');
		int fd = mkstemp(tmpl_buf.data());
		if (fd == -1) {
			throw IOException("read_ena_sequences: failed to create temp file for paired-end Aspera buffering");
		}
		temp_file_paths[run_idx] = std::string(tmpl_buf.data());

		// Stream-copy R1 from pipe to temp file in 1MB chunks (no full-file buffering)
		std::vector<char> copy_buf(COPY_BUF_SIZE);
		while (true) {
			int n = proc->ReadBounded(copy_buf.data(), COPY_BUF_SIZE);
			if (n == 0) {
				break;
			}
			if (n < 0) {
				close(fd);
				std::remove(temp_file_paths[run_idx].c_str());
				temp_file_paths[run_idx].clear();
				throw IOException("read_ena_sequences: pipe read error while buffering R1 for run '%s'",
				                  run.run_accession);
			}
			size_t written = 0;
			while (written < static_cast<size_t>(n)) {
				ssize_t w = write(fd, copy_buf.data() + written, static_cast<size_t>(n) - written);
				if (w < 0) {
					if (errno == EINTR) {
						continue;
					}
					int err = errno;
					close(fd);
					std::remove(temp_file_paths[run_idx].c_str());
					temp_file_paths[run_idx].clear();
					throw IOException("read_ena_sequences: failed to write temp file for paired-end Aspera buffering "
					                  "(run '%s', errno=%d: %s)",
					                  run.run_accession, err, strerror(err));
				}
				written += static_cast<size_t>(w);
			}
		}
		close(fd);

		// Advance tar stream to R2
		std::string filename2;
		size_t file_size2;
		if (!proc->NextFile(filename2, file_size2)) {
			std::remove(temp_file_paths[run_idx].c_str());
			temp_file_paths[run_idx].clear();
			throw IOException("read_ena_sequences: Aspera stream ended unexpectedly (expected R2 for run '%s')",
			                  run.run_accession);
		}

		std::string gz_hint2 = filename2;
		if (gz_hint2.empty() && run.aspera_paths.size() >= 2) {
			gz_hint2 = run.aspera_paths[1].remote_path;
		}
		bool is_gz2 = IsGzipped(gz_hint2);

		auto *s1 = CreateDuckDBSeqStream(fs, temp_file_paths[run_idx], is_gz);
		miint::AsperaSeqStream *s2 = nullptr;
		try {
			s2 = miint::CreateAsperaSeqStream(proc, is_gz2);
			readers[run_idx] = std::make_unique<miint::SequenceReader>(s1, s2, true);
		} catch (...) {
			delete s2;
			delete s1;
			throw;
		}
	}
}
#endif // MIINT_ASPERA_SUPPORTED

// ---- Bind ----

// Helper: get column value from a TSV row, or empty string if not present
static std::string GetCol(const std::vector<std::string> &row, int col) {
	return (col >= 0 && col < static_cast<int>(row.size())) ? row[col] : "";
}

// Resolve accessions to run info via ENA Portal API
static std::vector<ReadENASequencesTableFunction::RunInfo>
ResolveRuns(miint::ENAClient &client, const std::vector<std::string> &accessions, const std::string &prefer_format) {
	std::vector<ReadENASequencesTableFunction::RunInfo> runs;

	for (const auto &acc : accessions) {
		auto tsv = client.Search(acc, "read_run",
		                         "run_accession,sample_accession,experiment_accession,"
		                         "fastq_ftp,fastq_aspera,fastq_bytes,library_layout,"
		                         "submitted_ftp,submitted_aspera,submitted_bytes,submitted_format");
		auto parsed = miint::ENAParser::ParseTSV(tsv);

		// Find column indices by name
		int run_col = -1, sample_col = -1, exp_col = -1, ftp_col = -1, aspera_col = -1, bytes_col = -1, layout_col = -1;
		int sub_ftp_col = -1, sub_aspera_col = -1, sub_bytes_col = -1, sub_format_col = -1;
		for (size_t i = 0; i < parsed.column_names.size(); i++) {
			const auto &name = parsed.column_names[i];
			if (name == "run_accession")
				run_col = static_cast<int>(i);
			else if (name == "sample_accession")
				sample_col = static_cast<int>(i);
			else if (name == "experiment_accession")
				exp_col = static_cast<int>(i);
			else if (name == "fastq_ftp")
				ftp_col = static_cast<int>(i);
			else if (name == "fastq_aspera")
				aspera_col = static_cast<int>(i);
			else if (name == "fastq_bytes")
				bytes_col = static_cast<int>(i);
			else if (name == "library_layout")
				layout_col = static_cast<int>(i);
			else if (name == "submitted_ftp")
				sub_ftp_col = static_cast<int>(i);
			else if (name == "submitted_aspera")
				sub_aspera_col = static_cast<int>(i);
			else if (name == "submitted_bytes")
				sub_bytes_col = static_cast<int>(i);
			else if (name == "submitted_format")
				sub_format_col = static_cast<int>(i);
		}

		for (const auto &row : parsed.rows) {
			std::string run_accession = GetCol(row, run_col);
			std::string sample_accession = GetCol(row, sample_col);
			std::string experiment_accession = GetCol(row, exp_col);
			std::string ftp_field = GetCol(row, ftp_col);
			std::string aspera_field = GetCol(row, aspera_col);
			std::string bytes_field = GetCol(row, bytes_col);
			std::string layout = GetCol(row, layout_col);

			auto fastq_urls = miint::ENAParser::FTPtoHTTPS(ftp_field);
			bool has_fastq = !fastq_urls.empty();

			// Check for submitted SFF files
			std::string sub_ftp = GetCol(row, sub_ftp_col);
			std::string sub_aspera = GetCol(row, sub_aspera_col);
			std::string sub_bytes = GetCol(row, sub_bytes_col);
			std::string sub_format = GetCol(row, sub_format_col);
			auto sff_filter =
			    miint::ENAParser::FilterSubmittedByFormat(sub_ftp, sub_aspera, sub_format, sub_bytes, "SFF");
			bool has_sff = !sff_filter.urls.empty();

			// Format selection:
			// - 'sff':   use SFF if available, otherwise fall back to FASTQ silently
			// - 'auto':  use FASTQ if available, otherwise fall back to SFF
			// - 'fastq': FASTQ only, skip runs that lack FASTQ
			bool use_sff = false;
			if (prefer_format == "sff") {
				use_sff = has_sff;
			} else if (prefer_format == "auto") {
				use_sff = !has_fastq && has_sff;
			}

			if (use_sff) {
				// Flatten: one RunInfo per SFF file, all sharing the same accession metadata
				for (size_t si = 0; si < sff_filter.urls.size(); si++) {
					ReadENASequencesTableFunction::RunInfo info;
					info.run_accession = run_accession;
					info.sample_accession = sample_accession;
					info.experiment_accession = experiment_accession;
					info.format = ENASequenceFormat::SFF;
					info.sff_url = sff_filter.urls[si];
					info.is_paired = false; // SFF is always single-end
					info.total_bytes = sff_filter.total_bytes / sff_filter.urls.size();
					runs.push_back(std::move(info));
				}
			} else if (has_fastq) {
				// Standard FASTQ path (existing behavior)
				ReadENASequencesTableFunction::RunInfo info;
				info.run_accession = run_accession;
				info.sample_accession = sample_accession;
				info.experiment_accession = experiment_accession;
				info.fastq_urls = std::move(fastq_urls);
				info.aspera_paths = miint::AsperaUtils::ParseAsperaPaths(aspera_field);
				info.is_paired = (layout == "PAIRED");

				// Parse fastq_bytes
				std::vector<uint64_t> file_bytes;
				if (!bytes_field.empty()) {
					std::string::size_type bstart = 0;
					while (true) {
						auto bsemi = bytes_field.find(';', bstart);
						std::string bs = (bsemi == std::string::npos) ? bytes_field.substr(bstart)
						                                              : bytes_field.substr(bstart, bsemi - bstart);
						try {
							file_bytes.push_back(std::stoull(bs));
						} catch (...) {
							file_bytes.push_back(0);
						}
						if (bsemi == std::string::npos) {
							break;
						}
						bstart = bsemi + 1;
					}
				}

				// ENA 3-file paired-end filtering (existing logic)
				if (info.is_paired && info.fastq_urls.size() > 2) {
					std::vector<std::string> filtered;
					std::vector<uint64_t> filtered_bytes;
					for (size_t fi = 0; fi < info.fastq_urls.size(); fi++) {
						if (info.fastq_urls[fi].find("_1.fast") != std::string::npos ||
						    info.fastq_urls[fi].find("_2.fast") != std::string::npos) {
							filtered.push_back(info.fastq_urls[fi]);
							if (fi < file_bytes.size()) {
								filtered_bytes.push_back(file_bytes[fi]);
							}
						}
					}
					if (filtered.size() == 2) {
						info.fastq_urls = std::move(filtered);
						file_bytes = std::move(filtered_bytes);
					}
				}
				if (info.is_paired && info.aspera_paths.size() > 2) {
					std::vector<miint::AsperaPath> filtered;
					for (const auto &ap : info.aspera_paths) {
						if (ap.remote_path.find("_1.fast") != std::string::npos ||
						    ap.remote_path.find("_2.fast") != std::string::npos) {
							filtered.push_back(ap);
						}
					}
					if (filtered.size() == 2) {
						info.aspera_paths = std::move(filtered);
					}
				}

				info.total_bytes = 0;
				for (auto b : file_bytes) {
					info.total_bytes += b;
				}

				runs.push_back(std::move(info));
			}
			// else: no FASTQ and no SFF (or prefer_format='fastq' with no FASTQ) — skip run
		}
	}

	return runs;
}

unique_ptr<FunctionData> ReadENASequencesTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
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
		throw InvalidInputException("read_ena_sequences: first argument must be VARCHAR or VARCHAR[]");
	}

	if (accessions.empty()) {
		throw InvalidInputException("read_ena_sequences: at least one accession must be provided");
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
			throw InvalidInputException("read_ena_sequences: accession cannot be empty");
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
			    "read_ena_sequences: download_method must be 'auto', 'aspera', or 'http' (got '%s')", download_method);
		}
	}

	std::string prefer_format = "auto";
	auto pf_param = input.named_parameters.find("prefer_format");
	if (pf_param != input.named_parameters.end() && !pf_param->second.IsNull()) {
		prefer_format = pf_param->second.ToString();
		if (prefer_format != "auto" && prefer_format != "fastq" && prefer_format != "sff") {
			throw InvalidInputException(
			    "read_ena_sequences: prefer_format must be 'auto', 'fastq', or 'sff' (got '%s')", prefer_format);
		}
	}

	bool trim = true;
	auto trim_param = input.named_parameters.find("trim");
	if (trim_param != input.named_parameters.end() && !trim_param->second.IsNull()) {
		trim = trim_param->second.GetValue<bool>();
	}

	// Resolve accessions to run info via ENA Portal API
	auto &db = DatabaseInstance::GetDatabase(context);
	miint::ENAClient client(db);
	auto runs = ResolveRuns(client, accessions, prefer_format);

	if (runs.empty()) {
		throw IOException("read_ena_sequences: no sequence data found for the provided accession(s)");
	}

	auto data = make_uniq<Data>(std::move(runs), include_filepath, qual_offset, trim, prefer_format);

	// Determine whether to use Aspera (only relevant for FASTX runs; SFF uses HTTP only)
	data->use_aspera = false;
#if MIINT_ASPERA_SUPPORTED
	// Check if there are any FASTX runs that could benefit from Aspera
	bool has_fastx_runs = false;
	for (const auto &run : data->runs) {
		if (run.format == ENASequenceFormat::FASTX) {
			has_fastx_runs = true;
			break;
		}
	}
	if (has_fastx_runs && (download_method == "aspera" || download_method == "auto")) {
		bool aspera_available = miint::AsperaUtils::IsAvailable();
		if (download_method == "aspera" && !aspera_available) {
			throw IOException("read_ena_sequences: download_method='aspera' but ascp not found in PATH. "
			                  "Install IBM Aspera CLI or use download_method='http'");
		}
		if (aspera_available) {
			// Check that all FASTX runs have aspera paths
			bool all_have_aspera = true;
			for (const auto &run : data->runs) {
				if (run.format == ENASequenceFormat::FASTX && run.aspera_paths.empty()) {
					all_have_aspera = false;
					break;
				}
			}
			if (download_method == "aspera" && !all_have_aspera) {
				throw IOException("read_ena_sequences: download_method='aspera' but not all runs have Aspera paths");
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
		throw IOException("read_ena_sequences: Aspera is not supported on this platform");
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

unique_ptr<GlobalTableFunctionState> ReadENASequencesTableFunction::InitGlobal(ClientContext &context,
                                                                               TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	FileSystem &fs = FileSystem::GetFileSystem(context);
	auto state = make_uniq<GlobalState>(fs, data.runs, data.use_aspera, data.trim);

#if MIINT_ASPERA_SUPPORTED
	if (data.use_aspera) {
		state->aspera_config = data.aspera_config;
	}
#endif

	return state;
}

unique_ptr<LocalTableFunctionState> ReadENASequencesTableFunction::InitLocal(ExecutionContext &context,
                                                                             TableFunctionInitInput &input,
                                                                             GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// ---- Execute ----

void ReadENASequencesTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	miint::SequenceRecordBatch batch;
	size_t current_run_idx;

	// Helper: tear down a failed run's reader and process
	auto cleanup_run = [&](size_t idx) {
		auto &run = global_state.runs[idx];
		if (run.format == ENASequenceFormat::SFF) {
			global_state.sff_readers[idx].reset();
			if (!global_state.sff_temp_paths[idx].empty()) {
				std::remove(global_state.sff_temp_paths[idx].c_str());
				global_state.sff_temp_paths[idx].clear();
			}
		} else {
			global_state.readers[idx].reset();
		}
#if MIINT_ASPERA_SUPPORTED
		if (global_state.use_aspera && global_state.aspera_processes[idx]) {
			global_state.aspera_processes[idx].reset();
		}
		if (!global_state.temp_file_paths[idx].empty()) {
			std::remove(global_state.temp_file_paths[idx].c_str());
			global_state.temp_file_paths[idx].clear();
		}
#endif
	};

	// Helper: open reader for the current run
	auto open_run = [&](size_t idx) {
		auto &run = global_state.runs[idx];
		if (run.format == ENASequenceFormat::SFF) {
			global_state.OpenReaderSFF(idx);
		} else {
#if MIINT_ASPERA_SUPPORTED
			if (global_state.use_aspera) {
				global_state.OpenReaderAspera(idx);
			} else {
				global_state.OpenReader(idx);
			}
#else
			global_state.OpenReader(idx);
#endif
		}
	};

	// Claim-read-release loop (same pattern as read_fastx)
	while (true) {
		if (!local_state.has_run) {
			bool all_claimed = false;
			{
				lock_guard<mutex> read_lock(global_state.lock);
				if (global_state.next_run_idx >= global_state.runs.size()) {
					all_claimed = true;
				} else {
					local_state.current_run_idx = global_state.next_run_idx;
					global_state.next_run_idx++;
					local_state.has_run = true;
				}
			} // Lock released before I/O or warning output

			if (all_claimed) {
				// Print skipped-runs warning exactly once across all threads
				if (!global_state.skipped_warned.exchange(true, std::memory_order_acq_rel)) {
					lock_guard<mutex> skip_guard(global_state.skipped_lock);
					if (!global_state.skipped_runs.empty()) {
						std::string msg =
						    "read_ena_sequences: WARNING: " + std::to_string(global_state.skipped_runs.size()) +
						    " run(s) skipped due to download errors:";
						for (const auto &acc : global_state.skipped_runs) {
							msg += " " + acc;
						}
						Printer::Print(msg);
					}
				}
				output.SetCardinality(0);
				return;
			}

			// Open reader — retry once on failure before skipping
			try {
				open_run(local_state.current_run_idx);
			} catch (const std::exception &e) {
				auto &run = global_state.runs[local_state.current_run_idx];
				Printer::PrintF("read_ena_sequences: warning: run '%s' failed to open (%s), retrying...",
				                run.run_accession, e.what());
				cleanup_run(local_state.current_run_idx);
				global_state.run_sequence_counters[local_state.current_run_idx] = 1;
				try {
					open_run(local_state.current_run_idx);
				} catch (const std::exception &e2) {
					Printer::PrintF("read_ena_sequences: warning: run '%s' failed on retry (%s), skipping",
					                run.run_accession, e2.what());
					cleanup_run(local_state.current_run_idx);
					{
						lock_guard<mutex> skip_guard(global_state.skipped_lock);
						global_state.skipped_runs.push_back(run.run_accession);
					}
					global_state.bytes_completed.fetch_add(run.total_bytes, std::memory_order_relaxed);
					global_state.runs_completed.fetch_add(1, std::memory_order_relaxed);
					local_state.has_run = false;
					continue;
				}
			}
		}

		current_run_idx = local_state.current_run_idx;

		// Read batch — dispatch based on format
		auto read_batch = [&]() -> miint::SequenceRecordBatch {
			auto &run = global_state.runs[current_run_idx];
			if (run.format == ENASequenceFormat::SFF) {
				return global_state.sff_readers[current_run_idx]->read(STANDARD_VECTOR_SIZE);
			} else {
				return global_state.readers[current_run_idx]->read(STANDARD_VECTOR_SIZE);
			}
		};

		try {
			batch = read_batch();
		} catch (const std::exception &e) {
			// Mid-stream read failure. Tear down and retry from scratch.
			auto &run = global_state.runs[current_run_idx];
			Printer::PrintF("read_ena_sequences: warning: run '%s' failed mid-stream (%s), retrying...",
			                run.run_accession, e.what());
			cleanup_run(current_run_idx);
			global_state.run_sequence_counters[current_run_idx] = 1;
			try {
				open_run(current_run_idx);
				batch = read_batch();
			} catch (const std::exception &e2) {
				Printer::PrintF("read_ena_sequences: warning: run '%s' failed on retry (%s), skipping",
				                run.run_accession, e2.what());
				cleanup_run(current_run_idx);
				{
					lock_guard<mutex> skip_guard(global_state.skipped_lock);
					global_state.skipped_runs.push_back(run.run_accession);
				}
				global_state.bytes_completed.fetch_add(run.total_bytes, std::memory_order_relaxed);
				global_state.runs_completed.fetch_add(1, std::memory_order_relaxed);
				local_state.has_run = false;
				continue;
			}
		}

		if (batch.empty()) {
			// Run completed successfully — release resources
			auto &completed_run = global_state.runs[current_run_idx];
			if (completed_run.format == ENASequenceFormat::SFF) {
				global_state.sff_readers[current_run_idx].reset();
				if (!global_state.sff_temp_paths[current_run_idx].empty()) {
					std::remove(global_state.sff_temp_paths[current_run_idx].c_str());
					global_state.sff_temp_paths[current_run_idx].clear();
				}
			} else {
				global_state.readers[current_run_idx].reset();
			}
#if MIINT_ASPERA_SUPPORTED
			if (global_state.use_aspera && global_state.aspera_processes[current_run_idx]) {
				int exit_code = global_state.aspera_processes[current_run_idx]->WaitForExit();
				global_state.aspera_processes[current_run_idx].reset();
				if (exit_code != 0 && exit_code != -1) {
					throw IOException("read_ena_sequences: ascp exited with code %d for run '%s'", exit_code,
					                  completed_run.run_accession);
				}
			}
			if (!global_state.temp_file_paths[current_run_idx].empty()) {
				std::remove(global_state.temp_file_paths[current_run_idx].c_str());
				global_state.temp_file_paths[current_run_idx].clear();
			}
#endif
			global_state.bytes_completed.fetch_add(global_state.runs[current_run_idx].total_bytes,
			                                       std::memory_order_relaxed);
			global_state.runs_completed.fetch_add(1, std::memory_order_relaxed);
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
	// SFF stores raw Phred internally with offset 33; user's qual_offset only applies to FASTQ
	uint8_t effective_qual_offset = (run.format == ENASequenceFormat::SFF) ? 33 : bind_data.qual_offset;
	SetResultVectorListUInt8(output.data[field_idx++], batch.quals1, effective_qual_offset);

	// qual2 (column 6)
	if (batch.is_paired) {
		SetResultVectorListUInt8(output.data[field_idx++], batch.quals2, effective_qual_offset);
	} else {
		SetResultVectorNull(output.data[field_idx++]);
	}

	// run_accession (column 7) - constant for all rows in this chunk
	SetResultVectorFilepath(output.data[field_idx++], run.run_accession);

	// sample_accession (column 8) - constant
	SetResultVectorFilepath(output.data[field_idx++], run.sample_accession);

	// experiment_accession (column 9) - constant
	SetResultVectorFilepath(output.data[field_idx++], run.experiment_accession);

	// filepath (column 10) - optional
	if (bind_data.include_filepath) {
		std::string filepath;
		if (run.format == ENASequenceFormat::SFF) {
			filepath = run.sff_url;
		} else if (bind_data.use_aspera && !run.aspera_paths.empty()) {
			filepath = run.aspera_paths[0].host + ":" + run.aspera_paths[0].remote_path;
			for (size_t u = 1; u < run.aspera_paths.size(); u++) {
				filepath += ";" + run.aspera_paths[u].host + ":" + run.aspera_paths[u].remote_path;
			}
		} else if (!run.fastq_urls.empty()) {
			filepath = run.fastq_urls[0];
			for (size_t u = 1; u < run.fastq_urls.size(); u++) {
				filepath += ";" + run.fastq_urls[u];
			}
		}
		SetResultVectorFilepath(output.data[field_idx++], filepath);
	}

	output.SetCardinality(count);
}

// ---- Progress ----

double ReadENASequencesTableFunction::Progress(ClientContext &context, const FunctionData *bind_data,
                                               const GlobalTableFunctionState *global_state) {
	if (!global_state) {
		return -1.0;
	}
	auto &state = global_state->Cast<GlobalState>();
	if (state.total_runs == 0) {
		return 100.0;
	}

	// Prefer byte-based progress when ENA provided fastq_bytes metadata.
	// Falls back to run-count progress when byte sizes are unavailable.
	// Returning 0.0 when no work is done is safe — DuckDB's Kalman filter
	// defers initialization until progress > 0.
	if (state.total_bytes > 0) {
		uint64_t done = state.bytes_completed.load(std::memory_order_relaxed);
		return std::min(100.0, static_cast<double>(done) / static_cast<double>(state.total_bytes) * 100.0);
	}

	// Fallback: run-count based progress
	idx_t completed = state.runs_completed.load(std::memory_order_relaxed);
	return std::min(100.0, static_cast<double>(completed) / static_cast<double>(state.total_runs) * 100.0);
}

// ---- Registration ----

TableFunction ReadENASequencesTableFunction::GetFunction() {
	auto tf = TableFunction("read_ena_sequences", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	tf.named_parameters["qual_offset"] = LogicalType::BIGINT;
	tf.named_parameters["download_method"] = LogicalType::VARCHAR;
	tf.named_parameters["prefer_format"] = LogicalType::VARCHAR;
	tf.named_parameters["trim"] = LogicalType::BOOLEAN;
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	tf.table_scan_progress = Progress;
	return tf;
}

void ReadENASequencesTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
