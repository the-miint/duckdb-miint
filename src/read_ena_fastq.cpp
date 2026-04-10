#include "read_ena_fastq.hpp"
#include "SequenceRecord.hpp"
#include "duckdb/common/vector_size.hpp"

namespace duckdb {

// ---- Data ----

ReadENAFastqTableFunction::Data::Data(std::vector<RunInfo> runs, bool include_fp, uint8_t offset)
    : runs(std::move(runs)), include_filepath(include_fp), qual_offset(offset),
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

ReadENAFastqTableFunction::GlobalState::GlobalState(FileSystem &fs, const std::vector<RunInfo> &runs)
    : runs(runs), next_run_idx(0), fs(fs) {
	// Pre-allocate slots but don't open any connections yet.
	// Readers are created lazily in OpenReader() to avoid opening hundreds of
	// HTTP connections simultaneously for large projects.
	readers.resize(runs.size());
	run_sequence_counters.resize(runs.size(), 1);
}

void ReadENAFastqTableFunction::GlobalState::OpenReader(size_t run_idx) {
	if (readers[run_idx]) {
		return;
	}
	const auto &run = runs[run_idx];
	if (run.fastq_urls.empty()) {
		throw IOException("read_ena_fastq: no FASTQ URLs available for run '%s'", run.run_accession);
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
	readers[run_idx] = std::make_unique<miint::SequenceReader>(s1, s2);
}

// ---- Bind ----

// Resolve accessions to run info via ENA Portal API
static std::vector<ReadENAFastqTableFunction::RunInfo> ResolveRuns(miint::ENAClient &client,
                                                                   const std::vector<std::string> &accessions) {
	std::vector<ReadENAFastqTableFunction::RunInfo> runs;

	for (const auto &acc : accessions) {
		auto tsv = client.Search(acc, "read_run",
		                         "run_accession,sample_accession,experiment_accession,fastq_ftp,library_layout");
		auto parsed = miint::ENAParser::ParseTSV(tsv);

		// Find column indices by name
		int run_col = -1, sample_col = -1, exp_col = -1, ftp_col = -1, layout_col = -1;
		for (size_t i = 0; i < parsed.column_names.size(); i++) {
			if (parsed.column_names[i] == "run_accession")
				run_col = static_cast<int>(i);
			else if (parsed.column_names[i] == "sample_accession")
				sample_col = static_cast<int>(i);
			else if (parsed.column_names[i] == "experiment_accession")
				exp_col = static_cast<int>(i);
			else if (parsed.column_names[i] == "fastq_ftp")
				ftp_col = static_cast<int>(i);
			else if (parsed.column_names[i] == "library_layout")
				layout_col = static_cast<int>(i);
		}

		for (const auto &row : parsed.rows) {
			ReadENAFastqTableFunction::RunInfo info;
			info.run_accession = (run_col >= 0 && run_col < (int)row.size()) ? row[run_col] : "";
			info.sample_accession = (sample_col >= 0 && sample_col < (int)row.size()) ? row[sample_col] : "";
			info.experiment_accession = (exp_col >= 0 && exp_col < (int)row.size()) ? row[exp_col] : "";

			std::string ftp_field = (ftp_col >= 0 && ftp_col < (int)row.size()) ? row[ftp_col] : "";
			info.fastq_urls = miint::ENAParser::FTPtoHTTPS(ftp_field);

			std::string layout = (layout_col >= 0 && layout_col < (int)row.size()) ? row[layout_col] : "";
			info.is_paired = (layout == "PAIRED");

			if (!info.fastq_urls.empty()) {
				runs.push_back(std::move(info));
			}
		}
	}

	return runs;
}

unique_ptr<FunctionData> ReadENAFastqTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
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
		throw InvalidInputException("read_ena_fastq: first argument must be VARCHAR or VARCHAR[]");
	}

	if (accessions.empty()) {
		throw InvalidInputException("read_ena_fastq: at least one accession must be provided");
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
			throw InvalidInputException("read_ena_fastq: accession cannot be empty");
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

	// Resolve accessions to run info via ENA Portal API
	auto &db = DatabaseInstance::GetDatabase(context);
	miint::ENAClient client(db);
	auto runs = ResolveRuns(client, accessions);

	if (runs.empty()) {
		throw IOException("read_ena_fastq: no FASTQ data found for the provided accession(s)");
	}

	auto data = make_uniq<Data>(std::move(runs), include_filepath, qual_offset);

	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

// ---- InitGlobal / InitLocal ----

unique_ptr<GlobalTableFunctionState> ReadENAFastqTableFunction::InitGlobal(ClientContext &context,
                                                                           TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	FileSystem &fs = FileSystem::GetFileSystem(context);
	return make_uniq<GlobalState>(fs, data.runs);
}

unique_ptr<LocalTableFunctionState> ReadENAFastqTableFunction::InitLocal(ExecutionContext &context,
                                                                         TableFunctionInitInput &input,
                                                                         GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// ---- Execute ----

void ReadENAFastqTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
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
					output.SetCardinality(0);
					return;
				}
				local_state.current_run_idx = global_state.next_run_idx;
				global_state.next_run_idx++;
				local_state.has_run = true;
			} // Lock released before I/O

			// Open reader without holding the lock — safe because this thread
			// has exclusive ownership of the claimed run index
			global_state.OpenReader(local_state.current_run_idx);
		}

		current_run_idx = local_state.current_run_idx;
		batch = global_state.readers[current_run_idx]->read(STANDARD_VECTOR_SIZE);

		if (batch.empty()) {
			// Release the reader to free the HTTP connection for subsequent runs
			global_state.readers[current_run_idx].reset();
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
		std::string filepath = run.fastq_urls[0];
		for (size_t u = 1; u < run.fastq_urls.size(); u++) {
			filepath += ";" + run.fastq_urls[u];
		}
		SetResultVectorFilepath(output.data[field_idx++], filepath);
	}

	output.SetCardinality(count);
}

// ---- Registration ----

TableFunction ReadENAFastqTableFunction::GetFunction() {
	auto tf = TableFunction("read_ena_fastq", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	tf.named_parameters["qual_offset"] = LogicalType::BIGINT;
	return tf;
}

void ReadENAFastqTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
