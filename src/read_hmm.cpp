#include "read_hmm.hpp"

#include "HmmerProfile.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/common/vector/string_vector.hpp"

#include <cctype>
#include <optional>
#include <string_view>

namespace duckdb {

namespace {

enum Column : idx_t {
	NAME,
	ACCESSION,
	DESCRIPTION,
	LENGTH,
	ALPHABET,
	FORMAT,
	MAX_LENGTH,
	SEQUENCE_COUNT,
	EFFECTIVE_SEQUENCE_COUNT,
	CHECKSUM,
	GATHERING_SEQUENCE,
	GATHERING_DOMAIN,
	TRUSTED_SEQUENCE,
	TRUSTED_DOMAIN,
	NOISE_SEQUENCE,
	NOISE_DOMAIN,
	CONSENSUS,
	MSV_MU,
	MSV_LAMBDA,
	VITERBI_MU,
	VITERBI_LAMBDA,
	FORWARD_TAU,
	FORWARD_LAMBDA
};

// The first whitespace-delimited word of `line`, or "" for a blank line.
std::string_view FirstWord(std::string_view line) {
	size_t begin = 0;
	while (begin < line.size() && std::isspace(static_cast<unsigned char>(line[begin]))) {
		begin++;
	}
	size_t end = begin;
	while (end < line.size() && !std::isspace(static_cast<unsigned char>(line[end]))) {
		end++;
	}
	return line.substr(begin, end - begin);
}

void SetString(Vector &vec, idx_t row, const std::string &value) {
	FlatVector::GetDataMutable<string_t>(vec)[row] = StringVector::AddString(vec, value);
}

void SetString(Vector &vec, idx_t row, const std::optional<std::string> &value) {
	if (!value) {
		FlatVector::ValidityMutable(vec).SetInvalid(row);
		return;
	}
	SetString(vec, row, *value);
}

template <class T>
void SetNumber(Vector &vec, idx_t row, const std::optional<T> &value) {
	if (!value) {
		FlatVector::ValidityMutable(vec).SetInvalid(row);
		return;
	}
	FlatVector::GetDataMutable<T>(vec)[row] = *value;
}

void WriteRow(DataChunk &output, idx_t row, const miint::hmmer::ProfileHeader &h) {
	auto &out = output.data;
	SetString(out[NAME], row, h.name);
	SetString(out[ACCESSION], row, h.accession);
	SetString(out[DESCRIPTION], row, h.description);
	FlatVector::GetDataMutable<int64_t>(out[LENGTH])[row] = h.length;
	SetString(out[ALPHABET], row, h.alphabet);
	SetString(out[FORMAT], row, h.format);
	SetNumber(out[MAX_LENGTH], row, h.max_length);
	SetNumber(out[SEQUENCE_COUNT], row, h.sequence_count);
	SetNumber(out[EFFECTIVE_SEQUENCE_COUNT], row, h.effective_sequence_count);
	SetNumber(out[CHECKSUM], row, h.checksum);
	SetNumber(out[GATHERING_SEQUENCE], row, h.gathering_sequence);
	SetNumber(out[GATHERING_DOMAIN], row, h.gathering_domain);
	SetNumber(out[TRUSTED_SEQUENCE], row, h.trusted_sequence);
	SetNumber(out[TRUSTED_DOMAIN], row, h.trusted_domain);
	SetNumber(out[NOISE_SEQUENCE], row, h.noise_sequence);
	SetNumber(out[NOISE_DOMAIN], row, h.noise_domain);
	SetString(out[CONSENSUS], row, h.consensus);
	SetNumber(out[MSV_MU], row, h.msv_mu);
	SetNumber(out[MSV_LAMBDA], row, h.msv_lambda);
	SetNumber(out[VITERBI_MU], row, h.viterbi_mu);
	SetNumber(out[VITERBI_LAMBDA], row, h.viterbi_lambda);
	SetNumber(out[FORWARD_TAU], row, h.forward_tau);
	SetNumber(out[FORWARD_LAMBDA], row, h.forward_lambda);
}

// The smallest models HMMER ships, of one to six match states (testsuite/M1.hmm, and the
// HMMER2 ecori model in testsuite/i15-hmmconvert.pl), take 1.3 to 1.5 KB, and a Pfam model
// tens of KB, so 1000 bytes a model errs high, as EstimateFileScanCardinality wants.
unique_ptr<NodeStatistics> ReadHmmCardinality(ClientContext &context, const FunctionData *bind_data) {
	constexpr idx_t BYTES_PER_MODEL = 1000;
	return EstimateFileScanCardinality(bind_data->Cast<LineScanData>().file_paths, BYTES_PER_MODEL);
}

} // namespace

unique_ptr<FunctionData> ReadHmmTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                    vector<LogicalType> &return_types, vector<Identifier> &names) {
	names = {"hmm_name",         "hmm_accession",      "hmm_description",
	         "hmm_length",       "alphabet",           "format",
	         "max_length",       "sequence_count",     "effective_sequence_count",
	         "checksum",         "gathering_sequence", "gathering_domain",
	         "trusted_sequence", "trusted_domain",     "noise_sequence",
	         "noise_domain",     "consensus",          "msv_mu",
	         "msv_lambda",       "viterbi_mu",         "viterbi_lambda",
	         "forward_tau",      "forward_lambda"};
	return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR,  LogicalType::VARCHAR, LogicalType::BIGINT,
	                LogicalType::VARCHAR, LogicalType::VARCHAR,  LogicalType::BIGINT,  LogicalType::BIGINT,
	                LogicalType::FLOAT,   LogicalType::UINTEGER, LogicalType::FLOAT,   LogicalType::FLOAT,
	                LogicalType::FLOAT,   LogicalType::FLOAT,    LogicalType::FLOAT,   LogicalType::FLOAT,
	                LogicalType::VARCHAR, LogicalType::FLOAT,    LogicalType::FLOAT,   LogicalType::FLOAT,
	                LogicalType::FLOAT,   LogicalType::FLOAT,    LogicalType::FLOAT};
	return BindLineScan(context, input, names, return_types);
}

unique_ptr<LocalTableFunctionState> ReadHmmTableFunction::InitLocal(ExecutionContext &context,
                                                                    TableFunctionInitInput &input,
                                                                    GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// Each model is a record: the lines up to and including the next whose first word is "//",
// which HMMER reads as the model's end. Blank and comment lines between models belong to the
// next record, where HMMER skips them, so the line it reports counts from the previous "//".
// What follows the last "//" is a model only if it holds a line that is neither blank nor a
// comment: a model cut off before its "//" then fails in HMMER, as it should, rather than
// being dropped.
//
// As ScanLines does, each thread reads its files one at a time, start to end, and a chunk
// holds one file's rows only, so the filepath column is constant.
void ReadHmmTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &data = data_p.bind_data->Cast<LineScanData>();
	auto &global = data_p.global_state->Cast<LineScanGlobalState>();
	auto &local = data_p.local_state->Cast<LocalState>();

	idx_t row = 0;
	idx_t chunk_bytes = 0;
	std::string chunk_path = local.file ? local.file->path : std::string();

	auto read_model = [&]() {
		local.models_read++;
		miint::hmmer::ProfileHeader header;
		try {
			header = miint::hmmer::ParseProfile(local.record);
		} catch (const miint::hmmer::ProfileParseError &e) {
			const uint64_t line = local.record_start + static_cast<uint64_t>(e.Line()) - 1;
			throw InvalidInputException(data.function_name + ": " + local.file->path + ": model " +
			                            std::to_string(local.models_read) + " (line " + std::to_string(line) +
			                            "): " + e.what());
		}
		WriteRow(output, row, header);
		row++;
		chunk_bytes += local.record.size();
	};
	auto start_record = [&](uint64_t first_line) {
		local.record.clear();
		local.record_start = first_line;
		local.record_has_text = false;
	};

	while (row < STANDARD_VECTOR_SIZE && chunk_bytes < LINE_SCAN_MAX_CHUNK_BYTES) {
		if (!local.file) {
			if (row > 0) {
				break;
			}
			size_t idx;
			{
				std::lock_guard<std::mutex> guard(global.lock);
				if (global.next_file_idx >= data.file_paths.size()) {
					break;
				}
				idx = global.next_file_idx++;
			}
			local.file =
			    make_uniq<LineFile>(FileSystem::GetFileSystem(context), data.file_paths[idx], data.function_name);
			chunk_path = local.file->path;
			local.models_read = 0;
			start_record(1);
		}

		std::string_view line;
		if (!local.file->lines.Next(line)) {
			if (local.record_has_text) {
				read_model();
			}
			local.file.reset();
			continue;
		}
		local.record.append(line);
		local.record.push_back('\n');
		const auto word = FirstWord(line);
		if (!word.empty() && word[0] != '#') {
			local.record_has_text = true;
		}
		if (word == "//") {
			read_model();
			start_record(local.file->lines.LineNumber() + 1);
		}
	}
	if (data.include_filepath && row > 0) {
		SetResultVectorFilepath(output.data[output.ColumnCount() - 1], chunk_path);
	}
	output.SetChildCardinality(row);
}

void ReadHmmTableFunction::Register(ExtensionLoader &loader) {
	auto function = MakeLineScanFunction("read_hmm", Execute, Bind, InitLocal);
	function.cardinality = ReadHmmCardinality;
	loader.RegisterFunction(function);
}

} // namespace duckdb
