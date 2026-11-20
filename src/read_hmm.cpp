#include "read_hmm.hpp"

#include "HmmerProfile.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/common/vector/string_vector.hpp"

#include <optional>
#include <string_view>

namespace duckdb {

namespace {

namespace hmmer = miint::hmmer;

enum Column : idx_t {
	NAME,
	ACCESSION,
	DESCRIPTION,
	LENGTH,
	ALPHABET,
	FORMAT,
	MAX_INSTANCE_LENGTH,
	SEQUENCE_COUNT,
	EFFECTIVE_SEQUENCE_COUNT,
	CHECKSUM,
	GATHERING_THRESHOLD_SEQUENCE,
	GATHERING_THRESHOLD_DOMAIN,
	TRUSTED_CUTOFF_SEQUENCE,
	TRUSTED_CUTOFF_DOMAIN,
	NOISE_CUTOFF_SEQUENCE,
	NOISE_CUTOFF_DOMAIN,
	CONSENSUS,
	MSV_MU,
	MSV_LAMBDA,
	VITERBI_MU,
	VITERBI_LAMBDA,
	FORWARD_TAU,
	FORWARD_LAMBDA,
	COLUMN_COUNT
};

struct ColumnSpec {
	const char *name;
	LogicalTypeId type;
};

// In Column's order.
constexpr ColumnSpec COLUMNS[] = {
    {"hmm_name", LogicalTypeId::VARCHAR},
    {"hmm_accession", LogicalTypeId::VARCHAR},
    {"hmm_description", LogicalTypeId::VARCHAR},
    {"hmm_length", LogicalTypeId::BIGINT},
    {"alphabet", LogicalTypeId::VARCHAR},
    {"format", LogicalTypeId::VARCHAR},
    {"max_instance_length", LogicalTypeId::BIGINT},
    {"sequence_count", LogicalTypeId::BIGINT},
    {"effective_sequence_count", LogicalTypeId::FLOAT},
    {"checksum", LogicalTypeId::UINTEGER},
    {"gathering_threshold_sequence", LogicalTypeId::FLOAT},
    {"gathering_threshold_domain", LogicalTypeId::FLOAT},
    {"trusted_cutoff_sequence", LogicalTypeId::FLOAT},
    {"trusted_cutoff_domain", LogicalTypeId::FLOAT},
    {"noise_cutoff_sequence", LogicalTypeId::FLOAT},
    {"noise_cutoff_domain", LogicalTypeId::FLOAT},
    {"consensus", LogicalTypeId::VARCHAR},
    {"msv_mu", LogicalTypeId::FLOAT},
    {"msv_lambda", LogicalTypeId::FLOAT},
    {"viterbi_mu", LogicalTypeId::FLOAT},
    {"viterbi_lambda", LogicalTypeId::FLOAT},
    {"forward_tau", LogicalTypeId::FLOAT},
    {"forward_lambda", LogicalTypeId::FLOAT},
};
static_assert(sizeof(COLUMNS) / sizeof(COLUMNS[0]) == COLUMN_COUNT, "one ColumnSpec per Column");

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

void WriteRow(DataChunk &output, idx_t row, const hmmer::ProfileHeader &h) {
	auto &out = output.data;
	SetString(out[NAME], row, h.name);
	SetString(out[ACCESSION], row, h.accession);
	SetString(out[DESCRIPTION], row, h.description);
	FlatVector::GetDataMutable<int64_t>(out[LENGTH])[row] = h.length;
	SetString(out[ALPHABET], row, h.alphabet);
	SetString(out[FORMAT], row, h.format);
	SetNumber(out[MAX_INSTANCE_LENGTH], row, h.max_instance_length);
	SetNumber(out[SEQUENCE_COUNT], row, h.sequence_count);
	SetNumber(out[EFFECTIVE_SEQUENCE_COUNT], row, h.effective_sequence_count);
	SetNumber(out[CHECKSUM], row, h.checksum);
	SetNumber(out[GATHERING_THRESHOLD_SEQUENCE], row, h.gathering_threshold_sequence);
	SetNumber(out[GATHERING_THRESHOLD_DOMAIN], row, h.gathering_threshold_domain);
	SetNumber(out[TRUSTED_CUTOFF_SEQUENCE], row, h.trusted_cutoff_sequence);
	SetNumber(out[TRUSTED_CUTOFF_DOMAIN], row, h.trusted_cutoff_domain);
	SetNumber(out[NOISE_CUTOFF_SEQUENCE], row, h.noise_cutoff_sequence);
	SetNumber(out[NOISE_CUTOFF_DOMAIN], row, h.noise_cutoff_domain);
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
	for (const auto &column : COLUMNS) {
		names.emplace_back(column.name);
		return_types.emplace_back(column.type);
	}
	return BindLineScan(context, input, names, return_types);
}

unique_ptr<LocalTableFunctionState> ReadHmmTableFunction::InitLocal(ExecutionContext &context,
                                                                    TableFunctionInitInput &input,
                                                                    GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// Each model is a record: its lines from the first that HMMER doesn't skip, which holds its
// format tag, through the next whose first word is "//", which HMMER reads as its end. The blank
// and comment lines before a model are left out, so HMMER's line numbers count from the record's
// first line. Text after the last "//" that isn't blank or comment is a model too: one cut off
// before its "//" then fails in HMMER, as it should, rather than being dropped.
//
// A record goes to HMMER early when it can't be a model: when its first line holds no format tag
// (HMMER rejects it from that line alone), or when it has outgrown what HMMER can read. Reading on
// to a "//" would otherwise hold all of a file that isn't HMMER's, a FASTA file say, in memory.
//
// As ScanLines does, each thread reads its files one at a time, start to end, and a chunk holds
// one file's rows only, so the filepath column is constant.
void ReadHmmTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &data = data_p.bind_data->Cast<LineScanData>();
	auto &global = data_p.global_state->Cast<LineScanGlobalState>();
	auto &local = data_p.local_state->Cast<LocalState>();

	idx_t row = 0;
	idx_t chunk_bytes = 0;
	std::string chunk_path = local.file ? local.file->path : std::string();

	auto read_model = [&]() {
		// HMMER calibrates a HMMER2 model as it reads it, which takes milliseconds, and one chunk
		// can hold thousands: a cancel or max_execution_time is noticed here, between models.
		// The blocking check reads the clock each time; DuckDB's own reads it once in 256.
		context.BlockingInterruptCheck();
		local.models_read++;
		hmmer::ProfileHeader header;
		try {
			header = hmmer::ParseProfile(local.record);
		} catch (const hmmer::ProfileParseError &e) {
			const uint64_t line = local.record_start + static_cast<uint64_t>(e.Line()) - 1;
			throw InvalidInputException(data.function_name + ": " + local.file->path + ": model " +
			                            std::to_string(local.models_read) + " (line " + std::to_string(line) +
			                            "): " + e.what());
		}
		WriteRow(output, row, header);
		row++;
		chunk_bytes += local.record.size();
		local.record.clear();
	};

	while (row < STANDARD_VECTOR_SIZE && chunk_bytes < LINE_SCAN_MAX_CHUNK_BYTES) {
		if (!local.file) {
			if (row > 0 || !ClaimNextFile(context, data, global, local)) {
				break;
			}
			chunk_path = local.file->path;
			local.record.clear();
			local.models_read = 0;
		}

		std::string_view line;
		if (!local.file->lines.Next(line)) {
			if (!local.record.empty()) {
				read_model();
			}
			local.file.reset();
			continue;
		}
		const auto word = hmmer::FirstWord(line);
		const bool first_line = local.record.empty();
		if (first_line) {
			if (hmmer::IsBlankOrComment(word)) {
				continue;
			}
			local.record_start = local.file->lines.LineNumber();
		}
		local.record.append(line);
		local.record.push_back('\n');
		if ((first_line && !hmmer::IsFormatTag(word)) || hmmer::IsModelEnd(word) ||
		    local.record.size() > hmmer::MAX_RECORD_BYTES) {
			read_model();
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
