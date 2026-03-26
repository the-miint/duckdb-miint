#include "align_mafft.hpp"
#include "MafftAligner.hpp"
#include "SequenceReader.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"

namespace duckdb {

struct AlignMafftData : public TableFunctionData {
	std::string input_path;
};

struct AlignMafftGlobalState : public GlobalTableFunctionState {
	std::vector<std::string> names;
	std::vector<std::string> comments;
	std::vector<std::string> sequences;
	std::vector<int> original_lengths;
	int aligned_length = 0;
	idx_t current_row = 0;

	idx_t MaxThreads() const override {
		return 1;
	}
};

struct AlignMafftLocalState : public LocalTableFunctionState {};

static unique_ptr<FunctionData> AlignMafftBind(ClientContext &context, TableFunctionBindInput &input,
                                               vector<LogicalType> &return_types, vector<string> &names) {
	if (input.inputs.empty() || input.inputs[0].IsNull()) {
		throw InvalidInputException("align_mafft requires a file path argument");
	}
	auto path = input.inputs[0].GetValue<string>();

	auto &fs = FileSystem::GetFileSystem(context);
	if (!fs.FileExists(path)) {
		throw IOException("File not found: " + path);
	}

	auto data = make_uniq<AlignMafftData>();
	data->input_path = path;

	names = {"sequence_index", "read_id", "comment", "aligned_sequence", "original_length", "aligned_length"};
	return_types = {LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::VARCHAR,
	                LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER};

	return std::move(data);
}

static unique_ptr<GlobalTableFunctionState> AlignMafftInitGlobal(ClientContext &context,
                                                                 TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<AlignMafftData>();
	auto gstate = make_uniq<AlignMafftGlobalState>();

	// Read all sequences from input file
	miint::SequenceReader reader(data.input_path);
	std::vector<std::string> seq_names;
	std::vector<std::string> seq_comments;
	std::vector<std::string> seq_data;

	while (true) {
		auto batch = reader.read(STANDARD_VECTOR_SIZE);
		if (batch.size() == 0) {
			break;
		}
		for (idx_t i = 0; i < batch.size(); i++) {
			seq_names.push_back(batch.read_ids[i]);
			seq_comments.push_back(batch.comments[i]);
			seq_data.push_back(batch.sequences1[i]);
		}
	}

	if (seq_data.size() < 2) {
		throw InvalidInputException("align_mafft requires at least 2 sequences, got %d",
		                            static_cast<int>(seq_data.size()));
	}

	miint::MafftAligner aligner;
	auto result = aligner.align(seq_names, seq_comments, seq_data);

	gstate->names = std::move(result.names);
	gstate->comments = std::move(result.comments);
	gstate->sequences = std::move(result.sequences);
	gstate->original_lengths = std::move(result.original_lengths);
	gstate->aligned_length = result.aligned_length;

	return std::move(gstate);
}

static unique_ptr<LocalTableFunctionState> AlignMafftInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                               GlobalTableFunctionState *global_state) {
	return make_uniq<AlignMafftLocalState>();
}

static void AlignMafftExecute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<AlignMafftGlobalState>();

	idx_t total = gstate.names.size();
	if (gstate.current_row >= total) {
		output.SetCardinality(0);
		return;
	}

	idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - gstate.current_row);

	for (idx_t i = 0; i < count; i++) {
		idx_t row = gstate.current_row + i;

		FlatVector::GetData<int32_t>(output.data[0])[i] = static_cast<int32_t>(row);

		FlatVector::GetData<string_t>(output.data[1])[i] = StringVector::AddString(output.data[1], gstate.names[row]);

		if (gstate.comments[row].empty()) {
			FlatVector::Validity(output.data[2]).SetInvalid(i);
		} else {
			FlatVector::GetData<string_t>(output.data[2])[i] =
			    StringVector::AddString(output.data[2], gstate.comments[row]);
		}

		FlatVector::GetData<string_t>(output.data[3])[i] =
		    StringVector::AddString(output.data[3], gstate.sequences[row]);

		FlatVector::GetData<int32_t>(output.data[4])[i] = gstate.original_lengths[row];
		FlatVector::GetData<int32_t>(output.data[5])[i] = gstate.aligned_length;
	}

	gstate.current_row += count;
	output.SetCardinality(count);
}

TableFunction AlignMafftTableFunction::GetFunction() {
	auto tf = TableFunction("align_mafft", {LogicalType::VARCHAR}, AlignMafftExecute, AlignMafftBind,
	                        AlignMafftInitGlobal, AlignMafftInitLocal);
	return tf;
}

void AlignMafftTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
