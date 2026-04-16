#include "align_mafft.hpp"
#include "MafftAligner.hpp"
#include "sequence_table_reader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include <unordered_set>

namespace duckdb {

struct AlignMafftData : public TableFunctionData {
	std::string table_name;
};

struct AlignMafftGlobalState : public GlobalTableFunctionState {
	// Materialized alignment results. All sequences must be in memory
	// simultaneously because MAFFT's PartTree needs the complete set
	// before it can produce any aligned output.
	std::vector<std::string> names;
	std::vector<std::string> sequences;
	std::vector<int32_t> original_lengths;
	int32_t aligned_length = 0;
	idx_t current_row = 0;

	idx_t MaxThreads() const override {
		return 1;
	}
};

static unique_ptr<FunctionData> AlignMafftBind(ClientContext &context, TableFunctionBindInput &input,
                                               vector<LogicalType> &return_types, vector<string> &names) {
	if (input.inputs.empty() || input.inputs[0].IsNull()) {
		throw BinderException("align_mafft requires a sequence table name argument");
	}
	auto table_name = input.inputs[0].ToString();

	auto schema = ValidateSequenceTableSchema(context, table_name);
	if (schema.has_sequence2) {
		throw BinderException("align_mafft does not support paired-end data");
	}

	auto data = make_uniq<AlignMafftData>();
	data->table_name = std::move(table_name);

	names = {"sequence_index", "read_id", "aligned_sequence", "original_length", "aligned_length"};
	return_types = {LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::INTEGER,
	                LogicalType::INTEGER};

	return data;
}

static unique_ptr<GlobalTableFunctionState> AlignMafftInitGlobal(ClientContext &context,
                                                                 TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<AlignMafftData>();
	auto gstate = make_uniq<AlignMafftGlobalState>();

	auto loaded = LoadSingleEndSequences(context, data.table_name, "align_mafft", /*strict=*/true);

	if (loaded.sequences.size() < 2) {
		throw InvalidInputException("align_mafft requires at least 2 sequences, got %d",
		                            static_cast<int>(loaded.sequences.size()));
	}

	// MAFFT requires sequences of at least 6 characters; guard up-front so the
	// user gets a row-specific error instead of an opaque std::invalid_argument
	// out of the aligner.
	for (size_t i = 0; i < loaded.sequences.size(); i++) {
		if (loaded.sequences[i].size() < 6) {
			throw InvalidInputException("align_mafft: sequence for read_id '%s' is too short (%d chars, minimum 6)",
			                            loaded.labels[i], static_cast<int>(loaded.sequences[i].size()));
		}
	}

	// Duplicate read_ids make it impossible to join align_mafft output back to
	// the input table unambiguously (e.g., the deblur workflow pattern).
	std::unordered_set<std::string> seen;
	seen.reserve(loaded.labels.size());
	for (const auto &id : loaded.labels) {
		if (!seen.insert(id).second) {
			throw InvalidInputException("align_mafft: duplicate read_id '%s' in input table", id);
		}
	}

	// MafftAligner::align requires a comments vector sized to match names/sequences;
	// empty strings are ignored when building FASTA headers.
	std::vector<std::string> comments(loaded.labels.size(), "");

	miint::MafftAligner aligner;
	auto result = aligner.align(loaded.labels, comments, loaded.sequences);

	gstate->names = std::move(result.names);
	gstate->sequences = std::move(result.sequences);
	gstate->aligned_length = static_cast<int32_t>(result.aligned_length);
	gstate->original_lengths.reserve(result.original_lengths.size());
	for (auto len : result.original_lengths) {
		gstate->original_lengths.push_back(static_cast<int32_t>(len));
	}

	return gstate;
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

		FlatVector::GetData<int64_t>(output.data[0])[i] = static_cast<int64_t>(row);

		FlatVector::GetData<string_t>(output.data[1])[i] = StringVector::AddString(output.data[1], gstate.names[row]);

		FlatVector::GetData<string_t>(output.data[2])[i] =
		    StringVector::AddString(output.data[2], gstate.sequences[row]);

		FlatVector::GetData<int32_t>(output.data[3])[i] = gstate.original_lengths[row];
		FlatVector::GetData<int32_t>(output.data[4])[i] = gstate.aligned_length;
	}

	gstate.current_row += count;
	output.SetCardinality(count);
}

TableFunction AlignMafftTableFunction::GetFunction() {
	auto tf =
	    TableFunction("align_mafft", {LogicalType::VARCHAR}, AlignMafftExecute, AlignMafftBind, AlignMafftInitGlobal);
	// Match sibling aligners — allow downstream CTAS pipelines to parallelize.
	// Callers that want deterministic order should ORDER BY sequence_index.
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void AlignMafftTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
