#include "compute_pileup.hpp"

#include "PileupWalker.hpp"
#include "alignment_functions_internal.hpp"
#include "catalog_utils.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/main/materialized_query_result.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <algorithm>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace duckdb {

static constexpr const char *FN_NAME = "compute_pileup";

// ---------------------------------------------------------------------------
// Table-function state
// ---------------------------------------------------------------------------
struct PileupBindData : public TableFunctionData {
	std::string alignments_table;
	std::string reference_table;
};

struct PileupGlobalState : public GlobalTableFunctionState {
	std::unordered_map<std::string, std::string> ref;

	// Streaming alignment reader — conn must outlive alignment_stream.
	unique_ptr<Connection> conn;
	unique_ptr<QueryResult> alignment_stream;
	bool stream_exhausted = false;

	// Result buffer: pileup rows from the current alignment chunk.
	std::vector<miint::PileupRow> result_buffer;
	idx_t buffer_offset = 0;

	idx_t MaxThreads() const override {
		return 1;
	}
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
static void ValidateTableSchema(ClientContext &context, const std::string &table_name, const std::string &probe,
                                const char *role) {
	auto conn = MakeReadOnlyHelperConnection(context);
	std::string query = "SELECT " + probe + " FROM " + KeywordHelper::WriteOptionallyQuoted(table_name) + " LIMIT 0";
	auto result = conn.Query(query);
	if (result->HasError()) {
		throw BinderException("%s: %s table '%s' missing required column(s) — expected (%s) (%s)", FN_NAME, role,
		                      table_name, probe, result->GetError());
	}
}

// Load reference table → ref_id → sequence map.
//
// NOTE: deviation from the plan, which suggested `string_view`. We use owning
// `std::string` because the source `string_t` from MaterializedQueryResult is
// only valid for the lifetime of `result`, which is local to this function.
// A string_view map would dangle as soon as LoadReference returned. For
// human-scale references (chr1 = 250 MB) this doubles peak memory; for the
// Karst UMI use case (per-cluster centroids of a few kbp) it is negligible.
static std::unordered_map<std::string, std::string> LoadReference(ClientContext &context,
                                                                  const std::string &table_name) {
	std::unordered_map<std::string, std::string> ref;
	auto conn = MakeReadOnlyHelperConnection(context);
	std::string query = "SELECT ref_id, sequence FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);
	auto result = conn.Query(query);
	if (result->HasError()) {
		throw InvalidInputException("%s: failed to read reference table '%s': %s", FN_NAME, table_name,
		                            result->GetError());
	}
	auto &materialized = result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		UnifiedVectorFormat id_data, seq_data;
		chunk->data[0].ToUnifiedFormat(chunk->size(), id_data);
		chunk->data[1].ToUnifiedFormat(chunk->size(), seq_data);
		auto id_ptr = UnifiedVectorFormat::GetData<string_t>(id_data);
		auto seq_ptr = UnifiedVectorFormat::GetData<string_t>(seq_data);
		for (idx_t i = 0; i < chunk->size(); ++i) {
			auto idi = id_data.sel->get_index(i);
			auto seqi = seq_data.sel->get_index(i);
			if (!id_data.validity.RowIsValid(idi) || !seq_data.validity.RowIsValid(seqi)) {
				throw InvalidInputException("%s: NULL ref_id or sequence in reference table '%s'", FN_NAME, table_name);
			}
			std::string rid(id_ptr[idi].GetData(), id_ptr[idi].GetSize());
			std::string seq(seq_ptr[seqi].GetData(), seq_ptr[seqi].GetSize());
			ref.emplace(std::move(rid), std::move(seq));
		}
	}
	return ref;
}

// Process one chunk of alignments, appending pileup rows to the buffer.
static void ProcessAlignmentChunk(DataChunk &chunk, const std::unordered_map<std::string, std::string> &ref,
                                  std::vector<miint::PileupRow> &rows) {
	UnifiedVectorFormat read_id_data, ref_data, pos_data, cigar_data, seq_data, qual_data;
	chunk.data[0].ToUnifiedFormat(chunk.size(), read_id_data);
	chunk.data[1].ToUnifiedFormat(chunk.size(), ref_data);
	chunk.data[2].ToUnifiedFormat(chunk.size(), pos_data);
	chunk.data[3].ToUnifiedFormat(chunk.size(), cigar_data);
	chunk.data[4].ToUnifiedFormat(chunk.size(), seq_data);
	chunk.data[5].ToUnifiedFormat(chunk.size(), qual_data);

	auto read_id_ptr = UnifiedVectorFormat::GetData<string_t>(read_id_data);
	auto ref_ptr = UnifiedVectorFormat::GetData<string_t>(ref_data);
	auto pos_ptr = UnifiedVectorFormat::GetData<int64_t>(pos_data);
	auto cigar_ptr = UnifiedVectorFormat::GetData<string_t>(cigar_data);
	auto seq_ptr = UnifiedVectorFormat::GetData<string_t>(seq_data);
	auto qual_entries = UnifiedVectorFormat::GetData<list_entry_t>(qual_data);
	auto &qual_child_vec = ListVector::GetEntry(chunk.data[5]);
	auto qual_child_data = FlatVector::GetData<uint8_t>(qual_child_vec);

	for (idx_t i = 0; i < chunk.size(); ++i) {
		auto rdi = read_id_data.sel->get_index(i);
		auto rfi = ref_data.sel->get_index(i);
		auto pi = pos_data.sel->get_index(i);
		auto ci = cigar_data.sel->get_index(i);
		auto si = seq_data.sel->get_index(i);
		auto qi = qual_data.sel->get_index(i);

		if (!read_id_data.validity.RowIsValid(rdi) || !ref_data.validity.RowIsValid(rfi) ||
		    !pos_data.validity.RowIsValid(pi) || !cigar_data.validity.RowIsValid(ci) ||
		    !seq_data.validity.RowIsValid(si)) {
			throw InvalidInputException("%s: NULL in required alignment column", FN_NAME);
		}

		std::string read_id(read_id_ptr[rdi].GetData(), read_id_ptr[rdi].GetSize());
		std::string ref_name(ref_ptr[rfi].GetData(), ref_ptr[rfi].GetSize());
		int64_t pos = pos_ptr[pi];
		std::string cigar(cigar_ptr[ci].GetData(), cigar_ptr[ci].GetSize());
		std::string seq(seq_ptr[si].GetData(), seq_ptr[si].GetSize());

		auto ref_it = ref.find(ref_name);
		if (ref_it == ref.end()) {
			throw InvalidInputException("%s: alignment for read '%s' references '%s' but it is not present in the "
			                            "reference table",
			                            FN_NAME, read_id, ref_name);
		}

		bool qual_is_null = !qual_data.validity.RowIsValid(qi);
		const uint8_t *qual_bytes = nullptr;
		idx_t qual_len = 0;
		if (!qual_is_null) {
			const auto &entry = qual_entries[qi];
			qual_bytes = qual_child_data + entry.offset;
			qual_len = entry.length;
		}

		try {
			miint::PileupWalker::Walk(cigar, ref_name, ref_it->second, pos, read_id, seq, qual_bytes, qual_len,
			                          qual_is_null, rows);
		} catch (const miint::InvalidInputException &e) {
			throw InvalidInputException(e.what());
		}
	}
}

// ---------------------------------------------------------------------------
// Bind
// ---------------------------------------------------------------------------
static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
                                     vector<LogicalType> &return_types, vector<std::string> &names) {
	auto data = make_uniq<PileupBindData>();
	data->alignments_table = input.inputs[0].GetValue<std::string>();
	data->reference_table = input.inputs[1].GetValue<std::string>();

	ValidateTableSchema(context, data->alignments_table, "read_id, reference, position, cigar, sequence, qual",
	                    "alignments");
	ValidateTableSchema(context, data->reference_table, "ref_id, sequence", "reference");

	names = {"ref_id", "ref_pos", "read_id", "ref_base", "query_base", "query_qual", "insert_pos"};
	return_types = {LogicalType::VARCHAR, LogicalType::BIGINT,   LogicalType::VARCHAR, LogicalType::VARCHAR,
	                LogicalType::VARCHAR, LogicalType::UTINYINT, LogicalType::INTEGER};
	return data;
}

// ---------------------------------------------------------------------------
// InitGlobal — load reference (small) and open streaming alignment reader
// ---------------------------------------------------------------------------
static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<PileupBindData>();
	auto gstate = make_uniq<PileupGlobalState>();

	gstate->ref = LoadReference(context, data.reference_table);

	auto &db = DatabaseInstance::GetDatabase(context);
	gstate->conn = make_uniq<Connection>(db);
	InheritTempObjects(context, *gstate->conn);
	std::string query = "SELECT read_id, reference, position, cigar, sequence, qual FROM " +
	                    KeywordHelper::WriteOptionallyQuoted(data.alignments_table);
	gstate->alignment_stream = gstate->conn->SendQuery(query);
	if (gstate->alignment_stream->HasError()) {
		throw InvalidInputException("%s: failed to read alignments table '%s': %s", FN_NAME, data.alignments_table,
		                            gstate->alignment_stream->GetError());
	}

	return gstate;
}

// ---------------------------------------------------------------------------
// Execute — drain result buffer, refill from streaming alignment reader
// ---------------------------------------------------------------------------
static void EmitPileupRows(DataChunk &output, const std::vector<miint::PileupRow> &rows, idx_t offset, idx_t count) {
	auto &ref_id_vec = output.data[0];
	auto ref_pos_data = FlatVector::GetData<int64_t>(output.data[1]);
	auto &read_id_vec = output.data[2];
	auto &ref_base_vec = output.data[3];
	auto &query_base_vec = output.data[4];
	auto query_qual_data = FlatVector::GetData<uint8_t>(output.data[5]);
	auto insert_pos_data = FlatVector::GetData<int32_t>(output.data[6]);
	auto &ref_base_validity = FlatVector::Validity(output.data[3]);
	auto &query_base_validity = FlatVector::Validity(output.data[4]);
	auto &query_qual_validity = FlatVector::Validity(output.data[5]);

	ref_base_validity.SetAllValid(count);
	query_base_validity.SetAllValid(count);
	query_qual_validity.SetAllValid(count);

	for (idx_t i = 0; i < count; ++i) {
		const auto &r = rows[offset + i];
		FlatVector::GetData<string_t>(ref_id_vec)[i] = StringVector::AddString(ref_id_vec, r.ref_id);
		ref_pos_data[i] = r.ref_pos;
		FlatVector::GetData<string_t>(read_id_vec)[i] = StringVector::AddString(read_id_vec, r.read_id);
		insert_pos_data[i] = r.insert_pos;
		if (r.ref_base_is_null) {
			ref_base_validity.SetInvalid(i);
		} else {
			FlatVector::GetData<string_t>(ref_base_vec)[i] = StringVector::AddString(ref_base_vec, &r.ref_base, 1);
		}
		if (r.query_is_null) {
			query_base_validity.SetInvalid(i);
		} else {
			FlatVector::GetData<string_t>(query_base_vec)[i] =
			    StringVector::AddString(query_base_vec, &r.query_base, 1);
		}
		if (r.qual_is_null) {
			query_qual_validity.SetInvalid(i);
		} else {
			query_qual_data[i] = r.query_qual;
		}
	}

	output.SetCardinality(count);
}

static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	(void)context;
	auto &gstate = data_p.global_state->Cast<PileupGlobalState>();

	while (true) {
		idx_t available = gstate.result_buffer.size() - gstate.buffer_offset;
		if (available > 0) {
			idx_t count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			EmitPileupRows(output, gstate.result_buffer, gstate.buffer_offset, count);
			gstate.buffer_offset += count;
			return;
		}

		gstate.result_buffer.clear();
		gstate.buffer_offset = 0;

		if (gstate.stream_exhausted) {
			output.SetCardinality(0);
			return;
		}

		auto chunk = gstate.alignment_stream->Fetch();
		if (!chunk) {
			if (gstate.alignment_stream->HasError()) {
				throw InvalidInputException("%s: alignment stream error: %s", FN_NAME,
				                            gstate.alignment_stream->GetError());
			}
			gstate.stream_exhausted = true;
			output.SetCardinality(0);
			return;
		}
		if (chunk->size() == 0) {
			gstate.stream_exhausted = true;
			output.SetCardinality(0);
			return;
		}

		ProcessAlignmentChunk(*chunk, gstate.ref, gstate.result_buffer);
	}
}

// ---------------------------------------------------------------------------
// Register
// ---------------------------------------------------------------------------
void ComputePileupTableFunction::Register(ExtensionLoader &loader) {
	TableFunction tf("compute_pileup", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal);
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(tf);
}

} // namespace duckdb
