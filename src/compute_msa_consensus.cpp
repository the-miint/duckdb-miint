#include "compute_msa_consensus.hpp"

#include "MsaConsensus.hpp"
#include "alignment_functions_internal.hpp"

#include "duckdb/common/types/vector.hpp"
#include "duckdb/function/aggregate_function.hpp"
#include "duckdb/planner/expression/bound_aggregate_expression.hpp"

#include <cctype>
#include <cstdint>
#include <string>
#include <vector>

namespace duckdb {

namespace {

// Per-group accumulator. Every input row is materialised into the holder as
// (aligned_seq, qual). MAFFT places gaps inconsistently inside HP runs so
// HP-correction needs the full per-read aligned sequence; storing it here
// lets us compute both the column-vote consensus and the median per-read HP
// length without re-reading inputs.
//
// DuckDB requires the aggregate STATE struct to be trivially-move-constructible,
// so the actual storage lives on the heap behind a pointer. Mirrors the
// pointer-to-calculator pattern in compute_coverage_depth.cpp.
struct ConsensusHolder {
	std::vector<std::string> aligned_seqs;
	std::vector<std::vector<std::uint8_t>> quals;

	bool Empty() const {
		return aligned_seqs.empty();
	}
};

struct ConsensusState {
	ConsensusHolder *holder;

	ConsensusState() : holder(nullptr) {
	}

	bool Empty() const {
		return holder == nullptr || holder->Empty();
	}
};

// '-' is the only gap character we recognise. MAFFT can emit '.' for trailing
// unaligned ends in some modes — we treat anything not in {A,C,G,T,a,c,g,t}
// as a gap (no-information), matching the 5-state model in MsaConsensus.
char NormaliseBase(char c) {
	const char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
	if (u == 'A' || u == 'C' || u == 'G' || u == 'T' || u == '-') {
		return u;
	}
	return '-';
}

void IngestRow(ConsensusHolder &holder, const std::string &aligned_seq, const std::uint8_t *qual_data,
               std::size_t qual_len) {
	if (!holder.aligned_seqs.empty() && holder.aligned_seqs.front().size() != aligned_seq.size()) {
		throw InvalidInputException(
		    "compute_msa_consensus: inconsistent MSA row widths in group (%llu vs %llu) — every row in a "
		    "GROUP must share the same gapped alignment length",
		    holder.aligned_seqs.front().size(), aligned_seq.size());
	}

	// Validate qual matches the ungapped count.
	std::size_t nongap = 0;
	for (char c : aligned_seq) {
		if (NormaliseBase(c) != '-') {
			++nongap;
		}
	}
	if (qual_len != nongap) {
		throw InvalidInputException(
		    "compute_msa_consensus: qual list length %llu does not match ungapped sequence length %llu", qual_len,
		    nongap);
	}

	// Normalise once on ingest so we don't pay per-column toupper in finalize.
	std::string norm;
	norm.reserve(aligned_seq.size());
	for (char c : aligned_seq) {
		norm.push_back(NormaliseBase(c));
	}
	holder.aligned_seqs.push_back(std::move(norm));

	std::vector<std::uint8_t> q(qual_data, qual_data + qual_len);
	holder.quals.push_back(std::move(q));
}

struct ConsensusOperation {
	template <class STATE>
	static void Initialize(STATE &state) {
		state.holder = new ConsensusHolder();
	}

	template <class STATE>
	static void Destroy(STATE &state, AggregateInputData &) {
		delete state.holder;
		state.holder = nullptr;
	}

	static void Operation(Vector inputs[], AggregateInputData &, idx_t /*input_count*/, Vector &states, idx_t count) {
		auto &aligned_vec = inputs[0];
		auto &qual_vec = inputs[1];

		UnifiedVectorFormat aligned_data, qual_data, state_data;
		aligned_vec.ToUnifiedFormat(count, aligned_data);
		qual_vec.ToUnifiedFormat(count, qual_data);
		states.ToUnifiedFormat(count, state_data);

		auto aligned_ptr = UnifiedVectorFormat::GetData<string_t>(aligned_data);
		auto qual_entries = UnifiedVectorFormat::GetData<list_entry_t>(qual_data);
		auto state_ptr = UnifiedVectorFormat::GetData<ConsensusState *>(state_data);

		auto &qual_child_vec = ListVector::GetEntry(qual_vec);
		auto qual_child_data = FlatVector::GetData<std::uint8_t>(qual_child_vec);

		for (idx_t i = 0; i < count; ++i) {
			const auto si = state_data.sel->get_index(i);
			const auto ai = aligned_data.sel->get_index(i);
			const auto qi = qual_data.sel->get_index(i);

			if (!aligned_data.validity.RowIsValid(ai) || !qual_data.validity.RowIsValid(qi)) {
				continue;
			}

			auto *state = state_ptr[si];
			const std::string aligned(aligned_ptr[ai].GetData(), aligned_ptr[ai].GetSize());
			const auto entry = qual_entries[qi];
			IngestRow(*state->holder, aligned, qual_child_data + entry.offset, entry.length);
		}
	}

	template <class STATE, class OP>
	static void Combine(const STATE &source, STATE &target, AggregateInputData &) {
		if (source.Empty()) {
			return;
		}
		auto &src = *source.holder;
		auto &tgt = *target.holder;
		if (!tgt.aligned_seqs.empty() && tgt.aligned_seqs.front().size() != src.aligned_seqs.front().size()) {
			throw InvalidInputException(
			    "compute_msa_consensus: inconsistent MSA row widths across parallel partial states (%llu vs %llu)",
			    tgt.aligned_seqs.front().size(), src.aligned_seqs.front().size());
		}
		tgt.aligned_seqs.insert(tgt.aligned_seqs.end(), src.aligned_seqs.begin(), src.aligned_seqs.end());
		tgt.quals.insert(tgt.quals.end(), src.quals.begin(), src.quals.end());
	}

	static void Finalize(Vector &state_vector, AggregateInputData &, Vector &result, idx_t count, idx_t offset) {
		UnifiedVectorFormat state_data;
		state_vector.ToUnifiedFormat(count, state_data);
		auto states = UnifiedVectorFormat::GetData<ConsensusState *>(state_data);

		auto &result_validity = FlatVector::Validity(result);
		auto &entries = StructVector::GetEntries(result);
		auto &seq_vec = *entries[0];
		auto &qual_list_vec = *entries[1];
		auto qual_list_entries = FlatVector::GetData<list_entry_t>(qual_list_vec);

		for (idx_t i = 0; i < count; ++i) {
			const auto si = state_data.sel->get_index(i);
			auto &state = *states[si];

			if (state.Empty()) {
				result_validity.SetInvalid(i + offset);
				qual_list_entries[i + offset].offset = ListVector::GetListSize(qual_list_vec);
				qual_list_entries[i + offset].length = 0;
				continue;
			}

			const auto consensus = miint::BuildConsensus(state.holder->aligned_seqs, state.holder->quals);
			const auto &seq = consensus.first;
			const auto &qual = consensus.second;

			FlatVector::GetData<string_t>(seq_vec)[i + offset] =
			    StringVector::AddString(seq_vec, seq.data(), seq.size());

			const idx_t list_offset = ListVector::GetListSize(qual_list_vec);
			ListVector::Reserve(qual_list_vec, list_offset + qual.size());
			auto &qual_child = ListVector::GetEntry(qual_list_vec);
			auto qual_child_data = FlatVector::GetData<std::uint8_t>(qual_child);
			for (std::size_t k = 0; k < qual.size(); ++k) {
				qual_child_data[list_offset + k] = qual[k];
			}
			ListVector::SetListSize(qual_list_vec, list_offset + qual.size());
			qual_list_entries[i + offset].offset = list_offset;
			qual_list_entries[i + offset].length = qual.size();
		}
	}

	static bool IgnoreNull() {
		return true;
	}
};

} // namespace

void ComputeMsaConsensusFunction::Register(ExtensionLoader &loader) {
	const auto qual_type = LogicalType::LIST(LogicalType::UTINYINT);
	const auto return_type = LogicalType::STRUCT({{"seq", LogicalType::VARCHAR}, {"qual", qual_type}});

	auto fun = AggregateFunction(
	    "compute_msa_consensus", {LogicalType::VARCHAR, qual_type}, return_type,
	    AggregateFunction::StateSize<ConsensusState>,
	    AggregateFunction::StateInitialize<ConsensusState, ConsensusOperation>, ConsensusOperation::Operation,
	    AggregateFunction::StateCombine<ConsensusState, ConsensusOperation>, ConsensusOperation::Finalize, nullptr,
	    nullptr, AggregateFunction::StateDestroy<ConsensusState, ConsensusOperation>);

	loader.RegisterFunction(fun);
}

} // namespace duckdb
