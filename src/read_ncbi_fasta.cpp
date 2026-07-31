#include "read_ncbi_fasta.hpp"
#include "catalog_utils.hpp"
#include "miint_log.hpp"
#include "duckdb/common/vector_size.hpp"
#include <sstream>

namespace duckdb {

namespace {

// Build the fetch plan from the user's accession list. Walks input in order so a
// SQL ORDER BY sequence_index gives back the same ordering the user supplied.
// Consecutive SEQUENCE accessions accumulate into a SEQUENCE_BATCH unit up to
// batch_size; an ASSEMBLY accession flushes the in-flight batch and emits its
// own singleton unit. UNKNOWN accessions are treated like SEQUENCE (let NCBI
// surface the actual error) rather than rejected at Bind — that matches the
// pre-batch behavior of the single-accession path.
std::vector<ReadNCBIFastaTableFunction::WorkUnit> BuildWorkUnits(const std::vector<std::string> &accessions,
                                                                 int64_t batch_size) {
	using WorkUnit = ReadNCBIFastaTableFunction::WorkUnit;
	using Kind = ReadNCBIFastaTableFunction::WorkUnitKind;

	std::vector<WorkUnit> units;
	WorkUnit current_seq_batch {Kind::SEQUENCE_BATCH, {}};

	auto flush_seq = [&]() {
		if (!current_seq_batch.accessions.empty()) {
			units.push_back(std::move(current_seq_batch));
			current_seq_batch = WorkUnit {Kind::SEQUENCE_BATCH, {}};
		}
	};

	for (const auto &acc : accessions) {
		if (miint::NCBIParser::IsAssemblyAccession(acc)) {
			flush_seq();
			units.push_back(WorkUnit {Kind::ASSEMBLY, {acc}});
		} else {
			current_seq_batch.accessions.push_back(acc);
			if (static_cast<int64_t>(current_seq_batch.accessions.size()) >= batch_size) {
				flush_seq();
			}
		}
	}
	flush_seq();
	return units;
}

} // namespace

// Data constructor - sets up schema matching read_fastx
ReadNCBIFastaTableFunction::Data::Data(std::vector<std::string> accessions, std::vector<WorkUnit> work_units,
                                       const std::string &api_key, bool include_filepath, int64_t batch_size)
    : accessions(std::move(accessions)), work_units(std::move(work_units)), api_key(api_key),
      include_filepath(include_filepath), batch_size(batch_size) {
	// Schema matches read_fastx exactly
	names = {"sequence_index", "read_id", "comment", "sequence1", "sequence2", "qual1", "qual2"};
	types = {LogicalType::BIGINT,
	         LogicalType::VARCHAR,
	         LogicalType::VARCHAR,
	         LogicalType::VARCHAR,
	         LogicalType::VARCHAR,
	         LogicalType::LIST(LogicalType::UTINYINT),
	         LogicalType::LIST(LogicalType::UTINYINT)};

	if (include_filepath) {
		names.push_back("filepath");
		types.push_back(LogicalType::VARCHAR);
	}
}

// GlobalState constructor
ReadNCBIFastaTableFunction::GlobalState::GlobalState(DatabaseInstance &db, const std::string &api_key,
                                                     std::vector<WorkUnit> work_units)
    : work_units(std::move(work_units)), client(make_uniq<miint::NCBIClient>(db, api_key)), work_cursor(0),
      batch_offset(0), sequence_index(0) {
}

bool ReadNCBIFastaTableFunction::GlobalState::FetchNextBatch(ClientContext &context) {
	while (work_cursor < work_units.size()) {
		const auto &unit = work_units[work_cursor];
		work_cursor++;

		std::string fasta_text;
		std::vector<std::string> requested = unit.accessions;

		if (unit.kind == WorkUnitKind::ASSEMBLY) {
			fasta_text = client->FetchAssemblyFasta(unit.accessions.front());
			// Build the Datasets API download URL — that's the request the
			// fetch actually made, not an eutils URL.
			current_filepath_url = std::string(miint::NCBIClient::DATASETS_BASE) + "/genome/accession/" +
			                       unit.accessions.front() + "/download?include_annotation_type=GENOME_FASTA";
		} else {
			fasta_text = client->FetchFastaBatch(unit.accessions);
			// epost+efetch URLs need a fresh WebEnv handle that's already
			// expired by the time the user reads filepath, so surface the
			// equivalent `?id=A,B,C` form — same endpoint, reproducible.
			current_filepath_url = std::string(miint::NCBIClient::EUTILS_BASE) + "/efetch.fcgi?db=nuccore&id=" +
			                       miint::NCBIParser::JoinStrings(unit.accessions, ",") + "&rettype=fasta";
		}

		if (fasta_text.empty()) {
			// Empty response means NCBI dropped what we asked for. Record the
			// loss and continue rather than throw — one bad batch in a 10k-
			// accession run must not lose the remainder. The warning distinguishes
			// the two failure modes (sequence batch vs. single assembly) so the
			// user can tell whether it was epost+efetch or the Datasets API that
			// produced the empty body.
			for (const auto &acc : requested) {
				missing_accessions.push_back(acc);
			}
			if (unit.kind == WorkUnitKind::ASSEMBLY) {
				miint::EmitWarning(context,
				                   "read_ncbi_fasta: Datasets API returned empty FASTA for assembly accession '%s'",
				                   requested.front().c_str());
			} else {
				miint::EmitWarning(
				    context, "read_ncbi_fasta: NCBI returned empty response for sequence batch of %d accession(s): %s",
				    static_cast<int>(requested.size()), miint::NCBIParser::JoinStrings(requested, ",").c_str());
			}
			continue;
		}

		current_batch = miint::NCBIParser::ParseFasta(fasta_text);
		batch_offset = 0;

		// Diff requested vs returned for the sequence-batch path. Assembly FASTA
		// can legitimately return more or differently-named records (one ZIP
		// expands to many sequences), so diffing there would false-positive.
		if (unit.kind == WorkUnitKind::SEQUENCE_BATCH) {
			std::vector<std::string> returned_ids;
			returned_ids.reserve(current_batch.read_ids.size());
			for (const auto &id : current_batch.read_ids) {
				returned_ids.push_back(id);
			}
			auto missing = miint::NCBIParser::DiffMissingAccessions(requested, returned_ids);
			if (!missing.empty()) {
				miint::EmitWarning(context,
				                   "read_ncbi_fasta: NCBI omitted %d of %d requested accession(s) in batch: %s",
				                   static_cast<int>(missing.size()), static_cast<int>(requested.size()),
				                   miint::NCBIParser::JoinStrings(missing, ",").c_str());
				for (const auto &acc : missing) {
					missing_accessions.push_back(acc);
				}
			}
		}

		if (!current_batch.empty()) {
			return true;
		}
		// Empty parse on a non-empty response: continue to the next unit so a
		// single odd batch doesn't kill the whole scan.
	}

	// No more work. Emit end-of-scan summary exactly once. CAS guards against
	// re-emission if Execute is somehow called again after returning empty.
	bool expected = false;
	if (summary_emitted.compare_exchange_strong(expected, true)) {
		if (!missing_accessions.empty()) {
			miint::EmitWarning(context,
			                   "read_ncbi_fasta: SUMMARY: %d total accession(s) omitted by NCBI across this query: %s",
			                   static_cast<int>(missing_accessions.size()),
			                   miint::NCBIParser::JoinStrings(missing_accessions, ",").c_str());
		}
	}
	return false;
}

unique_ptr<FunctionData> ReadNCBIFastaTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                          vector<LogicalType> &return_types,
                                                          vector<std::string> &names) {
	// Parse accession(s) - can be VARCHAR or VARCHAR[]
	std::vector<std::string> accessions;

	if (input.inputs[0].IsNull()) {
		throw InvalidInputException(
		    "read_ncbi_fasta: accession cannot be NULL (use a VARCHAR literal or VARCHAR[] list)");
	}
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		accessions.push_back(input.inputs[0].ToString());
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			if (child.IsNull()) {
				throw InvalidInputException("read_ncbi_fasta: accession list cannot contain NULL");
			}
			accessions.push_back(child.ToString());
		}
	} else {
		throw InvalidInputException("read_ncbi_fasta: first argument must be VARCHAR or VARCHAR[]");
	}

	if (accessions.empty()) {
		throw InvalidInputException("read_ncbi_fasta: at least one accession must be provided");
	}

	// Validate that no accession is empty, and that none is actually a relation name
	// (#179 — passing a table of accessions by name used to be a silent no-op).
	for (const auto &acc : accessions) {
		if (acc.empty()) {
			throw InvalidInputException("read_ncbi_fasta: accession cannot be empty");
		}
		RejectRelationNameAsLiteral(context, "read_ncbi_fasta", acc);
	}

	// Parse api_key parameter (optional)
	std::string api_key;
	auto api_key_param = input.named_parameters.find("api_key");
	if (api_key_param != input.named_parameters.end() && !api_key_param->second.IsNull()) {
		api_key = api_key_param->second.ToString();
	}

	// Parse include_filepath parameter (optional)
	bool include_filepath = false;
	auto fp_param = input.named_parameters.find("include_filepath");
	if (fp_param != input.named_parameters.end() && !fp_param->second.IsNull()) {
		include_filepath = fp_param->second.GetValue<bool>();
	}

	// Parse batch_size parameter (optional). Default 500 matches NCBI's
	// practical POST size for epost; values <=0 are rejected loudly rather
	// than silently coerced — Rule 10.
	int64_t batch_size = 500;
	constexpr int64_t MAX_BATCH_SIZE = 10000; // Practical upper bound for NCBI epost POST bodies.
	auto bs_param = input.named_parameters.find("batch_size");
	if (bs_param != input.named_parameters.end() && !bs_param->second.IsNull()) {
		batch_size = bs_param->second.GetValue<int64_t>();
		if (batch_size <= 0) {
			throw InvalidInputException("read_ncbi_fasta: batch_size must be > 0 (got %lld)",
			                            static_cast<long long>(batch_size));
		}
		if (batch_size > MAX_BATCH_SIZE) {
			throw InvalidInputException("read_ncbi_fasta: batch_size %lld exceeds max %lld (NCBI epost limit)",
			                            static_cast<long long>(batch_size), static_cast<long long>(MAX_BATCH_SIZE));
		}
	}

	auto work_units = BuildWorkUnits(accessions, batch_size);
	auto data = make_uniq<Data>(std::move(accessions), std::move(work_units), api_key, include_filepath, batch_size);

	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> ReadNCBIFastaTableFunction::InitGlobal(ClientContext &context,
                                                                            TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto &db = DatabaseInstance::GetDatabase(context);
	return make_uniq<GlobalState>(db, data.api_key, data.work_units);
}

unique_ptr<LocalTableFunctionState> ReadNCBIFastaTableFunction::InitLocal(ExecutionContext &context,
                                                                          TableFunctionInitInput &input,
                                                                          GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

void ReadNCBIFastaTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	// Lock held across HTTP fetches inside FetchNextBatch. Intentional given
	// MaxThreads()==1 — the rate limiter already serializes calls anyway, and
	// keeping the lock simplifies the cursor/batch state machine. If this ever
	// goes multi-threaded, FetchNextBatch must move its HTTP work outside the
	// lock and use the rate limiter as the only serialization point.
	lock_guard<mutex> guard(global_state.lock);

	// If current batch is exhausted, fetch next work unit.
	while (global_state.current_batch.empty() || global_state.batch_offset >= global_state.current_batch.size()) {
		if (!global_state.FetchNextBatch(context)) {
			output.SetCardinality(0);
			return;
		}
	}

	// Determine how many records to output
	idx_t remaining = global_state.current_batch.size() - global_state.batch_offset;
	idx_t count = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);

	// Fill output vectors
	auto &batch = global_state.current_batch;
	size_t offset = global_state.batch_offset;

	// sequence_index (column 0)
	auto seq_idx_data = FlatVector::GetData<int64_t>(output.data[0]);
	for (idx_t i = 0; i < count; i++) {
		seq_idx_data[i] = global_state.sequence_index++;
	}

	// read_id (column 1)
	auto read_id_data = FlatVector::GetData<string_t>(output.data[1]);
	for (idx_t i = 0; i < count; i++) {
		read_id_data[i] = StringVector::AddString(output.data[1], batch.read_ids[offset + i]);
	}

	// comment (column 2) - nullable
	auto comment_data = FlatVector::GetData<string_t>(output.data[2]);
	auto &comment_validity = FlatVector::Validity(output.data[2]);
	for (idx_t i = 0; i < count; i++) {
		const auto &comment = batch.comments[offset + i];
		comment_data[i] = StringVector::AddString(output.data[2], comment);
		if (comment.empty()) {
			comment_validity.SetInvalid(i);
		}
	}

	// sequence1 (column 3)
	auto seq1_data = FlatVector::GetData<string_t>(output.data[3]);
	for (idx_t i = 0; i < count; i++) {
		seq1_data[i] = StringVector::AddString(output.data[3], batch.sequences1[offset + i]);
	}

	// sequence2 (column 4) - always NULL for NCBI FASTA
	output.data[4].SetVectorType(VectorType::CONSTANT_VECTOR);
	ConstantVector::SetNull(output.data[4], true);

	// qual1 (column 5) - always NULL for NCBI FASTA (no quality scores)
	output.data[5].SetVectorType(VectorType::CONSTANT_VECTOR);
	ConstantVector::SetNull(output.data[5], true);

	// qual2 (column 6) - always NULL
	output.data[6].SetVectorType(VectorType::CONSTANT_VECTOR);
	ConstantVector::SetNull(output.data[6], true);

	// filepath (column 7) - optional. Stores whichever URL the current work
	// unit actually used (eutils efetch for sequences, Datasets API for
	// assemblies) — set by FetchNextBatch when the unit was dispatched.
	if (bind_data.include_filepath) {
		output.data[7].SetVectorType(VectorType::CONSTANT_VECTOR);
		auto filepath_data = ConstantVector::GetData<string_t>(output.data[7]);
		*filepath_data = StringVector::AddString(output.data[7], global_state.current_filepath_url);
	}

	global_state.batch_offset += count;
	output.SetCardinality(count);
}

TableFunction ReadNCBIFastaTableFunction::GetFunction() {
	auto tf = TableFunction("read_ncbi_fasta", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["api_key"] = LogicalType::VARCHAR;
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	tf.named_parameters["batch_size"] = LogicalType::BIGINT;
	return tf;
}

void ReadNCBIFastaTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
