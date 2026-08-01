#include "read_ncbi.hpp"
#include "catalog_utils.hpp"
#include "miint_log.hpp"
#include "duckdb/common/vector_size.hpp"
#include <sstream>

namespace duckdb {

namespace {

// Chunk an in-order accession list into fixed-size batches. read_ncbi rejects
// assembly accessions at Bind, so unlike read_ncbi_fasta there's no per-type
// partition here — every input goes through epost+efetch as plain sequence IDs.
std::vector<std::vector<std::string>> ChunkBatches(const std::vector<std::string> &accessions, int64_t batch_size) {
	std::vector<std::vector<std::string>> batches;
	std::vector<std::string> current;
	for (const auto &acc : accessions) {
		current.push_back(acc);
		if (static_cast<int64_t>(current.size()) >= batch_size) {
			batches.push_back(std::move(current));
			current.clear();
		}
	}
	if (!current.empty()) {
		batches.push_back(std::move(current));
	}
	return batches;
}

} // namespace

ReadNCBITableFunction::Data::Data(std::vector<std::string> accessions, std::vector<std::vector<std::string>> batches,
                                  const std::string &api_key, int64_t batch_size)
    : accessions(std::move(accessions)), batches(std::move(batches)), api_key(api_key), batch_size(batch_size) {
	names = {"accession",   "version", "description",   "organism",
	         "taxonomy_id", "length",  "molecule_type", "update_date"};
	types = {LogicalType::VARCHAR, // accession
	         LogicalType::INTEGER, // version
	         LogicalType::VARCHAR, // description
	         LogicalType::VARCHAR, // organism
	         LogicalType::BIGINT,  // taxonomy_id
	         LogicalType::BIGINT,  // length
	         LogicalType::VARCHAR, // molecule_type
	         LogicalType::DATE};   // update_date
}

ReadNCBITableFunction::GlobalState::GlobalState(DatabaseInstance &db, const std::string &api_key,
                                                std::vector<std::vector<std::string>> batches)
    : batches(std::move(batches)), client(make_uniq<miint::NCBIClient>(db, api_key)), batch_cursor(0),
      result_offset(0) {
}

bool ReadNCBITableFunction::GlobalState::FetchNextBatch(ClientContext &context) {
	while (batch_cursor < batches.size()) {
		const auto &batch_ids = batches[batch_cursor];
		batch_cursor++;

		std::string xml = client->FetchGenBankXMLBatch(batch_ids);

		if (xml.empty()) {
			for (const auto &acc : batch_ids) {
				missing_accessions.push_back(acc);
			}
			miint::EmitWarning(context, "read_ncbi: NCBI returned empty response for batch of %d accession(s): %s",
			                   static_cast<int>(batch_ids.size()),
			                   miint::NCBIParser::JoinStrings(batch_ids, ",").c_str());
			continue;
		}

		auto parsed = miint::NCBIParser::ParseGenBankXMLBatch(xml);

		// DiffMissingAccessions compares base IDs after StripAccessionVersion,
		// so the version mismatch between user-supplied "NC_001416" and the
		// parser's "NC_001416.1" doesn't false-positive as missing.
		std::vector<std::string> returned_ids;
		returned_ids.reserve(parsed.size());
		for (const auto &meta : parsed) {
			if (!meta.accession.empty()) {
				returned_ids.push_back(meta.accession);
			}
		}
		auto missing = miint::NCBIParser::DiffMissingAccessions(batch_ids, returned_ids);
		if (!missing.empty()) {
			miint::EmitWarning(context, "read_ncbi: NCBI omitted %d of %d requested accession(s) in batch: %s",
			                   static_cast<int>(missing.size()), static_cast<int>(batch_ids.size()),
			                   miint::NCBIParser::JoinStrings(missing, ",").c_str());
			for (const auto &acc : missing) {
				missing_accessions.push_back(acc);
			}
		}

		bool any_appended = false;
		for (auto &meta : parsed) {
			if (!meta.accession.empty()) {
				metadata_results.push_back(std::move(meta));
				any_appended = true;
			}
		}
		if (any_appended) {
			return true;
		}
	}

	bool expected = false;
	if (summary_emitted.compare_exchange_strong(expected, true)) {
		if (!missing_accessions.empty()) {
			miint::EmitWarning(context,
			                   "read_ncbi: SUMMARY: %d total accession(s) omitted by NCBI across this query: %s",
			                   static_cast<int>(missing_accessions.size()),
			                   miint::NCBIParser::JoinStrings(missing_accessions, ",").c_str());
		}
	}
	return false;
}

unique_ptr<FunctionData> ReadNCBITableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                     vector<LogicalType> &return_types, vector<std::string> &names) {
	std::vector<std::string> accessions;

	if (input.inputs[0].IsNull()) {
		throw InvalidInputException("read_ncbi: accession cannot be NULL");
	}
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		accessions.push_back(input.inputs[0].ToString());
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			if (child.IsNull()) {
				throw InvalidInputException("read_ncbi: accession list cannot contain NULL");
			}
			accessions.push_back(child.ToString());
		}
	} else {
		throw InvalidInputException("read_ncbi: first argument must be VARCHAR or VARCHAR[]");
	}

	if (accessions.empty()) {
		throw InvalidInputException("read_ncbi: at least one accession must be provided");
	}

	for (const auto &acc : accessions) {
		if (acc.empty()) {
			throw InvalidInputException("read_ncbi: accession cannot be empty");
		}
		// #179 — passing a table of accessions by name used to be a silent no-op.
		RejectRelationNameAsLiteral(context, "read_ncbi", acc);
		if (miint::NCBIParser::IsAssemblyAccession(acc)) {
			throw InvalidInputException("read_ncbi: Assembly accession '%s' is not supported. "
			                            "Assembly accessions (GCF_/GCA_) represent collections of sequences. "
			                            "Use read_ncbi_fasta('%s') to retrieve sequences from this assembly, "
			                            "or find the component RefSeq accessions (e.g., NC_XXXXXX.X) for metadata.",
			                            acc, acc);
		}
	}

	std::string api_key;
	auto api_key_param = input.named_parameters.find("api_key");
	if (api_key_param != input.named_parameters.end() && !api_key_param->second.IsNull()) {
		api_key = api_key_param->second.ToString();
	}

	int64_t batch_size = 500;
	constexpr int64_t MAX_BATCH_SIZE = 10000; // Practical upper bound for NCBI epost POST bodies.
	auto bs_param = input.named_parameters.find("batch_size");
	if (bs_param != input.named_parameters.end() && !bs_param->second.IsNull()) {
		batch_size = bs_param->second.GetValue<int64_t>();
		if (batch_size <= 0) {
			throw InvalidInputException("read_ncbi: batch_size must be > 0 (got %lld)",
			                            static_cast<long long>(batch_size));
		}
		if (batch_size > MAX_BATCH_SIZE) {
			throw InvalidInputException("read_ncbi: batch_size %lld exceeds max %lld (NCBI epost limit)",
			                            static_cast<long long>(batch_size), static_cast<long long>(MAX_BATCH_SIZE));
		}
	}

	auto batches = ChunkBatches(accessions, batch_size);
	auto data = make_uniq<Data>(std::move(accessions), std::move(batches), api_key, batch_size);

	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> ReadNCBITableFunction::InitGlobal(ClientContext &context,
                                                                       TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto &db = DatabaseInstance::GetDatabase(context);
	return make_uniq<GlobalState>(db, data.api_key, data.batches);
}

unique_ptr<LocalTableFunctionState> ReadNCBITableFunction::InitLocal(ExecutionContext &context,
                                                                     TableFunctionInitInput &input,
                                                                     GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

void ReadNCBITableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	lock_guard<mutex> guard(global_state.lock);

	// Drain-then-refill: once we've emitted every record in the current
	// metadata_results buffer, clear it before fetching the next batch.
	// Without this the vector accumulates every batch's parsed metadata for
	// the lifetime of the scan — a 10k-accession query would hold all 10k
	// records in memory even though Execute is supposed to stream chunks.
	while (global_state.result_offset >= global_state.metadata_results.size()) {
		global_state.metadata_results.clear();
		global_state.result_offset = 0;
		if (!global_state.FetchNextBatch(context)) {
			output.SetCardinality(0);
			return;
		}
	}

	idx_t remaining = global_state.metadata_results.size() - global_state.result_offset;
	idx_t count = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);

	size_t offset = global_state.result_offset;

	for (idx_t i = 0; i < count; i++) {
		const auto &meta = global_state.metadata_results[offset + i];

		FlatVector::GetData<string_t>(output.data[0])[i] = StringVector::AddString(output.data[0], meta.accession);
		FlatVector::GetData<int32_t>(output.data[1])[i] = meta.version;
		FlatVector::GetData<string_t>(output.data[2])[i] = StringVector::AddString(output.data[2], meta.description);
		FlatVector::GetData<string_t>(output.data[3])[i] = StringVector::AddString(output.data[3], meta.organism);
		FlatVector::GetData<int64_t>(output.data[4])[i] = meta.taxonomy_id;
		FlatVector::GetData<int64_t>(output.data[5])[i] = meta.length;
		FlatVector::GetData<string_t>(output.data[6])[i] = StringVector::AddString(output.data[6], meta.molecule_type);

		if (!meta.update_date.empty()) {
			date_t date;
			idx_t pos;
			bool special;
			auto result = Date::TryConvertDate(meta.update_date.c_str(), meta.update_date.size(), pos, date, special);
			if (result == DateCastResult::SUCCESS) {
				FlatVector::GetData<date_t>(output.data[7])[i] = date;
			} else {
				FlatVector::Validity(output.data[7]).SetInvalid(i);
			}
		} else {
			FlatVector::Validity(output.data[7]).SetInvalid(i);
		}
	}

	global_state.result_offset += count;
	output.SetCardinality(count);
}

TableFunction ReadNCBITableFunction::GetFunction() {
	auto tf = TableFunction("read_ncbi", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["api_key"] = LogicalType::VARCHAR;
	tf.named_parameters["batch_size"] = LogicalType::BIGINT;
	return tf;
}

void ReadNCBITableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
