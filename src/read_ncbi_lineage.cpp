#include "read_ncbi_lineage.hpp"

#include "ensure_httpfs.hpp"
#include "miint_log.hpp"
#include "duckdb/common/vector_size.hpp"

#include <cstdlib>
#include <unordered_set>

namespace duckdb {

namespace {

std::vector<std::vector<std::string>> ChunkBatches(const std::vector<std::string> &taxids, int64_t batch_size) {
	std::vector<std::vector<std::string>> batches;
	std::vector<std::string> current;
	for (const auto &t : taxids) {
		current.push_back(t);
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

// A taxid argument is an integer; accept numeric values (normalized via int64) and
// VARCHAR (already the textual taxid).
std::string TaxidToString(const Value &v) {
	if (v.type().id() == LogicalTypeId::VARCHAR) {
		return v.ToString();
	}
	return std::to_string(v.GetValue<int64_t>());
}

// Emit a VARCHAR cell, or NULL when the collapsed rank is absent.
void SetRankCell(DataChunk &output, idx_t col, idx_t row, const std::string &value) {
	if (value.empty()) {
		FlatVector::Validity(output.data[col]).SetInvalid(row);
	} else {
		FlatVector::GetData<string_t>(output.data[col])[row] = StringVector::AddString(output.data[col], value);
	}
}

} // namespace

ReadNCBILineageTableFunction::Data::Data(std::vector<std::string> taxids, std::vector<std::vector<std::string>> batches,
                                         const std::string &api_key, int64_t batch_size)
    : taxids(std::move(taxids)), batches(std::move(batches)), api_key(api_key), batch_size(batch_size) {
	names = {"taxid", "name",   "rank",  "domain",  "phylum", "class",
	         "order", "family", "genus", "species", "strain", "lineage"};
	types = {LogicalType::BIGINT,  LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
	         LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
	         LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR};
}

ReadNCBILineageTableFunction::GlobalState::GlobalState(DatabaseInstance &db, const std::string &api_key,
                                                       std::vector<std::vector<std::string>> batches)
    : batches(std::move(batches)), client(make_uniq<miint::NCBIClient>(db, api_key)), batch_cursor(0),
      result_offset(0) {
}

bool ReadNCBILineageTableFunction::GlobalState::FetchNextBatch(ClientContext &context) {
	while (batch_cursor < batches.size()) {
		const auto &batch_ids = batches[batch_cursor];
		batch_cursor++;

		std::string xml = client->FetchTaxonomyXMLBatch(batch_ids);
		if (xml.empty()) {
			for (const auto &t : batch_ids) {
				missing_taxids.push_back(t);
			}
			miint::EmitWarning(context, "read_ncbi_lineage: NCBI returned empty response for batch of %d taxid(s): %s",
			                   static_cast<int>(batch_ids.size()),
			                   miint::NCBIParser::JoinStrings(batch_ids, ",").c_str());
			continue;
		}

		auto parsed = miint::ParseTaxonomyXML(xml);

		std::unordered_set<int64_t> returned;
		for (const auto &lin : parsed) {
			returned.insert(lin.taxid);
			// A queried taxid that NCBI has merged resolves to a taxon with a
			// different (current) TaxId; its old id is listed under <AkaTaxIds>.
			// Count those so a merged-but-resolved taxid isn't reported as omitted.
			for (int64_t aka : lin.aka_taxids) {
				returned.insert(aka);
			}
		}
		std::vector<std::string> missing;
		for (const auto &t : batch_ids) {
			if (returned.find(static_cast<int64_t>(std::strtoll(t.c_str(), nullptr, 10))) == returned.end()) {
				missing.push_back(t);
			}
		}
		if (!missing.empty()) {
			miint::EmitWarning(context, "read_ncbi_lineage: NCBI omitted %d of %d requested taxid(s) in batch: %s",
			                   static_cast<int>(missing.size()), static_cast<int>(batch_ids.size()),
			                   miint::NCBIParser::JoinStrings(missing, ",").c_str());
			for (const auto &t : missing) {
				missing_taxids.push_back(t);
			}
		}

		if (!parsed.empty()) {
			for (auto &lin : parsed) {
				results.push_back(std::move(lin));
			}
			return true;
		}
	}

	bool expected = false;
	if (summary_emitted.compare_exchange_strong(expected, true)) {
		if (!missing_taxids.empty()) {
			miint::EmitWarning(
			    context, "read_ncbi_lineage: SUMMARY: %d total taxid(s) omitted by NCBI across this query: %s",
			    static_cast<int>(missing_taxids.size()), miint::NCBIParser::JoinStrings(missing_taxids, ",").c_str());
		}
	}
	return false;
}

unique_ptr<FunctionData> ReadNCBILineageTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                            vector<LogicalType> &return_types,
                                                            vector<std::string> &names) {
	std::vector<std::string> taxids;

	if (input.inputs[0].IsNull()) {
		throw InvalidInputException("read_ncbi_lineage: taxid cannot be NULL");
	}
	if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : children) {
			if (child.IsNull()) {
				throw InvalidInputException("read_ncbi_lineage: taxid list cannot contain NULL");
			}
			taxids.push_back(TaxidToString(child));
		}
	} else {
		taxids.push_back(TaxidToString(input.inputs[0]));
	}

	if (taxids.empty()) {
		throw InvalidInputException("read_ncbi_lineage: at least one taxid must be provided");
	}

	std::string api_key;
	auto api_key_param = input.named_parameters.find("api_key");
	if (api_key_param != input.named_parameters.end() && !api_key_param->second.IsNull()) {
		api_key = api_key_param->second.ToString();
	}

	int64_t batch_size = 500;
	constexpr int64_t MAX_BATCH_SIZE = 10000; // NCBI epost limit
	auto bs_param = input.named_parameters.find("batch_size");
	if (bs_param != input.named_parameters.end() && !bs_param->second.IsNull()) {
		batch_size = bs_param->second.GetValue<int64_t>();
		if (batch_size <= 0) {
			throw InvalidInputException("read_ncbi_lineage: batch_size must be > 0 (got %lld)",
			                            static_cast<long long>(batch_size));
		}
		if (batch_size > MAX_BATCH_SIZE) {
			throw InvalidInputException("read_ncbi_lineage: batch_size %lld exceeds max %lld (NCBI epost limit)",
			                            static_cast<long long>(batch_size), static_cast<long long>(MAX_BATCH_SIZE));
		}
	}

	auto batches = ChunkBatches(taxids, batch_size);
	auto data = make_uniq<Data>(std::move(taxids), std::move(batches), api_key, batch_size);

	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ReadNCBILineageTableFunction::InitGlobal(ClientContext &context,
                                                                              TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	// The NCBI E-utilities client uses HTTPUtil, whose https support comes from httpfs.
	miint::EnsureHttpfsLoaded(context);
	auto &db = DatabaseInstance::GetDatabase(context);
	return make_uniq<GlobalState>(db, data.api_key, data.batches);
}

unique_ptr<LocalTableFunctionState>
ReadNCBILineageTableFunction::InitLocal(ExecutionContext &, TableFunctionInitInput &, GlobalTableFunctionState *) {
	return make_uniq<LocalState>();
}

void ReadNCBILineageTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	lock_guard<mutex> guard(global_state.lock);

	// Drain-then-refill so memory stays bounded to one batch (see read_ncbi).
	while (global_state.result_offset >= global_state.results.size()) {
		global_state.results.clear();
		global_state.result_offset = 0;
		if (!global_state.FetchNextBatch(context)) {
			output.SetCardinality(0);
			return;
		}
	}

	idx_t remaining = global_state.results.size() - global_state.result_offset;
	idx_t count = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);
	size_t offset = global_state.result_offset;

	for (idx_t i = 0; i < count; i++) {
		const auto &lin = global_state.results[offset + i];
		FlatVector::GetData<int64_t>(output.data[0])[i] = lin.taxid;
		FlatVector::GetData<string_t>(output.data[1])[i] = StringVector::AddString(output.data[1], lin.name);
		FlatVector::GetData<string_t>(output.data[2])[i] = StringVector::AddString(output.data[2], lin.rank);
		SetRankCell(output, 3, i, lin.domain);
		SetRankCell(output, 4, i, lin.phylum);
		SetRankCell(output, 5, i, lin.tax_class);
		SetRankCell(output, 6, i, lin.tax_order);
		SetRankCell(output, 7, i, lin.family);
		SetRankCell(output, 8, i, lin.genus);
		SetRankCell(output, 9, i, lin.species);
		SetRankCell(output, 10, i, lin.strain);
		FlatVector::GetData<string_t>(output.data[11])[i] = StringVector::AddString(output.data[11], lin.lineage);
	}

	global_state.result_offset += count;
	output.SetCardinality(count);
}

TableFunction ReadNCBILineageTableFunction::GetFunction() {
	auto tf = TableFunction("read_ncbi_lineage", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["api_key"] = LogicalType::VARCHAR;
	tf.named_parameters["batch_size"] = LogicalType::BIGINT;
	return tf;
}

void ReadNCBILineageTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
