// sylph_index_create() — build a sylph `.syldb` from a reference-sequence table.
// See sylph_index_create.hpp for the contract. Like rype_index_create, the build
// is a synchronous side effect in InitGlobal; Execute emits one status row.

#ifdef MIINT_HAS_SYLPH

#include "sylph_index_create.hpp"

#include "sequence_table_reader.hpp"
#include "sylph.h"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/materialized_query_result.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <string>
#include <vector>

namespace duckdb {

namespace {

// RAII guard so the FFI builder is freed on any throw between create and the
// normal end-of-InitGlobal free.
struct IndexBuilderGuard {
	::SylphIndexBuilder *builder = nullptr;
	~IndexBuilderGuard() {
		if (builder) {
			sylph_index_builder_free(builder);
		}
	}
};

[[noreturn]] void ThrowFFI(const char *prefix) {
	const char *err = sylph_get_last_error();
	throw IOException("sylph_index_create: %s: %s", prefix, err ? err : "<unknown>");
}

// Read one named INTEGER parameter into a bounded FFI field. Validates the
// value is in [1, hi] (0 is reserved for "use sylph default" and is left as the
// struct's zero, so callers never pass 0 explicitly here).
template <typename T>
void ApplyBoundedInt(TableFunctionBindInput &input, const std::string &param, int64_t hi, T &out) {
	auto it = input.named_parameters.find(param);
	if (it == input.named_parameters.end() || it->second.IsNull()) {
		return;
	}
	int64_t v = it->second.GetValue<int64_t>();
	if (v < 1 || v > hi) {
		throw BinderException("sylph_index_create: %s must be in [1, %lld] (got %lld)", param.c_str(), (long long)hi,
		                      (long long)v);
	}
	out = static_cast<T>(v);
}

} // namespace

// =============================================================================
// Bind
// =============================================================================
unique_ptr<FunctionData> SylphIndexCreateTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                             vector<LogicalType> &return_types,
                                                             vector<std::string> &names) {
	auto data = make_uniq<Data>();

	if (input.inputs.size() < 2) {
		throw BinderException("sylph_index_create requires source_table and output_path parameters");
	}
	data->source_table = input.inputs[0].ToString();
	data->output_path = input.inputs[1].ToString();

	// genome_id is required — it defines what a genome is.
	auto genome_param = input.named_parameters.find("genome_id");
	if (genome_param == input.named_parameters.end() || genome_param->second.IsNull() ||
	    genome_param->second.ToString().empty()) {
		throw BinderException("sylph_index_create: the genome_id parameter (grouping column) is required");
	}
	data->genome_id_col = genome_param->second.ToString();

	auto order_param = input.named_parameters.find("order_by");
	if (order_param != input.named_parameters.end() && !order_param->second.IsNull() &&
	    !order_param->second.ToString().empty()) {
		data->order_by_col = order_param->second.ToString();
	}

	// Seed sketch params from sylph's defaults — the single source of truth.
	if (sylph_genome_sketch_params_default(&data->sketch_params) != 0) {
		throw IOException("sylph_index_create: sylph_genome_sketch_params_default failed");
	}
	ApplyBoundedInt(input, "k", 255, data->sketch_params.k);
	ApplyBoundedInt(input, "c", 65535, data->sketch_params.c);
	ApplyBoundedInt(input, "min_spacing", 4294967295LL, data->sketch_params.min_spacing);
	// k is further constrained to sylph's supported sizes. The FFI re-checks
	// this, but validating here fails before the (non-atomic) build runs.
	if (data->sketch_params.k != 0 && data->sketch_params.k != 21 && data->sketch_params.k != 31) {
		throw BinderException("sylph_index_create: k must be 21 or 31 (got %d)", (int)data->sketch_params.k);
	}

	auto pseudotax_param = input.named_parameters.find("pseudotax");
	if (pseudotax_param != input.named_parameters.end() && !pseudotax_param->second.IsNull()) {
		data->sketch_params.pseudotax = pseudotax_param->second.GetValue<bool>() ? 1 : 0;
	}

	// Fail-fast schema validation (before the build side effect). read_id +
	// sequence1 are required; read_id may be VARCHAR or BIGINT.
	ValidateSequenceTableSchema(context, data->source_table, /*allow_bigint=*/true);

	// Validate the genome_id and order_by columns exist. A `LIMIT 0` probe works
	// for both tables and views and reports missing columns as a clean error.
	{
		Connection conn(*context.db);
		auto src = KeywordHelper::WriteOptionallyQuoted(data->source_table);
		auto gcol = KeywordHelper::WriteOptionallyQuoted(data->genome_id_col);
		auto ocol = KeywordHelper::WriteOptionallyQuoted(data->order_by_col);
		auto probe = conn.Query("SELECT " + gcol + ", " + ocol + " FROM " + src + " LIMIT 0");
		if (probe->HasError()) {
			throw BinderException("sylph_index_create: genome_id/order_by column check failed: %s", probe->GetError());
		}
		// Detect an optional `comment` column (present in read_fastx output). Its
		// presence switches on full-header contig-name reconstruction below.
		auto comment_probe = conn.Query("SELECT comment FROM " + src + " LIMIT 0");
		data->has_comment = !comment_probe->HasError();
	}

	names = data->names;
	return_types = data->types;
	return std::move(data);
}

// =============================================================================
// InitGlobal — performs the build (a synchronous side effect)
// =============================================================================
unique_ptr<GlobalTableFunctionState> SylphIndexCreateTableFunction::InitGlobal(ClientContext &context,
                                                                               TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	IndexBuilderGuard guard;
	guard.builder = sylph_index_builder_create(&data.sketch_params);
	if (guard.builder == nullptr) {
		ThrowFFI("index builder create failed");
	}

	Connection conn(*context.db);
	auto src = KeywordHelper::WriteOptionallyQuoted(data.source_table);
	auto gcol = KeywordHelper::WriteOptionallyQuoted(data.genome_id_col);
	auto ocol = KeywordHelper::WriteOptionallyQuoted(data.order_by_col);

	// Distinct genome keys (as VARCHAR), stably ordered so the resulting
	// genome_index ordering in the .syldb is reproducible.
	auto keys_res = conn.Query("SELECT DISTINCT CAST(" + gcol + " AS VARCHAR) AS gk FROM " + src + " WHERE " + gcol +
	                           " IS NOT NULL ORDER BY gk");
	if (keys_res->HasError()) {
		throw InvalidInputException("sylph_index_create: failed to enumerate genomes: %s", keys_res->GetError());
	}
	std::vector<std::string> genome_keys;
	for (idx_t r = 0; r < keys_res->RowCount(); r++) {
		auto v = keys_res->GetValue(0, r);
		if (v.IsNull()) {
			continue;
		}
		genome_keys.push_back(v.ToString());
	}
	if (genome_keys.empty()) {
		throw InvalidInputException("sylph_index_create: no non-NULL genome_id values found in '%s'",
		                            data.source_table.c_str());
	}

	// One genome at a time: begin → stream its contigs in order → end.
	for (const auto &key : genome_keys) {
		if (sylph_index_builder_begin_genome(guard.builder, key.c_str()) != 0) {
			ThrowFFI("begin_genome failed");
		}

		// Contig name = the full FASTA header. sylph (via needletail) stores the
		// whole header line as first_contig_name; read_fastx splits it into
		// read_id (first token) + comment (remainder). When a comment column is
		// present we rejoin them ("read_id comment") so a miint-built .syldb's
		// contig names match `sylph sketch`; otherwise the contig name is read_id.
		std::string contig_expr = data.has_comment
		                              ? "CAST(read_id AS VARCHAR) || CASE WHEN comment IS NOT NULL AND comment <> '' "
		                                "THEN ' ' || CAST(comment AS VARCHAR) ELSE '' END"
		                              : "CAST(read_id AS VARCHAR)";

		auto key_lit = Value(key).ToSQLString();
		auto sql = "SELECT " + contig_expr + " AS contig, sequence1 AS seq FROM " + src + " WHERE CAST(" + gcol +
		           " AS VARCHAR) = " + key_lit + " AND sequence1 IS NOT NULL ORDER BY " + ocol;
		auto stream = conn.SendQuery(sql);
		if (stream->HasError()) {
			throw InvalidInputException("sylph_index_create: failed to read genome '%s': %s", key.c_str(),
			                            stream->GetError());
		}

		while (true) {
			auto chunk = stream->Fetch();
			if (!chunk || chunk->size() == 0) {
				break;
			}
			auto &contig_vec = chunk->data[0];
			auto &seq_vec = chunk->data[1];
			contig_vec.Flatten(chunk->size());
			seq_vec.Flatten(chunk->size());
			auto contig_data = FlatVector::GetData<string_t>(contig_vec);
			auto &contig_valid = FlatVector::Validity(contig_vec);
			auto seq_data = FlatVector::GetData<string_t>(seq_vec);
			auto &seq_valid = FlatVector::Validity(seq_vec);

			for (idx_t i = 0; i < chunk->size(); i++) {
				if (!seq_valid.RowIsValid(i)) {
					continue; // SQL already filters NULL seq; defensive.
				}
				// contig_name needs NUL-termination for the C FFI; seq is passed
				// as (ptr, len) so its bytes need no terminator.
				std::string contig = contig_valid.RowIsValid(i) ? contig_data[i].GetString() : std::string();
				auto seq = seq_data[i];
				int rc = sylph_index_builder_add_contig(guard.builder, contig.c_str(),
				                                        reinterpret_cast<const unsigned char *>(seq.GetData()),
				                                        seq.GetSize());
				if (rc != 0) {
					ThrowFFI("add_contig failed");
				}
			}
		}

		if (sylph_index_builder_end_genome(guard.builder) != 0) {
			ThrowFFI("end_genome failed");
		}
	}

	if (sylph_index_builder_write(guard.builder, data.output_path.c_str()) != 0) {
		ThrowFFI("failed to write database");
	}
	gstate->num_genomes = sylph_index_builder_num_genomes(guard.builder);
	// guard frees the builder on return.
	return std::move(gstate);
}

unique_ptr<LocalTableFunctionState>
SylphIndexCreateTableFunction::InitLocal(ExecutionContext &, TableFunctionInitInput &, GlobalTableFunctionState *) {
	return make_uniq<LocalState>();
}

// =============================================================================
// Execute — emit the single status row
// =============================================================================
void SylphIndexCreateTableFunction::Execute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	if (gstate.done) {
		output.SetCardinality(0);
		return;
	}

	// Resolve the reported k/c from the sketch params (0 = sylph default).
	int32_t k = data.sketch_params.k != 0 ? data.sketch_params.k : 31;
	int32_t c = data.sketch_params.c != 0 ? data.sketch_params.c : 200;

	output.data[0].SetValue(0, Value(data.output_path));
	output.data[1].SetValue(0, Value::INTEGER(k));
	output.data[2].SetValue(0, Value::INTEGER(c));
	output.data[3].SetValue(0, Value::UBIGINT(gstate.num_genomes));
	output.data[4].SetValue(0, Value("ok"));

	output.SetCardinality(1);
	gstate.done = true;
}

// =============================================================================
// Registration
// =============================================================================
TableFunction SylphIndexCreateTableFunction::GetFunction() {
	TableFunction tf("sylph_index_create", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal,
	                 InitLocal);

	tf.named_parameters["genome_id"] = LogicalType::VARCHAR;
	tf.named_parameters["order_by"] = LogicalType::VARCHAR;
	tf.named_parameters["k"] = LogicalType::INTEGER;
	tf.named_parameters["c"] = LogicalType::INTEGER;
	tf.named_parameters["min_spacing"] = LogicalType::INTEGER;
	tf.named_parameters["pseudotax"] = LogicalType::BOOLEAN;

	return tf;
}

void SylphIndexCreateTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb

#endif // MIINT_HAS_SYLPH
