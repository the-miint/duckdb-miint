#include "miint_named_parameter.hpp"
// sylph_index_create() — build a sylph `.syldb` / `.syl2db` from a reference-sequence table.
// See sylph_index_create.hpp for the contract. The build is a synchronous side
// effect in InitGlobal: the distinct genome ids are enumerated, then N worker
// threads each claim genomes off a shared counter and sketch one genome at a time
// via a `WHERE genome_id = <id>` query (so only ~one genome is resident per
// thread, and a clustered source is pruned to that genome's row groups); each
// worker sketches into its own builder, the builders are merged, and the database
// is written once. Execute then emits the single status row.

#ifdef MIINT_HAS_SYLPH

#include "sylph_index_create.hpp"

#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "sequence_table_reader.hpp"
#include "sylph.h"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/parallel/task_scheduler.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <atomic>
#include <limits>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace duckdb {

namespace {

[[noreturn]] void ThrowFFI(const char *prefix) {
	const char *err = sylph_get_last_error();
	throw IOException("sylph_index_create: %s: %s", prefix, err ? err : "<unknown>");
}

// Read one named INTEGER parameter into a bounded FFI field. Validates the
// value is in [1, hi] (0 is reserved for "use sylph default" and is left as the
// struct's zero, so callers never pass 0 explicitly here).
template <typename T>
void ApplyBoundedInt(TableFunctionBindInput &input, const std::string &param, int64_t hi, T &out) {
	auto it = input.named_parameters.find(Identifier(param));
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
                                                             vector<Identifier> &names) {
	auto data = make_uniq<Data>();

	if (input.inputs.size() < 2) {
		throw BinderException("sylph_index_create requires source_table and output_path parameters");
	}
	data->source_table = input.inputs[0].ToString();
	RejectCTERelationName(input, data->source_table);
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

	// Genome-sketch parallelism. 0 = auto (DuckDB scheduler thread count).
	auto threads_param = input.named_parameters.find("threads");
	if (threads_param != input.named_parameters.end() && !threads_param->second.IsNull()) {
		auto t = threads_param->second.GetValue<int64_t>();
		if (t < 0) {
			throw BinderException("sylph_index_create: threads must be >= 0 (got %lld)", (long long)t);
		}
		data->user_threads = static_cast<uint32_t>(t);
	}

	// Output format: an explicit `two_stage` parameter decides; without it the
	// output_path suffix is the default (`.syl2db` = two-stage), as with
	// DuckDB's COPY and the sylph CLI's converter. The converter's own
	// preconditions are checked here so the (non-atomic) build never runs for
	// a database that could not be written; the FFI re-checks them.
	auto two_stage_param = input.named_parameters.find("two_stage");
	if (two_stage_param != input.named_parameters.end() && !two_stage_param->second.IsNull()) {
		data->two_stage = two_stage_param->second.GetValue<bool>();
	} else {
		data->two_stage = StringUtil::EndsWith(data->output_path, ".syl2db");
	}
	if (sylph_two_stage_params_default(&data->two_stage_params) != 0) {
		throw IOException("sylph_index_create: sylph_two_stage_params_default failed");
	}
	const bool has_two_stage_knobs =
	    input.named_parameters.contains("screen_c") || input.named_parameters.contains("min_sparse_kmers");
	if (data->two_stage) {
		ApplyBoundedInt(input, "screen_c", 4294967295LL, data->two_stage_params.screen_c);
		ApplyBoundedInt(input, "min_sparse_kmers", 4294967295LL, data->two_stage_params.min_sparse_kmers);
		if (data->sketch_params.pseudotax == 0) {
			throw BinderException("sylph_index_create: a two-stage (.syl2db) database needs profiling k-mers; "
			                      "pseudotax := false cannot be combined with a two-stage output");
		}
		// 0 = sylph's default c (200), the same resolution Execute reports.
		const uint32_t dense_c = data->sketch_params.c != 0 ? data->sketch_params.c : 200;
		if (data->two_stage_params.screen_c < dense_c) {
			throw BinderException("sylph_index_create: screen_c (%u) must be >= c (%u); the stage-1 screen can only "
			                      "be sparser than the dense sketch",
			                      data->two_stage_params.screen_c, dense_c);
		}
	} else if (has_two_stage_knobs) {
		throw BinderException("sylph_index_create: screen_c / min_sparse_kmers only apply to a two-stage output "
		                      "(two_stage := true or a .syl2db output_path)");
	}

	// Fail-fast schema validation (before the build side effect). read_id +
	// sequence1 are required; read_id may be VARCHAR or BIGINT.
	ValidateSequenceTableSchema(context, data->source_table, /*allow_bigint=*/true);

	// Validate the genome_id and order_by columns exist. A `LIMIT 0` probe works
	// for both tables and views and reports missing columns as a clean error.
	{
		auto conn = MakeReadOnlyHelperConnection(context);
		auto src = KeywordHelper::WriteOptionallyQuoted(data->source_table);
		auto gcol = KeywordHelper::WriteOptionallyQuoted(data->genome_id_col);
		auto ocol = KeywordHelper::WriteOptionallyQuoted(data->order_by_col);
		auto probe = conn.Query("SELECT " + gcol + ", " + ocol + " FROM " + src + " LIMIT 0");
		if (probe->HasError()) {
			throw BinderException("sylph_index_create: genome_id/order_by column check failed: %s", probe->GetError());
		}
		// genome_id is an identifier column — hold it to the same type contract
		// as read_id / the other tools (VARCHAR, BIGINT, or UUID). Its value is
		// persisted as the genome's file_name; the decimal/canonical string form
		// (what CAST/ToString produces) matches the id_column codec.
		if (!IsAllowedIdType(probe->GetTypes()[0])) {
			throw BinderException("sylph_index_create: genome_id column '%s' must be %s (got %s)",
			                      data->genome_id_col.c_str(), AllowedIdTypeList(),
			                      probe->GetTypes()[0].ToString().c_str());
		}
		// Detect an optional `comment` column (present in read_fastx output). Its
		// presence switches on full-header contig-name reconstruction below.
		auto comment_probe = conn.Query("SELECT comment FROM " + src + " LIMIT 0");
		data->has_comment = !comment_probe->HasError();
	}

	names = StringsToIdentifiers(data->names);
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

	auto src = KeywordHelper::WriteOptionallyQuoted(data.source_table);
	auto gcol = KeywordHelper::WriteOptionallyQuoted(data.genome_id_col);
	auto ocol = KeywordHelper::WriteOptionallyQuoted(data.order_by_col);

	// Contig name = the full FASTA header. sylph (via needletail) stores the whole
	// header line as first_contig_name; read_fastx splits it into read_id (first
	// token) + comment (remainder). When a comment column is present we rejoin
	// them so a miint-built .syldb's contig names match `sylph sketch`.
	const std::string contig_expr = data.has_comment
	                                    ? "CAST(read_id AS VARCHAR) || CASE WHEN comment IS NOT NULL AND comment <> '' "
	                                      "THEN ' ' || CAST(comment AS VARCHAR) ELSE '' END"
	                                    : "CAST(read_id AS VARCHAR)";

	// Enumerate the distinct genome ids up front (one lightweight scan of just the
	// id column). Each becomes an independent unit of work below.
	std::vector<Value> ids;
	{
		auto conn = MakeReadOnlyHelperConnection(context);
		auto res = conn.Query("SELECT DISTINCT " + gcol + " AS gk FROM " + src + " WHERE " + gcol +
		                      " IS NOT NULL ORDER BY gk");
		if (res->HasError()) {
			throw InvalidInputException("sylph_index_create: failed to enumerate genomes: %s", res->GetError());
		}
		// v2.0: Collection().GetValue rebuilds every row per call; build the rows once.
		auto rows = res->Collection().GetRows();
		for (idx_t r = 0; r < rows.size(); r++) {
			auto v = rows.GetValue(0, r);
			if (!v.IsNull()) {
				ids.push_back(std::move(v));
			}
		}
	}
	if (ids.empty()) {
		throw InvalidInputException("sylph_index_create: no non-NULL genome_id values found in '%s'",
		                            data.source_table.c_str());
	}

	// One builder per worker thread (mirrors `sylph sketch -t <cores>`; `threads`
	// overrides). Merged into the first before the single write. Freed on any throw.
	const idx_t db_threads = static_cast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	idx_t nthreads = data.user_threads != 0 ? static_cast<idx_t>(data.user_threads) : db_threads;
	nthreads = std::max<idx_t>(1, std::min<idx_t>(nthreads, ids.size()));

	std::vector<::SylphIndexBuilder *> builders(nthreads, nullptr);
	struct BuildersGuard {
		std::vector<::SylphIndexBuilder *> *v;
		~BuildersGuard() {
			if (v) {
				for (auto *b : *v) {
					if (b) {
						sylph_index_builder_free(b);
					}
				}
			}
		}
	} builders_guard {&builders};
	for (auto *&b : builders) {
		b = sylph_index_builder_create(&data.sketch_params);
		if (b == nullptr) {
			ThrowFFI("index builder create failed");
		}
	}

	std::atomic<size_t> next_genome {0};
	std::atomic<bool> failed {false};
	std::mutex err_mutex;
	std::string err_msg;

	// Each worker claims genome ids off the counter and sketches one at a time: a
	// `WHERE genome_id = <id>` query reads just that genome. Within it, consecutive
	// rows sharing the contig name are one contig (so a read_fastx whole-contig row
	// and Qiita's 64 KB chunk rows both work); first_contig_name is the lowest-
	// order_by contig.
	auto worker = [&](idx_t tid) {
		try {
			::SylphIndexBuilder *b = builders[tid];
			auto conn = MakeReadOnlyHelperConnection(context);
			size_t i;
			while (!failed.load(std::memory_order_relaxed) &&
			       (i = next_genome.fetch_add(1, std::memory_order_relaxed)) < ids.size()) {
				const std::string name = ids[i].ToString();   // genome id / file_name
				const std::string lit = ids[i].ToSQLString(); // typed literal → prune-friendly
				auto res = conn.Query("SELECT CAST(" + ocol + " AS BIGINT) AS ord, " + contig_expr +
				                      " AS nm, sequence1 AS seq FROM " + src + " WHERE " + gcol + " = " + lit +
				                      " AND sequence1 IS NOT NULL ORDER BY nm, ord");
				if (res->HasError()) {
					throw InvalidInputException("sylph_index_create: failed to read genome %s: %s", name.c_str(),
					                            res->GetError());
				}

				// Reassemble contigs (consecutive rows sharing nm) and sketch each.
				std::string cur_name, cur_seq;
				int64_t cur_order = 0;
				bool have_contig = false;
				auto flush = [&]() {
					if (sylph_index_builder_add_contig(b, name.c_str(), cur_order, cur_name.c_str(),
					                                   reinterpret_cast<const unsigned char *>(cur_seq.data()),
					                                   cur_seq.size()) != 0) {
						ThrowFFI("add_contig failed");
					}
				};
				auto rows = res->Collection().GetRows();
				for (idx_t r = 0; r < rows.size(); r++) {
					auto ordv = rows.GetValue(0, r);
					int64_t ord = ordv.IsNull() ? std::numeric_limits<int64_t>::max() : ordv.GetValue<int64_t>();
					std::string nm = rows.GetValue(1, r).ToString();
					std::string seq = rows.GetValue(2, r).ToString();
					if (!have_contig || nm != cur_name) {
						if (have_contig) {
							flush();
						}
						cur_name = std::move(nm);
						cur_seq = std::move(seq);
						cur_order = ord;
						have_contig = true;
					} else {
						cur_seq += seq;
						if (ord < cur_order) {
							cur_order = ord;
						}
					}
				}
				if (have_contig) {
					flush();
				}
				if (sylph_index_builder_end_genome(b, name.c_str()) != 0) {
					ThrowFFI("end_genome failed");
				}
			}
		} catch (const std::exception &e) {
			failed.store(true, std::memory_order_relaxed);
			std::lock_guard<std::mutex> lk(err_mutex);
			if (err_msg.empty()) {
				err_msg = e.what();
			}
		}
	};

	std::vector<std::thread> pool;
	pool.reserve(nthreads);
	for (idx_t t = 0; t < nthreads; t++) {
		pool.emplace_back(worker, t);
	}
	for (auto &th : pool) {
		th.join();
	}
	if (failed.load()) {
		throw IOException("sylph_index_create: %s", err_msg.c_str());
	}

	// Merge all per-thread builders into the first, then write once.
	for (idx_t t = 1; t < nthreads; t++) {
		if (sylph_index_builder_merge(builders[0], builders[t]) != 0) {
			ThrowFFI("merge failed");
		}
	}
	gstate->num_genomes = sylph_index_builder_num_genomes(builders[0]);
	if (data.two_stage) {
		if (sylph_index_builder_write_two_stage(builders[0], data.output_path.c_str(), &data.two_stage_params) != 0) {
			ThrowFFI("failed to write two-stage database");
		}
	} else if (sylph_index_builder_write(builders[0], data.output_path.c_str()) != 0) {
		ThrowFFI("failed to write database");
	}
	// builders_guard frees all builders on return.
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
		output.SetChildCardinality(0);
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

	output.SetChildCardinality(1);
	gstate.done = true;
}

// =============================================================================
// Registration
// =============================================================================
TableFunction SylphIndexCreateTableFunction::GetFunction() {
	TableFunction tf("sylph_index_create", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal,
	                 InitLocal);

	AddNamedParameter(tf, "genome_id", LogicalType::VARCHAR);
	AddNamedParameter(tf, "order_by", LogicalType::VARCHAR);
	AddNamedParameter(tf, "k", LogicalType::INTEGER);
	AddNamedParameter(tf, "c", LogicalType::INTEGER);
	AddNamedParameter(tf, "min_spacing", LogicalType::INTEGER);
	AddNamedParameter(tf, "pseudotax", LogicalType::BOOLEAN);
	AddNamedParameter(tf, "threads", LogicalType::INTEGER);
	AddNamedParameter(tf, "two_stage", LogicalType::BOOLEAN);
	AddNamedParameter(tf, "screen_c", LogicalType::INTEGER);
	AddNamedParameter(tf, "min_sparse_kmers", LogicalType::INTEGER);

	return tf;
}

void SylphIndexCreateTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb

#endif // MIINT_HAS_SYLPH
