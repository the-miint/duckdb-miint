#include "match_short_barcodes.hpp"
#include "catalog_utils.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/parser/keyword_helper.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/main/materialized_query_result.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

namespace duckdb {

// ---------------------------------------------------------------------------
// 2-bit packing
// ---------------------------------------------------------------------------
// A=00, C=01, G=10, T=11. 32 bases per uint64_t. Padding bits are 0; both
// query and ref share the same length and same all-zero padding, so XOR of
// padding bits is always 0 and does not contribute to Hamming distance.

struct PackedBarcode {
	std::string id;
	std::vector<uint64_t> words;
};

static inline uint8_t Encode(char c, const char *fn_name) {
	switch (c) {
	case 'A':
	case 'a':
		return 0;
	case 'C':
	case 'c':
		return 1;
	case 'G':
	case 'g':
		return 2;
	case 'T':
	case 't':
		return 3;
	default:
		throw InvalidInputException("%s: invalid base '%c' (only ACGT allowed)", fn_name, c);
	}
}

static void Pack(const std::string &seq, std::vector<uint64_t> &out, const char *fn_name) {
	const size_t n_words = (seq.size() + 31) / 32;
	out.assign(n_words, 0);
	for (size_t i = 0; i < seq.size(); ++i) {
		uint8_t code = Encode(seq[i], fn_name);
		size_t word = i / 32;
		size_t shift = (i % 32) * 2;
		out[word] |= static_cast<uint64_t>(code) << shift;
	}
}

// Number of 2-bit positions where a and b differ. Returns as soon as nm
// reaches `early_exit_at` so multi-word barcodes that blow the budget on the
// first word are rejected before scanning the rest. The returned value is
// only guaranteed accurate when ≤ early_exit_at; the caller rejects anything
// above the budget unconditionally.
static inline int32_t Hamming(const uint64_t *a, const uint64_t *b, size_t n_words, int32_t early_exit_at) {
	int32_t nm = 0;
	for (size_t i = 0; i < n_words; ++i) {
		uint64_t x = a[i] ^ b[i];
		// Collapse each 2-bit position to a single "differs?" bit
		uint64_t differs = (x | (x >> 1)) & 0x5555555555555555ULL;
		nm += __builtin_popcountll(differs);
		if (nm >= early_exit_at) {
			return nm;
		}
	}
	return nm;
}

// ---------------------------------------------------------------------------
// Match hit
// ---------------------------------------------------------------------------
struct MatchHit {
	uint32_t query_idx;
	uint32_t ref_idx;
	int32_t nm;
};

// ---------------------------------------------------------------------------
// Table-function state
// ---------------------------------------------------------------------------
struct MatchData : public TableFunctionData {
	std::string query_table;
	std::string ref_table;
	int32_t max_nm = 0;
	bool report_all = true;
};

struct MatchGlobalState : public GlobalTableFunctionState {
	std::vector<PackedBarcode> queries;
	std::vector<PackedBarcode> refs;
	std::vector<MatchHit> hits;
	idx_t hit_offset = 0;

	idx_t MaxThreads() const override {
		return 1;
	}
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
static constexpr const char *FN_NAME = "match_short_barcodes";

// Read (id, sequence) pairs from the named table via a separate connection
// (the `docs/internals/reading-tables-views.md` recipe — avoids context-locked
// deadlocks). Validates all sequences have the same length and packs them.
// Returns the per-barcode length; throws on mixed lengths or non-ACGT bases.
static size_t LoadAndPack(ClientContext &context, const std::string &table_name, const char *role,
                          std::vector<PackedBarcode> &out) {
	auto conn = MakeReadOnlyHelperConnection(context);
	std::string query = "SELECT id, sequence FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);
	auto result = conn.Query(query);
	if (result->HasError()) {
		throw InvalidInputException("%s: failed to read %s table '%s': %s", FN_NAME, role, table_name,
		                            result->GetError());
	}

	auto &materialized = result->Cast<MaterializedQueryResult>();
	size_t expected_len = 0;
	bool first = true;
	while (auto chunk = materialized.Fetch()) {
		auto &id_vec = chunk->data[0];
		auto &seq_vec = chunk->data[1];
		UnifiedVectorFormat id_data, seq_data;
		id_vec.ToUnifiedFormat(chunk->size(), id_data);
		seq_vec.ToUnifiedFormat(chunk->size(), seq_data);
		auto id_ptr = UnifiedVectorFormat::GetData<string_t>(id_data);
		auto seq_ptr = UnifiedVectorFormat::GetData<string_t>(seq_data);

		for (idx_t i = 0; i < chunk->size(); ++i) {
			auto idi = id_data.sel->get_index(i);
			auto seqi = seq_data.sel->get_index(i);
			if (!id_data.validity.RowIsValid(idi) || !seq_data.validity.RowIsValid(seqi)) {
				throw InvalidInputException("%s: NULL id or sequence in %s table '%s'", FN_NAME, role, table_name);
			}
			std::string id(id_ptr[idi].GetData(), id_ptr[idi].GetSize());
			std::string seq(seq_ptr[seqi].GetData(), seq_ptr[seqi].GetSize());
			if (seq.empty()) {
				throw InvalidInputException("%s: empty sequence in %s table '%s' (id=%s)", FN_NAME, role, table_name,
				                            id);
			}
			if (first) {
				expected_len = seq.size();
				first = false;
			} else if (seq.size() != expected_len) {
				throw InvalidInputException("%s: %s table '%s' has inconsistent barcode length (expected %zu, got %zu)",
				                            FN_NAME, role, table_name, expected_len, seq.size());
			}
			PackedBarcode pb;
			pb.id = std::move(id);
			Pack(seq, pb.words, FN_NAME);
			out.push_back(std::move(pb));
		}
	}
	return expected_len;
}

// Verify a table exists (separately from row reading) — gives a clean error
// before InitGlobal materialization. Both tables and views work via the
// catalog TABLE_ENTRY lookup.
static void ValidateTableExists(ClientContext &context, const std::string &table_name) {
	auto conn = MakeReadOnlyHelperConnection(context);
	std::string query = "SELECT id, sequence FROM " + KeywordHelper::WriteOptionallyQuoted(table_name) + " LIMIT 0";
	auto result = conn.Query(query);
	if (result->HasError()) {
		throw BinderException("%s: '%s' does not exist or is missing 'id'/'sequence' columns (%s)", FN_NAME, table_name,
		                      result->GetError());
	}
}

// ---------------------------------------------------------------------------
// Bind
// ---------------------------------------------------------------------------
static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
                                     vector<LogicalType> &return_types, vector<std::string> &names) {
	auto data = make_uniq<MatchData>();
	data->query_table = input.inputs[0].GetValue<std::string>();
	data->ref_table = input.inputs[1].GetValue<std::string>();

	ValidateTableExists(context, data->query_table);
	ValidateTableExists(context, data->ref_table);

	auto max_nm_it = input.named_parameters.find("max_nm");
	if (max_nm_it == input.named_parameters.end()) {
		throw BinderException("%s requires 'max_nm' parameter (max Hamming distance, e.g. max_nm:=2)", FN_NAME);
	}
	data->max_nm = max_nm_it->second.GetValue<int32_t>();
	if (data->max_nm < 0) {
		throw InvalidInputException("%s: max_nm must be >= 0 (got %d)", FN_NAME, data->max_nm);
	}

	auto report_it = input.named_parameters.find("report_all");
	if (report_it != input.named_parameters.end()) {
		data->report_all = report_it->second.GetValue<bool>();
	}

	names = {"query_id", "ref_id", "nm"};
	return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::INTEGER};

	return data;
}

// ---------------------------------------------------------------------------
// InitGlobal — load and pack both tables, compute all hits
// ---------------------------------------------------------------------------
static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<MatchData>();
	auto gstate = make_uniq<MatchGlobalState>();

	const size_t qlen = LoadAndPack(context, data.query_table, "query", gstate->queries);
	const size_t rlen = LoadAndPack(context, data.ref_table, "ref", gstate->refs);

	// Both tables may be empty; only fail if lengths differ when both are
	// non-empty.
	if (!gstate->queries.empty() && !gstate->refs.empty() && qlen != rlen) {
		throw InvalidInputException("%s: query barcode length (%zu) does not match ref barcode length (%zu)", FN_NAME,
		                            qlen, rlen);
	}

	if (gstate->queries.empty() || gstate->refs.empty()) {
		return gstate;
	}

	const size_t n_words = gstate->queries[0].words.size();
	const int32_t max_nm = data.max_nm;

	// Brute Q × R Hamming. Reference set is small (thousands), brute force is
	// fine and SIMD-friendly via popcount intrinsic.
	for (uint32_t qi = 0; qi < gstate->queries.size(); ++qi) {
		const uint64_t *qa = gstate->queries[qi].words.data();
		for (uint32_t ri = 0; ri < gstate->refs.size(); ++ri) {
			const uint64_t *ra = gstate->refs[ri].words.data();
			int32_t nm = Hamming(qa, ra, n_words, max_nm + 1);
			if (nm <= max_nm) {
				gstate->hits.push_back(MatchHit {qi, ri, nm});
			}
		}
	}

	// report_all := false → reduce to one best hit per query, tie-break by
	// ref_id ASC.
	if (!data.report_all) {
		// Sort by (query_idx ASC, nm ASC, ref_id ASC) — the first hit per
		// query group is therefore the best (lowest nm, lex-min ref_id on tie).
		std::sort(gstate->hits.begin(), gstate->hits.end(), [&](const MatchHit &a, const MatchHit &b) {
			if (a.query_idx != b.query_idx) {
				return a.query_idx < b.query_idx;
			}
			if (a.nm != b.nm) {
				return a.nm < b.nm;
			}
			return gstate->refs[a.ref_idx].id < gstate->refs[b.ref_idx].id;
		});
		std::vector<MatchHit> best;
		best.reserve(gstate->queries.size());
		uint32_t cur_q = std::numeric_limits<uint32_t>::max();
		for (auto &h : gstate->hits) {
			if (h.query_idx != cur_q) {
				cur_q = h.query_idx;
				best.push_back(h);
			}
		}
		gstate->hits = std::move(best);
	}

	return gstate;
}

// ---------------------------------------------------------------------------
// Execute — paginate hits into output chunks
// ---------------------------------------------------------------------------
static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	(void)context;
	auto &gstate = data_p.global_state->Cast<MatchGlobalState>();

	if (gstate.hit_offset >= gstate.hits.size()) {
		output.SetCardinality(0);
		return;
	}

	idx_t remaining = gstate.hits.size() - gstate.hit_offset;
	idx_t count = std::min(remaining, static_cast<idx_t>(STANDARD_VECTOR_SIZE));

	auto &query_id_vec = output.data[0];
	auto &ref_id_vec = output.data[1];
	auto nm_data = FlatVector::GetData<int32_t>(output.data[2]);

	for (idx_t i = 0; i < count; ++i) {
		const auto &hit = gstate.hits[gstate.hit_offset + i];
		FlatVector::GetData<string_t>(query_id_vec)[i] =
		    StringVector::AddString(query_id_vec, gstate.queries[hit.query_idx].id);
		FlatVector::GetData<string_t>(ref_id_vec)[i] = StringVector::AddString(ref_id_vec, gstate.refs[hit.ref_idx].id);
		nm_data[i] = hit.nm;
	}

	output.SetCardinality(count);
	gstate.hit_offset += count;
}

// ---------------------------------------------------------------------------
// Register
// ---------------------------------------------------------------------------
void MatchShortBarcodesTableFunction::Register(ExtensionLoader &loader) {
	TableFunction tf("match_short_barcodes", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal);
	tf.named_parameters["max_nm"] = LogicalType::INTEGER;
	tf.named_parameters["report_all"] = LogicalType::BOOLEAN;
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(tf);
}

} // namespace duckdb
