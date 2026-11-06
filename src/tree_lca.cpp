#include "tree_lca.hpp"

#include "TreeIndex.hpp"
#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "miint_named_parameter.hpp"
#include "tree_function_utils.hpp"
#include "tree_table_reader.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector/flat_vector.hpp"

#include <algorithm>
#include <atomic>
#include <functional>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace duckdb {

namespace {

constexpr const char *kFn = "tree_lca";

struct BindData : public TableFunctionData {
	std::string tree_table;
	std::string pairs_table;
	LogicalType name_type;
	bool ignore_missing = false;
};

// A name and the tip it names, as held in the global state's TipMatches: each name is
// stored once, however many pairs it is in.
using TipEntry = TipMatches::value_type;

struct Pair {
	const TipEntry *a;
	const TipEntry *b;
};

struct GlobalState : public GlobalTableFunctionState {
	// Node-based, so the entries the pairs point at never move.
	TipMatches tips;
	IndexedTree tree;
	std::vector<Pair> pairs;
	std::atomic<idx_t> next_pair {0};
	idx_t max_threads = 1;

	idx_t MaxThreads() const override {
		return max_threads;
	}
};

unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input, vector<LogicalType> &return_types,
                              vector<Identifier> &names) {
	auto data = make_uniq<BindData>();
	data->tree_table = input.inputs[0].ToString();
	data->pairs_table = input.inputs[1].ToString();
	for (const auto *relation : {&data->tree_table, &data->pairs_table}) {
		RejectCTERelationName(input, *relation);
	}

	// Schemas only: no rows are read until InitGlobal (docs/internals/architecture.md).
	ValidateTreeTableSchema(context, data->tree_table);
	auto tree_columns = GetTableOrViewColumns(context, data->tree_table, "Tree table");
	auto pair_columns = GetTableOrViewColumns(context, data->pairs_table, "Pairs table");

	const LogicalType tree_name = TipNameType(kFn, data->tree_table, tree_columns);
	RejectTreeAs(kFn, data->tree_table, data->pairs_table, "pass the pairs as a relation of their own");
	LogicalType a_type, b_type;
	if (!FindColumnType(pair_columns, "a", a_type)) {
		throw BinderException("%s: pairs table '%s' has no 'a' column", kFn, data->pairs_table);
	}
	if (!FindColumnType(pair_columns, "b", b_type)) {
		throw BinderException("%s: pairs table '%s' has no 'b' column", kFn, data->pairs_table);
	}
	if (a_type != tree_name || b_type != tree_name) {
		throw BinderException("%s: tip names must have one type in the tree and the pairs, but 'name' is %s in '%s', "
		                      "and 'a' and 'b' are %s and %s in '%s'; cast them to one of %s",
		                      kFn, tree_name.ToString(), data->tree_table, a_type.ToString(), b_type.ToString(),
		                      data->pairs_table, AllowedIdTypeList());
	}
	data->name_type = tree_name;
	data->ignore_missing = GetBoolParam(input, "ignore_missing", false);

	names.emplace_back("a");
	return_types.push_back(data->name_type);
	names.emplace_back("b");
	return_types.push_back(data->name_type);
	names.emplace_back("lca_index");
	return_types.push_back(LogicalType::BIGINT);
	names.emplace_back("distance");
	return_types.push_back(LogicalType::DOUBLE);
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind = input.bind_data->Cast<BindData>();
	auto state = make_uniq<GlobalState>();

	// The pairs first, once: their names are the only ones kept from the tree. A pair with a
	// NULL name is dropped, and a repeated pair kept once. The set that finds repeats is
	// released before the tree is read.
	{
		struct PairHash {
			size_t operator()(const std::pair<const TipEntry *, const TipEntry *> &p) const {
				return std::hash<const void *>()(p.first) ^
				       (std::hash<const void *>()(p.second) * 0x9e3779b97f4a7c15ULL);
			}
		};
		std::unordered_set<std::pair<const TipEntry *, const TipEntry *>, PairHash> seen;
		std::vector<std::string> a_names, b_names;
		std::vector<bool> a_nulls, b_nulls;
		ReadRelation(context, kFn, bind.pairs_table, {{"a", bind.name_type}, {"b", bind.name_type}},
		             [&](DataChunk &chunk) {
			             ExtractIdColumnAsStrings(chunk, 0, bind.name_type, a_names, a_nulls);
			             ExtractIdColumnAsStrings(chunk, 1, bind.name_type, b_names, b_nulls);
			             for (idx_t i = 0; i < chunk.size(); ++i) {
				             if (a_nulls[i] || b_nulls[i]) {
					             continue;
				             }
				             const TipEntry *a = &*state->tips.emplace(std::move(a_names[i]), TipMatch {}).first;
				             const TipEntry *b = &*state->tips.emplace(std::move(b_names[i]), TipMatch {}).first;
				             if (seen.emplace(a, b).second) {
					             state->pairs.push_back({a, b});
				             }
			             }
		             });
	}

	state->tree = ReadIndexedTree(context, kFn, bind.tree_table, state->tips, bind.name_type,
	                              miint::LengthPolicy::UnknownAllowed, "read");
	RequireUnambiguousTips(kFn, bind.tree_table, state->tips, "a name in a pair");
	if (!bind.ignore_missing) {
		std::vector<std::string> missing;
		for (const auto &tip : state->tips) {
			if (tip.second.count == 0) {
				missing.push_back(tip.first);
			}
		}
		if (!missing.empty()) {
			std::sort(missing.begin(), missing.end());
			throw InvalidInputException("%s: %d name(s) in the pairs are not tips in the tree '%s': %s; pass "
			                            "ignore_missing := true to drop their pairs",
			                            kFn, static_cast<int64_t>(missing.size()), bind.tree_table, QuoteSome(missing));
		}
	}
	// Under ignore_missing, a pair with a name that is no tip has no row.
	auto &pairs = state->pairs;
	pairs.erase(std::remove_if(pairs.begin(), pairs.end(),
	                           [](const Pair &pair) { return pair.a->second.count != 1 || pair.b->second.count != 1; }),
	            pairs.end());

	state->max_threads = ThreadsForClaims(context, state->pairs.size(), STANDARD_VECTOR_SIZE);
	return std::move(state);
}

// Each call answers the next block of pairs; threads claim blocks independently, so the
// rows come out in no particular order.
void Execute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind = data_p.bind_data->Cast<BindData>();
	auto &state = data_p.global_state->Cast<GlobalState>();
	const auto &index = state.tree.index;
	const idx_t start = state.next_pair.fetch_add(STANDARD_VECTOR_SIZE);
	if (start >= state.pairs.size()) {
		output.SetChildCardinality(0);
		return;
	}
	const idx_t count = std::min<idx_t>(STANDARD_VECTOR_SIZE, state.pairs.size() - start);

	auto &a_vec = output.data[0];
	auto &b_vec = output.data[1];
	auto lca = FlatVector::GetDataMutable<int64_t>(output.data[2]);
	auto distance = FlatVector::GetDataMutable<double>(output.data[3]);
	auto &distance_valid = FlatVector::ValidityMutable(output.data[3]);
	for (idx_t i = 0; i < count; ++i) {
		const auto &pair = state.pairs[start + i];
		EmitIdCell(a_vec, i, pair.a->first, bind.name_type);
		EmitIdCell(b_vec, i, pair.b->first, bind.name_type);
		const uint32_t a = pair.a->second.node;
		const uint32_t b = pair.b->second.node;
		const uint32_t l = index.lca(a, b);
		lca[i] = state.tree.node_ids[l];
		try {
			SetDistance(distance, distance_valid, i, index.path_length(a, b, l));
		} catch (const std::overflow_error &) {
			throw InvalidInputException("%s: the path between '%s' and '%s' is outside DOUBLE's range", kFn,
			                            pair.a->first, pair.b->first);
		}
	}
	output.SetChildCardinality(count);
}

} // namespace

TableFunction TreeLcaTableFunction::GetFunction() {
	TableFunction tf("tree_lca", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal);
	AddNamedParameter(tf, "ignore_missing", LogicalType::BOOLEAN);
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void TreeLcaTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
