#include "tree_nearest_tip.hpp"

#include "NearestTip.hpp"
#include "TreeIndex.hpp"
#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "miint_named_parameter.hpp"
#include "tree_table_reader.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/value_map.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/parallel/task_scheduler.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <algorithm>
#include <atomic>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace duckdb {

namespace {

constexpr const char *kFn = "tree_nearest_tip";
constexpr uint32_t kNoGroup = UINT32_MAX;

struct BindData : public TableFunctionData {
	std::string tree_table;
	std::string targets_table;
	std::string queries_table;
	miint::NearestTipOptions options;
	bool ignore_missing = false;
	LogicalType name_type;
	bool has_groups = false;
	LogicalType grp_type;
};

// A query to answer: its name, its group as given (for the output), its tip, and the
// group of targets it competes in (kNoGroup when it has none).
struct QueryEntry {
	std::string name;
	Value grp;
	uint32_t node;
	uint32_t group;
};

struct GlobalState : public GlobalTableFunctionState {
	miint::TreeIndex index;
	// The caller's node_index for each node of the index, for lca_index.
	std::vector<int64_t> node_ids;
	// One target set per group (a single one without groups); each points at `index`.
	std::vector<miint::NearestTipTargets> groups;
	std::unordered_map<uint32_t, std::string> target_names;
	std::vector<QueryEntry> queries;
	std::atomic<idx_t> next_query {0};
	idx_t max_threads = 1;

	idx_t MaxThreads() const override {
		return max_threads;
	}
};

// The type of `column` among a relation's columns, matched case-insensitively as DuckDB
// matches column names.
bool FindColumnType(const TableOrViewColumns &columns, const std::string &column, LogicalType &type) {
	for (idx_t i = 0; i < columns.names.size(); ++i) {
		if (StringUtil::CIEquals(columns.names[i], column)) {
			type = columns.types[i];
			return true;
		}
	}
	return false;
}

// Bind fixed the output types. A relation can change type without a catalog change (a
// view over a file that is rewritten, under a prepared statement), so a column read with
// another type is refused rather than emitted under the old one. The check uses the types
// of the read itself (docs/internals/reading-tables-views.md says why not the catalog's).
void RequireBoundType(const std::string &relation, const std::string &column, const LogicalType &now,
                      const LogicalType &bound) {
	if (now != bound) {
		throw InvalidInputException(
		    "%s: column '%s' of '%s' is %s now but was %s when the statement was bound; run or prepare it again", kFn,
		    column, relation, now.ToString(), bound.ToString());
	}
}

struct NamedRow {
	std::string name;
	Value grp;
};

// Read `name`, and `grp` with groups, from a targets or queries relation: once, at
// execution, over a helper connection (docs/internals/reading-tables-views.md). Rows
// with a NULL name are dropped.
std::vector<NamedRow> ReadNames(ClientContext &context, const std::string &relation, const BindData &bind) {
	auto conn = MakeReadOnlyHelperConnection(context);
	const std::string sql =
	    std::string("SELECT name") + (bind.has_groups ? ", grp" : "") + " FROM " + SQLIdentifier::ToString(relation);
	auto result = conn.Query(sql);
	if (result->HasError()) {
		throw InvalidInputException("%s: failed to read '%s': %s", kFn, relation, result->GetError());
	}
	auto &types = result->GetTypes();
	RequireBoundType(relation, "name", types[0], bind.name_type);
	if (bind.has_groups) {
		RequireBoundType(relation, "grp", types[1], bind.grp_type);
	}

	std::vector<NamedRow> rows;
	std::vector<std::string> names;
	std::vector<bool> nulls;
	while (auto chunk = result->Fetch()) {
		if (chunk->size() == 0) {
			break;
		}
		context.InterruptCheck();
		ExtractIdColumnAsStrings(*chunk, 0, bind.name_type, names, nulls);
		for (idx_t i = 0; i < chunk->size(); ++i) {
			if (!nulls[i]) {
				rows.push_back({std::move(names[i]), bind.has_groups ? chunk->GetValue(1, i) : Value()});
			}
		}
	}
	return rows;
}

// Sort ids so a tie goes to the smallest name in its type's order: VARCHAR by bytes,
// whatever its collation; BIGINT by value; UUID by its canonical text, which DuckDB orders
// the same way as the UUIDs themselves.
void SortInTypeOrder(std::vector<std::string> &ids, const LogicalType &type) {
	switch (type.id()) {
	case LogicalTypeId::VARCHAR:
	case LogicalTypeId::UUID:
		std::sort(ids.begin(), ids.end());
		return;
	case LogicalTypeId::BIGINT: {
		std::vector<std::pair<int64_t, std::string>> keyed;
		keyed.reserve(ids.size());
		for (auto &id : ids) {
			keyed.emplace_back(::miint::ParseIdAsInt64(id).value(), std::move(id));
		}
		std::sort(keyed.begin(), keyed.end());
		for (size_t i = 0; i < keyed.size(); ++i) {
			ids[i] = std::move(keyed[i].second);
		}
		return;
	}
	default:
		throw InternalException("%s: unsupported tip name type %s", kFn, type.ToString());
	}
}

// "'a', 'b', ... (n more)": the first ten names, quoted.
std::string QuoteSome(const std::vector<std::string> &names) {
	constexpr size_t cap = 10;
	std::string out;
	for (size_t i = 0; i < names.size() && i < cap; ++i) {
		out += (i ? ", '" : "'") + names[i] + "'";
	}
	if (names.size() > cap) {
		out += ", ... (" + std::to_string(names.size() - cap) + " more)";
	}
	return out;
}

unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input, vector<LogicalType> &return_types,
                              vector<Identifier> &names) {
	auto data = make_uniq<BindData>();
	data->tree_table = input.inputs[0].ToString();
	data->targets_table = input.inputs[1].ToString();
	data->queries_table = input.inputs[2].ToString();
	for (const auto *relation : {&data->tree_table, &data->targets_table, &data->queries_table}) {
		RejectCTERelationName(input, *relation);
	}

	// Schemas only: no rows are read until InitGlobal (docs/internals/architecture.md).
	ValidateTreeTableSchema(context, data->tree_table);
	auto tree_columns = GetTableOrViewColumns(context, data->tree_table, "Tree table");
	auto target_columns = GetTableOrViewColumns(context, data->targets_table, "Targets table");
	auto query_columns = GetTableOrViewColumns(context, data->queries_table, "Queries table");

	LogicalType tree_name, target_name, query_name, tree_length;
	if (!FindColumnType(tree_columns, "name", tree_name)) {
		throw BinderException("%s: tree table '%s' has no 'name' column", kFn, data->tree_table);
	}
	if (!FindColumnType(tree_columns, "branch_length", tree_length)) {
		throw BinderException("%s: tree table '%s' has no 'branch_length' column", kFn, data->tree_table);
	}
	// The tree is read in its own projection, so passing it as targets or queries would read
	// it twice; its internal labels are not tips either.
	for (const auto *names_table : {&data->targets_table, &data->queries_table}) {
		if (StringUtil::CIEquals(*names_table, data->tree_table)) {
			throw BinderException("%s: '%s' is the tree table; pass the tip names as a relation of their own, such "
			                      "as a view of the tree's tips",
			                      kFn, *names_table);
		}
	}
	if (!FindColumnType(target_columns, "name", target_name)) {
		throw BinderException("%s: targets table '%s' has no 'name' column", kFn, data->targets_table);
	}
	if (!FindColumnType(query_columns, "name", query_name)) {
		throw BinderException("%s: queries table '%s' has no 'name' column", kFn, data->queries_table);
	}
	if (target_name != tree_name || query_name != tree_name) {
		throw BinderException("%s: tip names must have the same type in all three tables, but 'name' is %s in '%s', "
		                      "%s in '%s' and %s in '%s'; cast them to one of %s",
		                      kFn, tree_name.ToString(), data->tree_table, target_name.ToString(), data->targets_table,
		                      query_name.ToString(), data->queries_table, AllowedIdTypeList());
	}
	if (!IsAllowedIdType(tree_name)) {
		throw BinderException("%s: tip names must be %s, but they are %s; cast them", kFn, AllowedIdTypeList(),
		                      tree_name.ToString());
	}
	data->name_type = tree_name;

	LogicalType target_grp, query_grp;
	const bool targets_grouped = FindColumnType(target_columns, "grp", target_grp);
	const bool queries_grouped = FindColumnType(query_columns, "grp", query_grp);
	if (targets_grouped != queries_grouped) {
		throw BinderException("%s: 'grp' must be in both the targets and the queries table, or in neither; it is only "
		                      "in '%s'",
		                      kFn, targets_grouped ? data->targets_table : data->queries_table);
	}
	if (targets_grouped && target_grp != query_grp) {
		throw BinderException("%s: 'grp' has different types in '%s' (%s) and '%s' (%s); cast them to one type", kFn,
		                      data->targets_table, target_grp.ToString(), data->queries_table, query_grp.ToString());
	}
	data->has_groups = targets_grouped;
	data->grp_type = target_grp;

	auto scope = input.named_parameters.find(Identifier("scope"));
	if (scope != input.named_parameters.end() && !scope->second.IsNull()) {
		const auto value = StringUtil::Lower(StringValue::Get(scope->second));
		if (value == "path") {
			data->options.scope = miint::NearestTipScope::Path;
		} else if (value == "clade") {
			data->options.scope = miint::NearestTipScope::Clade;
		} else {
			throw BinderException("%s: scope must be 'path' or 'clade', not '%s'", kFn,
			                      StringValue::Get(scope->second));
		}
	}
	data->options.exclude_self = GetBoolParam(input, "exclude_self", true);
	data->ignore_missing = GetBoolParam(input, "ignore_missing", false);

	names.emplace_back("query");
	return_types.push_back(data->name_type);
	if (data->has_groups) {
		names.emplace_back("grp");
		return_types.push_back(data->grp_type);
	}
	names.emplace_back("target");
	return_types.push_back(data->name_type);
	names.emplace_back("distance");
	return_types.push_back(LogicalType::DOUBLE);
	names.emplace_back("n_tied");
	return_types.push_back(LogicalType::BIGINT);
	names.emplace_back("lca_index");
	return_types.push_back(LogicalType::BIGINT);
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind = input.bind_data->Cast<BindData>();
	auto state = make_uniq<GlobalState>();

	// Each relation is read exactly once. Targets and queries are read in their own types,
	// which ReadNames checks against the bound output types. Tree names are only matched,
	// as text, and never emitted, so the tree's reader (ReadTreeArrays) reads them as text,
	// after checking their type against the bound one.
	// One relation as both targets and queries (identifiers are case-insensitive) is read
	// once for both, so both roles see the same rows.
	const auto targets = ReadNames(context, bind.targets_table, bind);
	const bool queries_are_targets = StringUtil::CIEquals(bind.targets_table, bind.queries_table);
	const auto queries_read =
	    queries_are_targets ? std::vector<NamedRow>() : ReadNames(context, bind.queries_table, bind);
	const auto &queries = queries_are_targets ? targets : queries_read;

	// The tips that queries and targets name: which node, and how many tips carry the name.
	struct TipMatch {
		uint32_t node = 0;
		uint32_t count = 0;
	};
	std::unordered_map<std::string, TipMatch> tips;
	for (const auto *rows : {&targets, &queries}) {
		for (const auto &row : *rows) {
			tips.emplace(row.name, TipMatch {});
		}
	}

	{
		// Only the names that queries and targets use are kept from the tree.
		std::unordered_set<std::string_view> keep_names;
		keep_names.reserve(tips.size());
		for (const auto &tip : tips) {
			keep_names.insert(tip.first);
		}
		TreeTableArrays tree;
		try {
			tree = ReadTreeArrays(context, bind.tree_table, keep_names, bind.name_type);
			state->index = miint::TreeIndex::Build(tree.tree.parent, tree.tree.branch_length, &tree.tree.node_ids);
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("%s: the tree '%s' cannot be searched: %s", kFn, bind.tree_table, e.what());
		}
		for (const auto &named : tree.named) {
			if (state->index.is_tip(named.first)) {
				auto &match = tips.at(named.second);
				match.node = named.first;
				++match.count;
			}
		}
		state->node_ids = std::move(tree.tree.node_ids);
		// The arrays are released here: the index and the names in use are all the search needs.
	}

	std::vector<std::string> ambiguous;
	for (const auto &tip : tips) {
		if (tip.second.count > 1) {
			ambiguous.push_back(tip.first);
		}
	}
	if (!ambiguous.empty()) {
		std::sort(ambiguous.begin(), ambiguous.end());
		throw InvalidInputException(
		    "%s: tip name '%s' names %d tips in the tree '%s'; a query or target must name exactly one tip", kFn,
		    ambiguous.front(), static_cast<int64_t>(tips[ambiguous.front()].count), bind.tree_table);
	}
	if (!bind.ignore_missing) {
		for (const auto &role : {std::make_pair("query", &queries), std::make_pair("target", &targets)}) {
			std::set<std::string> missing;
			for (const auto &row : *role.second) {
				if (tips[row.name].count == 0) {
					missing.insert(row.name);
				}
			}
			if (!missing.empty()) {
				throw InvalidInputException(
				    "%s: %d %s name(s) are not tips in the tree '%s': %s; pass ignore_missing := true to skip them",
				    kFn, static_cast<int64_t>(missing.size()), role.first, bind.tree_table,
				    QuoteSome(std::vector<std::string>(missing.begin(), missing.end())));
			}
		}
	}

	// Target names ranked in their type's order; a tie goes to the smallest rank.
	std::unordered_map<std::string, uint32_t> rank_of;
	{
		std::vector<std::string> ranked;
		for (const auto &row : targets) {
			if (tips[row.name].count == 1 && rank_of.emplace(row.name, 0).second) {
				ranked.push_back(row.name);
			}
		}
		SortInTypeOrder(ranked, bind.name_type);
		for (uint32_t r = 0; r < ranked.size(); ++r) {
			rank_of[ranked[r]] = r;
		}
	}

	// One target set per group. Without groups every target is in group 0; with them, a
	// target with a NULL group is in none.
	value_map_t<uint32_t> group_of;
	std::vector<std::vector<uint32_t>> group_nodes(bind.has_groups ? 0 : 1);
	std::vector<std::vector<uint32_t>> group_ranks(bind.has_groups ? 0 : 1);
	std::unordered_set<uint64_t> placed; // (group << 32) | node, so a repeated row counts once
	for (const auto &row : targets) {
		const auto &match = tips[row.name];
		if (match.count != 1) {
			continue; // missing, and skipped under ignore_missing
		}
		uint32_t group = 0;
		if (bind.has_groups) {
			if (row.grp.IsNull()) {
				continue;
			}
			auto found = group_of.find(row.grp);
			if (found == group_of.end()) {
				group = static_cast<uint32_t>(group_nodes.size());
				group_of.emplace(row.grp, group);
				group_nodes.emplace_back();
				group_ranks.emplace_back();
			} else {
				group = found->second;
			}
		}
		if (!placed.insert((static_cast<uint64_t>(group) << 32) | match.node).second) {
			continue;
		}
		group_nodes[group].push_back(match.node);
		group_ranks[group].push_back(rank_of.at(row.name));
		state->target_names.emplace(match.node, row.name);
	}
	state->groups.reserve(group_nodes.size());
	for (size_t g = 0; g < group_nodes.size(); ++g) {
		state->groups.emplace_back(state->index, group_nodes[g], group_ranks[g]);
	}

	// One entry per distinct (name, grp). Query groups that no target has still make
	// distinct rows, so they are told apart by value, not by target group.
	value_map_t<uint32_t> query_group_key;
	struct NameKeyHash {
		size_t operator()(const std::pair<std::string, uint32_t> &p) const {
			return std::hash<std::string>()(p.first) ^ (std::hash<uint32_t>()(p.second) * 0x9e3779b97f4a7c15ULL);
		}
	};
	std::unordered_set<std::pair<std::string, uint32_t>, NameKeyHash> seen;
	for (const auto &row : queries) {
		const auto &match = tips[row.name];
		if (match.count != 1) {
			continue; // missing, and dropped under ignore_missing
		}
		uint32_t key = 0;
		uint32_t group = 0;
		if (bind.has_groups) {
			key = kNoGroup;
			group = kNoGroup;
			if (!row.grp.IsNull()) {
				key = query_group_key.emplace(row.grp, static_cast<uint32_t>(query_group_key.size())).first->second;
				auto found = group_of.find(row.grp);
				if (found != group_of.end()) {
					group = found->second;
				}
			}
		}
		if (seen.emplace(row.name, key).second) {
			state->queries.push_back({row.name, row.grp, match.node, group});
		}
	}

	const idx_t blocks = (state->queries.size() + STANDARD_VECTOR_SIZE - 1) / STANDARD_VECTOR_SIZE;
	const idx_t threads = TaskScheduler::GetScheduler(context).NumberOfThreads();
	state->max_threads = std::max<idx_t>(1, std::min(blocks, threads));
	return std::move(state);
}

// Each call answers the next block of queries; threads claim blocks independently, so the
// rows come out in no particular order.
void Execute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind = data_p.bind_data->Cast<BindData>();
	auto &state = data_p.global_state->Cast<GlobalState>();
	const idx_t start = state.next_query.fetch_add(STANDARD_VECTOR_SIZE);
	if (start >= state.queries.size()) {
		output.SetChildCardinality(0);
		return;
	}
	const idx_t count = std::min<idx_t>(STANDARD_VECTOR_SIZE, state.queries.size() - start);

	idx_t col = 0;
	auto &query_vec = output.data[col++];
	Vector *grp_vec = bind.has_groups ? &output.data[col++] : nullptr;
	auto &target_vec = output.data[col++];
	auto &distance_vec = output.data[col++];
	auto &tied_vec = output.data[col++];
	auto &lca_vec = output.data[col++];
	auto distance = FlatVector::GetDataMutable<double>(distance_vec);
	auto tied = FlatVector::GetDataMutable<int64_t>(tied_vec);
	auto lca = FlatVector::GetDataMutable<int64_t>(lca_vec);
	auto &distance_valid = FlatVector::ValidityMutable(distance_vec);
	auto &lca_valid = FlatVector::ValidityMutable(lca_vec);

	for (idx_t i = 0; i < count; ++i) {
		const auto &entry = state.queries[start + i];
		EmitIdCell(query_vec, i, entry.name, bind.name_type);
		if (grp_vec) {
			grp_vec->SetValue(i, entry.grp);
		}
		miint::NearestTipResult found;
		if (entry.group != kNoGroup) {
			try {
				found = state.groups[entry.group].Search(entry.node, bind.options);
			} catch (const std::overflow_error &) {
				throw InvalidInputException("%s: the path from '%s' to its nearest target is outside DOUBLE's range",
				                            kFn, entry.name);
			}
		}
		tied[i] = static_cast<int64_t>(found.n_tied);
		if (found.target == miint::NearestTipResult::NO_TARGET) {
			FlatVector::ValidityMutable(target_vec).SetInvalid(i);
			distance_valid.SetInvalid(i);
			lca_valid.SetInvalid(i);
			continue;
		}
		EmitIdCell(target_vec, i, state.target_names.at(found.target), bind.name_type);
		distance[i] = found.distance;
		distance_valid.SetValid(i);
		lca[i] = state.node_ids[found.lca];
		lca_valid.SetValid(i);
	}
	output.SetChildCardinality(count);
}

} // namespace

TableFunction TreeNearestTipTableFunction::GetFunction() {
	TableFunction tf("tree_nearest_tip", {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute,
	                 Bind, InitGlobal);
	AddNamedParameter(tf, "scope", LogicalType::VARCHAR);
	AddNamedParameter(tf, "exclude_self", LogicalType::BOOLEAN);
	AddNamedParameter(tf, "ignore_missing", LogicalType::BOOLEAN);
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void TreeNearestTipTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
