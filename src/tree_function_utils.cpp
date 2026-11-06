#include "tree_function_utils.hpp"

#include "id_column_utils.hpp"
#include "miint_streaming_query.hpp"
#include "tree_table_reader.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/exception/binder_exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/parallel/task_scheduler.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <algorithm>
#include <stdexcept>
#include <string_view>
#include <unordered_set>

namespace duckdb {

bool FindColumnType(const TableOrViewColumns &columns, const std::string &column, LogicalType &type) {
	for (idx_t i = 0; i < columns.names.size(); ++i) {
		if (StringUtil::CIEquals(columns.names[i], column)) {
			type = columns.types[i];
			return true;
		}
	}
	return false;
}

LogicalType TipNameType(const char *fn, const std::string &tree_table, const TableOrViewColumns &tree_columns) {
	LogicalType type;
	if (!FindColumnType(tree_columns, "name", type)) {
		throw BinderException("%s: tree table '%s' has no 'name' column", fn, tree_table);
	}
	if (!IsAllowedIdType(type)) {
		throw BinderException("%s: tip names must be %s, but they are %s; cast them", fn, AllowedIdTypeList(),
		                      type.ToString());
	}
	return type;
}

void RejectTreeAs(const char *fn, const std::string &tree_table, const std::string &relation, const char *remedy) {
	if (StringUtil::CIEquals(relation, tree_table)) {
		throw BinderException("%s: '%s' is the tree table; %s", fn, relation, remedy);
	}
}

void ReadRelation(ClientContext &context, const char *fn, const std::string &relation,
                  const std::vector<std::pair<std::string, LogicalType>> &columns,
                  const std::function<void(DataChunk &)> &consume) {
	std::string sql = "SELECT ";
	for (size_t i = 0; i < columns.size(); ++i) {
		sql += (i ? ", " : "") + columns[i].first;
	}
	sql += " FROM " + SQLIdentifier::ToString(relation);

	auto conn = MakeReadOnlyHelperConnection(context);
	auto stream = SubmitStream(conn, sql);
	if (stream->HasError()) {
		throw InvalidInputException("%s: failed to read '%s': %s", fn, relation, stream->GetError());
	}
	auto &types = stream->GetTypes();
	for (size_t i = 0; i < columns.size(); ++i) {
		if (types[i] != columns[i].second) {
			throw InvalidInputException(
			    "%s: column '%s' of '%s' is %s now but was %s when the statement was bound; run or prepare it again",
			    fn, columns[i].first, relation, types[i].ToString(), columns[i].second.ToString());
		}
	}
	while (auto chunk = stream->Fetch()) {
		if (chunk->size() == 0) {
			break;
		}
		context.InterruptCheck();
		consume(*chunk);
	}
	// A read that fails part-way ends the stream like a finished one; only the error tells.
	if (stream->HasError()) {
		throw InvalidInputException("%s: failed to read '%s': %s", fn, relation, stream->GetError());
	}
}

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

IndexedTree ReadIndexedTree(ClientContext &context, const char *fn, const std::string &tree_table, TipMatches &tips,
                            const LogicalType &name_type, miint::LengthPolicy policy, const char *use) {
	std::unordered_set<std::string_view> keep_names;
	keep_names.reserve(tips.size());
	for (const auto &tip : tips) {
		keep_names.insert(tip.first);
	}
	IndexedTree out;
	TreeTableArrays tree;
	try {
		tree = ReadTreeArrays(context, tree_table, keep_names, name_type);
		out.index = miint::TreeIndex::Build(tree.tree.parent, tree.tree.branch_length, &tree.tree.node_ids, policy);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s: the tree '%s' cannot be %s: %s", fn, tree_table, use, e.what());
	}
	for (const auto &named : tree.named) {
		if (out.index.is_tip(named.first)) {
			auto &match = tips.at(named.second);
			match.node = named.first;
			++match.count;
		}
	}
	out.node_ids = std::move(tree.tree.node_ids);
	// The rest of the arrays are released here: the index, the ids and the matches are all
	// a caller needs.
	return out;
}

idx_t ThreadsForClaims(ClientContext &context, idx_t items, idx_t per_claim) {
	const idx_t claims = (items + per_claim - 1) / per_claim;
	const idx_t threads = TaskScheduler::GetScheduler(context).NumberOfThreads();
	return std::max<idx_t>(1, std::min(claims, threads));
}

void RequireUnambiguousTips(const char *fn, const std::string &tree_table, const TipMatches &tips, const char *who) {
	const std::string *first = nullptr;
	for (const auto &tip : tips) {
		if (tip.second.count > 1 && (!first || tip.first < *first)) {
			first = &tip.first;
		}
	}
	if (first) {
		throw InvalidInputException("%s: tip name '%s' names %d tips in the tree '%s'; %s must name exactly one tip",
		                            fn, *first, static_cast<int64_t>(tips.at(*first).count), tree_table, who);
	}
}

} // namespace duckdb
