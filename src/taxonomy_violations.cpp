#include "taxonomy_violations.hpp"

#include "TaxonomyStore.hpp"
#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "tree_function_utils.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/sql_identifier.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/vector/flat_vector.hpp"

#include <algorithm>
#include <array>
#include <string>
#include <utility>
#include <vector>

namespace duckdb {

namespace {

constexpr const char *kFn = "taxonomy_violations";
constexpr size_t kRanks = miint::kLineageRanks.size();

struct BindData : public TableFunctionData {
	std::string taxonomy;
	LogicalType id_type;
	// The rank columns' types as bound, in kLineageRanks order: all VARCHAR, collations kept.
	std::array<LogicalType, kRanks> rank_types;
};

// One output row. The pointers point into GlobalState::rows, or at static text; nullptr is NULL.
struct Violation {
	const char *kind;
	// The input row it names by feature_id.
	size_t row;
	const char *rank;
	const std::string *name;
	const std::string *parent;
	const char *detail;
};

struct GlobalState : public GlobalTableFunctionState {
	std::vector<miint::TaxonomyRow> rows;
	// Whether each row's feature_id is NULL; its text is then "".
	std::vector<bool> id_null;
	std::vector<Violation> violations;
	size_t next = 0;
};

// The type of feature_id, after checking that `columns` are exactly the shared lineage schema's:
// feature_id, of an id type, and the eight rank columns, VARCHAR. Matched by name, in any order and
// case, so a misspelled or extra column is an error rather than ignored.
LogicalType ValidateLineageSchema(const std::string &taxonomy, const TableOrViewColumns &columns,
                                  std::array<LogicalType, kRanks> &rank_types) {
	for (const auto &name : columns.names) {
		bool known = StringUtil::CIEquals(name, "feature_id");
		for (const auto &rank : miint::kLineageRanks) {
			known = known || StringUtil::CIEquals(name, rank.column);
		}
		if (!known) {
			throw BinderException("%s: taxonomy '%s' has a column '%s' outside the shared lineage schema; select "
			                      "feature_id and the eight rank columns, domain through strain",
			                      kFn, taxonomy, name);
		}
	}
	LogicalType id_type;
	if (!FindColumnType(columns, "feature_id", id_type)) {
		throw BinderException("%s: taxonomy '%s' has no 'feature_id' column", kFn, taxonomy);
	}
	if (!IsAllowedIdType(id_type)) {
		throw BinderException("%s: feature_id must be %s, but it is %s; cast it", kFn, AllowedIdTypeList(),
		                      id_type.ToString());
	}
	for (size_t r = 0; r < kRanks; r++) {
		const char *column = miint::kLineageRanks[r].column;
		if (!FindColumnType(columns, column, rank_types[r])) {
			throw BinderException("%s: taxonomy '%s' has no '%s' column", kFn, taxonomy, column);
		}
		if (rank_types[r].id() != LogicalTypeId::VARCHAR) {
			throw BinderException("%s: column '%s' of taxonomy '%s' is %s; the rank columns must be VARCHAR", kFn,
			                      column, taxonomy, rank_types[r].ToString());
		}
	}
	return id_type;
}

unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input, vector<LogicalType> &return_types,
                              vector<Identifier> &names) {
	auto data = make_uniq<BindData>();
	data->taxonomy = input.inputs[0].ToString();
	RejectCTERelationName(input, data->taxonomy);

	// Schemas only: no rows are read until InitGlobal (docs/internals/architecture.md).
	auto columns = GetTableOrViewColumns(context, data->taxonomy, std::string(kFn) + ": taxonomy table");
	data->id_type = ValidateLineageSchema(data->taxonomy, columns, data->rank_types);

	names.emplace_back("kind");
	return_types.push_back(LogicalType::VARCHAR);
	names.emplace_back("feature_id");
	return_types.push_back(data->id_type);
	for (const char *column : {"rank", "name", "parent", "detail"}) {
		names.emplace_back(column);
		return_types.push_back(LogicalType::VARCHAR);
	}
	return std::move(data);
}

const char *KindText(miint::TaxonomyViolationKind kind) {
	switch (kind) {
	case miint::TaxonomyViolationKind::InvalidId:
		return "invalid_id";
	case miint::TaxonomyViolationKind::InvalidValue:
		return "invalid_value";
	case miint::TaxonomyViolationKind::DuplicateId:
		return "duplicate_id";
	case miint::TaxonomyViolationKind::Gap:
		return "gap";
	}
	throw InternalException("%s: unknown violation kind", kFn);
}

// The taxonomy's rows, in the relation's order, ids as text. A NULL id is "" with id_null set,
// since its text cannot tell it from an empty VARCHAR id.
void ReadTaxonomy(ClientContext &context, const BindData &bind, std::vector<miint::TaxonomyRow> &rows,
                  std::vector<bool> &id_null) {
	std::vector<std::pair<std::string, LogicalType>> columns {{SQLIdentifier::ToString("feature_id"), bind.id_type}};
	for (size_t r = 0; r < kRanks; r++) {
		columns.emplace_back(SQLIdentifier::ToString(miint::kLineageRanks[r].column), bind.rank_types[r]);
	}
	std::vector<std::string> ids;
	std::vector<bool> nulls;
	ReadRelation(context, kFn, bind.taxonomy, columns, [&](DataChunk &chunk) {
		ExtractIdColumnAsStrings(chunk, 0, bind.id_type, ids, nulls);
		std::array<UnifiedVectorFormat, kRanks> ranks;
		for (size_t r = 0; r < kRanks; r++) {
			chunk.data[r + 1].ToUnifiedFormat(ranks[r]);
		}
		for (idx_t i = 0; i < chunk.size(); i++) {
			miint::TaxonomyRow row {std::move(ids[i]), {}};
			for (size_t r = 0; r < kRanks; r++) {
				const auto idx = ranks[r].sel->get_index(i);
				if (ranks[r].validity.RowIsValid(idx)) {
					row.ranks[r] = UnifiedVectorFormat::GetData<string_t>(ranks[r])[idx].GetString();
				}
			}
			rows.push_back(std::move(row));
			id_null.push_back(nulls[i]);
		}
	});
}

unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind = input.bind_data->Cast<BindData>();
	auto state = make_uniq<GlobalState>();
	ReadTaxonomy(context, bind, state->rows, state->id_null);

	const auto &rows = state->rows;
	auto &violations = state->violations;
	for (const auto &v : miint::FindRowViolations(rows)) {
		Violation out {KindText(v.kind), v.row, nullptr, nullptr, nullptr, nullptr};
		switch (v.kind) {
		case miint::TaxonomyViolationKind::InvalidId:
			out.detail = state->id_null[v.row] ? "is NULL" : "is empty";
			break;
		case miint::TaxonomyViolationKind::InvalidValue:
			out.rank = miint::kLineageRanks[v.rank].column;
			out.name = &*rows[v.row].ranks[v.rank];
			out.detail = v.reason;
			break;
		case miint::TaxonomyViolationKind::DuplicateId:
			break;
		case miint::TaxonomyViolationKind::Gap:
			out.rank = miint::kLineageRanks[v.rank].column;
			break;
		}
		violations.push_back(out);
	}
	// One row per input row under each parent of a label, so the rows at fault can be found;
	// an empty parent rank is a NULL parent.
	for (const auto &label : miint::FindMultipleParents(rows)) {
		for (const auto &under : label.rows_by_parent) {
			for (const size_t row : under) {
				const auto &ranks = rows[row].ranks;
				const auto &parent = ranks[label.rank - 1];
				violations.push_back({"multiple_parents", row, miint::kLineageRanks[label.rank].column,
				                      &*ranks[label.rank], parent ? &*parent : nullptr, nullptr});
			}
		}
	}
	return std::move(state);
}

void SetText(Vector &out, idx_t row, const char *text) {
	auto &validity = FlatVector::ValidityMutable(out);
	if (!text) {
		validity.SetInvalid(row);
		return;
	}
	FlatVector::GetDataMutable<string_t>(out)[row] = StringVector::AddString(out, text);
	validity.SetValid(row);
}

void SetText(Vector &out, idx_t row, const std::string *text) {
	auto &validity = FlatVector::ValidityMutable(out);
	if (!text) {
		validity.SetInvalid(row);
		return;
	}
	FlatVector::GetDataMutable<string_t>(out)[row] = StringVector::AddString(out, *text);
	validity.SetValid(row);
}

void Execute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind = data_p.bind_data->Cast<BindData>();
	auto &state = data_p.global_state->Cast<GlobalState>();
	const idx_t count = std::min<idx_t>(STANDARD_VECTOR_SIZE, state.violations.size() - state.next);
	for (idx_t i = 0; i < count; i++) {
		const auto &v = state.violations[state.next + i];
		SetText(output.data[0], i, v.kind);
		if (state.id_null[v.row]) {
			FlatVector::ValidityMutable(output.data[1]).SetInvalid(i);
		} else {
			EmitIdCell(output.data[1], i, state.rows[v.row].feature_id, bind.id_type);
		}
		SetText(output.data[2], i, v.rank);
		SetText(output.data[3], i, v.name);
		SetText(output.data[4], i, v.parent);
		SetText(output.data[5], i, v.detail);
	}
	state.next += count;
	output.SetChildCardinality(count);
}

} // namespace

TableFunction TaxonomyViolationsTableFunction::GetFunction() {
	TableFunction tf(kFn, {LogicalType::VARCHAR}, Execute, Bind, InitGlobal);
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void TaxonomyViolationsTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
