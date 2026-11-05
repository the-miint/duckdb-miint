#include "parse_lineage.hpp"

#include "table_function_common.hpp"
#include "taxonomy_lineage.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/common/vector/string_vector.hpp"
#include "duckdb/common/vector/struct_vector.hpp"
#include "duckdb/execution/expression_executor.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"

#include <stdexcept>
#include <string>
#include <string_view>

namespace duckdb {

static LogicalType ParseLineageResultType() {
	child_list_t<LogicalType> fields;
	for (const auto &rank : miint::kLineageRanks) {
		fields.emplace_back(rank.column, LogicalType::VARCHAR);
	}
	return LogicalType::STRUCT(std::move(fields));
}

struct ParseLineageBindData : public FunctionData {
	explicit ParseLineageBindData(char domain_letter) : domain_letter(domain_letter) {
	}

	unique_ptr<FunctionData> Copy() const override {
		return make_uniq<ParseLineageBindData>(domain_letter);
	}
	bool Equals(const FunctionData &other) const override {
		return domain_letter == other.Cast<ParseLineageBindData>().domain_letter;
	}

	char domain_letter;
};

// domain_letter must be one constant for the whole call: a per-row letter would let
// one taxonomy mix d__ and k__ rows, which the letter exists to prevent.
static unique_ptr<FunctionData> ParseLineageBind(BindScalarFunctionInput &input) {
	auto &args = input.GetArguments();
	if (!args[1]->IsFoldable()) {
		throw InvalidInputException("parse_lineage: domain_letter must be a constant 'd' or 'k'; a column or a "
		                            "prepared-statement parameter is not accepted");
	}
	Value letter = ExpressionExecutor::EvaluateScalar(input.GetClientContext(), *args[1]);
	if (letter.IsNull()) {
		throw InvalidInputException("parse_lineage: domain_letter must be 'd' or 'k', not NULL");
	}
	auto text = letter.GetValue<string>();
	if (text != "d" && text != "k") {
		throw InvalidInputException("parse_lineage: domain_letter must be 'd' or 'k', got '%s'", text);
	}
	return make_uniq<ParseLineageBindData>(text[0]);
}

static void ParseLineageExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &func_expr = state.expr.Cast<BoundFunctionExpression>();
	const char domain_letter = func_expr.BindInfo()->Cast<ParseLineageBindData>().domain_letter;

	UnifiedVectorFormat lineage_data;
	args.data[0].ToUnifiedFormat(lineage_data);
	auto lineages = UnifiedVectorFormat::GetData<string_t>(lineage_data);
	auto &columns = StructVector::GetEntries(result);

	for (idx_t row = 0; row < args.size(); ++row) {
		auto idx = lineage_data.sel->get_index(row);
		if (!lineage_data.validity.RowIsValid(idx)) {
			SetStructRowNull(result, row);
			continue;
		}
		miint::LineageRanks ranks;
		try {
			const auto &lineage = lineages[idx];
			ranks = miint::ParseLineageString(std::string_view(lineage.GetData(), lineage.GetSize()), domain_letter);
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("parse_lineage: %s", e.what());
		}
		for (idx_t col = 0; col < ranks.size(); ++col) {
			auto &column = columns[col];
			if (ranks[col]) {
				FlatVector::GetDataMutable<string_t>(column)[row] = StringVector::AddString(column, *ranks[col]);
			} else {
				FlatVector::SetNull(column, row, true);
			}
		}
	}
}

void ParseLineageFunction::Register(ExtensionLoader &loader) {
	ScalarFunction fn("parse_lineage",
	                  {FunctionParameter(Identifier("lineage"), LogicalType::VARCHAR),
	                   FunctionParameter(Identifier("domain_letter"), LogicalType::VARCHAR, Value("d"),
	                                     FunctionParameterKind::KEYWORD_ONLY)},
	                  ParseLineageResultType(), ParseLineageExecute, ParseLineageBind);
	fn.SetNullHandling(FunctionNullHandling::SPECIAL_HANDLING);
	fn.SetFallible();
	loader.RegisterFunction(fn);
}

} // namespace duckdb
