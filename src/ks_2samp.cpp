#include "ks_2samp.hpp"
#include "alignment_functions_internal.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/types/vector.hpp"

namespace duckdb {

// ks_2samp(a, b [, method]) -> STRUCT(statistic DOUBLE, pvalue DOUBLE)
//
// Two-sided two-sample Kolmogorov-Smirnov test (issue #218), micov's last reason to
// depend on SciPy. The algorithm lives in miint::KsTwoSample; this is only the
// DuckDB glue.
//
// NULL handling. A NULL row must be marked NULL on the struct AND on both of its
// fields, which is what SetStructRowNull (table_function_common.hpp) does; marking only
// the struct leaks the children's contents through struct_extract, and this function
// shipped with that bug once.
//
//   * a NULL list, or a NULL `method`, yields NULL -- ordinary NULL-in/NULL-out,
//     which is why this function keeps DEFAULT null handling rather than declaring
//     SPECIAL_HANDLING. Marking it special would tell the optimizer we might return
//     a real value from a NULL input, which is not true, and would forfeit the
//     NULL-propagation it can otherwise do.
//   * an EMPTY list also yields NULL. The statistic is undefined without at least
//     one observation on each side, and returning 0.0 there would read as "the two
//     samples are identical" -- the opposite of what an empty group means.
//   * NULL ELEMENTS INSIDE a list are dropped and the test runs on the remainder,
//     matching compress_intervals and the QC aggregates. A list whose every element
//     is NULL therefore reduces to the empty case and yields NULL.
//
// NaN is rejected outright rather than dropped: it has no position in the sort
// order, so it is a data error rather than an absence, and silently discarding it
// would change the sample size a p-value is computed against.
//
// `method` is validated BEFORE the NULL checks on purpose, so that a bogus method is
// reported even for rows whose data is NULL and whose result is therefore NULL
// anyway. That is load-bearing for a NULL column value -- verified: a table row of
// (NULL, [3.0]) with method 'bogus' raises, whereas checking validity first would
// return NULL and swallow the typo.
//
// It is NOT unconditional, and the limit is worth knowing: a LITERAL NULL argument
// is constant-folded away before the function runs at all. `EXPLAIN SELECT
// ks_2samp(NULL, [3.0], 'bogus')` shows a plain `NULL` projection over a dummy scan,
// so nothing validates the method there. That is the price of DEFAULT null handling
// and it is the right trade: the alternative, declaring SPECIAL_HANDLING, would
// forfeit NULL propagation on every real query to catch a typo in a query that has
// no data in it.

namespace {

LogicalType KsTwoSampleReturnType() {
	return LogicalType::STRUCT({{"statistic", LogicalType::DOUBLE}, {"pvalue", LogicalType::DOUBLE}});
}

// The non-NULL elements of one list row.
void CollectSample(UnifiedVectorFormat &list_fmt, UnifiedVectorFormat &child_fmt, idx_t row_idx,
                   std::vector<double> &out) {
	out.clear();
	auto list_entries = UnifiedVectorFormat::GetData<list_entry_t>(list_fmt);
	auto &entry = list_entries[list_fmt.sel->get_index(row_idx)];
	auto child_values = UnifiedVectorFormat::GetData<double>(child_fmt);
	for (idx_t k = 0; k < entry.length; k++) {
		const auto ci = child_fmt.sel->get_index(entry.offset + k);
		if (!child_fmt.validity.RowIsValid(ci)) {
			continue;
		}
		out.push_back(child_values[ci]);
	}
}

void KsTwoSampleScalarFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	const idx_t count = args.size();
	const bool has_method = args.ColumnCount() >= 3;

	UnifiedVectorFormat a_fmt;
	UnifiedVectorFormat b_fmt;
	args.data[0].ToUnifiedFormat(count, a_fmt);
	args.data[1].ToUnifiedFormat(count, b_fmt);

	UnifiedVectorFormat a_child_fmt;
	UnifiedVectorFormat b_child_fmt;
	ListVector::GetEntry(args.data[0]).ToUnifiedFormat(ListVector::GetListSize(args.data[0]), a_child_fmt);
	ListVector::GetEntry(args.data[1]).ToUnifiedFormat(ListVector::GetListSize(args.data[1]), b_child_fmt);

	UnifiedVectorFormat method_fmt;
	if (has_method) {
		args.data[2].ToUnifiedFormat(count, method_fmt);
	}

	auto &entries = StructVector::GetEntries(result);
	auto statistic_data = FlatVector::GetData<double>(*entries[0]);
	auto pvalue_data = FlatVector::GetData<double>(*entries[1]);

	// Hoisted: the payload pointer does not change per row.
	auto method_values = has_method ? UnifiedVectorFormat::GetData<string_t>(method_fmt) : nullptr;

	std::vector<double> a;
	std::vector<double> b;

	// One handler for the whole scan. Every miint::InvalidInputException raised below is a
	// hard error that aborts the query, so per-statement try blocks bought nothing except
	// three copies of the same rethrow.
	try {
		for (idx_t i = 0; i < count; i++) {
			bool method_is_null = false;
			if (has_method) {
				const auto mi = method_fmt.sel->get_index(i);
				if (!method_fmt.validity.RowIsValid(mi)) {
					method_is_null = true;
				} else {
					miint::ValidateKsMethod(method_values[mi].GetString());
				}
			}

			const auto ai = a_fmt.sel->get_index(i);
			const auto bi = b_fmt.sel->get_index(i);
			if (method_is_null || !a_fmt.validity.RowIsValid(ai) || !b_fmt.validity.RowIsValid(bi)) {
				SetStructRowNull(result, i);
				continue;
			}

			CollectSample(a_fmt, a_child_fmt, i, a);
			CollectSample(b_fmt, b_child_fmt, i, b);

			// NaN is checked BEFORE the empty short circuit, or a NaN in one sample would go
			// unreported whenever the other sample is empty -- ks_2samp([], ['nan']) is a
			// data error, not an absence, and the docs promise NaN always raises.
			miint::RejectKsNaN(a, b);

			// The SQL contract turns an empty sample into NULL; the algorithm treats it
			// as a precondition violation and throws. This is where the two meet.
			if (a.empty() || b.empty()) {
				SetStructRowNull(result, i);
				continue;
			}

			const auto r = miint::KsTwoSample(a, b);
			statistic_data[i] = r.statistic;
			pvalue_data[i] = r.pvalue;
		}
	} catch (const miint::InvalidInputException &e) {
		throw InvalidInputException("%s", e.what());
	}
}

} // namespace

void KsTwoSampleFunction::Register(ExtensionLoader &loader) {
	const auto list_double = LogicalType::LIST(LogicalType::DOUBLE);

	ScalarFunction two_arg("ks_2samp", {list_double, list_double}, KsTwoSampleReturnType(), KsTwoSampleScalarFunction);
	ScalarFunction three_arg("ks_2samp", {list_double, list_double, LogicalType::VARCHAR}, KsTwoSampleReturnType(),
	                         KsTwoSampleScalarFunction);

	ScalarFunctionSet function_set("ks_2samp");
	function_set.AddFunction(two_arg);
	function_set.AddFunction(three_arg);
	loader.RegisterFunction(function_set);
}

} // namespace duckdb
