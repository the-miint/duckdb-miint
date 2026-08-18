#include "alignment_functions.hpp"
#include "alignment_functions_internal.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/string_type.hpp"
#include "duckdb/function/scalar_function.hpp"
#include <string>

namespace duckdb {

// Main function implementation
static void AlignmentSeqIdentityScalarFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &cigar_vector = args.data[0];
	auto &nm_vector = args.data[1];
	auto &md_vector = args.data[2];
	auto &type_vector = args.data[3];

	// Manually handle 4 arguments since DuckDB doesn't have QuaternaryExecutor
	UnifiedVectorFormat cigar_data, nm_data, md_data, type_data;
	cigar_vector.ToUnifiedFormat(args.size(), cigar_data);
	nm_vector.ToUnifiedFormat(args.size(), nm_data);
	md_vector.ToUnifiedFormat(args.size(), md_data);
	type_vector.ToUnifiedFormat(args.size(), type_data);

	auto cigar_ptr = UnifiedVectorFormat::GetData<string_t>(cigar_data);
	auto nm_ptr = UnifiedVectorFormat::GetData<int64_t>(nm_data);
	auto md_ptr = UnifiedVectorFormat::GetData<string_t>(md_data);
	auto type_ptr = UnifiedVectorFormat::GetData<string_t>(type_data);

	auto result_data = FlatVector::GetData<double>(result);
	auto &result_validity = FlatVector::Validity(result);

	for (idx_t i = 0; i < args.size(); i++) {
		auto cigar_idx = cigar_data.sel->get_index(i);
		auto nm_idx = nm_data.sel->get_index(i);
		auto md_idx = md_data.sel->get_index(i);
		auto type_idx = type_data.sel->get_index(i);

		// Check validity for required parameters (cigar and type)
		if (!cigar_data.validity.RowIsValid(cigar_idx) || !type_data.validity.RowIsValid(type_idx)) {
			result_validity.SetInvalid(i);
			continue;
		}

		string_t cigar = cigar_ptr[cigar_idx];
		string_t type = type_ptr[type_idx];

		// Handle optional parameters - treat NULL as missing (-1 for nm, empty for md)
		int64_t nm = nm_data.validity.RowIsValid(nm_idx) ? nm_ptr[nm_idx] : -1;
		string_t md = md_data.validity.RowIsValid(md_idx) ? md_ptr[md_idx] : string_t("", 0);

		// Handle NULL or unmapped CIGAR
		if (cigar.GetSize() == 0 || (cigar.GetSize() == 1 && cigar.GetData()[0] == '*')) {
			result_validity.SetInvalid(i);
			continue;
		}

		// Get type string
		std::string type_str(type.GetData(), type.GetSize());

		// Parse CIGAR
		miint::CigarStats cigar_stats;
		try {
			std::string cigar_std(cigar.GetData(), cigar.GetSize());
			cigar_stats = miint::ParseCigar(cigar_std);
		} catch (const miint::InvalidInputException &) {
			result_validity.SetInvalid(i);
			continue;
		}

		double identity = 0.0;

		if (type_str == "gap_excluded") {
			// gap_excluded: #matches / (#matches + #mismatches)
			// Requires MD tag
			if (md.GetSize() == 0) {
				result_validity.SetInvalid(i);
				continue;
			}

			miint::MdStats md_stats;
			try {
				std::string md_std(md.GetData(), md.GetSize());
				md_stats = miint::ParseMd(md_std);
			} catch (const miint::InvalidInputException &) {
				result_validity.SetInvalid(i);
				continue;
			}

			int64_t total = md_stats.matches + md_stats.mismatches;
			if (total <= 0) {
				result_validity.SetInvalid(i);
				continue;
			}

			identity = static_cast<double>(md_stats.matches) / static_cast<double>(total);

		} else if (type_str == "blast") {
			// blast: #matches / alignment_columns
			// #matches = alignment_columns - NM
			// Requires NM tag
			if (nm < 0) {
				result_validity.SetInvalid(i);
				continue;
			}

			if (cigar_stats.alignment_columns <= 0) {
				result_validity.SetInvalid(i);
				continue;
			}

			// Validate NM tag doesn't exceed alignment length (per SAM spec)
			if (nm > cigar_stats.alignment_columns) {
				result_validity.SetInvalid(i);
				continue;
			}

			int64_t matches = cigar_stats.alignment_columns - nm;
			identity = static_cast<double>(matches) / static_cast<double>(cigar_stats.alignment_columns);

		} else if (type_str == "gap_compressed") {
			// gap_compressed: sequence identity with gap compression
			// Formula: 1 - (n - g + o) / (m + o) = (m - n + g) / (m + o)
			// where:
			//   m = sum(M/=/X) - match operations in CIGAR
			//   n = NM tag - edit distance
			//   g = sum(I+D) - total gap bases (insertions + deletions)
			//   o = gap_opens - number of gap opening events
			// This treats consecutive indel operations as a single event.
			// Reference: Heng Li's blog post "On the definition of sequence identity"
			// https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
			// Requires NM tag
			if (nm < 0) {
				result_validity.SetInvalid(i);
				continue;
			}

			int64_t m = cigar_stats.matches;
			int64_t n = nm;
			int64_t g = cigar_stats.insertions + cigar_stats.deletions;
			int64_t o = cigar_stats.gap_opens;

			// Validate NM tag is reasonable (shouldn't exceed matches + gaps)
			if (n > m + g) {
				result_validity.SetInvalid(i);
				continue;
			}

			int64_t denominator = m + o;
			if (denominator <= 0) {
				result_validity.SetInvalid(i);
				continue;
			}

			int64_t numerator = m - n + g;
			identity = static_cast<double>(numerator) / static_cast<double>(denominator);

		} else if (type_str == "cigar") {
			// cigar: identity from extended CIGAR ops (= and X) only.
			// Math lives in miint::ComputeCigarIdentity (internal helper); std::nullopt
			// means "can't compute from CIGAR alone" (M-only, mixed M+=/X, or degenerate).
			auto maybe_identity = miint::ComputeCigarIdentity(cigar_stats);
			if (!maybe_identity.has_value()) {
				result_validity.SetInvalid(i);
				continue;
			}
			identity = *maybe_identity;

		} else {
			throw InvalidInputException("Invalid type parameter for alignment_seq_identity: '%s'. "
			                            "Must be 'gap_excluded', 'blast', 'gap_compressed', or 'cigar'.",
			                            type_str);
		}

		result_data[i] = identity;
	}
}

ScalarFunction AlignmentSeqIdentityFunction::GetFunction() {
	ScalarFunction func("alignment_seq_identity",
	                    {LogicalType::VARCHAR, LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::VARCHAR},
	                    LogicalType::DOUBLE, AlignmentSeqIdentityScalarFunction);

	// Allow NULL values for optional parameters (nm and md)
	func.null_handling = FunctionNullHandling::SPECIAL_HANDLING;

	// Set default value for type parameter
	func.arguments[3] = LogicalType::VARCHAR;
	func.varargs = LogicalType::INVALID;

	return func;
}

void AlignmentSeqIdentityFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

// cigar_sequence_identity(cigar) — one-arg convenience over the type='cigar'
// branch of alignment_seq_identity. Same math (via miint::ComputeCigarIdentity).
// Returns NULL when identity can't be computed from CIGAR alone (M-only CIGAR,
// mixed M+=/X, degenerate CIGARs with no =/X ops).
static void CigarSequenceIdentityScalarFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &cigar_vector = args.data[0];

	UnifiedVectorFormat cigar_data;
	cigar_vector.ToUnifiedFormat(args.size(), cigar_data);
	auto cigar_ptr = UnifiedVectorFormat::GetData<string_t>(cigar_data);

	result.SetVectorType(VectorType::FLAT_VECTOR);
	auto result_data = FlatVector::GetData<double>(result);
	auto &result_validity = FlatVector::Validity(result);

	for (idx_t i = 0; i < args.size(); i++) {
		auto cigar_idx = cigar_data.sel->get_index(i);

		if (!cigar_data.validity.RowIsValid(cigar_idx)) {
			result_validity.SetInvalid(i);
			continue;
		}

		auto cigar = cigar_ptr[cigar_idx];
		if (cigar.GetSize() == 0 || (cigar.GetSize() == 1 && cigar.GetData()[0] == '*')) {
			result_validity.SetInvalid(i);
			continue;
		}

		try {
			std::string cigar_std(cigar.GetData(), cigar.GetSize());
			miint::CigarStats cigar_stats = miint::ParseCigar(cigar_std);

			auto maybe_identity = miint::ComputeCigarIdentity(cigar_stats);
			if (!maybe_identity.has_value()) {
				result_validity.SetInvalid(i);
				continue;
			}
			result_data[i] = *maybe_identity;
		} catch (const miint::InvalidInputException &e) {
			throw InvalidInputException(e.what());
		}
	}
}

ScalarFunction CigarSequenceIdentityFunction::GetFunction() {
	ScalarFunction func("cigar_sequence_identity", {LogicalType::VARCHAR}, LogicalType::DOUBLE,
	                    CigarSequenceIdentityScalarFunction);
	func.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	return func;
}

void CigarSequenceIdentityFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

// cigar_query_length implementation
static void CigarQueryLengthScalarFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &cigar_vector = args.data[0];
	auto &include_hard_clips_vector = args.data[1];

	BinaryExecutor::Execute<string_t, bool, int64_t>(
	    cigar_vector, include_hard_clips_vector, result, args.size(), [&](string_t cigar, bool include_hard_clips) {
		    // Handle NULL or unmapped CIGAR
		    if (cigar.GetSize() == 0 || (cigar.GetSize() == 1 && cigar.GetData()[0] == '*')) {
			    return int64_t(0);
		    }

		    try {
			    // Parse CIGAR
			    std::string cigar_std(cigar.GetData(), cigar.GetSize());
			    miint::CigarStats cigar_stats = miint::ParseCigar(cigar_std);

			    return miint::ComputeQueryLength(cigar_stats, include_hard_clips);
		    } catch (const miint::InvalidInputException &e) {
			    // Convert miint exceptions to DuckDB exceptions
			    throw InvalidInputException(e.what());
		    }
	    });
}

ScalarFunction CigarQueryLengthFunction::GetFunction() {
	ScalarFunction func("cigar_query_length", {LogicalType::VARCHAR, LogicalType::BOOLEAN}, LogicalType::BIGINT,
	                    CigarQueryLengthScalarFunction);

	// Allow NULL CIGAR (returns NULL)
	func.null_handling = FunctionNullHandling::SPECIAL_HANDLING;

	// Set default value for include_hard_clips parameter (defaults to true)
	func.arguments[1] = LogicalType::BOOLEAN;
	func.varargs = LogicalType::INVALID;

	return func;
}

void CigarQueryLengthFunction::Register(ExtensionLoader &loader) {
	// Register overload with both parameters
	ScalarFunction func_two_params = GetFunction();

	// Register overload with single parameter (include_hard_clips defaults to true)
	ScalarFunction func_one_param(
	    "cigar_query_length", {LogicalType::VARCHAR}, LogicalType::BIGINT,
	    [](DataChunk &args, ExpressionState &state, Vector &result) {
		    UnaryExecutor::Execute<string_t, int64_t>(args.data[0], result, args.size(), [&](string_t cigar) {
			    // Handle NULL or unmapped CIGAR
			    if (cigar.GetSize() == 0 || (cigar.GetSize() == 1 && cigar.GetData()[0] == '*')) {
				    return int64_t(0);
			    }

			    try {
				    // Parse CIGAR - let exceptions propagate for invalid input
				    std::string cigar_std(cigar.GetData(), cigar.GetSize());
				    miint::CigarStats cigar_stats = miint::ParseCigar(cigar_std);

				    // Default to include_hard_clips = true
				    return miint::ComputeQueryLength(cigar_stats, true);
			    } catch (const miint::InvalidInputException &e) {
				    // Convert miint exceptions to DuckDB exceptions
				    throw InvalidInputException(e.what());
			    }
		    });
	    });
	func_one_param.null_handling = FunctionNullHandling::SPECIAL_HANDLING;

	// Register both overloads as a function set
	ScalarFunctionSet function_set("cigar_query_length");
	function_set.AddFunction(func_one_param);
	function_set.AddFunction(func_two_params);
	loader.RegisterFunction(function_set);
}

// cigar_query_coverage implementation
static void CigarQueryCoverageScalarFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &cigar_vector = args.data[0];
	auto &type_vector = args.data[1];

	BinaryExecutor::Execute<string_t, string_t, double>(
	    cigar_vector, type_vector, result, args.size(), [&](string_t cigar, string_t type) {
		    // Handle NULL or unmapped CIGAR - return 0.0 for empty/unmapped
		    if (cigar.GetSize() == 0 || (cigar.GetSize() == 1 && cigar.GetData()[0] == '*')) {
			    return 0.0;
		    }

		    try {
			    // Parse CIGAR - let exceptions propagate for invalid input
			    std::string cigar_std(cigar.GetData(), cigar.GetSize());
			    miint::CigarStats cigar_stats = miint::ParseCigar(cigar_std);

			    // Get type string
			    std::string type_str(type.GetData(), type.GetSize());

			    // Compute coverage
			    return miint::ComputeQueryCoverage(cigar_stats, type_str);
		    } catch (const miint::InvalidInputException &e) {
			    // Convert miint exceptions to DuckDB exceptions
			    throw InvalidInputException(e.what());
		    }
	    });
}

ScalarFunction CigarQueryCoverageFunction::GetFunction() {
	ScalarFunction func("cigar_query_coverage", {LogicalType::VARCHAR, LogicalType::VARCHAR}, LogicalType::DOUBLE,
	                    CigarQueryCoverageScalarFunction);

	// Allow NULL values (returns NULL for NULL CIGAR, error for invalid type)
	func.null_handling = FunctionNullHandling::SPECIAL_HANDLING;

	// Set default value for type parameter (defaults to 'aligned')
	func.arguments[1] = LogicalType::VARCHAR;
	func.varargs = LogicalType::INVALID;

	return func;
}

void CigarQueryCoverageFunction::Register(ExtensionLoader &loader) {
	// Register overload with both parameters
	ScalarFunction func_two_params = GetFunction();

	// Register overload with single parameter (type defaults to 'aligned')
	ScalarFunction func_one_param(
	    "cigar_query_coverage", {LogicalType::VARCHAR}, LogicalType::DOUBLE,
	    [](DataChunk &args, ExpressionState &state, Vector &result) {
		    UnaryExecutor::Execute<string_t, double>(args.data[0], result, args.size(), [&](string_t cigar) {
			    // Handle NULL or unmapped CIGAR - return 0.0 for empty/unmapped
			    if (cigar.GetSize() == 0 || (cigar.GetSize() == 1 && cigar.GetData()[0] == '*')) {
				    return 0.0;
			    }

			    try {
				    // Parse CIGAR - let exceptions propagate for invalid input
				    std::string cigar_std(cigar.GetData(), cigar.GetSize());
				    miint::CigarStats cigar_stats = miint::ParseCigar(cigar_std);

				    // Default to type = 'aligned'
				    return miint::ComputeQueryCoverage(cigar_stats, "aligned");
			    } catch (const miint::InvalidInputException &e) {
				    // Convert miint exceptions to DuckDB exceptions
				    throw InvalidInputException(e.what());
			    }
		    });
	    });
	func_one_param.null_handling = FunctionNullHandling::SPECIAL_HANDLING;

	// Register both overloads as a function set
	ScalarFunctionSet function_set("cigar_query_coverage");
	function_set.AddFunction(func_one_param);
	function_set.AddFunction(func_two_params);
	loader.RegisterFunction(function_set);
}

// cigar_query_intervals(cigar, flags, [type='aligned']) implementation
//
// Returns the query positions this alignment covers as half-open [start, stop)
// intervals in the ORIGINAL READ's orientation. cigar_query_coverage answers "how much
// of the read does this one record explain" and cannot see sibling records; a read
// spanning the origin of a circular reference is emitted as several records, and
// pooling them requires their query footprints on a common axis. Hence intervals rather
// than a count.
//
// `start`/`stop` deliberately match compress_intervals' struct field names so the two
// compose without renaming.
static const LogicalType &CigarQueryIntervalsReturnType() {
	static const LogicalType type =
	    LogicalType::LIST(LogicalType::STRUCT({{"start", LogicalType::BIGINT}, {"stop", LogicalType::BIGINT}}));
	return type;
}

// One body serves both the two- and three-argument overloads, keyed off ColumnCount().
static void CigarQueryIntervalsScalarFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	const idx_t count = args.size();
	const bool has_type = args.ColumnCount() > 2;

	UnifiedVectorFormat cigar_fmt, flags_fmt, type_fmt;
	args.data[0].ToUnifiedFormat(count, cigar_fmt);
	args.data[1].ToUnifiedFormat(count, flags_fmt);
	if (has_type) {
		args.data[2].ToUnifiedFormat(count, type_fmt);
	}
	const auto cigar_data = UnifiedVectorFormat::GetData<string_t>(cigar_fmt);
	const auto flags_data = UnifiedVectorFormat::GetData<uint16_t>(flags_fmt);
	const auto type_data = has_type ? UnifiedVectorFormat::GetData<string_t>(type_fmt) : nullptr;

	result.SetVectorType(VectorType::FLAT_VECTOR);
	auto list_entries = FlatVector::GetData<list_entry_t>(result);
	auto &result_validity = FlatVector::Validity(result);

	// The interval count is only known by computing the intervals, so there is no cheap
	// sizing pass as in sequence_split. Rather than buffer the whole chunk and copy it in
	// afterwards -- which would double peak memory for no benefit -- grow the list child
	// per row and write straight into it, the way compute_coverage_depth's Finalize does.
	idx_t total = 0;
	for (idx_t row = 0; row < count; row++) {
		const auto cigar_idx = cigar_fmt.sel->get_index(row);
		const auto flags_idx = flags_fmt.sel->get_index(row);
		list_entries[row].offset = total;
		list_entries[row].length = 0;

		const idx_t type_idx = has_type ? type_fmt.sel->get_index(row) : 0;
		const bool valid = cigar_fmt.validity.RowIsValid(cigar_idx) && flags_fmt.validity.RowIsValid(flags_idx) &&
		                   (!has_type || type_fmt.validity.RowIsValid(type_idx));
		if (!valid) {
			result_validity.SetInvalid(row);
			continue;
		}

		// An unmapped or empty CIGAR needs no special case: ParseCigarOperations yields no
		// operations for it and ComputeQueryIntervals then covers nothing, which is the
		// empty list rather than an unknown one. Going through the normal path also means
		// the `type` argument is validated on unmapped rows, so a typo'd type fails on
		// every batch instead of only on batches that happen to hold a mapped record.
		// Deliberately stricter than cigar_query_coverage, which returns 0.0 for an
		// unmapped CIGAR without inspecting its type at all.
		try {
			auto ops = miint::ParseCigarOperations(cigar_data[cigar_idx].GetString());
			// 0x10 is the SAM reverse-strand bit; the CIGAR is written in reference
			// orientation, so it decides which end of the read the clips sit at.
			auto intervals = miint::ComputeQueryIntervals(ops, (flags_data[flags_idx] & 0x10) != 0,
			                                              has_type ? type_data[type_idx].GetString() : "aligned");
			if (intervals.empty()) {
				continue;
			}
			// Reserve may reallocate the child, so re-fetch its data pointers each time.
			ListVector::Reserve(result, total + intervals.size());
			auto &struct_children = StructVector::GetEntries(ListVector::GetEntry(result));
			auto start_data = FlatVector::GetData<int64_t>(*struct_children[0]);
			auto stop_data = FlatVector::GetData<int64_t>(*struct_children[1]);
			for (idx_t i = 0; i < intervals.size(); i++) {
				start_data[total + i] = intervals[i].start;
				stop_data[total + i] = intervals[i].stop;
			}
			list_entries[row].length = intervals.size();
			total += intervals.size();
		} catch (const miint::InvalidInputException &e) {
			throw InvalidInputException(e.what());
		}
	}

	ListVector::SetListSize(result, total);
}

ScalarFunction CigarQueryIntervalsFunction::GetFunction() {
	ScalarFunction func("cigar_query_intervals", {LogicalType::VARCHAR, LogicalType::USMALLINT, LogicalType::VARCHAR},
	                    CigarQueryIntervalsReturnType(), CigarQueryIntervalsScalarFunction);
	func.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	return func;
}

void CigarQueryIntervalsFunction::Register(ExtensionLoader &loader) {
	// Two-argument overload; type defaults to 'aligned' inside the shared body.
	ScalarFunction func_two_args("cigar_query_intervals", {LogicalType::VARCHAR, LogicalType::USMALLINT},
	                             CigarQueryIntervalsReturnType(), CigarQueryIntervalsScalarFunction);
	func_two_args.null_handling = FunctionNullHandling::SPECIAL_HANDLING;

	ScalarFunctionSet function_set("cigar_query_intervals");
	function_set.AddFunction(func_two_args);
	function_set.AddFunction(GetFunction());
	loader.RegisterFunction(function_set);
}

// cigar_pooled_identity(cigar) — sequence identity for a read the aligner reported as
// several alignment records, as sum(=) / sum(alignment columns) over the group.
//
// An aggregate rather than cigar_sequence_identity(string_agg(cigar, '')): concatenating
// CIGARs produces a string the SAM spec forbids, because clipping is only legal at an end,
// and its arithmetic came out right only because ParseCigar accumulates operation by
// operation without validating structure. Summing the counters needs no such indulgence,
// generalises to metrics that concatenation would get wrong (anything adjacency-sensitive,
// gap_opens above all), and drops both the per-group string copy and a whole reparse.
//
// The state is miint::IdentityCounts itself: four counters, no allocation, and Combine is
// four adds. On a single record the result is cigar_sequence_identity of that record,
// including every input where that is NULL — the rules live in one place and are applied
// once, to the sums.
struct CigarPooledIdentityOperation {
	template <class STATE>
	static void Initialize(STATE &state) {
		state = miint::IdentityCounts();
	}

	template <class INPUT_TYPE, class STATE, class OP>
	static void Operation(STATE &state, const INPUT_TYPE &input, AggregateUnaryInput &unary_input) {
		try {
			std::string cigar(input.GetData(), input.GetSize());
			state.Add(miint::ToIdentityCounts(miint::ParseCigar(cigar)));
		} catch (const miint::InvalidInputException &e) {
			throw InvalidInputException(e.what());
		}
	}

	template <class INPUT_TYPE, class STATE, class OP>
	static void ConstantOperation(STATE &state, const INPUT_TYPE &input, AggregateUnaryInput &unary_input,
	                              idx_t count) {
		for (idx_t i = 0; i < count; i++) {
			Operation<INPUT_TYPE, STATE, OP>(state, input, unary_input);
		}
	}

	template <class STATE, class OP>
	static void Combine(const STATE &source, STATE &target, AggregateInputData &aggr_input_data) {
		target.Add(source);
	}

	template <class T, class STATE>
	static void Finalize(STATE &state, T &target, AggregateFinalizeData &finalize_data) {
		auto identity = miint::ComputeCigarIdentity(state);
		if (!identity.has_value()) {
			finalize_data.ReturnNull();
			return;
		}
		target = *identity;
	}

	// A record with no CIGAR contributes no operations, so it contributes nothing. An
	// unmapped or empty CIGAR reaches ParseCigar and lands in the same place.
	static bool IgnoreNull() {
		return true;
	}
};

AggregateFunction CigarPooledIdentityFunction::GetFunction() {
	auto function =
	    AggregateFunction::UnaryAggregate<miint::IdentityCounts, string_t, double, CigarPooledIdentityOperation>(
	        LogicalType::VARCHAR, LogicalType::DOUBLE);
	function.name = "cigar_pooled_identity";
	// Four commutative adds; the result cannot depend on the order records arrive in. Stated
	// explicitly because DuckDB assumes the opposite by default, which costs the window path
	// its partial-aggregate reuse and wraps the aggregate in a sort it does not need -- and
	// windowing this over a mate pair is a documented use.
	function.SetOrderDependent(AggregateOrderDependent::NOT_ORDER_DEPENDENT);
	return function;
}

void CigarPooledIdentityFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
