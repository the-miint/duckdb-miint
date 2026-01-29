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

		} else {
			throw InvalidInputException("Invalid type parameter for alignment_seq_identity. "
			                            "Must be 'gap_excluded', 'blast', or 'gap_compressed'.");
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

// alignment_query_length implementation
static void AlignmentQueryLengthScalarFunction(DataChunk &args, ExpressionState &state, Vector &result) {
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

ScalarFunction AlignmentQueryLengthFunction::GetFunction() {
	ScalarFunction func("alignment_query_length", {LogicalType::VARCHAR, LogicalType::BOOLEAN}, LogicalType::BIGINT,
	                    AlignmentQueryLengthScalarFunction);

	// Allow NULL CIGAR (returns NULL)
	func.null_handling = FunctionNullHandling::SPECIAL_HANDLING;

	// Set default value for include_hard_clips parameter (defaults to true)
	func.arguments[1] = LogicalType::BOOLEAN;
	func.varargs = LogicalType::INVALID;

	return func;
}

void AlignmentQueryLengthFunction::Register(ExtensionLoader &loader) {
	// Register overload with both parameters
	ScalarFunction func_two_params = GetFunction();

	// Register overload with single parameter (include_hard_clips defaults to true)
	ScalarFunction func_one_param(
	    "alignment_query_length", {LogicalType::VARCHAR}, LogicalType::BIGINT,
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
	ScalarFunctionSet function_set("alignment_query_length");
	function_set.AddFunction(func_one_param);
	function_set.AddFunction(func_two_params);
	loader.RegisterFunction(function_set);
}

// alignment_query_coverage implementation
static void AlignmentQueryCoverageScalarFunction(DataChunk &args, ExpressionState &state, Vector &result) {
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

ScalarFunction AlignmentQueryCoverageFunction::GetFunction() {
	ScalarFunction func("alignment_query_coverage", {LogicalType::VARCHAR, LogicalType::VARCHAR}, LogicalType::DOUBLE,
	                    AlignmentQueryCoverageScalarFunction);

	// Allow NULL values (returns NULL for NULL CIGAR, error for invalid type)
	func.null_handling = FunctionNullHandling::SPECIAL_HANDLING;

	// Set default value for type parameter (defaults to 'aligned')
	func.arguments[1] = LogicalType::VARCHAR;
	func.varargs = LogicalType::INVALID;

	return func;
}

void AlignmentQueryCoverageFunction::Register(ExtensionLoader &loader) {
	// Register overload with both parameters
	ScalarFunction func_two_params = GetFunction();

	// Register overload with single parameter (type defaults to 'aligned')
	ScalarFunction func_one_param(
	    "alignment_query_coverage", {LogicalType::VARCHAR}, LogicalType::DOUBLE,
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
	ScalarFunctionSet function_set("alignment_query_coverage");
	function_set.AddFunction(func_one_param);
	function_set.AddFunction(func_two_params);
	loader.RegisterFunction(function_set);
}

} // namespace duckdb
