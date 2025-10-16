#include "alignment_functions.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/string_type.hpp"
#include "duckdb/function/scalar_function.hpp"
#include <cctype>
#include <string>

namespace miint {

// CIGAR parsing result structure
struct CigarStats {
	int64_t matches = 0;           // M, =, X operations
	int64_t match_ops = 0;         // = operations only
	int64_t mismatch_ops = 0;      // X operations only
	int64_t insertions = 0;        // I operations
	int64_t deletions = 0;         // D operations
	int64_t gap_opens = 0;         // Number of gap opening events
	int64_t alignment_columns = 0; // M + I + D operations
};

// Parse CIGAR string and extract statistics
static CigarStats ParseCigar(const std::string &cigar_str) {
	CigarStats stats;
	const char *cigar = cigar_str.data();
	size_t len = cigar_str.size();

	if (len == 0 || (len == 1 && cigar[0] == '*')) {
		return stats; // Empty or unmapped
	}

	int64_t op_len = 0;
	char prev_op_type = '\0'; // Track previous operation type for gap opens

	for (size_t i = 0; i < len; i++) {
		char c = cigar[i];

		if (std::isdigit(c)) {
			op_len = op_len * 10 + (c - '0');
		} else {
			// Operation character
			if (op_len == 0) {
				throw duckdb::InvalidInputException("Invalid CIGAR string: operation without length");
			}

			switch (c) {
			case 'M': // Match or mismatch (alignment match)
				stats.matches += op_len;
				stats.alignment_columns += op_len;
				break;
			case '=': // Sequence match
				stats.matches += op_len;
				stats.match_ops += op_len;
				stats.alignment_columns += op_len;
				break;
			case 'X': // Sequence mismatch
				stats.matches += op_len;
				stats.mismatch_ops += op_len;
				stats.alignment_columns += op_len;
				break;
			case 'I': // Insertion to the reference
				stats.insertions += op_len;
				stats.alignment_columns += op_len;
				// Count gap open if previous op was not I
				if (prev_op_type != 'I') {
					stats.gap_opens++;
				}
				break;
			case 'D': // Deletion from the reference
				stats.deletions += op_len;
				stats.alignment_columns += op_len;
				// Count gap open if previous op was not D
				if (prev_op_type != 'D') {
					stats.gap_opens++;
				}
				break;
			case 'N': // Skipped region (ref skip, e.g., intron)
			case 'S': // Soft clipping
			case 'H': // Hard clipping
			case 'P': // Padding
				// These don't contribute to alignment columns
				break;
			default:
				throw duckdb::InvalidInputException("Invalid CIGAR operation: " + std::string(1, c));
			}

			prev_op_type = c;
			op_len = 0;
		}
	}

	return stats;
}

// MD parsing result structure
struct MdStats {
	int64_t matches = 0;
	int64_t mismatches = 0;
};

// Parse MD tag string and extract match/mismatch counts
static MdStats ParseMd(const std::string &md_str) {
	MdStats stats;
	const char *md = md_str.data();
	size_t len = md_str.size();

	if (len == 0) {
		return stats; // Empty MD tag
	}

	int64_t match_len = 0;

	for (size_t i = 0; i < len; i++) {
		char c = md[i];

		if (std::isdigit(c)) {
			match_len = match_len * 10 + (c - '0');
		} else if (c == '^') {
			// Deletion marker: skip following deleted bases
			if (match_len > 0) {
				stats.matches += match_len;
				match_len = 0;
			}
			i++; // Skip the '^'
			// Skip deletion bases
			while (i < len && std::isalpha(md[i])) {
				i++;
			}
			i--; // Back up one since loop will increment
		} else if (std::isalpha(c)) {
			// Mismatch base
			if (match_len > 0) {
				stats.matches += match_len;
				match_len = 0;
			}
			stats.mismatches++;
		}
	}

	// Add remaining matches
	if (match_len > 0) {
		stats.matches += match_len;
	}

	return stats;
}

} // namespace miint

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
		} catch (const InvalidInputException &) {
			result_validity.SetInvalid(i);
			continue;
		}

		double identity = 0.0;
		bool valid_result = true;

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
			} catch (const InvalidInputException &) {
				result_validity.SetInvalid(i);
				continue;
			}

			int64_t total = md_stats.matches + md_stats.mismatches;
			if (total == 0) {
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

			if (cigar_stats.alignment_columns == 0) {
				result_validity.SetInvalid(i);
				continue;
			}

			int64_t matches = cigar_stats.alignment_columns - nm;
			identity = static_cast<double>(matches) / static_cast<double>(cigar_stats.alignment_columns);

		} else if (type_str == "gap_compressed") {
			// gap_compressed: 1 - (n - g + o) / (m + o) = (m - n + g) / (m + o)
			// where m = sum(M/=/X), n = NM, g = sum(I+D), o = gap_opens
			// Requires NM tag
			if (nm < 0) {
				result_validity.SetInvalid(i);
				continue;
			}

			int64_t m = cigar_stats.matches;
			int64_t n = nm;
			int64_t g = cigar_stats.insertions + cigar_stats.deletions;
			int64_t o = cigar_stats.gap_opens;

			int64_t denominator = m + o;
			if (denominator == 0) {
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

} // namespace duckdb
