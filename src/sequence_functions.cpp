#include "sequence_functions.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <array>

namespace duckdb {

// DNA complement lookup table for all 256 ASCII characters
// 0 means invalid character
static constexpr std::array<char, 256> CreateDnaComplementTable() {
	std::array<char, 256> table = {};
	// Initialize all to 0 (invalid)
	for (size_t i = 0; i < 256; i++) {
		table[i] = 0;
	}
	// Uppercase DNA bases
	table['A'] = 'T';
	table['T'] = 'A';
	table['G'] = 'C';
	table['C'] = 'G';
	// Uppercase IUPAC ambiguity codes
	table['R'] = 'Y'; // A or G -> T or C
	table['Y'] = 'R'; // C or T -> G or A
	table['S'] = 'S'; // G or C -> C or G
	table['W'] = 'W'; // A or T -> T or A
	table['K'] = 'M'; // G or T -> C or A
	table['M'] = 'K'; // A or C -> T or G
	table['B'] = 'V'; // not A -> not T
	table['D'] = 'H'; // not C -> not G
	table['H'] = 'D'; // not G -> not C
	table['V'] = 'B'; // not T -> not A
	table['N'] = 'N'; // any -> any
	// Lowercase DNA bases
	table['a'] = 't';
	table['t'] = 'a';
	table['g'] = 'c';
	table['c'] = 'g';
	// Lowercase IUPAC ambiguity codes
	table['r'] = 'y';
	table['y'] = 'r';
	table['s'] = 's';
	table['w'] = 'w';
	table['k'] = 'm';
	table['m'] = 'k';
	table['b'] = 'v';
	table['d'] = 'h';
	table['h'] = 'd';
	table['v'] = 'b';
	table['n'] = 'n';
	// Gap characters
	table['-'] = '-';
	table['.'] = '.';
	return table;
}

// RNA complement lookup table for all 256 ASCII characters
// 0 means invalid character
static constexpr std::array<char, 256> CreateRnaComplementTable() {
	std::array<char, 256> table = {};
	// Initialize all to 0 (invalid)
	for (size_t i = 0; i < 256; i++) {
		table[i] = 0;
	}
	// Uppercase RNA bases
	table['A'] = 'U';
	table['U'] = 'A';
	table['G'] = 'C';
	table['C'] = 'G';
	// Uppercase IUPAC ambiguity codes
	table['R'] = 'Y'; // A or G -> U or C
	table['Y'] = 'R'; // C or U -> G or A
	table['S'] = 'S'; // G or C -> C or G
	table['W'] = 'W'; // A or U -> U or A
	table['K'] = 'M'; // G or U -> C or A
	table['M'] = 'K'; // A or C -> U or G
	table['B'] = 'V'; // not A -> not U
	table['D'] = 'H'; // not C -> not G
	table['H'] = 'D'; // not G -> not C
	table['V'] = 'B'; // not U -> not A
	table['N'] = 'N'; // any -> any
	// Lowercase RNA bases
	table['a'] = 'u';
	table['u'] = 'a';
	table['g'] = 'c';
	table['c'] = 'g';
	// Lowercase IUPAC ambiguity codes
	table['r'] = 'y';
	table['y'] = 'r';
	table['s'] = 's';
	table['w'] = 'w';
	table['k'] = 'm';
	table['m'] = 'k';
	table['b'] = 'v';
	table['d'] = 'h';
	table['h'] = 'd';
	table['v'] = 'b';
	table['n'] = 'n';
	// Gap characters
	table['-'] = '-';
	table['.'] = '.';
	return table;
}

static constexpr auto DNA_COMPLEMENT_TABLE = CreateDnaComplementTable();
static constexpr auto RNA_COMPLEMENT_TABLE = CreateRnaComplementTable();

struct DnaReverseComplementOperator {
	template <class INPUT_TYPE, class RESULT_TYPE>
	static RESULT_TYPE Operation(INPUT_TYPE input, Vector &result) {
		auto input_data = input.GetData();
		auto input_len = input.GetSize();

		// Pre-allocate result string for performance
		auto result_str = StringVector::EmptyString(result, input_len);
		auto result_data = result_str.GetDataWriteable();

		// Reverse complement: iterate input in reverse, compute complement
		for (idx_t i = 0; i < input_len; i++) {
			unsigned char base = static_cast<unsigned char>(input_data[input_len - 1 - i]);
			char complement = DNA_COMPLEMENT_TABLE[base];

			if (complement == 0) {
				throw InvalidInputException("Invalid DNA base '%c' at position %llu", static_cast<char>(base),
				                            input_len - i);
			}

			result_data[i] = complement;
		}

		result_str.Finalize();
		return result_str;
	}
};

struct RnaReverseComplementOperator {
	template <class INPUT_TYPE, class RESULT_TYPE>
	static RESULT_TYPE Operation(INPUT_TYPE input, Vector &result) {
		auto input_data = input.GetData();
		auto input_len = input.GetSize();

		// Pre-allocate result string for performance
		auto result_str = StringVector::EmptyString(result, input_len);
		auto result_data = result_str.GetDataWriteable();

		// Reverse complement: iterate input in reverse, compute complement
		for (idx_t i = 0; i < input_len; i++) {
			unsigned char base = static_cast<unsigned char>(input_data[input_len - 1 - i]);
			char complement = RNA_COMPLEMENT_TABLE[base];

			if (complement == 0) {
				throw InvalidInputException("Invalid RNA base '%c' at position %llu", static_cast<char>(base),
				                            input_len - i);
			}

			result_data[i] = complement;
		}

		result_str.Finalize();
		return result_str;
	}
};

static void SequenceDnaReverseComplementFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::ExecuteString<string_t, string_t, DnaReverseComplementOperator>(args.data[0], result, args.size());
}

static void SequenceRnaReverseComplementFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::ExecuteString<string_t, string_t, RnaReverseComplementOperator>(args.data[0], result, args.size());
}

void SequenceFunctions::Register(ExtensionLoader &loader) {
	ScalarFunction sequence_dna_reverse_complement("sequence_dna_reverse_complement", {LogicalType::VARCHAR},
	                                               LogicalType::VARCHAR, SequenceDnaReverseComplementFunction);
	loader.RegisterFunction(sequence_dna_reverse_complement);

	ScalarFunction sequence_rna_reverse_complement("sequence_rna_reverse_complement", {LogicalType::VARCHAR},
	                                               LogicalType::VARCHAR, SequenceRnaReverseComplementFunction);
	loader.RegisterFunction(sequence_rna_reverse_complement);
}

} // namespace duckdb
