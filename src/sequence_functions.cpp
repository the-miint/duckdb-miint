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

// DNA regexp lookup table for all 256 ASCII characters
// nullptr means invalid character
static constexpr std::array<const char *, 256> CreateDnaRegexpTable() {
	std::array<const char *, 256> table = {};
	// Initialize all to nullptr (invalid)
	for (size_t i = 0; i < 256; i++) {
		table[i] = nullptr;
	}
	// Uppercase unambiguous bases
	table['A'] = "A";
	table['C'] = "C";
	table['G'] = "G";
	table['T'] = "T";
	// Uppercase IUPAC ambiguity codes
	table['R'] = "[AG]";   // A or G
	table['Y'] = "[CT]";   // C or T
	table['S'] = "[CG]";   // C or G
	table['W'] = "[AT]";   // A or T
	table['K'] = "[GT]";   // G or T
	table['M'] = "[AC]";   // A or C
	table['B'] = "[CGT]";  // not A
	table['D'] = "[AGT]";  // not C
	table['H'] = "[ACT]";  // not G
	table['V'] = "[ACG]";  // not T
	table['N'] = "[ACGT]"; // any
	// Lowercase unambiguous bases
	table['a'] = "a";
	table['c'] = "c";
	table['g'] = "g";
	table['t'] = "t";
	// Lowercase IUPAC ambiguity codes
	table['r'] = "[ag]";
	table['y'] = "[ct]";
	table['s'] = "[cg]";
	table['w'] = "[at]";
	table['k'] = "[gt]";
	table['m'] = "[ac]";
	table['b'] = "[cgt]";
	table['d'] = "[agt]";
	table['h'] = "[act]";
	table['v'] = "[acg]";
	table['n'] = "[acgt]";
	// Gap characters - match any character in regex
	table['-'] = ".";
	table['.'] = ".";
	return table;
}

// RNA regexp lookup table for all 256 ASCII characters
// nullptr means invalid character
static constexpr std::array<const char *, 256> CreateRnaRegexpTable() {
	std::array<const char *, 256> table = {};
	// Initialize all to nullptr (invalid)
	for (size_t i = 0; i < 256; i++) {
		table[i] = nullptr;
	}
	// Uppercase unambiguous bases
	table['A'] = "A";
	table['C'] = "C";
	table['G'] = "G";
	table['U'] = "U";
	// Uppercase IUPAC ambiguity codes
	table['R'] = "[AG]";   // A or G
	table['Y'] = "[CU]";   // C or U
	table['S'] = "[CG]";   // C or G
	table['W'] = "[AU]";   // A or U
	table['K'] = "[GU]";   // G or U
	table['M'] = "[AC]";   // A or C
	table['B'] = "[CGU]";  // not A
	table['D'] = "[AGU]";  // not C
	table['H'] = "[ACU]";  // not G
	table['V'] = "[ACG]";  // not U
	table['N'] = "[ACGU]"; // any
	// Lowercase unambiguous bases
	table['a'] = "a";
	table['c'] = "c";
	table['g'] = "g";
	table['u'] = "u";
	// Lowercase IUPAC ambiguity codes
	table['r'] = "[ag]";
	table['y'] = "[cu]";
	table['s'] = "[cg]";
	table['w'] = "[au]";
	table['k'] = "[gu]";
	table['m'] = "[ac]";
	table['b'] = "[cgu]";
	table['d'] = "[agu]";
	table['h'] = "[acu]";
	table['v'] = "[acg]";
	table['n'] = "[acgu]";
	// Gap characters - match any character in regex
	table['-'] = ".";
	table['.'] = ".";
	return table;
}

static constexpr auto DNA_REGEXP_TABLE = CreateDnaRegexpTable();
static constexpr auto RNA_REGEXP_TABLE = CreateRnaRegexpTable();

// Molecule type strings for error messages
static constexpr const char DNA_TYPE[] = "DNA";
static constexpr const char RNA_TYPE[] = "RNA";

// Templated reverse complement operator - works for both DNA and RNA
template <const std::array<char, 256> &COMPLEMENT_TABLE, const char *MOLECULE_TYPE>
struct ReverseComplementOperator {
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
			char complement = COMPLEMENT_TABLE[base];

			if (complement == 0) {
				throw InvalidInputException("Invalid %s base '%c' at position %llu", MOLECULE_TYPE,
				                            static_cast<char>(base), input_len - i);
			}

			result_data[i] = complement;
		}

		result_str.Finalize();
		return result_str;
	}
};

// Templated regexp conversion operator - works for both DNA and RNA
template <const std::array<const char *, 256> &REGEXP_TABLE, const char *MOLECULE_TYPE>
struct AsRegexpOperator {
	template <class INPUT_TYPE, class RESULT_TYPE>
	static RESULT_TYPE Operation(INPUT_TYPE input, Vector &result) {
		auto input_data = input.GetData();
		auto input_len = input.GetSize();

		// First pass: validate and calculate output length
		idx_t output_len = 0;
		for (idx_t i = 0; i < input_len; i++) {
			unsigned char base = static_cast<unsigned char>(input_data[i]);
			const char *regexp = REGEXP_TABLE[base];

			if (regexp == nullptr) {
				throw InvalidInputException("Invalid %s base '%c' at position %llu", MOLECULE_TYPE,
				                            static_cast<char>(base), i + 1);
			}

			// Calculate length of this regexp fragment
			const char *p = regexp;
			while (*p) {
				output_len++;
				p++;
			}
		}

		// Allocate result string with exact size
		auto result_str = StringVector::EmptyString(result, output_len);
		auto result_data = result_str.GetDataWriteable();

		// Second pass: build the regexp string
		idx_t output_pos = 0;
		for (idx_t i = 0; i < input_len; i++) {
			unsigned char base = static_cast<unsigned char>(input_data[i]);
			const char *regexp = REGEXP_TABLE[base];

			// Copy the regexp fragment
			while (*regexp) {
				result_data[output_pos++] = *regexp++;
			}
		}

		result_str.Finalize();
		return result_str;
	}
};

// Type aliases for DNA and RNA operators
using DnaReverseComplementOperator = ReverseComplementOperator<DNA_COMPLEMENT_TABLE, DNA_TYPE>;
using RnaReverseComplementOperator = ReverseComplementOperator<RNA_COMPLEMENT_TABLE, RNA_TYPE>;
using DnaAsRegexpOperator = AsRegexpOperator<DNA_REGEXP_TABLE, DNA_TYPE>;
using RnaAsRegexpOperator = AsRegexpOperator<RNA_REGEXP_TABLE, RNA_TYPE>;

static void SequenceDnaReverseComplementFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::ExecuteString<string_t, string_t, DnaReverseComplementOperator>(args.data[0], result, args.size());
}

static void SequenceRnaReverseComplementFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::ExecuteString<string_t, string_t, RnaReverseComplementOperator>(args.data[0], result, args.size());
}

static void SequenceDnaAsRegexpFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::ExecuteString<string_t, string_t, DnaAsRegexpOperator>(args.data[0], result, args.size());
}

static void SequenceRnaAsRegexpFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::ExecuteString<string_t, string_t, RnaAsRegexpOperator>(args.data[0], result, args.size());
}

void SequenceFunctions::Register(ExtensionLoader &loader) {
	ScalarFunction sequence_dna_reverse_complement("sequence_dna_reverse_complement", {LogicalType::VARCHAR},
	                                               LogicalType::VARCHAR, SequenceDnaReverseComplementFunction);
	loader.RegisterFunction(sequence_dna_reverse_complement);

	ScalarFunction sequence_rna_reverse_complement("sequence_rna_reverse_complement", {LogicalType::VARCHAR},
	                                               LogicalType::VARCHAR, SequenceRnaReverseComplementFunction);
	loader.RegisterFunction(sequence_rna_reverse_complement);

	ScalarFunction sequence_dna_as_regexp("sequence_dna_as_regexp", {LogicalType::VARCHAR}, LogicalType::VARCHAR,
	                                      SequenceDnaAsRegexpFunction);
	loader.RegisterFunction(sequence_dna_as_regexp);

	ScalarFunction sequence_rna_as_regexp("sequence_rna_as_regexp", {LogicalType::VARCHAR}, LogicalType::VARCHAR,
	                                      SequenceRnaAsRegexpFunction);
	loader.RegisterFunction(sequence_rna_as_regexp);
}

} // namespace duckdb
