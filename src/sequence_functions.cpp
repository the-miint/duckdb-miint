#include "sequence_functions.hpp"
#include "sequence_utils.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/limits.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <array>

namespace duckdb {

// Import shared complement tables from sequence_utils.hpp
static constexpr auto &DNA_COMPLEMENT_TABLE = miint::DNA_COMPLEMENT_TABLE;
static constexpr auto &RNA_COMPLEMENT_TABLE = miint::RNA_COMPLEMENT_TABLE;

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

// sequence_split(seq, chunk_size) -> LIST(STRUCT(chunk_index INTEGER, chunk_data VARCHAR))
//
// Fixed-width chunking in a single linear pass. Replaces the quadratic SQL macro
// (list_transform + substring inside a lambda, which loses the ASCII fast path and
// degrades to O(L^2)). Here each chunk is a direct byte slice -> O(L) total, with the
// same bounded ~O(L)/record memory profile as the list_transform shape it replaces.
//   - last chunk = remainder (no padding); chunk_index dense, 0-based, ascending.
//   - empty sequence -> empty list (0 chunks); NULL seq/chunk_size -> NULL.
//
// The O(L^2) blow-up is DuckDB issue #23229: a lambda-captured column loses the
// statistics that select substring's O(1) ASCII fast path, so substring falls back to
// the Unicode path and rescans from byte 0 on every element. If that is fixed upstream,
// the pure-SQL list_transform macro becomes linear again -- prefer that and retire this
// function. https://github.com/duckdb/duckdb/issues/23229
static const LogicalType &SequenceSplitReturnType() {
	static const LogicalType type = LogicalType::LIST(
	    LogicalType::STRUCT({{"chunk_index", LogicalType::INTEGER}, {"chunk_data", LogicalType::VARCHAR}}));
	return type;
}

static void SequenceSplitFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	const idx_t count = args.size();

	UnifiedVectorFormat seq_fmt, cs_fmt;
	args.data[0].ToUnifiedFormat(count, seq_fmt);
	args.data[1].ToUnifiedFormat(count, cs_fmt);
	const auto seq_data = UnifiedVectorFormat::GetData<string_t>(seq_fmt);
	const auto cs_data = UnifiedVectorFormat::GetData<int32_t>(cs_fmt);

	auto list_entries = FlatVector::GetData<list_entry_t>(result);
	auto &result_validity = FlatVector::Validity(result);

	// Pass 1: validate, NULL-propagate, and lay out per-row [offset, length) into the child.
	idx_t total_chunks = 0;
	for (idx_t row = 0; row < count; row++) {
		const auto si = seq_fmt.sel->get_index(row);
		const auto ci = cs_fmt.sel->get_index(row);
		if (!seq_fmt.validity.RowIsValid(si) || !cs_fmt.validity.RowIsValid(ci)) {
			result_validity.SetInvalid(row);
			list_entries[row].offset = total_chunks;
			list_entries[row].length = 0;
			continue;
		}
		const int32_t cs = cs_data[ci];
		if (cs <= 0) {
			throw InvalidInputException("sequence_split: chunk_size must be > 0 (got %d)", cs);
		}
		const idx_t len = seq_data[si].GetSize();
		const idx_t k = (len + static_cast<idx_t>(cs) - 1) / static_cast<idx_t>(cs); // ceil; len 0 -> 0
		// chunk_index is INTEGER; fail loud rather than silently wrap if the chunk count would
		// exceed it. Only reachable with a pathologically small chunk_size on a multi-GB record.
		if (k > static_cast<idx_t>(NumericLimits<int32_t>::Maximum())) {
			throw InvalidInputException("sequence_split: %llu-byte sequence at chunk_size %d yields %llu chunks, "
			                            "exceeding the INTEGER chunk_index limit (%d); use a larger chunk_size",
			                            static_cast<unsigned long long>(len), cs, static_cast<unsigned long long>(k),
			                            NumericLimits<int32_t>::Maximum());
		}
		list_entries[row].offset = total_chunks;
		list_entries[row].length = k;
		total_chunks += k;
	}

	// Reserve the flat child once (no per-row realloc), then fill.
	ListVector::Reserve(result, total_chunks);
	ListVector::SetListSize(result, total_chunks);
	auto &struct_vec = ListVector::GetEntry(result);
	auto &struct_children = StructVector::GetEntries(struct_vec);
	auto idx_data = FlatVector::GetData<int32_t>(*struct_children[0]); // chunk_index
	auto &data_vec = *struct_children[1];                              // chunk_data VARCHAR
	auto chunk_str = FlatVector::GetData<string_t>(data_vec);

	// Pass 2: slice. chunk_data is a copy into the result heap (must outlive the input).
	for (idx_t row = 0; row < count; row++) {
		if (!result_validity.RowIsValid(row)) {
			continue;
		}
		const auto si = seq_fmt.sel->get_index(row);
		const auto ci = cs_fmt.sel->get_index(row);
		const idx_t cs = static_cast<idx_t>(cs_data[ci]);
		const string_t &s = seq_data[si];
		const char *ptr = s.GetData();
		const idx_t len = s.GetSize();
		const idx_t base = list_entries[row].offset;
		const idx_t k = list_entries[row].length;
		for (idx_t i = 0; i < k; i++) {
			const idx_t start = i * cs;
			const idx_t this_len = (start + cs <= len) ? cs : (len - start);
			const idx_t child_idx = base + i;
			idx_data[child_idx] = static_cast<int32_t>(i);
			chunk_str[child_idx] = StringVector::AddString(data_vec, ptr + start, this_len);
		}
	}
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

	ScalarFunction sequence_split("sequence_split", {LogicalType::VARCHAR, LogicalType::INTEGER},
	                              SequenceSplitReturnType(), SequenceSplitFunction);
	sequence_split.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	loader.RegisterFunction(sequence_split);
}

} // namespace duckdb
