#include "alignment_functions.hpp"
#include "alignment_functions_internal.hpp"
#include "documented_function.hpp"
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
			// Requires CIGAR to use =/X (not M). Does not need NM or MD tags.
			// Formula: match_ops / alignment_columns
			// Returns NULL if:
			//   - CIGAR uses only M (ambiguous — can't distinguish matches from mismatches)
			//   - CIGAR mixes M with =/X (inconsistent)
			if (cigar_stats.match_ops + cigar_stats.mismatch_ops == 0) {
				// No = or X ops observed (M-only, or degenerate S/N/P-only CIGAR)
				result_validity.SetInvalid(i);
				continue;
			}

			// If M ops are present alongside =/X, the CIGAR is inconsistent — return NULL
			int64_t m_only = cigar_stats.matches - cigar_stats.match_ops - cigar_stats.mismatch_ops;
			if (m_only > 0) {
				result_validity.SetInvalid(i);
				continue;
			}

			// alignment_columns is guaranteed > 0 here because =/X ops exist
			identity = static_cast<double>(cigar_stats.match_ops) / static_cast<double>(cigar_stats.alignment_columns);

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
	static const std::string description = R"DOC(
Calculate sequence identity between read and reference using one of four
methods derived from Heng Li's [blog post on the definition of sequence
identity](https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity).

**IMPORTANT:** these calculations assume the alignment uses `=`/`X` operations
(e.g. `--xeq` with bowtie2, `-eqx` with minimap2, etc.) unless the `cigar` type
is selected.

### Identity types

The `type` parameter selects which formula to apply (default `'gap_compressed'`):

1. **`'gap_excluded'`** — ignore gaps, only consider match/mismatch positions.
   - Formula: `#matches / (#matches + #mismatches)`
   - Requires: CIGAR + MD tag
   - Use case: genetic divergence between species

2. **`'blast'`** — traditional BLAST-style identity.
   - Formula: `#matches / alignment_columns`
   - Requires: CIGAR + NM tag
   - Use case: general similarity measurement
   - Note: large indels significantly lower identity

3. **`'gap_compressed'`** (default) — count consecutive gaps as single events.
   - Formula: `(m - n + g) / (m + o)` where m=M_columns, n=NM, g=gap_bases, o=gap_opens
   - Equivalent to: `1 - (n - g + o) / (m + o)` from the blog post
   - Requires: CIGAR + NM tag
   - Use case: filtering alignments (recommended)
   - Note: more robust to structural variations

4. **`'cigar'`** — identity from extended CIGAR operations only (no tags needed).
   - Formula: `match_ops / alignment_columns` where `match_ops` = count of `=` ops
   - Requires: CIGAR with `=`/`X` ops (not `M`). NM and MD tags are ignored.
   - Returns NULL if CIGAR uses `M` (ambiguous) or mixes `M` with `=`/`X` (inconsistent).
   - Use case: computing identity on trimmed alignments from `alignment_slice`,
     where tags have been invalidated.

**Returns** a DOUBLE between 0.0 and 1.0, or NULL for unmapped reads or missing
required tags.
)DOC";
	RegisterDocumentedScalar(
	    loader, GetFunction(), description, {"cigar", "nm", "md", "type"},
	    {
	        // examples

	        "-- Calculate gap-compressed identity (default)\n"
	        "SELECT read_id, alignment_seq_identity(cigar, tag_nm, tag_md, 'gap_compressed') AS identity\n"
	        "FROM read_alignments('alignments.sam')\n"
	        "WHERE tag_nm IS NOT NULL;",

	        "-- Filter high-quality alignments\n"
	        "SELECT COUNT(*)\n"
	        "FROM read_alignments('alignments.bam')\n"
	        "WHERE alignment_seq_identity(cigar, tag_nm, tag_md) > 0.95;",

	        "-- Compare different identity methods\n"
	        "SELECT read_id,\n"
	        "  alignment_seq_identity(cigar, tag_nm, tag_md, 'gap_excluded') AS gap_excl,\n"
	        "  alignment_seq_identity(cigar, tag_nm, tag_md, 'blast') AS blast,\n"
	        "  alignment_seq_identity(cigar, tag_nm, tag_md, 'gap_compressed') AS gap_comp\n"
	        "FROM read_alignments('alignments.sam')\n"
	        "WHERE tag_nm IS NOT NULL AND tag_md IS NOT NULL;",

	        "-- Compute identity on sliced alignments (tags NULLed, cigar method still works)\n"
	        "SELECT read_id, alignment_seq_identity(cigar, tag_nm, tag_md, 'cigar') AS identity\n"
	        "FROM alignment_slice('my_alignments', 1000, 2000);",
	    },
	    /*alias_of=*/"", /*categories=*/ {"alignment-quality"});
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

	static const std::string description = R"DOC(
Calculate the total query length from a CIGAR string. Useful for understanding
read lengths and computing query coverage manually.

### Behavior

- When `include_hard_clips=true` (default): returns `M + I + S + = + X + H` —
  every query-consuming operation.
- When `include_hard_clips=false`: returns `M + I + S + = + X` — matches
  HTSlib's `bam_cigar2qlen`.
- Soft clips (`S`) are always included (they're present in the sequence field).
- Hard clips (`H`) are included only when the parameter is true.
- Deletions (`D`) and reference skips (`N`) don't consume query, so they're
  never counted.

**Returns** BIGINT — total query length, or 0 for NULL/unmapped reads.
)DOC";
	const std::vector<std::string> shared_examples = {
	    "-- Get query length including hard clips (default)\n"
	    "SELECT read_id, alignment_query_length(cigar) AS query_len\n"
	    "FROM read_alignments('alignments.sam');",

	    "-- Get query length excluding hard clips (matches bam_cigar2qlen)\n"
	    "SELECT read_id, alignment_query_length(cigar, false) AS query_len\n"
	    "FROM read_alignments('alignments.sam');",

	    "-- Compare lengths with and without hard clips\n"
	    "SELECT read_id, cigar,\n"
	    "  alignment_query_length(cigar, true) AS len_with_hard,\n"
	    "  alignment_query_length(cigar, false) AS len_without_hard\n"
	    "FROM read_alignments('alignments.sam')\n"
	    "WHERE cigar LIKE '%H%';",

	    "-- Calculate average query length per reference\n"
	    "SELECT reference, AVG(alignment_query_length(cigar)) AS avg_query_len\n"
	    "FROM read_alignments('alignments.bam')\n"
	    "WHERE NOT alignment_is_unmapped(flags)\n"
	    "GROUP BY reference;",
	};

	// Construct LogicalType from LogicalTypeId enum to avoid ODR-using the
	// LogicalType::VARCHAR / BOOLEAN static-const class members; binding them to
	// vector::emplace_back's forwarding reference emits a vague-linkage weak
	// symbol that collides with libduckdb_static.a under the loadable
	// extension's link settings.
	FunctionDescription one_param_desc;
	one_param_desc.description = description;
	one_param_desc.parameter_names.emplace_back("cigar");
	one_param_desc.parameter_types.emplace_back(LogicalTypeId::VARCHAR);
	one_param_desc.categories.emplace_back("alignment-quality");
	for (const auto &e : shared_examples) {
		one_param_desc.examples.emplace_back(e);
	}

	FunctionDescription two_param_desc;
	two_param_desc.description = description;
	two_param_desc.parameter_names.emplace_back("cigar");
	two_param_desc.parameter_names.emplace_back("include_hard_clips");
	two_param_desc.parameter_types.emplace_back(LogicalTypeId::VARCHAR);
	two_param_desc.parameter_types.emplace_back(LogicalTypeId::BOOLEAN);
	two_param_desc.categories.emplace_back("alignment-quality");
	for (const auto &e : shared_examples) {
		two_param_desc.examples.emplace_back(e);
	}

	RegisterDocumentedScalarSet(loader, function_set, {one_param_desc, two_param_desc});
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

	static const std::string description = R"DOC(
Calculate the proportion of query bases covered by the reference alignment.
Useful for assessing how much of a read actually aligns versus being clipped.

### Coverage types

The `type` parameter selects which formula to apply (default `'aligned'`):

1. **`'aligned'`** (default) — bases that align to the reference.
   - Formula: `(M + = + X) / (M + I + S + = + X + H)`
   - Only counts bases that match/mismatch the reference.
   - Insertions and clips reduce coverage.

2. **`'mapped'`** — bases that are mapped (not clipped).
   - Formula: `(M + I + = + X) / (M + I + S + = + X + H)`
   - Counts insertions as "mapped" even though they don't align.
   - Only clips reduce coverage.

### Behavior

- The denominator (query length) always includes hard clips.
- Returns 0.0 for reads with only clipping operations (no alignment).
- Deletions (`D`) and reference skips (`N`) don't affect coverage (they don't
  consume query).

**Returns** a DOUBLE between 0.0 and 1.0, or NULL for NULL CIGAR.

### Use cases

- **`'aligned'`**: assess alignment quality (how much of the read actually
  matches the reference).
- **`'mapped'`**: identify heavily clipped reads (adapters, chimeras,
  low-quality ends).
)DOC";
	const std::vector<std::string> shared_examples = {
	    "-- Get aligned coverage (default)\n"
	    "SELECT read_id, alignment_query_coverage(cigar) AS aligned_cov\n"
	    "FROM read_alignments('alignments.sam');",

	    "-- Get mapped coverage (includes insertions)\n"
	    "SELECT read_id, alignment_query_coverage(cigar, 'mapped') AS mapped_cov\n"
	    "FROM read_alignments('alignments.sam');",

	    "-- Compare aligned vs mapped coverage\n"
	    "SELECT read_id, cigar,\n"
	    "  alignment_query_coverage(cigar, 'aligned') AS aligned_cov,\n"
	    "  alignment_query_coverage(cigar, 'mapped') AS mapped_cov\n"
	    "FROM read_alignments('alignments.sam')\n"
	    "WHERE cigar LIKE '%I%';  -- Reads with insertions show the difference",

	    "-- Filter reads with high query coverage\n"
	    "SELECT COUNT(*)\n"
	    "FROM read_alignments('alignments.bam')\n"
	    "WHERE alignment_query_coverage(cigar, 'aligned') > 0.9;",

	    "-- Find heavily clipped reads\n"
	    "SELECT read_id, cigar, alignment_query_coverage(cigar) AS coverage\n"
	    "FROM read_alignments('alignments.sam')\n"
	    "WHERE alignment_query_coverage(cigar) < 0.5\n"
	    "ORDER BY coverage;",

	    "-- Calculate average coverage per reference\n"
	    "SELECT reference,\n"
	    "  AVG(alignment_query_coverage(cigar, 'aligned')) AS avg_aligned_cov,\n"
	    "  AVG(alignment_query_coverage(cigar, 'mapped')) AS avg_mapped_cov\n"
	    "FROM read_alignments('alignments.bam')\n"
	    "WHERE NOT alignment_is_unmapped(flags)\n"
	    "GROUP BY reference;",
	};

	// LogicalTypeId enum: see ODR comment in alignment_query_length above.
	FunctionDescription one_param_desc;
	one_param_desc.description = description;
	one_param_desc.parameter_names.emplace_back("cigar");
	one_param_desc.parameter_types.emplace_back(LogicalTypeId::VARCHAR);
	one_param_desc.categories.emplace_back("alignment-quality");
	for (const auto &e : shared_examples) {
		one_param_desc.examples.emplace_back(e);
	}

	FunctionDescription two_param_desc;
	two_param_desc.description = description;
	two_param_desc.parameter_names.emplace_back("cigar");
	two_param_desc.parameter_names.emplace_back("type");
	two_param_desc.parameter_types.emplace_back(LogicalTypeId::VARCHAR);
	two_param_desc.parameter_types.emplace_back(LogicalTypeId::VARCHAR);
	two_param_desc.categories.emplace_back("alignment-quality");
	for (const auto &e : shared_examples) {
		two_param_desc.examples.emplace_back(e);
	}

	RegisterDocumentedScalarSet(loader, function_set, {one_param_desc, two_param_desc});
}

} // namespace duckdb
