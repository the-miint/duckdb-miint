#include "alignment_flag_functions.hpp"
#include "documented_function.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

static void AlignmentIsPairedFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &flags_vector = args.data[0];
	UnaryExecutor::Execute<uint16_t, bool>(flags_vector, result, args.size(),
	                                       [&](uint16_t flags) { return (flags & 0x1) != 0; });
}

static void AlignmentIsProperPairFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &flags_vector = args.data[0];
	UnaryExecutor::Execute<uint16_t, bool>(flags_vector, result, args.size(),
	                                       [&](uint16_t flags) { return (flags & 0x2) != 0; });
}

static void AlignmentIsUnmappedFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &flags_vector = args.data[0];
	UnaryExecutor::Execute<uint16_t, bool>(flags_vector, result, args.size(),
	                                       [&](uint16_t flags) { return (flags & 0x4) != 0; });
}

static void AlignmentIsMateUnmappedFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &flags_vector = args.data[0];
	UnaryExecutor::Execute<uint16_t, bool>(flags_vector, result, args.size(),
	                                       [&](uint16_t flags) { return (flags & 0x8) != 0; });
}

static void AlignmentIsReverseFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &flags_vector = args.data[0];
	UnaryExecutor::Execute<uint16_t, bool>(flags_vector, result, args.size(),
	                                       [&](uint16_t flags) { return (flags & 0x10) != 0; });
}

static void AlignmentIsMateReverseFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &flags_vector = args.data[0];
	UnaryExecutor::Execute<uint16_t, bool>(flags_vector, result, args.size(),
	                                       [&](uint16_t flags) { return (flags & 0x20) != 0; });
}

static void AlignmentIsRead1Function(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &flags_vector = args.data[0];
	UnaryExecutor::Execute<uint16_t, bool>(flags_vector, result, args.size(),
	                                       [&](uint16_t flags) { return (flags & 0x40) != 0; });
}

static void AlignmentIsRead2Function(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &flags_vector = args.data[0];
	UnaryExecutor::Execute<uint16_t, bool>(flags_vector, result, args.size(),
	                                       [&](uint16_t flags) { return (flags & 0x80) != 0; });
}

static void AlignmentIsSecondaryFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &flags_vector = args.data[0];
	UnaryExecutor::Execute<uint16_t, bool>(flags_vector, result, args.size(),
	                                       [&](uint16_t flags) { return (flags & 0x100) != 0; });
}

static void AlignmentIsPrimaryFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &flags_vector = args.data[0];
	UnaryExecutor::Execute<uint16_t, bool>(flags_vector, result, args.size(), [&](uint16_t flags) {
		return (flags & 0x100) == 0 && (flags & 0x800) == 0;
	});
}

static void AlignmentIsQcFailedFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &flags_vector = args.data[0];
	UnaryExecutor::Execute<uint16_t, bool>(flags_vector, result, args.size(),
	                                       [&](uint16_t flags) { return (flags & 0x200) != 0; });
}

static void AlignmentIsDuplicateFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &flags_vector = args.data[0];
	UnaryExecutor::Execute<uint16_t, bool>(flags_vector, result, args.size(),
	                                       [&](uint16_t flags) { return (flags & 0x400) != 0; });
}

static void AlignmentIsSupplementaryFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &flags_vector = args.data[0];
	UnaryExecutor::Execute<uint16_t, bool>(flags_vector, result, args.size(),
	                                       [&](uint16_t flags) { return (flags & 0x800) != 0; });
}

namespace {

// Local convenience: every alignment-flag function has the same shape
// (USMALLINT flags -> BOOLEAN), so wrap the helper with a one-liner.
void RegisterFlag(ExtensionLoader &loader, const std::string &name, scalar_function_t fn,
                  const std::string &description, const std::string &alias_of = "") {
	const std::string filter_example =
	    "-- Filter pattern: combine flag tests with boolean operators\n"
	    "SELECT read_id, flags\n"
	    "FROM read_alignments('alignments.sam')\n"
	    "WHERE " + name + "(flags) AND NOT alignment_is_unmapped(flags);";
	const std::string single_example =
	    "SELECT read_id, flags, " + name + "(flags) AS flag\n"
	    "FROM read_alignments('alignments.sam') LIMIT 10;";
	RegisterDocumentedScalar(loader, ScalarFunction(name, {LogicalType::USMALLINT}, LogicalType::BOOLEAN, fn),
	                         description, {"flags"}, {single_example, filter_example}, alias_of, {"sam-flags"});
}

} // namespace

void AlignmentFlagFunctions::Register(ExtensionLoader &loader) {
	const std::string desc_paired = "Returns true if the alignment is part of a paired-end read (SAM flag 0x1).";
	RegisterFlag(loader, "alignment_is_paired", AlignmentIsPairedFunction, desc_paired);
	RegisterFlag(loader, "is_paired", AlignmentIsPairedFunction, desc_paired, "alignment_is_paired");

	const std::string desc_proper = "Returns true if both mates in a paired alignment are mapped in the expected "
	                                "orientation and within the expected insert size (SAM flag 0x2).";
	RegisterFlag(loader, "alignment_is_proper_pair", AlignmentIsProperPairFunction, desc_proper);
	RegisterFlag(loader, "is_proper_pair", AlignmentIsProperPairFunction, desc_proper, "alignment_is_proper_pair");

	const std::string desc_unmapped = "Returns true if the read is unmapped (SAM flag 0x4).";
	RegisterFlag(loader, "alignment_is_unmapped", AlignmentIsUnmappedFunction, desc_unmapped);
	RegisterFlag(loader, "is_unmapped", AlignmentIsUnmappedFunction, desc_unmapped, "alignment_is_unmapped");

	const std::string desc_mate_unmapped = "Returns true if the mate of this read is unmapped (SAM flag 0x8).";
	RegisterFlag(loader, "alignment_is_mate_unmapped", AlignmentIsMateUnmappedFunction, desc_mate_unmapped);
	RegisterFlag(loader, "is_munmap", AlignmentIsMateUnmappedFunction, desc_mate_unmapped, "alignment_is_mate_unmapped");

	const std::string desc_reverse = "Returns true if the read is mapped to the reverse strand (SAM flag 0x10).";
	RegisterFlag(loader, "alignment_is_reverse", AlignmentIsReverseFunction, desc_reverse);
	RegisterFlag(loader, "is_reverse", AlignmentIsReverseFunction, desc_reverse, "alignment_is_reverse");

	const std::string desc_mate_reverse =
	    "Returns true if the mate of this read is mapped to the reverse strand (SAM flag 0x20).";
	RegisterFlag(loader, "alignment_is_mate_reverse", AlignmentIsMateReverseFunction, desc_mate_reverse);
	RegisterFlag(loader, "is_mreverse", AlignmentIsMateReverseFunction, desc_mate_reverse, "alignment_is_mate_reverse");

	const std::string desc_read1 = "Returns true if this is the first read in a pair (SAM flag 0x40).";
	RegisterFlag(loader, "alignment_is_read1", AlignmentIsRead1Function, desc_read1);
	RegisterFlag(loader, "is_read1", AlignmentIsRead1Function, desc_read1, "alignment_is_read1");

	const std::string desc_read2 = "Returns true if this is the second read in a pair (SAM flag 0x80).";
	RegisterFlag(loader, "alignment_is_read2", AlignmentIsRead2Function, desc_read2);
	RegisterFlag(loader, "is_read2", AlignmentIsRead2Function, desc_read2, "alignment_is_read2");

	const std::string desc_secondary =
	    "Returns true if this alignment is a secondary alignment (SAM flag 0x100). A read with multiple possible "
	    "mappings has one primary and zero or more secondary alignments.";
	RegisterFlag(loader, "alignment_is_secondary", AlignmentIsSecondaryFunction, desc_secondary);
	RegisterFlag(loader, "is_secondary", AlignmentIsSecondaryFunction, desc_secondary, "alignment_is_secondary");

	// Primary is the negation of secondary OR supplementary; no SAM flag of its own.
	RegisterFlag(loader, "alignment_is_primary", AlignmentIsPrimaryFunction,
	             "Returns true if this is the primary alignment for the read — i.e. neither secondary "
	             "(SAM flag 0x100) nor supplementary (SAM flag 0x800).");

	const std::string desc_qc = "Returns true if the read failed quality checks (SAM flag 0x200).";
	RegisterFlag(loader, "alignment_is_qc_failed", AlignmentIsQcFailedFunction, desc_qc);
	RegisterFlag(loader, "is_qcfail", AlignmentIsQcFailedFunction, desc_qc, "alignment_is_qc_failed");

	const std::string desc_dup =
	    "Returns true if the read is a PCR or optical duplicate, as marked by an upstream tool (SAM flag 0x400).";
	RegisterFlag(loader, "alignment_is_duplicate", AlignmentIsDuplicateFunction, desc_dup);
	RegisterFlag(loader, "is_dup", AlignmentIsDuplicateFunction, desc_dup, "alignment_is_duplicate");

	const std::string desc_supp = "Returns true if this is a supplementary (chimeric) alignment (SAM flag 0x800). "
	                              "Supplementary alignments are non-linear parts of a chimeric alignment.";
	RegisterFlag(loader, "alignment_is_supplementary", AlignmentIsSupplementaryFunction, desc_supp);
	RegisterFlag(loader, "is_supplementary", AlignmentIsSupplementaryFunction, desc_supp, "alignment_is_supplementary");
}

} // namespace duckdb
