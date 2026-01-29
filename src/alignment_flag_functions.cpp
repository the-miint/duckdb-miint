#include "alignment_flag_functions.hpp"
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

void AlignmentFlagFunctions::Register(ExtensionLoader &loader) {
	ScalarFunction alignment_is_paired("alignment_is_paired", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                                   AlignmentIsPairedFunction);
	loader.RegisterFunction(alignment_is_paired);
	ScalarFunction is_paired("is_paired", {LogicalType::USMALLINT}, LogicalType::BOOLEAN, AlignmentIsPairedFunction);
	loader.RegisterFunction(is_paired);

	ScalarFunction alignment_is_proper_pair("alignment_is_proper_pair", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                                        AlignmentIsProperPairFunction);
	loader.RegisterFunction(alignment_is_proper_pair);
	ScalarFunction is_proper_pair("is_proper_pair", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                              AlignmentIsProperPairFunction);
	loader.RegisterFunction(is_proper_pair);

	ScalarFunction alignment_is_unmapped("alignment_is_unmapped", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                                     AlignmentIsUnmappedFunction);
	loader.RegisterFunction(alignment_is_unmapped);
	ScalarFunction is_unmapped("is_unmapped", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                           AlignmentIsUnmappedFunction);
	loader.RegisterFunction(is_unmapped);

	ScalarFunction alignment_is_mate_unmapped("alignment_is_mate_unmapped", {LogicalType::USMALLINT},
	                                          LogicalType::BOOLEAN, AlignmentIsMateUnmappedFunction);
	loader.RegisterFunction(alignment_is_mate_unmapped);
	ScalarFunction is_munmap("is_munmap", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                         AlignmentIsMateUnmappedFunction);
	loader.RegisterFunction(is_munmap);

	ScalarFunction alignment_is_reverse("alignment_is_reverse", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                                    AlignmentIsReverseFunction);
	loader.RegisterFunction(alignment_is_reverse);
	ScalarFunction is_reverse("is_reverse", {LogicalType::USMALLINT}, LogicalType::BOOLEAN, AlignmentIsReverseFunction);
	loader.RegisterFunction(is_reverse);

	ScalarFunction alignment_is_mate_reverse("alignment_is_mate_reverse", {LogicalType::USMALLINT},
	                                         LogicalType::BOOLEAN, AlignmentIsMateReverseFunction);
	loader.RegisterFunction(alignment_is_mate_reverse);
	ScalarFunction is_mreverse("is_mreverse", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                           AlignmentIsMateReverseFunction);
	loader.RegisterFunction(is_mreverse);

	ScalarFunction alignment_is_read1("alignment_is_read1", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                                  AlignmentIsRead1Function);
	loader.RegisterFunction(alignment_is_read1);
	ScalarFunction is_read1("is_read1", {LogicalType::USMALLINT}, LogicalType::BOOLEAN, AlignmentIsRead1Function);
	loader.RegisterFunction(is_read1);

	ScalarFunction alignment_is_read2("alignment_is_read2", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                                  AlignmentIsRead2Function);
	loader.RegisterFunction(alignment_is_read2);
	ScalarFunction is_read2("is_read2", {LogicalType::USMALLINT}, LogicalType::BOOLEAN, AlignmentIsRead2Function);
	loader.RegisterFunction(is_read2);

	ScalarFunction alignment_is_secondary("alignment_is_secondary", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                                      AlignmentIsSecondaryFunction);
	loader.RegisterFunction(alignment_is_secondary);

	ScalarFunction alignment_is_primary("alignment_is_primary", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                                    AlignmentIsPrimaryFunction);
	loader.RegisterFunction(alignment_is_primary);

	ScalarFunction is_secondary("is_secondary", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                            AlignmentIsSecondaryFunction);
	loader.RegisterFunction(is_secondary);

	ScalarFunction alignment_is_qc_failed("alignment_is_qc_failed", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                                      AlignmentIsQcFailedFunction);
	loader.RegisterFunction(alignment_is_qc_failed);
	ScalarFunction is_qcfail("is_qcfail", {LogicalType::USMALLINT}, LogicalType::BOOLEAN, AlignmentIsQcFailedFunction);
	loader.RegisterFunction(is_qcfail);

	ScalarFunction alignment_is_duplicate("alignment_is_duplicate", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                                      AlignmentIsDuplicateFunction);
	loader.RegisterFunction(alignment_is_duplicate);
	ScalarFunction is_dup("is_dup", {LogicalType::USMALLINT}, LogicalType::BOOLEAN, AlignmentIsDuplicateFunction);
	loader.RegisterFunction(is_dup);

	ScalarFunction alignment_is_supplementary("alignment_is_supplementary", {LogicalType::USMALLINT},
	                                          LogicalType::BOOLEAN, AlignmentIsSupplementaryFunction);
	loader.RegisterFunction(alignment_is_supplementary);
	ScalarFunction is_supplementary("is_supplementary", {LogicalType::USMALLINT}, LogicalType::BOOLEAN,
	                                AlignmentIsSupplementaryFunction);
	loader.RegisterFunction(is_supplementary);
}

} // namespace duckdb
