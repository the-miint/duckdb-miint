#pragma once

#include "SortMeRNAAligner.hpp"
#include "align_common.hpp"
#include "align_result_utils.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include <cstdint>
#include <string>
#include <vector>

namespace duckdb {

// Output schema for align_sortmerna_rrna — carries fields SAM cannot: e_value,
// identity, coverage, edit_distance, segment_idx. Paired-end produces two rows
// per input row with segment_idx 0 (fwd) and 1 (rev).
inline std::vector<std::string> GetSortMeRNARRNAOutputNames() {
	return {"read_id", "aligned", "strand",   "ref_name", "ref_start",     "ref_end",    "cigar",
	        "score",   "e_value", "identity", "coverage", "edit_distance", "segment_idx"};
}

// `query_id_type` drives `read_id`. Must be VARCHAR or BIGINT. No default —
// every caller must commit to a type so the audit can't slip. The other 12
// columns have fixed types (ref_name is free-form FASTA header text, never
// an id column).
inline std::vector<LogicalType> GetSortMeRNARRNAOutputTypes(const LogicalType &query_id_type) {
	return {query_id_type,        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::VARCHAR,
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::INTEGER,
	        LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::INTEGER,
	        LogicalType::INTEGER};
}

// Emit [offset, offset+count) rows of `batch` into the output DataChunk.
// `query_id_type` (VARCHAR or BIGINT) drives the read_id column; the other
// 12 columns have fixed types matching GetSortMeRNARRNAOutputTypes.
inline void OutputSortMeRNARRNABatch(DataChunk &output, const miint::SortMeRNAResultBatch &batch, idx_t offset,
                                     idx_t count, const LogicalType &query_id_type) {
	D_ASSERT(IsAllowedIdType(query_id_type));
	idx_t col = 0;
	EmitIdColumnFromStrings(output.data[col++], batch.read_ids, offset, count, query_id_type);
	SetAlignResultInt32(output.data[col++], batch.aligned, offset, count);
	SetAlignResultInt32(output.data[col++], batch.strands, offset, count);
	SetAlignResultString(output.data[col++], batch.ref_names, offset, count);
	SetAlignResultInt32(output.data[col++], batch.ref_starts, offset, count);
	SetAlignResultInt32(output.data[col++], batch.ref_ends, offset, count);
	SetAlignResultString(output.data[col++], batch.cigars, offset, count);
	SetAlignResultInt32(output.data[col++], batch.scores, offset, count);
	SetAlignResultDouble(output.data[col++], batch.e_values, offset, count);
	SetAlignResultDouble(output.data[col++], batch.identities, offset, count);
	SetAlignResultDouble(output.data[col++], batch.coverages, offset, count);
	SetAlignResultInt32(output.data[col++], batch.edit_distances, offset, count);
	SetAlignResultInt32(output.data[col++], batch.segment_indices, offset, count);
	output.SetCardinality(count);
}

// Project a SortMeRNAResultBatch slice onto the shared SAM output schema
// (see align_common.hpp::GetAlignmentOutputNames). Sortmerna does not compute
// MAPQ/MD/SA/XS/YT; those columns are NULL or set to conventional defaults
// (mapq=255 "unavailable"). Paired rows arrive consecutively from
// SortMeRNAAligner::align() as (segment=0, segment=1), so the mate fields for
// row j can be read from batch[j ^ 1] without cross-chunk lookups.
//
// `query_id_type` drives `read_id`; `subject_id_type` drives `reference` and
// `mate_reference`. Both must be VARCHAR, BIGINT, or UUID. The parameter mirrors the
// shape of OutputSAMRecordBatch (align_common.hpp) and lets the entry-point
// asserts catch a forgetful caller, but for sortmerna `subject_id_type` is
// always VARCHAR — subject names come from FASTA files on disk via
// `ref_paths`, never a user-provided table. The body therefore emits the SAM
// `=` mate-reference sentinel verbatim (VARCHAR carries `=` in-band); there
// is no BIGINT-subject resolution path because there is no BIGINT subject.
inline void OutputSortMeRNASamBatch(DataChunk &output, const miint::SortMeRNAResultBatch &batch, idx_t offset,
                                    idx_t count, bool is_paired, const LogicalType &query_id_type,
                                    const LogicalType &subject_id_type) {
	D_ASSERT(IsAllowedIdType(query_id_type));
	D_ASSERT(IsAllowedIdType(subject_id_type));

	// Pair-consecutive invariant: SortMeRNAAligner::align() emits segments as
	// (seg=0, seg=1, seg=0, seg=1, ...). Mate lookup below uses `j ^ 1`, which
	// requires both `offset` and `count` to be even so no pair straddles a
	// chunk boundary or runs off the end of the batch. Execute's
	// count = min(available, STANDARD_VECTOR_SIZE) preserves both invariants
	// as long as the total result-buffer size is even, which is guaranteed
	// for paired-end input. Assert here rather than relying on the caller alone.
	D_ASSERT(!is_paired || (offset % 2 == 0));
	D_ASSERT(!is_paired || (count % 2 == 0));

	idx_t col = 0;
	auto &read_id_vec = output.data[col++];
	auto &flags_vec = output.data[col++];
	auto &reference_vec = output.data[col++];
	auto &position_vec = output.data[col++];
	auto &stop_position_vec = output.data[col++];
	auto &mapq_vec = output.data[col++];
	auto &cigar_vec = output.data[col++];
	auto &mate_reference_vec = output.data[col++];
	auto &mate_position_vec = output.data[col++];
	auto &template_length_vec = output.data[col++];
	auto &tag_as_vec = output.data[col++];
	auto &tag_xs_vec = output.data[col++];
	auto &tag_ys_vec = output.data[col++];
	auto &tag_xn_vec = output.data[col++];
	auto &tag_xm_vec = output.data[col++];
	auto &tag_xo_vec = output.data[col++];
	auto &tag_xg_vec = output.data[col++];
	auto &tag_nm_vec = output.data[col++];
	auto &tag_yt_vec = output.data[col++];
	auto &tag_md_vec = output.data[col++];
	auto &tag_sa_vec = output.data[col++];

	auto flags_data = FlatVector::GetData<uint16_t>(flags_vec);
	auto position_data = FlatVector::GetData<int64_t>(position_vec);
	auto stop_position_data = FlatVector::GetData<int64_t>(stop_position_vec);
	auto mapq_data = FlatVector::GetData<uint8_t>(mapq_vec);
	auto cigar_data = FlatVector::GetData<string_t>(cigar_vec);
	auto mate_position_data = FlatVector::GetData<int64_t>(mate_position_vec);
	auto template_length_data = FlatVector::GetData<int64_t>(template_length_vec);
	auto tag_as_data = FlatVector::GetData<int64_t>(tag_as_vec);
	auto tag_nm_data = FlatVector::GetData<int64_t>(tag_nm_vec);

	// Match the SetAlignResult*Nullable pattern in align_result_utils.hpp:
	// mark validity explicitly before per-row writes rather than trusting the
	// DataChunk's recycled state. tag_as / tag_nm flip to valid on aligned
	// rows; the remaining nine tags are always NULL for sortmerna output.
	auto &tag_as_validity = FlatVector::Validity(tag_as_vec);
	auto &tag_nm_validity = FlatVector::Validity(tag_nm_vec);
	tag_as_validity.SetAllInvalid(count);
	tag_nm_validity.SetAllInvalid(count);
	FlatVector::Validity(tag_xs_vec).SetAllInvalid(count);
	FlatVector::Validity(tag_ys_vec).SetAllInvalid(count);
	FlatVector::Validity(tag_xn_vec).SetAllInvalid(count);
	FlatVector::Validity(tag_xm_vec).SetAllInvalid(count);
	FlatVector::Validity(tag_xo_vec).SetAllInvalid(count);
	FlatVector::Validity(tag_xg_vec).SetAllInvalid(count);
	FlatVector::Validity(tag_yt_vec).SetAllInvalid(count);
	FlatVector::Validity(tag_md_vec).SetAllInvalid(count);
	FlatVector::Validity(tag_sa_vec).SetAllInvalid(count);

	// Stage id-typed columns into vector<string> slices and bulk-emit via
	// EmitIdColumnFromStrings after the row loop — supports VARCHAR and BIGINT
	// uniformly. The "*" sentinel for unaligned rows encodes NULL through the
	// BIGINT codec (read_id only — subject side is always VARCHAR here).
	std::vector<std::string> emit_read_ids(count);
	std::vector<std::string> emit_references(count);
	std::vector<std::string> emit_mate_references(count);

	for (idx_t i = 0; i < count; ++i) {
		const idx_t j = offset + i;
		const bool aligned = batch.aligned[j] != 0;
		const bool reverse_strand = batch.strands[j] == 0;

		emit_read_ids[i] = batch.read_ids[j];

		uint16_t flags = 0;
		if (is_paired) {
			flags |= 0x1; // paired
			// segment_idx comes from SortMeRNAAligner::align: 0 = fwd, 1 = rev.
			if (batch.segment_indices[j] == 0) {
				flags |= 0x40; // first in pair
			} else {
				flags |= 0x80; // second in pair
			}
			const idx_t mate_j = j ^ 1;
			const bool mate_aligned = batch.aligned[mate_j] != 0;
			const bool mate_reverse = batch.strands[mate_j] == 0;
			// SAM spec 1.4: 0x2 is "each segment properly aligned according to
			// the aligner." Sortmerna is an rRNA classifier, not a paired-end
			// genome aligner — it aligns mates independently. We set 0x2 when
			// both mates aligned (to the same reference below), deliberately
			// weaker than the standard "concordant orientation within expected
			// insert size" meaning. Downstream consumers reading the bit as
			// "concordant" will be wrong, which is an intrinsic mismatch
			// between sortmerna's output and SAM's paired-end semantics.
			if (aligned && mate_aligned) {
				flags |= 0x2;
			}
			if (!mate_aligned) {
				flags |= 0x8; // mate unmapped
			}
			if (mate_aligned && mate_reverse) {
				flags |= 0x20; // mate reverse strand
			}
		}
		if (!aligned) {
			flags |= 0x4; // unmapped
		} else if (reverse_strand) {
			flags |= 0x10; // reverse strand
		}
		flags_data[i] = flags;

		if (aligned) {
			emit_references[i] = batch.ref_names[j];
			position_data[i] = batch.ref_starts[j];
			// sortmerna's ref_end is 1-based inclusive; the SAM schema uses
			// 1-based half-open stop_position (matches bowtie2/minimap2
			// conventions elsewhere in this codebase).
			stop_position_data[i] = batch.ref_ends[j] + 1;
			cigar_data[i] = StringVector::AddString(cigar_vec, batch.cigars[j]);
		} else {
			emit_references[i] = "*";
			position_data[i] = 0;
			stop_position_data[i] = 0;
			cigar_data[i] = StringVector::AddString(cigar_vec, "*");
		}

		// SAM convention for "MAPQ unavailable" — sortmerna does not compute
		// mapping quality, so we emit 255 for aligned and unaligned alike.
		mapq_data[i] = 255;

		if (is_paired) {
			const idx_t mate_j = j ^ 1;
			const bool mate_aligned = batch.aligned[mate_j] != 0;
			if (aligned && mate_aligned && batch.ref_names[j] == batch.ref_names[mate_j]) {
				// VARCHAR-only subject side (see helper docstring), so `=`
				// passes through verbatim — no BIGINT resolution needed.
				emit_mate_references[i] = "=";
				mate_position_data[i] = batch.ref_starts[mate_j];
			} else if (mate_aligned) {
				emit_mate_references[i] = batch.ref_names[mate_j];
				mate_position_data[i] = batch.ref_starts[mate_j];
			} else {
				emit_mate_references[i] = "*";
				mate_position_data[i] = 0;
			}
		} else {
			emit_mate_references[i] = "*";
			mate_position_data[i] = 0;
		}
		template_length_data[i] = 0;

		if (aligned) {
			tag_as_data[i] = batch.scores[j];
			tag_nm_data[i] = batch.edit_distances[j];
			tag_as_validity.SetValid(i);
			tag_nm_validity.SetValid(i);
		}
		// Sortmerna produces neither secondary-alignment scores nor bowtie2's
		// XS/XN/XM/XO/XG/YT/MD tags nor minimap2's SA chain; they were marked
		// NULL via SetAllInvalid(count) above.
	}

	EmitIdColumnFromStrings(read_id_vec, emit_read_ids, 0, count, query_id_type);
	EmitIdColumnFromStrings(reference_vec, emit_references, 0, count, subject_id_type);
	EmitIdColumnFromStrings(mate_reference_vec, emit_mate_references, 0, count, subject_id_type);

	output.SetCardinality(count);
}

} // namespace duckdb
