#include "align_bowtie2_daemon_common.hpp"

#include "gpl_boundary/arrow_ipc.hpp"

#include "nanoarrow/nanoarrow.h"

#include <cstring>

namespace duckdb {
namespace bt2_daemon {

namespace gb = ::duckdb::miint::gpl_boundary;

const char *const kOutputColumnNames[kNumOutputColumns] = {
    "read_id",        "flags",         "reference",       "position", "stop_position", "mapq",   "cigar",
    "mate_reference", "mate_position", "template_length", "tag_as",   "tag_xs",        "tag_ys", "tag_xn",
    "tag_xm",         "tag_xo",        "tag_xg",          "tag_nm",   "tag_yt",        "tag_md", "tag_sa"};

void PopulateOutputSchema(std::vector<std::string> &names, std::vector<LogicalType> &types) {
	names = {kOutputColumnNames, kOutputColumnNames + kNumOutputColumns};
	types = {LogicalType::VARCHAR,   // read_id
	         LogicalType::USMALLINT, // flags
	         LogicalType::VARCHAR,   // reference
	         LogicalType::BIGINT,    // position
	         LogicalType::BIGINT,    // stop_position
	         LogicalType::UTINYINT,  // mapq
	         LogicalType::VARCHAR,   // cigar
	         LogicalType::VARCHAR,   // mate_reference
	         LogicalType::BIGINT,    // mate_position
	         LogicalType::BIGINT,    // template_length
	         LogicalType::BIGINT,    // tag_as   — Int32 on wire, widened
	         LogicalType::BIGINT,    // tag_xs   — Int32 on wire, widened
	         LogicalType::BIGINT,    // tag_ys   — Int32 on wire, widened
	         LogicalType::BIGINT,    // tag_xn   — Int32 on wire, widened
	         LogicalType::BIGINT,    // tag_xm   — Int32 on wire, widened
	         LogicalType::BIGINT,    // tag_xo   — Int32 on wire, widened
	         LogicalType::BIGINT,    // tag_xg   — Int32 on wire, widened
	         LogicalType::BIGINT,    // tag_nm   — Int32 on wire, widened
	         LogicalType::VARCHAR,   // tag_yt
	         LogicalType::VARCHAR,   // tag_md
	         LogicalType::VARCHAR};  // tag_sa
}

void ValidateOutputSchema(const ArrowSchema &schema) {
	if (schema.n_children != kNumOutputColumns) {
		throw IOException("align_bowtie2: daemon returned unexpected schema (%lld columns, expected %d)",
		                  static_cast<long long>(schema.n_children), kNumOutputColumns);
	}
	for (int c = 0; c < kNumOutputColumns; ++c) {
		const auto *child = schema.children[c];
		if (!child || !child->name || std::strcmp(child->name, kOutputColumnNames[c]) != 0) {
			throw IOException("align_bowtie2: schema drift at column %d — expected '%s', got '%s'", c,
			                  kOutputColumnNames[c], (child && child->name) ? child->name : "(null)");
		}
	}
}

// =============================================================================
// BuildQueryIpc — nanoarrow encoder for {read_id, sequence1, sequence2?,
// qual1?, qual2?}. Matches gpl-boundary's bowtie2-align input schema.
// =============================================================================

std::vector<uint8_t> BuildQueryIpc(const QueryBatch &qb, const QueryArrowSchema &schema_flags) {
	const int n_cols = schema_flags.num_columns();
	ArrowSchema schema {};
	auto rc = ArrowSchemaInitFromType(&schema, NANOARROW_TYPE_STRUCT);
	if (rc != NANOARROW_OK) {
		throw InternalException("align_bowtie2: ArrowSchemaInit failed");
	}
	rc = ArrowSchemaAllocateChildren(&schema, n_cols);
	if (rc != NANOARROW_OK) {
		schema.release(&schema);
		throw InternalException("align_bowtie2: ArrowSchemaAllocateChildren failed");
	}
	int idx = 0;
	ArrowSchemaInitFromType(schema.children[idx], NANOARROW_TYPE_STRING);
	ArrowSchemaSetName(schema.children[idx], "read_id");
	++idx;
	ArrowSchemaInitFromType(schema.children[idx], NANOARROW_TYPE_STRING);
	ArrowSchemaSetName(schema.children[idx], "sequence1");
	++idx;
	if (schema_flags.has_sequence2) {
		ArrowSchemaInitFromType(schema.children[idx], NANOARROW_TYPE_STRING);
		ArrowSchemaSetName(schema.children[idx], "sequence2");
		++idx;
	}
	if (schema_flags.has_qual1) {
		ArrowSchemaInitFromType(schema.children[idx], NANOARROW_TYPE_STRING);
		ArrowSchemaSetName(schema.children[idx], "qual1");
		++idx;
	}
	if (schema_flags.has_qual2) {
		ArrowSchemaInitFromType(schema.children[idx], NANOARROW_TYPE_STRING);
		ArrowSchemaSetName(schema.children[idx], "qual2");
		++idx;
	}

	ArrowArray array {};
	ArrowError err {};
	if (ArrowArrayInitFromSchema(&array, &schema, &err) != NANOARROW_OK) {
		schema.release(&schema);
		throw InternalException("align_bowtie2: ArrowArrayInit failed: %s", err.message);
	}
	if (ArrowArrayStartAppending(&array) != NANOARROW_OK) {
		array.release(&array);
		schema.release(&schema);
		throw InternalException("align_bowtie2: ArrowArrayStartAppending failed");
	}

	const idx_t n_rows = qb.read_ids.size();
	for (idx_t row = 0; row < n_rows; ++row) {
		int c = 0;
		ArrowStringView rv {qb.read_ids[row].data(), static_cast<int64_t>(qb.read_ids[row].size())};
		if (ArrowArrayAppendString(array.children[c++], rv) != NANOARROW_OK) {
			array.release(&array);
			schema.release(&schema);
			throw InternalException("align_bowtie2: read_id append failed at row %lld", static_cast<long long>(row));
		}
		ArrowStringView s1v {qb.sequence1[row].data(), static_cast<int64_t>(qb.sequence1[row].size())};
		if (ArrowArrayAppendString(array.children[c++], s1v) != NANOARROW_OK) {
			array.release(&array);
			schema.release(&schema);
			throw InternalException("align_bowtie2: sequence1 append failed at row %lld", static_cast<long long>(row));
		}
		if (schema_flags.has_sequence2) {
			ArrowArray *child = array.children[c++];
			if (row < qb.sequence2_valid.size() && qb.sequence2_valid[row]) {
				ArrowStringView v {qb.sequence2[row].data(), static_cast<int64_t>(qb.sequence2[row].size())};
				if (ArrowArrayAppendString(child, v) != NANOARROW_OK) {
					array.release(&array);
					schema.release(&schema);
					throw InternalException("align_bowtie2: sequence2 append failed");
				}
			} else {
				if (ArrowArrayAppendNull(child, 1) != NANOARROW_OK) {
					array.release(&array);
					schema.release(&schema);
					throw InternalException("align_bowtie2: sequence2 null append failed");
				}
			}
		}
		if (schema_flags.has_qual1) {
			ArrowArray *child = array.children[c++];
			if (row < qb.qual1_valid.size() && qb.qual1_valid[row]) {
				ArrowStringView v {qb.qual1[row].data(), static_cast<int64_t>(qb.qual1[row].size())};
				if (ArrowArrayAppendString(child, v) != NANOARROW_OK) {
					array.release(&array);
					schema.release(&schema);
					throw InternalException("align_bowtie2: qual1 append failed");
				}
			} else {
				if (ArrowArrayAppendNull(child, 1) != NANOARROW_OK) {
					array.release(&array);
					schema.release(&schema);
					throw InternalException("align_bowtie2: qual1 null append failed");
				}
			}
		}
		if (schema_flags.has_qual2) {
			ArrowArray *child = array.children[c++];
			if (row < qb.qual2_valid.size() && qb.qual2_valid[row]) {
				ArrowStringView v {qb.qual2[row].data(), static_cast<int64_t>(qb.qual2[row].size())};
				if (ArrowArrayAppendString(child, v) != NANOARROW_OK) {
					array.release(&array);
					schema.release(&schema);
					throw InternalException("align_bowtie2: qual2 append failed");
				}
			} else {
				if (ArrowArrayAppendNull(child, 1) != NANOARROW_OK) {
					array.release(&array);
					schema.release(&schema);
					throw InternalException("align_bowtie2: qual2 null append failed");
				}
			}
		}
		if (ArrowArrayFinishElement(&array) != NANOARROW_OK) {
			array.release(&array);
			schema.release(&schema);
			throw InternalException("align_bowtie2: FinishElement failed at row %lld", static_cast<long long>(row));
		}
	}
	if (ArrowArrayFinishBuildingDefault(&array, &err) != NANOARROW_OK) {
		array.release(&array);
		schema.release(&schema);
		throw InternalException("align_bowtie2: ArrowArrayFinishBuilding failed: %s", err.message);
	}

	std::vector<uint8_t> bytes;
	try {
		bytes = gb::EncodeIpcStream(&schema, &array, 1);
	} catch (...) {
		array.release(&array);
		schema.release(&schema);
		throw;
	}
	if (array.release) {
		array.release(&array);
	}
	if (schema.release) {
		schema.release(&schema);
	}
	return bytes;
}

// =============================================================================
// EmitChunkRows — extract one DuckDB chunk's worth of rows from a decoded
// Arrow batch. Matches phylogeny_fasttree.cpp:711-737 decoder pattern.
// =============================================================================

namespace {

template <typename T>
T read_fixed(const ArrowArray &col, idx_t logical_index) {
	const auto *buf = reinterpret_cast<const T *>(col.buffers[1]);
	return buf[static_cast<idx_t>(col.offset) + logical_index];
}

bool is_null_at(const ArrowArray &arr, idx_t logical_index) {
	if (!arr.buffers[0]) {
		return false;
	}
	if (arr.null_count == 0) {
		return false;
	}
	const auto *bitmap = static_cast<const uint8_t *>(arr.buffers[0]);
	const idx_t abs = static_cast<idx_t>(arr.offset) + logical_index;
	return (bitmap[abs / 8] & (1u << (abs % 8))) == 0;
}

void emit_string(Vector &out, idx_t out_row, const ArrowArray &col, idx_t logical_index) {
	const auto *offsets = static_cast<const int32_t *>(col.buffers[1]);
	const auto *data = static_cast<const char *>(col.buffers[2]);
	const idx_t a = static_cast<idx_t>(col.offset) + logical_index;
	const int32_t start = offsets[a];
	const int32_t end = offsets[a + 1];
	const int32_t len = end - start;
	if (len < 0) {
		throw IOException("align_bowtie2: corrupt utf8 offsets at row %lld (start=%d end=%d)",
		                  static_cast<long long>(a), start, end);
	}
	FlatVector::GetData<string_t>(out)[out_row] = StringVector::AddString(out, data + start, static_cast<idx_t>(len));
}

} // namespace

void EmitChunkRows(DataChunk &output, idx_t to_emit, idx_t row_start, const ArrowArray &batch) {
	auto &v_read_id = output.data[0];
	auto *out_flags = FlatVector::GetData<uint16_t>(output.data[1]);
	auto &v_reference = output.data[2];
	auto *out_position = FlatVector::GetData<int64_t>(output.data[3]);
	auto *out_stop = FlatVector::GetData<int64_t>(output.data[4]);
	auto *out_mapq = FlatVector::GetData<uint8_t>(output.data[5]);
	auto &v_cigar = output.data[6];
	auto &v_mate_ref = output.data[7];
	auto *out_mate_pos = FlatVector::GetData<int64_t>(output.data[8]);
	auto *out_tlen = FlatVector::GetData<int64_t>(output.data[9]);
	auto *out_tag_as = FlatVector::GetData<int64_t>(output.data[10]);
	auto *out_tag_xs = FlatVector::GetData<int64_t>(output.data[11]);
	auto *out_tag_ys = FlatVector::GetData<int64_t>(output.data[12]);
	auto *out_tag_xn = FlatVector::GetData<int64_t>(output.data[13]);
	auto *out_tag_xm = FlatVector::GetData<int64_t>(output.data[14]);
	auto *out_tag_xo = FlatVector::GetData<int64_t>(output.data[15]);
	auto *out_tag_xg = FlatVector::GetData<int64_t>(output.data[16]);
	auto *out_tag_nm = FlatVector::GetData<int64_t>(output.data[17]);
	auto &v_tag_yt = output.data[18];
	auto &v_tag_md = output.data[19];
	auto &v_tag_sa = output.data[20];

	auto &mask_tag_as = FlatVector::Validity(output.data[10]);
	auto &mask_tag_xs = FlatVector::Validity(output.data[11]);
	auto &mask_tag_ys = FlatVector::Validity(output.data[12]);
	auto &mask_tag_xn = FlatVector::Validity(output.data[13]);
	auto &mask_tag_xm = FlatVector::Validity(output.data[14]);
	auto &mask_tag_xo = FlatVector::Validity(output.data[15]);
	auto &mask_tag_xg = FlatVector::Validity(output.data[16]);
	auto &mask_tag_nm = FlatVector::Validity(output.data[17]);
	auto &mask_tag_yt = FlatVector::Validity(v_tag_yt);
	auto &mask_tag_md = FlatVector::Validity(v_tag_md);
	auto &mask_tag_sa = FlatVector::Validity(v_tag_sa);

	const auto &col_read_id = *batch.children[0];
	const auto &col_flags = *batch.children[1];
	const auto &col_reference = *batch.children[2];
	const auto &col_position = *batch.children[3];
	const auto &col_stop = *batch.children[4];
	const auto &col_mapq = *batch.children[5];
	const auto &col_cigar = *batch.children[6];
	const auto &col_mate_ref = *batch.children[7];
	const auto &col_mate_pos = *batch.children[8];
	const auto &col_tlen = *batch.children[9];
	const auto &col_tag_as = *batch.children[10];
	const auto &col_tag_xs = *batch.children[11];
	const auto &col_tag_ys = *batch.children[12];
	const auto &col_tag_xn = *batch.children[13];
	const auto &col_tag_xm = *batch.children[14];
	const auto &col_tag_xo = *batch.children[15];
	const auto &col_tag_xg = *batch.children[16];
	const auto &col_tag_nm = *batch.children[17];
	const auto &col_tag_yt = *batch.children[18];
	const auto &col_tag_md = *batch.children[19];
	const auto &col_tag_sa = *batch.children[20];

	for (idx_t i = 0; i < to_emit; ++i) {
		const idx_t li = row_start + i;

		emit_string(v_read_id, i, col_read_id, li);
		out_flags[i] = read_fixed<uint16_t>(col_flags, li);
		emit_string(v_reference, i, col_reference, li);
		out_position[i] = read_fixed<int64_t>(col_position, li);
		out_stop[i] = read_fixed<int64_t>(col_stop, li);
		out_mapq[i] = read_fixed<uint8_t>(col_mapq, li);
		emit_string(v_cigar, i, col_cigar, li);
		emit_string(v_mate_ref, i, col_mate_ref, li);
		out_mate_pos[i] = read_fixed<int64_t>(col_mate_pos, li);
		out_tlen[i] = read_fixed<int64_t>(col_tlen, li);

		auto widen = [&](const ArrowArray &col, int64_t *out, ValidityMask &mask) {
			if (is_null_at(col, li)) {
				mask.SetInvalid(i);
			} else {
				out[i] = static_cast<int64_t>(read_fixed<int32_t>(col, li));
			}
		};
		widen(col_tag_as, out_tag_as, mask_tag_as);
		widen(col_tag_xs, out_tag_xs, mask_tag_xs);
		widen(col_tag_ys, out_tag_ys, mask_tag_ys);
		widen(col_tag_xn, out_tag_xn, mask_tag_xn);
		widen(col_tag_xm, out_tag_xm, mask_tag_xm);
		widen(col_tag_xo, out_tag_xo, mask_tag_xo);
		widen(col_tag_xg, out_tag_xg, mask_tag_xg);
		widen(col_tag_nm, out_tag_nm, mask_tag_nm);

		auto emit_nullable_str = [&](Vector &v, ValidityMask &mask, const ArrowArray &col) {
			if (is_null_at(col, li)) {
				mask.SetInvalid(i);
			} else {
				emit_string(v, i, col, li);
			}
		};
		emit_nullable_str(v_tag_yt, mask_tag_yt, col_tag_yt);
		emit_nullable_str(v_tag_md, mask_tag_md, col_tag_md);
		emit_nullable_str(v_tag_sa, mask_tag_sa, col_tag_sa);
	}
}

} // namespace bt2_daemon
} // namespace duckdb
