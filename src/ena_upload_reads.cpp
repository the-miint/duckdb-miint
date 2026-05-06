// SPDX-License-Identifier: MIT
//
// `ena_upload_reads` table function — see ena_upload_reads.hpp.

#include "ena_upload_reads.hpp"

#include "aspera_send.hpp"
#include "aspera_utils.hpp"
#ifdef MIINT_HAS_CURL
#include "curl_send.hpp"
#endif
#include "ena_upload_helpers.hpp"
#include "fastq_encoder.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/crypto/md5.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/main/secret/secret.hpp"
#include "duckdb/main/secret/secret_manager.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unistd.h> // unlink, rmdir, mkdtemp
#include <vector>
#include <zlib.h>

namespace duckdb {

namespace {

// =====================================================================
// Bind data
// =====================================================================
struct ENAUploadReadsBindData : public TableFunctionData {
	string relation_name;
	string secret_name; // empty → no secret (only valid for file://)
	string target_url;  // raw value as supplied by the user
	uint8_t qual_offset = 33;
	FastqLayoutMode layout_mode = FastqLayoutMode::AUTO;
	string aspera_max_rate; // empty → no -l flag
};

// =====================================================================
// Materialised input rows
// =====================================================================
struct UploadInputRows {
	vector<string> sample_ref;
	vector<string> read_id;
	vector<string> sequence1;
	vector<vector<uint8_t>> qual1;
	vector<bool> has_seq2;
	vector<string> sequence2;
	vector<vector<uint8_t>> qual2;
	vector<int64_t> sequence_index;
	vector<bool> sequence_index_valid;
	uint64_t total_seq_bytes = 0; // running total of sequence1 + sequence2 lengths
};

// One per (sample_ref) after grouping + layout resolution.
struct SampleGroup {
	string sample_ref;
	FastqLayoutMode resolved_layout;
	vector<size_t> row_indices; // already sorted by sequence_index NULLS LAST
};

// One per output FASTQ file (single sample → 1 or 2 files).
struct EmittedFile {
	string sample_ref;
	string filename;
	string filetype;
	string md5_hex;
	uint64_t bytes_written;
	string layout_name;
};

struct ENAUploadReadsGlobalState : public GlobalTableFunctionState {
	UploadInputRows input;
	vector<SampleGroup> samples;
	UploadTargetURL target;

	// Authenticated-transport fields (left empty when transport=LOCAL_FILE).
	// Both ASPERA and CURL paths consume `user` + `password`; the Aspera
	// path additionally needs the binary + key.
	string user;
	string password;
	string aspera_ascp_path;
	string aspera_key_path;
	string aspera_max_rate;

	vector<EmittedFile> emitted;
	idx_t emit_cursor = 0;

	idx_t MaxThreads() const override {
		return 1;
	}
};

// =====================================================================
// Input materialisation (separate connection, read-time)
// =====================================================================

// Find a column by case-insensitive name. Returns -1 when missing.
int FindColumn(const vector<string> &names, const char *target) {
	for (idx_t i = 0; i < names.size(); i++) {
		if (StringUtil::CIEquals(names[i], target)) {
			return static_cast<int>(i);
		}
	}
	return -1;
}

void RequireType(const string &col, const LogicalType &actual, LogicalTypeId expected) {
	if (actual.id() != expected) {
		throw InvalidInputException("ena_upload_reads: input column '%s' must be %s (got %s)", col,
		                            LogicalTypeIdToString(expected), actual.ToString());
	}
}

void RequireListUtinyint(const string &col, const LogicalType &actual) {
	if (actual.id() != LogicalTypeId::LIST) {
		throw InvalidInputException("ena_upload_reads: input column '%s' must be LIST(UTINYINT) (got %s)", col,
		                            actual.ToString());
	}
	auto &child = ListType::GetChildType(actual);
	if (child.id() != LogicalTypeId::UTINYINT) {
		throw InvalidInputException("ena_upload_reads: input column '%s' must be LIST(UTINYINT) (got LIST(%s))", col,
		                            child.ToString());
	}
}

vector<uint8_t> ExtractQualList(Vector &list_vec, UnifiedVectorFormat &list_data, idx_t row) {
	auto entries = UnifiedVectorFormat::GetData<list_entry_t>(list_data);
	auto sel_idx = list_data.sel->get_index(row);
	auto &child = ListVector::GetEntry(list_vec);
	auto child_data = FlatVector::GetData<uint8_t>(child);
	const idx_t len = entries[sel_idx].length;
	const idx_t offset = entries[sel_idx].offset;
	vector<uint8_t> out(len);
	std::memcpy(out.data(), child_data + offset, len);
	return out;
}

// Hard caps on the materialised input. We pull the entire relation into
// memory once (one std::string per sequence + one std::vector per quality
// list), which is fine for the typical "submit a few GB of FASTQ per batch"
// workflow but will OOM on truly large inputs. These bounds turn the silent
// OOM into an actionable error pointing the user at the right escape hatch
// (split-by-sample, or re-call per chunk). Streaming would remove the cap;
// it's deferred until a user actually trips this.
constexpr idx_t MAX_INPUT_ROWS = 50'000'000;
constexpr uint64_t MAX_INPUT_SEQUENCE_BYTES = 5ULL * 1024 * 1024 * 1024; // 5 GB

void MaterialiseInput(ClientContext &context, const string &relation_name, UploadInputRows &out) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const string query = "SELECT * FROM " + KeywordHelper::WriteOptionallyQuoted(relation_name);
	auto result = conn.Query(query);
	if (result->HasError()) {
		throw InvalidInputException("ena_upload_reads: failed to read relation '%s': %s", relation_name,
		                            result->GetError());
	}
	auto &materialized = result->Cast<MaterializedQueryResult>();
	auto &types = materialized.types;
	auto &names = materialized.names;

	const int sample_ref_idx = FindColumn(names, "sample_ref");
	const int read_id_idx = FindColumn(names, "read_id");
	const int sequence1_idx = FindColumn(names, "sequence1");
	const int qual1_idx = FindColumn(names, "qual1");
	const int sequence2_idx = FindColumn(names, "sequence2");
	const int qual2_idx = FindColumn(names, "qual2");
	const int sequence_index_idx = FindColumn(names, "sequence_index");

	if (sample_ref_idx < 0) {
		throw InvalidInputException("ena_upload_reads: input relation must include a 'sample_ref' column");
	}
	if (read_id_idx < 0) {
		throw InvalidInputException("ena_upload_reads: input relation must include a 'read_id' column");
	}
	if (sequence1_idx < 0) {
		throw InvalidInputException("ena_upload_reads: input relation must include a 'sequence1' column");
	}
	if (qual1_idx < 0) {
		throw InvalidInputException("ena_upload_reads: input relation must include a 'qual1' column");
	}
	const bool has_r2_columns = sequence2_idx >= 0 && qual2_idx >= 0;

	RequireType("sample_ref", types[sample_ref_idx], LogicalTypeId::VARCHAR);
	RequireType("read_id", types[read_id_idx], LogicalTypeId::VARCHAR);
	RequireType("sequence1", types[sequence1_idx], LogicalTypeId::VARCHAR);
	RequireListUtinyint("qual1", types[qual1_idx]);
	if (has_r2_columns) {
		RequireType("sequence2", types[sequence2_idx], LogicalTypeId::VARCHAR);
		RequireListUtinyint("qual2", types[qual2_idx]);
	}

	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}
		const idx_t n = chunk->size();

		UnifiedVectorFormat sample_ref_data, read_id_data, sequence1_data, qual1_data;
		chunk->data[sample_ref_idx].ToUnifiedFormat(n, sample_ref_data);
		chunk->data[read_id_idx].ToUnifiedFormat(n, read_id_data);
		chunk->data[sequence1_idx].ToUnifiedFormat(n, sequence1_data);
		chunk->data[qual1_idx].ToUnifiedFormat(n, qual1_data);

		auto sample_ref_strs = UnifiedVectorFormat::GetData<string_t>(sample_ref_data);
		auto read_id_strs = UnifiedVectorFormat::GetData<string_t>(read_id_data);
		auto sequence1_strs = UnifiedVectorFormat::GetData<string_t>(sequence1_data);

		UnifiedVectorFormat sequence2_data, qual2_data;
		const string_t *sequence2_strs = nullptr;
		if (has_r2_columns) {
			chunk->data[sequence2_idx].ToUnifiedFormat(n, sequence2_data);
			chunk->data[qual2_idx].ToUnifiedFormat(n, qual2_data);
			sequence2_strs = UnifiedVectorFormat::GetData<string_t>(sequence2_data);
		}

		UnifiedVectorFormat sequence_index_data;
		if (sequence_index_idx >= 0) {
			chunk->data[sequence_index_idx].ToUnifiedFormat(n, sequence_index_data);
		}

		for (idx_t row = 0; row < n; row++) {
			const idx_t global_row = out.sample_ref.size();
			if (global_row >= MAX_INPUT_ROWS) {
				throw InvalidInputException(
				    "ena_upload_reads: input relation exceeds %llu rows; split by sample_ref and call per-batch",
				    MAX_INPUT_ROWS);
			}
			const auto sr_idx = sample_ref_data.sel->get_index(row);
			if (!sample_ref_data.validity.RowIsValid(sr_idx)) {
				throw InvalidInputException("ena_upload_reads: NULL sample_ref at row %llu", global_row);
			}
			const auto sample_ref = sample_ref_strs[sr_idx].GetString();
			if (sample_ref.empty()) {
				throw InvalidInputException("ena_upload_reads: empty sample_ref at row %llu", global_row);
			}
			if (sample_ref.find('/') != string::npos) {
				throw InvalidInputException(
				    "ena_upload_reads: sample_ref '%s' contains '/' which is not allowed in a filename", sample_ref);
			}
			// Reject path-traversal components — the sample_ref is concatenated
			// directly into a destination path for `file://` writes, so
			// `..` could escape the upload directory.
			if (sample_ref == ".." || sample_ref == "." || sample_ref.find("..") != string::npos) {
				throw InvalidInputException(
				    "ena_upload_reads: sample_ref '%s' contains '..' (path-traversal not allowed)", sample_ref);
			}

			const auto rid_idx = read_id_data.sel->get_index(row);
			if (!read_id_data.validity.RowIsValid(rid_idx)) {
				throw InvalidInputException("ena_upload_reads: NULL read_id at row %llu", global_row);
			}
			const auto seq1_idx = sequence1_data.sel->get_index(row);
			if (!sequence1_data.validity.RowIsValid(seq1_idx)) {
				throw InvalidInputException("ena_upload_reads: NULL sequence1 at row %llu", global_row);
			}
			const auto q1_idx = qual1_data.sel->get_index(row);
			if (!qual1_data.validity.RowIsValid(q1_idx)) {
				throw InvalidInputException("ena_upload_reads: NULL qual1 at row %llu", global_row);
			}

			out.sample_ref.push_back(sample_ref);
			out.read_id.push_back(read_id_strs[rid_idx].GetString());
			out.sequence1.push_back(sequence1_strs[seq1_idx].GetString());
			auto qual1 = ExtractQualList(chunk->data[qual1_idx], qual1_data, row);
			if (qual1.size() != out.sequence1.back().size()) {
				throw InvalidInputException(
				    "ena_upload_reads: qual1 length (%llu) does not match sequence1 length (%llu) at row %llu",
				    qual1.size(), out.sequence1.back().size(), global_row);
			}
			out.qual1.push_back(std::move(qual1));
			out.total_seq_bytes += out.sequence1.back().size();

			if (has_r2_columns) {
				const auto s2 = sequence2_data.sel->get_index(row);
				const auto q2 = qual2_data.sel->get_index(row);
				const bool s2_valid = sequence2_data.validity.RowIsValid(s2);
				const bool q2_valid = qual2_data.validity.RowIsValid(q2);
				if (s2_valid != q2_valid) {
					throw InvalidInputException(
					    "ena_upload_reads: sequence2 and qual2 must both be NULL or both non-NULL at row %llu",
					    global_row);
				}
				out.has_seq2.push_back(s2_valid);
				if (s2_valid) {
					out.sequence2.push_back(sequence2_strs[s2].GetString());
					auto q2_vec = ExtractQualList(chunk->data[qual2_idx], qual2_data, row);
					if (q2_vec.size() != out.sequence2.back().size()) {
						throw InvalidInputException(
						    "ena_upload_reads: qual2 length (%llu) does not match sequence2 length (%llu) at row %llu",
						    q2_vec.size(), out.sequence2.back().size(), global_row);
					}
					out.qual2.push_back(std::move(q2_vec));
					out.total_seq_bytes += out.sequence2.back().size();
				} else {
					out.sequence2.emplace_back();
					out.qual2.emplace_back();
				}
			} else {
				out.has_seq2.push_back(false);
				out.sequence2.emplace_back();
				out.qual2.emplace_back();
			}

			if (out.total_seq_bytes > MAX_INPUT_SEQUENCE_BYTES) {
				throw InvalidInputException(
				    "ena_upload_reads: input sequences exceed %llu bytes (cap is to keep the in-memory pipeline "
				    "bounded; split by sample_ref and call per-batch)",
				    MAX_INPUT_SEQUENCE_BYTES);
			}

			if (sequence_index_idx >= 0) {
				const auto si = sequence_index_data.sel->get_index(row);
				if (sequence_index_data.validity.RowIsValid(si)) {
					out.sequence_index.push_back(UnifiedVectorFormat::GetData<int64_t>(sequence_index_data)[si]);
					out.sequence_index_valid.push_back(true);
				} else {
					out.sequence_index.push_back(0);
					out.sequence_index_valid.push_back(false);
				}
			} else {
				out.sequence_index.push_back(0);
				out.sequence_index_valid.push_back(false);
			}
		}
	}
}

// =====================================================================
// Grouping + layout resolution
// =====================================================================

void GroupAndResolveLayout(const UploadInputRows &input, FastqLayoutMode requested, vector<SampleGroup> &out) {
	std::map<string, vector<size_t>> by_sample;
	for (size_t i = 0; i < input.sample_ref.size(); i++) {
		by_sample[input.sample_ref[i]].push_back(i);
	}
	for (auto &kv : by_sample) {
		SampleGroup g;
		g.sample_ref = kv.first;
		g.row_indices = std::move(kv.second);
		// Sort by sequence_index, NULLS LAST. Stable so equal/null rows keep
		// original input order, which matches reader expectations.
		std::stable_sort(g.row_indices.begin(), g.row_indices.end(), [&input](size_t a, size_t b) {
			const bool av = input.sequence_index_valid[a];
			const bool bv = input.sequence_index_valid[b];
			if (av && bv) {
				return input.sequence_index[a] < input.sequence_index[b];
			}
			if (av) {
				return true; // valid before null
			}
			return false;
		});
		vector<bool> has_r2;
		has_r2.reserve(g.row_indices.size());
		for (auto idx : g.row_indices) {
			has_r2.push_back(input.has_seq2[idx]);
		}
		g.resolved_layout = ResolveLayout(g.sample_ref, requested, has_r2);
		out.push_back(std::move(g));
	}
}

// =====================================================================
// Encoder → gzip → MD5 stream
// =====================================================================
//
// Shared core for the file-sink (push, eager file write) and the libcurl
// streaming producer (pull, in-memory buffering). Both run identical zlib +
// MD5 plumbing — only the destination of the gzipped output differs. The
// destination is supplied via a `ChunkSink` callback so consumers can route
// bytes wherever they need.
//
// MD5 is computed over the bytes handed to the sink, which is exactly what
// the consumer ultimately writes / uploads — the digest matches.

class GzipMd5Stream {
public:
	using ChunkSink = std::function<void(const uint8_t *data, std::size_t size)>;

	explicit GzipMd5Stream(ChunkSink sink_fn) : sink(std::move(sink_fn)) {
		std::memset(&zs, 0, sizeof(zs));
		// 16 + MAX_WBITS selects gzip framing (vs raw deflate or zlib wrapper).
		if (deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 16 + MAX_WBITS, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
			throw IOException("ena_upload_reads: deflateInit2 failed");
		}
		zs_initialized = true;
		work_buf.resize(64 * 1024);
	}

	~GzipMd5Stream() {
		if (zs_initialized) {
			deflateEnd(&zs);
		}
	}

	GzipMd5Stream(const GzipMd5Stream &) = delete;
	GzipMd5Stream &operator=(const GzipMd5Stream &) = delete;

	// Feed `size` bytes of plain input. Gzipped output is forwarded to `sink`
	// in 64 KB chunks as deflate produces it.
	void Write(const uint8_t *data, std::size_t size) {
		if (finished) {
			throw IOException("ena_upload_reads: GzipMd5Stream::Write after Finish");
		}
		// zlib's avail_in is `uInt` (32-bit). Loop in case a caller passes a
		// larger buffer than fits in one zlib call.
		const uint8_t *p = data;
		std::size_t remaining = size;
		while (remaining > 0) {
			const std::size_t step = std::min<std::size_t>(remaining, std::numeric_limits<uInt>::max());
			zs.next_in = const_cast<unsigned char *>(p);
			zs.avail_in = static_cast<uInt>(step);
			while (zs.avail_in > 0) {
				DrainOnce(Z_NO_FLUSH);
			}
			p += step;
			remaining -= step;
		}
	}

	// Flush remaining gzip state with Z_FINISH until deflate returns
	// Z_STREAM_END. Returns (md5_hex, bytes_emitted_to_sink). Idempotent —
	// calling twice throws.
	std::pair<string, uint64_t> Finish() {
		if (finished) {
			throw IOException("ena_upload_reads: GzipMd5Stream::Finish called twice");
		}
		finished = true;
		zs.next_in = nullptr;
		zs.avail_in = 0;
		while (true) {
			zs.next_out = work_buf.data();
			zs.avail_out = static_cast<uInt>(work_buf.size());
			int rc = deflate(&zs, Z_FINISH);
			if (rc < 0 && rc != Z_BUF_ERROR) {
				throw IOException("ena_upload_reads: deflate(Z_FINISH) failed: %s", zs.msg ? zs.msg : "<no msg>");
			}
			const std::size_t produced = work_buf.size() - zs.avail_out;
			if (produced > 0) {
				EmitChunk(produced);
			}
			if (rc == Z_STREAM_END) {
				break;
			}
			if (produced == 0) {
				// deflate(Z_FINISH) must either emit output or return
				// Z_STREAM_END given a non-empty output buffer; otherwise
				// we'd spin forever.
				throw IOException("ena_upload_reads: deflate(Z_FINISH) made no progress");
			}
		}
		return {md5_ctx.FinishHex(), bytes_emitted};
	}

	bool Finished() const {
		return finished;
	}

private:
	void DrainOnce(int flush_mode) {
		zs.next_out = work_buf.data();
		zs.avail_out = static_cast<uInt>(work_buf.size());
		int rc = deflate(&zs, flush_mode);
		if (rc < 0) {
			throw IOException("ena_upload_reads: deflate failed: %s", zs.msg ? zs.msg : "<no msg>");
		}
		const std::size_t produced = work_buf.size() - zs.avail_out;
		if (produced > 0) {
			EmitChunk(produced);
		}
	}

	void EmitChunk(std::size_t n) {
		md5_ctx.Add(work_buf.data(), n);
		sink(work_buf.data(), n);
		bytes_emitted += n;
	}

	ChunkSink sink;
	z_stream zs {};
	bool zs_initialized = false;
	bool finished = false;
	vector<uint8_t> work_buf;
	uint64_t bytes_emitted = 0;
	MD5Context md5_ctx;
};

// File-backed wrapper. Sink callback writes each gzipped chunk straight to
// the open FileHandle — no extra buffering between deflate and disk.
class GzipMd5FileSink {
public:
	GzipMd5FileSink(FileSystem &fs, const string &path) : fs(fs), path(path) {
		file = fs.OpenFile(path, FileFlags::FILE_FLAGS_WRITE | FileFlags::FILE_FLAGS_FILE_CREATE);
		// Lambda captures `this`; `stream` is destroyed before `file` (declared
		// after it), so the lambda's `*this->file` access is always valid.
		stream = make_uniq<GzipMd5Stream>([this](const uint8_t *p, std::size_t n) {
			this->fs.Write(*this->file, const_cast<uint8_t *>(p), static_cast<int64_t>(n));
		});
	}

	GzipMd5FileSink(const GzipMd5FileSink &) = delete;
	GzipMd5FileSink &operator=(const GzipMd5FileSink &) = delete;

	void Write(const char *data, std::size_t size) {
		stream->Write(reinterpret_cast<const uint8_t *>(data), size);
	}

	// Returns (md5_hex, bytes_written). Closes the file.
	std::pair<string, uint64_t> Finish() {
		auto result = stream->Finish();
		file->Close();
		return result;
	}

private:
	FileSystem &fs;
	string path;
	// Order matters: `file` is destroyed AFTER `stream`, so the sink lambda
	// stored inside `stream` can safely reference `*file` for the lifetime of
	// the stream.
	unique_ptr<FileHandle> file;
	unique_ptr<GzipMd5Stream> stream;
};

// =====================================================================
// Streaming encode → gzip → MD5 producer for libcurl
// =====================================================================
//
// libcurl pulls bytes via CURLOPT_READFUNCTION; this producer satisfies that
// contract by encoding records on demand. Wraps a GzipMd5Stream whose sink
// appends gzipped chunks into `out_buf`; `Read` drains `out_buf` to libcurl
// and produces more on demand.
//
// State machine:
//   1) Encoding: fetch next row, encode FASTQ via the encoder, feed the
//      record bytes to `stream.Write` — gzip output lands in `out_buf`.
//   2) Finishing: when no more rows, call `stream.Finish()` (Z_FINISH)
//      which drains the trailer into `out_buf` and yields the final
//      md5 + bytes_emitted, which we cache for `Finish()`.
//   3) EOF: Read returns 0 once `out_buf` is drained AND `stream.Finished()`.
//
// `Finish()` is valid only after Read has returned 0.

class StreamingGzipMd5Producer {
public:
	StreamingGzipMd5Producer(const UploadInputRows &input, const SampleGroup &group, FastqLayoutMode layout,
	                         int which_file, uint8_t qual_offset)
	    : input(input), group(group), layout(layout), which_file(which_file), encoder(qual_offset),
	      stream([this](const uint8_t *p, std::size_t n) { out_buf.insert(out_buf.end(), p, p + n); }) {
	}

	StreamingGzipMd5Producer(const StreamingGzipMd5Producer &) = delete;
	StreamingGzipMd5Producer &operator=(const StreamingGzipMd5Producer &) = delete;

	std::size_t Read(char *buf, std::size_t max_bytes) {
		while (true) {
			if (out_buf_pos < out_buf.size()) {
				const std::size_t avail = out_buf.size() - out_buf_pos;
				const std::size_t to_copy = std::min(avail, max_bytes);
				std::memcpy(buf, out_buf.data() + out_buf_pos, to_copy);
				out_buf_pos += to_copy;
				return to_copy;
			}
			out_buf.clear();
			out_buf_pos = 0;
			if (stream.Finished()) {
				return 0;
			}
			ProduceMore();
		}
	}

	std::pair<string, uint64_t> Finish() {
		if (!stream.Finished()) {
			throw IOException("ena_upload_reads: StreamingGzipMd5Producer::Finish before EOF");
		}
		return cached_result;
	}

private:
	void ProduceMore() {
		std::vector<uint8_t> record_buf;
		auto record_sink = [&record_buf](const char *data, std::size_t size) {
			record_buf.insert(record_buf.end(), reinterpret_cast<const uint8_t *>(data),
			                  reinterpret_cast<const uint8_t *>(data) + size);
		};
		if (!EncodeOneRecord(record_sink)) {
			cached_result = stream.Finish();
			return;
		}
		stream.Write(record_buf.data(), record_buf.size());
	}

	bool EncodeOneRecord(const FastqEncoder::Sink &sink) {
		if (row_cursor >= group.row_indices.size()) {
			return false;
		}
		const auto rid = group.row_indices[row_cursor++];
		const char *id = input.read_id[rid].data();
		const std::size_t id_len = input.read_id[rid].size();

		if (layout == FastqLayoutMode::SINGLE || (layout == FastqLayoutMode::PAIRED && which_file == 0)) {
			const auto &seq = input.sequence1[rid];
			const auto &q = input.qual1[rid];
			encoder.Encode(sink, id, id_len, nullptr, 0, seq.data(), seq.size(), q.data(), q.size());
		} else if (layout == FastqLayoutMode::PAIRED && which_file == 1) {
			const auto &seq = input.sequence2[rid];
			const auto &q = input.qual2[rid];
			encoder.Encode(sink, id, id_len, nullptr, 0, seq.data(), seq.size(), q.data(), q.size());
		} else if (layout == FastqLayoutMode::PAIRED_INTERLEAVED) {
			const auto &s1 = input.sequence1[rid];
			const auto &q1 = input.qual1[rid];
			encoder.Encode(sink, id, id_len, nullptr, 0, s1.data(), s1.size(), q1.data(), q1.size());
			const auto &s2 = input.sequence2[rid];
			const auto &q2 = input.qual2[rid];
			encoder.Encode(sink, id, id_len, nullptr, 0, s2.data(), s2.size(), q2.data(), q2.size());
		} else {
			// AUTO must be resolved by ResolveLayout before this class is
			// instantiated. A future enum addition would land here too.
			throw std::logic_error("StreamingGzipMd5Producer::EncodeOneRecord: unhandled layout");
		}
		return true;
	}

	const UploadInputRows &input;
	const SampleGroup &group;
	const FastqLayoutMode layout;
	const int which_file;
	FastqEncoder encoder;

	// `out_buf` MUST be declared before `stream` so that the sink lambda
	// stored inside `stream` is destroyed (with `stream`) before `out_buf` —
	// otherwise the lambda's last access could touch a destroyed vector.
	std::vector<std::uint8_t> out_buf;
	std::size_t out_buf_pos = 0;
	std::size_t row_cursor = 0;
	std::pair<string, uint64_t> cached_result;
	GzipMd5Stream stream;
};

// =====================================================================
// Sample → file write
// =====================================================================

// Encode the rows of one sample into the supplied sink. `which_file` controls
// which side of the layout we emit (0 = R1 / single / interleaved; 1 = R2
// only, used by the split-paired layout's second pass).
//
// `qual_offset` is wired through the encoder constructor by the caller; this
// helper only forwards reads, so the active offset is implicit in `encoder`.
void EncodeSampleRows(const UploadInputRows &in, const SampleGroup &group, FastqLayoutMode layout, int which_file,
                      FastqEncoder &encoder, GzipMd5FileSink &sink) {
	auto sink_fn = [&sink](const char *data, std::size_t size) {
		sink.Write(data, size);
	};

	for (auto rid : group.row_indices) {
		const char *id = in.read_id[rid].data();
		const std::size_t id_len = in.read_id[rid].size();

		if (layout == FastqLayoutMode::SINGLE || (layout == FastqLayoutMode::PAIRED && which_file == 0)) {
			const auto &seq = in.sequence1[rid];
			const auto &q = in.qual1[rid];
			encoder.Encode(sink_fn, id, id_len, nullptr, 0, seq.data(), seq.size(), q.data(), q.size());
		} else if (layout == FastqLayoutMode::PAIRED && which_file == 1) {
			const auto &seq = in.sequence2[rid];
			const auto &q = in.qual2[rid];
			encoder.Encode(sink_fn, id, id_len, nullptr, 0, seq.data(), seq.size(), q.data(), q.size());
		} else if (layout == FastqLayoutMode::PAIRED_INTERLEAVED) {
			const auto &s1 = in.sequence1[rid];
			const auto &q1 = in.qual1[rid];
			encoder.Encode(sink_fn, id, id_len, nullptr, 0, s1.data(), s1.size(), q1.data(), q1.size());
			const auto &s2 = in.sequence2[rid];
			const auto &q2 = in.qual2[rid];
			encoder.Encode(sink_fn, id, id_len, nullptr, 0, s2.data(), s2.size(), q2.data(), q2.size());
		} else {
			// AUTO must be resolved before this point; see ResolveLayout.
			throw std::logic_error("EncodeSampleRows: unhandled FastqLayoutMode");
		}
	}
}

void RunUpload(ClientContext &context, const ENAUploadReadsBindData &bind, ENAUploadReadsGlobalState &gs) {
	auto &fs = FileSystem::GetFileSystem(context);

	if (gs.target.transport == UploadTransport::LOCAL_FILE) {
		fs.CreateDirectoriesRecursive(gs.target.remote_dir);
	}

	for (auto &group : gs.samples) {
		auto filenames = OutputFilenames(group.sample_ref, group.resolved_layout);
		for (size_t f = 0; f < filenames.size(); f++) {
			const auto &fname = filenames[f];
			const int which_file = (group.resolved_layout == FastqLayoutMode::PAIRED && f == 1) ? 1 : 0;

			string md5_hex;
			uint64_t bytes_written = 0;

			if (gs.target.transport == UploadTransport::CURL) {
#ifdef MIINT_HAS_CURL
				// Streaming path — no temp file. libcurl pulls bytes from
				// the producer until EOF; the producer encodes records on
				// demand and computes MD5 over the gzipped output as it
				// flows through.
				StreamingGzipMd5Producer producer(gs.input, group, group.resolved_layout, which_file, bind.qual_offset);
				miint::CurlUploadOptions opts;
				opts.url = gs.target.url_for_curl + fname;
				opts.user = gs.user;
				opts.password = gs.password;
				opts.create_dirs = true;
				// `MIINT_CURL_VERBOSE=1` enables CURLOPT_VERBOSE wire trace on
				// stderr — invaluable for diagnosing FTPS/HTTPS upload hangs
				// against the live ENA endpoint.
				if (const char *v = std::getenv("MIINT_CURL_VERBOSE"); v && std::string(v) == "1") {
					opts.verbose = true;
				}
				auto curl_producer = [&producer](char *buf, std::size_t max_bytes) -> std::size_t {
					return producer.Read(buf, max_bytes);
				};
				auto result = miint::RunCurlUpload(opts, curl_producer);
				if (!result.error_message.empty()) {
					throw IOException("ena_upload_reads: %s upload failed for sample '%s' file '%s': %s",
					                  gs.target.scheme, group.sample_ref, fname, result.error_message);
				}
				auto finished = producer.Finish();
				md5_hex = std::move(finished.first);
				bytes_written = finished.second;
#else
				throw IOException("ena_upload_reads: %s:// transport requires libcurl, which is disabled in this build "
				                  "(use aspera:// or file:// instead, or rebuild with -DMIINT_ENABLE_CURL=ON)",
				                  gs.target.scheme);
#endif
			} else {
				// Encode → gzip → MD5 → file. For Aspera we then ascp the
				// file. For LOCAL_FILE the destination IS the local path.
				string write_path;
				string temp_dir;
				if (gs.target.transport == UploadTransport::LOCAL_FILE) {
					write_path = gs.target.remote_dir + fname;
				} else {
					// Aspera: write to a private temp directory so the
					// basename matches the desired remote name. mkdtemp
					// gives us a unique path; we unlink + rmdir after
					// ascp finishes.
					string tmpl = miint::GetTempDir() + "/ena-upload-XXXXXX";
					vector<char> buf(tmpl.begin(), tmpl.end());
					buf.push_back('\0');
					if (mkdtemp(buf.data()) == nullptr) {
						throw IOException("ena_upload_reads: mkdtemp failed: %s", std::strerror(errno));
					}
					temp_dir.assign(buf.data());
					write_path = temp_dir + "/" + fname;
				}

				// Cleanup helper — runs whether the encode/gzip/transport
				// block succeeds or throws. For LOCAL_FILE we leave the
				// destination in place regardless (matches COPY behaviour).
				auto cleanup_temp = [&]() noexcept {
					if (!temp_dir.empty()) {
						(void)unlink(write_path.c_str());
						(void)rmdir(temp_dir.c_str());
					}
				};

				try {
					FastqEncoder encoder(bind.qual_offset);
					GzipMd5FileSink sink(fs, write_path);
					EncodeSampleRows(gs.input, group, group.resolved_layout, which_file, encoder, sink);
					auto finished = sink.Finish();
					md5_hex = std::move(finished.first);
					bytes_written = finished.second;

#if MIINT_ASPERA_SUPPORTED && defined(MIINT_STATIC_BUILD)
					if (gs.target.transport == UploadTransport::ASPERA) {
						AsperaSendOptions opts;
						opts.ascp_path = gs.aspera_ascp_path;
						opts.key_path = gs.aspera_key_path;
						opts.user = gs.user;
						opts.host = gs.target.host;
						opts.port = 33001;
						opts.local_path = write_path;
						opts.remote_dir = gs.target.remote_dir;
						opts.max_rate = gs.aspera_max_rate;
						auto argv = BuildAscpSendArgv(opts);

						auto result = miint::RunAsperaSend(argv, gs.password);
						if (result.exit_code != 0) {
							throw IOException("ena_upload_reads: ascp failed (exit %d) for sample '%s' file '%s': %s",
							                  result.exit_code, group.sample_ref, fname, result.stderr_output);
						}
					}
#else
					if (gs.target.transport == UploadTransport::ASPERA) {
						throw IOException("ena_upload_reads: Aspera transport is not supported on this build/platform");
					}
#endif
				} catch (...) {
					cleanup_temp();
					throw;
				}
				cleanup_temp();
			}

			EmittedFile e;
			e.sample_ref = group.sample_ref;
			e.filename = fname;
			e.filetype = "fastq";
			e.md5_hex = std::move(md5_hex);
			e.bytes_written = bytes_written;
			e.layout_name = FastqLayoutModeName(group.resolved_layout);
			gs.emitted.push_back(std::move(e));
		}
	}
}

// =====================================================================
// Secret resolution (authenticated transports — Aspera, libcurl)
// =====================================================================

void ResolveUploadCredentials(ClientContext &context, const string &transport_label, const string &secret_name,
                              string &user_out, string &password_out) {
	if (secret_name.empty()) {
		throw BinderException("ena_upload_reads: %s transport requires a SECRET — pass `secret := 'name'`",
		                      transport_label);
	}
	auto &mgr = SecretManager::Get(context);
	auto txn = CatalogTransaction::GetSystemCatalogTransaction(context);
	auto entry = mgr.GetSecretByName(txn, secret_name);
	if (!entry) {
		throw BinderException("ena_upload_reads: secret '%s' not found", secret_name);
	}
	auto kv = dynamic_cast<const KeyValueSecret *>(entry->secret.get());
	if (!kv) {
		throw InvalidInputException("ena_upload_reads: secret '%s' is not a KeyValueSecret", secret_name);
	}
	auto user_val = kv->TryGetValue("user");
	auto password_val = kv->TryGetValue("password");
	if (user_val.IsNull() || password_val.IsNull()) {
		throw BinderException("ena_upload_reads: secret '%s' is missing user or password", secret_name);
	}
	user_out = user_val.ToString();
	password_out = password_val.ToString();
}

// =====================================================================
// Bind / InitGlobal / Execute
// =====================================================================

unique_ptr<FunctionData> Bind(ClientContext &, TableFunctionBindInput &input, vector<LogicalType> &return_types,
                              vector<string> &names) {
	auto bind = make_uniq<ENAUploadReadsBindData>();

	auto get_string_param = [&](const char *key, string &out, bool required) {
		auto it = input.named_parameters.find(key);
		if (it == input.named_parameters.end() || it->second.IsNull()) {
			if (required) {
				throw BinderException("ena_upload_reads: required named parameter '%s' is missing", key);
			}
			return;
		}
		out = it->second.ToString();
	};

	get_string_param("relation", bind->relation_name, /*required=*/true);
	get_string_param("secret", bind->secret_name, /*required=*/false);
	string target_url = "aspera://webin2.ebi.ac.uk/";
	get_string_param("target_url", target_url, /*required=*/false);
	bind->target_url = target_url;

	auto qo_it = input.named_parameters.find("qual_offset");
	if (qo_it != input.named_parameters.end() && !qo_it->second.IsNull()) {
		const auto v = qo_it->second.GetValue<int64_t>();
		if (v != 33 && v != 64) {
			throw BinderException("ena_upload_reads: qual_offset must be 33 or 64 (got %lld)", v);
		}
		bind->qual_offset = static_cast<uint8_t>(v);
	}

	string layout_name = "auto";
	get_string_param("layout", layout_name, /*required=*/false);
	try {
		bind->layout_mode = ParseFastqLayoutMode(layout_name);
	} catch (const std::exception &e) {
		throw BinderException("ena_upload_reads: %s", e.what());
	}

	auto rate_it = input.named_parameters.find("aspera_rate_limit_mbps");
	if (rate_it != input.named_parameters.end() && !rate_it->second.IsNull()) {
		const auto v = rate_it->second.GetValue<int64_t>();
		if (v <= 0) {
			throw BinderException("ena_upload_reads: aspera_rate_limit_mbps must be a positive integer");
		}
		bind->aspera_max_rate = std::to_string(v) + "m";
	}

	names = {"sample_ref", "filename", "filetype", "md5", "bytes_written", "layout"};
	return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
	                LogicalType::VARCHAR, LogicalType::UBIGINT, LogicalType::VARCHAR};
	return std::move(bind);
}

unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind = input.bind_data->Cast<ENAUploadReadsBindData>();
	auto gs = make_uniq<ENAUploadReadsGlobalState>();

	UploadTargetURL parsed;
	try {
		parsed = ParseUploadTargetURL(bind.target_url);
	} catch (const std::exception &e) {
		throw BinderException("ena_upload_reads: invalid target_url '%s': %s", bind.target_url, e.what());
	}
	gs->target = parsed;
	gs->aspera_max_rate = bind.aspera_max_rate;

	if (gs->target.transport == UploadTransport::ASPERA) {
		ResolveUploadCredentials(context, "Aspera", bind.secret_name, gs->user, gs->password);
		auto &db = DatabaseInstance::GetDatabase(context);
		gs->aspera_ascp_path = miint::AsperaUtils::FindAscp();
		if (gs->aspera_ascp_path.empty()) {
			throw IOException("ena_upload_reads: 'ascp' not found on PATH; install Aspera CLI or use file:// "
			                  "transport for local testing");
		}
		gs->aspera_key_path = miint::AsperaUtils::ResolveKey(db, gs->aspera_ascp_path, /*required=*/true);
	} else if (gs->target.transport == UploadTransport::CURL) {
#ifdef MIINT_HAS_CURL
		// Plain user/password auth via Basic (HTTP) or USERPWD (FTP). Most
		// hosts will require this; we always pass them through if a secret
		// was supplied. A future "anonymous FTP" use case could relax this
		// to make `secret` optional for ftp:// targets.
		ResolveUploadCredentials(context, "libcurl", bind.secret_name, gs->user, gs->password);
#else
		throw BinderException("ena_upload_reads: %s:// transport requires libcurl, which is disabled in this build "
		                      "(use aspera:// or file:// instead, or rebuild with -DMIINT_ENABLE_CURL=ON)",
		                      gs->target.scheme);
#endif
	}

	MaterialiseInput(context, bind.relation_name, gs->input);
	GroupAndResolveLayout(gs->input, bind.layout_mode, gs->samples);

	// Run the entire upload eagerly during init. The execution side just
	// drains the per-file rows from `emitted`. This keeps the error model
	// simple (any failure aborts the statement) and matches the typical
	// per-batch upload size — usually under a few GB of FASTQ.
	RunUpload(context, bind, *gs);
	return std::move(gs);
}

void Execute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &gs = data_p.global_state->Cast<ENAUploadReadsGlobalState>();
	const idx_t remaining = gs.emitted.size() - gs.emit_cursor;
	if (remaining == 0) {
		output.SetCardinality(0);
		return;
	}
	const idx_t to_emit = std::min<idx_t>(remaining, STANDARD_VECTOR_SIZE);

	auto sample_ref = FlatVector::GetData<string_t>(output.data[0]);
	auto filename = FlatVector::GetData<string_t>(output.data[1]);
	auto filetype = FlatVector::GetData<string_t>(output.data[2]);
	auto md5_v = FlatVector::GetData<string_t>(output.data[3]);
	auto bytes_v = FlatVector::GetData<uint64_t>(output.data[4]);
	auto layout_v = FlatVector::GetData<string_t>(output.data[5]);

	for (idx_t i = 0; i < to_emit; i++) {
		const auto &row = gs.emitted[gs.emit_cursor + i];
		sample_ref[i] = StringVector::AddString(output.data[0], row.sample_ref);
		filename[i] = StringVector::AddString(output.data[1], row.filename);
		filetype[i] = StringVector::AddString(output.data[2], row.filetype);
		md5_v[i] = StringVector::AddString(output.data[3], row.md5_hex);
		bytes_v[i] = row.bytes_written;
		layout_v[i] = StringVector::AddString(output.data[5], row.layout_name);
	}
	gs.emit_cursor += to_emit;
	output.SetCardinality(to_emit);
}

} // namespace

TableFunction ENAUploadReadsTableFunction::GetFunction() {
	TableFunction tf("ena_upload_reads", {}, Execute, Bind, InitGlobal);
	tf.named_parameters["relation"] = LogicalType::VARCHAR;
	tf.named_parameters["secret"] = LogicalType::VARCHAR;
	tf.named_parameters["target_url"] = LogicalType::VARCHAR;
	tf.named_parameters["qual_offset"] = LogicalType::BIGINT;
	tf.named_parameters["layout"] = LogicalType::VARCHAR;
	tf.named_parameters["aspera_rate_limit_mbps"] = LogicalType::BIGINT;
	return tf;
}

void ENAUploadReadsTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
