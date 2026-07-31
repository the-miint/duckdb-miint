// SPDX-License-Identifier: MIT
//
// `ena_upload_reads` table function — see ena_upload_reads.hpp.

#include "ena_upload_reads.hpp"

#include "aspera_send.hpp"
#include "aspera_utils.hpp"
#ifdef MIINT_HAS_CURL
#include "catalog_utils.hpp"
#include "curl_send.hpp"
#endif
#include "ena_upload_helpers.hpp"
#include "fastq_encoder.hpp"
#include "remote_file_helper.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/random_engine.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/crypto/md5.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/main/secret/secret.hpp"
#include "duckdb/main/secret/secret_manager.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <unistd.h> // unlink, rmdir for temp-staging cleanup (available on MinGW)
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
// Per-sample upload plan
// =====================================================================
// One per distinct sample_ref, produced by the aggregate pre-pass. Carries only
// the metadata needed to drive that sample's streaming data pass — no row
// payload is held, which is what keeps peak memory bounded regardless of how
// large the dataset (or any single sample) is.
struct SamplePlan {
	string sample_ref;
	FastqLayoutMode resolved_layout;
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
	vector<SamplePlan> samples;
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
// Schema validation, sample planning, and row extraction (read-time)
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

// Validate the input relation's required columns + types against a zero-row
// probe, so we never materialise data just to check the schema. Returns whether
// the optional R2 columns (sequence2 + qual2) are present.
bool ValidateSchemaDetectR2(Connection &conn, const string &relation_name) {
	const string query = "SELECT * FROM " + KeywordHelper::WriteOptionallyQuoted(relation_name) + " LIMIT 0";
	auto result = conn.Query(query);
	if (result->HasError()) {
		throw InvalidInputException("ena_upload_reads: failed to read relation '%s': %s", relation_name,
		                            result->GetError());
	}
	auto &names = result->names;
	auto &types = result->types;

	const int sample_ref_idx = FindColumn(names, "sample_ref");
	const int read_id_idx = FindColumn(names, "read_id");
	const int sequence1_idx = FindColumn(names, "sequence1");
	const int qual1_idx = FindColumn(names, "qual1");
	const int sequence2_idx = FindColumn(names, "sequence2");
	const int qual2_idx = FindColumn(names, "qual2");

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
	RequireType("sample_ref", types[sample_ref_idx], LogicalTypeId::VARCHAR);
	RequireType("read_id", types[read_id_idx], LogicalTypeId::VARCHAR);
	RequireType("sequence1", types[sequence1_idx], LogicalTypeId::VARCHAR);
	RequireListUtinyint("qual1", types[qual1_idx]);

	// R2 is opt-in but all-or-nothing: a relation with only one of the two
	// columns is a malformed schema, not a single-end relation. Treating it as
	// single-end would silently drop the present mate, so reject it loudly.
	const bool has_sequence2 = sequence2_idx >= 0;
	const bool has_qual2 = qual2_idx >= 0;
	if (has_sequence2 != has_qual2) {
		throw InvalidInputException(
		    "ena_upload_reads: input relation has '%s' but not '%s' — paired input requires both (or neither)",
		    has_sequence2 ? "sequence2" : "qual2", has_sequence2 ? "qual2" : "sequence2");
	}
	const bool has_r2_columns = has_sequence2; // == has_qual2
	if (has_r2_columns) {
		RequireType("sequence2", types[sequence2_idx], LogicalTypeId::VARCHAR);
		RequireListUtinyint("qual2", types[qual2_idx]);
	}
	return has_r2_columns;
}

// A sample_ref becomes a filename component (and, for file://, a path segment),
// so reject anything that could break the name or escape the upload directory.
void ValidateSampleRef(const string &sample_ref) {
	if (sample_ref.empty()) {
		throw InvalidInputException("ena_upload_reads: empty sample_ref");
	}
	if (sample_ref.find('/') != string::npos) {
		throw InvalidInputException("ena_upload_reads: sample_ref '%s' contains '/' which is not allowed in a filename",
		                            sample_ref);
	}
	if (sample_ref == "." || sample_ref.find("..") != string::npos) {
		throw InvalidInputException("ena_upload_reads: sample_ref '%s' is '.' or contains '..' (path-traversal not "
		                            "allowed)",
		                            sample_ref);
	}
}

// Cheap aggregate pre-pass: a single projected GROUP BY scan that learns each
// sample's R2 pattern — `all_paired` (bool_and) and `any_paired` (bool_or) of
// "sequence2 IS NOT NULL" — without ever touching the sequence payload. Resolves
// every sample's layout up front, failing fast on a mixed sample or a
// layout/mode conflict BEFORE any upload begins (so an invalid sample can't
// leave earlier samples half-uploaded), and yields the list of samples to
// stream. O(#samples) memory.
void PlanSamples(Connection &conn, const string &relation_name, FastqLayoutMode requested, bool has_r2_columns,
                 vector<SamplePlan> &out) {
	const string quoted = KeywordHelper::WriteOptionallyQuoted(relation_name);
	const string r2_expr = has_r2_columns ? "sequence2 IS NOT NULL" : "false";
	const string query = "SELECT sample_ref, bool_and(" + r2_expr + ") AS all_paired, bool_or(" + r2_expr +
	                     ") AS any_paired FROM " + quoted + " GROUP BY sample_ref";
	auto result = conn.Query(query);
	if (result->HasError()) {
		throw InvalidInputException("ena_upload_reads: failed to scan relation '%s': %s", relation_name,
		                            result->GetError());
	}
	while (auto chunk = result->Fetch()) {
		const idx_t n = chunk->size();
		for (idx_t row = 0; row < n; row++) {
			auto sr_val = chunk->data[0].GetValue(row);
			if (sr_val.IsNull()) {
				throw InvalidInputException("ena_upload_reads: NULL sample_ref");
			}
			const string sample_ref = sr_val.ToString();
			ValidateSampleRef(sample_ref);
			// bool_and / bool_or over a non-empty group are never NULL.
			const bool all_paired = chunk->data[1].GetValue(row).GetValue<bool>();
			const bool any_paired = chunk->data[2].GetValue(row).GetValue<bool>();
			SamplePlan plan;
			plan.sample_ref = sample_ref;
			try {
				plan.resolved_layout = ResolveLayoutFromCounts(sample_ref, requested, all_paired, any_paired);
			} catch (const std::exception &e) {
				throw InvalidInputException("ena_upload_reads: %s", e.what());
			}
			out.push_back(std::move(plan));
		}
	}
}

// =====================================================================
// Encoder → gzip → MD5 stream
// =====================================================================
//
// Streaming zlib + MD5 plumbing. The gzipped output is handed to a `ChunkSink`
// callback (here always a file write — every transport stages to a file). MD5 is
// computed over the bytes handed to the sink, which is exactly what the consumer
// ultimately writes / uploads, so the digest matches the on-disk file.

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
// Per-chunk FASTQ encoding
// =====================================================================
//
// Encode the rows of one streamed DataChunk into the sample's gzip sink(s). The
// data query selects a fixed column order — 0=read_id, 1=sequence1, 2=qual1,
// (3=sequence2, 4=qual2) — so indices are positional. R1 always goes to sink0;
// for split PAIRED, R2 goes to sink1; for PAIRED_INTERLEAVED, R2 follows R1 into
// sink0 (one pass over the row → positional pairing by construction). The R2
// columns are guaranteed present whenever the resolved layout needs them
// (PlanSamples would have rejected the sample otherwise).
void EncodeChunk(DataChunk &chunk, FastqLayoutMode layout, FastqEncoder &encoder, GzipMd5FileSink &sink0,
                 GzipMd5FileSink *sink1, const string &sample_ref) {
	const idx_t n = chunk.size();

	UnifiedVectorFormat read_id_data, sequence1_data, qual1_data;
	chunk.data[0].ToUnifiedFormat(n, read_id_data);
	chunk.data[1].ToUnifiedFormat(n, sequence1_data);
	chunk.data[2].ToUnifiedFormat(n, qual1_data);
	auto read_id_strs = UnifiedVectorFormat::GetData<string_t>(read_id_data);
	auto sequence1_strs = UnifiedVectorFormat::GetData<string_t>(sequence1_data);

	const bool need_r2 = layout == FastqLayoutMode::PAIRED || layout == FastqLayoutMode::PAIRED_INTERLEAVED;
	UnifiedVectorFormat sequence2_data, qual2_data;
	const string_t *sequence2_strs = nullptr;
	if (need_r2) {
		chunk.data[3].ToUnifiedFormat(n, sequence2_data);
		chunk.data[4].ToUnifiedFormat(n, qual2_data);
		sequence2_strs = UnifiedVectorFormat::GetData<string_t>(sequence2_data);
	}

	auto write0 = [&sink0](const char *data, std::size_t size) {
		sink0.Write(data, size);
	};
	// R2 destination: split PAIRED → sink1; PAIRED_INTERLEAVED → sink0 (R2 right
	// after R1 in the same file). `r2_sink` is non-null exactly when `need_r2`, so
	// `write_r2` is only ever invoked through a valid pointer.
	GzipMd5FileSink *r2_sink = need_r2 ? (layout == FastqLayoutMode::PAIRED ? sink1 : &sink0) : nullptr;
	auto write_r2 = [r2_sink](const char *data, std::size_t size) {
		r2_sink->Write(data, size);
	};

	for (idx_t row = 0; row < n; row++) {
		const auto rid_i = read_id_data.sel->get_index(row);
		if (!read_id_data.validity.RowIsValid(rid_i)) {
			throw InvalidInputException("ena_upload_reads: NULL read_id in sample '%s'", sample_ref);
		}
		const auto s1_i = sequence1_data.sel->get_index(row);
		if (!sequence1_data.validity.RowIsValid(s1_i)) {
			throw InvalidInputException("ena_upload_reads: NULL sequence1 in sample '%s'", sample_ref);
		}
		const auto q1_i = qual1_data.sel->get_index(row);
		if (!qual1_data.validity.RowIsValid(q1_i)) {
			throw InvalidInputException("ena_upload_reads: NULL qual1 in sample '%s'", sample_ref);
		}

		const char *id = read_id_strs[rid_i].GetData();
		const std::size_t id_len = read_id_strs[rid_i].GetSize();
		const char *s1 = sequence1_strs[s1_i].GetData();
		const std::size_t s1_len = sequence1_strs[s1_i].GetSize();
		auto q1 = ExtractQualList(chunk.data[2], qual1_data, row);
		if (q1.size() != s1_len) {
			throw InvalidInputException(
			    "ena_upload_reads: qual1 length (%llu) does not match sequence1 length (%llu) in sample '%s'",
			    (uint64_t)q1.size(), (uint64_t)s1_len, sample_ref);
		}
		encoder.Encode(write0, id, id_len, nullptr, 0, s1, s1_len, q1.data(), q1.size());

		if (!need_r2) {
			continue;
		}
		const auto s2_i = sequence2_data.sel->get_index(row);
		const auto q2_i = qual2_data.sel->get_index(row);
		if (!sequence2_data.validity.RowIsValid(s2_i) || !qual2_data.validity.RowIsValid(q2_i)) {
			// PlanSamples resolved this sample as paired; a NULL R2 here means the
			// relation returned different rows between the pre-pass and this scan
			// (a non-deterministic view). Refuse rather than emit a broken pair.
			throw InvalidInputException("ena_upload_reads: sample '%s' resolved as paired but a row has NULL R2 "
			                            "(is the input relation non-deterministic across scans?)",
			                            sample_ref);
		}
		const char *s2 = sequence2_strs[s2_i].GetData();
		const std::size_t s2_len = sequence2_strs[s2_i].GetSize();
		auto q2 = ExtractQualList(chunk.data[4], qual2_data, row);
		if (q2.size() != s2_len) {
			throw InvalidInputException(
			    "ena_upload_reads: qual2 length (%llu) does not match sequence2 length (%llu) in sample '%s'",
			    (uint64_t)q2.size(), (uint64_t)s2_len, sample_ref);
		}
		encoder.Encode(write_r2, id, id_len, nullptr, 0, s2, s2_len, q2.data(), q2.size());
	}
}

#ifdef MIINT_HAS_CURL
// Upload an already-staged local file to a curl URL (ftp/ftps/http/https). The
// encode path stages every transport's output to a file (uniform); for curl we
// then stream that file to the socket via a file-backed read callback.
miint::CurlUploadResult UploadFileViaCurl(FileSystem &fs, const string &path, miint::CurlUploadOptions opts) {
	auto handle = fs.OpenFile(path, FileFlags::FILE_FLAGS_READ);
	opts.expected_size = static_cast<long long>(fs.GetFileSize(*handle));
	auto producer = [&fs, &handle](char *buf, std::size_t max_bytes) -> std::size_t {
		return static_cast<std::size_t>(fs.Read(*handle, buf, static_cast<int64_t>(max_bytes)));
	};
	return miint::RunCurlUpload(opts, producer);
}
#endif

// Encode + stage + transport one sample's file(s). Every transport stages the
// gzipped output to a file first — the destination for file://, a private temp
// dir for aspera/curl — then ships it: ascp for Aspera, a file-backed upload for
// curl, nothing more for file://. The sample's rows are streamed once (both
// sinks stay open across the single pass for split-paired), so peak memory is
// one DataChunk + one zlib window regardless of how large the sample is.
void UploadOneSample(ClientContext &context, const ENAUploadReadsBindData &bind, ENAUploadReadsGlobalState &gs,
                     Connection &conn, const string &data_query_prefix, const SamplePlan &plan) {
	auto &fs = FileSystem::GetFileSystem(context);
	const auto filenames = OutputFilenames(plan.sample_ref, plan.resolved_layout); // 1 (single/interleaved) or 2

	// On builds without Aspera support, fail-fast if Aspera transport was
	// requested — the ascp binary path and its config aren't compiled in.
#if !(MIINT_ASPERA_SUPPORTED && defined(MIINT_STATIC_BUILD))
	if (gs.target.transport == UploadTransport::ASPERA) {
		throw IOException("ena_upload_reads: Aspera transport is not supported on this build/platform");
	}
#endif

	const bool stage = gs.target.transport != UploadTransport::LOCAL_FILE;
	string temp_dir;
	vector<string> write_paths;
	if (!stage) {
		for (auto &fname : filenames) {
			write_paths.push_back(gs.target.remote_dir + fname);
		}
	} else {
		// Private temp dir so each staged file's basename matches its remote
		// name; removed after transport. mkdtemp(3) is POSIX-only (MinGW's libc
		// omits it), so build a uniquely-named directory through DuckDB's
		// FileSystem instead — portable across POSIX, Windows and Emscripten.
		// GetTempDirectory honours DuckDB's configured temp_directory and the
		// platform temp path; CreateDirectory is non-recursive, so ensure that
		// root exists first. A 64-bit random suffix makes collision with an
		// existing directory negligible; the bounded retry covers the off chance.
		const string temp_root = miint::RemoteFileHelper::GetTempDirectory(context);
		fs.CreateDirectoriesRecursive(temp_root);
		RandomEngine rng;
		for (idx_t attempt = 0;; attempt++) {
			string candidate = fs.JoinPath(
			    temp_root, StringUtil::Format("ena-upload-%016llx", (unsigned long long)rng.NextRandomInteger64()));
			if (!fs.DirectoryExists(candidate)) {
				fs.CreateDirectory(candidate);
				temp_dir = candidate;
				break;
			}
			if (attempt >= 64) {
				throw IOException("ena_upload_reads: could not create a unique staging directory under '%s'",
				                  temp_root);
			}
		}
		for (auto &fname : filenames) {
			write_paths.push_back(fs.JoinPath(temp_dir, fname));
		}
	}

	// Cleanup runs whether the block below succeeds or throws. For LOCAL_FILE
	// temp_dir is empty, so the destination is left in place (matches COPY).
	auto cleanup_temp = [&]() noexcept {
		if (!temp_dir.empty()) {
			for (auto &wp : write_paths) {
				(void)unlink(wp.c_str());
			}
			(void)rmdir(temp_dir.c_str());
		}
	};

	try {
		vector<unique_ptr<GzipMd5FileSink>> sinks;
		sinks.reserve(write_paths.size());
		for (auto &wp : write_paths) {
			sinks.push_back(make_uniq<GzipMd5FileSink>(fs, wp));
		}
		GzipMd5FileSink *sink1 = sinks.size() > 1 ? sinks[1].get() : nullptr;

		// Stream this sample's rows once and encode straight into the sink(s).
		// SendQuery (not a prepared statement) is deliberate: it defaults to a
		// streaming result, so Fetch pulls one DataChunk at a time and peak memory
		// stays bounded. A prepared statement's result output_type defaults to
		// FORCE_MATERIALIZED — it would buffer the entire sample in RAM, defeating
		// the whole refactor. Draining to exhaustion closes the stream before the
		// next sample's query. If EncodeChunk throws mid-stream the result is
		// abandoned un-closed, but the exception aborts the statement: `conn` is
		// torn down (its dtor runs cleanup) and never reused, so the dangling
		// active query never matters here.
		FastqEncoder encoder(bind.qual_offset);
		auto result = conn.SendQuery(data_query_prefix + KeywordHelper::WriteQuoted(plan.sample_ref, '\''));
		if (result->HasError()) {
			throw InvalidInputException("ena_upload_reads: failed to read sample '%s': %s", plan.sample_ref,
			                            result->GetError());
		}
		while (auto chunk = result->Fetch()) {
			if (chunk->size() == 0) {
				continue;
			}
			EncodeChunk(*chunk, plan.resolved_layout, encoder, *sinks[0], sink1, plan.sample_ref);
		}

		// Finish each file (flushes the gzip trailer, closes the handle, yields
		// the md5 over exactly the on-disk bytes), then transport it.
		for (size_t f = 0; f < filenames.size(); f++) {
			auto finished = sinks[f]->Finish();
			string md5_hex = std::move(finished.first);
			const uint64_t bytes_written = finished.second;

			if (gs.target.transport == UploadTransport::CURL) {
#ifdef MIINT_HAS_CURL
				miint::CurlUploadOptions opts;
				opts.url = gs.target.url_for_curl + filenames[f];
				opts.user = gs.user;
				opts.password = gs.password;
				opts.create_dirs = true;
				// `MIINT_CURL_VERBOSE=1` enables CURLOPT_VERBOSE wire trace on
				// stderr — invaluable for diagnosing FTPS/HTTPS upload hangs.
				if (const char *v = std::getenv("MIINT_CURL_VERBOSE"); v && std::string(v) == "1") {
					opts.verbose = true;
				}
				auto curl_result = UploadFileViaCurl(fs, write_paths[f], opts);
				if (!curl_result.error_message.empty()) {
					throw IOException("ena_upload_reads: %s upload failed for sample '%s' file '%s': %s",
					                  gs.target.scheme, plan.sample_ref, filenames[f], curl_result.error_message);
				}
#else
				throw IOException("ena_upload_reads: %s:// transport requires libcurl, which is disabled in this "
				                  "build (use aspera:// or file:// instead, or rebuild with -DMIINT_ENABLE_CURL=ON)",
				                  gs.target.scheme);
#endif
			}
#if MIINT_ASPERA_SUPPORTED && defined(MIINT_STATIC_BUILD)
			else if (gs.target.transport == UploadTransport::ASPERA) {
				AsperaSendOptions opts;
				opts.ascp_path = gs.aspera_ascp_path;
				opts.key_path = gs.aspera_key_path;
				opts.user = gs.user;
				opts.host = gs.target.host;
				opts.port = 33001;
				opts.local_path = write_paths[f];
				opts.remote_dir = gs.target.remote_dir;
				opts.max_rate = gs.aspera_max_rate;
				auto argv = BuildAscpSendArgv(opts);
				auto ascp_result = miint::RunAsperaSend(argv, gs.password);
				if (ascp_result.exit_code != 0) {
					throw IOException("ena_upload_reads: ascp failed (exit %d) for sample '%s' file '%s': %s",
					                  ascp_result.exit_code, plan.sample_ref, filenames[f], ascp_result.stderr_output);
				}
			}
#endif
			// LOCAL_FILE: the destination IS write_paths[f]; nothing more to do.

			EmittedFile e;
			e.sample_ref = plan.sample_ref;
			e.filename = filenames[f];
			e.filetype = "fastq";
			e.md5_hex = std::move(md5_hex);
			e.bytes_written = bytes_written;
			e.layout_name = FastqLayoutModeName(plan.resolved_layout);
			gs.emitted.push_back(std::move(e));
		}
	} catch (...) {
		cleanup_temp();
		throw;
	}
	cleanup_temp();
}

// Plan the samples (cheap aggregate pre-pass), then stream each one to its
// destination, one at a time. A single Connection is reused across samples; each
// sample's streaming result is fully drained before the next query (DuckDB
// allows only one active stream per connection).
void RunStreamingUpload(ClientContext &context, const ENAUploadReadsBindData &bind, ENAUploadReadsGlobalState &gs) {
	auto conn = MakeReadOnlyHelperConnection(context);

	const bool has_r2_columns = ValidateSchemaDetectR2(conn, bind.relation_name);
	PlanSamples(conn, bind.relation_name, bind.layout_mode, has_r2_columns, gs.samples);

	if (gs.target.transport == UploadTransport::LOCAL_FILE) {
		auto &fs = FileSystem::GetFileSystem(context);
		fs.CreateDirectoriesRecursive(gs.target.remote_dir);
	}

	// Per-sample data query, completed in UploadOneSample by appending a quoted
	// sample_ref literal. WriteQuoted escapes the value, so a sample_ref with an
	// embedded quote stays safe (the filename-level checks in ValidateSampleRef
	// reject '/' and '..' but allow quotes).
	const string data_query_prefix = "SELECT read_id, sequence1, qual1" +
	                                 string(has_r2_columns ? ", sequence2, qual2" : "") + " FROM " +
	                                 KeywordHelper::WriteOptionallyQuoted(bind.relation_name) + " WHERE sample_ref = ";
	for (auto &plan : gs.samples) {
		UploadOneSample(context, bind, gs, conn, data_query_prefix, plan);
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

	// Plan the samples (cheap aggregate pre-pass, fail-fast layout resolution)
	// and run the entire upload eagerly during init, streaming one sample at a
	// time. The execution side just drains the per-file rows from `emitted`. This
	// keeps the error model simple (any failure aborts the statement) and keeps
	// peak memory bounded to one DataChunk + one zlib window regardless of size.
	RunStreamingUpload(context, bind, *gs);
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
