#pragma once
#include <climits>
#include <cstddef>
#include <optional>
#include <variant>
#include <vector>
#include <string>
#include <memory>
#include <kseq++/seqio.hpp>
#include "SequenceRecord.hpp"

#ifdef MIINT_STATIC_BUILD
#include "duckdb_seq_stream.hpp"
#include "aspera_stream.hpp"
#endif

namespace miint {

// A local-file sequence stream that checks gzopen() actually succeeded.
//
// klibpp::SeqStreamIn hands gzopen()'s result straight to KStream without looking at it.
// A NULL handle does not fail loudly: it produces a stream that yields zero records, which
// SequenceReader cannot tell apart from a genuinely empty file, so every failed open --
// EACCES, EMFILE/ENFILE on descriptor exhaustion, ENOMEM -- was reported as
// "Empty file: <path>" naming a file that is fine. ext/kseq++ is vendored unpatched, so the
// check lives here instead. Deriving from SeqStreamIn's own base lets us open the file
// ourselves and read errno with nothing running in between; SeqStreamIn adds only
// constructors, so behaviour is otherwise identical.
class CheckedSeqStreamIn : public klibpp::SeqStreamIn::base_type {
public:
	explicit CheckedSeqStreamIn(const std::string &path)
	    : klibpp::SeqStreamIn::base_type(OpenOrThrow(path), gzread, gzclose), path_(path) {
	}

	// zlib's account of how the input ended, empty if it ended cleanly. A gzip file that is
	// cut short does not make gzread fail: it returns what it could inflate and then 0, as at
	// a real end of file, and only gzerror says the data stopped early.
	std::string EndOfInputError() const;

private:
	static gzFile OpenOrThrow(const std::string &path);

	std::string path_;
};

class SequenceReader {
public:
	//! With `read_ahead`, each file is read ahead of read() on a thread of its own, so that
	//! inflating and parsing it overlap the caller's work. It stops reading ahead once it holds
	//! this many bytes of records per file, or read()'s own byte budget if that is smaller. Each
	//! poll asks for as many records as the last one read says fit in the room left (at least
	//! one, at most sixteen), so a poll passes it by about a record, or by up to sixteen where
	//! record lengths jump. Only regular files are read ahead: a read from a pipe can wait on
	//! its writer forever, and the reader could then never be dropped.
	static constexpr size_t READ_AHEAD_BYTES = 4 * 1024 * 1024;

	//! Whether a scan of `files` files (pairs of them, if `paired`), all read at once, can read
	//! each ahead and stay within `threads`: a file being read takes one thread, and its
	//! read-ahead one more per file of the pair.
	static bool ReadAheadFits(size_t files, bool paired, size_t threads) {
		return files * (paired ? 3 : 2) <= threads;
	}

	// An empty input is an error unless `allow_empty`, which reads it as no records. Either
	// way, a pair with exactly one empty mate is an error: its other mate's records have no mates.
	explicit SequenceReader(const std::string &path1, const std::optional<std::string> &path2 = std::nullopt,
	                        bool allow_format_mismatch = false, bool allow_empty = false, bool read_ahead = false);
	~SequenceReader();
	// A read-ahead holds a reference to the stream it reads, which lives in this object: a copy
	// or a move would leave it reading the old one.
	SequenceReader(const SequenceReader &) = delete;
	SequenceReader &operator=(const SequenceReader &) = delete;
	SequenceReader(SequenceReader &&) = delete;
	SequenceReader &operator=(SequenceReader &&) = delete;

#ifdef MIINT_STATIC_BUILD
	// Streaming constructor for remote files via DuckDB FileHandle. Ownership of
	// the streams is transferred via DuckDBSeqStreamHandle: pass a nullptr-valued
	// handle (still constructed with the duckdb_seq_close deleter) for the
	// second slot when the run is single-end.
	SequenceReader(DuckDBSeqStreamHandle stream1, DuckDBSeqStreamHandle stream2_or_null,
	               bool allow_format_mismatch = false, bool allow_empty = false);

#if MIINT_ASPERA_SUPPORTED
	// Streaming constructor for Aspera pipe-backed streams.
	SequenceReader(AsperaSeqStreamHandle stream1, AsperaSeqStreamHandle stream2_or_null,
	               bool allow_format_mismatch = false);
#endif // MIINT_ASPERA_SUPPORTED
#endif // MIINT_STATIC_BUILD

	// Read up to `n` records or until cumulative bytes across the accumulated batch exceed
	// `max_bytes`, whichever comes first. Byte accounting sums `name + comment + seq + qual`
	// across both mates in paired-end mode. At least one record is always returned if the
	// stream has data, even when a single record alone exceeds `max_bytes` (starvation guard).
	// Defaults to SIZE_MAX, preserving the row-only behavior for existing callers.
	SequenceRecordBatch read(const int n, const size_t max_bytes = SIZE_MAX);

	// Whether sequence1 (`segment` 0) or sequence2 (1) is FASTQ: its first header starts with '@'.
	// A FASTQ record with no bases has an empty quality, as a FASTA record does, so a record's
	// quality cannot tell the formats apart.
	bool IsFastq(int segment) const {
		return segment == 0 ? state1_.fastq : state2_.fastq;
	}

	// Test seam: the largest number of records materialized by any single underlying stream
	// poll over this reader's lifetime, by read() or by a read-ahead. Lets tests assert that
	// large records do not trigger oversized prefetches (see dynamic poll sizing in
	// read_se/read_pe and in ReadAhead::ReadUntilStopped).
	size_t MaxPollCount() const;

	// Test seam: the bytes of records read ahead of read() and not yet taken, for sequence1
	// (segment 0) or sequence2 (segment 1).
	size_t ReadAheadBytes(int segment) const;

	// Test seam: whether the read-ahead of `segment` has stopped, and only taking records from
	// it will start it again: it holds READ_AHEAD_BYTES (or read()'s budget), or it reached the
	// end of the file or failed. A segment that reads nothing ahead is idle.
	bool ReadAheadIdle(int segment) const;

	// Test seam: whether `segment` is read ahead, on the thread its first poll started.
	bool ReadAheadStarted(int segment) const;

private:
	using SeqStreamIn = CheckedSeqStreamIn;

#if defined(MIINT_STATIC_BUILD) && MIINT_ASPERA_SUPPORTED
	using StreamVar = std::variant<std::unique_ptr<SeqStreamIn>, std::unique_ptr<DuckDBSeqStreamIn>,
	                               std::unique_ptr<AsperaSeqStreamIn>>;
#elif defined(MIINT_STATIC_BUILD)
	using StreamVar = std::variant<std::unique_ptr<SeqStreamIn>, std::unique_ptr<DuckDBSeqStreamIn>>;
#else
	using StreamVar = std::variant<std::unique_ptr<SeqStreamIn>>;
#endif

	// Where the record parser is in one stream.
	struct StreamState {
		uint64_t records_read = 0;
		// The next header's '>' or '@' when the previous record's sequence lines ended by
		// consuming it, otherwise 0.
		char pending_header = 0;
		// Nothing parsed yet. Text before the first header is skipped for FASTA, which can open
		// with ';' comments, and is an error for FASTQ.
		bool at_start = true;
		// The first header is '@'. Set by the first parse, which the constructor runs, so before
		// any read-ahead thread touches this state.
		bool fastq = false;
	};

	enum class ParseResult { RECORD, END, MALFORMED };

	// Parse one record. On MALFORMED, `problem` says what is wrong with it.
	template <class Stream>
	static ParseResult parse_record(Stream &stream, klibpp::KSeq &rec, StreamState &state, std::string &problem);

	// Read up to `n` records. Throws, naming `source` and the record, when the stream fails
	// rather than ends or a record is malformed.
	static std::vector<klibpp::KSeq> read_stream(StreamVar &var, int n, const std::string &source, StreamState &state);

	class ReadAhead;

	StreamVar sequence1_reader_;
	std::optional<StreamVar> sequence2_reader_;
	// What each stream is called in errors: its path where there is one.
	std::string source1_ = "sequence1";
	std::string source2_ = "sequence2";
	StreamState state1_;
	StreamState state2_;
	// Declared after the streams, names and states they read, so they stop before those go.
	bool read_ahead_enabled_ = false;
	size_t peeked_record_bytes_[2] = {0, 0};
	bool read_ahead_unavailable_[2] = {false, false}; // its thread could not be started
	std::unique_ptr<ReadAhead> read_ahead_[2];

	bool paired_;
	bool first_read_; // Track if we need to return buffered data
	std::vector<klibpp::KSeq> buffered_read1_;
	std::vector<klibpp::KSeq> buffered_read2_;

	// Running estimate of bytes per record (per pair in PE), used to size each poll so a
	// single read_stream call never materializes far more than the byte budget. 0 until the
	// first record is accepted.
	size_t observed_record_bytes_ = 0;
	size_t max_poll_count_ = 0; // test seam; see MaxPollCount()

	SequenceRecordBatch read_se(const int n, const size_t max_bytes);
	SequenceRecordBatch read_pe(const int n, const size_t max_bytes);
	// Up to `n` records of `segment`, fewer only at its end: read here, or taken from its
	// read-ahead, which the first poll starts.
	std::vector<klibpp::KSeq> poll(int segment, int n, size_t max_bytes);
};
}; // namespace miint
