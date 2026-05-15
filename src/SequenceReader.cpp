#include <SequenceReader.hpp>
#include <algorithm>

namespace miint {
// Helper function to extract base read ID by stripping /[1-9] suffix and comments
static std::string base_read_id(const std::string &id) {
	// First, strip comments (everything after first space)
	size_t space_pos = id.find(' ');
	std::string base = (space_pos == std::string::npos) ? id : id.substr(0, space_pos);

	// Check if it ends with /[1-9] (single digit 1-9 only)
	// Need at least 3 chars for pattern "x/1"
	if (base.length() >= 3) {
		size_t len = base.length();
		char last_char = base[len - 1];
		char second_last = base[len - 2];

		// If ends with /[1-9], strip that suffix
		if (second_last == '/' && last_char >= '1' && last_char <= '9') {
			base = base.substr(0, len - 2);
		}
	}

	return base;
}

// Helper function to check if two read IDs match after normalization
static void check_ids(const std::string &name1, const std::string &name2) {
	std::string base_id1 = base_read_id(name1);
	std::string base_id2 = base_read_id(name2);

	if (base_id1 != base_id2) {
		throw std::runtime_error("Mismatched read IDs: " + name1 + " vs " + name2);
	}
}

// Validate that paired-end files use the same format (both FASTA or both FASTQ)
static void validate_format_consistency(bool is_fasta1, bool is_fasta2, bool allow_format_mismatch) {
	if (!allow_format_mismatch && is_fasta1 != is_fasta2) {
		throw std::runtime_error("Cannot mix FASTA and FASTQ formats: sequence1 is " +
		                         std::string(is_fasta1 ? "FASTA" : "FASTQ") + ", sequence2 is " +
		                         std::string(is_fasta2 ? "FASTA" : "FASTQ"));
	}
}

// Read from whichever stream variant is active
std::vector<klibpp::KSeq> SequenceReader::read_stream(StreamVar &var, int n) {
	return std::visit([n](auto &reader) { return reader->read(static_cast<std::vector<klibpp::KSeq>::size_type>(n)); },
	                  var);
}

SequenceReader::SequenceReader(const std::string &path1, const std::optional<std::string> &path2,
                               bool allow_format_mismatch)
    : first_read_(true) {
	sequence1_reader_ = std::make_unique<SeqStreamIn>(path1.c_str());

	// Check if first file is empty by attempting to peek at first record
	buffered_read1_ = read_stream(sequence1_reader_, 1);
	bool is_empty1 = buffered_read1_.empty();
	bool is_fasta1 = !is_empty1 && buffered_read1_[0].qual.empty();

	if (is_empty1) {
		throw std::runtime_error("Empty file: " + path1);
	}

	// Keep buffered_read1_ to return on first read() call instead of recreating reader

	paired_ = path2.has_value() && (path2->length() > 0);
	if (paired_) {
		sequence2_reader_.emplace(std::make_unique<SeqStreamIn>(path2->c_str()));

		// Check if second file is empty and detect format
		buffered_read2_ = read_stream(sequence2_reader_.value(), 1);
		bool is_empty2 = buffered_read2_.empty();
		bool is_fasta2 = !is_empty2 && buffered_read2_[0].qual.empty();

		if (is_empty2) {
			throw std::runtime_error("Empty file: " + path2.value());
		}

		validate_format_consistency(is_fasta1, is_fasta2, allow_format_mismatch);

		// Keep buffered_read2_ to return on first read() call instead of recreating reader
	}
}

#ifdef MIINT_STATIC_BUILD
SequenceReader::SequenceReader(DuckDBSeqStreamHandle stream1, DuckDBSeqStreamHandle stream2_or_null,
                               bool allow_format_mismatch)
    : first_read_(true) {
	// Ownership protocol: each input handle owns its raw pointer until the
	// matching make_unique<KStreamIn>() succeeds. Holding the kstream in a
	// local before `release()` and the variant move-assign closes any window
	// where two owners could exist: if make_unique throws, the input handle
	// still owns the stream and frees it on unwind; if make_unique succeeds,
	// `release()` runs before any variant op, and the variant move-assign
	// from rvalue unique_ptr is noexcept. On any later peek throwing, the
	// kstream wrapper (now inside the variant) is the sole owner.
	auto kstream1 = std::make_unique<DuckDBSeqStreamIn>(stream1.get(), duckdb_seq_read, duckdb_seq_close);
	stream1.release();
	sequence1_reader_ = std::move(kstream1);

	buffered_read1_ = read_stream(sequence1_reader_, 1);
	bool is_empty1 = buffered_read1_.empty();
	bool is_fasta1 = !is_empty1 && buffered_read1_[0].qual.empty();

	if (is_empty1) {
		throw std::runtime_error("Empty stream");
	}

	paired_ = (stream2_or_null.get() != nullptr);
	if (paired_) {
		auto kstream2 = std::make_unique<DuckDBSeqStreamIn>(stream2_or_null.get(), duckdb_seq_read, duckdb_seq_close);
		stream2_or_null.release();
		sequence2_reader_.emplace(std::move(kstream2));

		buffered_read2_ = read_stream(sequence2_reader_.value(), 1);
		bool is_empty2 = buffered_read2_.empty();
		bool is_fasta2 = !is_empty2 && buffered_read2_[0].qual.empty();

		if (is_empty2) {
			throw std::runtime_error("Empty stream (sequence2)");
		}

		validate_format_consistency(is_fasta1, is_fasta2, allow_format_mismatch);
	}
}

#if MIINT_ASPERA_SUPPORTED
SequenceReader::SequenceReader(AsperaSeqStreamHandle stream1, AsperaSeqStreamHandle stream2_or_null,
                               bool allow_format_mismatch)
    : first_read_(true) {
	// See ownership protocol on the DuckDB+DuckDB ctor above. Unlike that
	// ctor, both kstream wrappers are constructed up front before either
	// peek: Aspera streams are backed by live ascp pipes, and ordering the
	// transfers/peeks per-stream would interleave reads from the two
	// processes in ways the original (pre-handle) code took care to avoid.
	auto kstream1 = std::make_unique<AsperaSeqStreamIn>(stream1.get(), aspera_seq_read, aspera_seq_close);
	stream1.release();
	sequence1_reader_ = std::move(kstream1);

	paired_ = (stream2_or_null.get() != nullptr);
	if (paired_) {
		auto kstream2 = std::make_unique<AsperaSeqStreamIn>(stream2_or_null.get(), aspera_seq_read, aspera_seq_close);
		stream2_or_null.release();
		sequence2_reader_.emplace(std::move(kstream2));
	}

	buffered_read1_ = read_stream(sequence1_reader_, 1);
	if (buffered_read1_.empty()) {
		throw std::runtime_error("Empty stream");
	}
	bool is_fasta1 = buffered_read1_[0].qual.empty();

	if (paired_) {
		buffered_read2_ = read_stream(sequence2_reader_.value(), 1);
		if (buffered_read2_.empty()) {
			throw std::runtime_error("Empty stream (sequence2)");
		}
		bool is_fasta2 = buffered_read2_[0].qual.empty();

		validate_format_consistency(is_fasta1, is_fasta2, allow_format_mismatch);
	}
}
#endif // MIINT_ASPERA_SUPPORTED
#endif // MIINT_STATIC_BUILD

// Records polled per stream call. Larger values amortize per-call overhead (vector
// allocation) across more records; smaller values cap the transient memory overshoot
// when a single poll pulls in records that together exceed max_bytes. 16 is a
// compromise: for short-read FASTQ a 2048-record chunk needs ~128 polls instead of
// 2048 (cuts ~94% of per-record allocation overhead), while for bacterial genomes
// (~5 MB/record) the worst-case overshoot is 16 * 5 MB = 80 MB -- negligible next to
// a 512 MB default budget. Any tail records beyond the stop point are pushed back
// into buffered_read1_ / buffered_read2_ for the next call, so nothing is dropped.
static constexpr int POLL_CHUNK_SIZE = 16;

// Bytes contributed by a single record. Counted against the per-chunk budget so that
// FASTQ quality (and paired mates) are accounted for alongside the sequence itself.
static size_t record_bytes(const klibpp::KSeq &rec) {
	return rec.name.size() + rec.comment.size() + rec.seq.size() + rec.qual.size();
}

// Append a KSeq into an SE batch, stripping sequence whitespace in place.
static void append_se(SequenceRecordBatch &batch, klibpp::KSeq &rec) {
	rec.seq.erase(std::remove_if(rec.seq.begin(), rec.seq.end(), [](unsigned char c) { return std::isspace(c); }),
	              rec.seq.end());
	batch.read_ids.emplace_back(base_read_id(rec.name));
	batch.comments.emplace_back(rec.comment);
	batch.sequences1.emplace_back(std::move(rec.seq));
	batch.quals1.emplace_back(rec.qual);
}

SequenceRecordBatch SequenceReader::read_se(const int n, const size_t max_bytes) {
	SequenceRecordBatch batch(false);
	batch.reserve(n);

	size_t cumulative_bytes = 0;

	// Drain a polled chunk into `batch`, pushing any unused tail into buffered_read1_ when
	// a stop condition (budget hit or batch full) fires mid-chunk. Returns true when the
	// whole chunk was consumed and caller should poll again, false otherwise.
	auto drain_chunk = [&](std::vector<klibpp::KSeq> &chunk) -> bool {
		for (size_t i = 0; i < chunk.size(); i++) {
			auto &rec = chunk[i];
			const size_t rec_bytes = record_bytes(rec);
			// Starvation guard: if batch is still empty, accept at least one record even
			// when it alone exceeds max_bytes -- otherwise the pipeline stalls.
			if (!batch.empty() && cumulative_bytes + rec_bytes > max_bytes) {
				buffered_read1_.clear();
				buffered_read1_.reserve(chunk.size() - i);
				for (size_t j = i; j < chunk.size(); j++) {
					buffered_read1_.push_back(std::move(chunk[j]));
				}
				first_read_ = true;
				return false;
			}
			append_se(batch, rec);
			cumulative_bytes += rec_bytes;
			const bool done = (int)batch.size() >= n || cumulative_bytes >= max_bytes;
			if (done) {
				if (i + 1 < chunk.size()) {
					buffered_read1_.clear();
					buffered_read1_.reserve(chunk.size() - (i + 1));
					for (size_t j = i + 1; j < chunk.size(); j++) {
						buffered_read1_.push_back(std::move(chunk[j]));
					}
					first_read_ = true;
				}
				return false;
			}
		}
		return true;
	};

	// Replay pre-buffered records first (peeked at construction, or tail from a prior call).
	if (first_read_) {
		first_read_ = false;
		auto buffered = std::move(buffered_read1_);
		buffered_read1_.clear();
		if (!drain_chunk(buffered)) {
			return batch;
		}
	}

	// Poll the stream in chunks so the budget can short-circuit before the next I/O round.
	while ((int)batch.size() < n && cumulative_bytes < max_bytes) {
		const int remaining = n - (int)batch.size();
		const int poll_n = remaining < POLL_CHUNK_SIZE ? remaining : POLL_CHUNK_SIZE;
		auto more = read_stream(sequence1_reader_, poll_n);
		if (more.empty()) {
			break;
		}
		if (!drain_chunk(more)) {
			break;
		}
	}

	return batch;
}

// Append a paired KSeq into a PE batch, stripping sequence whitespace on both mates.
static void append_pe(SequenceRecordBatch &batch, klibpp::KSeq &rec1, klibpp::KSeq &rec2) {
	rec1.seq.erase(std::remove_if(rec1.seq.begin(), rec1.seq.end(), [](unsigned char c) { return std::isspace(c); }),
	               rec1.seq.end());
	rec2.seq.erase(std::remove_if(rec2.seq.begin(), rec2.seq.end(), [](unsigned char c) { return std::isspace(c); }),
	               rec2.seq.end());

	batch.read_ids.emplace_back(base_read_id(rec1.name));
	batch.comments.emplace_back(rec1.comment);
	batch.sequences1.emplace_back(std::move(rec1.seq));
	batch.sequences2.emplace_back(std::move(rec2.seq));
	batch.quals1.emplace_back(rec1.qual);
	batch.quals2.emplace_back(rec2.qual);
}

SequenceRecordBatch SequenceReader::read_pe(const int n, const size_t max_bytes) {
	SequenceRecordBatch batch(true);
	batch.reserve(n);

	size_t cumulative_bytes = 0;

	// Push [i..end) of both chunks into buffered_read1_ / buffered_read2_ so the next call
	// replays them. Kept in a helper to keep mate buffers in lockstep.
	auto push_back_tail = [&](std::vector<klibpp::KSeq> &c1, std::vector<klibpp::KSeq> &c2, size_t i) {
		buffered_read1_.clear();
		buffered_read2_.clear();
		const size_t tail = c1.size() - i;
		buffered_read1_.reserve(tail);
		buffered_read2_.reserve(tail);
		for (size_t j = i; j < c1.size(); j++) {
			buffered_read1_.push_back(std::move(c1[j]));
			buffered_read2_.push_back(std::move(c2[j]));
		}
		first_read_ = true;
	};

	// Drain paired chunks into `batch`, validating mate-id parity and enforcing the budget.
	auto drain_chunks = [&](std::vector<klibpp::KSeq> &c1, std::vector<klibpp::KSeq> &c2) -> bool {
		if (c1.size() != c2.size()) {
			const std::string &id = !c1.empty() ? c1.back().name : c2.back().name;
			throw std::runtime_error("Mismatched number of records: missing mate for " + id);
		}
		for (size_t i = 0; i < c1.size(); i++) {
			auto &rec1 = c1[i];
			auto &rec2 = c2[i];
			check_ids(rec1.name, rec2.name);
			const size_t pair_bytes = record_bytes(rec1) + record_bytes(rec2);
			if (!batch.empty() && cumulative_bytes + pair_bytes > max_bytes) {
				push_back_tail(c1, c2, i);
				return false;
			}
			append_pe(batch, rec1, rec2);
			cumulative_bytes += pair_bytes;
			const bool done = (int)batch.size() >= n || cumulative_bytes >= max_bytes;
			if (done) {
				if (i + 1 < c1.size()) {
					push_back_tail(c1, c2, i + 1);
				}
				return false;
			}
		}
		return true;
	};

	if (first_read_) {
		first_read_ = false;
		auto buf1 = std::move(buffered_read1_);
		auto buf2 = std::move(buffered_read2_);
		buffered_read1_.clear();
		buffered_read2_.clear();

		if (buf1.size() != buf2.size()) {
			// Should not happen -- peek-at-construction and push_back_tail keep them in sync.
			throw std::runtime_error("Buffered paired records out of sync (internal error)");
		}
		if (!drain_chunks(buf1, buf2)) {
			return batch;
		}
	}

	while ((int)batch.size() < n && cumulative_bytes < max_bytes) {
		const int remaining = n - (int)batch.size();
		const int poll_n = remaining < POLL_CHUNK_SIZE ? remaining : POLL_CHUNK_SIZE;
		auto more1 = read_stream(sequence1_reader_, poll_n);
		auto more2 = read_stream(sequence2_reader_.value(), poll_n);

		if (more1.empty() && more2.empty()) {
			break;
		}
		if (!drain_chunks(more1, more2)) {
			break;
		}
	}

	return batch;
}

SequenceRecordBatch SequenceReader::read(const int n, const size_t max_bytes) {
	if (paired_) {
		return read_pe(n, max_bytes);
	} else {
		return read_se(n, max_bytes);
	}
}
}; // namespace miint
