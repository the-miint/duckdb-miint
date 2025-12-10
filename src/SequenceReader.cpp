#include <SequenceReader.hpp>

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

SequenceReader::SequenceReader(const std::string &path1, const std::optional<std::string> &path2) : first_read_(true) {
	sequence1_reader_ = std::make_unique<SeqStreamIn>(path1.c_str());

	// Check if first file is empty by attempting to peek at first record
	buffered_read1_ = sequence1_reader_->read(1);
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
		buffered_read2_ = sequence2_reader_.value()->read(1);
		bool is_empty2 = buffered_read2_.empty();
		bool is_fasta2 = !is_empty2 && buffered_read2_[0].qual.empty();

		if (is_empty2) {
			throw std::runtime_error("Empty file: " + path2.value());
		}

		// Validate format consistency
		if (is_fasta1 != is_fasta2) {
			throw std::runtime_error("Cannot mix FASTA and FASTQ formats: sequence1 is " +
			                         std::string(is_fasta1 ? "FASTA" : "FASTQ") + ", sequence2 is " +
			                         std::string(is_fasta2 ? "FASTA" : "FASTQ"));
		}

		// Keep buffered_read2_ to return on first read() call instead of recreating reader
	}
}

SequenceRecordBatch SequenceReader::read_se(const int n) {
	SequenceRecordBatch batch(false);
	batch.reserve(n);

	std::vector<klibpp::KSeq> reads;

	// On first read, return buffered records first
	if (first_read_) {
		first_read_ = false;
		reads = std::move(buffered_read1_);

		// If we got fewer than n records from buffer, read more
		int remaining = n - (int)reads.size();
		if (remaining > 0) {
			auto more = sequence1_reader_->read(remaining);
			reads.insert(reads.end(), more.begin(), more.end());
		}
	} else {
		reads = sequence1_reader_->read(n);
	}

	// Populate batch directly from KSeq records
	for (auto &rec : reads) {
		batch.read_ids.emplace_back(base_read_id(rec.name));
		batch.comments.emplace_back(rec.comment);
		batch.sequences1.emplace_back(rec.seq);
		batch.quals1.emplace_back(rec.qual);
	}

	return batch;
}

SequenceRecordBatch SequenceReader::read_pe(const int n) {
	SequenceRecordBatch batch(true);
	batch.reserve(n);

	// Read sequence1 and sequence2 sequentially
	// DuckDB's parallel execution (MaxThreads) provides better parallelism
	// than internal threading with small batch sizes
	std::vector<klibpp::KSeq> read1s, read2s;

	// On first read, start with buffered records
	if (first_read_) {
		first_read_ = false;
		read1s = std::move(buffered_read1_);
		read2s = std::move(buffered_read2_);

		// If we got fewer than n records from buffer, read more
		int remaining = n - (int)read1s.size();
		if (remaining > 0) {
			auto more1 = sequence1_reader_->read(remaining);
			auto more2 = sequence2_reader_.value()->read(remaining);
			read1s.insert(read1s.end(), more1.begin(), more1.end());
			read2s.insert(read2s.end(), more2.begin(), more2.end());
		}
	} else {
		read1s = sequence1_reader_->read(n);
		read2s = sequence2_reader_.value()->read(n);
	}

	if (read1s.size() != read2s.size()) {
		auto id1 = read1s.back().name;
		throw std::runtime_error("Mismatched number of records: missing mate for " + id1);
	}

	// Populate batch directly from KSeq records
	for (size_t i = 0; i < read1s.size(); i++) {
		auto &rec1 = read1s[i];
		auto &rec2 = read2s[i];

		// Validate that read IDs match
		check_ids(rec1.name, rec2.name);

		batch.read_ids.emplace_back(base_read_id(rec1.name));
		batch.comments.emplace_back(rec1.comment);
		batch.sequences1.emplace_back(rec1.seq);
		batch.sequences2.emplace_back(rec2.seq);
		batch.quals1.emplace_back(rec1.qual);
		batch.quals2.emplace_back(rec2.qual);
	}

	return batch;
}

SequenceRecordBatch SequenceReader::read(const int n) {
	if (paired_) {
		return read_pe(n);
	} else {
		return read_se(n);
	}
}
}; // namespace miint
