#include <SequenceReader.hpp>

namespace miint {
SequenceReader::SequenceReader(const std::string &path1, const std::optional<std::string> &path2)
    : first_read_(true) {
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

std::vector<SequenceRecord> SequenceReader::read_se(const int n) {
	std::vector<SequenceRecord> out;
	out.reserve(n);

	// On first read, return buffered records first
	if (first_read_) {
		first_read_ = false;
		for (auto &rec : buffered_read1_) {
			out.emplace_back(rec);
		}
		buffered_read1_.clear();

		// If we got fewer than n records from buffer, read more
		int remaining = n - out.size();
		if (remaining > 0) {
			for (auto rec : sequence1_reader_->read(remaining)) {
				out.emplace_back(rec);
			}
		}
	} else {
		for (auto rec : sequence1_reader_->read(n)) {
			out.emplace_back(rec);
		}
	}

	return out;
}

std::vector<SequenceRecord> SequenceReader::read_pe(const int n) {
	std::vector<SequenceRecord> out;
	out.reserve(n);

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
		int remaining = n - read1s.size();
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

	for (size_t i = 0; i < read1s.size(); i++) {
		auto &rec1 = read1s[i];
		auto &rec2 = read2s[i];
		out.emplace_back(rec1, rec2);
	}

	return out;
}

std::vector<SequenceRecord> SequenceReader::read(const int n) {
	if (paired_) {
		return read_pe(n);
	} else {
		return read_se(n);
	}
}
}; // namespace miint
