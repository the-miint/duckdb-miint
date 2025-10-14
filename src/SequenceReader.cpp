#include <SequenceReader.hpp>

namespace miint {
SequenceReader::SequenceReader(const std::string &path1, const std::optional<std::string> &path2) {
	sequence1_reader_ = std::make_unique<SeqStreamIn>(path1.c_str());

	// Check if first file is empty by attempting to peek at first record
	auto peek1 = sequence1_reader_->read(1);
	bool is_empty1 = peek1.empty();
	bool is_fasta1 = !is_empty1 && peek1[0].qual.empty();

	if (is_empty1) {
		throw std::runtime_error("Empty file: " + path1);
	}

	// Recreate reader after peek (kseq++ consumed the first record)
	sequence1_reader_ = std::make_unique<SeqStreamIn>(path1.c_str());

	paired_ = path2.has_value() && (path2->length() > 0);
	if (paired_) {
		sequence2_reader_.emplace(std::make_unique<SeqStreamIn>(path2->c_str()));

		// Check if second file is empty and detect format
		auto peek2 = sequence2_reader_.value()->read(1);
		bool is_empty2 = peek2.empty();
		bool is_fasta2 = !is_empty2 && peek2[0].qual.empty();

		if (is_empty2) {
			throw std::runtime_error("Empty file: " + path2.value());
		}

		// Validate format consistency
		if (is_fasta1 != is_fasta2) {
			throw std::runtime_error("Cannot mix FASTA and FASTQ formats: sequence1 is " +
			                         std::string(is_fasta1 ? "FASTA" : "FASTQ") + ", sequence2 is " +
			                         std::string(is_fasta2 ? "FASTA" : "FASTQ"));
		}

		// Recreate reader after peek
		sequence2_reader_.emplace(std::make_unique<SeqStreamIn>(path2->c_str()));
	}
}

std::vector<SequenceRecord> SequenceReader::read_se(const int n) const {
	std::vector<SequenceRecord> out;
	out.reserve(n);

	for (auto rec : sequence1_reader_->read(n)) {
		out.emplace_back(rec);
	}

	return std::move(out);
}

std::vector<SequenceRecord> SequenceReader::read_pe(const int n) const {
	std::vector<SequenceRecord> out;
	out.reserve(n);

	// Read sequence1 and sequence2 sequentially
	// DuckDB's parallel execution (MaxThreads) provides better parallelism
	// than internal threading with small batch sizes
	std::vector<klibpp::KSeq> read1s, read2s;

	read1s = sequence1_reader_->read(n);
	read2s = sequence2_reader_.value()->read(n);

	if (read1s.size() != read2s.size()) {
		auto id1 = read1s.back().name;
		throw std::runtime_error("Mismatched number of records: missing mate for " + id1);
	}

	for (size_t i = 0; i < read1s.size(); i++) {
		auto &rec1 = read1s[i];
		auto &rec2 = read2s[i];
		out.emplace_back(rec1, rec2);
	}

	return std::move(out);
}

std::vector<SequenceRecord> SequenceReader::read(const int n) const {
	if (paired_) {
		return std::move(read_pe(n));
	} else {
		return std::move(read_se(n));
	}
}
}; // namespace miint
