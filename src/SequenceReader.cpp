#include <SequenceReader.hpp>
#include <thread>

namespace miint {
SequenceReader::SequenceReader(const std::string &path1, const std::optional<std::string> &path2) {
	reader1_ = std::make_unique<SeqStreamIn>(path1.c_str());

	paired_ = path2.has_value() && (path2->length() > 0);
	if (paired_) {
		reader2_.emplace(std::make_unique<SeqStreamIn>(path2->c_str()));
	}
}

std::vector<SequenceRecord> SequenceReader::read_se(const int n) const {
	std::vector<SequenceRecord> out;
	out.reserve(n);

	for (auto rec : reader1_->read(n)) {
		out.emplace_back(rec);
	}

	return std::move(out);
}

std::vector<SequenceRecord> SequenceReader::read_pe(const int n) const {
	std::vector<SequenceRecord> out;
	out.reserve(n);

	// this _works_ but unclear if we really get benefit
	// it may be as "n" tends to be low, should consider buffering
	// TODO: read large blocks, then "read" from there
	std::vector<klibpp::KSeq> read1s, read2s;

	std::thread t1([&] { read1s = reader1_->read(n); });
	std::thread t2([&] { read2s = reader2_.value()->read(n); });

	// Block until both threads complete
	t1.join();
	t2.join();

	// read1s = reader1_->read(n);
	// read2s = reader2_.value()->read(n);

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
