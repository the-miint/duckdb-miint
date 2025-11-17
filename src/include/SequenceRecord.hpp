#pragma once
#include <optional>
#include <string>
#include "QualScore.hpp"
#include "kseq++/kseq++.hpp"
#include <stdexcept>
#include <kseq++/seqio.hpp>

namespace miint {
class UnpairedError : public std::runtime_error {
public:
	explicit UnpairedError(const std::string &msg) : std::runtime_error(msg) {
	}
};

enum class SequenceRecordField { READ_ID = 0, COMMENT, SEQUENCE1, SEQUENCE2, QUAL1, QUAL2 };

class SequenceRecord {
public:
	std::string read_id;
	std::string comment;
	std::string read1;
	std::optional<std::string> read2;
	QualScore qual1;
	std::optional<QualScore> qual2;
	bool is_paired;

	// qual as string
	// note: we cannot assert the read IDs of r1 and r2 are the same through this
	// method of construction
	SequenceRecord(const std::string &id, const std::string &cmt, const std::string &r1, const std::string &q1,
	               std::optional<std::string> r2 = std::nullopt, std::optional<std::string> q2 = std::nullopt) noexcept
	    : read_id(id), comment(cmt), read1(r1), qual1(q1), is_paired(r2.has_value() && q2.has_value()),
	      read2(std::move(r2)), qual2(std::move(q2)) {
	}

	// qual as QualScore
	// note: we cannot assert the read IDs of r1 and r2 are the same through this
	// method of construction
	SequenceRecord(const std::string &id, const std::string &cmt, const std::string &r1, const QualScore &q1,
	               std::optional<std::string> r2 = std::nullopt, std::optional<QualScore> q2 = std::nullopt) noexcept
	    : read_id(id), comment(cmt), read1(r1), qual1(q1), is_paired(r2.has_value() && q2.has_value()),
	      read2(std::move(r2)), qual2(std::move(q2)) {
	}

	explicit SequenceRecord(klibpp::KSeq &r1) noexcept
	    : comment(r1.comment), read1(r1.seq), qual1(r1.qual), is_paired(false) {
		read_id = base_read_id(r1.name);
	}

	SequenceRecord(klibpp::KSeq &r1, klibpp::KSeq &r2)
	    : comment(r1.comment), read1(r1.seq), qual1(r1.qual), read2(r2.seq), qual2(r2.qual), is_paired(true) {
		check_ids(r1, r2);
		read_id = base_read_id(r1.name);
	}

	const std::string &GetString(const SequenceRecordField field) const {
		switch (field) {
		case SequenceRecordField::READ_ID:
			return read_id;
		case SequenceRecordField::COMMENT:
			return comment;
		case SequenceRecordField::SEQUENCE1:
			return read1;
		case SequenceRecordField::SEQUENCE2:
			if (is_paired) {
				return read2.value();
			}
			throw UnpairedError("Attempted to access paired read with an unpaired record");
		case SequenceRecordField::QUAL1:
		case SequenceRecordField::QUAL2:
		default:
			throw std::invalid_argument("Invalid field");
		}
	}

	const QualScore &GetQual(const SequenceRecordField field) const {
		switch (field) {
		case SequenceRecordField::QUAL1:
			return qual1;
		case SequenceRecordField::QUAL2:
			if (is_paired) {
				return qual2.value();
			}
			throw UnpairedError("Attempted to access paired read with an unpaired record");
		case SequenceRecordField::READ_ID:
		case SequenceRecordField::COMMENT:
		case SequenceRecordField::SEQUENCE1:
		case SequenceRecordField::SEQUENCE2:
		default:
			throw std::invalid_argument("Invalid field");
		}
	}

	const size_t GetLength(const SequenceRecordField field) const {
		switch (field) {
		case SequenceRecordField::SEQUENCE1:
		case SequenceRecordField::QUAL1:
			return length_read1();
		case SequenceRecordField::SEQUENCE2:
		case SequenceRecordField::QUAL2:
			return length_read2();
		default:
			throw std::invalid_argument("Invalid field");
		}
	}
	const size_t length_read1() const {
		return read1.length();
	}

	const size_t length_read2() const {
		if (is_paired) {
			return read2.value().length();
		}
		throw UnpairedError("Attempted to access paired read with an unpaired record");
	}

private:
	inline std::string base_read_id(const std::string &id) const {
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

	inline void check_ids(const klibpp::KSeq &r1, const klibpp::KSeq &r2) const {
		std::string base_id1 = base_read_id(r1.name);
		std::string base_id2 = base_read_id(r2.name);

		if (base_id1 != base_id2) {
			throw std::runtime_error("Mismatched read IDs: " + r1.name + " vs " + r2.name);
		}
	}
}; // SequenceRecord
}; // namespace miint
