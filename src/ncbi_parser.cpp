#include "ncbi_parser.hpp"
#include <algorithm>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>

namespace miint {

// Default warning callback writes to stderr
void DefaultWarningCallback(const std::string &msg) {
	std::cerr << "[NCBI Warning] " << msg << std::endl;
}

// Helper to check if string starts with prefix (safe, no substr)
static bool StartsWith(const std::string &str, const std::string &prefix) {
	return str.length() >= prefix.length() && str.compare(0, prefix.length(), prefix) == 0;
}

AccessionType NCBIParser::DetectAccessionType(const std::string &accession) {
	if (accession.empty()) {
		return AccessionType::UNKNOWN;
	}

	// Assembly accessions: GCF_ or GCA_ prefix (safe - uses StartsWith)
	if (StartsWith(accession, "GCF_") || StartsWith(accession, "GCA_")) {
		return AccessionType::ASSEMBLY;
	}

	// Sequence accessions: Various RefSeq/GenBank prefixes
	static const std::vector<std::string> sequence_prefixes = {"NC_", "NM_", "NP_", "NR_", "XM_", "XP_",
	                                                           "XR_", "NG_", "NT_", "NW_", "NZ_", "AC_"};

	for (const auto &prefix : sequence_prefixes) {
		if (StartsWith(accession, prefix)) {
			return AccessionType::SEQUENCE;
		}
	}

	return AccessionType::UNKNOWN;
}

bool NCBIParser::IsAssemblyAccession(const std::string &accession) {
	return DetectAccessionType(accession) == AccessionType::ASSEMBLY;
}

// Extract accession from various FASTA ID formats
// Handles: "NC_001416.1", "ref|NC_001416.1|", "gi|123|ref|NC_001416.1|description"
std::string NCBIParser::ExtractAccessionFromFastaId(const std::string &fasta_id, std::string &remainder) {
	remainder.clear();

	if (fasta_id.empty()) {
		return "";
	}

	// Check for pipe-delimited format (e.g., "gi|123|ref|NC_001416.1|description")
	if (fasta_id.find('|') != std::string::npos) {
		std::vector<std::string> parts;
		std::istringstream ss(fasta_id);
		std::string part;
		while (std::getline(ss, part, '|')) {
			parts.push_back(part);
		}

		// Look for ref|accession| pattern
		for (size_t i = 0; i + 1 < parts.size(); i++) {
			if (parts[i] == "ref" || parts[i] == "gb" || parts[i] == "emb" || parts[i] == "dbj") {
				// Next part is the accession
				std::string accession = parts[i + 1];
				// Everything after is remainder (join remaining parts)
				if (i + 2 < parts.size()) {
					for (size_t j = i + 2; j < parts.size(); j++) {
						if (!remainder.empty()) {
							remainder += " ";
						}
						remainder += parts[j];
					}
				}
				return accession;
			}
		}

		// Fallback: use first non-empty part as accession
		for (const auto &p : parts) {
			if (!p.empty() && p != "gi" && p != "ref" && p != "gb") {
				// Check if it looks like an accession (has underscore or dot)
				if (p.find('_') != std::string::npos || p.find('.') != std::string::npos) {
					return p;
				}
			}
		}
	}

	// Simple format: just split on first space
	size_t space_pos = fasta_id.find(' ');
	if (space_pos != std::string::npos) {
		remainder = fasta_id.substr(space_pos + 1);
		return fasta_id.substr(0, space_pos);
	}

	return fasta_id;
}

SequenceRecordBatch NCBIParser::ParseFasta(const std::string &fasta_text) {
	SequenceRecordBatch batch(false); // Not paired

	if (fasta_text.empty()) {
		return batch;
	}

	std::istringstream stream(fasta_text);
	std::string line;
	std::string current_id;
	std::string current_comment;
	std::string current_sequence;

	auto flush_record = [&]() {
		if (!current_id.empty()) {
			batch.read_ids.push_back(current_id);
			batch.comments.push_back(current_comment);
			batch.sequences1.push_back(current_sequence);
			batch.quals1.push_back(QualScore("")); // Empty quality for FASTA
		}
		current_id.clear();
		current_comment.clear();
		current_sequence.clear();
	};

	while (std::getline(stream, line)) {
		// Remove trailing whitespace/CR
		while (!line.empty() && (line.back() == '\r' || line.back() == ' ' || line.back() == '\t')) {
			line.pop_back();
		}

		if (line.empty()) {
			continue;
		}

		if (!line.empty() && line[0] == '>') {
			// New record - flush previous
			flush_record();

			// Parse header: >id comment or >gi|123|ref|accession|description
			std::string header = line.substr(1);

			// Use smart extraction that handles pipe-delimited format
			current_id = ExtractAccessionFromFastaId(header, current_comment);
		} else {
			// Sequence line
			current_sequence += line;
		}
	}

	// Flush last record
	flush_record();

	return batch;
}

// Helper to extract text between XML tags (improved robustness)
// Finds the first occurrence of <tag>...</tag> at the top level
static std::string ExtractXMLValue(const std::string &xml, const std::string &tag) {
	std::string open_tag = "<" + tag + ">";
	std::string close_tag = "</" + tag + ">";

	// Find first opening tag
	size_t start = xml.find(open_tag);
	if (start == std::string::npos) {
		// Try with attributes: <tag attr="value">
		std::string open_tag_prefix = "<" + tag + " ";
		start = xml.find(open_tag_prefix);
		if (start == std::string::npos) {
			return "";
		}
		// Find the end of the opening tag
		size_t tag_end = xml.find('>', start);
		if (tag_end == std::string::npos) {
			return "";
		}
		start = tag_end + 1;
	} else {
		start += open_tag.length();
	}

	// Find matching closing tag (handle nested tags by counting)
	int depth = 1;
	size_t pos = start;
	while (depth > 0 && pos < xml.length()) {
		size_t next_open = xml.find("<" + tag, pos);
		size_t next_close = xml.find(close_tag, pos);

		if (next_close == std::string::npos) {
			return ""; // Malformed XML
		}

		if (next_open != std::string::npos && next_open < next_close) {
			// Nested opening tag
			depth++;
			pos = next_open + 1;
		} else {
			// Closing tag
			depth--;
			if (depth == 0) {
				return xml.substr(start, next_close - start);
			}
			pos = next_close + close_tag.length();
		}
	}

	return "";
}

// Extract taxonomy ID from db_xref qualifier
static int64_t ExtractTaxonomyId(const std::string &xml) {
	std::regex taxon_regex(R"(taxon:(\d+))");
	std::smatch match;
	if (std::regex_search(xml, match, taxon_regex)) {
		try {
			return std::stoll(match[1].str());
		} catch (...) {
			return 0;
		}
	}
	return 0;
}

// Parse version from accession.version format
static int32_t ExtractVersion(const std::string &accession_version) {
	size_t dot = accession_version.rfind('.');
	if (dot != std::string::npos && dot + 1 < accession_version.length()) {
		try {
			return std::stoi(accession_version.substr(dot + 1));
		} catch (...) {
			return 0;
		}
	}
	return 0;
}

// Validate date components
static bool ValidateDate(int day, int month, int year) {
	if (month < 1 || month > 12)
		return false;
	if (day < 1 || day > 31)
		return false;
	if (year < 1900 || year > 2100)
		return false;

	// Check days per month
	static const int days_in_month[] = {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
	if (day > days_in_month[month])
		return false;

	return true;
}

GenBankMetadata NCBIParser::ParseGenBankXML(const std::string &xml) {
	GenBankMetadata metadata;

	if (xml.empty()) {
		return metadata;
	}

	// Extract fields from GenBank XML format
	metadata.accession = ExtractXMLValue(xml, "GBSeq_accession-version");
	if (metadata.accession.empty()) {
		metadata.accession = ExtractXMLValue(xml, "GBSeq_primary-accession");
	}

	metadata.version = ExtractVersion(metadata.accession);
	metadata.description = ExtractXMLValue(xml, "GBSeq_definition");
	metadata.organism = ExtractXMLValue(xml, "GBSeq_organism");
	metadata.molecule_type = ExtractXMLValue(xml, "GBSeq_moltype");

	// Parse length
	std::string length_str = ExtractXMLValue(xml, "GBSeq_length");
	if (!length_str.empty()) {
		try {
			metadata.length = std::stoll(length_str);
		} catch (...) {
			metadata.length = 0;
		}
	}

	// Extract taxonomy ID from feature table
	metadata.taxonomy_id = ExtractTaxonomyId(xml);

	// Parse update date (format: DD-MON-YYYY -> YYYY-MM-DD)
	std::string update_date = ExtractXMLValue(xml, "GBSeq_update-date");
	if (!update_date.empty()) {
		static const std::map<std::string, int> months = {{"JAN", 1}, {"FEB", 2},  {"MAR", 3},  {"APR", 4},
		                                                  {"MAY", 5}, {"JUN", 6},  {"JUL", 7},  {"AUG", 8},
		                                                  {"SEP", 9}, {"OCT", 10}, {"NOV", 11}, {"DEC", 12}};

		std::regex date_regex(R"((\d{1,2})-([A-Za-z]{3})-(\d{4}))");
		std::smatch match;
		if (std::regex_match(update_date, match, date_regex)) {
			int day = std::stoi(match[1].str());
			std::string month_str = match[2].str();
			// Convert to uppercase for lookup
			std::transform(month_str.begin(), month_str.end(), month_str.begin(), ::toupper);
			int year = std::stoi(match[3].str());

			auto it = months.find(month_str);
			if (it != months.end()) {
				int month = it->second;
				if (ValidateDate(day, month, year)) {
					char buf[11];
					snprintf(buf, sizeof(buf), "%04d-%02d-%02d", year, month, day);
					metadata.update_date = buf;
				}
			}
		}
	}

	return metadata;
}

// Detect source type from feature table header or accession
static std::string DetectSource(const std::string &seqid) {
	if (seqid.empty()) {
		return "unknown";
	}

	// RefSeq accessions start with specific prefixes
	if (StartsWith(seqid, "NC_") || StartsWith(seqid, "NM_") || StartsWith(seqid, "NP_") || StartsWith(seqid, "NR_") ||
	    StartsWith(seqid, "XM_") || StartsWith(seqid, "XP_") || StartsWith(seqid, "XR_") || StartsWith(seqid, "NG_") ||
	    StartsWith(seqid, "NT_") || StartsWith(seqid, "NW_") || StartsWith(seqid, "NZ_")) {
		return "RefSeq";
	}

	// GenBank/INSDC accessions (typically 1-2 letters followed by digits)
	if (seqid.length() >= 5) {
		bool is_genbank = true;
		size_t letter_count = 0;
		for (size_t i = 0; i < seqid.length() && seqid[i] != '.'; i++) {
			if (std::isalpha(seqid[i])) {
				letter_count++;
				if (letter_count > 2)
					is_genbank = false;
			} else if (std::isdigit(seqid[i])) {
				break;
			} else {
				is_genbank = false;
				break;
			}
		}
		if (is_genbank && letter_count > 0) {
			return "GenBank";
		}
	}

	return "NCBI";
}

// Check if a position string contains complex location (join, complement, order)
static bool IsComplexLocation(const std::string &pos_str) {
	return pos_str.find("join") != std::string::npos || pos_str.find("order") != std::string::npos ||
	       pos_str.find("complement") != std::string::npos || pos_str.find("..") != std::string::npos ||
	       pos_str.find(',') != std::string::npos;
}

// Parse a simple position (handling partial indicators < and >)
static int64_t ParseSimplePosition(const std::string &pos_str) {
	std::string clean = pos_str;
	// Remove partial indicators
	if (!clean.empty() && (clean[0] == '<' || clean[0] == '>')) {
		clean = clean.substr(1);
	}
	// Remove any non-digit suffix
	size_t end = clean.find_first_not_of("0123456789");
	if (end != std::string::npos) {
		clean = clean.substr(0, end);
	}
	if (clean.empty()) {
		return 0;
	}
	try {
		return std::stoll(clean);
	} catch (...) {
		return 0;
	}
}

FeatureAnnotationBatch NCBIParser::ParseFeatureTable(const std::string &feature_table, WarningCallback warn_callback) {
	FeatureAnnotationBatch batch;

	if (feature_table.empty()) {
		return batch;
	}

	std::istringstream stream(feature_table);
	std::string line;
	std::string current_seqid;
	std::string current_source;
	FeatureAnnotation *current_feature = nullptr;
	bool warned_complex_location = false;

	// Lambda for emitting warnings
	auto warn = [&](const std::string &msg) {
		if (warn_callback) {
			warn_callback(msg);
		}
	};

	while (std::getline(stream, line)) {
		// Remove trailing whitespace/CR
		while (!line.empty() && (line.back() == '\r' || line.back() == ' ')) {
			line.pop_back();
		}

		if (line.empty()) {
			continue;
		}

		// Header line: >Feature ref|NC_001416.1|
		if (!line.empty() && line[0] == '>') {
			static const std::string feature_prefix = ">Feature";
			if (StartsWith(line, feature_prefix)) {
				// Extract sequence ID from header
				std::string header =
				    line.length() > feature_prefix.length() + 1 ? line.substr(feature_prefix.length() + 1) : "";

				// Remove leading whitespace
				size_t start = header.find_first_not_of(" \t");
				if (start != std::string::npos) {
					header = header.substr(start);
				}

				// Handle ref|NC_001416.1| format (safe with StartsWith)
				if (StartsWith(header, "ref|")) {
					header = header.substr(4);
					size_t end = header.find('|');
					if (end != std::string::npos) {
						header = header.substr(0, end);
					}
				} else if (StartsWith(header, "gb|") || StartsWith(header, "emb|") || StartsWith(header, "dbj|")) {
					header = header.substr(3);
					size_t end = header.find('|');
					if (end != std::string::npos) {
						header = header.substr(0, end);
					}
				}

				current_seqid = header;
				current_source = DetectSource(current_seqid);
			}
			continue;
		}

		// Check if this is a qualifier line (starts with tabs) or feature line
		if (!line.empty() && line[0] == '\t') {
			// Qualifier line - add to current feature
			if (current_feature != nullptr) {
				// Count leading tabs
				size_t tab_count = 0;
				while (tab_count < line.size() && line[tab_count] == '\t') {
					tab_count++;
				}

				// Validate tab count (should be 3 for standard INSDC format)
				if (tab_count != 3) {
					warn("Unexpected indentation (" + std::to_string(tab_count) +
					     " tabs) in feature table, expected 3");
				}

				std::string qualifier_line = line.substr(tab_count);

				// Split on tab: key\tvalue
				size_t tab_pos = qualifier_line.find('\t');
				std::string key, value;
				if (tab_pos != std::string::npos) {
					key = qualifier_line.substr(0, tab_pos);
					value = qualifier_line.substr(tab_pos + 1);
				} else {
					key = qualifier_line;
					value = "";
				}

				current_feature->attrs.emplace_back(key, value);

				// Handle codon_start for CDS phase calculation (#18)
				if (current_feature->type == "CDS" && key == "codon_start") {
					try {
						int codon_start = std::stoi(value);
						// codon_start 1 -> phase 0, codon_start 2 -> phase 2, codon_start 3 -> phase 1
						if (codon_start == 1) {
							current_feature->phase = 0;
						} else if (codon_start == 2) {
							current_feature->phase = 2;
						} else if (codon_start == 3) {
							current_feature->phase = 1;
						}
					} catch (...) {
						// Keep default phase
					}
				}
			}
		} else {
			// Feature line: start\tstop\ttype
			std::istringstream feature_stream(line);
			std::string start_str, stop_str, type;

			if (!(feature_stream >> start_str >> stop_str >> type)) {
				continue; // Invalid line
			}

			// Check for complex locations (#17)
			if (!warned_complex_location && (IsComplexLocation(start_str) || IsComplexLocation(stop_str))) {
				warn("Complex feature location detected (join/complement/order). "
				     "Using outer bounds only. For full location support, use GenBank flat file format.");
				warned_complex_location = true;
			}

			// Create new feature
			batch.features.emplace_back();
			current_feature = &batch.features.back();

			current_feature->seqid = current_seqid;
			current_feature->source = current_source;
			current_feature->type = type;
			current_feature->has_score = false;

			// Parse positions using int64_t (#5)
			int64_t pos1 = ParseSimplePosition(start_str);
			int64_t pos2 = ParseSimplePosition(stop_str);

			// Determine strand based on position order
			if (pos1 > pos2) {
				// Complement strand (positions are reversed in feature table)
				current_feature->position = pos2;
				current_feature->stop_position = pos1;
				current_feature->strand = "-";
			} else {
				current_feature->position = pos1;
				current_feature->stop_position = pos2;
				current_feature->strand = "+";
			}

			// Default phase - will be overridden by codon_start qualifier for CDS
			if (type == "CDS") {
				current_feature->phase = 0;
			} else {
				current_feature->phase = -1;
			}
		}
	}

	return batch;
}

} // namespace miint
