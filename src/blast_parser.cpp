#include "blast_parser.hpp"
#include <sstream>
#include <stdexcept>

namespace miint {

BlastSubmitResult BlastParser::ParseSubmitResponse(const std::string &html) {
	BlastSubmitResult result;

	auto extract_qblast_value = [&](const std::string &key) -> std::string {
		auto pos = html.find(key);
		if (pos == std::string::npos) {
			return "";
		}
		pos += key.size();
		while (pos < html.size() && (html[pos] == ' ' || html[pos] == '=')) {
			pos++;
		}
		auto end = pos;
		while (end < html.size() && html[end] != '\n' && html[end] != '\r' && html[end] != ' ') {
			end++;
		}
		return html.substr(pos, end - pos);
	};

	result.rid = extract_qblast_value("RID ");
	std::string rtoe_str = extract_qblast_value("RTOE ");
	if (!rtoe_str.empty()) {
		try {
			result.rtoe = std::stoi(rtoe_str);
		} catch (...) {
			result.rtoe = 0;
		}
	}

	auto error_marker = html.find("class=\"error\">");
	if (error_marker != std::string::npos) {
		auto tag_close = html.find('>', error_marker);
		if (tag_close != std::string::npos) {
			auto content_start = tag_close + 1;
			auto content_end = html.find('<', content_start);
			if (content_end != std::string::npos) {
				result.error = html.substr(content_start, content_end - content_start);
			}
		}
	}

	return result;
}

BlastStatus BlastParser::ParseStatusResponse(const std::string &text) {
	if (text.find("Status=WAITING") != std::string::npos) {
		return BlastStatus::WAITING;
	}
	if (text.find("Status=READY") != std::string::npos) {
		return BlastStatus::READY;
	}
	return BlastStatus::UNKNOWN;
}

std::vector<BlastHit> BlastParser::ParseTabularResults(const std::string &tabular) {
	std::vector<BlastHit> hits;
	if (tabular.empty()) {
		return hits;
	}

	std::istringstream stream(tabular);
	std::string line;
	while (std::getline(stream, line)) {
		if (!line.empty() && line.back() == '\r') {
			line.pop_back();
		}
		if (line.empty() || line[0] == '#') {
			continue;
		}

		std::vector<std::string> fields;
		std::istringstream line_stream(line);
		std::string field;
		while (std::getline(line_stream, field, '\t')) {
			fields.push_back(field);
		}

		if (fields.size() < 12) {
			continue;
		}

		try {
			BlastHit hit;
			hit.query_id = fields[0];
			hit.subject_id = fields[1];
			hit.pct_identity = std::stod(fields[2]);
			hit.alignment_length = std::stoi(fields[3]);
			hit.mismatches = std::stoi(fields[4]);
			hit.gap_opens = std::stoi(fields[5]);
			hit.query_start = std::stoll(fields[6]);
			hit.query_end = std::stoll(fields[7]);
			hit.subject_start = std::stoll(fields[8]);
			hit.subject_end = std::stoll(fields[9]);
			hit.evalue = std::stod(fields[10]);
			hit.bit_score = std::stod(fields[11]);
			hits.push_back(std::move(hit));
		} catch (...) {
			continue;
		}
	}

	return hits;
}

bool BlastParser::ValidateProgram(const std::string &program) {
	return program == "blastn" || program == "blastp" || program == "blastx" || program == "tblastn" ||
	       program == "tblastx";
}

std::string BlastParser::BuildFastaPayload(const std::vector<std::string> &ids,
                                           const std::vector<std::string> &sequences) {
	if (ids.size() != sequences.size()) {
		throw std::invalid_argument("ids and sequences must have the same size");
	}
	std::string result;
	size_t total = 0;
	for (size_t i = 0; i < ids.size(); i++) {
		total += 1 + ids[i].size() + 1 + sequences[i].size() + 1;
	}
	result.reserve(total);
	for (size_t i = 0; i < ids.size(); i++) {
		result += '>';
		result += ids[i];
		result += '\n';
		result += sequences[i];
		result += '\n';
	}
	return result;
}

std::string BlastParser::UrlEncode(const std::string &value) {
	std::string encoded;
	encoded.reserve(value.size() + value.size() / 4);
	for (unsigned char c : value) {
		if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '-' || c == '_' ||
		    c == '.' || c == '~') {
			encoded += static_cast<char>(c);
		} else {
			char hex[4];
			snprintf(hex, sizeof(hex), "%%%02X", c);
			encoded += hex;
		}
	}
	return encoded;
}

std::string BlastParser::BuildSubmitBody(const std::string &program, const std::string &database,
                                         const std::string &fasta_query, double evalue, int max_targets,
                                         bool megablast) {
	std::ostringstream body;
	body << "CMD=Put";
	body << "&PROGRAM=" << UrlEncode(program);
	body << "&DATABASE=" << UrlEncode(database);

	std::ostringstream evalue_ss;
	evalue_ss << evalue;
	body << "&EXPECT=" << evalue_ss.str();

	body << "&HITLIST_SIZE=" << max_targets;
	if (megablast) {
		body << "&MEGABLAST=on";
	}
	body << "&QUERY=" << UrlEncode(fasta_query);
	return body.str();
}

std::vector<BlastParser::SequenceBatch> BlastParser::SplitIntoBatches(const std::vector<std::string> &ids,
                                                                      const std::vector<std::string> &sequences,
                                                                      size_t max_bytes) {
	std::vector<SequenceBatch> batches;
	if (ids.empty()) {
		return batches;
	}

	SequenceBatch current;
	size_t current_size = 0;

	for (size_t i = 0; i < ids.size(); i++) {
		// FASTA record: ">id\nsequence\n"
		size_t record_size = 1 + ids[i].size() + 1 + sequences[i].size() + 1;

		if (!current.ids.empty() && current_size + record_size > max_bytes) {
			batches.push_back(std::move(current));
			current = SequenceBatch {};
			current_size = 0;
		}

		current.ids.push_back(ids[i]);
		current.sequences.push_back(sequences[i]);
		current_size += record_size;
	}

	if (!current.ids.empty()) {
		batches.push_back(std::move(current));
	}

	return batches;
}

} // namespace miint
