#include "ena_parser.hpp"
#include <expat.h>
#include <sstream>
#include <stdexcept>

namespace miint {

// ---- Accession type detection ----

ENAAccessionType ENAParser::DetectAccessionType(const std::string &accession) {
	if (accession.empty()) {
		return ENAAccessionType::UNKNOWN;
	}

	// Study: PRJNA*, PRJEB*, PRJDB*, ERP*, SRP*, DRP*
	if (accession.rfind("PRJNA", 0) == 0 || accession.rfind("PRJEB", 0) == 0 || accession.rfind("PRJDB", 0) == 0) {
		return ENAAccessionType::STUDY;
	}
	if (accession.rfind("ERP", 0) == 0 || accession.rfind("SRP", 0) == 0 || accession.rfind("DRP", 0) == 0) {
		return ENAAccessionType::STUDY;
	}

	// Sample: SAMN*, SAME*, SAMD*
	if (accession.rfind("SAMN", 0) == 0 || accession.rfind("SAME", 0) == 0 || accession.rfind("SAMD", 0) == 0) {
		return ENAAccessionType::SAMPLE;
	}

	// Run: SRR*, ERR*, DRR*
	if (accession.rfind("SRR", 0) == 0 || accession.rfind("ERR", 0) == 0 || accession.rfind("DRR", 0) == 0) {
		return ENAAccessionType::RUN;
	}

	// Experiment: SRX*, ERX*, DRX*
	if (accession.rfind("SRX", 0) == 0 || accession.rfind("ERX", 0) == 0 || accession.rfind("DRX", 0) == 0) {
		return ENAAccessionType::EXPERIMENT;
	}

	return ENAAccessionType::UNKNOWN;
}

std::string ENAParser::QueryParamForAccessionType(ENAAccessionType type) {
	switch (type) {
	case ENAAccessionType::STUDY:
		return "study_accession";
	case ENAAccessionType::SAMPLE:
		return "sample_accession";
	case ENAAccessionType::RUN:
		return "run_accession";
	case ENAAccessionType::EXPERIMENT:
		return "experiment_accession";
	default:
		return "accession";
	}
}

// ---- Validation ----

void ENAParser::ValidateAccession(const std::string &accession) {
	for (char c : accession) {
		if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' && c != '-' && c != '.') {
			throw std::invalid_argument("Invalid character in accession: '" + accession + "'");
		}
	}
}

void ENAParser::ValidateFields(const std::string &fields) {
	for (char c : fields) {
		if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' && c != ',') {
			throw std::invalid_argument("Invalid character in fields parameter: '" + fields + "'");
		}
	}
}

// ---- URL construction ----

std::string ENAParser::BuildSearchURL(const std::string &accession, const std::string &result_type,
                                      const std::string &fields) {
	ValidateAccession(accession);
	ValidateFields(fields);

	auto acc_type = DetectAccessionType(accession);
	auto query_param = QueryParamForAccessionType(acc_type);

	std::ostringstream url;
	url << PORTAL_BASE << "/search?"
	    << "result=" << result_type << "&query=" << query_param << "%3D%22" << accession << "%22"
	    << "&fields=" << fields << "&limit=0&format=tsv";
	return url.str();
}

std::string ENAParser::BuildSearchURLBatch(const std::vector<std::string> &accessions, ENAAccessionType accession_type,
                                           const std::string &result_type, const std::string &fields) {
	if (accessions.empty()) {
		throw std::invalid_argument("BuildSearchURLBatch: accessions vector must not be empty");
	}
	if (accession_type == ENAAccessionType::UNKNOWN) {
		throw std::invalid_argument("BuildSearchURLBatch: cannot build a batched query for UNKNOWN accession type");
	}
	for (const auto &acc : accessions) {
		ValidateAccession(acc);
	}
	ValidateFields(fields);

	auto query_param = QueryParamForAccessionType(accession_type);

	// URL-encoded form of: <query_param> IN ("acc1","acc2",...)
	//   space = %20   ( = %28   ) = %29   , = %2C   " = %22
	std::ostringstream url;
	url << PORTAL_BASE << "/search?"
	    << "result=" << result_type << "&query=" << query_param << "%20IN%20%28";
	for (size_t i = 0; i < accessions.size(); i++) {
		if (i > 0) {
			url << "%2C";
		}
		url << "%22" << accessions[i] << "%22";
	}
	url << "%29&fields=" << fields << "&limit=0&format=tsv";
	return url.str();
}

std::string ENAParser::BuildXMLURL(const std::vector<std::string> &accessions) {
	for (const auto &acc : accessions) {
		ValidateAccession(acc);
	}

	std::ostringstream url;
	url << BROWSER_BASE << "/xml/";
	for (size_t i = 0; i < accessions.size(); i++) {
		if (i > 0) {
			url << ",";
		}
		url << accessions[i];
	}
	return url.str();
}

// ---- TSV parsing ----

static std::vector<std::string> SplitTabs(const std::string &line) {
	std::vector<std::string> fields;
	std::string::size_type start = 0;
	while (true) {
		auto pos = line.find('\t', start);
		if (pos == std::string::npos) {
			fields.push_back(line.substr(start));
			break;
		}
		fields.push_back(line.substr(start, pos - start));
		start = pos + 1;
	}
	return fields;
}

ENATSVResult ENAParser::ParseTSV(const std::string &tsv) {
	ENATSVResult result;
	if (tsv.empty()) {
		return result;
	}

	std::istringstream stream(tsv);
	std::string line;

	// First line is header
	if (!std::getline(stream, line)) {
		return result;
	}
	if (!line.empty() && line.back() == '\r') {
		line.pop_back();
	}
	result.column_names = SplitTabs(line);

	// Remaining lines are data
	while (std::getline(stream, line)) {
		if (!line.empty() && line.back() == '\r') {
			line.pop_back();
		}
		if (line.empty()) {
			continue;
		}
		result.rows.push_back(SplitTabs(line));
	}

	return result;
}

// ---- XML parsing for sample attributes ----

struct XMLParserState {
	std::string current_sample_accession;
	std::string current_tag;
	std::string current_value;
	std::string char_buffer;
	bool in_tag = false;
	bool in_value = false;
	std::vector<SampleAttribute> attributes;
};

static void XMLCALL xml_start_element(void *user_data, const char *name, const char **attrs) {
	auto *state = static_cast<XMLParserState *>(user_data);
	std::string elem(name);

	if (elem == "SAMPLE") {
		for (int i = 0; attrs[i]; i += 2) {
			if (std::string(attrs[i]) == "accession") {
				state->current_sample_accession = attrs[i + 1];
				break;
			}
		}
	} else if (elem == "SAMPLE_ATTRIBUTE") {
		state->current_tag.clear();
		state->current_value.clear();
	} else if (elem == "TAG") {
		state->in_tag = true;
		state->char_buffer.clear();
	} else if (elem == "VALUE") {
		state->in_value = true;
		state->char_buffer.clear();
	}
}

static void XMLCALL xml_end_element(void *user_data, const char *name) {
	auto *state = static_cast<XMLParserState *>(user_data);
	std::string elem(name);

	if (elem == "TAG") {
		state->current_tag = state->char_buffer;
		state->in_tag = false;
	} else if (elem == "VALUE") {
		state->current_value = state->char_buffer;
		state->in_value = false;
	} else if (elem == "SAMPLE_ATTRIBUTE") {
		if (!state->current_tag.empty()) {
			state->attributes.push_back({state->current_sample_accession, state->current_tag, state->current_value});
		}
	}
}

static void XMLCALL xml_char_data(void *user_data, const char *s, int len) {
	auto *state = static_cast<XMLParserState *>(user_data);
	if (state->in_tag || state->in_value) {
		state->char_buffer.append(s, static_cast<size_t>(len));
	}
}

std::vector<SampleAttribute> ENAParser::ParseSampleAttributesXML(const std::string &xml) {
	XMLParserState state;

	XML_Parser parser = XML_ParserCreate(nullptr);
	if (!parser) {
		throw std::runtime_error("Failed to create XML parser");
	}

	// RAII guard ensures XML_ParserFree on all exit paths (including exceptions from callbacks)
	struct ParserGuard {
		XML_Parser p;
		~ParserGuard() {
			XML_ParserFree(p);
		}
	} guard {parser};

	XML_SetUserData(parser, &state);
	XML_SetElementHandler(parser, xml_start_element, xml_end_element);
	XML_SetCharacterDataHandler(parser, xml_char_data);

	if (XML_Parse(parser, xml.data(), static_cast<int>(xml.size()), XML_TRUE) == XML_STATUS_ERROR) {
		auto error_msg = std::string(XML_ErrorString(XML_GetErrorCode(parser)));
		throw std::runtime_error("XML parse error: " + error_msg);
	}

	return state.attributes;
}

// ---- Filter submitted files by format ----

// Split a semicolon-delimited string into a vector of strings.
static std::vector<std::string> SplitSemicolon(const std::string &field) {
	std::vector<std::string> parts;
	if (field.empty()) {
		return parts;
	}
	std::string::size_type start = 0;
	while (true) {
		auto pos = field.find(';', start);
		if (pos == std::string::npos) {
			parts.push_back(field.substr(start));
			break;
		}
		parts.push_back(field.substr(start, pos - start));
		start = pos + 1;
	}
	return parts;
}

ENAParser::SubmittedFilterResult ENAParser::FilterSubmittedByFormat(const std::string &submitted_ftp,
                                                                    const std::string &submitted_aspera,
                                                                    const std::string &submitted_format,
                                                                    const std::string &submitted_bytes,
                                                                    const std::string &target_format) {
	SubmittedFilterResult result;

	auto ftp_parts = SplitSemicolon(submitted_ftp);
	auto aspera_parts = SplitSemicolon(submitted_aspera);
	auto format_parts = SplitSemicolon(submitted_format);
	auto bytes_parts = SplitSemicolon(submitted_bytes);

	// Use minimum count across ftp and format (both are required)
	size_t count = std::min(ftp_parts.size(), format_parts.size());

	for (size_t i = 0; i < count; i++) {
		if (format_parts[i] != target_format) {
			continue;
		}

		// Convert FTP path to HTTPS URL
		auto urls = FTPtoHTTPS(ftp_parts[i]);
		if (!urls.empty()) {
			result.urls.push_back(std::move(urls[0]));
		}

		// Preserve raw aspera entry if available at this index
		if (i < aspera_parts.size() && !aspera_parts[i].empty()) {
			result.aspera_raw.push_back(aspera_parts[i]);
		}

		// Accumulate bytes
		if (i < bytes_parts.size()) {
			try {
				result.total_bytes += std::stoull(bytes_parts[i]);
			} catch (...) {
				// Invalid number — treat as 0
			}
		}
	}

	return result;
}

// ---- FTP to HTTPS ----

std::vector<std::string> ENAParser::FTPtoHTTPS(const std::string &ftp_field) {
	std::vector<std::string> urls;
	if (ftp_field.empty()) {
		return urls;
	}

	std::string::size_type start = 0;
	while (true) {
		auto pos = ftp_field.find(';', start);
		std::string path;
		if (pos == std::string::npos) {
			path = ftp_field.substr(start);
		} else {
			path = ftp_field.substr(start, pos - start);
		}
		if (!path.empty()) {
			// Strip ftp:// prefix if present before prepending https://
			if (path.rfind("ftp://", 0) == 0) {
				path = path.substr(6);
			}
			urls.push_back("https://" + path);
		}
		if (pos == std::string::npos) {
			break;
		}
		start = pos + 1;
	}

	return urls;
}

// ---- Default fields ----

std::string ENAParser::DefaultFields(const std::string &result_type) {
	if (result_type == "read_run") {
		return "run_accession,experiment_accession,sample_accession,study_accession,"
		       "fastq_ftp,fastq_aspera,fastq_bytes,fastq_md5,"
		       "library_strategy,library_source,library_selection,library_layout,library_name,"
		       "instrument_model,instrument_platform,"
		       "read_count,base_count,"
		       "sample_alias,sample_title,scientific_name,tax_id,"
		       "collection_date,country,lat,lon,depth";
	}
	if (result_type == "sample") {
		return "sample_accession,secondary_sample_accession,sample_alias,sample_title,"
		       "description,scientific_name,tax_id,tax_lineage,"
		       "collection_date,country,lat,lon,depth,"
		       "environment_biome,environment_feature,environment_material,"
		       "host,isolation_source";
	}
	if (result_type == "study") {
		return "study_accession,secondary_study_accession,study_title,study_description,"
		       "center_name,first_public,last_updated,scientific_name,tax_id";
	}
	return "";
}

} // namespace miint
