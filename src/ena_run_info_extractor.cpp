#include "ena_run_info_extractor.hpp"

namespace miint {

ENAColumnIndexMap ENAColumnIndexMap::FromHeader(const std::vector<std::string> &header) {
	ENAColumnIndexMap cols;
	for (size_t i = 0; i < header.size(); i++) {
		const auto &name = header[i];
		auto idx = static_cast<int>(i);
		if (name == "run_accession") {
			cols.run_accession = idx;
		} else if (name == "sample_accession") {
			cols.sample_accession = idx;
		} else if (name == "experiment_accession") {
			cols.experiment_accession = idx;
		} else if (name == "study_accession") {
			cols.study_accession = idx;
		} else if (name == "fastq_ftp") {
			cols.fastq_ftp = idx;
		} else if (name == "fastq_aspera") {
			cols.fastq_aspera = idx;
		} else if (name == "fastq_bytes") {
			cols.fastq_bytes = idx;
		} else if (name == "fastq_md5") {
			cols.fastq_md5 = idx;
		} else if (name == "library_layout") {
			cols.library_layout = idx;
		} else if (name == "submitted_ftp") {
			cols.submitted_ftp = idx;
		} else if (name == "submitted_aspera") {
			cols.submitted_aspera = idx;
		} else if (name == "submitted_bytes") {
			cols.submitted_bytes = idx;
		} else if (name == "submitted_format") {
			cols.submitted_format = idx;
		}
	}
	return cols;
}

std::string ENAColumnIndexMap::Get(const std::vector<std::string> &row, int col) {
	if (col < 0 || col >= static_cast<int>(row.size())) {
		return "";
	}
	return row[col];
}

static std::vector<uint64_t> ParseBytesField(const std::string &bytes_field) {
	std::vector<uint64_t> out;
	if (bytes_field.empty()) {
		return out;
	}
	std::string::size_type start = 0;
	while (true) {
		auto semi = bytes_field.find(';', start);
		std::string token =
		    (semi == std::string::npos) ? bytes_field.substr(start) : bytes_field.substr(start, semi - start);
		try {
			out.push_back(std::stoull(token));
		} catch (...) {
			out.push_back(0);
		}
		if (semi == std::string::npos) {
			break;
		}
		start = semi + 1;
	}
	return out;
}

// Splits a `;`-separated field into its raw string tokens. Used for fastq_md5,
// which -- unlike fastq_bytes -- has no numeric interpretation and no
// "malformed token" case worth defaulting (a malformed/empty token is simply
// kept as-is; ENARunInfoExtractor only cares whether the resulting count
// lines up with fastq_urls, not the content of any one token).
static std::vector<std::string> SplitSemicolonField(const std::string &field) {
	std::vector<std::string> out;
	if (field.empty()) {
		return out;
	}
	std::string::size_type start = 0;
	while (true) {
		auto semi = field.find(';', start);
		out.push_back((semi == std::string::npos) ? field.substr(start) : field.substr(start, semi - start));
		if (semi == std::string::npos) {
			break;
		}
		start = semi + 1;
	}
	return out;
}

std::vector<ENARunInfo> ENARunInfoExtractor::FromTSVRow(const std::vector<std::string> &row,
                                                        const ENAColumnIndexMap &cols,
                                                        const std::string &prefer_format) {
	std::vector<ENARunInfo> out;

	const std::string run_accession = ENAColumnIndexMap::Get(row, cols.run_accession);
	const std::string sample_accession = ENAColumnIndexMap::Get(row, cols.sample_accession);
	const std::string experiment_accession = ENAColumnIndexMap::Get(row, cols.experiment_accession);
	const std::string ftp_field = ENAColumnIndexMap::Get(row, cols.fastq_ftp);
	const std::string aspera_field = ENAColumnIndexMap::Get(row, cols.fastq_aspera);
	const std::string bytes_field = ENAColumnIndexMap::Get(row, cols.fastq_bytes);
	const std::string md5_field = ENAColumnIndexMap::Get(row, cols.fastq_md5);
	const std::string layout = ENAColumnIndexMap::Get(row, cols.library_layout);
	const std::string sub_ftp = ENAColumnIndexMap::Get(row, cols.submitted_ftp);
	const std::string sub_aspera = ENAColumnIndexMap::Get(row, cols.submitted_aspera);
	const std::string sub_bytes = ENAColumnIndexMap::Get(row, cols.submitted_bytes);
	const std::string sub_format = ENAColumnIndexMap::Get(row, cols.submitted_format);

	auto fastq_urls = ENAParser::FTPtoHTTPS(ftp_field);
	const bool has_fastq = !fastq_urls.empty();

	auto sff_filter = ENAParser::FilterSubmittedByFormat(sub_ftp, sub_aspera, sub_format, sub_bytes, "SFF");
	const bool has_sff = !sff_filter.urls.empty();

	bool use_sff = false;
	if (prefer_format == "sff") {
		use_sff = has_sff;
	} else if (prefer_format == "auto") {
		use_sff = !has_fastq && has_sff;
	}
	// prefer_format == "fastq" → use_sff stays false, and if no FASTQ we skip the row entirely

	if (use_sff) {
		// Flatten: one RunInfo per SFF file. use_sff is only set when has_sff is true,
		// so sff_filter.urls is guaranteed non-empty here.
		const uint64_t per_file_bytes = sff_filter.total_bytes / sff_filter.urls.size();
		for (const auto &url : sff_filter.urls) {
			ENARunInfo info;
			info.run_accession = run_accession;
			info.sample_accession = sample_accession;
			info.experiment_accession = experiment_accession;
			info.format = ENASequenceFormat::SFF;
			info.sff_url = url;
			info.is_paired = false;
			info.total_bytes = per_file_bytes;
			out.push_back(std::move(info));
		}
		return out;
	}

	if (!has_fastq) {
		return out; // prefer_format='fastq' with no FASTQ, or no formats at all
	}

	ENARunInfo info;
	info.run_accession = run_accession;
	info.sample_accession = sample_accession;
	info.experiment_accession = experiment_accession;
	info.fastq_urls = std::move(fastq_urls);
	info.aspera_paths = AsperaUtils::ParseAsperaPaths(aspera_field);
	info.is_paired = (layout == "PAIRED");

	auto file_bytes = ParseBytesField(bytes_field);
	auto file_md5 = SplitSemicolonField(md5_field);

	// ENA 3-file paired-end filtering: drop the non-_1/_2 file when present
	if (info.is_paired && info.fastq_urls.size() > 2) {
		std::vector<std::string> filtered;
		std::vector<uint64_t> filtered_bytes;
		std::vector<std::string> filtered_md5;
		for (size_t fi = 0; fi < info.fastq_urls.size(); fi++) {
			const auto &url = info.fastq_urls[fi];
			if (url.find("_1.fast") != std::string::npos || url.find("_2.fast") != std::string::npos) {
				filtered.push_back(url);
				if (fi < file_bytes.size()) {
					filtered_bytes.push_back(file_bytes[fi]);
				}
				if (fi < file_md5.size()) {
					filtered_md5.push_back(file_md5[fi]);
				}
			}
		}
		if (filtered.size() == 2) {
			info.fastq_urls = std::move(filtered);
			file_bytes = std::move(filtered_bytes);
			file_md5 = std::move(filtered_md5);
		}
	}
	if (info.is_paired && info.aspera_paths.size() > 2) {
		std::vector<AsperaPath> filtered;
		for (const auto &ap : info.aspera_paths) {
			if (ap.remote_path.find("_1.fast") != std::string::npos ||
			    ap.remote_path.find("_2.fast") != std::string::npos) {
				filtered.push_back(ap);
			}
		}
		if (filtered.size() == 2) {
			info.aspera_paths = std::move(filtered);
		}
	}

	info.total_bytes = 0;
	for (auto b : file_bytes) {
		info.total_bytes += b;
	}

	// fastq_md5 must align 1:1 with the (possibly filtered) fastq_urls. ENA
	// occasionally omits the field or its count doesn't match (rare metadata
	// inconsistency); a mismatch degrades to "no md5 for this run" rather than
	// throwing or misaligning indices -- md5 verification is best-effort
	// integrity metadata, not required for the run to be usable.
	if (file_md5.size() == info.fastq_urls.size()) {
		info.fastq_md5 = std::move(file_md5);
	}

	out.push_back(std::move(info));
	return out;
}

} // namespace miint
