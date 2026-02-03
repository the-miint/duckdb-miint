#include "zip_utils.hpp"
#include "miniz.hpp"
#include <cstring>

namespace miint {

std::string ExtractFromZip(const std::string &zip_data, const std::string &substring) {
	duckdb_miniz::mz_zip_archive zip;
	memset(&zip, 0, sizeof(zip));

	// Initialize ZIP reader from memory
	if (!duckdb_miniz::mz_zip_reader_init_mem(&zip, zip_data.data(), zip_data.size(), 0)) {
		throw std::runtime_error("Failed to read ZIP archive");
	}

	// Proper RAII guard to ensure ZIP reader cleanup on exception or return
	struct ZipReaderCleanup {
		duckdb_miniz::mz_zip_archive *zip_ptr;
		~ZipReaderCleanup() {
			duckdb_miniz::mz_zip_reader_end(zip_ptr);
		}
	} cleanup_guard {&zip};

	std::string result;
	bool found = false;

	// Iterate through files in ZIP
	duckdb_miniz::mz_uint num_files = duckdb_miniz::mz_zip_reader_get_num_files(&zip);
	for (duckdb_miniz::mz_uint i = 0; i < num_files; i++) {
		duckdb_miniz::mz_zip_archive_file_stat file_stat;
		if (!duckdb_miniz::mz_zip_reader_file_stat(&zip, i, &file_stat)) {
			continue;
		}

		std::string filename = file_stat.m_filename;

		// Check if filename contains substring (e.g., "_genomic.fna")
		// Note: This is substring matching, not glob/regex pattern matching
		if (filename.find(substring) != std::string::npos) {
			found = true;

			// Get the reliable uncompressed size from the file_stat we already have.
			// This is read from the ZIP central directory header and is authoritative.
			// miniz's extract_to_heap may allocate a larger buffer for alignment/padding,
			// so we must use this size to avoid including extra null bytes.
			size_t expected_size = static_cast<size_t>(file_stat.m_uncomp_size);

			// Pre-allocate result string and extract directly into it.
			// This avoids malloc/free assumptions about miniz's allocator.
			result.resize(expected_size);
			if (expected_size > 0) {
				if (!duckdb_miniz::mz_zip_reader_extract_to_mem(&zip, i, result.data(), expected_size, 0)) {
					throw std::runtime_error(std::string("Failed to extract '") + filename +
					                         "' from ZIP: " + duckdb_miniz::mz_zip_get_error_string(zip.m_last_error));
				}
			}
			break;
		}
	}

	if (!found) {
		throw std::runtime_error("No file containing '" + substring + "' found in ZIP archive (" +
		                         std::to_string(num_files) + " files checked)");
	}

	return result;
}

} // namespace miint
