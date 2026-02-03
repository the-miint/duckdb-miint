#pragma once

#include <stdexcept>
#include <string>

namespace miint {

// ZIP extraction utility
// Extracts the first file whose name contains 'substring' from a ZIP archive in memory.
// For example, substring="_genomic.fna" matches "ncbi_dataset/data/GCF_000001/GCF_000001_genomic.fna"
//
// Returns the extracted file content as a string (may be empty if the file is empty).
// Throws std::runtime_error if:
//   - ZIP archive cannot be read
//   - No file containing the substring is found
//   - Extraction fails
//
// Note: This performs substring matching, not glob or regex pattern matching.
std::string ExtractFromZip(const std::string &zip_data, const std::string &substring);

} // namespace miint
