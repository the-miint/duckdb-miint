#pragma once

#include "duckdb/main/client_context.hpp"
#include <string>
#include <unordered_map>
#include <vector>

namespace duckdb {

// A single read's sequence and quality data (forward strand, as stored in FASTQ).
struct SequenceDataEntry {
	std::string sequence1;
	std::vector<uint8_t> qual1; // raw Phred (0-93), no ASCII offset
	std::string sequence2;      // empty if single-end
	std::vector<uint8_t> qual2;
};

using SequenceDataMap = std::unordered_map<std::string, SequenceDataEntry>;

// Schema info for a sequence data table.
struct SequenceDataSchema {
	bool has_sequence2 = false;
	bool has_qual2 = false;
};

// Validate that a table/view has the required read_fastx schema columns:
// Required: read_id (VARCHAR), sequence1 (VARCHAR), qual1 (LIST(UTINYINT)).
// Optional: sequence2 (VARCHAR), qual2 (LIST(UTINYINT)) — only needed for paired-end.
// Throws BinderException on missing required columns or wrong types.
SequenceDataSchema ValidateSequenceDataSchema(ClientContext &context, const std::string &table_name);

// Read all sequence data from a table/view into memory, keyed by read_id.
// Uses a separate Connection to avoid deadlocking the caller's context.
// Throws on duplicate read_id, NULL read_id, or NULL sequence1.
SequenceDataMap ReadSequenceDataTable(ClientContext &context, const std::string &table_name);

} // namespace duckdb
