#pragma once

#include "SequenceRecord.hpp"
#include <functional>
#include <string>
#include <vector>

namespace miint {

// Warning callback type for parser warnings
using WarningCallback = std::function<void(const std::string &)>;

// Default warning callback that writes to stderr
void DefaultWarningCallback(const std::string &msg);

// Accession type for routing to appropriate API
enum class AccessionType {
	ASSEMBLY, // GCF_/GCA_ -> Datasets API for metadata
	SEQUENCE, // NC_/NM_/NP_/etc -> E-utilities
	UNKNOWN
};

// Parsed GenBank metadata for read_ncbi
struct GenBankMetadata {
	std::string accession;
	int32_t version = 0;
	std::string description;
	std::string organism;
	int64_t taxonomy_id = 0;
	int64_t length = 0;
	std::string molecule_type;
	std::string update_date; // YYYY-MM-DD format
	std::string extra_json;  // JSON string for additional fields
};

// Parsed feature annotation for read_ncbi_annotation (GFF3-compatible)
struct FeatureAnnotation {
	std::string seqid;                                      // Sequence ID
	std::string source;                                     // Source (e.g., "RefSeq", "GenBank")
	std::string type;                                       // Feature type (gene, CDS, etc.)
	int64_t position = 0;                                   // Start position (1-based)
	int64_t stop_position = 0;                              // End position
	double score = 0.0;                                     // Score (or NaN if not applicable)
	bool has_score = false;                                 // Whether score is valid
	std::string strand;                                     // + or - or .
	int32_t phase = -1;                                     // 0, 1, 2, or -1 if not applicable
	std::vector<std::pair<std::string, std::string>> attrs; // Key-value attributes
};

// Batch of feature annotations
struct FeatureAnnotationBatch {
	std::vector<FeatureAnnotation> features;

	bool empty() const {
		return features.empty();
	}
	size_t size() const {
		return features.size();
	}
};

// NCBI data parsing utilities (no DuckDB dependencies)
class NCBIParser {
public:
	// Accession type detection
	static AccessionType DetectAccessionType(const std::string &accession);
	static bool IsAssemblyAccession(const std::string &accession);

	// Parse GenBank XML into metadata struct
	// Returns empty metadata if parsing fails
	static GenBankMetadata ParseGenBankXML(const std::string &xml);

	// Parse FASTA text into SequenceRecordBatch (reuses existing struct)
	// Returns batch with is_paired=false, empty quals
	// For NCBI FASTA, extracts accession as read_id (handles pipe-delimited format)
	static SequenceRecordBatch ParseFasta(const std::string &fasta_text);

	// Parse INSDC feature table into annotation batch (GFF3-compatible)
	// Limitations:
	// - Complex locations (join, complement of join) emit warnings and use outer bounds
	// - Phase is derived from codon_start qualifier for CDS features
	// Set warn_callback to nullptr or empty function for quiet mode
	static FeatureAnnotationBatch ParseFeatureTable(const std::string &feature_table,
	                                                WarningCallback warn_callback = DefaultWarningCallback);

	// Helper to extract accession from various FASTA ID formats
	// Handles: "NC_001416.1", "ref|NC_001416.1|", "gi|123|ref|NC_001416.1|description"
	static std::string ExtractAccessionFromFastaId(const std::string &fasta_id, std::string &remainder);
};

} // namespace miint
