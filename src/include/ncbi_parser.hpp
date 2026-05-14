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

// epost handshake result: server-issued tokens that authorize a follow-up efetch.
// Empty fields signal a parse failure; callers (NCBIClient::EPostIds) wrap that
// into an IOException with the raw server response embedded for debuggability.
struct EPostResult {
	std::string webenv;
	std::string query_key;
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

	// Strip trailing ".N" version suffix from an accession, returning the base accession.
	// NCBI accepts unversioned IDs (e.g. "NC_000913") but its FASTA/XML responses include
	// versioned headers (e.g. ">NC_000913.3 ..."). Missing-accession detection has to
	// compare base IDs or it will report every unversioned request as missing.
	static std::string StripAccessionVersion(const std::string &accession);

	// Parse a multi-record GenBank XML response (<GBSet> wrapping multiple <GBSeq>) into
	// one GenBankMetadata per record, in submission order. The single-record ParseGenBankXML
	// silently returns only the first record when fed a batch response, so callers handling
	// batched efetch must use this entrypoint.
	static std::vector<GenBankMetadata> ParseGenBankXMLBatch(const std::string &xml);

	// Diff requested vs. returned accessions, comparing base IDs (after StripAccessionVersion).
	// Returns the subset of `requested` not present in `returned`, in original requested order.
	// NCBI silently omits invalid IDs from batch responses — without this diff the caller
	// has no way to surface lost data (Rule 10: fail loud).
	static std::vector<std::string> DiffMissingAccessions(const std::vector<std::string> &requested,
	                                                      const std::vector<std::string> &returned);

	// Extract content of a single XML tag (first occurrence). Handles attribute-bearing
	// open tags and nested same-name tags via depth counting. Exposed publicly so the
	// EPostIds handshake parser can reuse it for <WebEnv> and <QueryKey>.
	static std::string ExtractXMLTagValue(const std::string &xml, const std::string &tag);

	// Parse <WebEnv> and <QueryKey> from an NCBI epost response. Returns an
	// EPostResult with empty fields if either tag is absent — caller (NCBIClient::EPostIds)
	// is responsible for surfacing the failure with the raw response embedded.
	// Kept here, separate from NCBIClient, so the unit-test binary (which doesn't
	// link the DuckDB HTTP stack) can exercise the parsing in isolation.
	static EPostResult ParseEPostResponse(const std::string &xml);

	// Join a vector of strings with `sep`. Tiny utility, but the alternative
	// (ostringstream + i==0 guard) was duplicated half a dozen times across
	// the batched-fetch code paths; collapse to one site so a future format
	// change happens in one place.
	static std::string JoinStrings(const std::vector<std::string> &parts, const std::string &sep);
};

} // namespace miint
