#include "ena_search_field_registry.hpp"

#include <algorithm>
#include <cctype>

namespace miint {

// Curated from ENA Portal API /returnFields?result=sample on 2026-04-23.
// Subset chosen to cover the microbiome-analysis sample metadata that users
// of read_ena_attributes typically filter on (geography, taxonomy, host,
// environment, sampling context, physical properties, strain info).
//
// REGENERATION
// ------------
// 1. List the live field set:
//      SELECT field_name FROM ena_searchable_fields('sample') ORDER BY field_name;
// 2. Cross-reference with the plan's curation criteria in
//    localdocs/PLAN-ena-predicate-maxseqs.md (decision #8): keep only fields
//    useful for microbiome analysis (host_*, environment_*, collection_date,
//    country, lat, lon, depth, tax_id, scientific_name, sample_alias,
//    description, sampling_site, isolation_source, and related).
// 3. Preserve all entries already present here that still exist upstream —
//    removing one silently breaks user queries that pushed down successfully.
// 4. Paste the new lowercase names into the set below, update the date in
//    this comment, and bump any tests that enumerate the set size.
//
// Note: read_run / experiment / study field sets will get their own registries
// when pushdown is extended to other result types.
static const std::set<std::string> kSampleFields = {
    // Sample identifiers
    "sample_accession",
    "secondary_sample_accession",
    "sample_alias",
    "sample_title",
    "sample_description",
    "description",
    // Taxonomy
    "scientific_name",
    "tax_id",
    "tax_lineage",
    "target_gene",
    // Time
    "collection_date",
    "collection_date_start",
    "collection_date_end",
    "first_public",
    "last_updated",
    // Geography
    "country",
    "lat",
    "lon",
    "depth",
    "altitude",
    "elevation",
    // Environment (MIxS-style)
    "environment_biome",
    "environment_feature",
    "environment_material",
    "environmental_medium",
    "environmental_sample",
    "broad_scale_environmental_context",
    "local_environmental_context",
    // Host
    "host",
    "host_body_site",
    "host_scientific_name",
    "host_tax_id",
    "host_sex",
    "host_status",
    "host_phenotype",
    "host_genotype",
    // Sampling context
    "isolation_source",
    "sampling_site",
    "sampling_platform",
    "sampling_campaign",
    "collected_by",
    "identified_by",
    // Project / study linkage
    "project_name",
    "study_accession",
    "center_name",
    "broker_name",
    "investigation_type",
    "checklist",
    // Strain / variety
    "strain",
    "sub_species",
    "sub_strain",
    "cultivar",
    "ecotype",
    "serotype",
    "serovar",
    // Physical properties
    "ph",
    "salinity",
    "temperature",
    // Disease / status
    "disease",
    "sample_capture_status",
    // Misc descriptive
    "keywords",
};

static std::string CanonicalizeTag(const std::string &tag) {
	auto first = tag.find_first_not_of(" \t\n\r");
	if (first == std::string::npos) {
		return {};
	}
	auto last = tag.find_last_not_of(" \t\n\r");
	std::string trimmed = tag.substr(first, last - first + 1);
	std::transform(trimmed.begin(), trimmed.end(), trimmed.begin(),
	               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	return trimmed;
}

bool ENASearchFieldRegistry::IsSearchableSampleField(const std::string &tag) {
	auto canonical = CanonicalizeTag(tag);
	if (canonical.empty()) {
		return false;
	}
	return kSampleFields.count(canonical) > 0;
}

const std::set<std::string> &ENASearchFieldRegistry::SearchableSampleFields() {
	return kSampleFields;
}

} // namespace miint
