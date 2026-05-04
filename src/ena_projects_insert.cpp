// SPDX-License-Identifier: MIT
//
// Pure-data helper for INSERT INTO ena.projects. See ena_projects_insert.hpp.

#include "ena_projects_insert.hpp"

#include "ena_receipt_parser.hpp"
#include "ena_submit_outcome.hpp"

#include <stdexcept>

namespace miint {

namespace {

struct ProjectSubmitTraits {
	static void SetEnvelopeArray(SubmissionSpec &env, const std::vector<ProjectSpec> &specs) {
		env.projects = specs;
	}
	static const char *ReceiptObjectType() {
		return "PROJECT";
	}
	static std::string BuildEnvelope(const SubmissionSpec &env) {
		return BuildEnvelopeJSON(env);
	}
	static const char *ContentType() {
		return "application/json";
	}
	static ENAProjectInsertResult BuildRow(const ProjectSpec &spec, const ENAObjectReceipt &obj) {
		ENAProjectInsertResult row;
		row.alias = spec.alias;
		row.prjeb_accession = obj.accession;
		row.status = obj.status;
		row.hold_until_date = obj.hold_until_date;
		for (const auto &ext : obj.ext_ids) {
			// Case-insensitive: the receipt XSD's EXT_ID/@type enumeration is a
			// free-form string and observed casings vary across endpoints.
			if (EqualsIgnoreCase(ext.type, "study")) {
				row.erp_accession = ext.accession;
				break;
			}
		}
		return row;
	}
};

} // namespace

ENASubmissionOutcome SubmitProjectInsertOutcome(const std::vector<ProjectSpec> &projects,
                                                const ENAProjectInsertOptions &opts, const ENAPostFn &post_fn) {
	return SubmitENAObjectOutcome<ProjectSubmitTraits, ProjectSpec, ENAProjectInsertOptions, ENASubmissionOutcome>(
	    projects, opts, post_fn);
}

std::vector<ENAProjectInsertResult> SubmitProjectInsert(const std::vector<ProjectSpec> &projects,
                                                        const ENAProjectInsertOptions &opts, const ENAPostFn &post_fn) {
	auto outcome = SubmitProjectInsertOutcome(projects, opts, post_fn);
	if (!outcome.success) {
		const std::string detail =
		    outcome.error_messages.empty() ? "no error detail" : FlattenENAErrors(outcome.error_messages);
		throw std::runtime_error("ENA submission failed: " + detail);
	}
	return outcome.rows;
}

} // namespace miint
