// SPDX-License-Identifier: MIT
// Pure-data helper for INSERT INTO ena.samples. See ena_samples_insert.hpp.

#include "ena_samples_insert.hpp"

#include "ena_receipt_parser.hpp"
#include "ena_submit_outcome.hpp"

#include <stdexcept>

namespace miint {

namespace {

struct SampleSubmitTraits {
	static void SetEnvelopeArray(SubmissionSpec &env, const std::vector<SampleSpec> &specs) {
		env.samples = specs;
	}
	static const char *ReceiptObjectType() {
		return "SAMPLE";
	}
	static std::string BuildEnvelope(const SubmissionSpec &env) {
		return BuildEnvelopeJSON(env);
	}
	static const char *ContentType() {
		return "application/json";
	}
	static ENASampleInsertResult BuildRow(const SampleSpec &spec, const ENAObjectReceipt &obj) {
		ENASampleInsertResult row;
		row.alias = spec.alias;
		row.ers_accession = obj.accession;
		row.status = obj.status;
		row.hold_until_date = obj.hold_until_date;
		for (const auto &ext : obj.ext_ids) {
			// Case-insensitive: the ENA receipt XSD enumerates this attribute
			// as a free-form string; observed casings include "biosample" (V2),
			// "BioSample" (some legacy paths), and "BIOSAMPLE" (rare).
			if (EqualsIgnoreCase(ext.type, "biosample")) {
				row.samea_accession = ext.accession;
				break;
			}
		}
		return row;
	}
};

} // namespace

ENASamplesSubmissionOutcome SubmitSampleInsertOutcome(const std::vector<SampleSpec> &samples,
                                                      const ENASampleInsertOptions &opts, const ENAPostFn &post_fn) {
	return SubmitENAObjectOutcome<SampleSubmitTraits, SampleSpec, ENASampleInsertOptions, ENASamplesSubmissionOutcome>(
	    samples, opts, post_fn);
}

std::vector<ENASampleInsertResult> SubmitSampleInsert(const std::vector<SampleSpec> &samples,
                                                      const ENASampleInsertOptions &opts, const ENAPostFn &post_fn) {
	auto outcome = SubmitSampleInsertOutcome(samples, opts, post_fn);
	if (!outcome.success) {
		const std::string detail =
		    outcome.error_messages.empty() ? "no error detail" : FlattenENAErrors(outcome.error_messages);
		throw std::runtime_error("ENA submission failed: " + detail);
	}
	return outcome.rows;
}

} // namespace miint
