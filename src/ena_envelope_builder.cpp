// ENA Webin envelope builder. See ena_envelope_builder.hpp.

#include "ena_envelope_builder.hpp"

#include <algorithm>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace miint {

namespace {

// Compile-time enum sets for client-side validation. Sourced from the SRA
// XSDs in localdocs/ena-research-webin-v2-deep.md §4.4 and §6.1. Membership
// is rough — the goal is to reject typos before the network round-trip, not
// to be an exhaustive compatibility matrix. ENA's server is the authority.
const std::unordered_set<std::string> &LibraryStrategies() {
	static const std::unordered_set<std::string> kSet = {"WGS",
	                                                     "WGA",
	                                                     "WXS",
	                                                     "RNA-Seq",
	                                                     "ssRNA-seq",
	                                                     "miRNA-Seq",
	                                                     "ncRNA-Seq",
	                                                     "FL-cDNA",
	                                                     "EST",
	                                                     "Hi-C",
	                                                     "ATAC-seq",
	                                                     "WCS",
	                                                     "RAD-Seq",
	                                                     "CLONE",
	                                                     "POOLCLONE",
	                                                     "AMPLICON",
	                                                     "CLONEEND",
	                                                     "FINISHING",
	                                                     "ChIP-Seq",
	                                                     "MNase-Seq",
	                                                     "DNase-Hypersensitivity",
	                                                     "Bisulfite-Seq",
	                                                     "CTS",
	                                                     "MRE-Seq",
	                                                     "MeDIP-Seq",
	                                                     "MBD-Seq",
	                                                     "Tn-Seq",
	                                                     "VALIDATION",
	                                                     "FAIRE-seq",
	                                                     "SELEX",
	                                                     "RIP-Seq",
	                                                     "ChIA-PET",
	                                                     "Synthetic-Long-Read",
	                                                     "Targeted-Capture",
	                                                     "Tethered Chromatin Conformation Capture",
	                                                     "NOMe-Seq",
	                                                     "ChM-Seq",
	                                                     "GBS",
	                                                     "Ribo-Seq",
	                                                     "OTHER"};
	return kSet;
}

const std::unordered_set<std::string> &LibrarySources() {
	static const std::unordered_set<std::string> kSet = {
	    "GENOMIC",     "GENOMIC SINGLE CELL", "TRANSCRIPTOMIC", "TRANSCRIPTOMIC SINGLE CELL",
	    "METAGENOMIC", "METATRANSCRIPTOMIC",  "SYNTHETIC",      "VIRAL RNA",
	    "OTHER"};
	return kSet;
}

const std::unordered_set<std::string> &LibrarySelections() {
	static const std::unordered_set<std::string> kSet = {"RANDOM",
	                                                     "PCR",
	                                                     "RANDOM PCR",
	                                                     "RT-PCR",
	                                                     "HMPR",
	                                                     "MF",
	                                                     "repeat fractionation",
	                                                     "size fractionation",
	                                                     "MSLL",
	                                                     "cDNA",
	                                                     "cDNA_randomPriming",
	                                                     "cDNA_oligo_dT",
	                                                     "PolyA",
	                                                     "Oligo-dT",
	                                                     "Inverse rRNA",
	                                                     "Inverse rRNA selection",
	                                                     "ChIP",
	                                                     "ChIP-Seq",
	                                                     "MNase",
	                                                     "DNase",
	                                                     "Hybrid Selection",
	                                                     "Reduced Representation",
	                                                     "Restriction Digest",
	                                                     "5-methylcytidine antibody",
	                                                     "MBD2 protein methyl-CpG binding domain",
	                                                     "CAGE",
	                                                     "RACE",
	                                                     "MDA",
	                                                     "padlock probes capture method",
	                                                     "other",
	                                                     "unspecified"};
	return kSet;
}

const std::unordered_set<std::string> &Platforms() {
	static const std::unordered_set<std::string> kSet = {
	    "LS454",     "ILLUMINA",          "OXFORD_NANOPORE", "PACBIO_SMRT", "ION_TORRENT", "BGISEQ",
	    "ABI_SOLID", "COMPLETE_GENOMICS", "HELICOS",         "CAPILLARY",   "DNBSEQ"};
	return kSet;
}

const std::unordered_set<std::string> &RunFiletypes() {
	static const std::unordered_set<std::string> kSet = {"sra",
	                                                     "srf",
	                                                     "sff",
	                                                     "fastq",
	                                                     "fasta",
	                                                     "tab",
	                                                     "454_native",
	                                                     "454_native_seq",
	                                                     "454_native_qual",
	                                                     "Helicos_native",
	                                                     "Illumina_native",
	                                                     "Illumina_native_seq",
	                                                     "Illumina_native_prb",
	                                                     "Illumina_native_int",
	                                                     "Illumina_native_qseq",
	                                                     "Illumina_native_scarf",
	                                                     "SOLiD_native",
	                                                     "SOLiD_native_csfasta",
	                                                     "SOLiD_native_qual",
	                                                     "PacBio_HDF5",
	                                                     "bam",
	                                                     "cram",
	                                                     "CompleteGenomics_native",
	                                                     "OxfordNanopore_native"};
	return kSet;
}

void RequireMembership(const std::unordered_set<std::string> &set, const std::string &value, const char *field_name,
                       const std::string &alias) {
	if (set.find(value) == set.end()) {
		throw std::runtime_error("ENA envelope: unknown " + std::string(field_name) + " '" + value + "' for alias '" +
		                         alias + "'");
	}
}

// JSON string escaping per RFC 8259 §7. Treats input as opaque bytes; multibyte
// UTF-8 characters pass through unchanged. Control bytes (< 0x20) are escaped
// as \uXXXX; named escapes used where defined.
void AppendJsonString(std::string &out, const std::string &s) {
	out.push_back('"');
	for (auto c : s) {
		const auto u = static_cast<unsigned char>(c);
		switch (u) {
		case '"':
			out.append("\\\"");
			break;
		case '\\':
			out.append("\\\\");
			break;
		case '\b':
			out.append("\\b");
			break;
		case '\f':
			out.append("\\f");
			break;
		case '\n':
			out.append("\\n");
			break;
		case '\r':
			out.append("\\r");
			break;
		case '\t':
			out.append("\\t");
			break;
		default:
			if (u < 0x20) {
				char buf[8];
				std::snprintf(buf, sizeof(buf), "\\u%04x", u);
				out.append(buf);
			} else {
				out.push_back(c);
			}
		}
	}
	out.push_back('"');
}

const char *ActionName(ENAAction a) {
	switch (a) {
	case ENAAction::ADD:
		return "ADD";
	case ENAAction::MODIFY:
		return "MODIFY";
	case ENAAction::CANCEL:
		return "CANCEL";
	case ENAAction::HOLD:
		return "HOLD";
	case ENAAction::RELEASE:
		return "RELEASE";
	case ENAAction::VALIDATE:
		return "VALIDATE";
	}
	throw std::logic_error("ENA envelope: unhandled ENAAction value");
}

// True iff `s` contains at least one non-whitespace byte. Used to reject
// targets like "   " that would otherwise pass an `.empty()` check and emit
// `<CANCEL target="   "/>` to the server.
bool HasNonWhitespace(const std::string &s) {
	for (char c : s) {
		if (c != ' ' && c != '\t' && c != '\r' && c != '\n') {
			return true;
		}
	}
	return false;
}

// Single source of truth for "does this spec carry a target accession or
// refname?". Empty / whitespace-only fields don't count as targets.
bool IsTargetedAction(const SubmissionSpec &env) {
	return HasNonWhitespace(env.target_accession) || HasNonWhitespace(env.target_refname);
}

// True iff the spec carries any object-body content.
bool HasBodyContent(const SubmissionSpec &env) {
	return !env.projects.empty() || !env.samples.empty() || !env.experiments.empty() || !env.runs.empty();
}

// Pure-data validators shared by the JSON and XML emitters. Each Append*
// pair (Append<X> for JSON, AppendXml<X> for XML) calls the matching
// ValidateXxx at the top so the precondition string lives in one place.
void ValidateActions(const SubmissionSpec &env) {
	const bool has_target = IsTargetedAction(env);

	// Reject whitespace-only targets that pass the `.empty()` check but would
	// emit garbage to the server. Caught here rather than at the call site so
	// callers get a clear message instead of an opaque server-side rejection.
	if (!env.target_accession.empty() && !HasNonWhitespace(env.target_accession)) {
		throw std::runtime_error("ENA envelope: target_accession must contain non-whitespace characters");
	}
	if (!env.target_refname.empty() && !HasNonWhitespace(env.target_refname)) {
		throw std::runtime_error("ENA envelope: target_refname must contain non-whitespace characters");
	}

	// Lifecycle actions that operate on an existing accession need a target.
	if ((env.action == ENAAction::CANCEL || env.action == ENAAction::RELEASE) && !has_target) {
		throw std::runtime_error(std::string("ENA envelope: ") + ActionName(env.action) +
		                         " requires target_accession or target_refname");
	}

	// HOLD has two distinct shapes:
	//   1. Body-pattern (forward-dated submission): action=ADD with
	//      hold_until_date set; the emitter pairs ADD with a sibling
	//      <HOLD HoldUntilDate=.../> action automatically.
	//   2. Targeted-pattern (post-hoc embargo extension): action=HOLD with a
	//      target accession plus hold_until_date.
	if (env.action == ENAAction::HOLD) {
		if (!has_target && env.hold_until_date.empty()) {
			throw std::runtime_error("ENA envelope: HOLD action requires hold_until_date");
		}
		if (!has_target && !env.hold_until_date.empty()) {
			throw std::runtime_error(
			    "ENA envelope: with hold_until_date, use action=ADD; the HOLD entry is added automatically");
		}
		if (has_target && env.hold_until_date.empty()) {
			throw std::runtime_error("ENA envelope: HOLD with a target requires hold_until_date "
			                         "(the new embargo end date, e.g. '2027-12-31')");
		}
	}

	// ADD/MODIFY/VALIDATE identify objects via the body sets, not via the
	// action element. A target on these is almost certainly a programming
	// mistake; reject it loudly.
	if ((env.action == ENAAction::ADD || env.action == ENAAction::MODIFY || env.action == ENAAction::VALIDATE) &&
	    has_target) {
		throw std::runtime_error(std::string("ENA envelope: ") + ActionName(env.action) +
		                         " action does not take a target accession or refname");
	}

	// Targeted lifecycle actions carry no body — the target is the entire
	// payload. Reject body content rather than silently dropping it: a caller
	// reusing a populated SubmissionSpec for a CANCEL would otherwise look
	// like a successful submit-then-cancel, when really the body was discarded.
	if (has_target && HasBodyContent(env)) {
		throw std::runtime_error(std::string("ENA envelope: ") + ActionName(env.action) +
		                         " with target_accession/target_refname must not carry body content "
		                         "(projects/samples/experiments/runs); the target is the entire payload");
	}

	// hold_until_date is only meaningful with ADD (forward-dated submission)
	// or HOLD (post-hoc embargo). Setting it on CANCEL/RELEASE/MODIFY/VALIDATE
	// is silently ignored by the wire format, which would otherwise let a
	// caller mistakenly populate `outcome.hold_until_date` and then write
	// it to submission_log even though it had no effect on the server.
	if (!env.hold_until_date.empty() && env.action != ENAAction::ADD && env.action != ENAAction::HOLD) {
		throw std::runtime_error(std::string("ENA envelope: ") + ActionName(env.action) +
		                         " action does not take a hold_until_date "
		                         "(only ADD and HOLD do)");
	}
}

void ValidateExperimentSpec(const ExperimentSpec &e) {
	if (e.alias.empty()) {
		throw std::runtime_error("ENA envelope: experiment alias must be non-empty");
	}
	RequireMembership(LibraryStrategies(), e.library_strategy, "library_strategy", e.alias);
	RequireMembership(LibrarySources(), e.library_source, "library_source", e.alias);
	RequireMembership(LibrarySelections(), e.library_selection, "library_selection", e.alias);
	RequireMembership(Platforms(), e.platform, "platform", e.alias);
	if (e.instrument_model.empty()) {
		throw std::runtime_error("ENA envelope: experiment '" + e.alias + "' instrument_model must be non-empty");
	}
}

void ValidateRunFile(const RunFile &f, const std::string &run_alias) {
	if (f.filename.empty()) {
		throw std::runtime_error("ENA envelope: run '" + run_alias + "' file.filename must be non-empty");
	}
	if (f.checksum.empty()) {
		throw std::runtime_error("ENA envelope: run '" + run_alias + "' file '" + f.filename +
		                         "' checksum must be non-empty");
	}
	if (f.filetype.empty()) {
		throw std::runtime_error("ENA envelope: run '" + run_alias + "' file '" + f.filename +
		                         "' filetype must be non-empty (e.g. 'fastq', 'bam', 'cram')");
	}
	RequireMembership(RunFiletypes(), f.filetype, "run.file.filetype", run_alias);
}

void ValidateRunSpec(const RunSpec &r) {
	if (r.alias.empty()) {
		throw std::runtime_error("ENA envelope: run alias must be non-empty");
	}
	if (r.files.empty()) {
		throw std::runtime_error("ENA envelope: run '" + r.alias + "' must have at least one file");
	}
}

void AppendActions(std::string &out, const SubmissionSpec &env) {
	ValidateActions(env);
	// Targeted lifecycle actions are XML-only in the current build. The V2
	// JSON dispatcher's behaviour for submission-level actions with `target=`
	// hasn't been verified live, and we don't want a silently-malformed
	// envelope (no `target` field) to slip out via the JSON path.
	if (IsTargetedAction(env)) {
		throw std::runtime_error("ENA envelope: targeted lifecycle actions (CANCEL, RELEASE, "
		                         "targeted HOLD) must be built via BuildEnvelopeXML, not JSON");
	}
	out.append("\"actions\":[");
	out.append("{\"type\":");
	AppendJsonString(out, ActionName(env.action));
	out.push_back('}');
	if (!env.hold_until_date.empty()) {
		out.append(",{\"type\":\"HOLD\",\"holdUntilDate\":");
		AppendJsonString(out, env.hold_until_date);
		out.push_back('}');
	}
	out.push_back(']');
}

void AppendProject(std::string &out, const ProjectSpec &p) {
	if (p.alias.empty()) {
		throw std::runtime_error("ENA envelope: project alias must be non-empty");
	}
	out.push_back('{');
	out.append("\"alias\":");
	AppendJsonString(out, p.alias);
	out.append(",\"title\":");
	AppendJsonString(out, p.title);
	// `description` is always emitted: wwwdev's XSD validator intermittently
	// rejects PROJECT documents missing a populated <DESCRIPTION> element
	// ("Expected element 'DESCRIPTION' before the end of the content in
	// element PROJECT") even though SRA.project.xsd declares it minOccurs=0.
	// Behaviour is non-deterministic across runs (one submission accepted,
	// the next rejected, same alias shape), so we always send it. When the
	// user didn't provide one we reuse `title` so the resulting element has
	// non-empty content.
	out.append(",\"description\":");
	AppendJsonString(out, p.description.empty() ? p.title : p.description);
	out.append(p.is_umbrella ? ",\"umbrellaProject\":{}" : ",\"sequencingProject\":{}");
	out.push_back('}');
}

void AppendSample(std::string &out, const SampleSpec &s) {
	if (s.alias.empty()) {
		throw std::runtime_error("ENA envelope: sample alias must be non-empty");
	}
	if (s.taxon_id <= 0) {
		throw std::runtime_error("ENA envelope: sample.taxon_id must be > 0 (got " + std::to_string(s.taxon_id) +
		                         " for alias '" + s.alias + "')");
	}
	out.push_back('{');
	out.append("\"alias\":");
	AppendJsonString(out, s.alias);
	if (!s.title.empty()) {
		out.append(",\"title\":");
		AppendJsonString(out, s.title);
	}
	if (!s.description.empty()) {
		out.append(",\"description\":");
		AppendJsonString(out, s.description);
	}
	out.append(",\"organism\":{\"taxonId\":");
	AppendJsonString(out, std::to_string(s.taxon_id));
	if (!s.scientific_name.empty()) {
		out.append(",\"scientificName\":");
		AppendJsonString(out, s.scientific_name);
	}
	out.push_back('}');

	const bool any_attrs = !s.checklist.empty() || !s.attributes.empty();
	if (any_attrs) {
		out.append(",\"attributes\":[");
		bool first = true;
		if (!s.checklist.empty()) {
			out.append("{\"tag\":\"ENA-CHECKLIST\",\"value\":");
			AppendJsonString(out, s.checklist);
			out.push_back('}');
			first = false;
		}
		for (const auto &kv : s.attributes) {
			if (!first) {
				out.push_back(',');
			}
			out.append("{\"tag\":");
			AppendJsonString(out, kv.first);
			out.append(",\"value\":");
			AppendJsonString(out, kv.second);
			// Optional `unit` from the sparse attribute_units vector — some
			// checklist attributes (e.g. ERC000015 lat/lon → `DD`) are
			// rejected by the server without it. JSON key is singular
			// (`unit`) even though the SRA XML element is `<UNITS>` and the
			// SQL column is plural (`attribute_units`); this is what the V2
			// JSON validator expects. Linear lookup is fine; the vector is
			// small per sample (typically 0–3 entries) and tag strings are
			// short.
			for (const auto &u : s.attribute_units) {
				if (u.first == kv.first) {
					if (!u.second.empty()) {
						out.append(",\"unit\":");
						AppendJsonString(out, u.second);
					}
					break;
				}
			}
			out.push_back('}');
			first = false;
		}
		out.push_back(']');
	}
	out.push_back('}');
}

// Emit a {"refname":"..."} or {"accession":"..."} cross-reference. accession
// wins when both are set. `field_name` and `alias` are only used for the
// "missing both" error message.
void AppendRef(std::string &out, const RefDescriptor &ref, const char *field_name, const std::string &alias) {
	if (ref.accession.empty() && ref.refname.empty()) {
		throw std::runtime_error("ENA envelope: " + std::string(field_name) + " required for alias '" + alias + "'");
	}
	out.push_back('{');
	if (!ref.accession.empty()) {
		out.append("\"accession\":");
		AppendJsonString(out, ref.accession);
	} else {
		out.append("\"refname\":");
		AppendJsonString(out, ref.refname);
	}
	out.push_back('}');
}

void AppendExperiment(std::string &out, const ExperimentSpec &e) {
	ValidateExperimentSpec(e);

	out.push_back('{');
	out.append("\"alias\":");
	AppendJsonString(out, e.alias);
	if (!e.title.empty()) {
		out.append(",\"title\":");
		AppendJsonString(out, e.title);
	}
	out.append(",\"studyRef\":");
	AppendRef(out, e.study_ref, "study_ref", e.alias);
	out.append(",\"design\":{");
	if (!e.design_description.empty()) {
		out.append("\"designDescription\":");
		AppendJsonString(out, e.design_description);
		out.push_back(',');
	}
	out.append("\"sampleDescriptor\":");
	AppendRef(out, e.sample_ref, "sample_ref", e.alias);
	out.append(",\"libraryDescriptor\":{");
	if (!e.library_name.empty()) {
		out.append("\"libraryName\":");
		AppendJsonString(out, e.library_name);
		out.push_back(',');
	}
	out.append("\"libraryStrategy\":");
	AppendJsonString(out, e.library_strategy);
	out.append(",\"librarySource\":");
	AppendJsonString(out, e.library_source);
	out.append(",\"librarySelection\":");
	AppendJsonString(out, e.library_selection);
	out.append(",\"libraryLayout\":");
	out.append(e.library_layout == ENALibraryLayout::PAIRED ? "{\"paired\":{}}" : "{\"single\":{}}");
	out.append("}}"); // close libraryDescriptor + design

	out.append(",\"platform\":{");
	AppendJsonString(out, e.platform);
	out.append(":{\"instrumentModel\":");
	AppendJsonString(out, e.instrument_model);
	out.append("}}"); // close platform.<NAME> + platform

	out.push_back('}');
}

void AppendRunFile(std::string &out, const RunFile &f, const std::string &run_alias) {
	ValidateRunFile(f, run_alias);

	out.push_back('{');
	out.append("\"filename\":");
	AppendJsonString(out, f.filename);
	out.append(",\"filetype\":");
	AppendJsonString(out, f.filetype);
	out.append(",\"checksumMethod\":\"MD5\",\"checksum\":");
	AppendJsonString(out, f.checksum);
	out.push_back('}');
}

void AppendRun(std::string &out, const RunSpec &r) {
	ValidateRunSpec(r);

	out.push_back('{');
	out.append("\"alias\":");
	AppendJsonString(out, r.alias);
	if (!r.title.empty()) {
		out.append(",\"title\":");
		AppendJsonString(out, r.title);
	}
	out.append(",\"experimentRef\":");
	AppendRef(out, r.experiment_ref, "experiment_ref", r.alias);
	out.append(",\"files\":[");
	bool first = true;
	for (const auto &f : r.files) {
		if (!first) {
			out.push_back(',');
		}
		AppendRunFile(out, f, r.alias);
		first = false;
	}
	out.push_back(']');
	out.push_back('}');
}

template <typename T, typename Appender>
void AppendArray(std::string &out, bool &needs_comma, const char *key, const std::vector<T> &items, Appender appender) {
	if (items.empty()) {
		return;
	}
	if (needs_comma) {
		out.push_back(',');
	}
	out.push_back('"');
	out.append(key);
	out.append("\":[");
	bool first = true;
	for (const auto &item : items) {
		if (!first) {
			out.push_back(',');
		}
		appender(out, item);
		first = false;
	}
	out.push_back(']');
	needs_comma = true;
}

} // namespace

// =====================================================================
// XML envelope (experiments + runs path)
// =====================================================================
//
// The V2 server's JSON dispatcher is implemented for project + sample only;
// SRA-side objects (experiment, run, analysis) require XML. We emit the
// canonical SRA-XSD shape wrapped in a `<WEBIN>` document so a single POST
// carries both `<SUBMISSION>` (action) and the requested object set.

namespace {

// XML attribute / element-text escaping per https://www.w3.org/TR/xml/#syntax
// (single quote only needed when the surrounding attribute is single-quoted,
// but we escape all five for safety regardless of context).
void AppendXmlEscaped(std::string &out, const std::string &s) {
	for (char c : s) {
		switch (c) {
		case '&':
			out.append("&amp;");
			break;
		case '<':
			out.append("&lt;");
			break;
		case '>':
			out.append("&gt;");
			break;
		case '"':
			out.append("&quot;");
			break;
		case '\'':
			out.append("&apos;");
			break;
		default:
			out.push_back(c);
		}
	}
}

void AppendXmlElement(std::string &out, const char *tag, const std::string &value) {
	out.push_back('<');
	out.append(tag);
	out.push_back('>');
	AppendXmlEscaped(out, value);
	out.append("</");
	out.append(tag);
	out.push_back('>');
}

// Emit `<TAG attr="value"/>` for a study/sample/experiment cross-reference.
// Uses `accession` when the descriptor has one, otherwise `refname`. Throws
// when both are empty (matches the JSON path's invariant).
void AppendXmlRef(std::string &out, const char *element, const RefDescriptor &ref, const char *field_name,
                  const std::string &alias) {
	if (ref.accession.empty() && ref.refname.empty()) {
		throw std::runtime_error("ENA envelope: " + std::string(field_name) + " required for alias '" + alias + "'");
	}
	out.push_back('<');
	out.append(element);
	out.push_back(' ');
	out.append(ref.accession.empty() ? "refname" : "accession");
	out.append("=\"");
	AppendXmlEscaped(out, ref.accession.empty() ? ref.refname : ref.accession);
	out.append("\"/>");
}

void AppendXmlActions(std::string &out, const SubmissionSpec &env) {
	ValidateActions(env);
	if (IsTargetedAction(env)) {
		// Targeted lifecycle action — accession wins over refname when both
		// are set, matching the RefDescriptor convention used elsewhere.
		const std::string &target = HasNonWhitespace(env.target_accession) ? env.target_accession : env.target_refname;
		out.append("<ACTIONS><ACTION><");
		out.append(ActionName(env.action));
		out.append(" target=\"");
		AppendXmlEscaped(out, target);
		if (env.action == ENAAction::HOLD) {
			out.append("\" HoldUntilDate=\"");
			AppendXmlEscaped(out, env.hold_until_date);
		}
		out.append("\"/></ACTION></ACTIONS>");
		return;
	}

	out.append("<ACTIONS><ACTION><");
	out.append(ActionName(env.action));
	out.append("/></ACTION>");
	if (!env.hold_until_date.empty()) {
		out.append("<ACTION><HOLD HoldUntilDate=\"");
		AppendXmlEscaped(out, env.hold_until_date);
		out.append("\"/></ACTION>");
	}
	out.append("</ACTIONS>");
}

void AppendXmlExperiment(std::string &out, const ExperimentSpec &e) {
	ValidateExperimentSpec(e);

	out.append("<EXPERIMENT alias=\"");
	AppendXmlEscaped(out, e.alias);
	out.append("\">");
	if (!e.title.empty()) {
		AppendXmlElement(out, "TITLE", e.title);
	}
	AppendXmlRef(out, "STUDY_REF", e.study_ref, "study_ref", e.alias);
	out.append("<DESIGN>");
	// DESIGN_DESCRIPTION is XSD-mandatory but may be empty; emit empty when
	// the user didn't provide one (matches webin-cli behaviour).
	out.append("<DESIGN_DESCRIPTION>");
	AppendXmlEscaped(out, e.design_description);
	out.append("</DESIGN_DESCRIPTION>");
	AppendXmlRef(out, "SAMPLE_DESCRIPTOR", e.sample_ref, "sample_ref", e.alias);
	out.append("<LIBRARY_DESCRIPTOR>");
	if (!e.library_name.empty()) {
		AppendXmlElement(out, "LIBRARY_NAME", e.library_name);
	}
	AppendXmlElement(out, "LIBRARY_STRATEGY", e.library_strategy);
	AppendXmlElement(out, "LIBRARY_SOURCE", e.library_source);
	AppendXmlElement(out, "LIBRARY_SELECTION", e.library_selection);
	out.append("<LIBRARY_LAYOUT>");
	out.append(e.library_layout == ENALibraryLayout::PAIRED ? "<PAIRED/>" : "<SINGLE/>");
	out.append("</LIBRARY_LAYOUT>");
	out.append("</LIBRARY_DESCRIPTOR>");
	out.append("</DESIGN>");
	out.append("<PLATFORM><");
	out.append(e.platform);
	out.append("><INSTRUMENT_MODEL>");
	AppendXmlEscaped(out, e.instrument_model);
	out.append("</INSTRUMENT_MODEL></");
	out.append(e.platform);
	out.append("></PLATFORM>");
	out.append("</EXPERIMENT>");
}

void AppendXmlRunFile(std::string &out, const RunFile &f, const std::string &run_alias) {
	ValidateRunFile(f, run_alias);

	out.append("<FILE filename=\"");
	AppendXmlEscaped(out, f.filename);
	out.append("\" filetype=\"");
	AppendXmlEscaped(out, f.filetype);
	out.append("\" checksum_method=\"MD5\" checksum=\"");
	AppendXmlEscaped(out, f.checksum);
	out.append("\"/>");
}

void AppendXmlRun(std::string &out, const RunSpec &r) {
	ValidateRunSpec(r);
	out.append("<RUN alias=\"");
	AppendXmlEscaped(out, r.alias);
	out.append("\">");
	if (!r.title.empty()) {
		AppendXmlElement(out, "TITLE", r.title);
	}
	AppendXmlRef(out, "EXPERIMENT_REF", r.experiment_ref, "experiment_ref", r.alias);
	out.append("<DATA_BLOCK><FILES>");
	for (const auto &f : r.files) {
		AppendXmlRunFile(out, f, r.alias);
	}
	out.append("</FILES></DATA_BLOCK>");
	out.append("</RUN>");
}

} // namespace

std::string BuildEnvelopeXML(const SubmissionSpec &env) {
	std::string out;
	out.reserve(512);
	out.append(R"(<?xml version="1.0" encoding="UTF-8"?>)");
	out.append("<WEBIN>");
	out.append("<SUBMISSION>");
	AppendXmlActions(out, env);
	out.append("</SUBMISSION>");
	// ValidateActions (called from AppendXmlActions above) rejects the
	// targeted-action + body-content combination, so reaching here with body
	// content implies an untargeted action.
	if (!env.experiments.empty()) {
		out.append("<EXPERIMENT_SET>");
		for (const auto &e : env.experiments) {
			AppendXmlExperiment(out, e);
		}
		out.append("</EXPERIMENT_SET>");
	}
	if (!env.runs.empty()) {
		out.append("<RUN_SET>");
		for (const auto &r : env.runs) {
			AppendXmlRun(out, r);
		}
		out.append("</RUN_SET>");
	}
	out.append("</WEBIN>");
	return out;
}

std::string BuildEnvelopeJSON(const SubmissionSpec &env) {
	std::string out;
	out.reserve(256);
	out.push_back('{');
	out.append("\"submission\":{");
	AppendActions(out, env);
	out.push_back('}');

	bool needs_comma = true;
	AppendArray(out, needs_comma, "projects", env.projects, AppendProject);
	AppendArray(out, needs_comma, "samples", env.samples, AppendSample);
	AppendArray(out, needs_comma, "experiments", env.experiments, AppendExperiment);
	AppendArray(out, needs_comma, "runs", env.runs, AppendRun);

	out.push_back('}');
	return out;
}

} // namespace miint
