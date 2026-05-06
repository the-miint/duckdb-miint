// SPDX-License-Identifier: MIT
//
// Shared scaffolding for the four test_ena_*_insert.cpp files. CapturedPost
// records the last POST observed by the mock post functor; MakeReceiptXML
// builds a canned receipt XML body from a list of (element, alias, accession,
// ext_id?) entries.

#pragma once

#include "ena_post_fn.hpp"

#include <string>
#include <vector>

namespace miint_test {

// Records the last POST so tests can inspect URL / body / auth / content-type.
struct CapturedPost {
	std::string url;
	std::string body;
	std::string user;
	std::string password;
	std::string content_type;
};

// Build a post functor that captures the request into `captured` and
// returns a fixed `response`. Tests inspect `captured` after the call.
// Use this for tests where the response body is fixed independent of
// request shape; the existing per-table tests build their response inline
// because they want it to depend on the captured aliases.
inline miint::ENAPostFn StubPost(CapturedPost &captured, const std::string &response) {
	return [&captured, response](const std::string &url, const std::string &body, const std::string &user,
	                             const std::string &password, const std::string &content_type) {
		captured = {url, body, user, password, content_type};
		return response;
	};
}

struct ReceiptObjectFixture {
	std::string element; // PROJECT / SAMPLE / EXPERIMENT / RUN
	std::string alias;
	std::string accession;        // server-assigned primary (PRJEBxxx, ERSxxx, …)
	std::string ext_id_type;      // optional; "" → no <EXT_ID/> emitted
	std::string ext_id_accession; // required when ext_id_type non-empty
};

inline std::string MakeReceiptXML(const std::vector<ReceiptObjectFixture> &objects, bool success = true,
                                  const std::string &error = "") {
	std::string out = "<?xml version=\"1.0\"?><RECEIPT receiptDate=\"2026-05-03T12:00:00Z\" "
	                  "submissionFile=\"mock\" success=\"";
	out += (success ? "true" : "false");
	out += "\">";
	for (const auto &o : objects) {
		out += "<" + o.element + " accession=\"" + o.accession + "\" alias=\"" + o.alias + "\" status=\"PRIVATE\"";
		if (o.ext_id_type.empty()) {
			out += "/>";
		} else {
			out += "><EXT_ID accession=\"" + o.ext_id_accession + "\" type=\"" + o.ext_id_type + "\"/></" + o.element +
			       ">";
		}
	}
	out += "<SUBMISSION accession=\"ERA999\" alias=\"mock\"/><ACTIONS>ADD</ACTIONS>";
	if (!success && !error.empty()) {
		out += "<MESSAGES><ERROR>" + error + "</ERROR></MESSAGES>";
	}
	out += "</RECEIPT>";
	return out;
}

} // namespace miint_test
