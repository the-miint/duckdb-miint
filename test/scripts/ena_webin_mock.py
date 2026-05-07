#!/usr/bin/env python3
"""Tiny in-process HTTP mock of the ENA Webin V2 /submit endpoint.

Used by SQL tests in the duckdb-miint extension that need to exercise the
INSERT INTO ena.{...} pipeline end-to-end without hitting wwwdev.ebi.ac.uk.

Behavior:
  - Requires HTTP Basic auth with user starting with "Webin-" (any password).
    Returns 401 otherwise.
  - On POST to /submit (or /submit/queue), inspects the JSON envelope and
    returns an XML receipt with one <PROJECT>/<SAMPLE>/... element per object,
    each carrying a deterministic accession derived from its alias (so tests
    can assert specific values without parsing).
  - Accession scheme: project alias "p1" -> PRJEB1000 + ERP1000;
    sample alias "s1" -> ERS2000 + SAMEA2000. The trailing number is a hash
    of the alias mod 10000, offset by object-type base.
  - Every receipt declares success="true" unless the alias contains the
    literal substring "FAIL", in which case success="false" and an error
    message is emitted.
  - On GET to /portal/api/search (Phase 8 alias collision check), parses the
    `query=<kind>_alias IN ("a1","a2",...)` filter and returns a TSV with
    one row per alias starting with the literal prefix "EXISTS_". All other
    aliases are reported as not present.
  - On GET to /browser/api/xml/<accession> (Phase 8 checklist registry):
      * "TST000001" -> strict fixture (mandatory + units + CV) for the
        Phase 8 Step 8b validation tests.
      * Anything else -> permissive fixture (every label the existing
        sample-insert tests use, all marked optional) so older mock tests
        keep round-tripping without rewriting their attribute MAPs.

Invoked by run_tests.sh on a free port; URL is exported as ENA_WEBIN_MOCK_URL,
MIINT_ENA_PORTAL_URL_BASE = "$ENA_WEBIN_MOCK_URL/portal/api", and
MIINT_ENA_CHECKLIST_URL_BASE = "$ENA_WEBIN_MOCK_URL/browser/api/xml".
"""

from __future__ import annotations

import base64
import hashlib
import json
import os
import re
import sys
from http.server import BaseHTTPRequestHandler, HTTPServer

_CHECKLIST_TST000001 = """<?xml version="1.0" encoding="UTF-8"?>
<CHECKLIST_SET>
<CHECKLIST accession="TST000001" checklistType="Sample">
  <IDENTIFIERS><PRIMARY_ID>TST000001</PRIMARY_ID></IDENTIFIERS>
  <DESCRIPTOR>
    <LABEL>miint Phase 8 test checklist</LABEL>
    <NAME>miint_phase8_test</NAME>
    <DESCRIPTION>Synthetic checklist used by ena_checklist_validation_mock.test.</DESCRIPTION>
    <FIELD_GROUP restrictionType="Any number or none of the fields">
      <FIELD>
        <LABEL>project name</LABEL>
        <NAME>project_name</NAME>
        <FIELD_TYPE><TEXT_FIELD/></FIELD_TYPE>
        <MANDATORY>mandatory</MANDATORY>
      </FIELD>
      <FIELD>
        <LABEL>collection date</LABEL>
        <NAME>collection_date</NAME>
        <FIELD_TYPE><TEXT_FIELD/></FIELD_TYPE>
        <MANDATORY>mandatory</MANDATORY>
      </FIELD>
      <FIELD>
        <LABEL>geographic location (latitude)</LABEL>
        <NAME>geographic_location_latitude</NAME>
        <FIELD_TYPE><TEXT_FIELD/></FIELD_TYPE>
        <UNITS><UNIT>DD</UNIT></UNITS>
        <MANDATORY>mandatory</MANDATORY>
      </FIELD>
      <FIELD>
        <LABEL>country</LABEL>
        <NAME>country</NAME>
        <FIELD_TYPE>
          <TEXT_CHOICE_FIELD>
            <TEXT_VALUE><VALUE>USA</VALUE></TEXT_VALUE>
            <TEXT_VALUE><VALUE>Canada</VALUE></TEXT_VALUE>
            <TEXT_VALUE><VALUE>Mexico</VALUE></TEXT_VALUE>
          </TEXT_CHOICE_FIELD>
        </FIELD_TYPE>
        <MANDATORY>mandatory</MANDATORY>
      </FIELD>
      <FIELD>
        <LABEL>ploidy</LABEL>
        <NAME>ploidy</NAME>
        <FIELD_TYPE><TEXT_FIELD/></FIELD_TYPE>
        <MANDATORY>optional</MANDATORY>
      </FIELD>
    </FIELD_GROUP>
  </DESCRIPTOR>
</CHECKLIST>
</CHECKLIST_SET>
"""

_PERMISSIVE_LABELS = (
    "project name",
    "collection date",
    "geographic location (country and/or sea)",
    "geographic location (latitude)",
    "geographic location (longitude)",
    "broad-scale environmental context",
    "local environmental context",
    "environmental medium",
)
_CHECKLIST_PERMISSIVE = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<CHECKLIST_SET>\n'
    '<CHECKLIST accession="ERC000000-mock" checklistType="Sample">\n'
    '  <IDENTIFIERS><PRIMARY_ID>ERC000000-mock</PRIMARY_ID></IDENTIFIERS>\n'
    '  <DESCRIPTOR>\n'
    '    <LABEL>miint mock permissive checklist</LABEL>\n'
    '    <NAME>miint_mock_permissive</NAME>\n'
    '    <FIELD_GROUP restrictionType="Any number or none of the fields">\n'
    + "".join(
        f"      <FIELD><LABEL>{label}</LABEL><NAME>{label.replace(' ', '_')}</NAME>"
        f"<FIELD_TYPE><TEXT_FIELD/></FIELD_TYPE><MANDATORY>optional</MANDATORY></FIELD>\n"
        for label in _PERMISSIVE_LABELS
    )
    + '    </FIELD_GROUP>\n'
    '  </DESCRIPTOR>\n'
    '</CHECKLIST>\n'
    '</CHECKLIST_SET>\n'
)


ACCESSION_BASE = {
    # (base_offset, primary_prefix_on_element, ext_id_prefix)
    # Primary prefix is what ENA returns as the <ELEMENT accession="..."/>;
    # ext_id_prefix is the EXT_ID accession (study ERP for projects,
    # BioSample SAMEA for samples). See ena-research-webin-v2-deep.md §1.
    "PROJECT": (1000, "PRJEB", "ERP"),
    "SAMPLE": (2000, "ERS", "SAMEA"),
    "EXPERIMENT": (3000, "ERX", "ERX"),
    "RUN": (4000, "ERR", "ERR"),
    "ANALYSIS": (5000, "ERZ", "ERZ"),
}


def _accession_for(alias: str, kind: str) -> tuple[str, str]:
    base, primary_prefix, ext_prefix = ACCESSION_BASE[kind]
    digest = int(hashlib.sha256(alias.encode("utf-8")).hexdigest(), 16)
    n = base + (digest % 10000)
    return f"{primary_prefix}{n}", f"{ext_prefix}{n}"


def _xml_escape(s: str) -> str:
    return (
        s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;").replace("'", "&apos;")
    )


def _parse_xml_envelope(body: str) -> tuple[dict, str, str, str]:
    """Best-effort regex extraction of action / hold / target / object aliases
    from a `<WEBIN>...</WEBIN>` document. Mirrors the JSON envelope shape so
    the receipt builder doesn't care which form the request used. This mock
    doesn't run a real XML parser — the duckdb-miint envelope format is
    constrained and stable enough for a regex pass.

    Returns (envelope, action_name, hold_until, target). `target` is set when
    the action element carries `target="..."` (CANCEL / RELEASE / targeted HOLD)."""
    envelope: dict = {}
    action_name = "ADD"
    hold_until = ""
    target = ""
    # Action element can be:
    #   <CANCEL/>                                       (untargeted)
    #   <CANCEL target="ERS123"/>                       (targeted CANCEL/RELEASE)
    #   <HOLD HoldUntilDate="2027-01-01"/>              (body-pattern HOLD)
    #   <HOLD target="ERS123" HoldUntilDate="2027-..."/> (targeted HOLD)
    # Capture target= and HoldUntilDate= attributes in any order.
    actions = re.findall(
        r'<ACTION>\s*<(\w+)((?:\s+\w+="[^"]*")*)\s*/>\s*</ACTION>',
        body,
    )
    for tag, attrs in actions:
        attr_dict = dict(re.findall(r'(\w+)="([^"]*)"', attrs))
        hold_attr = attr_dict.get("HoldUntilDate", "")
        target_attr = attr_dict.get("target", "")
        if tag == "HOLD":
            if hold_attr:
                hold_until = hold_attr
            if target_attr:
                target = target_attr
                action_name = "HOLD"
        elif tag in ("ADD", "MODIFY", "CANCEL", "RELEASE", "VALIDATE"):
            action_name = tag
            if target_attr:
                target = target_attr
    for kind_plural, set_tag, item_tag in (
        ("experiments", "EXPERIMENT_SET", "EXPERIMENT"),
        ("runs", "RUN_SET", "RUN"),
        ("analyses", "ANALYSIS_SET", "ANALYSIS"),
    ):
        block = re.search(rf"<{set_tag}>(.*?)</{set_tag}>", body, re.DOTALL)
        if not block:
            continue
        items = []
        for m in re.finditer(rf"<{item_tag}\s+alias=\"([^\"]+)\"", block.group(1)):
            items.append({"alias": m.group(1)})
        envelope[kind_plural] = items
    return envelope, action_name, hold_until, target


def _build_receipt(
    envelope: dict, action_name: str = "ADD", hold_until: str = "", target: str = ""
) -> tuple[str, bool]:
    if not action_name and not hold_until:
        actions = envelope.get("submission", {}).get("actions", [])
        action_name = "ADD"
        for a in actions:
            t = a.get("type", "")
            if t == "HOLD":
                hold_until = a.get("holdUntilDate", "")
            elif t in ("ADD", "MODIFY", "CANCEL", "RELEASE", "VALIDATE"):
                action_name = t

    objects: list[tuple[str, str]] = []
    for kind_plural, kind_singular in (
        ("projects", "PROJECT"),
        ("samples", "SAMPLE"),
        ("experiments", "EXPERIMENT"),
        ("runs", "RUN"),
        ("analyses", "ANALYSIS"),
    ):
        for obj in envelope.get(kind_plural, []) or []:
            objects.append((kind_singular, obj.get("alias", "")))

    # Target-based failure trigger for lifecycle ops (which have no body
    # objects to encode FAIL into): a target containing "FAIL" produces a
    # failure receipt. Mirrors the alias-based trigger for ADD-style
    # submissions.
    target_fail = "FAIL" in target if target else False
    is_lifecycle = bool(target)
    body_fail = any("FAIL" in alias for _, alias in objects)
    # `success` here is the boolean returned to callers as the round-trip
    # verdict (and matches `outcome.success` in the C++ lifecycle code:
    # success iff no <ERROR> elements were emitted). The XML attribute on
    # <RECEIPT> follows real ENA semantics — see comments below.
    success = not target_fail and not body_fail
    # On wwwdev cross-submission lifecycle responses (CANCEL / RELEASE /
    # HOLD on an existing accession), ENA sets RECEIPT @success="false"
    # even when the action took effect — the attribute means "did this
    # submission produce new accessions", which is always false for
    # lifecycle ops. The actual outcome is in <INFO>/<ERROR> children.
    # Match that behaviour here so SubmitLifecycle's
    # `success = errors.empty()` interpretation gets exercised by tests.
    receipt_success_attr = "false" if is_lifecycle else ("true" if success else "false")
    parts: list[str] = []
    parts.append('<?xml version="1.0" encoding="UTF-8"?>')
    parts.append(
        f'<RECEIPT receiptDate="2026-05-03T12:00:00.000Z" ' f'submissionFile="mock" success="{receipt_success_attr}">'
    )
    for kind_singular, alias in objects:
        if not alias:
            continue
        primary, ext = _accession_for(alias, kind_singular)
        parts.append(
            f'<{kind_singular} accession="{_xml_escape(primary)}" '
            f'alias="{_xml_escape(alias)}" status="PRIVATE"'
            f'{" holdUntilDate=" + chr(34) + _xml_escape(hold_until) + chr(34) if hold_until else ""}>'
        )
        ext_type = (
            "study"
            if kind_singular == "PROJECT"
            else ("biosample" if kind_singular == "SAMPLE" else kind_singular.lower())
        )
        parts.append(f'<EXT_ID accession="{_xml_escape(ext)}" type="{ext_type}"/>')
        parts.append(f'</{kind_singular}>')
    parts.append('<SUBMISSION accession="ERA1234567" alias="mock"/>')
    parts.append(f'<ACTIONS>{action_name}</ACTIONS>')
    # Failure-path messages (always <ERROR>, mirrors ENA).
    if not success:
        parts.append('<MESSAGES>')
        for _, alias in objects:
            if "FAIL" in alias:
                parts.append(f'<ERROR>mock validation failure for alias ' f'{_xml_escape(alias)}</ERROR>')
        if target_fail:
            parts.append(f'<ERROR>mock validation failure for target {_xml_escape(target)}</ERROR>')
        parts.append('</MESSAGES>')
    elif is_lifecycle:
        # Success-path lifecycle: emit <INFO> mirroring ENA's wwwdev format
        # ("...is set to cancelled/released/private status."). SubmitLifecycle
        # treats this as success because errors.empty() is true.
        status_word = {"CANCEL": "cancelled", "RELEASE": "public", "HOLD": "private"}.get(action_name, "updated")
        parts.append('<MESSAGES>')
        parts.append(f'<INFO>accession "{_xml_escape(target)}" is set to {status_word} status.</INFO>')
        parts.append('</MESSAGES>')
    parts.append('</RECEIPT>')
    return "".join(parts), success


class Handler(BaseHTTPRequestHandler):
    def log_message(self, format, *args):  # silence stderr
        return

    def _check_auth(self) -> bool:
        header = self.headers.get("Authorization", "")
        if not header.startswith("Basic "):
            return False
        try:
            user_pw = base64.b64decode(header[6:]).decode("utf-8", "replace")
        except Exception:
            return False
        return user_pw.startswith("Webin-") and ":" in user_pw

    def do_GET(self):
        # Checklist XML (Phase 8 Step 8b). Anonymous GET — public reference
        # data on the real EBI browser API.
        if self.path.startswith("/browser/api/xml/"):
            xml = self._handle_checklist_fetch(self.path)
            payload = xml.encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/xml; charset=utf-8")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
            return
        # Portal API alias collision lookup. Anonymous — the real ENA
        # portal at /ena/portal/api/search has no authenticated mode (Basic
        # auth → HTTP 500), and ena-upload-cli's check_remote.py likewise
        # queries it without credentials.
        if self.path.startswith("/portal/api/search"):
            tsv = self._handle_portal_search(self.path)
            if tsv is None:
                self.send_response(400)
                self.end_headers()
                return
            payload = tsv.encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "text/tab-separated-values; charset=utf-8")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
            return
        # Default: health check.
        self.send_response(200)
        self.send_header("Content-Type", "text/plain")
        self.end_headers()
        self.wfile.write(b"ena-webin-mock ok\n")

    def _handle_checklist_fetch(self, path: str) -> str:
        accession = path.rsplit("/", 1)[-1]
        if accession == "TST000001":
            return _CHECKLIST_TST000001
        # Default permissive fixture: a synthetic checklist whose fields are
        # the union of labels used by the legacy ena_samples_insert_mock.test
        # and ena_full_pipeline.test, all marked optional with no units / no
        # CV. Pre-Phase-8 mock tests don't carry exhaustive attribute sets;
        # this lets the validator pass through those tests unchanged while
        # still being strict for new tests that opt into TST000001.
        return _CHECKLIST_PERMISSIVE

    def _handle_portal_search(self, path: str) -> str | None:
        """Parse `?query=<kind>_alias IN ("a1","a2")&fields=<kind>_alias&...`
        out of the URL and return a TSV body. Returns None on a malformed
        request so the caller can emit HTTP 400."""
        from urllib.parse import urlparse, parse_qs, unquote

        qs = parse_qs(urlparse(path).query)
        fields = qs.get("fields", [""])[0]
        query = unquote(qs.get("query", [""])[0])
        if not fields or not query:
            return None
        # Accept both "<field> IN (\"a\",\"b\")" and the single
        # "<field> = \"a\"" forms. Extract bare aliases.
        aliases: list[str] = []
        if " IN " in query:
            inner = query.split(" IN ", 1)[1].strip()
            if inner.startswith("(") and inner.endswith(")"):
                inner = inner[1:-1]
            for part in inner.split(","):
                part = part.strip()
                if part.startswith('"') and part.endswith('"'):
                    aliases.append(part[1:-1])
        elif "=" in query:
            _, rhs = query.split("=", 1)
            rhs = rhs.strip()
            if rhs.startswith('"') and rhs.endswith('"'):
                aliases.append(rhs[1:-1])
        # Only aliases that start with the literal prefix "EXISTS_" are
        # treated as already-present in the submission account. SQL tests
        # use this convention to drive the collision path without a stateful
        # mock. (Production users won't see this behaviour because the real
        # ENA portal API replies with their actual submission history.)
        present = [a for a in aliases if a.startswith("EXISTS_")]
        lines = [fields] + present
        return "\n".join(lines) + "\n"

    def do_POST(self):
        if not self._check_auth():
            self.send_response(401)
            self.send_header("WWW-Authenticate", 'Basic realm="webin"')
            self.end_headers()
            return
        if self.path not in ("/submit", "/submit/queue"):
            self.send_response(404)
            self.end_headers()
            return

        length = int(self.headers.get("Content-Length", "0") or "0")
        body = self.rfile.read(length) if length else b""
        ctype = self.headers.get("Content-Type", "")
        action_name = ""
        hold_until = ""
        target = ""
        if "json" in ctype:
            try:
                envelope = json.loads(body.decode("utf-8"))
            except Exception:
                self.send_response(400)
                self.end_headers()
                return
        elif "xml" in ctype:
            envelope, action_name, hold_until, target = _parse_xml_envelope(body.decode("utf-8", "replace"))
        else:
            self.send_response(415)
            self.end_headers()
            return

        receipt, _ok = _build_receipt(envelope, action_name, hold_until, target)
        payload = receipt.encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/xml; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)


def main() -> int:
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    server = HTTPServer(("127.0.0.1", port), Handler)
    bound_port = server.server_address[1]
    pid_file = os.environ.get("ENA_WEBIN_MOCK_PID_FILE")
    if pid_file:
        with open(pid_file, "w") as f:
            f.write(str(os.getpid()))
    print(f"ENA_WEBIN_MOCK_URL=http://127.0.0.1:{bound_port}", flush=True)
    server.serve_forever()
    return 0


if __name__ == "__main__":
    sys.exit(main())
