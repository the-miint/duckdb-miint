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

Invoked by run_tests.sh on a free port; URL is exported as ENA_WEBIN_MOCK_URL
and MIINT_ENA_PORTAL_URL_BASE = "$ENA_WEBIN_MOCK_URL/portal/api".
"""

from __future__ import annotations

import base64
import hashlib
import json
import os
import re
import sys
from http.server import BaseHTTPRequestHandler, HTTPServer

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


def _parse_xml_envelope(body: str) -> tuple[dict, str, str]:
    """Best-effort regex extraction of action / hold / object aliases from a
    `<WEBIN>...</WEBIN>` document. Mirrors the JSON envelope shape so the
    receipt builder doesn't care which form the request used. This mock
    doesn't run a real XML parser — the duckdb-miint envelope format is
    constrained and stable enough for a regex pass."""
    envelope: dict = {}
    action_name = "ADD"
    hold_until = ""
    actions = re.findall(r"<ACTION>\s*<(\w+)(?:\s+HoldUntilDate=\"([^\"]*)\")?\s*/>\s*</ACTION>", body)
    for tag, hold in actions:
        if tag == "HOLD":
            hold_until = hold or ""
        elif tag in ("ADD", "MODIFY", "CANCEL", "RELEASE", "VALIDATE"):
            action_name = tag
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
    return envelope, action_name, hold_until


def _build_receipt(envelope: dict, action_name: str = "ADD", hold_until: str = "") -> tuple[str, bool]:
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

    success = not any("FAIL" in alias for _, alias in objects)
    parts: list[str] = []
    parts.append('<?xml version="1.0" encoding="UTF-8"?>')
    parts.append(
        f'<RECEIPT receiptDate="2026-05-03T12:00:00.000Z" '
        f'submissionFile="mock" success="{"true" if success else "false"}">'
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
    if not success:
        parts.append('<MESSAGES>')
        for _, alias in objects:
            if "FAIL" in alias:
                parts.append(f'<ERROR>mock validation failure for alias ' f'{_xml_escape(alias)}</ERROR>')
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
        # Portal API alias collision lookup (Phase 8 Step 8a). Authenticated.
        if self.path.startswith("/portal/api/search"):
            if not self._check_auth():
                self.send_response(401)
                self.send_header("WWW-Authenticate", 'Basic realm="webin"')
                self.end_headers()
                return
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
        if "json" in ctype:
            try:
                envelope = json.loads(body.decode("utf-8"))
            except Exception:
                self.send_response(400)
                self.end_headers()
                return
        elif "xml" in ctype:
            envelope, action_name, hold_until = _parse_xml_envelope(body.decode("utf-8", "replace"))
        else:
            self.send_response(415)
            self.end_headers()
            return

        receipt, _ok = _build_receipt(envelope, action_name, hold_until)
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
