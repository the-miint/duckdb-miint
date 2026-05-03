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

Invoked by run_tests.sh on a free port; URL is exported as ENA_WEBIN_MOCK_URL.
"""

from __future__ import annotations

import base64
import hashlib
import json
import os
import sys
from http.server import BaseHTTPRequestHandler, HTTPServer

ACCESSION_BASE = {
    "PROJECT": (1000, "PRJEB", "ERP"),
    "SAMPLE": (2000, "SAMEA", "ERS"),
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


def _build_receipt(envelope: dict) -> tuple[str, bool]:
    actions = envelope.get("submission", {}).get("actions", [])
    action_name = "ADD"
    hold_until = ""
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
        # health check
        self.send_response(200)
        self.send_header("Content-Type", "text/plain")
        self.end_headers()
        self.wfile.write(b"ena-webin-mock ok\n")

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
        if "json" not in ctype:
            self.send_response(415)
            self.end_headers()
            return
        try:
            envelope = json.loads(body.decode("utf-8"))
        except Exception:
            self.send_response(400)
            self.end_headers()
            return

        receipt, _ok = _build_receipt(envelope)
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
