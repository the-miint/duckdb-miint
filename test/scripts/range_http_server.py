#!/usr/bin/env python3
"""Serve a directory over HTTP, with byte ranges, for the HTTP reader tests.

`python3 -m http.server` ignores Range headers and answers every GET with the
whole file. httpfs then falls back to downloading the full file and logs a
WARNING, which DuckDB 2.0's shell prints with the query result, so the BAM
checks in test/shell/read_https.sh fail; and the range-request path that real
servers (S3, ENA, NCBI) put httpfs on is never exercised. This server answers
`Range: bytes=FIRST-LAST`, the one form httpfs sends, with 206 Partial
Content as those servers do. Anything else is served as http.server serves it.

Usage: range_http_server.py PORT DIRECTORY   (listens on 127.0.0.1 only)
"""

import functools
import http.server
import io
import os
import re
import sys
from http import HTTPStatus

RANGE = re.compile(r"bytes=(\d+)-(\d+)")


class RangeRequestHandler(http.server.SimpleHTTPRequestHandler):
    def send_head(self):
        path = self.translate_path(self.path)
        if not os.path.isfile(path):
            # Directory listings, redirects and 404s are unchanged.
            return super().send_head()
        with open(path, "rb") as f:
            data = f.read()
        size = len(data)
        m = RANGE.fullmatch(self.headers.get("Range", ""))
        if not m or int(m[1]) > int(m[2]):
            # A server may ignore a Range it does not handle (RFC 9110 14.2).
            body = data
            self.send_response(HTTPStatus.OK)
        elif int(m[1]) >= size:
            # What real servers answer for a range past the end.
            self.send_response(HTTPStatus.REQUESTED_RANGE_NOT_SATISFIABLE)
            self.send_header("Content-Range", f"bytes */{size}")
            self.send_header("Content-Length", "0")
            self.end_headers()
            return None
        else:
            first, last = int(m[1]), min(int(m[2]), size - 1)
            body = data[first : last + 1]
            self.send_response(HTTPStatus.PARTIAL_CONTENT)
            self.send_header("Content-Range", f"bytes {first}-{last}/{size}")
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Content-Type", self.guess_type(path))
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Last-Modified", self.date_time_string(os.stat(path).st_mtime))
        self.end_headers()
        return io.BytesIO(body)


def main():
    port, directory = int(sys.argv[1]), sys.argv[2]
    handler = functools.partial(RangeRequestHandler, directory=directory)
    http.server.ThreadingHTTPServer(("127.0.0.1", port), handler).serve_forever()


if __name__ == "__main__":
    main()
