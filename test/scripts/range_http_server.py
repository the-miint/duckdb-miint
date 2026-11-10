#!/usr/bin/env python3
"""Serve a directory over HTTP, with byte ranges, for the HTTP reader tests.

`python3 -m http.server` ignores Range headers and answers every GET with the
whole file. httpfs then falls back to downloading the full file and logs a
WARNING, which DuckDB 2.0's shell prints with the query result, so the BAM
checks in test/shell/read_https.sh fail; and the range-request path that real
servers (S3, ENA, NCBI) put httpfs on is never exercised. This server answers
`Range: bytes=FIRST-LAST`, the one form httpfs sends, with 206 Partial
Content as those servers do. Anything else is served as http.server serves it.

Two path prefixes make a GET of a file fail as downloads do:
- /cut/N/<path>: the first N GETs of <path> send three quarters of the body, then
  reset the connection, as a dropped download does.
- /slow/<path>: sends the first 600 bytes of the body, then a byte every 100 ms.

Usage: range_http_server.py PORT DIRECTORY   (listens on 127.0.0.1 only)
"""

import collections
import functools
import http.server
import io
import os
import re
import socket
import struct
import sys
import threading
import time
from http import HTTPStatus

RANGE = re.compile(r"bytes=(\d+)-(\d+)")
CUT = re.compile(r"/cut/(\d+)(/.*)")
SLOW = "/slow"
SLOW_FIRST_BYTES = 600
SLOW_INTERVAL_SECONDS = 0.1

_cut_gets = collections.Counter()
_cut_lock = threading.Lock()


class RangeRequestHandler(http.server.SimpleHTTPRequestHandler):
    _was_reset = False

    def do_GET(self):
        cut = CUT.fullmatch(self.path)
        if cut:
            with _cut_lock:
                _cut_gets[self.path] += 1
                due = _cut_gets[self.path] <= int(cut[1])
            self.path = cut[2]
            if due:
                self._send_cut()
                return
        elif self.path.startswith(SLOW + "/"):
            self.path = self.path[len(SLOW) :]
            self._send_slow()
            return
        super().do_GET()

    def _send_cut(self):
        body = self.send_head()
        if body is None:
            return
        data = body.read()
        self.wfile.write(data[: len(data) * 3 // 4])
        self.wfile.flush()
        # SO_LINGER 0 makes close() send RST: the client sees a transport error, not a clean EOF.
        self.connection.setsockopt(socket.SOL_SOCKET, socket.SO_LINGER, struct.pack("ii", 1, 0))
        self.connection.close()
        self.close_connection = True
        self._was_reset = True

    def _send_slow(self):
        body = self.send_head()
        if body is None:
            return
        data = body.read()
        try:
            self.wfile.write(data[:SLOW_FIRST_BYTES])
            self.wfile.flush()
            for i in range(SLOW_FIRST_BYTES, len(data)):
                time.sleep(SLOW_INTERVAL_SECONDS)
                self.wfile.write(data[i : i + 1])
                self.wfile.flush()
        except (BrokenPipeError, ConnectionResetError):
            # The client stopped reading, which is what the test waits for.
            self.close_connection = True

    def finish(self):
        try:
            super().finish()
        except OSError:
            if not self._was_reset:
                raise

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
