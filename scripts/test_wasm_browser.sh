#!/usr/bin/env bash
# Browser check for the miint WASM extension: load it into duckdb-wasm the way
# the playground does, in headless Chrome, and run scripts/wasm_browser/steps.json.
#
# This, not test_wasm_load.sh, is the check for Rust code paths: the load test
# links a main module and JS glue built from the duckdb/ submodule, while this
# runs the duckdb-wasm release users get (docs/wasm-testing.md).
#
# Prerequisites:
#   - Node.js 22.12 or newer and npm (puppeteer-core is installed into scripts/wasm_browser/)
#   - Google Chrome or Chromium; set CHROME to its binary if it is not found
#   - Network access: duckdb-wasm is fetched from cdn.jsdelivr.net, as in the playground
#   - An extension build: ./scripts/build_wasm.sh
#
# Usage:
#   ./scripts/test_wasm_browser.sh [extension.wasm]   # default: build/wasm_eh's, as test_wasm_load.sh
#   ./scripts/test_wasm_browser.sh --native [duckdb]  # same steps natively (default: build/release/duckdb)
#
# --native is how expected values are produced: a new step's `expect` is what the
# native build returns.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DRIVER_DIR="$SCRIPT_DIR/wasm_browser"

if ! command -v node &>/dev/null; then
    echo "ERROR: node not found. Install Node.js."
    exit 1
fi

if [ "${1:-}" = "--native" ]; then
    DUCKDB="${2:-$REPO_DIR/build/release/duckdb}"
    [ -x "$DUCKDB" ] || { echo "ERROR: $DUCKDB not found. Build it with: bash build.sh"; exit 1; }
    exec node "$DRIVER_DIR/driver.mjs" --native "$DUCKDB"
fi

WASM_FILE="${1:-$REPO_DIR/build/wasm_eh/extension/miint/miint.duckdb_extension.wasm}"
if [ ! -f "$WASM_FILE" ]; then
    echo "ERROR: no extension .wasm at $WASM_FILE. Build it with: ./scripts/build_wasm.sh"
    exit 1
fi

# puppeteer-core 25 needs Node 22.12 or newer.
NODE_VERSION="$(node --version)"
IFS=. read -r NODE_MAJOR NODE_MINOR _ <<<"${NODE_VERSION#v}"
if [ "$NODE_MAJOR" -lt 22 ] || { [ "$NODE_MAJOR" -eq 22 ] && [ "$NODE_MINOR" -lt 12 ]; }; then
    echo "ERROR: Node $NODE_VERSION is too old; the browser check needs Node 22.12 or newer."
    exit 1
fi
if ! command -v npm &>/dev/null; then
    echo "ERROR: npm not found. It comes with Node.js."
    exit 1
fi

if [ -z "${CHROME:-}" ]; then
    for c in "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
             google-chrome google-chrome-stable chromium chromium-browser; do
        if [ -x "$c" ] || command -v "$c" &>/dev/null; then
            CHROME="$(command -v "$c" || echo "$c")"
            break
        fi
    done
fi
if [ -z "${CHROME:-}" ]; then
    echo "ERROR: Chrome or Chromium not found. Set CHROME to its binary."
    exit 1
fi

# Install, or reinstall when node_modules no longer matches package-lock.json.
# npm records what it installed in node_modules/.package-lock.json: the lockfile's
# packages without the root entry.
if ! node -e '
    const fs = require("fs");
    const lock = JSON.parse(fs.readFileSync(process.argv[1] + "/package-lock.json", "utf8")).packages;
    const installed = JSON.parse(fs.readFileSync(process.argv[1] + "/node_modules/.package-lock.json", "utf8")).packages;
    delete lock[""];
    process.exit(JSON.stringify(lock) === JSON.stringify(installed) ? 0 : 1);
' "$DRIVER_DIR" 2>/dev/null; then
    echo "Installing puppeteer-core into $DRIVER_DIR ..."
    npm ci --prefix "$DRIVER_DIR" --ignore-scripts --no-audit --no-fund
fi

# The DuckDB version the extension targets: the one scripts/check-duckdb-version.sh treats as
# canonical, which build_wasm.sh builds by default.
DUCKDB_VERSION="$("$REPO_DIR/scripts/check-duckdb-version.sh" --print)"

echo "=== Browser check: $WASM_FILE ==="
exec node "$DRIVER_DIR/driver.mjs" --wasm "$WASM_FILE" --chrome "$CHROME" --duckdb-version "$DUCKDB_VERSION"
