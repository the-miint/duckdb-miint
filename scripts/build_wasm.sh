#!/usr/bin/env bash
# Build the miint WASM extension locally and verify it has no unresolved symbols.
#
# Prerequisites:
#   - emscripten (emsdk) installed and activated
#   - VCPKG_TOOLCHAIN_PATH set (or vcpkg/ symlink in repo root)
#   - Rust with wasm32-unknown-emscripten target: rustup target add wasm32-unknown-emscripten
#   - Node.js (for import verification)
#
# Usage:
#   ./scripts/build_wasm.sh              # Build + verify wasm_eh (default)
#   ./scripts/build_wasm.sh --build-only # Build without verification
#   ./scripts/build_wasm.sh --all        # Build all three WASM variants
#   ./scripts/build_wasm.sh --clean      # Clean WASM build and ext/ artifacts, then rebuild
#
# The CI builds wasm_eh for shell.duckdb.org, so that's the default target.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

VERIFY=true
ALL_TARGETS=false
TARGET="wasm_eh"
CLEAN=false

for arg in "$@"; do
    case "$arg" in
        --build-only) VERIFY=false ;;
        --all)        ALL_TARGETS=true ;;
        --mvp)        TARGET="wasm_mvp" ;;
        --eh)         TARGET="wasm_eh" ;;
        --threads)    TARGET="wasm_threads" ;;
        --clean)      CLEAN=true ;;
        -h|--help)
            echo "Usage: $0 [--build-only] [--all] [--mvp|--eh|--threads] [--clean]"
            echo ""
            echo "  Default: build wasm_eh and verify no unresolved library imports."
            echo ""
            echo "  --build-only  Skip import verification after building"
            echo "  --all         Build all three WASM variants (mvp, eh, threads)"
            echo "  --mvp/--eh/--threads  Build a specific variant (default: eh)"
            echo "  --clean       Remove WASM build dir and ext/ artifacts before building"
            exit 0
            ;;
        *) echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

# Check prerequisites
if ! command -v emcc &>/dev/null; then
    echo "ERROR: emcc not found. Install and activate emsdk first."
    echo "  git clone https://github.com/emscripten-core/emsdk.git"
    echo "  cd emsdk && ./emsdk install latest && ./emsdk activate latest && source emsdk_env.sh"
    exit 1
fi

if [ -z "${VCPKG_TOOLCHAIN_PATH:-}" ]; then
    if [ -d "$REPO_DIR/vcpkg" ]; then
        export VCPKG_TOOLCHAIN_PATH="$REPO_DIR/vcpkg/scripts/buildsystems/vcpkg.cmake"
        echo "Auto-detected VCPKG_TOOLCHAIN_PATH=$VCPKG_TOOLCHAIN_PATH"
    else
        echo "ERROR: VCPKG_TOOLCHAIN_PATH not set and no vcpkg/ directory found."
        exit 1
    fi
fi

# WASM builds use the wasm32-emscripten vcpkg triplet (set by CI and the DuckDB
# Makefile via VCPKG_EMSDK_FLAGS). Export it here so vcpkg can find the right
# packages if the Makefile doesn't set it from the toolchain file.
export VCPKG_TARGET_TRIPLET="${VCPKG_TARGET_TRIPLET:-wasm32-emscripten}"

if ! rustup target list --installed 2>/dev/null | grep -q wasm32-unknown-emscripten; then
    echo "Installing Rust wasm32-unknown-emscripten target..."
    rustup target add wasm32-unknown-emscripten
fi

EMCC_VERSION=$(emcc --version 2>&1 | head -1 | grep -oP '\d+\.\d+\.\d+' || echo "unknown")
echo "=== WASM build (emcc $EMCC_VERSION) ==="
echo "    CI uses emsdk 3.1.71"
echo ""

cd "$REPO_DIR"

if $CLEAN; then
    echo "Cleaning WASM build artifacts..."
    for t in wasm_mvp wasm_eh wasm_threads; do
        rm -rf "build/$t"
    done
    # External libs build in-source; clean their objects so they rebuild with correct flags
    make -C ext/minimap2 clean 2>/dev/null || true
    make -C ext/WFA2-lib clean 2>/dev/null || true
    make -C ext/htslib-1.22.1 clean 2>/dev/null || true
    make -C ext/mafft/core clean 2>/dev/null || true
    echo ""
fi

if $ALL_TARGETS; then
    TARGETS="wasm_mvp wasm_eh wasm_threads"
else
    TARGETS="$TARGET"
fi

BUILD_OK=true
for t in $TARGETS; do
    echo "--- Building $t ---"
    mkdir -p "build"
    if ! make "$t" 2>&1 | tee "build/${t}_build.log"; then
        echo ""
        echo "ERROR: Build failed for $t. Check build/${t}_build.log"
        BUILD_OK=false
        continue
    fi
    WASM_FILE=$(find "build/$t" -name "miint.duckdb_extension.wasm" 2>/dev/null | head -1)
    if [ -z "$WASM_FILE" ]; then
        echo "ERROR: miint.duckdb_extension.wasm not found after building $t"
        BUILD_OK=false
        continue
    fi
    SIZE=$(du -h "$WASM_FILE" | cut -f1)
    echo "SUCCESS: $WASM_FILE ($SIZE)"
    echo ""
done

if ! $BUILD_OK; then
    echo "FAILED: One or more WASM builds failed."
    exit 1
fi

if $VERIFY; then
    echo "=== Verifying WASM extension imports ==="
    echo ""

    if ! command -v node &>/dev/null; then
        echo "WARNING: Node.js not found, skipping import verification."
        echo "Install Node.js to enable WASM import checking."
        exit 0
    fi

    VERIFY_OK=true
    for t in $TARGETS; do
        WASM_FILE=$(find "build/$t" -name "miint.duckdb_extension.wasm" 2>/dev/null | head -1)
        if [ -z "$WASM_FILE" ]; then
            continue
        fi

        echo "--- Checking $t: $WASM_FILE ---"

        # Use Node.js to compile the WASM module and inspect imports/exports.
        #
        # In PIC WASM side modules, library symbols appear as both:
        #   - GOT.func.X / GOT.mem.X imports (PIC mechanism, resolved at load time)
        #   - Direct exports of X (the actual symbol definition)
        #
        # A GOT import WITH a matching export is correct (self-resolved PIC).
        # A GOT import WITHOUT a matching export means the symbol wasn't linked.
        # A direct env import of a library symbol means it wasn't linked at all.
        RESULT=$(node -e "
            const fs = require('fs');
            const buf = fs.readFileSync('$WASM_FILE');
            WebAssembly.compile(buf).then(mod => {
                const imports = WebAssembly.Module.imports(mod);
                const exports = WebAssembly.Module.exports(mod);
                const exportNames = new Set(exports.map(e => e.name));

                // Patterns for symbols that must be statically linked into the extension
                const libPatterns = [
                    /^gz(close|open|d?open|read|write|gets|puts|eof|rewind|seek|tell|flush)/,
                    /^(deflate|inflate|compress|uncompress|adler32|crc32)/,
                    /^hts_/,  /^sam_/,  /^bam_/,  /^bcf_/,  /^bgzf_/,
                    /^mm_/,  /^mm2_/,
                    /^wavefront/,
                    /^XML_/,
                    /^rype_/,
                ];
                const isLib = name => libPatterns.some(p => p.test(name));

                // Check 1: Direct (env) imports of library symbols = not linked at all
                const badDirect = imports.filter(i => i.module === 'env' && isLib(i.name));

                // Check 2: GOT imports without matching exports = linked but not exported
                const gotImports = imports.filter(i =>
                    (i.module === 'GOT.func' || i.module === 'GOT.mem') && isLib(i.name)
                );
                const badGot = gotImports.filter(i => !exportNames.has(i.name));

                // Check 3: Extension init function must be exported
                const hasInit = exportNames.has('miint_duckdb_cpp_init');

                const problems = [];
                if (!hasInit) problems.push('Missing miint_duckdb_cpp_init export');
                badDirect.forEach(i => problems.push('Unlinked: env.' + i.name));
                badGot.forEach(i => problems.push('GOT without export: ' + i.module + '.' + i.name));

                if (problems.length > 0) {
                    console.log('FAIL');
                    problems.forEach(p => console.log('  ' + p));
                } else {
                    const resolvedGot = gotImports.length - badGot.length;
                    console.log('OK');
                    console.log('  Imports: ' + imports.length + ' total, ' + gotImports.length + ' library GOT (all self-resolved)');
                    console.log('  Exports: ' + exports.length + ' total, miint_duckdb_cpp_init present');
                }
            }).catch(err => {
                console.log('COMPILE_ERROR ' + err.message);
            });
        " 2>&1)

        STATUS=$(echo "$RESULT" | head -1)
        case "$STATUS" in
            OK)
                echo "  PASS"
                echo "$RESULT" | tail -n+2
                ;;
            FAIL)
                echo "  FAIL:"
                echo "$RESULT" | tail -n+2
                VERIFY_OK=false
                ;;
            COMPILE_ERROR*)
                echo "  FAIL: WebAssembly.compile error: $(echo "$STATUS" | cut -d' ' -f2-)"
                VERIFY_OK=false
                ;;
            *)
                echo "  FAIL: Unexpected output: $RESULT"
                VERIFY_OK=false
                ;;
        esac
        echo ""
    done

    if ! $VERIFY_OK; then
        echo "FAILED: Import verification found issues."
        exit 1
    fi

    echo "=== All checks passed ==="
fi

echo "Done."
