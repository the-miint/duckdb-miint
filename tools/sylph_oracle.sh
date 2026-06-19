#!/usr/bin/env bash
# Regenerate data/sylph/tiny.syldb and data/sylph/expected_profile.tsv from
# the embedded sylph fork's CLI.
#
# Run from the repo root:
#     bash tools/sylph_oracle.sh
#
# Mirrors the sortmerna oracle pattern (see data/sortmerna/README.md).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

# Locate the sylph CLI. After Phase 1 wires CMake, the ExternalProject build
# output is the canonical path. As a developer convenience for Phase 0/1, we
# also accept a cargo build inside the submodule itself (faster iteration than
# rebuilding the whole extension).
EMBEDDED_CLI="$REPO_ROOT/build/release/extension/miint/sylph_build-prefix/src/sylph_build/target/release/sylph"
SUBMODULE_CLI="$REPO_ROOT/ext/sylph/target/release/sylph"

if [[ -x "$EMBEDDED_CLI" ]]; then
    SYLPH="$EMBEDDED_CLI"
elif [[ -x "$SUBMODULE_CLI" ]]; then
    SYLPH="$SUBMODULE_CLI"
    echo "Note: using submodule sylph build at $SUBMODULE_CLI (Phase 0/1 mode)" >&2
else
    echo "ERROR: no sylph CLI found. Build first via 'bash build.sh', or build the" >&2
    echo "       submodule directly via 'cd ext/sylph && cargo build --release'." >&2
    echo "       (For the latter, you may need to strip conda from PATH:" >&2
    echo "        env -i HOME=\$HOME PATH=/usr/bin:/usr/local/bin:\$HOME/.cargo/bin cargo build --release)" >&2
    exit 1
fi

cd data/sylph

echo "Building tiny.syldb from 3 E. coli reference genomes..."
"$SYLPH" sketch \
    tiny_refs/e.coli-EC590.fasta.gz \
    tiny_refs/e.coli-K12.fasta.gz \
    tiny_refs/e.coli-o157.fasta.gz \
    -o tiny

echo "Profiling tiny K12 paired reads against tiny.syldb..."
"$SYLPH" profile tiny.syldb \
    -1 tiny_reads_R1.fq.gz \
    -2 tiny_reads_R2.fq.gz \
    -o expected_profile.tsv

echo "Done. Updated:"
echo "  data/sylph/tiny.syldb"
echo "  data/sylph/expected_profile.tsv"
echo
echo "Update the SHA pin (after committing the embedded fork):"
echo "  ( cd ext/sylph && git rev-parse HEAD ) > data/sylph/tiny_oracle.submodule.sha"
