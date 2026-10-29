set -e

# Guard against the make-clean-build-target race (issue #83). A single `make`
# invocation with both `clean` and a build goal lets -j run them concurrently,
# and clean ends up rm-ing .o files mid-archive.
if grep -rnE --include='*.cmake' --include=CMakeLists.txt '\bmake\b.*\bclean\b.*(\blib[A-Za-z0-9_]+\.a\b|\ball\b)' CMakeLists.txt cmake/; then
    echo "ERROR: BUILD_COMMAND has 'make ... clean <target>' on one line (race per #83)." >&2
    echo "Split into two sequential COMMAND steps." >&2
    exit 1
fi

# MIINT_TEST_JOBS=N runs the SQL files and the C++ tests on N parallel workers
# (test/scripts/run_tests_parallel.py) instead of one after another. The same tests run
# either way. Unset, as in CI, it is 1: the serial pass. Checked before the setup below.
MIINT_TEST_JOBS="${MIINT_TEST_JOBS:-1}"
case "$MIINT_TEST_JOBS" in
    *[!0-9]* | 0*)
        echo "ERROR: MIINT_TEST_JOBS must be a positive integer, got '$MIINT_TEST_JOBS'" >&2
        exit 1
        ;;
esac

# mmvec compiles its objective once per instruction set and picks at load time, so
# a fit's result depends on the CPU running it -- the same property scikit-bio has
# through OpenBLAS, and an accepted one. Expected values in the SQL tests are
# carved against the BASELINE kernel, so pin it here: without this the suite would
# assert different numbers on an AVX-512 CI runner than on an AVX2 one, and a real
# regression would be indistinguishable from a change of machine.
#
# The wide variants are not left untested by this. test/cpp/test_MMvec.cpp calls
# each one directly -- no environment variable involved, since DetectIsa() memoizes
# and one process can only ever observe one dispatch decision -- and checks them
# against baseline within the carved kIsa*Tol bands.
#
# Honours an externally-set value so `MIINT_SIMD=avx512 bash run_tests.sh` still
# works for a deliberate cross-check.
export MIINT_SIMD="${MIINT_SIMD:-baseline}"
echo "MIINT_SIMD=$MIINT_SIMD (mmvec kernel; expected values are carved against 'baseline')"

# Under this config every error fails its test; by default unittest skips some
# (docs/testing.md, "Guards and skipped tests").
export DUCKDB_TEST_CONFIG="${DUCKDB_TEST_CONFIG:-$PWD/test/configs/no_skipped_errors.json}"

# Start local HTTP server for HTTPS reader tests (unless already set externally)
HTTP_SERVER_PID=""
if [ -z "$MIINT_HTTPS_TEST_URL" ]; then
    # Find a free port: bind to port 0, read assigned port, close immediately
    HTTPS_TEST_PORT=$(python3 -c "import socket; s=socket.socket(); s.bind(('',0)); print(s.getsockname()[1]); s.close()")
    # Not `python3 -m http.server`: it ignores Range headers, see the script.
    python3 test/scripts/range_http_server.py "$HTTPS_TEST_PORT" data/ > /dev/null 2>&1 &
    HTTP_SERVER_PID=$!
    trap "kill $HTTP_SERVER_PID 2>/dev/null || true" EXIT

    # Wait for server to be ready (up to 3 seconds)
    HTTP_SERVER_READY=0
    for i in $(seq 1 15); do
        if curl -s "http://127.0.0.1:$HTTPS_TEST_PORT/" > /dev/null 2>&1; then
            HTTP_SERVER_READY=1
            break
        fi
        sleep 0.2
    done

    if [ "$HTTP_SERVER_READY" -ne 1 ]; then
        echo "ERROR: HTTP test server failed to start on port $HTTPS_TEST_PORT"
        exit 1
    fi

    export MIINT_HTTPS_TEST_URL="http://127.0.0.1:$HTTPS_TEST_PORT"
fi

# Detect optional external tools
if command -v ascp &> /dev/null; then
    export ASPERA_AVAILABLE=1
fi

# gpl-boundary daemon (process-isolation host for GPL tools incl. fasttree).
# When present we exercise the full phylogeny_fasttree round trip; when not,
# parity tests skip gracefully.
if command -v gpl-boundary &> /dev/null; then
    export GPL_BOUNDARY_AVAILABLE=1
fi

# miint looks for gpl-boundary in its install cache before PATH, but DuckDB 2.0's unittest overrides $HOME, so
# the tests can't see ~/.cache/miint/bin and fall back to whatever (possibly older) copy is on PATH. Point at the
# cached binary explicitly (same lookup as MiintGplBoundaryCacheDir) unless the caller already chose one.
GPL_BOUNDARY_CACHED="${XDG_CACHE_HOME:-$HOME/.cache}/miint/bin/gpl-boundary"
if [ -z "${MIINT_GPL_BOUNDARY_PATH:-}" ] && [ -x "$GPL_BOUNDARY_CACHED" ]; then
    export MIINT_GPL_BOUNDARY_PATH="$GPL_BOUNDARY_CACHED"
fi

# skb-arrow (process host for scikit-bio). Same $HOME problem: point at the binary install_skb_arrow() puts in
# miint's cache unless the caller chose one, then guard on whichever binary miint will run.
SKB_ARROW_CACHED="${XDG_CACHE_HOME:-$HOME/.cache}/miint/bin/skb-arrow"
if [ -z "${MIINT_SKB_ARROW_PATH:-}" ] && [ -x "$SKB_ARROW_CACHED" ]; then
    export MIINT_SKB_ARROW_PATH="$SKB_ARROW_CACHED"
fi
if [ -x "${MIINT_SKB_ARROW_PATH:-$(command -v skb-arrow)}" ]; then
    export SKB_ARROW_AVAILABLE=1
fi

# Independent native FastTree binary (bioconda or built from source), used only
# to regenerate the phylogeny_fasttree parity oracles below. Distinct from
# gpl-boundary's vendored ext/fasttree library, which is what the test drives.
FASTTREE_AVAILABLE=""
if command -v FastTree &> /dev/null; then
    FASTTREE_AVAILABLE=1
fi

# phylogeny_fasttree parity oracles.
#
# data/fasttree/{tiny,moderate}.fa are committed input fixtures; their
# *.golden.nwk siblings are the bioconda FastTree binary's outputs for those
# fixtures with -seed 42 -nosupport (regenerated by setting
# MIINT_FASTTREE_REGENERATE=1 below). Two sidecars per fixture:
#   *.golden.fixture.sha       — sha256 of the input fixture. The parity test
#                                compares it with the fixture and fails when the
#                                fasta was edited without an oracle refresh.
#   *.golden.fasttree.version  — `FastTree | grep Version` of the binary that
#                                generated the oracle.
#
# To regenerate: rm data/fasttree/*.golden.{nwk,fixture.sha,fasttree.version}
# and rerun run_tests.sh with MIINT_FASTTREE_REGENERATE=1.
ft_get_version() {
    # FastTree writes its banner to stderr ("FastTree Version 2.2.0 ..."). The
    # 2>&1 fold + first-line grep is enough; if FastTree isn't on PATH the
    # grep returns empty and the caller can detect the empty string.
    FastTree 2>&1 < /dev/null | head -1 | grep -oE 'Version [0-9][0-9.]*' | head -1
}

ft_regenerate_one() {
    # Args: fixture, golden, sha_file, version_file
    # Captures FastTree stderr to a sibling .regen.log on failure so the
    # diagnostic survives — /dev/null was throwing it away.
    local fixture="$1"
    local golden="$2"
    local sha_file="$3"
    local version_file="$4"
    if [ "${MIINT_FASTTREE_REGENERATE:-0}" != "1" ]; then
        return 0
    fi
    if [ -z "${FASTTREE_AVAILABLE:-}" ]; then
        echo "Warning: MIINT_FASTTREE_REGENERATE=1 but FastTree not on PATH; oracles not refreshed"
        return 0
    fi
    if [ ! -f "$fixture" ]; then
        echo "Warning: fixture $fixture missing; cannot regenerate oracle"
        return 0
    fi
    local logfile="${golden%.nwk}.regen.log"
    echo "Regenerating FastTree parity oracle from $fixture..."
    if FastTree -nt -seed 42 -nosupport "$fixture" > "$golden" 2> "$logfile"; then
        sha256sum "$fixture" | awk '{print $1}' > "$sha_file"
        local v
        v=$(ft_get_version)
        if [ -n "$v" ]; then
            echo "$v" > "$version_file"
        fi
        rm -f "$logfile"
        echo "FastTree oracle regenerated at $golden ($v)"
    else
        rm -f "$golden"
        echo "Warning: FastTree oracle regeneration failed for $fixture; diagnostic at $logfile"
    fi
}

ft_regenerate_one "data/fasttree/tiny.fa"     "data/fasttree/tiny.golden.nwk"     "data/fasttree/tiny.golden.fixture.sha"     "data/fasttree/tiny.golden.fasttree.version"
ft_regenerate_one "data/fasttree/moderate.fa" "data/fasttree/moderate.golden.nwk" "data/fasttree/moderate.golden.fixture.sha" "data/fasttree/moderate.golden.fasttree.version"

# ASR (ancestral-state) parity goldens: independent ground-truth reconstructions
# (Brownian-motion GLS via the phylogenetic VCV; ape::ace for the discrete Mk models).
# The NUMERIC CSVs under data/asr/ ARE the fixed expected output -- COMMITTED, so they
# are always present and need no availability gate. The offline generator is a dev-only
# R + ape (GPL) tool deliberately kept out of this BSD tree
# (../duckdb-miint-localdocs/gen_asr_oracle.R); ape's code is never committed or
# distributed -- only its numeric output is.
#
# These previously sat behind a require-env MIINT_ASR_PARITY_OK gate driven by
# `sha256sum -c` against a data/asr/goldens.sha256 manifest. Both the gate and the
# manifest were removed for the same reason as the data/simsurvey/ ones below:
# `sha256sum -c` also fails on a MISMATCH, so a corrupted or edited golden SKIPPED the
# parity test and left CI green -- the exact opposite of the intended protection, and a
# Rule 10 (fail loud) violation. The parity tests themselves are the integrity check,
# and they now always run: a corrupted golden fails them, and a missing one fails the
# read_csv_auto that loads it.

# NOTE: the Kuczynski-2010 oracle bands, community_distances distance goldens and
# cluster_kmeans/cluster_upgma parity goldens under data/simsurvey/ are COMMITTED,
# so they are always present and need no availability gate. They previously sat
# behind sha256 gates; those were removed because `sha256sum -c` also fails on a
# MISMATCH, which turned an edited golden into a silently SKIPPED parity test --
# the opposite of the intended protection. The parity tests themselves are the
# integrity check, and they now always run.

# Procrustes parity goldens: independent numeric ground truth from SciPy
# (scipy.spatial.procrustes full case) + a NumPy port of the q2#338 partial
# technique, cross-checked against the shipped q2_diversity._partial_procrustes.
# The expected CSVs under data/procrustes/ are the committed, SHA-pinned output
# (data/procrustes/goldens.sha256); the generator (test/scripts/gen_procrustes_oracle.py)
# needs only SciPy/NumPy (BSD, in-tree-compatible). procrustes_parity.test checks
# the goldens against that manifest and fails on a mismatch. To refresh: rerun
# the generator under a SciPy env.

# Live NCBI BLAST tests, opt-in. miint's BLAST client waits up to about two hours per
# submission (src/blast_client.cpp), and that the endpoint answers says nothing about
# NCBI's queue (one submission was still queued after 28 minutes on 2026-10-15). That is
# past the httpfs pass's default limit below, so raise it too:
# `MIINT_BLAST_LIVE_TEST=1 MIINT_HTTPFS_TIMEOUT_SECONDS=10800 bash run_tests.sh`. With the
# opt-in set, the probe still skips them when BLAST does not answer.
if [ "${MIINT_BLAST_LIVE_TEST:-0}" = "1" ] &&
    curl -sSf --max-time 5 -o /dev/null \
        "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi" 2>/dev/null; then
    export BLAST_LIVE_AVAILABLE=1
fi

# Live ENA tests, opt-in: `MIINT_ENA_LIVE_TEST=1 bash run_tests.sh` (why: docs/testing.md,
# "Guards and skipped tests"). They run under with_live_ena: without the mock's URLs, which
# this script exports below, and with ENA_AVAILABLE, which is set nowhere else. With the
# opt-in set, the probe still skips them when the Portal does not answer.
source test/scripts/ena_mock_env.sh
unset ENA_AVAILABLE
ENA_LIVE=0
with_live_ena() {
    (
        unset_ena_mock_env
        export ENA_AVAILABLE=1
        "$@"
    )
}
if [ "${MIINT_ENA_LIVE_TEST:-0}" = "1" ]; then
    if curl -sSf --max-time 10 -o /dev/null \
            "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession%3D%22ERR1074767%22&fields=run_accession&limit=1&format=tsv" \
            2>/dev/null; then
        ENA_LIVE=1
    else
        echo "Warning: MIINT_ENA_LIVE_TEST=1, but ENA's Portal API did not answer; the live ENA tests will skip"
    fi
fi

# NCBI network reachability, for the live read_ncbi* tests: they use both E-utilities and
# the Datasets API, so both must answer.
if curl -sSf --max-time 5 -o /dev/null "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi" 2>/dev/null &&
    curl -sSf --max-time 5 -o /dev/null "https://api.ncbi.nlm.nih.gov/datasets/v2/version" 2>/dev/null; then
    export NCBI_AVAILABLE=1
fi

# NCBI taxonomy dump reachability (opt-in integration test; the test downloads
# ~60MB). Cache into a contained repo dir (not the user's ~/.cache/miint) so the
# test is self-contained and re-runs hit the cache instead of re-downloading.
if curl -sSf --max-time 10 -I -o /dev/null \
        "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" 2>/dev/null; then
    export MIINT_TAXDUMP_REMOTE=1
    export MIINT_TAXONOMY_CACHE_DIR="$PWD/data/cache/taxonomy"
    mkdir -p "$MIINT_TAXONOMY_CACHE_DIR"
fi

# ena_upload_reads streaming memory proof (opt-in, expensive). Uploads a
# multi-GB lazy view to file:// under a tight memory_limit; the streaming
# implementation completes in bounded memory whereas the old materialise-the-
# whole-relation path would OOM. Off by default (slow, multi-GB encode); run
# with `MIINT_ENA_BIG_STREAM_TEST=1 bash run_tests.sh` to include it.
if [ "${MIINT_ENA_BIG_STREAM_TEST:-0}" = "1" ]; then
    export MIINT_ENA_BIG_STREAM_TEST=1
fi

# ENA Webin V2 submission test endpoint. Live-integration tests for
# INSERT INTO ena.{projects,samples,...} require a real Webin account
# (free, registered at https://www.ebi.ac.uk/ena/submit/webin/). When the
# user provides credentials AND the dev/test server responds, we export
# ENA_WEBIN_TEST_AVAILABLE so guarded SQL tests can run.
if [ -n "$ENA_WEBIN_TEST_USER" ] && [ -n "$ENA_WEBIN_TEST_PASSWORD" ]; then
    export ENA_WEBIN_TEST_USER
    export ENA_WEBIN_TEST_PASSWORD
    if curl -sSf --max-time 3 -o /dev/null \
            "https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/swagger-ui/index.html" \
            2>/dev/null; then
        export ENA_WEBIN_TEST_AVAILABLE=1
    fi
fi

# The ENA mock (test/scripts/ena_webin_mock.py): canned receipts for INSERT INTO
# ena.{projects,samples,...} with no live credentials, and the Portal answers and files
# for the runs in data/ena/. Portal, checklist and Reports requests go to it, except under
# with_live_ena.
ENA_MOCK_PID=""
if [ -z "$ENA_WEBIN_MOCK_URL" ] && python3 -c "import http.server" 2>/dev/null; then
    # Inherit the existing trap (HTTP_SERVER_PID) and add the mock kill, set before the
    # mock starts so an exit while it boots still stops it.
    trap 'kill $HTTP_SERVER_PID $ENA_MOCK_PID 2>/dev/null || true' EXIT
    if ! start_ena_mock; then
        echo "Warning: the ENA mock failed to start; its tests will skip, which fails the httpfs pass"
    fi
fi

if conda run -n massql python3 -c "from massql import msql_engine" 2>/dev/null; then
    export MASSQL_PYTHON_AVAILABLE=1
fi

# Concurrent-LOAD test, opt-in because it takes minutes: `MIINT_LOAD_RACE_TEST=1 bash run_tests.sh`
# (CI sets it nightly). The extension loads only into the DuckDB it was built for, so Python's
# duckdb must match.
if [ "${MIINT_LOAD_RACE_TEST:-0}" = "1" ]; then
    export MIINT_TEST_PYTHON="${MIINT_TEST_PYTHON:-python3}"
    if [ "v$("$MIINT_TEST_PYTHON" -c 'import duckdb; print(duckdb.__version__)' 2>/dev/null)" = \
        "$(./build/release/duckdb -csv -noheader -c 'SELECT version();' 2>/dev/null)" ]; then
        export DUCKDB_PYTHON_AVAILABLE=1
    else
        echo "ERROR: MIINT_LOAD_RACE_TEST=1 but $MIINT_TEST_PYTHON has no duckdb matching this build (set MIINT_TEST_PYTHON)" >&2
        exit 1
    fi
fi

# Download and cache MassQL test data if not present
CACHE_DIR="data/cache"
mkdir -p "$CACHE_DIR"
MZML_FILE="$CACHE_DIR/bld_plt1_07_120_1.mzML"
EXPECTED_SHA="9be69d251cbf7d9b3c439ccc76dd65f055ddac96b946fc2f565e9a95bfd4ca46"

if [ ! -f "$MZML_FILE" ]; then
    echo "Downloading MassQL test data..."
    if ! wget -q -O "$MZML_FILE" "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?forceDownload=true&file=f.MSV000085944/ccms_peak/raw_data/bld_plt1_07_120_1.mzML"; then
        rm -f "$MZML_FILE"
        echo "Warning: MassQL test data download failed, skipping MassQL tests"
    fi
fi
if echo "$EXPECTED_SHA  $MZML_FILE" | sha256sum -c --quiet 2>/dev/null; then
    export MASSQL_TEST_DATA="$MZML_FILE"
fi

# Download and cache GNPS test data if not present
GNPS_FILE="$CACHE_DIR/GNPS00002_A3_p.mzML"
if [ ! -f "$GNPS_FILE" ]; then
    echo "Downloading GNPS test data..."
    if ! wget -q -O "$GNPS_FILE" "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?forceDownload=true&file=f.MSV000084494/ccms_peak/raw/GNPS00002_A3_p.mzML"; then
        rm -f "$GNPS_FILE"
        echo "Warning: GNPS test data download failed, skipping GNPS tests"
    fi
fi
if [ -f "$GNPS_FILE" ]; then
    export MASSQL_GNPS_DATA=1
fi

# Download and cache real mzXML test data if not present
MZXML_FILE="$CACHE_DIR/10125_P1_RA12_01_49.mzXML"
MZXML_EXPECTED_SHA="d026938157e5f640568c371ed6a0bca0de34e8c181da1443298671005e7ffe5a"
if [ ! -f "$MZXML_FILE" ]; then
    echo "Downloading mzXML test data..."
    if ! wget -q -O "$MZXML_FILE" "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000080179%2Fccms_peak%2FmzXMLs%2F10125_P1_RA12_01_49.mzXML&forceDownload=true"; then
        rm -f "$MZXML_FILE"
        echo "Warning: mzXML test data download failed, skipping mzXML real data tests"
    fi
fi
if echo "$MZXML_EXPECTED_SHA  $MZXML_FILE" | sha256sum -c --quiet 2>/dev/null; then
    export MZXML_REAL_DATA="$MZXML_FILE"
fi

# Detect compile-time optional features by querying the built extension
if echo "SELECT * FROM miint_versions() WHERE library = 'HDF5';" | ./build/release/duckdb -csv 2>/dev/null | grep -q HDF5; then
    export HDF5_AVAILABLE=1
fi
if echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'detect_chimera_uchime';" | ./build/release/duckdb -csv -noheader 2>/dev/null | grep -q 1; then
    export VSEARCH_AVAILABLE=1
fi
if echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'align_mafft';" | ./build/release/duckdb -csv -noheader 2>/dev/null | grep -q 1; then
    export MAFFT_AVAILABLE=1
fi
if echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'sourcetracker';" | ./build/release/duckdb -csv -noheader 2>/dev/null | grep -q 1; then
    export ST3_AVAILABLE=1
fi
# krepp's index regions are only compiled in when an OpenMP runtime was found at
# configure time; without them krepp_index_create refuses threads > 1. Both
# branches are exported, because the refusal is worth asserting too - a build
# that took the parameter and silently ran on one core is the failure this is
# guarding against.
# Unset first: these two are mutually exclusive assertions, not one availability
# flag, so a value inherited from the caller's environment does not merely skip a
# file - it runs the wrong one and fails it. Verified: KREPP_OPENMP_ABSENT=1
# pre-set on an OpenMP build makes krepp_index_create_no_threads.test fail.
unset KREPP_OPENMP_AVAILABLE KREPP_OPENMP_ABSENT
if echo "SELECT * FROM miint_versions() WHERE library = 'krepp';" | ./build/release/duckdb -csv 2>/dev/null | grep -q krepp; then
    if echo "SELECT * FROM miint_versions() WHERE library = 'krepp-openmp';" | ./build/release/duckdb -csv 2>/dev/null | grep -q krepp-openmp; then
        export KREPP_OPENMP_AVAILABLE=1
    else
        export KREPP_OPENMP_ABSENT=1
    fi
fi
if echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'align_abpoa';" | ./build/release/duckdb -csv 2>/dev/null | grep -q 1; then
    export ABPOA_AVAILABLE=1
fi
if echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'align_sortmerna_rrna';" | ./build/release/duckdb -csv -noheader 2>/dev/null | grep -q 1; then
    export SORTMERNA_AVAILABLE=1
fi
# The skb-arrow table functions (MIINT_HAS_GPL_BOUNDARY), for tests that need them bound but not run.
if echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'skb_mantel' AND function_type = 'table';" | ./build/release/duckdb -csv -noheader 2>/dev/null | grep -q 1; then
    export SKB_ARROW_COMPILED=1
fi
if echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'unifrac_pcoa' AND function_type = 'table';" | ./build/release/duckdb -csv -noheader 2>/dev/null | grep -q 1; then
    export UNIFRAC_AVAILABLE=1
fi
if echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'sylph_profile';" | ./build/release/duckdb -csv -noheader 2>/dev/null | grep -q 1; then
    export SYLPH_AVAILABLE=1
fi
if echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'place_krepp';" | ./build/release/duckdb -csv -noheader 2>/dev/null | grep -q 1; then
    export KREPP_AVAILABLE=1
fi

# Toy krepp index for the place_krepp end-to-end test.
#
# krepp ships no test suite. What it ships is the tutorial data in
# ext/krepp/test/ - 25 genomes, a rooted guide tree and 100 reads - plus a
# README quickstart that builds an index from them. We rebuild that index here
# rather than committing one: it is ~69 MB of derived data (the cmer table is
# 8 bytes per reference k-mer, so no realistic genome set makes it small), and
# its on-disk layout is tied to the krepp that wrote it - v0.9.0 swapped the
# recursion branches in Node::generate_tree, so an index written by a different
# version is misread rather than rejected.
#
# Needs the krepp binary and xz to unpack the references. Without either,
# test/sql/place_krepp_toy.test skips; the validation tests in
# test/sql/place_krepp.test still run.
#
# NOTE ON THE BINARY. As of 2026-09-10 bioconda's newest krepp is 0.9.1, while
# the submodule pins 31205033, whose VERSION string is v0.10.0 and which is not
# tagged - so `conda install bioconda::krepp` gives a binary older than the
# linked library. Build it from source at the submodule's commit
# (`git clone https://github.com/bo1929/krepp && cd krepp &&
# git checkout <ext/krepp's commit> && git submodule update --init --recursive
# && make`). The stamp below records the CLI banner precisely so a mismatch
# forces a rebuild rather than being silently reused, but nothing compares it
# against the linked version - see MIINT_KREPP_TOY_INDEX in the docs.
#
# The stamp pins three things, so the ~20 s build happens once and repeats only
# when one of them moves. Everything lands in data/krepp/, which is gitignored.
#   - ext/krepp's submodule HEAD, because that is what the extension links.
#   - the CLI's version banner, because that is what writes the index, and it is
#     a different build from the linked one (a separately compiled CLI).
#   - the sha of the reference tarball, i.e. the input data itself.
# The banner alone is not enough: it expands PRINT_VERSION, a hardcoded string
# in common.hpp that does not change with every commit. Two commits
# either side of the generate_tree change both report v0.9.0 while writing
# incompatible indexes - precisely the drift described above. Same shape as
# SORTMERNA_REAL_ORACLE, which records the submodule sha for the same reason.
#
# The index directory is NAMED for the stamp rather than stamped alongside it.
# `krepp index` only creates the directory and writes files whose names encode
# the resolved config (cmer-m<m>r<r>-{frac,no_frac}); it never clears what is
# already there. Rebuilding in place after a krepp change that moves m or r
# would leave the old partials beside the new ones, and DiscoverPartials would
# find two complete suffix groups and load both - exactly the mixed-version
# index this stamp exists to prevent. A fresh directory per stamp sidesteps it
# without deleting anything (CLAUDE.md: never `rm` without permission), at the
# cost of leaving old indexes on disk; data/krepp/ is gitignored, so remove them
# by hand when you care.
KREPP_TOY_DIR=data/krepp

# krepp's tutorial references, unpacked once.
#
# Hoisted out of the CLI-gated build below because krepp_index_create needs the
# references and nothing else - it builds the index in-process, so its
# round-trip test runs on a machine with no krepp binary at all. The CLI build
# further down reuses what this unpacks.
krepp_toy_ref_count() { ls "$KREPP_TOY_DIR"/references_toy/*.fna 2>/dev/null | wc -l | tr -d ' '; }
if [ -n "$KREPP_AVAILABLE" ] && command -v xz &> /dev/null \
   && [ -f ext/krepp/test/references_toy.tar.gz ] && [ -s ext/krepp/test/input_map.tsv ]; then
    # How many there are supposed to be, from the corpus's own manifest, so a
    # half-finished extraction is retried instead of being mistaken for a
    # complete one by "at least one .fna exists".
    #
    # `grep -c` exits 1 when it counts zero lines and 2 when the file is
    # missing. Here it is the whole right-hand side of an assignment, so set -e
    # checks it and either would abort the entire test run. The `-s` test above
    # rules both out; `|| true` is belt for whoever edits that gate next.
    KREPP_TOY_WANT_REFS="$(grep -c . ext/krepp/test/input_map.tsv || true)"
    if [ "$(krepp_toy_ref_count)" != "$KREPP_TOY_WANT_REFS" ]; then
        mkdir -p "$KREPP_TOY_DIR"
        # `|| true` is load-bearing. set -e exempts every command in an AND-OR
        # list EXCEPT the last one, so the last is exactly the position it does
        # check - and xz exits 1 if any member fails to decode. Without the
        # trailing `|| true` that failure aborts the entire test run, every SQL
        # and C++ test after this point, instead of skipping one optional test.
        # The count check below is what decides whether the unpack worked.
        tar -xzf ext/krepp/test/references_toy.tar.gz -C "$KREPP_TOY_DIR" \
            && xz -df "$KREPP_TOY_DIR"/references_toy/*.fna.xz 2>/dev/null || true
    fi
    if [ "$(krepp_toy_ref_count)" = "$KREPP_TOY_WANT_REFS" ]; then
        export MIINT_KREPP_TOY_REFS="$KREPP_TOY_DIR/references_toy"
    else
        echo "Warning: unpacked $(krepp_toy_ref_count) of $KREPP_TOY_WANT_REFS krepp toy references,"
        echo "         so krepp_index_create's round-trip test is being skipped."
    fi
elif [ -n "$KREPP_AVAILABLE" ]; then
    # Say which of the two it was. A silent skip here looks identical to a
    # passing round-trip test.
    if ! command -v xz &> /dev/null; then
        echo "Note: xz not on PATH, so krepp's toy references cannot be unpacked and"
        echo "      krepp_index_create's round-trip test is being skipped."
    else
        echo "Note: krepp's toy reference corpus (ext/krepp/test/references_toy.tar.gz"
        echo "      and input_map.tsv) is missing or empty, so krepp_index_create's"
        echo "      round-trip test is being skipped."
    fi
fi
if [ -n "$KREPP_AVAILABLE" ] && command -v krepp &> /dev/null && command -v xz &> /dev/null; then
    # Every component is checked, because an empty one would silently match an
    # empty one on the other side and turn the pin into a no-op.
    KREPP_TOY_WANT=""
    KREPP_TOY_SHA="$(git -C ext/krepp rev-parse HEAD 2>/dev/null || true)"
    KREPP_TOY_BANNER="$(krepp --help 2>&1 | head -1 || true)"
    KREPP_TOY_TARBALL="$(sha256sum ext/krepp/test/references_toy.tar.gz 2>/dev/null | awk '{print $1}')"
    if [ -z "$KREPP_TOY_SHA" ] || [ -z "$KREPP_TOY_BANNER" ] || [ -z "$KREPP_TOY_TARBALL" ]; then
        echo "Warning: cannot stamp krepp toy index (submodule sha, krepp banner or tarball sha unavailable);"
        echo "         skipping place_krepp end-to-end tests"
    else
        KREPP_TOY_WANT="$KREPP_TOY_SHA|$KREPP_TOY_BANNER|$KREPP_TOY_TARBALL"
    fi
fi
if [ -n "$KREPP_TOY_WANT" ]; then
    KREPP_TOY_KEY="$(printf '%s' "$KREPP_TOY_WANT" | sha256sum | cut -c1-12)"
    KREPP_TOY_INDEX="$KREPP_TOY_DIR/index_toy-$KREPP_TOY_KEY"
    # Built under .partial and moved into place only on success. `krepp index`
    # creates its output directory in the CLI callback, before it indexes
    # anything, so a failure at any point leaves a directory behind that exists
    # and is incomplete. Testing `-d` on the final path would then both export
    # that partial index and, because the same `-d` guards the build, never
    # retry it. The mv is atomic, and nothing is deleted.
    KREPP_TOY_PARTIAL="$KREPP_TOY_INDEX.partial"
    if [ ! -d "$KREPP_TOY_INDEX" ] && [ ! -d "$KREPP_TOY_PARTIAL" ]; then
        echo "Building krepp toy index (once per krepp version) ..."
        mkdir -p "$KREPP_TOY_DIR"
        if [ -n "$MIINT_KREPP_TOY_REFS" ] \
           && awk -v d="$PWD/$KREPP_TOY_DIR" -F'\t' '{print $1 "\t" d "/references_toy/" $1 ".fna"}' \
                ext/krepp/test/input_map.tsv > "$KREPP_TOY_DIR/input_map.tsv" \
           && krepp index -h 11 -k 27 -w 35 -o "$KREPP_TOY_PARTIAL" \
                -i "$KREPP_TOY_DIR/input_map.tsv" -t ext/krepp/test/tree_toy.nwk \
                --num-threads 4 > "$KREPP_TOY_DIR/index_toy-$KREPP_TOY_KEY.log" 2>&1 \
           && mv "$KREPP_TOY_PARTIAL" "$KREPP_TOY_INDEX"; then
            :
        else
            echo "Warning: krepp toy index build failed, skipping place_krepp end-to-end tests"
            echo "         (see $KREPP_TOY_DIR/index_toy-$KREPP_TOY_KEY.log; the incomplete build is"
            echo "          left at $KREPP_TOY_PARTIAL — remove it to retry)"
        fi
    elif [ ! -d "$KREPP_TOY_INDEX" ]; then
        # Reached when a previous run failed and parked a .partial. Without this
        # branch the guard above is simply false, nothing builds, nothing warns,
        # and place_krepp_toy.test skips forever with CI green.
        echo "Warning: an incomplete krepp toy index is parked at $KREPP_TOY_PARTIAL,"
        echo "         so place_krepp end-to-end tests are being skipped."
        echo "         Remove that directory to retry the build."
    fi
    if [ -d "$KREPP_TOY_INDEX" ]; then
        export MIINT_KREPP_TOY_INDEX="$KREPP_TOY_INDEX"
    fi
elif [ -n "$KREPP_AVAILABLE" ]; then
    # Say so. bioconda's newest krepp is older than the pinned commit, so "CLI not
    # found" is the expected case rather than the exotic one, and a silent skip
    # here looks identical to a passing end-to-end suite.
    if ! command -v krepp &> /dev/null; then
        echo "Note: krepp CLI not on PATH, so the toy index cannot be built and"
        echo "      place_krepp end-to-end tests are being skipped. Build krepp from"
        echo "      source at ext/krepp's commit (https://github.com/bo1929/krepp) to"
        echo "      run them; bioconda's build is older than the linked library."
    elif ! command -v xz &> /dev/null; then
        echo "Note: xz not on PATH, so krepp's reference tarball cannot be unpacked"
        echo "      and place_krepp end-to-end tests are being skipped."
    fi
fi
# libcurl streaming-upload transport (off on macOS — vsearch/OpenSSL symbol clash).
if echo "SELECT 1 FROM miint_versions() WHERE library = 'libcurl';" | ./build/release/duckdb -csv -noheader 2>/dev/null | grep -q 1; then
    export MIINT_HAS_CURL=1
fi
# Compile-time gate (separate from runtime GPL_BOUNDARY_AVAILABLE which
# tracks whether the binary is on PATH). Set when the table function is
# registered — i.e., MIINT_HAS_GPL_BOUNDARY was on at build time.
if echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'phylogeny_fasttree' AND function_type = 'table';" | ./build/release/duckdb -csv -noheader 2>/dev/null | grep -q 1; then
    export PHYLOGENY_FASTTREE_AVAILABLE=1
fi

# Multi-part minimap2 index fixture for align_minimap2's streaming path
# (align_minimap2_multipart.test). minimap2 has no SQL-level way to build a
# multi-part .mmi (save_minimap2_index always builds single-part, via
# mm_idx_str); the CLI's own `-I <batch>` flag produces one by writing
# multiple mm_idx_dump blocks into the same file. We reproduce that exact
# on-disk layout without depending on a minimap2 CLI binary: build two tiny
# single-part indexes via save_minimap2_index, then concatenate their bytes
# (each is already a self-contained MM_IDX_MAGIC-prefixed dump).
#
# align_minimap2_multipart.test hardcodes query sequences that are PREFIXES of
# the two reference sequences below. Change one without the other and the test
# quietly stops matching instead of failing loudly — edit them together.
MINIMAP2_MULTIPART_DIR="data/shards"
MINIMAP2_MULTIPART_PART1="$MINIMAP2_MULTIPART_DIR/multipart_fixture_part1.mmi"
MINIMAP2_MULTIPART_PART2="$MINIMAP2_MULTIPART_DIR/multipart_fixture_part2.mmi"
MINIMAP2_MULTIPART_MMI="$MINIMAP2_MULTIPART_DIR/multipart_fixture.mmi"
mkdir -p "$MINIMAP2_MULTIPART_DIR"
MULTIPART_GEN_SQL="
CREATE TABLE _multipart_fixture_p1 AS SELECT * FROM (VALUES
    ('part1_ref', 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC')
) AS t(read_id, sequence1);
CREATE TABLE _multipart_fixture_p2 AS SELECT * FROM (VALUES
    ('part2_ref', 'TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAA')
) AS t(read_id, sequence1);
SELECT success FROM save_minimap2_index('_multipart_fixture_p1', '$MINIMAP2_MULTIPART_PART1', k := 5);
SELECT success FROM save_minimap2_index('_multipart_fixture_p2', '$MINIMAP2_MULTIPART_PART2', k := 5);
"
# The two single-part files are kept (not just the concatenated multi-part
# file) so the .test can build a UNION ALL of two single-index_path calls as
# an oracle that is correct by construction and independent of the streaming
# code path being tested.
if echo "$MULTIPART_GEN_SQL" | ./build/release/duckdb -csv -noheader > /dev/null 2>&1 \
    && cat "$MINIMAP2_MULTIPART_PART1" "$MINIMAP2_MULTIPART_PART2" > "$MINIMAP2_MULTIPART_MMI" 2>/dev/null; then
    export MIINT_MINIMAP2_MULTIPART_FIXTURE="$MINIMAP2_MULTIPART_MMI"
else
    echo "Warning: failed to generate multi-part minimap2 index fixture; align_minimap2_multipart.test will skip"
fi

# Phase 5: real-data regression oracle.
#
# The BLAST-format output of `sortmerna --blast '1 cigar qcov' -e 1` on the
# in-tree 100k-read environmental amplicons vs. gg_13_8_ref_set.fasta is
# checked in at data/sortmerna/real_oracle.blast.gz (1.9 MB gzipped). The
# submodule SHA that produced it is recorded alongside so we can detect when
# the oracle has gone stale.
#
# `-e 1` is deliberate: the library reports every hit and applies no e-value
# filter, while the CLI's default threshold became 1e-5 in sortmerna 7.0.
# Pinning the CLI to its old default of 1 keeps this a comparison of
# alignments rather than of filtering.
#
# To regenerate after bumping ext/sortmerna: remove both files and rerun
# run_tests.sh with MIINT_SORTMERNA_REAL_DATA=1. That rebuilds the oracle
# using the submodule-built native CLI (the exact binary our library links),
# so the oracle and our embedded library always run the same version.
SMR_NATIVE_CLI="./build/release/extension/miint/sortmerna_build/src/sortmerna/sortmerna"
SMR_REAL_QUERY="ext/sortmerna/test_data/set2_environmental_study_550_amplicon.fasta"
SMR_REAL_REF="ext/sortmerna/test_data/gg_13_8_ref_set.fasta"
SMR_ORACLE_GZ="data/sortmerna/real_oracle.blast.gz"
SMR_ORACLE_SHA="data/sortmerna/real_oracle.submodule.sha"
SMR_SUBMODULE_HEAD=$(git -C ext/sortmerna rev-parse HEAD 2>/dev/null || echo "unknown")

if [ "${MIINT_SORTMERNA_REAL_DATA:-0}" = "1" ] && [ -x "$SMR_NATIVE_CLI" ] \
        && [ -f "$SMR_REAL_QUERY" ] && [ -f "$SMR_REAL_REF" ]; then
    echo "Regenerating sortmerna real-data oracle via the native CLI..."
    SMR_TMP=$(mktemp -d)
    if "$SMR_NATIVE_CLI" \
            --ref "$SMR_REAL_REF" --reads "$SMR_REAL_QUERY" \
            --workdir "$SMR_TMP" --aligned "$SMR_TMP/aligned" \
            --blast '1 cigar qcov' -e 1 --num_alignments 1 --threads 4 > "$SMR_TMP/smr.log" 2>&1 \
            && [ -f "$SMR_TMP/aligned.blast" ]; then
        mkdir -p data/sortmerna
        gzip -c "$SMR_TMP/aligned.blast" > "$SMR_ORACLE_GZ"
        echo "$SMR_SUBMODULE_HEAD" > "$SMR_ORACLE_SHA"
        echo "Oracle regenerated at $SMR_ORACLE_GZ (submodule $SMR_SUBMODULE_HEAD)"
        rm -rf "$SMR_TMP"
    else
        # Preserve the CLI log for debugging — the tmpdir is about to go.
        SMR_FAIL_LOG="$(pwd)/smr_oracle_gen_fail.log"
        cp "$SMR_TMP/smr.log" "$SMR_FAIL_LOG" 2>/dev/null || true
        echo "Warning: sortmerna oracle generation failed; log saved to $SMR_FAIL_LOG"
        rm -rf "$SMR_TMP"
    fi
fi

# Export the oracle path only if both the gzipped file and its SHA record
# exist, AND the recorded SHA matches the current submodule HEAD. A mismatch
# means the submodule was bumped without regenerating — skip the test rather
# than silently validate against stale goldens.
if [ -f "$SMR_ORACLE_GZ" ] && [ -f "$SMR_ORACLE_SHA" ]; then
    SMR_RECORDED_SHA=$(cat "$SMR_ORACLE_SHA")
    if [ "$SMR_RECORDED_SHA" = "$SMR_SUBMODULE_HEAD" ]; then
        export SORTMERNA_REAL_ORACLE="$SMR_ORACLE_GZ"
    else
        echo "Warning: sortmerna oracle ($SMR_RECORDED_SHA) does not match submodule HEAD ($SMR_SUBMODULE_HEAD); skipping real-data test"
    fi
fi

# The main pass, on MIINT_TEST_JOBS workers (checked at the top).
if [ "$MIINT_TEST_JOBS" -gt 1 ]; then
    python3 test/scripts/run_tests_parallel.py --jobs "$MIINT_TEST_JOBS" \
        --unittest ./build/release/test/unittest --cpp-tests ./build/release/extension/miint/tests
else
    make test
    ./build/release/extension/miint/tests
fi

# read_fastx_fd_release.test must run under a low RLIMIT_NOFILE: it asserts that one scan
# over many paths holds descriptors proportional to thread count, not to path count, and
# a process with the usual limit cannot tell those two apart. The main pass above has no way
# to set a limit, so the file holds itself back with `require-env MIINT_LOW_FD_TEST` and
# runs here instead.
#
# Lowering the soft limit needs no privilege. If a hard limit already sits below the target
# we run anyway rather than skipping -- a tighter limit only makes the assertion stronger.
FD_RELEASE_TEST="test/sql/read_fastx_fd_release.test"
FD_RELEASE_LIMIT=256
echo "Running $FD_RELEASE_TEST with RLIMIT_NOFILE<=$FD_RELEASE_LIMIT..."
(
    ulimit -n "$FD_RELEASE_LIMIT" 2>/dev/null || true
    echo "  RLIMIT_NOFILE soft=$(ulimit -Sn) hard=$(ulimit -Hn)"
    # The || true above swallows the expected failure (hard limit already below target) along
    # with any unanticipated one. Since the soft limit can never exceed the hard limit, the
    # expected case always leaves us at or under the target -- so a soft limit still above it
    # here means the lower silently didn't take, and the test below would run with no real
    # descriptor pressure and pass without exercising anything.
    if [ "$(ulimit -Sn)" -gt "$FD_RELEASE_LIMIT" ]; then
        echo "  ERROR: RLIMIT_NOFILE is still $(ulimit -Sn), above the $FD_RELEASE_LIMIT target;" \
             "refusing to run $FD_RELEASE_TEST without real descriptor pressure" >&2
        exit 1
    fi
    MIINT_LOW_FD_TEST=1 ./build/release/test/unittest -w NoTests "$FD_RELEASE_TEST"
) || exit 1

# sylph_index_create_thread_limit.test must run where no new thread can start: it asserts that
# sylph_index_create then sketches every genome on the calling thread. The main pass cannot set
# RLIMIT_NPROC, so the file holds itself back with `require-env MIINT_NPROC_LIMIT_TEST` and runs
# here. The subshell lowers the limit and then execs the runner, since a fork would be refused
# too. Only Linux counts threads against RLIMIT_NPROC, and root is not held to it; elsewhere the
# file would pass without the limit doing anything, so it is skipped with a note.
NPROC_LIMIT_TEST="test/sql/sylph_index_create_thread_limit.test"
if [ "$(uname -s)" != "Linux" ]; then
    echo "Skipping $NPROC_LIMIT_TEST: only Linux counts threads against RLIMIT_NPROC"
elif [ "$(id -u)" = "0" ]; then
    echo "Skipping $NPROC_LIMIT_TEST: root is not held to RLIMIT_NPROC"
else
    echo "Running $NPROC_LIMIT_TEST with RLIMIT_NPROC=1..."
    (
        ulimit -u 1 || exit 1
        MIINT_NPROC_LIMIT_TEST=1 exec ./build/release/test/unittest -w NoTests "$NPROC_LIMIT_TEST"
    ) || exit 1
fi

# The httpfs pass: every SQL file that needs httpfs, all of which `make test` above skips.
# docs/testing.md ("Guards and skipped tests") says why, and how this pass loads httpfs and
# this build's miint. The files' require-env guards still decide whether each runs: a live
# test runs only when its service answered the probe above.
HTTPFS_TESTS=$(grep -rl --include='*.test' '^require httpfs' test/sql | paste -sd, -)
HTTPFS_TEST_HOME="$PWD/build/release/httpfs_test_home"
mkdir -p "$HTTPFS_TEST_HOME"
# INSTALL httpfs downloads it on the first run only. -unsigned because this build is unsigned.
if ! HOME="$HTTPFS_TEST_HOME" ./build/release/duckdb -unsigned -c "INSTALL httpfs; LOAD httpfs;
        FORCE INSTALL '$PWD/build/release/extension/miint/miint.duckdb_extension';" > /dev/null; then
    echo "ERROR: could not install httpfs and this build's miint into $HTTPFS_TEST_HOME; the tests that need httpfs cannot run" >&2
    exit 1
fi
# A stalled remote transfer would otherwise hang the run instead of failing it.
HTTPFS_TIMEOUT_SECONDS="${MIINT_HTTPFS_TIMEOUT_SECONDS:-1800}"
TIMEOUT_CMD=$(command -v timeout || command -v gtimeout || true)
if [ -z "$TIMEOUT_CMD" ]; then
    echo "WARNING: no timeout or gtimeout command, so a stalled transfer hangs this run instead of failing it" >&2
fi
# Which skips are allowed: docs/testing.md, "Guards and skipped tests". A skip does not
# change unittest's exit status, so a disallowed one, such as a skip for the mock server's
# URL that this script exports, is caught here.
HTTPFS_ALLOWED_SKIPS='ENA_AVAILABLE|NCBI_AVAILABLE|BLAST_LIVE_AVAILABLE|MIINT_TAXDUMP_REMOTE|ENA_WEBIN_TEST_(USER|PASSWORD|AVAILABLE)|ASPERA_AVAILABLE|HDF5_AVAILABLE'

# run_httpfs_files <comma-separated files> <log>: run them as the httpfs pass does, and
# exit the script when one fails or skips for a reason this run should have provided.
run_httpfs_files() {
    local files="$1" log="$2" rc=0 start=$SECONDS elapsed bad_skips
    # An empty filter would select every test, not none.
    if [ -z "$files" ]; then
        echo "ERROR: no test files for $log" >&2
        exit 1
    fi
    echo "Running the $(echo "$files" | tr ',' '\n' | wc -l) test files that need httpfs (log: $log)..."
    # on_init's LOAD fails every file when httpfs cannot load (under autoloading, `require
    # httpfs` alone passes without loading it). --keep-home keeps this HOME, which DuckDB 2.0's
    # unittest otherwise swaps for an empty sandbox of its own. --require miint fails a file
    # whose miint cannot load, where it would otherwise be skipped; -w NoTests fails the run
    # when any file in the list matches no test.
    HOME="$HTTPFS_TEST_HOME" DUCKDB_TEST_AUTOLOADING=all DUCKDB_TEST_ON_INIT='LOAD httpfs;' \
        ${TIMEOUT_CMD:+"$TIMEOUT_CMD" --kill-after=60 "$HTTPFS_TIMEOUT_SECONDS"} \
        ./build/release/test/unittest --keep-home --require miint -w NoTests "$files" \
        > "$log" 2>&1 || rc=$?
    elapsed=$((SECONDS - start))
    if [ "$rc" -ne 0 ]; then
        # In full: CI keeps the console output, not the log file.
        cat "$log"
        # timeout exits 124, or 137 when it has to kill; 137 is also any other SIGKILL, such
        # as the out-of-memory killer's, so only a run that lasted the limit counts as timed out.
        if { [ "$rc" -eq 124 ] || [ "$rc" -eq 137 ]; } && [ "$elapsed" -ge "$HTTPFS_TIMEOUT_SECONDS" ]; then
            echo "ERROR: tests that need httpfs did not finish within ${HTTPFS_TIMEOUT_SECONDS}s" >&2
        else
            echo "ERROR: tests that need httpfs failed (exit $rc after ${elapsed}s)" >&2
        fi
        exit 1
    fi
    tail -n 20 "$log"
    bad_skips=$(sed -n '/^Skipped tests for the following reasons:/,$p' "$log" | tail -n +2 |
        grep -vE "^(require-env ($HTTPFS_ALLOWED_SKIPS): [0-9]+)?$" || true)
    if [ -n "$bad_skips" ]; then
        echo "ERROR: the httpfs pass skipped files for reasons this run should have provided:" >&2
        echo "$bad_skips" >&2
        exit 1
    fi
}

run_httpfs_files "$HTTPFS_TESTS" "$PWD/build/release/httpfs_tests.log"

# The live ENA tests, which the pass above skipped (see MIINT_ENA_LIVE_TEST above).
if [ "$ENA_LIVE" = "1" ]; then
    ENA_LIVE_TESTS=$(grep -rl --include='*.test' '^require-env ENA_AVAILABLE' test/sql | paste -sd, -)
    with_live_ena run_httpfs_files "$ENA_LIVE_TESTS" "$PWD/build/release/httpfs_ena_live_tests.log" || exit 1
fi

# Run shell script tests. A script that loads httpfs runs with the httpfs pass's HOME
# (docs/testing.md, "Shell Script Tests").
echo "Running shell script tests..."
for test_script in test/shell/*.sh; do
    if [ -f "$test_script" ]; then
        echo "Running $test_script..."
        script_home="$HOME"
        if grep -qi 'load httpfs' "$test_script"; then
            script_home="$HTTPFS_TEST_HOME"
        fi
        # A script that checks ENA_AVAILABLE is a live ENA test, and skips itself without it.
        live=""
        if [ "$ENA_LIVE" = "1" ] && grep -q 'ENA_AVAILABLE' "$test_script"; then
            live=with_live_ena
        fi
        if ! $live env HOME="$script_home" bash "$test_script"; then
            echo "Shell test failed: $test_script"
            exit 1
        fi
    fi
done
echo "All shell script tests passed!"
