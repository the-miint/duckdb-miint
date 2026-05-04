# Plan: ENA Submission Support — TDD Implementation

## Resuming this work (read this first)

If you're picking this up in a fresh session:

- **Repo**: `/home/dtmcdonald/dev/duckdb-miint-alt`. Branch: `feat/ena-submission` (cut from `v1.5-variegata`).
- **Status as of 2026-05-03**: Phases 0–6 complete and committed (Phase 6 main pending commit). **Phase 7 (`ena.experiments` + `ena.runs`) is next.**
  - `a8e8ec0` Phase 0 — branch + scaffolding (no behavior change)
  - `100a2c7` Phase 1 — `CREATE SECRET (TYPE ENA, ...)` with password redaction + literal/env/file indirection
  - `f88e756` Phase 2 — `ENAClient::PostJSON`/`PostXML` with HTTP Basic auth + retry on 429/5xx (incl. 504)
  - `bcd10fe` Phase 3 — `BuildEnvelopeJSON` (V2 wire format) + `ParseReceiptXML` (canonical XSD-governed shape).
  - `321045d` Phase 4 — `ena.projects` virtual table INSERT (StorageExtension scaffolding + projects insert operator + `ena.submission_log`)
  - `bb94ea4` Phase 5 — `ena.samples` virtual table INSERT — **MILESTONE 1 (BioSample-first workflow)**. Shared helpers in `ena_insert_common.{cpp,hpp}`.
  - `bb59c26` Phase 6 prep — `FastqEncoder` extracted from `copy_fastq.cpp` into a reusable class with a `void(const char *, size_t)` sink callback (DuckDB-free, unit-testable)
  - (pending) Phase 6 main — `ena_upload_reads` table function + Aspera write-side helper.
- **Test status**: 797 SQL cases / 6799 assertions / 0 failed; Catch2 binary 6698 assertions in 114 cases / 0 failed. Full `bash run_tests.sh` clean.
- **Architecture-of-record**: `localdocs/ena-submission-design-v2.md`. Open decisions in §11 are resolved (default endpoint=test, naming=`ena.projects`, persistent `ena.submission_log`, error on overwrite, flat upload layout, system CA store).
- **Research artifacts** in `localdocs/`:
  - `ena-research-webin-v2-deep.md` — Webin V2 endpoints, receipt structure, XSDs, auth (HTTP Basic only)
  - `ena-research-microbiome-metadata.md` — MIxS/checklists, taxon IDs, library/instrument enumerations
  - `ena-research-tools-and-auth.md` — webin-cli/ena-upload-cli architecture, DuckDB Secrets pattern, real-world failure modes
  - `ena-submission-design.md` (v1) — original 5-approach survey, still readable but superseded by v2
- **Conventions** (from CLAUDE.md and stored memory):
  - Never use `git add -A` or `git add .` — stage files explicitly by path.
  - Never `rm` without permission.
  - Pre-commit hook runs `make format-check`. Auto-fix via `conda run -n duckdb-143 make format-fix`.
  - Build: `bash build.sh` (initial). Incremental: `cmake --build build/release --config Release`. Test: `bash run_tests.sh` (full) or `./build/release/extension/miint/tests "[tag]"` (targeted Catch2).
  - When adding library deps, update `miint_versions()` in `src/miint_extension.cpp` lines 114–168.
- **Workflow contract**:
  - Strict red/green/refactor TDD per phase.
  - `linus-code-reviewer` runs at end of each phase's refactor; address blocking issues before approval gate.
  - **STOP after each phase's review fixes** and wait for user "approved, continue" before starting the next phase.
  - At end of each phase, append a "Phase N done" subsection to this plan (`## Implementation log` section) with files touched, test deltas, deviations, follow-ups.
  - User asked previously to commit each phase before starting the next (matches the existing per-phase commit cadence; see Phases 0–5 plus Phase 6 prep).

**Up next**: Phase 7 (`ena.experiments` + `ena.runs`). Phase 6 main is committed (see Implementation log). Pre-Linus plan additions, retained:
- **Phase 5 review surfaced a Phase-7-prep ask**: template the per-table `Finalize` body and extract a CRTP Sink/Source base from `ena_projects_insert_op.cpp` + `ena_samples_insert_op.cpp` before the third concrete operator (experiments) lands. Plan this as the first refactor in Phase 7's Red phase, not the last.
- **Reports API SELECT path** (`SELECT * FROM ena.{projects,samples}` via `/ena/submit/report/...`) was deferred from Phase 5. Pick this up after Phase 6 if no consumer surfaces sooner.
- **Transport-for-uploads decision (2026-05-03)**: Phase 6 ships **Aspera-only**. FTP/FTPS is parked — add it later only if a real user environment can't run `ascp`. This swaps the original Phase 6 (FTPS) and the original Phase 8 (Aspera write-side) — Phase 8 now drops the Aspera write-side bullet and gains nothing else from Phase 6's deferral.

---

## Context

duckdb-miint currently has read-only ENA support (`read_ena_*` table functions) but no way to submit data to the European Nucleotide Archive. Microbiome users routinely need to register projects (`PRJEB…`) and BioSamples (`SAMEA…`) **before** any sequencing data exists — for grant reporting, manuscript references, and consortium coordination — then later upload reads and register experiments/runs once data lands. The submission also needs to honor MIxS/GSC checklists for downstream MGnify analysis.

The architectural design is locked in `localdocs/ena-submission-design-v2.md` (hybrid approach: virtual tables for the SRA object graph + one table function for upload + scalar functions for lifecycle, backed by HTTP Basic auth via DuckDB Secrets Manager). All six open questions from §11 of that doc are resolved.

This plan implements that design phase-by-phase using strict red/green/refactor TDD, with a code review and explicit user approval gate at every phase boundary. The plan file itself is updated at end of each phase to reflect actual progress.

## Branch & workflow

- **Branch**: `feat/ena-submission` cut from `v1.5-variegata`.
- **Commits**: one logical commit per red/green/refactor step within a phase. Final phase commit merges via PR or is left for user-driven merge.
- **Code review**: after each phase's refactor step, invoke `linus-code-reviewer` agent over the phase's diff. Surface findings; address blocking issues before requesting approval.
- **Approval gate**: after review feedback is addressed, **STOP** and wait for user "approved, continue" before starting next phase.
- **Plan update**: at end of each phase, append a "Phase N done" subsection here recording actual files touched, test counts, deviations, and any follow-up issues. Don't start phase N+1 until that subsection is written.
- **Documentation**: deferred to phase 10 (single doc pass at the end), per user direction.

## TDD cadence per phase

1. **Red** — write the failing tests first (SQL `.test`, Catch2 C++, or both), run the suite, confirm new tests fail and existing tests pass.
2. **Green** — minimal implementation to pass. No premature abstraction. Existing tests stay green.
3. **Refactor** — extract duplication, tighten naming, remove dead code. Tests stay green throughout.

### Test pyramid (credential-aware)

EBI does not publish shared test credentials — each user registers a free Webin account at <https://www.ebi.ac.uk/ena/submit/webin/> (one-time, ~5 min). To minimize the burden, the test suite is layered so the **majority of tests need no credentials at all**:

1. **Unit (no creds, no network)** — Catch2 tests against canned fixtures: envelope builder byte-equality vs official ENA tutorial XML/JSON, receipt parser over a fixture corpus, FASTQ encoder bit-exactness, alias collision logic with a mock fetcher, checklist field validation. **Bulk of the suite.** Runs everywhere.
2. **Integration with mock servers (no creds)** — SQL `.test` files using:
   - Local HTTP mock with canned receipts for envelope→POST→parse round-trips. Validates that `INSERT INTO ena.samples` correctly invokes `BuildEnvelope`, hits an HTTP endpoint, and parses a fixture receipt — no real Webin account needed. Implemented for Phases 4–5 (`test/scripts/ena_webin_mock.py`, spawned by `run_tests.sh`, exports `ENA_WEBIN_MOCK_URL`).
   - **No upload-side mock**: Phase 6 uploads use Aspera, which has no local mock. The Aspera transport is exercised only in the credentialed live tier; the encode→gzip→MD5 pipeline is exercised via a `file://` test bypass (`test/sql/ena_upload_reads_local.test`) without any network. (If FTP gets resurrected later, `pyftpdlib` would back it.)
3. **Live integration (creds required)** — guarded with `require-env ENA_WEBIN_TEST_USER` / `require-env ENA_WEBIN_TEST_PASSWORD` (pattern from `test/sql/read_ena_sequences.test`). One file per phase, asserting the full round-trip against `wwwdev.ebi.ac.uk`. Skipped silently in CI without secrets; runs locally for the maintainer.

Maintainer setup: register a Webin account once, set the two env vars in shell init or store as GitHub Actions secrets. Contributors hacking on internals never need them.

## Reference: existing patterns to reuse

From the exploration reports:

- **Separate-connection pattern** (`src/reference_table_reader.cpp` lines 31–80, also `src/placement_table_reader.cpp` lines 53–76, plus `docs/internals/reading-tables-views.md`) — for `ena_upload_reads` reading user-supplied relations.
- **HTTP retry/rate-limit/backoff** (`src/ena_client.cpp` lines 11–65) — extend with POST + Basic auth.
- **Aspera binary discovery** (`src/aspera_utils.cpp` lines 35–82) — extend to support `ASPERA_SCP_PASS` env var for write-side.
- **FASTQ encoder** — extracted in commit `bb59c26` (Phase 6 prep) as `src/fastq_encoder.{cpp,hpp}` with a `void(const char *, size_t)` sink callback. `copy_fastq.cpp` now wraps it; `ena_upload_reads` will plug a gzip+MD5 sink in front of it.
- **Test idioms** (`test/sql/read_ena_sequences.test`, `test_ENAClient.cpp` mock fetcher pattern around lines 692–726) — for new ENA tests.
- **CMake source-list registration** — every new `.cpp` must be added to the explicit list in `CMakeLists.txt`.
- **Function registration** — every new SQL surface registered in `src/miint_extension.cpp::LoadInternal()`.
- **`miint_versions()` table function** (`src/miint_extension.cpp` lines 114–168) — update if any new library is embedded.

## Phases

### Phase 0 — Branch + scaffolding (no behavior change)

**Goal**: clean foundation; nothing user-visible.

**Red**: N/A (scaffolding only).

**Green**:
- Cut branch `feat/ena-submission` from `v1.5-variegata`.
- Add stubs to `CMakeLists.txt` source list (will fill in later phases): `src/ena_secret.cpp`, `src/ena_envelope_builder.cpp`, `src/ena_receipt_parser.cpp`, `src/ena_submission_log.cpp`, `src/ena_storage.cpp`, `src/ena_upload_reads.cpp`, `src/ena_lifecycle.cpp` — each with a single `// stub` comment.
- Extend `run_tests.sh` to:
  - Detect & export `ENA_WEBIN_TEST_USER` / `ENA_WEBIN_TEST_PASSWORD` if already set.
  - Probe `https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/swagger-ui/index.html` with a 3s timeout; if creds present and probe succeeds, export `ENA_WEBIN_TEST_AVAILABLE=1`.
  - ~~Spawn a local `pyftpdlib` server on a free port for the upload-mock test tier; export `ENA_FTP_MOCK_URL`.~~ **Dropped — Phase 6 went Aspera-only on 2026-05-03; no FTP mock is spawned. The Phase 0 done-log records that this server was deferred to its consumer, and the consumer no longer needs it. Resurrect only if FTP comes back.**
  - Spawn a local HTTP mock that returns canned ENA receipts on `/submit`; export `ENA_WEBIN_MOCK_URL`.

**Refactor**: N/A.

**Verification**: `bash build.sh && bash run_tests.sh` green; no new tests added; existing test count unchanged.

**Code review**: skip (scaffolding only).

**Approval gate**: STOP.

---

### Phase 1 — DuckDB ENA Secret type

**Goal**: `CREATE SECRET (TYPE ENA, ...)` works; passwords redacted in `duckdb_secrets()`.

**Red**:
- New `test/sql/ena_secret.test`: `CREATE SECRET test_ena (TYPE ENA, USER 'Webin-x', PASSWORD 'shh', ENDPOINT 'test')` followed by `SELECT name, type, secret_string FROM duckdb_secrets() WHERE name='test_ena'` — expects `password=redacted`.
- New `test/cpp/test_ena_secret.cpp`: Catch2 tests for `ENASecret::Resolve` covering: literal `password`, `password_env` indirection, `password_file` indirection, missing env var error, missing file error.

**Green**:
- `src/include/ena_secret.hpp` — `ENASecret` class wrapping `KeyValueSecret`.
- `src/ena_secret.cpp` — `RegisterENASecretType(ExtensionLoader&)` and `CreateENASecretFunction`. Pattern from `ena-research-tools-and-auth.md` §4 (mirroring duckdb-postgres). `redact_keys = {"password", "bearer_token"}`.
- Wire into `src/miint_extension.cpp::LoadInternal()`.
- Add `test/cpp/test_ena_secret.cpp` to `CMakeLists.txt`.

**Refactor**: extract password indirection helper (`ResolveIndirectSecret`) if used in more than one place.

**Verification**: `./build/release/test/unittest "test/sql/ena_secret.test"` passes; `./build/release/extension/miint/tests "[ena_secret]"` passes; `bash run_tests.sh` overall green.

**Code review**: yes, full diff.

**Approval gate**: STOP.

---

### Phase 2 — HTTP POST + Basic auth in `ENAClient`

**Goal**: `ENAClient::PostJSON(url, body, user, pw)` works with retry/backoff identical to existing GET path.

**Red**:
- `test/cpp/test_ENAClient.cpp` additions: `TEST_CASE("ENAClient POST")` covering success, 401, 403, 5xx → retry → success, network error, body MIME headers, Basic auth header construction. Use mock-fetcher pattern around lines 692–726 of the existing test file (extend the mock to record method + headers + body).

**Green**:
- Refactor `ENAClient::MakeRequest` into `MakeRequest(method, url, body, content_type, basic_auth_user, basic_auth_pw)` keeping the rate-limit + retry loop. Add `GetText`, `PostJSON` convenience methods. Use `duckdb::PostRequestInfo`.
- Update `src/include/ena_client.hpp`.

**Refactor**: any common request-setup code lifted out.

**Verification**: existing GET-based ENA tests still pass; new POST tests pass.

**Code review**: yes.

**Approval gate**: STOP.

---

### Phase 3 — JSON envelope builder + receipt parser + `ena.submission_log`

**Goal**: pure-data plumbing for the wire format. No live-server dependency.

**Red**:
- `test/cpp/test_ena_envelope.cpp` — given a struct describing one project / one sample / one experiment / one run, builds the V2 JSON envelope and asserts byte-equality against canned fixtures. Fixtures sourced from `localdocs/ena-research-webin-v2-deep.md` §1.
- `test/cpp/test_ena_receipt.cpp` — parses canned XML and JSON receipts (success, validation failure, partial-success edge cases) and asserts the parsed struct.
- `test/sql/ena_submission_log.test` — verifies `ena.submission_log` table exists with the schema from `ena-submission-design-v2.md` §7.7 (no rows yet — populating it requires phases 4+).

**Green**:
- `src/include/ena_envelope_builder.hpp` + `src/ena_envelope_builder.cpp` — `BuildEnvelope(Action, Vec<ProjectSpec>...)` etc. Use `yyjson` (already a DuckDB dependency).
- `src/include/ena_receipt_parser.hpp` + `src/ena_receipt_parser.cpp` — `ParseReceipt(content_type, body)` returns `ENAReceipt { success, accessions[], errors[], era_accession, ... }`. Both XML (via `expat`, already vendored) and JSON paths.
- `src/ena_submission_log.cpp` — declares the `ena.submission_log` table; populated as a regular DuckDB table (not virtual) created at first ENA-related connection. Append-only `LogSubmission(ENALogRow)` helper.

**Refactor**: builder + parser share enums/structs (`ENAObjectType`, `ENAAction`).

**Verification**: all new Catch2 tests pass; `ena_submission_log.test` passes.

**Code review**: yes — heavy on edge cases for the parser.

**Approval gate**: STOP.

---

### Phase 4 — `ena.projects` virtual table (INSERT only)

**Goal**: `INSERT INTO ena.projects (...) VALUES (...) RETURNING prjeb_accession, erp_accession` works against the ENA test endpoint.

**Red** (three tiers):
- Catch2: StorageExtension scaffolding (load with `ATTACH 'ena:test' AS ena`; `LookupSchema(...)` returns the predefined schema; `CreateTable` errors with "fixed schema").
- Catch2: `ENAProjectsInsert::Sink` + `Finalize` against the mock-fetcher pattern (extended from `test_ENAClient.cpp`'s mock around lines 692–726): given a one-row chunk, it builds the right envelope, makes the right POST, and populates the RETURNING vector with the canned receipt's accessions.
- `test/sql/ena_projects_insert_mock.test` — uses `ENA_WEBIN_MOCK_URL` from phase 0; full SQL round-trip against the local mock, no live creds required.
- `test/sql/ena_projects_insert_live.test` — `require-env ENA_WEBIN_TEST_USER`, `ENA_WEBIN_TEST_PASSWORD`, `ENA_WEBIN_TEST_AVAILABLE`. Same query against `wwwdev.ebi.ac.uk`, asserts `RETURNING` produces a non-null `prjeb_accession` matching `^PRJEB[0-9]+$`. Skipped in CI without secrets.

**Green** — this phase introduces the most novel infrastructure; expect highest LOC:
- Fetch reference impl from `github.com/duckdb/duckdb-postgres` (clone to scratch, read, don't link). Mirror the StorageExtension layout: `attach_function_t`, `Catalog`, `TransactionManager`, `Transaction`, `SchemaCatalogEntry`, `TableCatalogEntry`, `PhysicalInsert`.
- `src/include/ena_storage.hpp` — class declarations.
- `src/ena_storage.cpp` — `ENAStorageExtension`, `ENACatalog`, `ENATransactionManager`, `ENATransaction`, `ENASchemaEntry`. Schemas hard-coded (5 tables: projects, samples, experiments, runs, analyses), but only `projects` is wired in this phase; the other four return "not implemented yet" if poked.
- `src/ena_projects_insert.cpp` — `ENAProjectsInsert` PhysicalOperator subclass. Sink collects rows; Finalize builds envelope (phase 3), POSTs (phase 2), parses receipt (phase 3), populates RETURNING vector, writes to `ena.submission_log` (phase 3).
- Wire into `LoadInternal()`.

**Refactor**: factor out a generic `ENAObjectInsert` base class to make samples/experiments/runs trivially reuse it in later phases.

**Verification**: SQL test passes against `wwwdev.ebi.ac.uk` (skipped if creds absent). Existing tests green.

**Code review**: yes — biggest review of the project, focus on the StorageExtension correctness.

**Approval gate**: STOP.

---

### Phase 5 — `ena.samples` virtual table (BioSample-first milestone)

**Goal**: `INSERT INTO ena.samples (...) RETURNING ers_accession, samea_accession` works. **This is the user's headline requirement** — register samples and get BioSample accessions before any sequencing data exists.

**Red**:
- Catch2: `ENASamplesInsert` against the mock-fetcher pattern, plus a fixture corpus of receipts (success, single-sample, multi-sample, validation failure with `success=false` + `messages.error[]`).
- `test/sql/ena_samples_insert_mock.test` — three subtests against the local mock:
  - Insert a single sample with `taxon_id=408170, checklist='ERC000015'` and required attributes (collection_date, lat/lon, country, ENVO env_broad/local/medium).
  - Insert two samples in one statement; both get accessions.
  - Insert with missing checklist-mandatory attribute; expect a clear error from receipt parser.
- `test/sql/ena_samples_insert_live.test` — same shape, against `wwwdev.ebi.ac.uk`, gated on credentials.
- Catch2 test for `SELECT * FROM ena.samples` against canned Reports API fixture (read-back path).

**Green**:
- `src/ena_samples_insert.cpp` — `ENASamplesInsert` derived from the base extracted in phase 4's refactor.
- Wire `ena.samples` table-entry to its insert operator.
- Implement read-only `SELECT * FROM ena.samples`/`ena.projects` via the Reports API (`/ena/submit/report/...`).

**Refactor**: extract attribute MAP serialization helper if duplicated with phase 4.

**Verification**: SQL test passes; running phases 4+5 together demonstrates the end-to-end BioSample-first workflow:
```sql
CREATE SECRET (TYPE ENA, USER 'Webin-x', PASSWORD_ENV 'PW');
INSERT INTO ena.projects (alias, title, project_type) VALUES ('p1', 't', 'METAGENOMIC') RETURNING prjeb_accession;
INSERT INTO ena.samples (alias, taxon_id, checklist, attributes)
  SELECT 'sample_' || i, 408170, 'ERC000015', MAP { ... } FROM range(1,5) RETURNING samea_accession;
```

**Code review**: yes.

**Approval gate**: STOP. **MILESTONE 1 reached** (BioSample-first workflow).

---

### Phase 6 — `ena_upload_reads` table function (Aspera write-side)

**Goal**: stream a multi-sample reads relation directly to the Webin upload area via Aspera (`ascp --mode=send`); return per-file metadata (sample_ref, filename, filetype, md5, bytes_written, layout).

**Phase 6 prep — done, committed `bb59c26`**: FASTQ encoder extracted into `src/fastq_encoder.{cpp,hpp}` with a `void(const char *, size_t)` sink callback. 7 Catch2 cases / 9 assertions covering Phred+33 byte-exact output, comment formatting, paired-end interleaving, quality-overflow rejection (Phred+33 and +64), encoder reuse without state leak. Existing `test/sql/copy_fastq.test` (103 assertions) still green. `copy_fastq.cpp::FastqCopySink` now wraps the encoder with a sink lambda over `MemoryStream::WriteData` — no behaviour change.

**Transport decision (resolved 2026-05-03): Aspera only.**
- Uses miint's existing Aspera infrastructure: `src/aspera_utils.{cpp,hpp}` (binary discovery, key resolution) + `src/aspera_stream.cpp` (`fork()` + `pipe()` + `execvp("ascp", ...)` pattern, currently read-side / receive).
- Auth: `ASPERA_SCP_PASS` env var set in the child process before `execvp` so the password never lands on a command line. The user/password come from the existing ENA secret (already wired for Phase 4/5 INSERT paths).
- Performance: Aspera handles TLS, congestion control, and resume natively — no need to layer our own retry logic.
- FTP/FTPS is **parked**, not deferred to a numbered phase. Add it later only if a real user environment can't run `ascp`. The FTP-via-subprocess-`curl` design discussed earlier is preserved at the bottom of this section as a reference if needed.
- libcurl in-process is parked indefinitely (it's in `vcpkg.json` from earlier work but no consumer needs it; either prune from `vcpkg.json` later or keep for future use).

**Remaining Red**:
- `test/cpp/test_ena_upload_reads.cpp` — layout detection per sample (single / paired / mixed = error). Pure-data, no DuckDB linkage. Plus an Aspera-args-builder unit test (verify `--mode=send`, `-d` for directory creation, `-Q` QoS, `-l` rate limit, `-i` key path, target URL shape) — captures the argv vector without forking.
- `test/sql/ena_upload_reads.test` — `require-env ASPERA_AVAILABLE` (run_tests.sh detects `ascp` on PATH and exports this) plus `require-env ENA_WEBIN_TEST_USER` + `ENA_WEBIN_TEST_PASSWORD`. Skipped silently when either is absent. Uploads to the live `webin2.ebi.ac.uk` Aspera endpoint — there is **no local Aspera mock**, so this is a credentialed integration test. Write the test to be resumable (use a timestamped sample alias) so re-runs against the shared test account don't collide.
- `test/sql/ena_upload_reads_local.test` — **unguarded**, exercises the encode → gzip → MD5 pipeline with a `target_url := 'file:///tmp/...'` (or temp dir) bypass. This proves the upstream pipeline byte-by-byte without Aspera or the network. The Aspera transport is the only piece this test doesn't cover.

**Remaining Green**:
- `src/include/ena_upload_reads.hpp` + `src/ena_upload_reads.cpp` — table function. Named params: `relation` (required), `secret` (required), `target_url` (default `aspera://webin2.ebi.ac.uk/`; supports `file://` for the local test bypass), optional `qual_offset` (default 33), optional `layout` (default `auto`; values `auto|single|paired|paired_interleaved`), optional `aspera_rate_limit_mbps` (default unset, passes through to `ascp -l`).
  - **Bind**: validate input schema (`sample_ref VARCHAR`, `read_id VARCHAR`, `sequence1 VARCHAR`, `qual1 LIST(UTINYINT)`, `sequence2 VARCHAR` nullable, `qual2 LIST(UTINYINT)` nullable, `sequence_index BIGINT` nullable). Output schema: `sample_ref VARCHAR`, `filename VARCHAR`, `filetype VARCHAR`, `md5 VARCHAR`, `bytes_written UBIGINT`, `layout VARCHAR`.
  - **InitGlobal**: read the relation via the separate-connection pattern (`docs/internals/reading-tables-views.md`; `src/reference_table_reader.cpp` template). Group rows by `sample_ref`; sort each group by `sequence_index` (NULLS LAST → original row order). Detect layout per sample: all R2 NULL = `single`; all R2 non-NULL = `paired`; mixed = error with the alias and offending row. Resolve Aspera config once via `AsperaUtils::FindAscp()` + `AsperaUtils::ResolveKey()`; bail out if `target_url` is `aspera://` and `ascp` isn't on PATH.
  - **Execute**: claim a sample atomically (mutex + queue). Encode reads through `FastqEncoder` → gzip (zlib `deflate` with gzip header, `Z_DEFAULT_COMPRESSION`) → `duckdb::MD5Context::Add`. Two write paths:
    1. **Aspera**: write the gzipped FASTQ to a temp file under `/tmp/ena-upload-XXXXXX.fastq.gz` (ascp doesn't read stdin reliably per existing `src/aspera_stream.cpp` comment), then fork+execvp `ascp --mode=send --user=$WEBIN_USER -d -Q -i $key tempfile target_url:/path`. Set `ASPERA_SCP_PASS` in `envp[]` for the child. Wait, capture stderr, surface non-zero exit as `IOException`. Unlink the temp file on success and on error (best-effort).
    2. **`file://`**: write straight to the resolved local path. Same FASTQ byte stream, no transport.
  - On success emit one row per uploaded file with the MD5 the encoder→gzip→MD5 pipeline computed locally.
  - **Files emitted per layout**:
    - single: `{sample_ref}.fastq.gz`
    - paired-split (default): `{sample_ref}_1.fastq.gz`, `{sample_ref}_2.fastq.gz`
    - paired-interleaved (`layout := 'paired_interleaved'`): `{sample_ref}.fastq.gz` (one file, R1+R2 alternating)
- **Aspera write-side helper**: extend `src/aspera_stream.cpp` with a write/send process class, or add `src/aspera_send.cpp` that builds the argv from `(AsperaConfig, target_url, password, local_path)` and forks. Keep the read-side path untouched.
- `run_tests.sh` already detects `ASPERA_AVAILABLE` (line ~10–28 area, mirrors `BOWTIE2_AVAILABLE`). Verify; add `ENA_WEBIN_TEST_USER` / `ENA_WEBIN_TEST_PASSWORD` env passthrough if not already present. No new test server to spawn.
- Wire the function into `LoadInternal()` in `src/miint_extension.cpp`. No `miint_versions()` change (Aspera is a runtime binary, not a linked library).

**Refactor**: per-sample workspace tightening (one zlib stream + one MD5 context + one encoder per sample; reuse buffers across reads in the same sample).

**Verification**:
- `./build/release/extension/miint/tests "[fastq_encoder]"` — already passes (committed in `bb59c26`).
- `./build/release/extension/miint/tests "[ena_upload_reads]"` — new layout-detection + ascp-args-builder coverage.
- `./build/release/test/unittest "test/sql/ena_upload_reads_local.test"` — local file:// round-trip; exercises encoder → gzip → MD5 deterministically without any external dependency.
- `ASPERA_AVAILABLE=1 ENA_WEBIN_TEST_USER=... ENA_WEBIN_TEST_PASSWORD=... ./build/release/test/unittest "test/sql/ena_upload_reads.test"` — live Aspera round-trip when creds + `ascp` both present (skipped otherwise).
- Existing `test/sql/copy_fastq.test` and `test/sql/ena_storage_attach.test` etc. unchanged.

**Code review**: yes — table function correctness, separate-connection pattern correctness, fork/pipe error handling, layout-detection edge cases, temp-file cleanup on error path, `ASPERA_SCP_PASS` not leaking via `getenv` reads after the child exits.

**Approval gate**: STOP.

**Out of scope for Phase 6** (parked):
- **FTP/FTPS upload via subprocess `curl`** — keep the design (in this paragraph) on file in case real users hit Aspera-blocked environments. The plan: `fork()` + `pipe()` + `execvp("curl", "--upload-file", "-", "--user", "X:Y", "ftp://webin2.ebi.ac.uk/path")`, password via env not argv. The encoder→gzip→MD5 pipeline doesn't change; only the transport-side write target swaps. `pyftpdlib` would back the SQL test mock. **Add only when there's a documented user request.**
- **In-process libcurl** — only worth the build-system churn if Aspera proves inadequate **and** `curl` subprocess proves inadequate.
- **Resume / partial-upload retry** — Aspera handles resume natively when configured; we don't add anything until a user reports a failure mode.

---

### Phase 7 — `ena.experiments` + `ena.runs` (full submission milestone)

**Goal**: complete the project → samples → upload → experiments → runs pipeline end-to-end.

**Red**:
- `test/sql/ena_full_pipeline.test` — guarded by Webin test creds + `ENA_WEBIN_TEST_AVAILABLE`. Builds 2 projects × 4 samples × 4 runs in one test, asserts all server-assigned accessions are returned and stored in `ena.submission_log`.
- Catch2 tests for envelope builder additions (experiments, runs).

**Green**:
- `src/ena_experiments_insert.cpp`, `src/ena_runs_insert.cpp` — derive from phase 4's base. Library/platform/instrument enumerations validated client-side (compile-time constants from `enasequence/webin-xml`).
- Wire `ena.experiments`, `ena.runs` tables.

**Refactor**: nothing major expected — base class should absorb most of it.

**Verification**: full-pipeline test passes against `wwwdev.ebi.ac.uk`.

**Code review**: yes.

**Approval gate**: STOP. **MILESTONE 2 reached** (full submission pipeline).

---

### Phase 8 — Alias collision check + checklist registry

**Goal**: production-quality polish before lifecycle work. (The original "Aspera write-side" bullet moved to Phase 6 when uploads went Aspera-only on 2026-05-03.)

**Red**:
- `test/cpp/test_ena_alias_check.cpp` — pre-INSERT alias lookup against a mock portal-API fetcher.
- `test/cpp/test_ena_checklist.cpp` — parse one bundled checklist XML (e.g. ERC000015) and assert mandatory fields enumerated correctly.

**Green**:
- Alias collision: `ENAClient::PortalSearch(object_type, alias)` reusing `MakeRequest`. Called pre-INSERT in each `ENA*Insert::PreSubmit`; throws on collision unless `if_not_exists` mode.
- Checklist registry: at extension load, parse bundled XMLs (`ena/checklists/ERC000NN.xml`) into static map. Insert path validates user-supplied attribute MAP keys against the chosen checklist's mandatory + optional fields.
- Update `miint_versions()` to include the bundled checklist commit hash.

**Refactor**: any duplication between the two sub-features.

**Verification**: all new tests pass; existing tests still pass.

**Code review**: yes.

**Approval gate**: STOP.

---

### Phase 9 — Lifecycle scalars (validate, hold, release, modify, cancel)

**Goal**: complete the SQL surface from the design doc §7.3.

**Red**:
- `test/sql/ena_lifecycle.test` — one subtest per scalar against `wwwdev`. `ena_validate(...)` against a deliberately invalid sample; expects `success=false` row. `ena_hold(accession, '2027-01-01')` then `SELECT hold_until_date FROM ena.samples WHERE alias=...`. `ena_release` then re-query. `ena_modify_sample(accession, attrs)` then re-query. `ena_cancel(accession)` then expect status='CANCELLED'.

**Green**:
- `src/ena_lifecycle.cpp` — five scalar / table functions. All thin wrappers over `BuildEnvelope` + `PostJSON` + `ParseReceipt` + log.

**Refactor**: nothing major.

**Verification**: lifecycle test passes.

**Code review**: yes.

**Approval gate**: STOP.

---

### Phase 10 — Documentation

**Goal**: reference docs match shipped behavior; design doc gets a "as-built" addendum.

**Updates** (no test changes):
- `docs/table-functions.md` — add `ena_upload_reads`, `read_ena_*` if changed.
- `docs/scalar-functions.md` — add `ena_validate`, `ena_hold`, `ena_release`, `ena_modify_*`, `ena_cancel`.
- `docs/ena-submission.md` — new file covering the canonical workflow with worked examples (mirrors §3 of the design v2 doc, with copy-pasteable SQL).
- `docs/internals/embedded-tools.md` — note checklist XML bundle.
- `docs/installation.md` — mention `ENA_WEBIN_TEST_USER`/`PASSWORD` env vars for users running the test suite locally.
- `localdocs/ena-submission-design-v2.md` — append "Implementation status" section linking to commit shas.
- `README.md` — one-line mention.

**Verification**: docs render in the doc viewer; no broken links; CI doc-link checker (if any) green.

**Code review**: skim only.

**Approval gate**: STOP. **PROJECT COMPLETE**.

---

## Critical files (modified or created)

Created:
- `src/ena_secret.cpp` + `.hpp`
- `src/ena_envelope_builder.cpp` + `.hpp`
- `src/ena_receipt_parser.cpp` + `.hpp`
- `src/ena_submission_log.cpp` + `.hpp`
- `src/ena_storage.cpp` + `.hpp`
- `src/ena_projects_insert.cpp`
- `src/ena_samples_insert.cpp`
- `src/ena_experiments_insert.cpp`
- `src/ena_runs_insert.cpp`
- `src/ena_upload_reads.cpp` + `.hpp`
- `src/ena_lifecycle.cpp` + `.hpp`
- `src/fastq_encoder.cpp` + `.hpp` (refactored out of `copy_fastq.cpp`)
- `ena/checklists/ERC000{11,13–25,28,29}.xml` (bundled)
- `test/sql/ena_secret.test`, `ena_projects_insert.test`, `ena_samples_insert.test`, `ena_submission_log.test`, `ena_upload_reads.test`, `ena_upload_reads_aspera.test`, `ena_full_pipeline.test`, `ena_lifecycle.test`
- `test/cpp/test_ena_secret.cpp`, `test_ena_envelope.cpp`, `test_ena_receipt.cpp`, `test_fastq_encoder.cpp`, `test_ena_upload_reads.cpp`, `test_ena_alias_check.cpp`, `test_ena_checklist.cpp`
- `docs/ena-submission.md`

Modified:
- `src/ena_client.cpp` + `.hpp` (POST + Basic auth)
- `src/aspera_utils.cpp` + `.hpp`, `src/aspera_stream.cpp` + `.hpp` (write-side support)
- `src/copy_fastq.cpp` (extract encoder)
- `src/miint_extension.cpp` (registrations, `miint_versions()`)
- `CMakeLists.txt` (source list)
- `run_tests.sh` (env guards — Aspera + Webin creds; HTTP receipt mock spawn; no FTP test server now that Phase 6 is Aspera-only)
- `test/cpp/test_ENAClient.cpp` (POST tests)
- `docs/table-functions.md`, `docs/scalar-functions.md`, `docs/internals/embedded-tools.md`, `docs/installation.md`, `README.md`
- `localdocs/ena-submission-design-v2.md` (as-built addendum)

## End-to-end verification (after Phase 10)

```bash
# Setup
git checkout feat/ena-submission
bash build.sh

# Smoke: existing tests still pass
bash run_tests.sh

# Smoke: ENA-specific (requires test creds)
export ENA_WEBIN_TEST_USER=Webin-XXXXX
export ENA_WEBIN_TEST_PASSWORD=...
bash run_tests.sh

# Manual canonical workflow
./build/release/duckdb <<EOF
CREATE SECRET my_ena (TYPE ENA, USER 'Webin-XXXXX', PASSWORD_ENV 'ENA_WEBIN_TEST_PASSWORD');
INSERT INTO ena.projects (alias, title, project_type)
  VALUES ('demo-' || strftime(now(), '%Y%m%d-%H%M%S'), 'Demo project', 'METAGENOMIC')
  RETURNING prjeb_accession, erp_accession;
INSERT INTO ena.samples (alias, taxon_id, checklist, attributes)
  VALUES ('demo-sample-1', 408170, 'ERC000015',
    MAP {'collection date': '2026-05-01',
         'geographic location (country and/or sea)': 'United States',
         'broad-scale environmental context': 'human-associated habitat [ENVO:00009003]',
         ...})
  RETURNING ers_accession, samea_accession;
SELECT * FROM ena.submission_log;
EOF
```

---

## Implementation log

(Updated at the end of each phase. Each entry: files actually touched, test count delta, deviations from plan, follow-up issues.)

### Phase 0 — done (commit a8e8ec0)

Files touched (9): `CMakeLists.txt`, `run_tests.sh`, plus 7 new stub `.cpp` files (`src/ena_secret.cpp`, `src/ena_envelope_builder.cpp`, `src/ena_receipt_parser.cpp`, `src/ena_submission_log.cpp`, `src/ena_storage.cpp`, `src/ena_upload_reads.cpp`, `src/ena_lifecycle.cpp`).

Tests delta: 0 new, 0 regressions. C++ tests `712 cases | 6636 assertions | 0 failures`. `bash run_tests.sh` end-to-end green.

Deviations from plan:
- **Mock servers (pyftpdlib FTP, HTTP receipt mock) deferred** to the phases that consume them (Phase 6 and Phase 4 respectively). Setting them up now would be dead infrastructure, and the canned receipts they return depend on fixtures that come from Phase 3.
- Skipped header (`.hpp`) stubs — added them as needed when the corresponding `.cpp` is implemented. CMake source list only enumerates `.cpp`.

Follow-ups: none.
### Phase 1 — done (commit 100a2c7)

Files touched (6): `CMakeLists.txt` (+1), `src/ena_secret.cpp` (+85), `src/miint_extension.cpp` (+2), `src/include/ena_secret.hpp` (NEW, 72 lines), `test/cpp/test_ena_secret.cpp` (NEW, 116 lines, 10 Catch2 cases), `test/sql/ena_secret.test` (NEW, 99 lines, 21 assertions).

Tests delta: +10 C++ cases (722 total), +1 SQL test file (21 assertions). 0 regressions.

Deviations from plan:
- **`ResolvePasswordIndirection` lives inline in the header**, not in `ena_secret.cpp`. Reason: miint's existing test executable doesn't link against duckdb (see `ena_parser.cpp` using `std::runtime_error`), and pulling the helper out as inline let it stay unit-testable. The DuckDB-glue `CreateSecret` callback wraps `runtime_error → InvalidInputException`.
- **`bearer_token` parameter dropped** entirely (was in design doc §6 as "deferred post-MVP"). Reviewer flagged it as a ghost parameter — registered but never stored. Re-add it in a future phase only if Webin V2 adds Bearer auth (currently Basic only per `ena-research-webin-v2-deep.md` §3).
- **Refactor step skipped** (no obvious extraction; `GetOption` is a 6-line local helper used twice in one file).

Code review (linus-code-reviewer): one blocking issue (bearer_token ghost parameter — fixed), two nits (comment accuracy + unnecessary `std::move` — both fixed). Verdict cleared after fixes.

Follow-ups: none.
### Phase 2 — done (commit f88e756)

Files touched (3): `src/ena_client.cpp` (+59), `src/include/ena_client.hpp` (+93/-5), `test/cpp/test_ENAClient.cpp` (+89). Total ~236 lines added.

Tests delta: +9 C++ cases (731 total), +27 assertions (6675 total). 0 regressions.

Deviations from plan:
- **Round-trip POST tests deferred to Phase 4** as anticipated. The plan said "extend the mock to record method + headers + body" but the existing MockFetcher (`test_ENAClient.cpp:692–726`) is a `string -> string` function decoupled from the HTTP layer; it does not represent the API surface for HTTP-level concerns. Phase 4's mock server is the right place.
- **All retryable behaviour is unit-tested via `IsRetryableStatus`** (extended to 504 per reviewer) rather than running through the retry loop. The retry loop is duplicated between GET and POST; refactoring to a generic helper deferred until Phase 3 lands more call sites.
- **`Bearer` auth still not exposed** — Webin V2 OpenAPI is Basic-only.

Code review (linus-code-reviewer):
- **Blocking — Accept/Content-Type coupling**: original code unconditionally set `Accept = Content-Type`. Refactored `PostBody` to take an explicit `accept_type` parameter (mismatch allowed), with `PostJSON`/`PostXML` defaulting to matching values.
- **Should-fix**: 504 added to retryable status set. Colon-in-userid validation in `BuildBasicAuthHeader` (RFC 7617 §2.1) added. Empty-string Base64 test added.
- **Nit**: dead-code post-loop throw documented as "// unreachable". Comments added explaining `PostRequestInfo::buffer_out` semantics + rebuild-per-iteration rationale.

Follow-ups (deferred):
- DRY the GET and POST retry loops in Phase 3 (will then be three call sites, the third tipping the cost-benefit toward extraction).
- Consider scrubbing credential `std::string` allocations after use if/when miint runs in a multi-tenant context.
### Phase 3 — done (commit bcd10fe)

Files touched (7): `CMakeLists.txt` (+4), `src/ena_envelope_builder.cpp` (+201), `src/ena_receipt_parser.cpp` (+156), `src/include/ena_envelope_builder.hpp` (NEW, 61 lines), `src/include/ena_receipt_parser.hpp` (NEW, 47 lines), `test/cpp/test_ena_envelope.cpp` (NEW, 219 lines, 14 cases / 31 assertions), `test/cpp/test_ena_receipt.cpp` (NEW, 173 lines, 6 cases / 53 assertions). Total ~859 lines added.

Tests delta: +20 C++ cases (751 total), +84 assertions (6759 total). 0 regressions.

Deviations from plan:
- **`ena.submission_log` table deferred to Phase 4.** The plan placed log scaffolding here, but the table has no consumer until Phase 4 wires `INSERT INTO ena.projects`. Land it with its first consumer rather than as an empty placeholder.
- **JSON receipt parsing deferred** — only XML supported. Phase 4+ will POST with `Accept: application/xml` so all receipts arrive in the canonical XSD-governed shape. JSON parser can be added if anyone needs it.
- **`SampleSpec.scientific_name` not exposed.** Was scaffolded then removed during refactor (dead code; no consumer until Phase 5+). Will reintroduce when needed.

Code review (linus-code-reviewer):
- **Blocking — HOLD invariant**: header documented validation (HOLD requires date) that the code didn't enforce, plus a silent double-HOLD bug when both `action=HOLD` and `hold_until_date` were set. Now: throws on either misconfiguration, with two new tests pinning the behavior.
- **Blocking — dead `elem_stack`** in receipt parser: removed; replaced with comment explaining why a single boolean suffices for the current XSD.
- **Should-fix — `int` truncation on >2GB receipts**: bounds check added before `XML_Parse`.
- **Should-fix — `ActionName` silent fallback**: replaced `return "ADD"` with `throw std::logic_error("unhandled ENAAction")` so a future enum addition fails loudly at runtime instead of silently emitting ADD.
- **Should-fix — UTF-8 precondition**: documented in header (caller responsible for valid UTF-8; no validation/U+FFFD substitution).
- **Nits not addressed**: `AppendArray` `needs_comma` naming, `IsObjectElement` linear scan over 11 strings — both judged not worth churn for current scale.

Follow-ups: none.
### Phase 4 — done (commit `321045d`)

Files touched (15):
- New code (8): `src/include/ena_storage.hpp` (227), `src/ena_storage.cpp` (~590), `src/include/ena_projects_insert.hpp` (70), `src/ena_projects_insert.cpp` (~115), `src/include/ena_projects_insert_op.hpp` (58), `src/ena_projects_insert_op.cpp` (~395)
- New tests (4): `test/cpp/test_ena_projects_insert.cpp` (190 lines, 6 cases / 25 assertions), `test/sql/ena_storage_attach.test` (47 assertions), `test/sql/ena_projects_insert_mock.test` (skipped without httpfs in test harness; consistent with `read_ncbi.test`), `test/sql/ena_projects_insert_live.test` (skipped without ENA_WEBIN_TEST_USER/PASSWORD/AVAILABLE)
- New mock infrastructure (1): `test/scripts/ena_webin_mock.py` (180 lines) — Python HTTPServer that returns canned XML receipts; spawned by `run_tests.sh`, exports `ENA_WEBIN_MOCK_URL`.
- Modified (5): `CMakeLists.txt` (+source list entries), `src/miint_extension.cpp` (registers `ENAStorageExtension` for `ena` storage type), `src/ena_secret.cpp` (adds optional `ENDPOINT_URL` parameter for mock/test override), `run_tests.sh` (mock server lifecycle), test SQL files.

Tests delta: +6 C++ cases (757 total), +25 assertions in those new cases. SQL tests: +1 file passing (`ena_storage_attach.test` = 47 assertions). Mock + live tests skip when their preconditions aren't met. Final: `bash run_tests.sh` green — 757 cases / 6684 assertions / 0 failed for the SQL tier; 113 cases / 6626 assertions / 0 failed for the C++ tier.

Manual round-trip (mock server + duckdb shell with `LOAD httpfs`) verifies:
- `INSERT INTO ena.projects (alias, title, project_type) VALUES (...) RETURNING ...` — server-assigned `prjeb_accession`/`erp_accession` round-trip correctly.
- `SELECT * FROM ena.submission_log` — both successful and failed submissions are logged, with raw receipt preserved on logical failure.

Deviations from plan:
- **Mock SQL test infrastructure landed here, not Phase 0.** The plan deferred mock-server setup, and indeed the receipt fixture only existed once Phase 3 landed; building it now means it has a real consumer. Skipped automatically when `httpfs` isn't pre-installed in the test harness (POST is implemented by httpfs, not the built-in DuckDB http client).
- **`umbrella_children` and `attributes` columns reject non-empty values with NotImplementedException** rather than being silently honoured. Acceptable for Phase 4; revisit when needed.
- **`hold_until_date` per-row agreement requirement**: ENA Webin V2 applies HOLD at the submission level, not per-object. We surface a clear error if rows in one INSERT disagree. Documented as known limitation.
- **`submission_id` column type is VARCHAR not UUID.** UUID's physical storage is `hugeint_t`, but our row already carries the string form. Going through UUID would mean converting on every emit for no user-visible gain. Reverted from initial UUID design after a debug-build assertion fired during Linus review verification.
- **Catalog refactor for shared submission_log column emission**: extracted `AddSubmissionLogColumns(ColumnList&)` so the catalog and the scan function share a single source of truth — caught during Linus review.
- **`ena_projects_insert.cpp` SubmitProjectInsertOutcome no longer throws on logical failure.** Server-reported `success=false` and parse failures populate `outcome.success=false` + `outcome.error_messages` while preserving `raw_receipt`; the Finalize side throws after logging so the submission_log row keeps the server response on failed XML parses.

Code review (linus-code-reviewer):
- Two raised "blockers" were misdiagnoses (RETURNING column count and `sink_state` null guard) — both stem from the reference `PhysicalInsert` in DuckDB which has the same patterns; verified `op.types = table.GetTypes()` when `return_chunk=true` regardless of what columns the user RETURNs (the executor wraps with PhysicalProjection). Replaced silent guards with `D_ASSERT` for clarity.
- Should-fix #4 (secret lookup): switched to `GetSecretByName` without storage filter so secrets in any backend resolve.
- Should-fix #5 (`dynamic_cast` of `KeyValueSecret`): replaced with pointer cast + null check returning a clean `InvalidInputException`.
- Should-fix #6 (raw_receipt loss on parse failure): refactored `SubmitProjectInsertOutcome` to never throw on logical failures; the wrapper `SubmitProjectInsert` still throws to preserve unit-test contract.
- Should-fix #8 (submission_log schema duplication): factored `AddSubmissionLogColumns(ColumnList&)` shared between catalog `BuildENATableInfo` and scan-time `EmitSubmissionLogColumns`.
- Should-fix #10 (NotImplementedScanFn dead code): clarified with comments + a defensive throw.
- Nit N1 (unused `emitted` variable): removed.
- Nit N2 (O(n) alias lookup): replaced with positional index since `outcome.rows` is built in `projects` order.
- Other nits (mock auth tightening, trap accumulation, content_type parameter): judged not worth churn for current scale.
- Issue #1 (column_index_map bounds): added an explicit bounds check returning `INVALID_INDEX` for undersized maps.

Bug found by manual round-trip during review fixes: the `ena.submission_log` SCAN path was untested by the existing attach SQL test (which only verifies schema). Original UUID column type triggered a debug-build assertion when the scan called `StringVector::AddString` on a UUID-typed Vector. Fixed by switching to VARCHAR.

Follow-ups (deferred):
- Mock SQL test currently requires `httpfs` to be auto-installable in the test harness (matches existing `read_ncbi.test` pattern). To run locally, either build with `ENABLE_EXTENSION_AUTOINSTALL=1` or pre-install httpfs in the active extension cache.
- Phase 5 should reuse the StorageExtension scaffolding (catalog, transaction manager, schema) without modification; only the per-table insert operator should grow. The `ENAProjectsInsert` class was deliberately specific rather than generic — extracting a base class is deferred to Phase 5 once the second concrete user (samples) lands.
- `umbrella_children` / `attributes` MAP serialization to be wired in Phase 5 (samples need the MAP path) and back-ported to projects.

### Phase 5 — done (commit `bb94ea4`, on top of Phase 4 commit `321045d`)

**MILESTONE 1 reached** — BioSample-first workflow works end-to-end (ATTACH → CREATE PROJECT → CREATE SAMPLES with ERS/SAMEA accessions returned).

Files touched (15):
- New code (6): `src/include/ena_samples_insert.hpp` (50), `src/ena_samples_insert.cpp` (~115), `src/include/ena_samples_insert_op.hpp` (60), `src/ena_samples_insert_op.cpp` (~310), `src/include/ena_insert_common.hpp` (70), `src/ena_insert_common.cpp` (~115), `src/include/ena_post_fn.hpp` (20).
- New tests (3): `test/cpp/test_ena_samples_insert.cpp` (200 lines, 7 cases / 33 assertions), `test/sql/ena_samples_insert_mock.test`, `test/sql/ena_samples_insert_live.test`.
- Modified (5): `src/include/ena_envelope_builder.hpp` (added `description`, `scientific_name` to SampleSpec), `src/ena_envelope_builder.cpp` (emits both fields), `src/include/ena_storage.hpp` (`ENASubmissionLogRow.duration_ms` int32→int64), `src/ena_storage.cpp` (PlanInsert dispatcher routes to ENAProjectsInsert OR ENASamplesInsert; `submission_log.duration_ms` schema column INTEGER→BIGINT), `src/ena_projects_insert_op.cpp` (refactored to consume shared helpers in `ena_insert_common.{hpp,cpp}`; ResolveInputColumn calls hoisted out of row loop), `src/ena_projects_insert.cpp` (case-insensitive ext-id type comparison via shared `EqualsIgnoreCase`), `src/include/ena_receipt_parser.hpp` (added `EqualsIgnoreCase`), `src/include/ena_projects_insert.hpp` (delegated `ENAPostFn` to `ena_post_fn.hpp`), `src/ena_insert_common.cpp` (RecordSubmissionLog stops truncating duration_ms to int32), `test/scripts/ena_webin_mock.py` (SAMPLE accession layout fixed: ERS primary, SAMEA EXT_ID).
- CMakeLists.txt: source list updates for the 4 new .cpp files in EXTENSION_SOURCES + 2 in TEST_SOURCES.

Tests delta: +7 C++ cases (764 total), +33 assertions (6696 total). +2 SQL test files (skipped without httpfs/creds, matching read_ncbi.test pattern). 0 regressions. `bash run_tests.sh`: 764 SQL cases / 6717 assertions / 0 failed.

Manual round-trip via mock + duckdb shell with `LOAD httpfs`:
- `INSERT INTO ena.samples (alias, taxon_id, checklist, attributes) VALUES (...) RETURNING alias, ers_accession, samea_accession` — server-assigned ERS/SAMEA round-trip; user-supplied attribute MAP echoed in submission_log.request_payload.
- Multi-sample insert in one statement — both get distinct accessions.
- Failure path (alias contains FAIL) — InvalidInputException with submission_log.success=false + raw_receipt preserved.
- Invariant rejection — `taxon_id <= 0` and empty `alias` rejected client-side before any HTTP traffic.

Deviations from plan:
- **Reports API SELECT path (`SELECT * FROM ena.{projects,samples}` via `/ena/submit/report/...`) deferred** to a future phase. The BioSample-first INSERT workflow doesn't need read-back, the GET endpoint requires a different mock-server harness than the submit-side, and the design doc §7.4 read-back doesn't have a Phase 5 consumer. Currently `SELECT * FROM ena.samples` throws `NotImplementedException` via the existing scan stub.
- **`ena.samples` lacks a `hold_until_date` column** (intentional — ENA Webin V2 applies HOLD at the submission level, not per-object; the projects table got a per-row column for ergonomics, but samples don't have a clear analogue without a session-scoped option).
- **`attributes` column emits NULL on RETURNING** for samples (same trade-off as Phase 4 projects: user-supplied MAP is preserved verbatim in `ena.submission_log.request_payload`).
- **Shared base class for the operators NOT extracted yet.** Pulled out free helpers (ResolveENACredentials, ResolveInputColumn, ValueToVarchar, ValueToDateString, GenerateSubmissionId, RecordSubmissionLog) into `ena_insert_common.{hpp,cpp}` but the Sink/Source machinery still duplicates ~250 lines per operator. Linus review flagged this — see follow-up.
- **`SampleSpec` re-introduced** with `description` and `scientific_name` (Phase 3 had scaffolded then removed them as dead code).
- **`taxon_id` schema kept as INTEGER** — NCBI taxon IDs trivially fit in int32 today (max ~3M; INT32_MAX is 2.1B); changing to BIGINT adds nothing.

Code review (linus-code-reviewer):
- **Blocking — dead `attribute_struct_values` loop** in samples `AppendReturningRows`: built `Value::STRUCT` objects for every row that were silently discarded with `(void)`. Fixed: deleted the loop and the redundant first SetValue, leaving the single NULL-emitting line.
- **Blocking — case-sensitive `ext.type == "biosample"` comparison** would silently return empty `samea_accession` if EBI ever returned `"BioSample"` or `"BIOSAMPLE"`. Fixed: added `EqualsIgnoreCase` in `ena_receipt_parser.hpp` (inline, pure-data layer can use it without DuckDB linkage). Applied to both samples (biosample) and projects (study) ext-id matching.
- **Should-fix — `duration_ms` int32 truncation**: int32 caps at ~35 minutes; large submissions over saturated wwwdev could plausibly exceed it. Fixed: row struct + schema both bumped to int64/BIGINT.
- **Should-fix — `ena_samples_insert.hpp` pulled in `ena_projects_insert.hpp`** transitively for `ENAPostFn`. Fixed: extracted to `src/include/ena_post_fn.hpp` (tiny dedicated header).
- **Should-fix — `ResolveInputColumn` called per row** in both projects and samples build loops. Fixed: hoisted optional column resolutions before the chunk loop.
- **Should-fix — misleading HOLD comment** in `BuildSamplesFromBuffer` (architecture discussion in a code comment with no logic). Fixed: trimmed to one line stating the Phase-5 limitation.
- Nits not addressed (mock auth permissiveness, alias-collision risk in the live test, namespace string conventions): judged not worth churn for current scale.

Linus follow-up flagged for Phase 7: **template the `Finalize` body** so adding experiments/runs/analyses doesn't copy 70+ more lines per operator. The right shape is a function template parameterised on `SpecT`/`OutcomeT`/the submit functor / build-from-buffer / append-returning callbacks; Sink state and Source state can also be a CRTP base or a single helper struct. Do this before Phase 7's third concrete operator lands, not after.

Follow-ups (deferred):
- Reports API `SELECT * FROM ena.{projects,samples}` (Phase 6+).
- Real MAP-value emission in RETURNING for `attributes` (Phase 6+).
- Linus follow-up: extract the templated Finalize / CRTP Sink-Source base before Phase 7.
- HOLD on samples (session option or per-statement option), if anyone asks.

### Phase 6 — done (encoder commit `bb59c26` + main commit pending)

**Phase 6 prep — done (commit `bb59c26`)**:
- New: `src/include/fastq_encoder.hpp` (50 lines), `src/fastq_encoder.cpp` (40 lines), `test/cpp/test_fastq_encoder.cpp` (160 lines, 7 cases / 9 assertions).
- Modified: `src/copy_fastq.cpp` (replaces static `WriteFastqRecordToBuffer` with a `FastqEncoder` instance + sink lambdas wrapping `MemoryStream::WriteData`), `CMakeLists.txt` (added encoder + test to source lists).
- Tests delta: +7 C++ cases (771 total). 0 regressions; existing `test/sql/copy_fastq.test` (103 assertions) unchanged.
- Design choice: encoder takes a `void(const char *, size_t)` sink callback rather than `MemoryStream &` — keeps the encoder DuckDB-free so the unit-test target can link the .cpp without `libduckdb` (matches the existing `test/cpp/test_ena_envelope.cpp` pattern). The COPY hot path wraps `WriteData` in a one-line lambda; same number of writes per record, no extra allocation.

**Phase 6 main — done** (Aspera-only after the 2026-05-03 transport decision; FTP/FTPS parked).

Files touched (10):
- New code (5): `src/include/ena_upload_helpers.hpp` (99 lines), `src/ena_upload_helpers.cpp` (247 lines), `src/include/ena_upload_reads.hpp` (26 lines), `src/ena_upload_reads.cpp` (780 lines), `src/include/aspera_send.hpp` (38 lines), `src/aspera_send.cpp` (126 lines).
- New tests (3): `test/cpp/test_ena_upload_reads.cpp` (306 lines, 26 cases / 75 assertions), `test/sql/ena_upload_reads_local.test` (149 lines, unguarded file:// round-trip, 72 assertions), `test/sql/ena_upload_reads.test` (44 lines, gated live Aspera test, skipped without creds).
- Modified (2): `CMakeLists.txt` (`src/ena_upload_helpers.cpp`, `src/ena_upload_reads.cpp`, `src/aspera_send.cpp` added to EXTENSION_SOURCES; helpers + test added to TEST_SOURCES), `src/miint_extension.cpp` (`#include <ena_upload_reads.hpp>` + `ENAUploadReadsTableFunction::Register(loader)`).

Tests delta: +26 C++ cases (797 SQL, 6799 SQL assertions, 0 failed; 114 C++ test cases, 6698 assertions, 0 failed in `bash run_tests.sh`). The new `test/sql/ena_upload_reads.test` is skipped silently without `ASPERA_AVAILABLE` + `ENA_WEBIN_TEST_USER` + `ENA_WEBIN_TEST_PASSWORD`; `test/sql/ena_upload_reads_local.test` runs unguarded.

Round-trip (local file:// transport, no creds):
- `SELECT * FROM ena_upload_reads(relation := 'tbl', target_url := 'file:///tmp/uploads/')` — emits one row per output file with sample_ref, filename, filetype='fastq', md5 (32-char hex), bytes_written (UBIGINT), layout. Files are byte-exact gzipped FASTQ; on-disk md5 matches the emitted column.
- Single, paired-split (default), and paired-interleaved layouts all verified.
- Mixed-layout in one sample → InvalidInputException naming the sample_ref and offending row.
- Empty input → 0 emitted rows.
- Unknown layout name (e.g. `'triple'`) → BinderException.

Deviations from plan:
- **Materialised input + hard cap, not streaming.** The implementation reads the user's relation via the separate-connection pattern into in-memory vectors, groups by sample_ref, sorts by sequence_index NULLS LAST, then encodes/uploads sequentially. Hard caps of 50M rows AND 5GB total sequence bytes turn the silent OOM into an actionable error pointing the user at split-by-sample. Streaming would remove the cap entirely; deferred until a user trips it. (See the user follow-up: libcurl's `CURLOPT_READFUNCTION` would let us drop temp files for FTPS/HTTPS-PUT — Aspera still needs temp files because `ascp` doesn't reliably read stdin.)
- **Aspera write-side via temp directory, not a temp file.** mkdtemp gives us a unique directory; we write `{sample_ref}.fastq.gz` (or `_1`/`_2`) inside it so the basename matches the desired remote name. Cleanup unlinks the file then rmdirs the directory after `ascp --mode=send` completes (success or failure). The cleanup runs through a try/catch wrapping the entire encode → finish → ascp block, not just the ascp call — initially only ascp failures cleaned up, but encode/gzip failures also leave a temp file behind that needs unlinking.
- **`MaxThreads()=1`** — single-threaded. Per-sample concurrency is feasible later (mutex + claim queue) but not needed for the typical batch size and would obscure error attribution.
- **`secret` is required only for Aspera transport.** `file://` skips secret resolution entirely (the local test runs unguarded).

Helpers (DuckDB-free, in `src/include/ena_upload_helpers.hpp`):
- `FastqLayoutMode` enum + `FastqLayoutModeName` / `ParseFastqLayoutMode` (case-insensitive)
- `ResolveLayout(sample_ref, requested, has_r2[])` — collapses AUTO based on per-row R2 presence; throws naming the sample + the first offending row index when the requested mode is violated. Three separate first-mismatch indices (auto / first-with-r2 / first-without-r2) so error messages always point to the right row.
- `OutputFilenames(sample_ref, layout)` — emits `{ref}.fastq.gz` (single, paired-interleaved) or `{ref}_1.fastq.gz` + `{ref}_2.fastq.gz` (paired-split).
- `ParseUploadTargetURL(url)` — parses `aspera://host[/path]` and `file://path`. Forgiving on `file://` to accept relative paths (DuckDB's per-test `__TEST_DIR__` substitutes a relative path).
- `BuildAscpSendArgv(opts)` — `[ascp, --mode=send, --user=, --host=, -P, port, -i, key, -Q, -d, (-l, max_rate)?, local_path, remote_dir]`. Password is NOT in argv — it's injected via `ASPERA_SCP_PASS` in the child env.

`src/aspera_send.cpp::RunAsperaSend(argv, password)` — fork+execvp `ascp`, child env has `ASPERA_SCP_PASS` set via `setenv` (which copies the string, so the parent's std::string lifetime is fine), `prctl(PR_SET_PDEATHSIG, SIGTERM)` on Linux, stderr captured via pipe (drained EINTR-safely), waitpid loop is EINTR-safe, `_exit(127)` if `/dev/null` open fails (rather than silently inheriting parent stdout/stdin into ascp).

Code review (linus-code-reviewer): five blockers + seven should-fixes + six nits. Addressed in this commit:
- **Blocker — unbounded materialisation.** Fixed: hard caps of `MAX_INPUT_ROWS = 50M` AND `MAX_INPUT_SEQUENCE_BYTES = 5GB`, both enforced with clear `InvalidInputException` messages pointing the user at split-by-sample.
- **Blocker — temp-file leak when GzipMd5FileSink::Finish throws.** Fixed: per-file work in `RunUpload` is wrapped in a try/catch around the entire encode → finish → ascp block; cleanup_temp lambda always runs (success or failure) for the Aspera path.
- **Blocker — fd leak in child if `/dev/null` open fails.** Fixed: child `_exit(127)`s rather than inheriting the parent's stdout/stdin into ascp; opens `/dev/null` `O_RDWR` and dup2s into both fd 0 and 1.
- **Blocker — `uInt` cast truncates large inputs in zlib.** Fixed: `GzipMd5FileSink::Write` chunks input in `UINT_MAX`-sized pieces.
- **Blocker (mis-diagnosis) — SQL injection via WriteOptionallyQuoted.** Verified false: `WriteOptionallyQuoted` quotes any text containing non-identifier characters (`;`, space, `"`, etc.) and escapes embedded quotes; the unquoted return path only fires for safe identifiers. Matches the convention used by `reference_table_reader.cpp` and `placement_table_reader.cpp`. No change.
- **Should-fix — EINTR-unsafe waitpid + read in aspera_send.** Fixed: both wrapped in `while(true)` loops handling EINTR.
- **Should-fix — ResolveLayout SINGLE/PAIRED first-row reporting.** Fixed: three separate "first offending row" counters (auto-mismatch vs first-with-r2 vs first-without-r2) so each requested mode reports the correct row. Two new tests pin the regression.
- **Should-fix — global row counter in MaterialiseInput error messages.** Fixed: `global_row = out.sample_ref.size()` is used instead of the per-chunk index.
- **Should-fix — `..` path-traversal in `sample_ref`.** Fixed: rejected client-side before any path is built.
- Should-fix — `aspera_send.hpp` API disappears on non-static builds. Acknowledged but left as-is; the matching `#if` guard pattern is consistent with `aspera_stream.{cpp,hpp}` and a future use site would also need the guard.
- Should-fix — temp-dir cleanup logging when unlink fails. Acknowledged but kept the `(void)` cast; the comment now explains why (best-effort by design).
- Should-fix — `file://` relative-path footgun. Acknowledged; documented in the URL parser comment that we accept this for the `__TEST_DIR__` test fixture pattern.
- Nits not addressed: `aspera_send.hpp` API placement, `prctl` per-thread caveat, "uncompressed bytes" output column, BIGINT-vs-UTINYINT type for `qual_offset`, live-test round-trip read-back, additional `sample_ref` charset validation. None blocked review approval; defer if/when needed.

Follow-ups (deferred):
- libcurl-based streaming transport for FTPS/HTTPS-PUT (would drop temp files for those paths but not for Aspera).
- Per-sample concurrency in the encode/upload loop (would need `MaxThreads()` > 1 + claim queue + per-sample workspace).
- Round-trip md5 verification in the live Aspera test (download via `ascp --mode=recv` and compare).
- Reports API `SELECT * FROM ena.{projects,samples}` (still owed from Phase 5).

Phase 8 has been trimmed accordingly (the old "Aspera write-side" bullet moved here, now done).


### Phase 7 — pending
### Phase 8 — pending
### Phase 9 — pending
### Phase 10 — pending
