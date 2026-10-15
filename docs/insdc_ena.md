# Working with EBI / ENA

MIINT provides two complementary integrations with EMBL-EBI's [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home):

1. **Reading** — query metadata and stream sequence data directly from ENA's public APIs as if they were DuckDB tables.
2. **Submitting** — register projects, samples, experiments, runs, and analyses with ENA's [Webin V2](https://www.ebi.ac.uk/ena/submit/webin/login) submission service, including FASTQ upload via Aspera or HTTPS/FTP.

Both surfaces share the same accession types, the same authentication model (Webin account credentials for submission, anonymous for reading), and route through the same `httpfs`-backed network layer. They can be mixed in a single SQL session — for example, validating that a sample's metadata round-trips correctly between a planned submission and the published catalog.

## Table of contents

- [Reading](#reading) - Detail on reading data from ENA.
  - [Metadata](#read-ena) - Run / sample / study etc level metadata.
  - [Attributes](#read-ena-attributes) - Submitter-defined sample attributes.
  - [Sequences](#read-ena-sequences) - Stream FASTA / FASTQ / SFF data into a `read_fastx` compatible schema.
  - [Searchable fields](#ena-searchable-fields) - List the fields ENA's Portal API can return for a result type.
- [Writing](#writing) - Detail on submitting data to ENA.
  - [Register](#register-secret) - Register your submission credentials. 
  - [Attach](#attach-the-catalog) - Attach the ENA catalog.
  - [Register the project](#register-the-project) - Register a project with ENA.
  - [Register samples](#register-samples) - Register samples with ENA.
  - [Upload sequence data](#upload-sequence-data) - Upload sequence data to ENA.
  - [Register experiments and runs](#register-experiments-and-runs) - Register the experiments and runs with ENA.
  - [Audit log](#audit-log) - The recorded audit trail from the submission process.
  - [Alias targeting](#alias-targeting) - Targeting submissions by alias. 
  - [Alias lookup](#alias-lookup) - Looking up submission information by alias.
  - [Status transition rules](#status-transition-rules) - Rules for how ENA statuses change.
  - [End to end example](#end-to-end-example) - An end-to-end example of the submission process.
  - [Caveats and limitations](#caveats-and-limitations) - Known caveats and limitations to be aware of.
  
### Reading

Read paths target ENA's [Portal API](https://www.ebi.ac.uk/ena/portal/api/), [Browser XML API](https://www.ebi.ac.uk/ena/browser/api/), and the read-data FTP/HTTPS endpoints. No credentials are required.

Every metadata request miint sends to ENA's Portal, Browser XML and Reports APIs, whether from a read or from the checks an `INSERT` makes first, is retried with backoff on HTTP 429, 500, 502, 503 and 504. The Portal API sometimes answers HTTP 200 with `ERROR occurred. Not all results may have been written.` in place of the results; that response is retried the same way. If the last retry gets it too, the query fails instead of returning a partial result, with two exceptions: a `LATERAL` call to `read_ena_sequences` skips that accession with a warning in `miint_warnings()`, and `INSERT INTO ena.samples` skips its client-side checklist check, leaving validation to Webin.

A an example read pattern joins a study's metadata with its sequence data:

```sql
WITH runs AS (
    SELECT run_accession, sample_accession, library_layout
    FROM read_ena('PRJEB11419')
    WHERE library_layout = 'PAIRED'
)
-- A lateral call: options such as max_sequences can't be set here (see Read ENA Sequences).
SELECT r.sample_accession, s.read_id, length(s.sequence1) AS r1_len
FROM runs r, read_ena_sequences(r.run_accession) s;
```

### Read ENA

Query metadata from the EBI European Nucleotide Archive (ENA) Portal API. Returns run, sample, or study metadata as a table.

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to EBI servers (www.ebi.ac.uk)

**Parameters:**
- `accession` (VARCHAR or VARCHAR[]): ENA/SRA accession(s). Supports study (PRJNA\*, PRJEB\*, ERP\*, SRP\*), sample (SAMN\*, SAME\*), run (SRR\*, ERR\*), and experiment (SRX\*, ERX\*) accessions.
- `result` (VARCHAR, optional, default 'read_run'): ENA result type. One of: `read_run`, `sample`, `study`. When the accession type doesn't match the result type (e.g., run accession with `result='study'`), the accession is automatically resolved.
- `fields` (VARCHAR, optional): Comma-separated list of ENA field names to return. If omitted, a sensible default set is used per result type.

**Output schema:**
- Columns are named after the requested fields. `read_count`, `base_count` and `tax_id` are BIGINT; `fastq_ftp`, `fastq_aspera` and `fastq_md5` are VARCHAR[] and `fastq_bytes` is BIGINT[] (ENA's `;`-separated lists, split); every other field is VARCHAR
- Default `read_run` fields include: `run_accession`, `experiment_accession`, `sample_accession`, `study_accession`, `fastq_ftp`, `fastq_bytes`, `fastq_md5`, `library_strategy`, `library_layout`, `instrument_model`, `read_count`, `base_count`, and more
- Default `sample` fields include: `sample_accession`, `scientific_name`, `tax_id`, `collection_date`, `country`, `lat`, `lon`, `depth`, and more
- Default `study` fields include: `study_accession`, `study_title`, `study_description`, `center_name`, and more

**Behavior:**
- Queries the ENA Portal API (`/search` endpoint) with `format=tsv`
- Accession type is auto-detected from prefix and mapped to the appropriate query parameter
- Cross-type resolution: e.g., a run accession with `result='study'` first resolves to the study accession, then queries the study
- Rate-limited to ~3 requests/second with retry on 429/500/502/503/504 (see [Reading](#reading) for an HTTP 200 that is retried too)

**Examples:**
```sql
-- Get all run metadata for a single run
SELECT * FROM read_ena('ERR1074767');

-- Get specific fields
SELECT run_accession, library_layout, read_count
FROM read_ena('ERR1074767', fields='run_accession,library_layout,read_count');

-- Get sample metadata (auto-resolves run -> sample)
SELECT * FROM read_ena('ERR1074767', result='sample');

-- Get study info
SELECT study_title FROM read_ena('PRJEB11419', result='study');

-- Get all runs for a BioProject
CREATE TABLE project_runs AS SELECT * FROM read_ena('PRJNA555783');
```

### Read ENA Attributes

Fetch sample attributes from EBI/ENA via the Browser XML API. Returns every attribute in a sample's XML: the submitter-defined key-value attributes, named as the submitter wrote them and including those the Portal API does not carry, plus ENA's own (such as `ENA-LAST-UPDATE`).

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to EBI servers (www.ebi.ac.uk)

**Parameters:**
- `accession` (VARCHAR or VARCHAR[]): ENA/SRA accession(s). Supports study, sample, run, and experiment accessions. Non-sample accessions are automatically resolved to their associated sample accession(s).

**Output schema:**
- `sample_accession` (VARCHAR): BioSample accession (SAMN\*/SAME\*)
- `tag` (VARCHAR): Attribute name as it appears in the sample XML (e.g., 'collection date', 'geographic location (country and/or sea)', custom fields, 'ENA-LAST-UPDATE')
- `value` (VARCHAR): Attribute value

**Behavior:**
- Resolves non-sample accessions to sample accessions via the Portal API
- Fetches XML from the Browser API in batches of 200 samples
- Parses `<SAMPLE_ATTRIBUTE>` `<TAG>`/`<VALUE>` pairs
- Returns ALL attributes, including custom ones (e.g., primer sequences, custom identifiers) that `read_ena` does not carry

**Tag names and filtering:**
Tags are the names in the sample XML, so the same attribute can be spelled differently from one sample to the next, even within one study (`collection_date` in one, `collection date` in another). A `WHERE tag = ...` filter is ordinary, case-sensitive SQL equality on those names, and every sample's XML is still fetched. The Portal API's sample fields carry many of the same attributes under one name each (`collection_date`, `country`, `host_body_site`); read them with `read_ena(..., result='sample', fields='...')`, which does not fetch each sample's XML, and list them with [`ena_searchable_fields('sample')`](#ena-searchable-fields).

**Examples:**
```sql
-- Get all attributes for a run's sample
SELECT * FROM read_ena_attributes('ERR1074767');

-- Get attributes for a specific sample
SELECT * FROM read_ena_attributes('SAMEA3610311');

-- Find primer information
SELECT tag, value FROM read_ena_attributes('ERR1074767')
WHERE tag LIKE '%primer%' OR tag LIKE '%Primer%';

-- Pivot attributes to wide format for a study. Samples in one study can
-- spell the same attribute differently, so list each spelling the study uses.
SELECT sample_accession,
       MAX(CASE WHEN tag IN ('collection_date', 'collection date') THEN value END) AS collection_date,
       MAX(CASE WHEN tag IN ('country', 'geographic location (country and/or sea)') THEN value END) AS country
FROM read_ena_attributes('PRJEB11419')
GROUP BY sample_accession;

-- The Portal API's sample fields come from read_ena, not read_ena_attributes:
-- one request for the study, no per-sample XML.
SELECT sample_accession
FROM read_ena('PRJEB11419', result='sample', fields='sample_accession,host_body_site')
WHERE host_body_site = 'UBERON:feces';
```

### ENA Searchable Fields

Enumerate the fields that ENA's Portal API can return for a result type (its `/returnFields` list). Use this to discover what field names are available for a given result type (sample, read_run, study, experiment, etc.) before constructing a query — a `fields=` list for `read_ena` (result types `read_run`, `sample` and `study`), or a direct Portal API `/search` URL.

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to EBI servers (www.ebi.ac.uk)

**Parameters:**
- `result_type` (VARCHAR, required): An ENA Portal API result type (e.g., `'sample'`, `'read_run'`, `'study'`, `'experiment'`)

**Output schema:**
- `field_name` (VARCHAR): Canonical field identifier accepted by ENA's `/search?fields=...` query parameter
- `type` (VARCHAR): ENA's declared type for the field (e.g., `text`, `number`, `date`, `controlled value`, `taxonomy`)
- `description` (VARCHAR, nullable): ENA's human-readable description, when available

**Behavior:**
- Issues exactly one HTTP call to `/returnFields?result=<result_type>&format=tsv` on first use and caches the parsed result for the remainder of the scan
- Validates `result_type` as alphanumeric + `_-.` to prevent URL injection

**Examples:**
```sql
-- All fields the sample result type can return
SELECT field_name, type FROM ena_searchable_fields('sample') ORDER BY field_name;

-- Check whether a field exists before requesting it
SELECT COUNT(*) > 0 AS exists
FROM ena_searchable_fields('sample')
WHERE field_name = 'host_body_site';
```

### Read ENA Sequences

Stream FASTA/FASTQ sequence data from EBI/ENA with run, sample, and experiment accession columns. Returns data in the same schema as `read_fastx` plus accession metadata columns, enabling direct association of sequence data with project metadata. Supports mixed FASTA/FASTQ paired-end data (unlike `read_fastx` which rejects format mismatches).

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to EBI servers (ftp.sra.ebi.ac.uk via HTTPS)

**Parameters:**
- `accession` (VARCHAR or VARCHAR[]): ENA/SRA accession(s). Supports study (bulk download all runs), sample, run, and experiment accessions.
- `include_filepath` (BOOLEAN, optional, default false): Add a filepath column naming where each run was read from: its HTTPS URL(s), or `host:/path` for a run read over Aspera. For paired-end runs, the two are semicolon-separated.
- `qual_offset` (BIGINT, optional, default 33): Quality score offset (33 for Phred+33/Sanger, 64 for Phred+64/Illumina 1.3+)
- `max_sequences` (BIGINT, optional, default 0): If `> 0`, stop emitting from each run after this many sequences. `0` (or NULL / absent) means unlimited. For paired-end runs the cap counts **pairs** (one output row per pair), not underlying FASTQ records — `max_sequences => N` yields at most N rows and corresponds to 2N downloaded reads. When downloading via Aspera the cap tears down the `ascp` transfer early, saving real bandwidth. For SFF runs the cap applies but the full file is downloaded before any record is parsed; a loud warning is printed in that case.
- `trim_sff` (BOOLEAN, optional, default true): For SFF runs, apply the quality and adapter clip positions from the SFF header to trim sequences and quality scores. Ignored for FASTQ runs. Named `trim_sff` rather than `trim` because `TRIM` is a SQL function keyword and this function is dual-path (supports both scalar and lateral invocation), which together prevent DuckDB's binder from accepting `trim=...`.
- `download_method` (VARCHAR, optional, default `'auto'`): `'http'`, `'aspera'` or `'auto'`. `'auto'` uses Aspera when `ascp` is on `PATH`, an Aspera key is found or downloaded, and every FASTQ run has Aspera paths; otherwise HTTP. When `ascp` exits before sending any data (a rejected key, an unreachable host), `'auto'` warns and reads that run over HTTP, and runs of the scan that start after it go straight to HTTP; `'aspera'` raises the error with `ascp`'s message. Lateral calls always use HTTP.
- `verify_md5` (BOOLEAN, optional, default **true**): Verify each downloaded FASTQ file's bytes against ENA's reported `fastq_md5` once the file has been read to completion, raising a hard error (`IOException`) on mismatch. On by default — it's cheap (the digest is computed incrementally over bytes already in flight) and catches silent truncation/corruption that byte- or row-count checks miss. Verification is automatically skipped (with a loud warning, never silently) when it doesn't apply:
  - ENA reported no `fastq_md5` for a given file.
  - The file is not gzip-compressed (verification hashes the raw pre-decompression bytes, matching the basis of ENA's reported digest; a non-gzip file has no such basis here).
  - The run routes through SFF (`OpenSFF`) or Aspera (`download_method => 'aspera'`) — neither transport is wired to the verification tap.
  - `max_sequences` capped the run before it reached the file's true end — verification requires having seen every byte, so a deliberately partial read is never treated as a failure.

Pass named parameters with arrow syntax (`name => value`). Because `read_ena_sequences` also supports lateral / correlated invocation, `name = value` is read as a comparison and fails with a Binder Error. For example: `read_ena_sequences('X', prefer_format => 'sff', trim_sff => false)`. This works only when the accession argument is a literal (or a list of literals). In a lateral / correlated call, where the accession comes from another table's column, DuckDB doesn't bind named parameters in any syntax (passing one fails with a Binder Error), so those calls run with the defaults ([#295](https://github.com/the-miint/duckdb-miint/issues/295)).

**Output schema:**
- `sequence_index` (BIGINT): 1-based sequence index (per run)
- `read_id` (VARCHAR): FASTQ read identifier
- `comment` (VARCHAR): FASTQ comment line (nullable)
- `sequence1` (VARCHAR): Forward read sequence
- `sequence2` (VARCHAR): Reverse read sequence (NULL for single-end)
- `qual1` (UTINYINT[]): Forward quality scores
- `qual2` (UTINYINT[]): Reverse quality scores (NULL for single-end)
- `run_accession` (VARCHAR): SRR/ERR run accession
- `sample_accession` (VARCHAR): SAMN/SAME BioSample accession
- `experiment_accession` (VARCHAR): SRX/ERX experiment accession
- `filepath` (VARCHAR, optional): HTTPS URL(s), or Aspera `host:/path`(s), the run was read from

**Behavior:**
- At bind time, queries ENA Portal API to discover runs, FASTQ URLs, and library layout (paired vs single-end)
- Constructs HTTPS URLs from ENA's FTP paths
- Streams FASTQ data through the same `SequenceReader`/`DuckDBSeqStream` infrastructure as `read_fastx`
- Paired-end detection is automatic from the `library_layout` field
- Supports parallel streaming across multiple runs (up to 8 concurrent)

**Examples:**
```sql
-- Stream reads from a single run
SELECT * FROM read_ena_sequences('ERR1074767') LIMIT 100;

-- Stream all reads for a BioProject directly to Parquet
COPY (SELECT * FROM read_ena_sequences('PRJNA555783'))
  TO 'project_sequences.parquet' (FORMAT PARQUET);

-- Join sequence data with metadata
CREATE TABLE runs AS SELECT * FROM read_ena('PRJNA555783');
CREATE TABLE seqs AS SELECT * FROM read_ena_sequences('PRJNA555783');

SELECT s.run_accession, r.library_strategy, COUNT(*) AS read_count
FROM seqs s
JOIN runs r USING (run_accession)
GROUP BY ALL;

-- Stream with filepath for provenance tracking
SELECT run_accession, read_id, filepath
FROM read_ena_sequences('ERR1074767', include_filepath => true) LIMIT 5;
```

**Notes:**
- For large projects, the FASTQ download can take significant time and bandwidth
- The `run_accession` column enables easy JOIN back to metadata from `read_ena`
- Quality scores are converted to numeric values using the specified offset (default Phred+33)
- **Failure handling**: on a transient open failure, the run is retried once.
  If the retry also fails, a single-run literal scan raises the error; a scan
  of several runs skips the run with a loud warning and continues. If every run
  of such a scan fails, it raises `read_ena_sequences: every run failed (<n>):
  <accessions>; first error: <message>` at the end rather than return no rows;
  so does a single-run scan whose run fails mid-stream, rather than return its
  partial rows. A lateral call skips a failed run with a warning, even when
  every run fails. On
  a mid-stream failure (connection dropped after some rows have been emitted),
  the run is NOT retried — re-reading from scratch would emit the same reads
  again with duplicate `sequence_index` values downstream. Instead, the run is
  recorded as skipped, a loud warning reports the run accession and the
  number of reads emitted before the failure, and the scan continues with
  other runs. If you see such a warning, re-run the query (the metadata
  lookup is cached) or use a smaller per-run selection to recover the
  truncated data.
- **md5 verification** (see `verify_md5` above, on by default): verification
  runs only after a file is read to its true end, so a mismatch between ENA's
  reported `fastq_md5` and the downloaded bytes means exactly one thing: the
  complete byte stream read for that file does not hash to the digest ENA
  published. It does not establish *why*, and the two causes are
  indistinguishable in the raised error: (a) a transfer cut on a gzip member
  boundary decompresses cleanly and looks like a normal end of stream, so it
  reaches verification with no transport error, and a re-run may well
  succeed; (b) the bytes genuinely disagree with the published digest, and a
  re-run reproduces the failure. A caller that retries should bound its
  retries rather than treat the error as always-permanent or
  always-transient. A truncation *inside* a gzip member is a different
  failure — a mid-stream read error, skipped with a warning, never an md5
  mismatch.

  Who receives the raised exception: a scan of a **single run** raises the
  mismatch directly to the caller. A scan of **several runs** (a `varchar[]`
  of accessions, or a project accession that expands to many runs) and **any
  lateral** call instead skip the run with a loud warning and continue, the
  same treatment as any other completion-time integrity failure, since one
  corrupt run should not discard rows already emitted for its siblings. If
  every run of a scan of several runs fails, the scan raises the "every run
  failed" error above, which quotes the first failure recorded.

  **Message contract**: the raised message is `read_ena_sequences: md5
  mismatch for '<label>': ENA reported <expected> but downloaded bytes hash
  to <actual>`, and never contains transport-failure wording — `connection`,
  `timed out`, `timeout`, `network`, `reset`, `refused`, `unreachable`,
  `temporarily`, `curl` — so a caller can classify on the text alone. This is
  pinned by `test/cpp/test_ena_md5.cpp`; rewording the message is a breaking
  change. (This describes the exception's message content —
  `duckdb::Exception` stores a JSON envelope, so a caller logging `e.what()`
  sees the message inside that JSON, not the bare string.)

**Lateral invocation (correlated arguments):**

Accessions can come from a correlated column via `LATERAL`. Each outer row opens its
own run; a `LIMIT` inside the lateral short-circuits that run's download (the reader
and any HTTP connection are torn down as soon as the limit is reached). Metadata
lookups share an in-memory LRU cache scoped to the query, so repeated outer-row
values and within-query batches avoid redundant ENA Portal API calls.

```sql
-- Find runs containing a probe sequence, downloading only until each match is seen.
SELECT r.run_accession
FROM read_ena('PRJEB11419') AS r,
     LATERAL (SELECT 1 FROM read_ena_sequences(r.run_accession)
              WHERE sequence1.contains('AACGTAGGTCACAAGCGTTGTCCGGA')
              LIMIT 1);
```

Limitations in lateral mode:
- Named parameters (`max_sequences`, `verify_md5`, `include_filepath`, ...) can't be passed in any syntax: the call
  fails to bind. Every option keeps its default ([#295](https://github.com/the-miint/duckdb-miint/issues/295)).
- Downloads always use HTTP, even when `ascp` is available (the lateral use case is
  short-circuit-driven, not throughput-driven).
- `table_scan_progress` reports `-1.0` (indeterminate) because the total work is
  driven by the outer side and not known at bind time.

## Writing

Submissions go through [ENA Webin V2](https://docs.ena.ebi.ac.uk/), which expects an XML envelope per object type and authenticates with a Webin-* account. duckdb-miint exposes this as a virtual catalog: register a Webin secret, `ATTACH` an `ena` database, and `INSERT INTO` the fixed submission tables.

The full submission topology of one cohort looks like:

```
project (1) ── studies the cohort
   │
   ├── sample (N) ── biological materials, with metadata attributes
   │       │
   │       └── experiment (M) ── library prep + platform
   │               │
   │               └── run (R) ── the actual sequence files (uploaded separately)
   │
   └── analysis (optional) ── derived results (assemblies, variant calls, …)
```

Every step uses the same SQL idiom: `INSERT INTO ena.<table> … RETURNING <server_assigned_columns>`. The `RETURNING` clause is how you capture ENA's accession (PRJEB / SAMEA / ERX / ERR / ERZ) and feed it into the next step.

### Register secret

```sql
-- Test endpoint (wwwdev.ebi.ac.uk) — use this for dry runs.
CREATE SECRET my_ena (
    TYPE ENA,
    USER 'Webin-12345',
    PASSWORD_ENV 'WEBIN_PASSWORD',     -- or PASSWORD '...' or PASSWORD_FILE '/path'
    ENDPOINT 'test'                    -- 'test' (default) or 'production'
);
```

| Parameter | Required | Notes |
|---|---|---|
| `USER` | yes | Your Webin submission account, e.g. `Webin-12345`. |
| `PASSWORD` | one of three | Literal string. Stored redacted; `duckdb_secrets()` shows `password=redacted`. |
| `PASSWORD_ENV` | one of three | Env var name (resolved at create time). Preferred for CI. |
| `PASSWORD_FILE` | one of three | Path to a file whose first line is the password. CRLF is stripped. |
| `ENDPOINT` | no | `'test'` (default; `wwwdev.ebi.ac.uk`) or `'production'` (`www.ebi.ac.uk`). |

Exactly one of `PASSWORD` / `PASSWORD_ENV` / `PASSWORD_FILE` must be provided. Specifying zero or more than one is an error.

### Attach the catalog

```sql
ATTACH '' AS ena (TYPE ENA, SECRET 'my_ena');
```

This registers a fixed-schema virtual catalog with one schema (`main`) and six tables. `CREATE TABLE` and other DDL against `ena.*` are rejected — the schema is hard-coded to match the Webin V2 object model.

| Table | Purpose | Columns |
|---|---|---|
| `ena.projects` | BioProject (`PRJEB`) registration. | `alias`, `title`, `description`, `project_type`, `umbrella_children`, `attributes`, `prjeb_accession`, `erp_accession`, `hold_until_date` |
| `ena.samples` | BioSample (`SAMEA`) registration with checklist-driven attributes. | `alias`, `title`, `description`, `taxon_id`, `scientific_name`, `checklist`, `attributes`, `attribute_units`, `ers_accession`, `samea_accession` |
| `ena.experiments` | Library / platform metadata (`ERX`). | `alias`, `title`, `study_ref`, `sample_descriptor`, `design_description`, `library_name`, `library_strategy`, `library_source`, `library_selection`, `library_layout`, `platform`, `instrument_model`, `erx_accession` |
| `ena.runs` | Sequence-file registration (`ERR`). | `alias`, `experiment_ref`, `title`, `files` (`LIST<STRUCT(filename, filetype, md5)>`), `err_accession` |
| `ena.analyses` | Derived results (`ERZ`). INSERT not yet implemented. | `alias`, `study_ref`, `analysis_type`, `accession` |
| `ena.submission_log` | Append-only in-memory bookkeeping. | `submission_id`, `submitted_at`, `endpoint`, `secret_name`, `action`, `object_type`, `n_objects`, `success`, `era_accession`, `request_payload`, `receipt`, `error_messages`, `duration_ms`, `target`, `object_aliases`, `object_accessions` |

The `submission_log` is the only table you `SELECT` from directly. The other five are submit-only — `SELECT` against them is not implemented (the Reports API powers the alias→accession lookup behind the lifecycle table fns, but doesn't surface a `SELECT *` view of the user's submission account; see Caveats). `INSERT INTO ena.analyses` is reserved for a future build and currently raises a binder exception. `DELETE FROM ena.<projects/samples/experiments/runs>` issues a `CANCEL` against Webin V2 (see Lifecycle).

### Register the project

```sql
INSERT INTO ena.projects (alias, title, project_type)
VALUES ('my-cohort-2026', 'Gut microbiome cohort 2026', 'METAGENOMIC')
RETURNING prjeb_accession, erp_accession;
-- ┌──────────────────┬────────────────┐
-- │ prjeb_accession  │ erp_accession  │
-- │     varchar      │    varchar     │
-- ├──────────────────┼────────────────┤
-- │ PRJEB12345       │ ERP00009876    │
-- └──────────────────┴────────────────┘
```

Aliases are submitter-chosen, must be unique within your account, and act as the foreign-key handle for downstream `_ref` columns.

### Register samples

Samples are validated against an ENA [checklist](https://www.ebi.ac.uk/ena/browser/checklists) (e.g. `ERC000015` for human gastrointestinal). Mandatory attributes for the checklist are checked client-side **before** the request is sent — missing fields surface as a clear error rather than an opaque server rejection.

```sql
INSERT INTO ena.samples (alias, taxon_id, checklist, attributes, attribute_units)
VALUES
  ('s1', 408170, 'ERC000015',
   MAP {
     'collection date': '2026-04-01',
     'project name': 'cohort-2026',
     'geographic location (country and/or sea)': 'USA',
     'geographic location (latitude)': '32.7157',
     'geographic location (longitude)': '-117.1611',
     'broad-scale environmental context': 'human-associated habitat [ENVO:00009003]',
     'local environmental context': 'gastrointestinal tract [ENVO:00002041]',
     'environmental medium': 'feces [UBERON:0001988]'
   },
   MAP {
     'geographic location (latitude)': 'DD',
     'geographic location (longitude)': 'DD'
   })
RETURNING samea_accession;
```

`attribute_units` is a sibling map keyed by the same attribute name. ENA's checklist defines which attributes require a unit (latitude/longitude need `'DD'`, depth needs `'m'`, etc.). Attributes that don't take a unit are simply absent from the map.

### Upload sequence data

Run registration is two steps: (a) upload the FASTQ files to your Webin "drop box" and capture the per-file MD5s, then (b) `INSERT INTO ena.runs` referencing those filenames + MD5s.

`ena_upload_reads` drives the upload from an in-DuckDB relation, encoding sequences and quality scores into FASTQ format on the fly, gzipping, computing MD5 in a single streaming pass, and shipping the bytes to ENA.

```sql
-- Required input columns: sample_ref, read_id, sequence1, qual1, sequence_index.
-- Optional: sequence2, qual2 (paired-end).
CREATE TABLE my_reads AS SELECT … ;

CREATE TABLE upload_results AS
SELECT * FROM ena_upload_reads(
    relation   := 'my_reads',
    secret     := 'my_ena',
    target_url := 'aspera://webin2.ebi.ac.uk/'   -- aspera://, file://, ftp(s)://, http(s)://
);
-- ┌────────────┬─────────────────┬──────────┬────────────────────┬───────────────┬────────┐
-- │ sample_ref │    filename     │ filetype │        md5         │ bytes_written │ layout │
-- │  varchar   │     varchar     │ varchar  │      varchar       │    ubigint    │ varchar│
-- ├────────────┼─────────────────┼──────────┼────────────────────┼───────────────┼────────┤
-- │ s1         │ s1_1.fastq.gz   │ fastq    │ 9b8932f85caa54e6…  │      87432    │ paired │
-- │ s1         │ s1_2.fastq.gz   │ fastq    │ 183d6a24e0c3704e…  │      88105    │ paired │
-- └────────────┴─────────────────┴──────────┴────────────────────┴───────────────┴────────┘
```

Layout is auto-detected from the input rows: presence of `sequence2` / `qual2` per sample-group decides single vs. paired. `layout := 'paired_interleaved'` overrides this to emit one interleaved file per sample. Mixed single/paired rows under the same `sample_ref` are rejected at scan time.

`target_url` schemes:
- `aspera://webin2.ebi.ac.uk/` — production Aspera transport, requires `ascp` on `PATH` (`MIINT_ENABLE_ASPERA` runtime check).
- `ftp://` / `ftps://` / `http://` / `https://` — libcurl streaming-upload transport.
- `file:///some/dir/` — write to a local directory; no `secret` parameter required. Useful for testing the encode→gzip→md5 pipeline without round-tripping through ENA.

### Register experiments and runs

```sql
INSERT INTO ena.experiments (
    alias, study_ref, sample_descriptor,
    library_strategy, library_source, library_selection, library_layout,
    platform, instrument_model)
VALUES ('e-s1', 'my-cohort-2026', 's1',
        'WGS', 'METAGENOMIC', 'RANDOM', 'PAIRED',
        'ILLUMINA', 'Illumina NovaSeq 6000')
RETURNING erx_accession;

INSERT INTO ena.runs (alias, experiment_ref, files)
SELECT 'r-' || sample_ref,
       'e-' || sample_ref,
       LIST({'filename': filename, 'filetype': filetype, 'md5': md5} ORDER BY filename)
FROM upload_results
GROUP BY sample_ref
RETURNING alias, err_accession;
```

The `files` column expects a `LIST<STRUCT(filename, filetype, md5)>` — exactly the shape `ena_upload_reads` emits, modulo the per-sample `LIST(...)` aggregation. Empty lists, NULL entries, empty `filetype`, and empty `md5` are rejected at bind time (the server would reject them anyway, this just surfaces the error earlier).

### Audit log

Every `INSERT` (and every lifecycle action — see below) appends a row to `ena.submission_log`, including the request payload, the raw receipt XML, and a stable `submission_id`. This is the source of truth for what was sent and when:

```sql
SELECT submitted_at, object_type, n_objects, success, era_accession, duration_ms
FROM ena.submission_log
ORDER BY submitted_at DESC
LIMIT 20;
```

Failed submissions are also logged (`success = false`, `error_messages` populated), and the `INSERT` itself raises so the surrounding transaction sees the failure.

`ena_upload_reads` is **not** logged — it's a transport-only operation (encode + gzip + ship FASTQ to webin2.ebi.ac.uk) and doesn't hit the Webin V2 submission endpoint. The server-side audit for uploads lives in your Webin account dashboard. The subsequent `INSERT INTO ena.runs` that references the uploaded files IS logged.

## Lifecycle: cancel, release, embargo, modify

Webin V2 supports five lifecycle actions on already-registered objects. The miint extension exposes them as table functions, plus a SQL `DELETE` shortcut for CANCEL:

| Action | Surface |
|---|---|
| Cancel | `ena_cancel(secret, accession \| (refname, kind), catalog?)` &nbsp;or&nbsp; `DELETE FROM ena.<table> WHERE <accession_col> \| alias = '…'` |
| Release | `ena_release(secret, accession \| (refname, kind), catalog?)` |
| Hold (embargo until a date) | `ena_hold(secret, accession \| (refname, kind), until, catalog?)` |
| Modify (project) | `ena_modify_project(secret, accession, alias, title, …)` |
| Modify (sample) | `ena_modify_sample(secret, accession, alias, taxon_id, attributes, …)` |
| Modify (experiment) | `ena_modify_experiment(secret, accession, alias, study_ref, sample_descriptor, …)` |
| Modify (run) | `ena_modify_run(secret, accession, alias, experiment_ref, files, …)` |
| Validate (dry-run) | `SET miint_ena_validate_only = true;` (any `INSERT INTO ena.X (…)`) |

All four lifecycle table functions return a single row with `action`, `target`, `success`, `era_accession`, `hold_until_date`, `error_messages`, `duration_ms` and append a row to `ena.submission_log`. The MODIFY family returns the same row shape as the corresponding `INSERT`'s `RETURNING`.

```sql
-- Cancel a previously-registered project by accession.
SELECT * FROM ena_cancel(secret => 'my_ena', accession => 'PRJEB12345');
-- Equivalent SQL DELETE:
DELETE FROM ena.projects WHERE prjeb_accession = 'PRJEB12345';

-- Release a held sample by accession.
SELECT * FROM ena_release(secret => 'my_ena', accession => 'ERS999999');

-- Embargo a sample until a future date.
SELECT * FROM ena_hold(secret => 'my_ena', accession => 'ERS999999', until => '2027-12-31');
```

**lifecycle rows:**

Lifecycle calls log the action plus the `target` accession (post-translation if the call was alias-targeted). `object_aliases` / `object_accessions` are empty for lifecycle rows — those parallel arrays only carry per-object children of an `INSERT` receipt, and lifecycle ops don't register new objects. `era_accession` is also empty: Webin V2 doesn't assign a submission accession for cross-submission lifecycle ops. The `success` column reflects whether the action took effect (the receipt's `<RECEIPT @success>` attribute is unreliable for lifecycle and is parsed semantically here from `<INFO>` vs `<ERROR>` children).

`object_type` is populated from the `kind` named parameter on alias-targeted lifecycle calls, and from the table name on `DELETE FROM ena.X` calls. Accession-targeted lifecycle calls (where the user passes only `accession`, no `kind` to disambiguate) leave `object_type` empty — the function is kind-agnostic on that path.

```sql
SELECT submitted_at, object_type, action, target, success, duration_ms
FROM ena.submission_log
WHERE action IN ('CANCEL', 'RELEASE', 'HOLD')
ORDER BY submitted_at DESC LIMIT 20;
```

### Alias targeting

The lifecycle table functions also accept the alias the object was originally registered with. The extension translates the alias to its server-assigned accession via the Webin Reports API automatically before issuing the lifecycle action — one HTTP GET to resolve, then the lifecycle POST. Pass exactly one of `accession` or `refname` (both at once is rejected at bind time).

When `refname` is set, the sibling `kind` named parameter is **required**. Aliases are unique per-account-per-kind on Webin; the kind disambiguates a reused alias. `kind` is one of `'projects'` / `'samples'` / `'experiments'` / `'runs'` (matches the `ena.<table>` names). The accession-targeted path is kind-agnostic (the prefix encodes it).

```sql
-- Cancel by alias — Reports API resolves to PRJEB12345 transparently.
SELECT * FROM ena_cancel(secret => 'my_ena',
                         refname => 'my-cohort-2026',
                         kind => 'projects');

-- DELETE-by-alias works the same way; kind comes from the table.
DELETE FROM ena.projects WHERE alias = 'my-cohort-2026';

-- Hold a sample by alias.
SELECT * FROM ena_hold(secret => 'my_ena',
                       refname => 's-pilot-001',
                       kind => 'samples',
                       until => '2027-12-31');
```

If the alias isn't registered in the submission account, the call surfaces a friendly `refname '…' not found in this submission account (kind=…)` error from the Reports lookup before any lifecycle POST is attempted. Verify the kind first — a wrong-kind lookup also returns "not found".

### Alias lookup

`ena.submission_log` retains every alias and its server-assigned accession from the `INSERT` receipt in two parallel `LIST(VARCHAR)` columns: `object_aliases` and `object_accessions`. Within the same session you can recover an accession without hitting the Reports API:

```sql
-- Inspect first, then act — useful when you want to confirm the accession
-- (e.g. log it, sanity-check the kind) before the lifecycle call.
SELECT object_accessions[list_position(object_aliases, 'my-cohort-2026')] AS prjeb
FROM ena.submission_log
WHERE list_contains(object_aliases, 'my-cohort-2026')
  AND object_type = 'projects'
  AND success
ORDER BY submitted_at DESC LIMIT 1;
-- → 'PRJEB12345'

DELETE FROM ena.projects WHERE prjeb_accession = 'PRJEB12345';
```

`submission_log` is **in-memory per attached catalog** — `DETACH ena` empties it. For the common alias-targeted lifecycle case, the table-function `refname` + `kind` form (above) handles cross-session lookup automatically. For long-term audit / forensic use, persist the log:

```sql
-- Persist the alias↔accession map across sessions (audit only — the
-- refname+kind form already handles cross-session lookup transparently):
COPY (
    SELECT submitted_at, object_type, object_aliases, object_accessions
    FROM ena.submission_log WHERE success
) TO 'ena_objects.parquet' (FORMAT PARQUET);
```

### Status-transition rules

ENA's submission status machine forbids certain transitions. The most common surprise is **PUBLIC → CANCELLED** — once an object has been `RELEASE`d (made public), it cannot subsequently be cancelled. The lifecycle table functions surface this as a server-side error message when you try, with no client-side gate (every other rule that the extension knows about gates at bind time). Plan accordingly: if you want to clean an object up, don't `RELEASE` it first.

Other transitions to be aware of:

- `HOLD → CANCEL`: works (held objects can be withdrawn).
- `CANCEL` is terminal — cancelled objects can't be revived.
- `RELEASE` is one-way for the submitter (the object goes PUBLIC and can be browsed by anyone; only ENA staff can suppress).

## End-to-end example

```sql
LOAD miint;

CREATE SECRET my_ena (TYPE ENA, USER 'Webin-12345',
                     PASSWORD_ENV 'WEBIN_PASSWORD', ENDPOINT 'test');
ATTACH '' AS ena (TYPE ENA, SECRET 'my_ena');

-- 1. Project
INSERT INTO ena.projects (alias, title, project_type)
VALUES ('cohort-2026', 'Gut microbiome cohort 2026', 'METAGENOMIC');

-- 2. Samples (truncated; full attribute set per checklist)
INSERT INTO ena.samples (alias, taxon_id, checklist, attributes, attribute_units)
SELECT 's' || row_id, 408170, 'ERC000015', attrs, units
FROM my_local_metadata;

-- 3. Upload FASTQs
CREATE TABLE uploads AS
SELECT * FROM ena_upload_reads(relation := 'my_reads', secret := 'my_ena',
                               target_url := 'aspera://webin2.ebi.ac.uk/');

-- 4. Experiments
INSERT INTO ena.experiments (alias, study_ref, sample_descriptor, library_strategy,
                             library_source, library_selection, library_layout,
                             platform, instrument_model)
SELECT 'e-' || alias, 'cohort-2026', alias,
       'WGS', 'METAGENOMIC', 'RANDOM', 'PAIRED',
       'ILLUMINA', 'Illumina NovaSeq 6000'
FROM ena_sample_aliases;

-- 5. Runs
INSERT INTO ena.runs (alias, experiment_ref, files)
SELECT 'r-' || sample_ref, 'e-' || sample_ref,
       LIST({'filename': filename, 'filetype': filetype, 'md5': md5} ORDER BY filename)
FROM uploads
GROUP BY sample_ref;

DETACH ena;
```

## Caveats and limitations

- `SELECT * FROM ena.<submission table>` is **not** supported (other than `submission_log`). The Webin V2 API doesn't expose registered objects through the same endpoint they were submitted to. The extension uses the [Reports API](https://www.ebi.ac.uk/ena/submit/report) for the alias→accession lookup that powers `refname` + `kind` lifecycle calls and `DELETE WHERE alias = '…'`, but doesn't surface a `SELECT *` view of the user's submission account through the catalog.
- The `test` endpoint (`wwwdev.ebi.ac.uk`) is a sandbox: receipts come back with valid-looking accessions, but those accessions are not addressable through the public Browser. Use it for dry runs and CI; switch to `production` for real submissions.
- File uploads require either `ascp` (Aspera) on `PATH` or libcurl-backed streaming. Aspera is faster for large files but requires the IBM Aspera client; HTTPS/FTP works without extra dependencies. The local `file://` transport is for testing the encode pipeline only — it does not push anything to ENA.
- WASM builds: ENA *reading* works in all WASM variants. ENA *submission* requires an upload transport — Aspera is unavailable on any WASM (no `fork`/`exec`) and libcurl is gated off on `wasm_threads` (vcpkg port lacks `-pthread`); only `wasm_eh` and `wasm_mvp` can submit, and only via curl-based HTTPS/FTP transports.
