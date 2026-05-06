# Working with EBI / ENA

duckdb-miint provides two complementary integrations with EMBL-EBI's [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home):

1. **Reading** — query metadata and stream sequence data directly from ENA's public APIs as if they were DuckDB tables.
2. **Submitting** — register projects, samples, experiments, runs, and analyses with ENA's [Webin V2](https://www.ebi.ac.uk/ena/submit/webin/login) submission service, including FASTQ upload via Aspera or HTTPS/FTP.

Both surfaces share the same accession types, the same authentication model (Webin account credentials for submission, anonymous for reading), and route through the same `httpfs`-backed network layer. They can be mixed in a single SQL session — for example, validating that a sample's metadata round-trips correctly between a planned submission and the published catalog.

## Reading

Read paths target ENA's [Portal API](https://www.ebi.ac.uk/ena/portal/api/), [Browser XML API](https://www.ebi.ac.uk/ena/browser/api/), and the read-data FTP/HTTPS endpoints. No credentials are required.

| Function | Purpose |
|---|---|
| [`read_ena`](table-functions.md#read_enaaccession-resultread_run-fields) | Run / sample / study metadata via the Portal `/search` API. |
| [`read_ena_attributes`](table-functions.md#read_ena_attributesaccession) | Submitter-defined sample attributes via the Browser XML API, with structured-search pushdown when filters target known fields. |
| [`ena_searchable_fields`](table-functions.md#ena_searchable_fieldsresult_type) | Enumerate the field allowlist that ENA's Portal `/search?result=…` endpoint accepts as structured filters. |
| [`read_ena_sequences`](table-functions.md#read_ena_sequencesaccession-include_filepathfalse-qual_offset33-max_sequences0) | Stream FASTA / FASTQ sequence data with run / sample / experiment columns attached, in `read_fastx`-compatible schema. |

A typical read pattern joins a study's metadata with its sequence data:

```sql
WITH runs AS (
    SELECT run_accession, sample_accession, library_layout
    FROM read_ena('PRJEB11419')
    WHERE library_layout = 'PAIRED'
)
SELECT r.sample_accession, s.read_id, length(s.sequence1) AS r1_len
FROM runs r, read_ena_sequences(r.run_accession, max_sequences => 1000) s;
```

## Submitting

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

### 1. Register a Webin secret

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

### 2. Attach the ENA catalog

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
| `ena.submission_log` | Append-only in-memory bookkeeping. | `submission_id`, `submitted_at`, `endpoint`, `secret_name`, `action`, `object_type`, `n_objects`, `success`, `era_accession`, `request_payload`, `receipt`, `error_messages`, `duration_ms` |

The `submission_log` is the only table you `SELECT` from. The other five are write-only — `SELECT` against them is reserved for a future Reports-API integration. `INSERT INTO ena.analyses` is reserved for a future build and currently raises a binder exception.

### 3. Register the project

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

### 4. Register samples

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

### 5. Upload sequence files

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

### 6. Register experiments and runs

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

Every INSERT appends a row to `ena.submission_log`, including the request payload, the raw receipt XML, and a stable `submission_id`. This is the source of truth for what was sent and when:

```sql
SELECT submitted_at, object_type, n_objects, success, era_accession, duration_ms
FROM ena.submission_log
ORDER BY submitted_at DESC
LIMIT 20;
```

Failed submissions are also logged (`success = false`, `error_messages` populated), and the `INSERT` itself raises so the surrounding transaction sees the failure.

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

## Caveats and limits

- `SELECT * FROM ena.<submission table>` is **not** supported (other than `submission_log`). The Webin V2 API doesn't expose registered objects through the same endpoint they were submitted to; round-tripping requires the [Reports API](https://www.ebi.ac.uk/ena/submit/report) (future work).
- The `test` endpoint (`wwwdev.ebi.ac.uk`) is a sandbox: receipts come back with valid-looking accessions, but those accessions are not addressable through the public Browser. Use it for dry runs and CI; switch to `production` for real submissions.
- File uploads require either `ascp` (Aspera) on `PATH` or libcurl-backed streaming. Aspera is faster for large files but requires the IBM Aspera client; HTTPS/FTP works without extra dependencies. The local `file://` transport is for testing the encode pipeline only — it does not push anything to ENA.
- WASM builds: ENA *reading* works in all WASM variants. ENA *submission* requires an upload transport — Aspera is unavailable on any WASM (no `fork`/`exec`) and libcurl is gated off on `wasm_threads` (vcpkg port lacks `-pthread`); only `wasm_eh` and `wasm_mvp` can submit, and only via curl-based HTTPS/FTP transports.
