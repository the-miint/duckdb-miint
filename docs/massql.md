# Mass Spectrometry & MassQL

MIINT brings [MassQL](https://github.com/mwang87/MassQueryLanguage) to DuckDB, letting you search mass spectrometry data using the same query syntax you know from the Python MassQL library — but running on top of a high-performance analytical database.

## Getting Started

If you're new to DuckDB, it's a fast SQL database that runs locally (no server needed). Install it from [duckdb.org](https://duckdb.org/docs/installation/), then:

```sql
-- One-time setup: install and load the MIINT extension
INSTALL miint FROM community;
LOAD miint;

-- Run a MassQL query directly against a file
SELECT * FROM massql(
  'QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857',
  'sample.mzML'
);
```

That's it. The `massql()` function takes your MassQL query string and a data source (file path or table name), and returns results as a table.

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `query` | VARCHAR | MassQL query string (positional, required) |
| `source` | VARCHAR | File path or table/view name (positional, required) |
| `sample_id` | VARCHAR | Column name to iterate over per-sample (named, optional) |

When `sample_id` is provided, the function queries distinct values of that column, runs the MassQL pipeline independently for each value, and streams results back with the sample identifier prepended as the first output column. The output column preserves the original column's type (VARCHAR, INTEGER, etc.).

## Quick Comparison: Python MassQL vs MIINT

| Python MassQL | MIINT (DuckDB) |
|---------------|----------------|
| `results_df = msql_engine.process_query("QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857", "sample.mzML")` | `SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857', 'sample.mzML');` |
| Returns a pandas DataFrame | Returns a table you can further query, filter, join, or export |
| Runs single-threaded in Python | Runs in DuckDB's parallel execution engine |
| Must fit in memory | Streams results; handles large files efficiently |

The MassQL query syntax inside the quotes is identical. The only difference is how you call it.

## MassQL Syntax Reference

### Basic Structure

```
QUERY [aggregation()] datatype [WHERE conditions] [FILTER conditions]
```

- **datatype**: `MS1DATA` or `MS2DATA`
- **WHERE**: Conditions that spectra must match (AND logic between conditions)
- **FILTER**: Post-filter on individual peaks within already-matched spectra

### Condition Fields

| Field | Description | Used with |
|-------|-------------|-----------|
| `MS1MZ` | MS1 peak m/z | MS1DATA |
| `MS2PROD` | MS2 product ion m/z | MS2DATA |
| `MS2PREC` | MS2 precursor m/z | MS2DATA |
| `MS2NL` | Neutral loss (precursor minus product) | MS2DATA |
| `RTMIN` | Minimum retention time (minutes) | Both |
| `RTMAX` | Maximum retention time (minutes) | Both |
| `SCANMIN` / `SCANMAX` | Scan number range | Both |
| `CHARGE` | Precursor charge state | MS2DATA |
| `POLARITY` | `Positive` or `Negative` | Both |

### Qualifiers

Qualifiers modify a condition. Append them with `:QUALIFIER=value`:

```sql
-- m/z tolerance (default is 0.1 Da)
SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857:TOLERANCEMZ=0.5', 'sample.mzML');
SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857:TOLERANCEPPM=5', 'sample.mzML');

-- Intensity thresholds
-- INTENSITYPERCENT: peak must be >= N% of base peak intensity
SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857:INTENSITYPERCENT=50', 'sample.mzML');

-- INTENSITYVALUE: peak must have >= N absolute intensity
SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857:INTENSITYVALUE=1000', 'sample.mzML');

-- INTENSITYTICPERCENT: peak must be >= N% of total ion current
SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857:INTENSITYTICPERCENT=10', 'sample.mzML');
```

**Additional qualifiers:**

| Qualifier | Description | Example |
|-----------|-------------|---------|
| `TOLERANCEMZ` | Absolute m/z tolerance (Da) | `:TOLERANCEMZ=0.5` |
| `TOLERANCEPPM` | Relative m/z tolerance (ppm) | `:TOLERANCEPPM=5` |
| `INTENSITYPERCENT` | Minimum % of base peak | `:INTENSITYPERCENT=50` |
| `INTENSITYVALUE` | Minimum absolute intensity | `:INTENSITYVALUE=1000` |
| `INTENSITYTICPERCENT` | Minimum % of TIC | `:INTENSITYTICPERCENT=10` |
| `MASSDEFECT` | Mass defect range | `:MASSDEFECT=range(min=0.1,max=0.3)` |
| `EXCLUDED` | Exclude matching spectra | `:EXCLUDED` |
| `CARDINALITY` | Min/max matches for OR lists | `:CARDINALITY=range(min=2,max=3)` |
| `OTHERSCAN` | Expand to neighboring scans by RT | `:OTHERSCAN=rtrange(left=0.5,right=0.5)` |

### Aggregation Functions

Aggregation controls what the query returns. Without aggregation, you get individual peaks.

```sql
-- No aggregation: returns individual peaks (m/z + intensity for each match)
SELECT * FROM massql('QUERY MS2DATA WHERE MS2PROD=167.0857', 'sample.mzML');

-- scaninfo: one row per matching spectrum, with metadata
SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857', 'sample.mzML');

-- scannum: just the spectrum indices (useful for counting)
SELECT * FROM massql('QUERY scannum(MS2DATA) WHERE MS2PROD=167.0857', 'sample.mzML');

-- scansum: sum of peak intensities per spectrum
SELECT * FROM massql('QUERY scansum(MS2DATA) WHERE MS2PROD=167.0857', 'sample.mzML');
```

| Aggregation | What you get |
|-------------|-------------|
| *(none)* | Every matching peak: m/z, intensity, normalized intensity |
| `scaninfo()` | One row per spectrum: scan number, RT, precursor info, TIC |
| `scannum()` | Just spectrum indices (for counting matches) |
| `scansum()` | Total intensity per spectrum |
| `scanmz()` | Distinct precursor m/z values |
| `scanmaxint()` | Max intensity per spectrum |
| `scanrangesum()` | Binned intensity sums per spectrum |

### X-Variable Patterns

X-variables let you find spectra where peaks occur at consistent m/z offsets, without specifying the absolute m/z:

```sql
-- Find spectra with two peaks separated by 14 Da (CH2 mass shift)
SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=X AND MS2PROD=X+14', 'sample.mzML');

-- Using chemical formulas for the offset
SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=X AND MS2PROD=X+formula(Fe)', 'sample.mzML');

-- Iron isotope pattern: coefficient on X
SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=X AND MS2PROD=2*(X-formula(Fe))', 'sample.mzML');

-- Cross-level: find MS2 spectra whose precursor matches an MS1 peak
SELECT * FROM massql('QUERY scannum(MS2DATA) WHERE MS1MZ=X AND MS2PREC=X', 'sample.mzML');

-- Constrain the X search range
SELECT * FROM massql(
  'QUERY scaninfo(MS2DATA) WHERE MS2PROD=X AND MS2PROD=X+14 AND X=range(min=100,max=500)',
  'sample.mzML'
);
```

### Y-Variable Patterns (Intensity Matching)

Y-variables find spectra where peaks have a specific intensity ratio:

```sql
-- Peak at X is the reference; peak at X+28 must be ~50% of reference intensity (within 20%)
SELECT * FROM massql(
  'QUERY scaninfo(MS2DATA) WHERE MS2PROD=X:INTENSITYMATCHREFERENCE AND MS2PROD=X+28:INTENSITYMATCH=Y*0.5:INTENSITYMATCHPERCENT=20',
  'sample.mzML'
);
```

### Chemical Formula Functions

Use `formula()` to compute monoisotopic mass from element symbols:

```sql
-- In MassQL queries
SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=X AND MS2PROD=X+formula(Fe)', 'sample.mzML');

-- As a standalone function (returns mass in Daltons)
SELECT formula('Fe');        -- 55.934936
SELECT formula('C6H12O6');   -- 180.063388
SELECT formula('H2O');       -- 18.010565
```

Additional mass functions:
- `aminoaciddelta(sequence)` — residue mass for amino acid sequences (e.g., `aminoaciddelta(GK)`)
- `peptide(sequence, charge=N, ion=X)` — fragment ion mass (e.g., `peptide(PEPTIDEK,charge=2,ion=b)`)

## Practical Examples

### Find spectra containing a product ion

```sql
SELECT * FROM massql(
  'QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857:TOLERANCEPPM=5',
  'sample.mzML'
);
```

### Precursor + product search

```sql
SELECT * FROM massql(
  'QUERY scannum(MS2DATA) WHERE MS2PREC=349.18 AND MS2PROD=167.09',
  'sample.mzML'
);
```

### Neutral loss search

```sql
SELECT * FROM massql(
  'QUERY scaninfo(MS2DATA) WHERE MS2NL=176.032',
  'sample.mzML'
);
```

### Retention time window

```sql
SELECT * FROM massql(
  'QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857 AND RTMIN=1.0 AND RTMAX=5.0',
  'sample.mzML'
);
```

### Exclude spectra with a contaminant peak

```sql
SELECT * FROM massql(
  'QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857 AND MS2PROD=391.28:EXCLUDED',
  'sample.mzML'
);
```

### Count matching spectra

```sql
-- DuckDB SQL tip: COUNT(*) counts rows in the result
SELECT COUNT(*) FROM massql(
  'QUERY scannum(MS2DATA) WHERE MS2PROD=167.0857:TOLERANCEPPM=5',
  'sample.mzML'
);
```

### Query per sample

When your table contains data from multiple samples, use `sample_id` to run the MassQL query independently for each sample:

```sql
-- Load multiple files with filepath as sample identifier
CREATE TABLE all_spectra AS
  SELECT * FROM read_mzml('data/*.mzML', include_filepath=true);

-- Run the query per file — results include 'filepath' as the first column
SELECT * FROM massql(
  'QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857',
  'all_spectra',
  sample_id := 'filepath'
);

-- Count matching spectra per sample
SELECT filepath, COUNT(*) FROM massql(
  'QUERY scannum(MS2DATA) WHERE MS2PROD=167.0857',
  'all_spectra',
  sample_id := 'filepath'
) GROUP BY filepath;
```

The `sample_id` column must not contain NULL values. It works with any column type (VARCHAR, INTEGER, etc.) and cannot be used with file path sources — load the file into a table first.

### Work with multiple files

```sql
-- Load multiple files into a table, then query
CREATE TABLE all_spectra AS
  SELECT * FROM read_mzml('data/*.mzML', include_filepath=true);

SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857', 'all_spectra');
```

### Export results to CSV

```sql
-- DuckDB SQL tip: COPY ... TO exports results to a file
COPY (
  SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=167.0857', 'sample.mzML')
) TO 'results.csv' (HEADER, DELIMITER ',');
```

---

## For SQL Users: Reading Mass Spectrometry Files Directly

If you prefer to write SQL rather than MassQL, you can query mzML/mzXML files directly:

```sql
-- Read an mzML file (returns one row per spectrum with 27 columns)
SELECT * FROM read_mzml('sample.mzML');

-- Read an mzXML file (identical output schema, enables UNION ALL across formats)
SELECT * FROM read_mzxml('sample.mzXML');

-- Read chromatograms from mzML
SELECT * FROM read_mzml_chromatograms('sample.mzML');
```

See [Reading files](reading.md#mzml-and-mzxml) for the full column schema.

### Useful SQL Macros for Mass Spec

These macros are used internally by the MassQL engine but are available for direct SQL use:

```sql
-- Unnest peak arrays into individual rows with normalized intensity
SELECT * FROM mzml_peaks('my_table') WHERE ms_level = 2;

-- m/z tolerance checks
SELECT mz_within(200.05, 200.0, 0.1);     -- true (within 0.1 Da)
SELECT mz_within_ppm(200.05, 200.0, 10);  -- true (within 10 ppm)

-- Scan-level summaries
SELECT * FROM mzml_scaninfo('my_table');   -- one row per spectrum
SELECT * FROM mzml_scansum('my_table');    -- intensity sums
SELECT * FROM mzml_scannum('my_table');    -- distinct scan indices
```

#### `mzml_peak_pair(relation, formula_str)`

Finds MS2 spectra that contain a **peak pair** — one peak at m/z = X and another at m/z = `2*X - formula(formula_str)`, within 0.1 Da — and returns all peaks from the matching spectra (composable with `mzml_scaninfo`). Useful for isotope/adduct-offset patterns keyed to an element or group.

- `relation` (VARCHAR): a table/view name holding [`read_mzml`](reading.md#mzml-and-mzxml) output.
- `formula_str` (VARCHAR): a chemical formula (e.g. `'Fe'`); its monoisotopic mass is computed via [`formula()`](#chemical-formula-functions).

Output columns mirror the input relation's schema.

```sql
-- MS2 spectra with an iron-offset peak pair
SELECT * FROM mzml_peak_pair('spectra', 'Fe');
```

> **Note:** `massql.md` documents the most commonly used mass-spec helper macros. The full set (`mzml_ms`, `mzml_filter_mz`/`_nl`, `mzml_x_*`, `mzml_isotope_pattern`, `massdefect`, `mz_massdefect_within`, `mzml_i_norm`, `mzml_or_cardinality`, `mzml_excluded_ms`, …) is registered and callable; a complete helper-macro reference is a tracked follow-up.
