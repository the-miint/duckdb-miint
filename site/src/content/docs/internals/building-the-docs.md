---
title: Building the docs
description: How to configure, generate, and serve the duckdb-miint documentation site locally.
---

The docs site is an Astro Starlight project under `site/`. Reference pages
for each registered function are generated at build time by introspecting the
compiled extension via `duckdb_functions()` and rendering one MD page per
function. Hand-authored prose (this page, the architecture doc, etc.) lives in
`site/src/content/docs/`.

The end-to-end flow is **build the extension once, then build the site**.

## Prerequisites

| Tool | Why | Where it comes from |
|---|---|---|
| The standard miint build env | To produce `build/release/extension/miint/miint.duckdb_extension` | See `docs/installation.md` (recommends `conda create -n duckdb-build -c conda-forge gcc gxx cmake ninja make perl`) |
| Stock DuckDB CLI matching the extension's version | Used by the introspector to get a baseline `duckdb_functions()` snapshot before `LOAD miint` | Default path is `/home/lpatel/software/miniconda3/envs/miint/bin/duckdb`. Override with `DUCKDB_BIN`. |
| Node 20+ | Astro toolchain | Bundled with [nvm](https://github.com/nvm-sh/nvm) or any Node 20+ install. `npm` ships with it. |

The introspector cannot use the repo's `build/release/duckdb` because that
binary auto-loads miint at startup, so a "before LOAD" snapshot is impossible.
Use a stock build of the matching DuckDB version.

## Steps

```sh
# 1. Build the extension (once per C++ change)
bash build.sh

# 2. Install site dependencies (once per package-lock.json change)
cd site && npm install

# 3. Generate reference + sidebar manifest from the built extension
npm run prebuild

# 4a. Static build to site/dist/
npm run build

# 4b. Or dev server with live reload (re-runs prebuild on each restart)
npm run dev
```

`npm run prebuild` runs two scripts in sequence:

- `scripts/introspect-functions.mjs` — diffs `duckdb_functions()` before vs after
  `LOAD`, emits `site/build/functions.json`.
- `scripts/generate-reference.mjs` — reads `functions.json`, writes one MD per
  function under `src/content/docs/reference/<type>/<category>/<name>.md`,
  one `index.md` per category, and `site/build/sidebar.json` consumed by
  `astro.config.mjs` to render the Reference sidebar.

Both scripts are idempotent — they wipe their output dirs first, so removing a
function from the C++ removes its page on the next run.

## Environment variables

| Variable | Default | Purpose |
|---|---|---|
| `DUCKDB_BIN` | `/home/lpatel/software/miniconda3/envs/miint/bin/duckdb` | Stock duckdb client of the matching version |
| `MIINT_EXT` | `build/release/extension/miint/miint.duckdb_extension` | Path to the built extension to introspect |

## Viewing the site over SSH

The dev server binds to `127.0.0.1:4321`. To view it from a workstation:

```sh
# On the workstation
ssh -N -L 4321:localhost:4321 <remote-host>
# Then open http://localhost:4321/
```

## Adding a new documented function

1. In the C++ registration site, switch from
   `loader.RegisterFunction(myFunction)` to
   `RegisterDocumentedScalar(...)` /
   `RegisterDocumentedScalarSet(...)` /
   `RegisterDocumentedTableFunction(...)`
   from `src/include/documented_function.hpp`. Pass description, parameter
   names, examples, optional `alias_of`, optional category.
2. Rebuild the extension.
3. Add the function name to `POC_SLICE` in
   `site/scripts/generate-reference.mjs` (this allowlist exists so the POC
   only generates pages for functions whose registrations have been
   converted; remove the allowlist when every registration is converted).
4. Re-run `npm run prebuild`. The new page appears under the matching
   category in the sidebar.

## What the generator does NOT (yet) read from the catalog

`duckdb_functions()` exposes description, parameter names, parameter types,
return type, examples, categories, and `alias_of`. It does **not** expose
named-parameter defaults for table functions — those live inside each Bind
callback. As a workaround, document defaults in the description prose (the
field accepts full markdown, so a `### Named parameters` table renders fine).
A small `miint_function_metadata()` debug table function backed by an
in-process registry is the proposed Phase B fix.
