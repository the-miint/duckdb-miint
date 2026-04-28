---
title: Building the docs
description: How to configure, generate, and serve the duckdb-miint documentation site locally.
---

The docs site is an Astro Starlight project under `site/`. Reference pages
for each registered function are generated at build time by introspecting
DuckDB's function catalog (`duckdb_functions()`) and rendering one MD page
per function. Hand-authored prose (this page, the architecture doc, etc.)
lives in `site/src/content/docs/`.

## Two introspection modes

By default, the docs build introspects the **released miint extension**
served from the public repository at <https://ftp.microbio.me/pub/miint>.
This means the generated reference reflects what users actually install,
not whatever happens to be in your local `build/` directory. The build
needs network access to that URL but no local C++ build.

For active development on a registration (you've edited a `RegisterDocumented*`
call and want to preview the new docs before releasing), set `MIINT_EXT`
to the path of a locally-built extension. The introspector then loads
that artifact directly and skips the released-repo path.

```sh
# Default: docs reflect the released extension.
npm run prebuild

# Local build: preview a registration you just edited.
MIINT_EXT=/path/to/build/release/extension/miint/miint.duckdb_extension \
  npm run prebuild
```

The local-build path uses `duckdb -unsigned` because locally-built
extensions aren't signed; the released path doesn't need it.

## Prerequisites

| Tool | Why | Where it comes from |
|---|---|---|
| Stock DuckDB CLI ≥ the version the released miint was built against | Used by the introspector to get a baseline `duckdb_functions()` snapshot, then diff against the post-LOAD catalog | Install from your package manager or download from <https://duckdb.org/docs/installation/>. The released miint repo only carries artifacts for specific DuckDB versions (currently v1.5.2+); a mismatched stock client gets a 404. Override the binary with `DUCKDB_BIN`. |
| Node 20+ | Astro toolchain | Bundled with [nvm](https://github.com/nvm-sh/nvm) or any Node 20+ install. `npm` ships with it. |
| (local-build mode only) The standard miint build env | To produce `build/release/extension/miint/miint.duckdb_extension` | See [Installation](/getting-started/installation/) — recommends `conda create -n duckdb-build -c conda-forge gcc gxx cmake ninja make perl` |

The introspector cannot use the repo's `build/release/duckdb` because that
binary auto-loads miint at startup, so a "before LOAD" snapshot is
impossible. Use a stock build of the matching DuckDB version.

The introspector always invokes `duckdb -unsigned`. The local-build path
needs it because locally-built extensions aren't signed; the released
path needs it because the dedicated miint repository
(`ftp.microbio.me/pub/miint`) is not yet signed by a key that stock
DuckDB recognizes for that URL — only the canonical
`community-extensions.duckdb.org` path is. Once miint releases are
mirrored into a signed channel this can be dropped.

## Steps (released-extension mode)

```sh
# 1. Install site dependencies (once per package-lock.json change)
cd site && npm install

# 2. Generate reference + sidebar manifest from the released extension
npm run prebuild

# 3a. Static build to site/dist/
npm run build

# 3b. Or dev server with live reload (re-runs prebuild on each restart)
npm run dev
```

## Steps (local-build mode)

```sh
# 1. Build the extension (once per C++ change)
bash build.sh

# 2. Install site dependencies (once per package-lock.json change)
cd site && npm install

# 3. Generate against the local build
MIINT_EXT=$(pwd)/../build/release/extension/miint/miint.duckdb_extension \
  npm run prebuild

# 4. Build or serve as above
npm run build   # or: npm run dev
```

`npm run prebuild` runs two scripts in sequence:

- `scripts/introspect-functions.mjs` — diffs `duckdb_functions()` before
  vs after the extension is loaded, emits `site/build/functions.json`.
- `scripts/generate-reference.mjs` — reads `functions.json`, writes one
  MD per function under `src/content/docs/reference/<type>/<category>/<name>.md`,
  one `index.md` per category, and `site/build/sidebar.json` consumed by
  `astro.config.mjs` to render the Reference sidebar.

Both scripts are idempotent — they wipe their output dirs first, so
removing a function from the C++ removes its page on the next run.

## Environment variables

| Variable | Default | Purpose |
|---|---|---|
| `DUCKDB_BIN` | `duckdb` (PATH lookup) | Stock duckdb client of the matching version |
| `MIINT_REPO` | `https://ftp.microbio.me/pub/miint` | Custom miint extension repository to install from |
| `MIINT_EXT` | _unset_ (use `MIINT_REPO`) | Path to a local `.duckdb_extension`. When set, `MIINT_REPO` is ignored and `duckdb` is invoked with `-unsigned`. |

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
   `site/scripts/generate-reference.mjs` (this allowlist exists so the
   POC only generates pages for functions whose registrations have been
   converted; remove the allowlist when every registration is converted).
4. Run `MIINT_EXT=... npm run prebuild` (local-build mode) to preview the
   new page. Once the change ships in a release, the default build picks
   it up automatically.

## What the generator does NOT (yet) read from the catalog

`duckdb_functions()` exposes description, parameter names, parameter
types, return type, examples, categories, and `alias_of`. It does **not**
expose named-parameter defaults for table functions — those live inside
each Bind callback. As a workaround, document defaults in the description
prose (the field accepts full markdown, so a `### Named parameters`
table renders fine). A small `miint_function_metadata()` debug table
function backed by an in-process registry is the proposed Phase B fix.
