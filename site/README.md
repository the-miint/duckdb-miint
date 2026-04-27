# duckdb-miint docs

The Astro Starlight project that builds the duckdb-miint documentation site.

Reference pages for each registered miint function are auto-generated at
build time by introspecting the compiled extension. Narrative pages
(architecture, embedded tools, this README, etc.) are hand-authored.

## Quick start

```sh
bash ../build.sh                          # build the extension once
cd site && npm install                    # install deps once
npm run dev                               # opens http://127.0.0.1:4321/
```

Full procedure, environment variables, and the rationale for the two-source
pipeline: see [`src/content/docs/internals/building-the-docs.md`](src/content/docs/internals/building-the-docs.md).
