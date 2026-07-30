# Extension updating 
When cloning this template, the target version of DuckDB should be the latest stable release of DuckDB. However, there 
will inevitably come a time when a new DuckDB is released and the extension repository needs updating. This process goes
as follows:

**The mechanical part is scripted.** `./scripts/bump-duckdb-version.sh vX.Y.Z` moves both submodules and rewrites every in-repo spelling below, then verifies with `./scripts/check-duckdb-version.sh` (which also runs in CI as the `DuckDB version consistency` job). It deliberately does not commit — a version bump can break the C++ build, so build and test before committing. The rest of this page is what those scripts do, plus the cross-repo steps they cannot do for you.

- Bump submodules
  - `./duckdb` should be set to latest tagged release
  - `./extension-ci-tools` should be set to updated branch corresponding to latest DuckDB release. We track the release-series branch (e.g. `v1.5-variegata`), which upstream keeps identical to the per-patch branch (`v1.5.5`), so check out that branch's head.
- Bump versions in `./github/workflows/MainDistributionPipeline.yml`
  - `duckdb_version` input in the `duckdb-stable-build` job
  - `duckdb_version` input in the `code-quality-check` job
  - the commented-out artifact name in the disabled `verify-wasm` job (`miint-<ver>-extension-...`), so it isn't stale when that job is re-enabled
  - `ci_tools_version` and the `uses:` refs need NO change — they point at the moving release-series branch, not a per-patch tag
- Bump the target version everywhere else it is spelled out
  - `Dockerfile` `ARG DUCKDB_VERSION` (local default; `docker.yml` resolves the newest Docker Hub tag at runtime and passes `--build-arg`)
  - `DUCKDB_VERSION` default in `scripts/cron-publish-extension.sh` — it drives the `miint-<ver>-extension-*` artifact prefix and the published paths, so a stale value makes every publish reject
  - `duckdb==` pin in `python/pyproject.toml` — the CLI `LOAD`s the extension, so its DuckDB must match exactly

Outside this repo (do these after CI is green and a release tag is pushed):

- The cron host's env file (e.g. `/etc/miint-publish.env`) overrides `DUCKDB_VERSION`; editing the script alone changes nothing there. The live FTP paths become `<ver>/` and `tagged/<ver>/`, leaving the previous version's trees serving until retired.
- The playground console in the umbrella site repo pins both `@duckdb/duckdb-wasm` and `EXT_VERSION`. A C++ extension only loads into a duckdb-wasm of the *same* core version, so bump both together and verify the npm build actually reports the new version before flipping.
- `duckdb/community-extensions` `extensions/miint/description.yml` carries the `ref` that `INSTALL miint FROM community` serves; open a PR bumping it to the new tag.

# API changes
DuckDB extensions built with this extension template are built against the internal C++ API of DuckDB. This API is not guaranteed to be stable.
What this means for extension development is that when updating your extensions DuckDB target version using the above steps, you may run into the fact that your extension no longer builds properly.

Currently, DuckDB does not (yet) provide a specific change log for these API changes, but it is generally not too hard to figure out what has changed.

For figuring out how and why the C++ API changed, we recommend using the following resources:
- DuckDB's [Release Notes](https://github.com/duckdb/duckdb/releases)
- DuckDB's history of [Core extension patches](https://github.com/duckdb/duckdb/commits/main/.github/patches/extensions)
- The git history of the relevant C++ Header file of the API that has changed