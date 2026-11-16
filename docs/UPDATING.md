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
  - `ci_tools_version` and the `uses:` refs need NO change within a release series — they point at the moving release-series branch, not a per-patch tag. Moving to a new series (e.g. `v1.5-variegata` → `v2.0-cyanoptera`) changes both, together with the `./extension-ci-tools` branch above, and the `BRANCH` default and example env file in `scripts/cron-publish-extension.sh` when the new series' branch becomes the one published, and the branch table and descriptor keys in [`internals/releasing.md`](internals/releasing.md#branches)
  - Before DuckDB tags a release, the workflow can build against the duckdb submodule's commit instead: add a job that calls extension-ci-tools' `_submodule_version.yml` and pass its `version` output as both `duckdb_version` inputs. `check-duckdb-version.sh` then reports PRE-RELEASE and skips the other spellings, and the artifacts carry the commit hash, which no released DuckDB loads and the publish cron rejects. On release day put the literal version back in both inputs, drop the job, and run the bump script
- Bump the target version everywhere else it is spelled out
  - `Dockerfile` `ARG DUCKDB_VERSION` (local default; `docker.yml` resolves the newest Docker Hub tag at runtime and passes `--build-arg`)
  - `DUCKDB_VERSION` default in `scripts/cron-publish-extension.sh` — it drives the `miint-<ver>-extension-*` artifact prefix and the published paths, so a stale value makes every publish reject
  - `duckdb==` pin in `python/pyproject.toml` — the CLI `LOAD`s the extension, so its DuckDB must match exactly
  - `@duckdb/duckdb-wasm` pin in `scripts/wasm_browser/index.html`, the runtime of the browser check (`scripts/test_wasm_browser.sh`). Its npm version does not follow DuckDB's numbering: pick the release that reports the new core in `pragma_version()`, the same one the playground moves to (its pin is in [`docs/playground.md`](https://github.com/the-miint/the-miint.github.io/blob/main/docs/playground.md) of the-miint.github.io). The scripts do not rewrite it; the browser check fails, naming the file, until it matches.

Outside this repo (do these after CI is green and a release tag is pushed):

- The cron host's env file (e.g. `/etc/miint-publish.env`) overrides `DUCKDB_VERSION`, and `BRANCH` on a series move; editing the script alone changes nothing there. The live FTP paths become `<ver>/` and `tagged/<ver>/`, leaving the previous version's trees serving until retired.
- The playground console in the umbrella site repo pins both `@duckdb/duckdb-wasm` and `EXT_VERSION`. A C++ extension only loads into a duckdb-wasm of the *same* core version, so bump both together and verify the npm build actually reports the new version before flipping.
- `duckdb/community-extensions` `extensions/miint/description.yml` carries the `ref` that `INSTALL miint FROM community` serves; it moves only with a release, by [`internals/releasing.md`](internals/releasing.md).

Traps when verifying a bump:

- `scripts/test_wasm_load.sh` reuses `build/wasm_shell/` whenever it holds a core, and stops when that core's recorded DuckDB version is not the extension's or none is recorded. After a bump, run it with `--clean`.
- CI runs no WASM check on a pull request or branch push: run Levels 2–4 of [`wasm-testing.md`](wasm-testing.md) locally before tagging.
- DuckDB release tags are lightweight, so `git ls-remote <url> 'refs/tags/vX.Y.Z^{}'` matches nothing; query the plain ref as well, as `scripts/check-duckdb-version.sh` does.

# API changes
DuckDB extensions built with this extension template are built against the internal C++ API of DuckDB. This API is not guaranteed to be stable.
What this means for extension development is that when updating your extensions DuckDB target version using the above steps, you may run into the fact that your extension no longer builds properly.

Currently, DuckDB does not (yet) provide a specific change log for these API changes, but it is generally not too hard to figure out what has changed.

For figuring out how and why the C++ API changed, we recommend using the following resources:
- DuckDB's [Release Notes](https://github.com/duckdb/duckdb/releases)
- DuckDB's history of [Core extension patches](https://github.com/duckdb/duckdb/commits/main/.github/patches/extensions)
- The git history of the relevant C++ Header file of the API that has changed