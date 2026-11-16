# Releasing

How a duckdb-miint release is numbered, cut, and delivered to `INSTALL miint FROM community`. Moving to a new DuckDB version is a separate procedure: [`../UPDATING.md`](../UPDATING.md).

## Versions and tags

Versions are semver on the SQL surface, independent of DuckDB's numbering:

| Bump | When |
|---|---|
| Major | A released function is removed or renamed, its parameters, output columns or types change, or a changed default changes results. |
| Minor | New functions or features, and moving to a new DuckDB minor version. |
| Patch | Fixes only, build fixes included. |

A release is an annotated tag, `vX.Y.Z`, preceded by one or more `vX.Y.Z-rc.N` on the commit it will name. The tag message is the release's summary.

A DuckDB 2.0 build of a tagged commit reports the tag as its version (`miint_version()`, `extension_version` in `duckdb_extensions()`), found with `git describe --tags --match "v*.*.*"`. It reports a commit hash otherwise. `git describe` prefers an annotated tag over a lightweight one on the same commit, and the newer of two annotated tags, so tag the final after the candidate, and never lightweight. DuckDB 1.5 builds always report a hash: its build passes the pattern with literal quotes, fixed for 2.0 by [duckdb/duckdb#24083](https://github.com/duckdb/duckdb/pull/24083).

## Branches

| Branch | Releases | Descriptor key |
|---|---|---|
| The default branch, for DuckDB's latest stable release (`v2.0-cyanoptera`) | every release | `ref` |
| The previous DuckDB line (`v1.5-variegata`) | frozen at `v1.0.1`; a patch release only on demand | the line's codename (`variegata`) |

When DuckDB opens its next release's feature freeze (a `vX.Y-<codename>` branch upstream), create the matching branch and set `ref_next` in the descriptor by hand, as community-extensions' [`UPDATING.md`](https://github.com/duckdb/community-extensions/blob/main/UPDATING.md) describes. At DuckDB's release that branch becomes the default, and the old one freezes under its codename key.

## How community-extensions builds miint

- The descriptor, `extensions/miint/description.yml` in [duckdb/community-extensions](https://github.com/duckdb/community-extensions), names one commit per DuckDB line (`ref`, `variegata`, ...). `extension.version` is a single field shared by all of them, so it follows `ref`.
- A merged descriptor change rebuilds only the keys it changed, on every platform, and publishes the result over the previous build. There are no versioned installs.
- At every DuckDB release, community-extensions rebuilds every extension from `ref` without asking. So `ref` only ever names a release tag.
- DuckDB staff review every descriptor PR, so release when there is something to ship, at most about monthly, and not per merge.

## Cutting a release

1. Pick a commit on the release branch whose push CI passed. Before a final tag, compare its SQL surface with the previous release to check the version bump: function names, parameters and return types, from `duckdb_functions()` after `LOAD` minus before. Table functions' output columns are not in that view.
2. Tag the candidate and push the tag:
   ```bash
   git tag -a v1.2.0-rc.1 -m "Release v1.2.0-rc.1 ..." <commit>
   git push origin v1.2.0-rc.1
   ```
   The tag push builds every platform, wasm_eh included, which branch pushes and pull requests do not.
3. When that run passes, propose it:
   ```bash
   ./scripts/release-community-extension.sh v1.2.0-rc.1
   ```
   This opens a draft PR from your fork. A fix means a new commit, `-rc.2` on it, and the same command, which pushes to the same PR.
4. When the tag's CI and the PR's CI have both passed, tag the final on the same commit, push it, and wait for its CI:
   ```bash
   git tag -a v1.2.0 -m "Release v1.2.0 ..." <same commit>
   git push origin v1.2.0
   ./scripts/release-community-extension.sh v1.2.0
   ```
   The script marks the PR ready without pushing to it: the build is the one already tested.
5. Publish the release notes:
   ```bash
   gh release create v1.2.0 --verify-tag --generate-notes --notes-start-tag <previous release>
   ```
6. The publish cron (`scripts/cron-publish-extension.sh`) puts the tag's builds, wasm_eh included, under `tagged/<duckdb version>/` on its next run.

A patch release on a frozen line names its key: `./scripts/release-community-extension.sh v1.0.2 variegata`.

## What the script refuses

`scripts/release-community-extension.sh` changes nothing unless the tag is an annotated `vX.Y.Z` or `vX.Y.Z-rc.N`, is on origin, names a commit on the key's release branch, and its tag push's CI run passed. `--dry-run` runs every check and prints the descriptor change. It never switches or edits your community-extensions checkout: it commits in a temporary worktree. It expects that checkout at `../community-extensions` (or `COMMUNITY_EXTENSIONS_DIR`), with your fork as `origin` and duckdb/community-extensions as `upstream`. `test/shell/release_community_extension.sh` covers each refusal and each PR state.
