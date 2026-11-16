# Releasing

How a duckdb-miint release is numbered, cut, and delivered to `INSTALL miint FROM community`. Moving to a new DuckDB release series is a separate procedure: [`docs/UPDATING.md`](https://github.com/the-miint/duckdb-miint/blob/HEAD/docs/UPDATING.md).

## How community-extensions builds miint

- The descriptor, `extensions/miint/description.yml` in [duckdb/community-extensions](https://github.com/duckdb/community-extensions), names one miint commit per DuckDB release series, under a key: `ref` for DuckDB's latest stable release, the series' codename (`variegata`, ...) for an older one, `ref_next` for an upcoming one. `extension.version` is a single field shared by all of them, so it follows `ref`.
- A merged descriptor change rebuilds only the keys it changed, on every platform the descriptor allows (wasm_mvp and wasm_threads among them, which miint's CI never builds), and publishes the result over the previous build. There are no versioned installs.
- At every DuckDB release, community-extensions rebuilds every extension from `ref` without asking. So `ref` only ever names a release tag.
- DuckDB staff review every descriptor PR, so release when there is something to ship, at most about monthly, and not per merge.

## Versions and tags

Versions are semver on the SQL surface, independent of DuckDB's numbering:

| Bump | When |
|---|---|
| Major | A released function is removed or renamed, its parameters, output columns or types change, or a changed default changes results. |
| Minor | New functions or features, and moving to a new DuckDB release series (1.5 → 2.0 included). |
| Patch | Fixes only, build fixes included. |

A release is an annotated tag, `vX.Y.Z`, preceded by one or more `vX.Y.Z-rc.N` on the commit it will name. The tag message is the release's summary.

A DuckDB 2.0 build of a tagged commit reports the tag as its version (`miint_version()`, `extension_version` in `duckdb_extensions()`), found with `git describe --tags --match "v*.*.*"`. It reports a commit hash otherwise. Among the release tags on the commit, `git describe` prefers an annotated tag to a lightweight one, then the newest, so tag the final after the candidate. DuckDB 1.5 builds always report a hash: its build passes the pattern with literal quotes, fixed for 2.0 by [duckdb/duckdb#24083](https://github.com/duckdb/duckdb/pull/24083).

## Branches

| Branch | Releases | Descriptor key |
|---|---|---|
| The default branch, for DuckDB's latest stable release (`v2.0-cyanoptera`) | every release | `ref` |
| An older series' branch (`v1.5-variegata`) | frozen: a patch release only on demand | the series' codename (`variegata`) |

When DuckDB opens its next release's feature freeze (a `vX.Y-<codename>` branch upstream), create the matching branch and set `ref_next` in the descriptor, by hand, to a release candidate's commit, as community-extensions' [`UPDATING.md`](https://github.com/duckdb/community-extensions/blob/main/UPDATING.md) describes. At DuckDB's release that branch becomes the default and its candidate's final tag goes to `ref`; add the old series' codename key, by hand, naming the commit of its last release tag. The script then maintains both.

## Cutting a release

1. Pick a commit on the release branch whose push CI passed. Before a final tag, compare its SQL surface with the previous release to check the version bump: function names, parameters and return types, from `duckdb_functions()` after `LOAD` minus before. Table functions' output columns are not in that view.
2. Tag the candidate and push the tag:
   ```bash
   git tag -a v1.2.0-rc.1 -m "Release v1.2.0-rc.1 ..." <commit>
   git push origin v1.2.0-rc.1
   ```
   The tag push adds wasm_eh to the platforms CI builds, which branch pushes and pull requests skip.
3. When that run passes, propose it:
   ```bash
   ./scripts/release-community-extension.sh v1.2.0-rc.1
   ```
   This opens a draft PR from your fork. A fix means a new commit, `-rc.2` on it, and the same command, which pushes to the same PR.
4. When the PR's own CI has passed, tag the final on the same commit and push it:
   ```bash
   git tag -a v1.2.0 -m "Release v1.2.0 ..." <same commit>
   git push origin v1.2.0
   ```
5. When the final tag's CI run has passed:
   ```bash
   ./scripts/release-community-extension.sh v1.2.0
   ```
   The script marks the PR ready, without pushing to it, once community-extensions has built the commit the PR names and passed (`build.yml` for `ref`, `build_<codename>.yml` for a frozen series); until then it says so and changes nothing.
6. Publish the release notes:
   ```bash
   gh release create v1.2.0 --verify-tag --generate-notes --notes-start-tag <previous release>
   ```
7. The publish cron (`scripts/cron-publish-extension.sh`) puts the tag's builds, wasm_eh included, under `tagged/<duckdb version>/` on its next run.

A patch release on a frozen series names its key, `./scripts/release-community-extension.sh v1.0.2 variegata`, and needs no candidate: the PR opens as a draft, and running the same command again once its CI has passed marks it ready.

## What the script refuses, and why

`scripts/release-community-extension.sh` proposes nothing when:

| Refused | Why |
|---|---|
| The argument is not `vX.Y.Z` or `vX.Y.Z-rc.N` | The descriptor names only tagged commits, the only ones whose build reports a version. |
| The tag is lightweight | The tag message is the release's summary, and `git describe` ranks any annotated tag on the commit above it. |
| A build of the commit would report another tag (`git describe` over origin's release tags on the commit, as the build runs it) | The binary would report that tag as its version, while the PR (and, for `ref`, `extension.version`) names this one. |
| The tag is not on origin, or origin is not the descriptor's `repo.github` | community-extensions checks the commit out from `repo.github`. |
| A local copy of the tag, or of another release tag on the commit, differs from origin's | The checks would read another commit, or `git describe` another set of tags, than community-extensions sees. |
| The commit is not on the key's release branch | Only merged code ships. |
| The tag push's CI run is missing, has not passed, or tested another commit | Only a tag push builds wasm_eh. |
| The key does not already hold a commit | So `repo.github` or an `extension:` field is never rewritten. |
| `origin` is not your fork of duckdb/community-extensions, or `upstream` is not duckdb/community-extensions | The change is proposed against upstream's descriptor, from a branch on your fork. |
| Your fork has the release's branch but no open PR from it | It holds an unfinished attempt on an old upstream; the script prints the command that deletes it. |

On the PR it proposes, every new PR opens as a draft, a ready PR goes back to draft before anything is pushed to it, a push rewrites the PR's description to name the new tag, and only a final tag, once community-extensions has built the commit the fork branch names and passed, marks it ready. `--dry-run` runs every check up to the descriptor change and prints the change (or that there is nothing to push) and the PR's state, pushing nothing and changing no PR; like every run, it fetches origin's release tags on the commit and the remotes' branches. The script never switches or edits your community-extensions checkout: it commits in a temporary worktree. It expects that checkout at `../community-extensions` (or `COMMUNITY_EXTENSIONS_DIR`). `test/shell/release_community_extension.sh` has a scenario for each refusal in this table and each PR transition.
