export VCPKG_TOOLCHAIN_PATH=`pwd`/vcpkg/scripts/buildsystems/vcpkg.cmake

# Keep the local vcpkg clone in sync with the baseline pinned in vcpkg.json.
# vcpkg here is a plain clone (not a submodule), so a bumped builtin-baseline
# won't resolve unless this checkout's version database contains that commit.
# CI pins the same SHA via lukka/run-vcpkg (vcpkgGitCommitId); this mirrors that
# for local builds. We only move the checkout when the baseline isn't already
# reachable from HEAD, so a newer vcpkg is never downgraded.
if [ -d vcpkg/.git ]; then
    BASELINE=$(grep builtin-baseline vcpkg.json 2>/dev/null | grep -oE '[0-9a-f]{40}' | head -1)
    if [ -n "$BASELINE" ] && ! git -C vcpkg merge-base --is-ancestor "$BASELINE" HEAD 2>/dev/null; then
        echo "vcpkg: syncing checkout to pinned baseline $BASELINE"
        git -C vcpkg cat-file -e "${BASELINE}^{commit}" 2>/dev/null \
            || git -C vcpkg fetch --quiet origin "$BASELINE" 2>/dev/null \
            || git -C vcpkg fetch --quiet origin
        if [ -n "$(git -C vcpkg status --porcelain)" ]; then
            echo "ERROR: vcpkg working tree has local changes; refusing to move it to $BASELINE." >&2
            echo "       Commit/stash inside vcpkg/ or check out the baseline manually, then rerun." >&2
            exit 1
        fi
        git -C vcpkg checkout --quiet "$BASELINE" \
            || { echo "ERROR: could not check out vcpkg baseline $BASELINE" >&2; exit 1; }
    fi
fi

GEN=ninja make
