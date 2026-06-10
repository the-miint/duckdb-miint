# Overlay port: build zlib-ng in ZLIB_COMPAT mode and install it as the `zlib`
# package. Because it is an overlay port named "zlib", it shadows the upstream
# zlib port for every consumer in the dependency graph (the manifest, plus
# transitive deps like HDF5 and RocksDB) and for kseq++/htslib via
# find_package(ZLIB). Result: one zlib in the build, and it is zlib-ng's
# SIMD-accelerated inflate -- the hot path for reading gzipped FASTQ/BAM.
#
# REF/SHA512 copied verbatim from vcpkg's upstream zlib-ng port (2.2.5). The
# package version (1.3.1, see vcpkg.json) is intentionally decoupled from the
# zlib-ng source tag: it advertises the zlib API level zlib-ng-compat emulates
# so dependents that pin "zlib >= 1.x" still resolve.
vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO zlib-ng/zlib-ng
    REF 2.2.5
    SHA512 b599ea24375d08fa098ed7c3b14548e0d9731a155a024a0904b0ae4a6d3491a69f0c0574d66b6e4af1e40f10e38b6b555d4c4b1fe3589ca83a5f97fbd92f635f
    HEAD_REF develop
)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        "-DZLIB_FULL_VERSION=1.3.1"
        -DZLIB_ENABLE_TESTS=OFF
        -DWITH_NEW_STRATEGIES=ON
        -DZLIB_COMPAT=ON
    OPTIONS_RELEASE
        -DWITH_OPTIM=ON
)
vcpkg_cmake_install()
vcpkg_copy_pdbs()

vcpkg_fixup_pkgconfig()

# In ZLIB_COMPAT mode zlib-ng installs its CMake config under lib/cmake/ZLIB,
# providing the standard ZLIB::ZLIB target that find_package(ZLIB) consumes.
vcpkg_cmake_config_fixup(CONFIG_PATH lib/cmake/ZLIB)

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/share"
                    "${CURRENT_PACKAGES_DIR}/debug/include")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE.md")
