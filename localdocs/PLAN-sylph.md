# Sylph Profile — TDD Implementation Plan

> **Status:** PROPOSED. Not yet started. This is the working spec; update as decisions land.

## Context

Embed [sylph](https://github.com/bluenote-1577/sylph) into duckdb-miint as a table function `sylph_profile(...)` that performs FracMinHash-based relative-abundance profiling of paired-end shotgun metagenomic reads against a pre-built `.syldb` reference database.

The sketch & profile compute runs **in-process, in-memory** against DuckDB-resident sequence data. The `.syldb` is loaded once per session and reused across many calls. Output is delivered zero-copy via Arrow C Data Interface.

User clone of upstream sylph 0.9.0 lives at `../sylph/` (sibling of `duckdb-miint/`). 0.9.0 is the latest upstream release as of writing. It will be set up as a git submodule at `ext/sylph/` pointing at `the-miint/sylph` (a fork to be created), branch `v0.9.0-miint`.

### Goals

1. SQL table function `sylph_profile(query_table, syldb_path)` returning per-genome abundance rows.
2. **Zero file I/O for sample reads** — sequences come from a DuckDB table/view, not FASTQ paths.
3. `.syldb` loaded once per session, mmap-backed, shared across queries.
4. Per-sample partitioning via the existing `per_sample_table_function` helper (mirrors `align_mafft`, `detect_chimera_uchime_denovo`).
5. Bit-identical numerical output to upstream `sylph profile --reads` on the same data, modulo documented divergence.
6. Cleanly disabled on WASM and Windows (sylph fork must compile elsewhere; integration is Linux/macOS first).

### Non-goals (v1)

- Sketch *building* (`sylph sketch`) from FFI — keep that as the upstream CLI; users prepare `.syldb` offline.
- `sylph contain` (containment-only without abundance). Will be a follow-up `sylph_contain` table function once the FFI shape is settled.
- Sketching *reference genomes* in DuckDB. References are already-sketched `.syldb` only.
- Long-read profiling tweaks. Short reads only in v1; long-read params can be added later.
- Any networking. `.syldb` must be a local file path.

---

## Architecture

```
ext/sylph/                                  — Rust fork of bluenote-1577/sylph
├── src/
│   ├── lib.rs                              — exports existing modules + new `profile_api`, `c_api`
│   ├── sketch.rs                           — UNCHANGED logic; new `sketch_pair_slices()` added
│   ├── contain.rs                          — refactored: I/O split from compute
│   ├── profile_api.rs                      — NEW: pure data-in/data-out profile entry point
│   └── c_api.rs                            — NEW: C/Arrow FFI surface (~600 LOC)
├── sylph.h                                  — NEW: hand-written C header (the contract)
├── Cargo.toml                              — add `[lib] crate-type = ["staticlib","rlib"]`,
│                                              `arrow` dep behind `arrow-ffi` feature,
│                                              move `clap` behind `cli` feature,
│                                              drop `tikv-jemallocator`,
│                                              drop `panic = "abort"` from release profile
└── Cargo.lock

duckdb-miint/
├── ext/sylph/                              — moved from ../sylph/, then forked
├── src/
│   ├── include/
│   │   ├── SylphDatabase.hpp               — RAII wrapper around `SylphDatabase*`
│   │   ├── SylphSketch.hpp                 — RAII wrapper around `SylphSketch*`
│   │   └── sylph_profile.hpp               — table function declaration
│   ├── SylphDatabase.cpp                   — load/cache singleton, error mapping
│   ├── SylphSketch.cpp                     — Arrow-input sketch builder
│   └── sylph_profile.cpp                   — table function (per-sample wired)
├── test/
│   ├── cpp/test_SylphDatabase.cpp          — load/free, version, num_genomes, error paths
│   ├── cpp/test_SylphSketch.cpp            — sketch from in-memory bytes vs. CLI sketch
│   └── sql/
│       ├── sylph_profile.test              — golden-data SQL tests
│       ├── sylph_profile_realworld.test    — large data, gated by env var
│       └── sylph_profile_sample_id.test    — per-sample partitioning
├── data/sylph/
│   ├── tiny_reads_R1.fq.gz                 — small paired reads, 2 known genomes
│   ├── tiny_reads_R2.fq.gz
│   ├── tiny.syldb                          — built from 5 small ref genomes
│   ├── tiny_refs/G*.fa                     — the 5 small ref genomes (for regen)
│   ├── expected_profile.tsv                — golden output from upstream `sylph profile`
│   ├── tiny_oracle.submodule.sha           — sylph fork SHA that produced expected_profile.tsv
│   └── README.md                           — regen instructions
└── tools/sylph_oracle.sh                   — script to regenerate expected_profile.tsv
                                              from the embedded fork's `sylph` binary
```

### Layered design

1. **Pure-Rust core (`profile_api`)** — `fn profile(genomes: &[GenomeSketch], sample: &SequencesSketch, args: ProfileArgs) -> Vec<OwnedAniResult>`. No I/O, no globals, no `process::exit`. This is the testable seam.
2. **C FFI (`c_api`)** — opaque handles, thread-local error string, `catch_unwind` at every entry, Arrow C Data Interface for output. Calls into `profile_api`.
3. **C++ RAII wrappers (`SylphDatabase`, `SylphSketch`)** — translate FFI errors into `std::runtime_error`, manage handle lifetimes, hide the `sylph.h` header from the rest of duckdb-miint.
4. **Table function (`sylph_profile`)** — wired through `per_sample_table_function`, reads sample sequences from a user-supplied table/view, calls into the C++ wrappers, materializes the returned Arrow batch into a `DataChunk`.

### Sample-input contract

Required input columns on the sample table:

| Column | Type | Meaning |
|---|---|---|
| `read_id` | VARCHAR | Read identifier (used only for read-pairing dedup; not surfaced) |
| `sequence1` | VARCHAR | R1 sequence (uppercase ACGTN) |
| `sequence2` | VARCHAR | R2 sequence; NULL or empty allowed for SE samples |

Optional column (when partitioning):

| Column | Type | Meaning |
|---|---|---|
| `sample_id` | any | Per-sample partition key (only when `sample_id:=` parameter set) |

This matches the schema produced by `read_fastx(r1:=..., r2:=...)`. Naming `sequence1`/`sequence2` keeps consistency with `detect_chimera_uchime` and friends — see the SQL examples in `PLAN-uchime.md`.

### Output schema

Always 9 columns, in this order. Adding columns later is non-breaking; reordering is.

| Column | Type | Meaning |
|---|---|---|
| `genome_index` | UINTEGER | Index into the `.syldb` (0..num_genomes), useful for joins |
| `genome_name` | VARCHAR | Genome filename / FASTA basename from `GenomeSketch.file_name` |
| `contig_name` | VARCHAR | First contig name from `GenomeSketch.first_contig_name` |
| `sequence_abundance` | DOUBLE | Fraction of effective coverage; sums to ~1.0 across rows |
| `taxonomic_abundance` | DOUBLE | Fraction of reads winner-take-all-assigned; sums to ~1.0 |
| `adjusted_ani` | DOUBLE | Coverage-corrected containment ANI |
| `eff_cov` | DOUBLE | Effective coverage estimate |
| `naive_ani` | DOUBLE | Raw containment ANI (no correction) |
| `kmers_reassigned` | UBIGINT | k-mers winner-take-all assigned to this genome |

When per-sample partitioning is on, a `sample_id` column is prepended, matching the source column's type — same convention as `align_mafft`/`uchime_*`.

### Named parameters

| Parameter | Type | Default | Notes |
|---|---|---|---|
| `syldb_path` | VARCHAR | (required, positional 2) | Local path to a `.syldb` |
| `sample_id` | VARCHAR | auto-detect | Column name for per-sample partition; auto-detects `sample_id` if present in the source schema (matches `read_fastx`'s convention) |
| `min_ani` | DOUBLE | 0.95 | Drop genomes below this adjusted ANI |
| `min_number_kmers` | UINTEGER | 50 | Drop genomes with fewer matched k-mers |
| `min_count_correct` | DOUBLE | 3.0 | λ-correction floor |
| `estimate_unknown` | BOOLEAN | false | Renormalize to fraction-of-reads-explained |
| `dedup_paired_reads` | BOOLEAN | true | Toggle sylph's read-pair dedup |
| `dedup_fpr` | DOUBLE | 0.001 | Cuckoo-filter FPR; 0 = exact dedup |
| `min_read_length` | UINTEGER | 0 | 0 = use sylph default; reads shorter dropped |
| `threads` | UINTEGER | 0 | 0 = use DuckDB scheduler thread count |

The defaults match upstream sylph 0.9.0's `--reads` profile defaults exactly. Tested in the SQL availability tests.

---

## The Rust fork: what to add, file by file

### `Cargo.toml` changes

```toml
[package]
# ... existing ...

[lib]
crate-type = ["staticlib", "rlib"]

[features]
default = ["cli", "fastx"]
cli = ["dep:clap", "dep:simple_logger", "dep:env_logger"]
fastx = ["dep:needletail"]                    # path-based sketching needs needletail
arrow-ffi = ["dep:arrow"]                     # FFI output

[dependencies]
# ... existing, but with these changes:
clap = { version = "3", features = ["derive"], optional = true }
simple_logger = { version = "3", features = ["stderr"], optional = true }
needletail = { version = "0.5.0", optional = true }
arrow = { version = "57", default-features = false, features = ["ffi"], optional = true }

# REMOVE entirely:
# [target.'cfg(target_env = "musl")'.dependencies]
# tikv-jemallocator = "0"

[profile.release]
# REMOVE: panic = "abort"           — incompatible with FFI catch_unwind
lto = true
```

`fastx` and `cli` go behind features so the duckdb-miint static link doesn't drag in the CLI surface.

### `src/lib.rs`

Add the two new modules:

```rust
pub mod sketch;
pub mod constants;
pub mod types;
pub mod seeding;
#[cfg(feature = "cli")]
pub mod cmdline;
pub mod contain;
pub mod inference;
pub mod inspect;

pub mod profile_api;     // NEW
#[cfg(feature = "arrow-ffi")]
pub mod c_api;           // NEW

#[cfg(target_arch = "x86_64")]
pub mod avx2_seeding;
```

### `src/sketch.rs` — minimal additions

Add a paired-end byte-slice sketcher that mirrors `sketch_pair_sequences` (sketch.rs:773) but takes an iterator yielding `(&[u8], Option<&[u8]>)` instead of opening files. Existing functions stay byte-for-byte identical.

```rust
/// Sketch paired-end reads from in-memory byte slices.
///
/// `reads` yields `(r1, Some(r2))` for paired-end or `(r1, None)` for SE.
/// All other parameters match `sketch_pair_sequences`.
pub fn sketch_pair_slices<'a, I>(
    reads: I,
    sample_name: Option<String>,
    c: usize,
    k: usize,
    no_dedup: bool,
    dedup_fpr: f64,
) -> SequencesSketch
where
    I: IntoIterator<Item = (&'a [u8], Option<&'a [u8]>)>,
{
    // Body is structurally identical to lines 790-895 of sketch_pair_sequences,
    // but operates on slices and never touches the filesystem.
    // Helpers reused as-is: extract_markers, pair_kmer, pair_kmer_single,
    //                      dup_removal_lsh_full_exact, dup_removal_lsh_full
    todo!()
}
```

This is the only user-visible behavioural addition to `sketch.rs`. Estimated ~80 lines.

### `src/contain.rs` — surgical refactor

The compute kernel is currently entangled with I/O and printing. Split it.

Functions to keep where they are: `get_stats` (line 613), `winner_table` (422), `derep_if_reassign_threshold` (365), `estimate_true_cov` (389), `estimate_covered_bases` (403), `ani_from_lambda` (829), `bootstrap_interval` (861), `get_kmer_identity` (913). All already pure compute over already-loaded data. Make them `pub(crate)` if not already.

Function to split: `contain` (line 115). Extract its post-load compute into a new pure function `run_profile_compute` that lives in `profile_api` (below). The existing `contain` becomes a thin wrapper: load files → call `run_profile_compute` → format output to TSV. The CLI behaviour is unchanged.

`get_seq_sketch` and `get_genome_sketches` (lines 556 and 494) — leave as-is; they're file-loading helpers used only by the CLI path.

### `src/profile_api.rs` — NEW

```rust
//! Pure data-in/data-out profile entry point.
//!
//! Takes already-loaded reference genomes and a sample sketch, returns
//! owned profile results with no borrows back into the inputs. This is
//! the seam tested by the C FFI and by the duckdb-miint integration.

use crate::types::{GenomeSketch, SequencesSketch, AdjustStatus};

/// Owned (non-borrowing) variant of `AniResult`. Strings are copied out
/// of the input GenomeSketch references so callers can free the input
/// independently of the results.
#[derive(Debug, Clone)]
pub struct OwnedAniResult {
    pub genome_index: u32,
    pub genome_name: String,
    pub contig_name: String,
    pub sequence_abundance: f64,
    pub taxonomic_abundance: f64,
    pub adjusted_ani: f64,
    pub eff_cov: f64,
    pub naive_ani: f64,
    pub kmers_reassigned: u64,
    pub lambda: AdjustStatus,
}

#[derive(Debug, Clone)]
pub struct ProfileArgs {
    pub min_ani: f64,
    pub min_count_correct: f64,
    pub min_number_kmers: u32,
    pub estimate_unknown: bool,
    pub num_threads: usize,
}

impl Default for ProfileArgs {
    fn default() -> Self {
        // Mirror `sylph profile --reads` defaults. Verified against
        // sylph 0.9.0 `cmdline.rs`.
        Self { min_ani: 0.95, min_count_correct: 3.0,
               min_number_kmers: 50, estimate_unknown: false, num_threads: 0 }
    }
}

pub fn run_profile_compute(
    genomes: &[GenomeSketch],
    sample: &SequencesSketch,
    args: &ProfileArgs,
) -> Vec<OwnedAniResult> {
    // Adapter that calls the existing private compute path:
    //   1. get_stats over (genomes, sample) -> Vec<AniResult>
    //   2. derep_if_reassign_threshold
    //   3. winner_table -> taxonomic abundance
    //   4. copy borrowed strings into OwnedAniResult
    // Threading: configure rayon thread pool if num_threads != 0.
    todo!()
}
```

### `src/c_api.rs` — NEW (~600 LOC)

Implements `sylph.h`. Follow rype's `c_api.rs` for shape: opaque pointer types via `Box::into_raw`/`Box::from_raw`, thread-local last-error via `thread_local! { static LAST_ERROR: RefCell<Option<CString>> }`, `catch_unwind` wrapping every public function, Arrow C Data Interface output via `arrow::ffi::FFI_ArrowArray::new` / `FFI_ArrowSchema::try_from`.

Concrete Rust signatures (matching `sylph.h`):

```rust
#[no_mangle]
pub extern "C" fn sylph_database_load(path: *const c_char) -> *mut SylphDatabase;
#[no_mangle]
pub extern "C" fn sylph_database_free(db: *mut SylphDatabase);
#[no_mangle]
pub extern "C" fn sylph_database_num_genomes(db: *const SylphDatabase) -> usize;

#[no_mangle]
pub extern "C" fn sylph_sketch_paired_arrow(
    r1_array: *const ffi::FFI_ArrowArray,
    r2_array: *const ffi::FFI_ArrowArray,         // may be null
    r1_schema: *const ffi::FFI_ArrowSchema,
    r2_schema: *const ffi::FFI_ArrowSchema,        // may be null
    params: *const SylphSketchParams,
) -> *mut SylphSketch;

#[no_mangle]
pub extern "C" fn sylph_sketch_builder_create(params: *const SylphSketchParams) -> *mut SylphSketch;
#[no_mangle]
pub extern "C" fn sylph_sketch_builder_add_pair(
    builder: *mut SylphSketch,
    r1: *const u8, r1_len: usize,
    r2: *const u8, r2_len: usize,                  // r2_len 0 = SE
) -> i32;
#[no_mangle]
pub extern "C" fn sylph_sketch_builder_finalize(builder: *mut SylphSketch) -> i32;
#[no_mangle]
pub extern "C" fn sylph_sketch_free(sketch: *mut SylphSketch);

#[no_mangle]
pub extern "C" fn sylph_profile(
    db: *const SylphDatabase,
    sample: *const SylphSketch,
    params: *const SylphProfileParams,
    out_array: *mut ffi::FFI_ArrowArray,
    out_schema: *mut ffi::FFI_ArrowSchema,
) -> i32;

#[no_mangle]
pub extern "C" fn sylph_get_last_error() -> *const c_char;
#[no_mangle]
pub extern "C" fn sylph_version() -> *const c_char;
#[no_mangle]
pub extern "C" fn sylph_miint_fork_version() -> *const c_char;
```

Internal types:

```rust
pub struct SylphDatabase {
    genomes: Vec<GenomeSketch>,
    syldb_kind: SyldbKind,        // params from the .syldb header (k, c)
}

pub enum SylphSketchInner {
    Building { ctx: BuilderCtx },  // accumulator for sketch_builder_*
    Finalized(Box<SequencesSketch>),
}
pub struct SylphSketch { inner: SylphSketchInner, params: SylphSketchParams }
```

Arrow input handling: `FFI_ArrowArray::from_raw` + `FFI_ArrowSchema::from_raw` → reconstitute `LargeStringArray` (DuckDB exports VARCHAR as `LargeUtf8`). Validate the schema is `LargeUtf8` (or `LargeBinary`); reject anything else with a clear error.

Arrow output handling: build a 9-column `RecordBatch` once per `sylph_profile` call, then extract its arrays into the caller-provided `FFI_ArrowArray` / `FFI_ArrowSchema` slots. arrow-rs handles the FFI export.

### `src/main.rs` — gate behind `cli` feature

```rust
#[cfg(feature = "cli")]
fn main() { /* existing body */ }

#[cfg(not(feature = "cli"))]
fn main() { /* empty — staticlib build doesn't need a binary */ }
```

### `sylph.h` — final FFI contract

The header drafted in our prior turn becomes the source of truth, lightly edited to match the Rust signatures above (notably `r2_len=0` instead of `NULL` for SE in the streaming builder, since byte-pointer-with-zero-length is more idiomatic in Rust FFI than null pointers).

---

## CMake integration (mirrors rype)

Add to `CMakeLists.txt` (around line 619, where rype starts):

```cmake
# ----- sylph (Rust, Arrow FFI) -----
option(MIINT_ENABLE_SYLPH "Build with sylph profiling support" ON)
if(EMSCRIPTEN OR WIN32)
    set(MIINT_ENABLE_SYLPH OFF)
    message(STATUS "sylph disabled on this platform")
endif()

if(MIINT_ENABLE_SYLPH)
    set(SYLPH_LIB_NAME "libsylph.a")
    set(SYLPH_TARGET_DIR ${CMAKE_CURRENT_BINARY_DIR}/sylph_target)
    set(SYLPH_CARGO_TARGET_FLAG "")
    set(SYLPH_LIB_SUBDIR "release")
    set(SYLPH_CONFIGURE_CMD "")

    # Cross-compilation matrix copied from rype block (apple cross only — no
    # WASM/MinGW because sylph is gated off above)
    if(APPLE AND OSX_BUILD_ARCH STREQUAL "x86_64")
        set(SYLPH_CARGO_TARGET_FLAG "--target" "x86_64-apple-darwin")
        set(SYLPH_LIB_SUBDIR "x86_64-apple-darwin/release")
        set(SYLPH_CONFIGURE_CMD ${RUSTUP} target add x86_64-apple-darwin)
    elseif(APPLE AND OSX_BUILD_ARCH STREQUAL "arm64")
        set(SYLPH_CARGO_TARGET_FLAG "--target" "aarch64-apple-darwin")
        set(SYLPH_LIB_SUBDIR "aarch64-apple-darwin/release")
        set(SYLPH_CONFIGURE_CMD ${RUSTUP} target add aarch64-apple-darwin)
    endif()

    set(SYLPH_LIB_PATH ${SYLPH_TARGET_DIR}/${SYLPH_LIB_SUBDIR}/${SYLPH_LIB_NAME})

    ExternalProject_Add(
        sylph_build
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/sylph
        CONFIGURE_COMMAND "${SYLPH_CONFIGURE_CMD}"
        BUILD_COMMAND ${CARGO} rustc --release --no-default-features
            --features arrow-ffi --features fastx
            --lib --crate-type=staticlib
            --target-dir ${SYLPH_TARGET_DIR}
            --manifest-path ${CMAKE_CURRENT_SOURCE_DIR}/ext/sylph/Cargo.toml
            ${SYLPH_CARGO_TARGET_FLAG}
        INSTALL_COMMAND ""
        BUILD_ALWAYS TRUE
        BUILD_BYPRODUCTS ${SYLPH_LIB_PATH}
    )

    add_library(sylph STATIC IMPORTED GLOBAL)
    set_target_properties(sylph PROPERTIES IMPORTED_LOCATION ${SYLPH_LIB_PATH})
    add_dependencies(sylph sylph_build)

    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_SYLPH
        SYLPH_GIT_VERSION="${SYLPH_GIT_VERSION}")
    target_link_libraries(${EXTENSION_NAME} sylph)
endif()
```

Register `SYLPH_GIT_VERSION` next to `RYPE_GIT_VERSION` (around line 137-157).

Source globbing: gate `sylph_profile.cpp`, `SylphDatabase.cpp`, `SylphSketch.cpp` with `if(MIINT_ENABLE_SYLPH)`. Mirror the vsearch pattern at line 802.

WASM `LINKED_LIBS` handling: **not needed** — sylph is gated off on WASM.

`run_tests.sh` guard: add a `SYLPH_AVAILABLE` env var detected by querying `duckdb_functions()` for `sylph_profile`, mirroring the `VSEARCH_AVAILABLE` block at line 100-101.

`miint_versions()` entry: add a Sylph row to the `miint_versions` table function so users can verify the embedded fork version (mirrors how HDF5/sortmerna show up).

---

## C++ side

### `src/include/SylphDatabase.hpp`

```cpp
#pragma once
#include <memory>
#include <string>

struct SylphDatabase;  // opaque from sylph.h

namespace miint {

class SylphDatabaseHandle {
public:
    explicit SylphDatabaseHandle(const std::string &path);  // throws on load failure
    ~SylphDatabaseHandle();
    SylphDatabaseHandle(const SylphDatabaseHandle &) = delete;
    SylphDatabaseHandle &operator=(const SylphDatabaseHandle &) = delete;

    const ::SylphDatabase *raw() const { return db_; }
    size_t num_genomes() const;

private:
    ::SylphDatabase *db_;
};

}  // namespace miint
```

### `src/include/SylphSketch.hpp`

Mirror the database wrapper. Add a static factory method that takes Arrow C Data Interface pointers (from a DuckDB-side `ArrowConverter`). Avoid making the C++ side depend on the `arrow-rs` crate — use the C ABI structs from `nanoarrow` (already available in the duckdb-miint dep tree via DuckDB) or hand-rolled `ArrowArray`/`ArrowSchema` C structs.

### `src/sylph_profile.cpp`

Wired via `per_sample_table_function.hpp` (see `site/src/content/docs/internals/per-sample-pattern.md` lines 1-50). Per-sample mode is the default expected workflow when users have multiple samples in one table.

Bind: validate input table has `read_id`, `sequence1`, `sequence2` columns; resolve `syldb_path`; load `.syldb` via a session-scoped cache (key: realpath; value: `shared_ptr<SylphDatabaseHandle>`).

Init: `MaxThreads()` returns `min(db_threads, num_samples)`. Per-sample parallelism takes precedence (DuckDB scheduler runs N samples concurrently). Within each sample, set rayon to `max(1, db_threads / num_samples_in_flight)` to fill leftover cores when num_samples < db_threads. Single-sample case becomes `rayon = db_threads, duckdb = 1` automatically; many-samples case becomes `rayon = 1, duckdb = N`. This matches the user's "encourage batching by sample" direction.

Implementation: each sample's exec call configures rayon via a scoped `rayon::ThreadPoolBuilder::new().num_threads(...).build_scoped(...)` so the thread count is per-call, not process-global. Test the math: 8 db threads + 1 sample → rayon=8; 8 db threads + 16 samples → rayon=1, duckdb=8 in parallel.

Execute: per sample, stream the sample's rows through `sylph_sketch_builder_add_pair`, finalize, call `sylph_profile`, materialize the returned Arrow batch into the output `DataChunk`.

Materialization: for the 9-column output, use the existing `arrow-zero-copy.md` recipe — DuckDB Vectors point into Arrow buffers with reference-counted release callbacks. The Arrow batch is small (one row per surviving genome), so cost is negligible regardless.

### Session-scoped database cache

Single-process map `unordered_map<string, weak_ptr<SylphDatabaseHandle>>` guarded by mutex, keyed by canonicalized path. Load returns a strong pointer, drops back to weak when no table function holds it. mtime-invalidation: if the file's mtime changes between loads, reload (defensive — `.syldb` is supposed to be immutable but users sometimes regenerate).

---

## Test-first plan (Red → Green → Refactor)

Each task lists the *test* that goes red first, then *what makes it green*. Same shape as `PLAN-uchime.md`.

### Phase 0 — fixture preparation (no code, only data) — **DONE 2026-04-30**

**Goal**: create the smallest possible reproducible test corpus that exercises the full compute path.

Decision deviating from initial plan: instead of synthesizing reads from 5 RefSeq genomes (which would need wgsim/ART, neither installed locally), we adopt sylph's own bundled `test_files/` fixtures verbatim. They're already redistributable under MIT/Apache-2.0 and small (~5 MB compressed).

Committed at `data/sylph/`:

- `tiny_refs/e.coli-{EC590,K12,o157}.fasta.gz` — 3 reference genomes (~1.5 MB each compressed)
- `tiny_reads_R{1,2}.fq.gz` — sylph's bundled K12 paired reads (gzipped, ~96 KB each)
- `tiny.syldb` — sketched database over the 3 refs
- `expected_profile.tsv` — golden output (md5 deterministic)
- `tiny_oracle.submodule.sha` — upstream sylph commit SHA pin (`cf6ee06...`)
- `README.md` — provenance + regen instructions
- `tools/sylph_oracle.sh` — regen script (uses embedded fork's CLI when available, falls back to sibling `../sylph/` clone in Phase 0)

Expected golden: only K12 passes the default `min_ani=0.95` threshold; sequence_abundance and taxonomic_abundance both 100.0000. EC590 and o157 are E. coli sister strains that don't pass the cutoff — the fixture is intentionally exercising the "containment ANI cutoff drops siblings" path.

**Caveat for later phases**: this single-genome golden doesn't test multi-genome abundance estimation. When Phase 5 lands, we'll add a separate gated corpus with mixed reads from multiple genomes for the multi-row golden. For Phases 1–4, single-genome is sufficient to exercise schema, FFI, and SQL correctness.

### Phase 1 — fork plumbing (Rust, no FFI yet)

**RED test 1.1** — once clap/needletail/simple_logger are made `optional` in Cargo.toml without gating their imports, `cargo check --no-default-features --lib` fails with "unresolved import `clap`" / "unresolved import `needletail`".

**GREEN — DONE 2026-04-30 (sylph commit `3467148`):** Three features added: `cli` (clap+simple_logger), `fastx` (needletail), and an implicit "pure compute" mode when neither is on. `cli` requires `fastx` (the orchestrators use needletail). `[lib] crate-type = ["staticlib", "rlib"]` added. Dropped `tikv-jemallocator` and `panic="abort"`. Gated `cmdline` and `inspect` modules in lib.rs; gated each `SketchArgs`/`ContainArgs`-taking function and each needletail-using function. Stub `fn main()` for the no-cli build.

Verification:
- `cargo build --release` (default features) → green, 2m07s, produces working `target/release/sylph`.
- `cargo check --no-default-features --lib` → green.
- `cargo rustc --release --no-default-features --lib --crate-type=staticlib` → green, 6.9 MB `libsylph.a`.
- Oracle script reproduces md5 `d8344148...` for `expected_profile.tsv` from the new feature-gated build — bit-identical to upstream behaviour.

Known cosmetic warnings (not blockers): ~30 `unused import` / `dead code` warnings in the no-default-features build are dormant compute helpers that Phases 1.2 (sketch_pair_slices) and 1.3 (profile_api) will exercise.

`arrow-ffi` feature is **not** added in 1.1 — it lands with the `arrow` dep in Phase 2.1 when `c_api` is introduced.

**RED test 1.2** — `sketch_pair_slices` produces the same `kmer_counts` HashMap as `sketch_pair_sequences` for a fixed pair of small FASTQ files vs. their byte contents.

**GREEN — DONE 2026-04-30 (sylph commit `7aceef1`):** Added `SketchPairBuilder` (streaming primitive) and `sketch_pair_slices` (one-shot iterator wrapper) at the bottom of `src/sketch.rs`. Always available (no feature gate). 3 tests in `src/sketch.rs::tests` all green: paired-end equivalence, single-end path, and streaming-builder vs. one-shot equivalence. The streaming builder is what the FFI's `sylph_sketch_builder_*` entry points will wrap directly.

**RED test 1.3** — `run_profile_compute(genomes, sample, default_args)` against in-memory fixtures returns rows matching `expected_profile.tsv` to high precision.

**GREEN — DONE 2026-04-30 (sylph commit `40a5233`):** Created `src/profile_api.rs` with `ProfileArgs`, `LambdaEstimator`, `OwnedAniResult`, and `run_profile_compute`. Refactored `get_stats` and `bootstrap_interval` to take `&ProfileArgs` instead of `&ContainArgs` and ungated them (the cli gate previously blocked their use from FFI builds). Exposed `winner_table`/`derep_if_reassign_threshold`/`estimate_true_cov`/`estimate_covered_bases`/`get_kmer_identity` as `pub(crate)`. `contain::contain` now builds a `ProfileArgs` from its `ContainArgs` via `profile_args_from_contain_args`. CLI behaviour unchanged: oracle TSV md5 still `d8344148...`.

`tests/profile_api_integration.rs` (2 tests, both green): in-process sketch of 3 E. coli refs + the bundled K12 paired reads, run profile, assert only K12 passes (sister strains drop below default 95% ANI cutoff), tax/seq abundance ≈ 100%, ANI fields within 0.5 pp of expected_profile.tsv. Query mode (pseudotax=false) keeps all 3 genomes.

### Phase 2 — C FFI (still no DuckDB)

**RED test 2.1** — `SylphDatabaseHandle("data/sylph/tiny.syldb")` constructs successfully, `num_genomes() == 3`, destructor frees without crash. Loading a nonexistent path throws with a useful error message.

**GREEN — DONE 2026-04-30 (sylph commit `47e02a5`, miint commit `170f2da`):**
- sylph fork: `src/c_api.rs` (~270 LOC) with thread-local last-error, panic-safe `catch_unwind` wrappers, opaque `SylphDatabase`, version helpers. `sylph.h` is the contract.
- duckdb-miint: CMake `MIINT_ENABLE_SYLPH` option + ExternalProject_Add building libsylph.a; `MIINT_HAS_SYLPH` define + sylph linked into ${EXTENSION_NAME} / ${LOADABLE_EXTENSION_NAME} / `tests`. `src/SylphDatabase.{hpp,cpp}` RAII wrapper. `test/cpp/test_SylphDatabase.cpp` Catch2 test with 4 cases.
- 6 sylph c_api unit tests + 1 integration test pass; full duckdb-miint build verification deferred until vcpkg deps finish rebuilding.

**RED test 2.2** — sylph_sketch_builder_create/add_pair/finalize/free FFI, plus tests for round-trip add_pair, NULL-handling, post-finalize behaviour.

**GREEN — DONE 2026-04-30 (sylph commit `b1d9fe0`):** Added `SylphSketch` opaque (Building/Finalized two-state), `SylphSketchParams`, and 4 FFI entries. Internal `sketch_borrow_finalized` helper for Phase 2.4. 5 c_api unit tests pass. Updated sylph.h. Staticlib now exports 10 sylph_* symbols.

**RED test 2.3** — round-trip: build a sketch via `sylph_sketch_builder_*` and compare against `sketch_pair_sequences` on the same data. Equal `kmer_counts`.

**GREEN — DONE 2026-04-30 (sylph commit `bbadaa7`):** `c_api::tests::ffi_builder_matches_path_sketcher` (gated by `fastx`) feeds the bundled K12 paired reads through both sketchers with `dedup_fpr=0` and asserts `kmer_counts` HashMap equality. Test passes — confirms the FFI is a faithful wrapper of the streaming builder.

**RED test 2.4** — `sylph_profile` against an in-memory syldb with the FFI sketch produces results matching `expected_profile.tsv`, returned via Arrow C Data Interface.

**GREEN — DONE 2026-04-30 (sylph commit `dd2b49f`):** Added `arrow = 57` (default-features off, ffi feature) behind a new `arrow-ffi` Cargo feature. `SylphProfileParams` C struct (layout-stable with reserved fields, 0-default-encoded). `sylph_profile(db, sample, params, out_array, out_schema)` populates caller-owned FFI_ArrowArray + FFI_ArrowSchema slots with the 9-column RecordBatch. `tests/c_api_profile.rs` integration test (gated `arrow-ffi+fastx`) builds an in-memory syldb, sketches K12 reads via the FFI builder, runs sylph_profile, decodes the FFI batch back into arrow types, and asserts: 1 row, 9 columns in the documented order, K12 in genome_name, abundances ≈ 100%, adjusted_ani within 0.5pp of 98.89. Test passes (cargo runtime 19m37s on first compile, then 0.57s for the test itself). The arrow-ffi staticlib is 8.7 MB and exports 11 sylph_* symbols (10 from Phase 2.2 + sylph_profile).

duckdb-miint CMakeLists.txt updated: `cargo rustc` invocation now uses `--no-default-features --features arrow-ffi` so the linked libsylph.a includes sylph_profile (needed by Phase 3's table function).

### Phase 3 — DuckDB table function

**RED test 3.1** — `SELECT * FROM duckdb_functions() WHERE function_name = 'sylph_profile'` returns 1 row.

**GREEN — IMPLEMENTED 2026-04-30 (miint commit `479f78f`):** registered `sylph_profile` via `RegisterDocumentedTableFunction`. Output schema fixed at 9 columns. SQL test asserts function appears in `duckdb_functions()` with `function_type = table` and the schema is visible via DESCRIBE.

**RED test 3.2/3.3** — golden K12-only single-row recovery from `data/sylph/tiny.syldb`, with `taxonomic_abundance ≈ 100%`, `adjusted_ani ≈ 98.89` (matching `expected_profile.tsv`).

**GREEN — IMPLEMENTED 2026-04-30 (miint commit `479f78f`):** Bind validates source schema via `ValidateSequenceTableSchema`. InitGlobal loads the syldb via `SylphDatabaseHandle` and creates a `QuerySequenceStream` + empty sketch. Execute streams sequence rows through `sylph_sketch_builder_add_pair`, finalizes, calls `sylph_profile`, decodes the Arrow C Data Interface output into a 9-column DataChunk via the manual decoder in `IngestProfileBatch`. SQL test asserts row count and golden values to 0.5 percentage-points.

**RED test 3.4** — error paths: nonexistent table, missing `sequence1` column, missing `.syldb`.

**GREEN — IMPLEMENTED 2026-04-30:** bind-time errors flow through `ValidateSequenceTableSchema` → `BinderException("does not exist" / "missing required column 'sequence1'")`. `.syldb` errors flow through `SylphDatabaseHandle` → `IOException("failed to open ...")`. SQL test asserts each substring.

**RED test 3.5** — view support.

**GREEN — IMPLEMENTED 2026-04-30:** uses `QuerySequenceStream` which transparently reads from views or tables. SQL test creates a VIEW over the base table and asserts identical golden output.

**RED test 3.6** — named parameter validation.

**GREEN — IMPLEMENTED 2026-04-30:** all 8 named parameters (`min_ani`, `min_number_kmers`, `min_count_correct`, `estimate_unknown`, `dedup_paired_reads`, `dedup_fpr`, `threads`, `sample_id`) wired through to `SylphProfileParams`. SQL test asserts `min_ani := 50` → range error, `min_ani := 0.999` → drops K12 (no rows).

**VERIFIED 2026-05-01 (miint commit `9eed9e7`):** all Phase 3 SQL tests pass end-to-end against the built extension.
- `./build/release/extension/miint/tests "[sylph]"` → 10 assertions in 4 cases (Phase 2.1 Catch2).
- `./build/release/test/unittest "test/sql/sylph_profile.test"` → 21 assertions, all green.
- `./build/release/test/unittest "test/sql/sylph_profile_sample_id.test"` → 14 assertions, all green (Phase 4.1 per-sample).
- duckdb shell: `SELECT * FROM sylph_profile('reads', 'data/sylph/tiny.syldb')` returns K12 at 100% taxonomic abundance, 98.89% adjusted_ani — exact match to upstream sylph CLI's `expected_profile.tsv`.

**Full regression suite (2026-05-01):**
- C++ Catch2 (`./build/release/extension/miint/tests`) → 714 cases / 6502 assertions, 0 failed (34 skipped — bowtie2 not installed).
- SQL (`bash run_tests.sh`) → 169 / 170 cases pass. The 1 failure is `test/sql/read_jplace.test:75` (a glob test parsing multiple jplace files), which OOMs on this host under concurrent test execution — pre-existing environmental flake, **unrelated to sylph integration**.

### Phase 4 — per-sample partitioning

**RED test 4.1** — partition a multi-sample table via `sample_id:='sample_id'`. Output prepends a `sample_id` column.

**GREEN — IMPLEMENTED 2026-04-30 (miint commits `479f78f`, `7e92225`):** `Data` extends with `has_sample_id` + `sample_info`; `GlobalState` extends `PerSampleGlobalState`; new `LocalState` owns a per-thread `Connection`. Bind calls `DiscoverSamples` when `sample_id` is non-empty. Execute dispatches: per-sample mode loops through `ClaimNextSample` → build per-sample TEMP VIEW → sketch + profile → emit. Empty-buffer-skip-to-next-sample fix in `7e92225`. `max_threads_hint=1` in v1 (sortmerna precedent).

**RED test 4.2** — concurrency benchmark. **DEFERRED → planned 2026-05-01 PM.**

v1 used `max_threads_hint=1` copied from the sortmerna precedent. That precedent does NOT apply to sylph: sortmerna has a process-wide `g_run_mutex` that serializes everything; sylph's profile path is mutex-free (read-only mmap'd DB, per-call local sketch, rayon worker pool spawned per invocation).

**Performance impact of the bad default**, measured on PRJNA798058 (164 samples × ~7.7M paired reads avg, GTDB-r232 0.95-derep syldb on a 16-core/96G `long`-partition node):

- **Serial baseline (current `max_threads_hint=1`)** — sbatch 20176, in flight 2026-05-01 ~10:25, ETA ~2.3 h. Per-sample ~50 s sketch+profile after one-time DB load amortization.
- **Cross-sample parallel target** — change `max_threads_hint=0` (use DuckDB scheduler thread count) + pass `pp.num_threads=1` to sylph_profile so the inner rayon pool is single-threaded. Predicted wall: **~10–20 min** for the same workload. Net ~10-15× speedup.

**Why cross-sample > rayon-inside for this workload:** `extract_markers` (AVX2 minimizer extraction) is memory-bandwidth-bound; rayon-inside-one-sample plateaus around 6–10× single-thread on 16 cores. Cross-sample parallelism scales near-linearly because samples are fully independent — no cross-sample memory contention, no sync overhead.

**Code change (single-line):**
```cpp
// src/sylph_profile.cpp, InitGlobal:
//   was: InitPerSampleGlobal(context, *gstate, ..., /*max_threads_hint=*/1);
//   now:
InitPerSampleGlobal(context, *gstate, data.sample_info.sample_values.size(), /*max_threads_hint=*/0);
// And in BuildProfileParams (per-sample call):
pp.num_threads = (data.threads == 0) ? 1 : data.threads;  // was: data.threads
```

**Test additions:**
- A new SQL test in `test/sql/sylph_profile_sample_id.test` that asserts a multi-sample run completes correctly with multi-thread scheduling (no race/clobber between threads on the shared DB).
- Optional perf check (gated by `MIINT_SYLPH_PERF`) that times the all-samples run and asserts wall < 30 min on a 16-core box.

**Plan:** apply the change after job 20176 finishes (so we have a clean serial baseline to compare against), rebuild, resubmit the same all-samples sbatch, and record both wall times in this plan doc. Bonus: 32-core nodes should drop wall to ~5-10 min — sylph's per-sample work is pure CPU, scales linearly with cores until memory bandwidth saturates.

**Updated 2026-05-01 PM (post-baseline):** baseline measurement landed: serial all-samples = **2 h 55 min** (job 20176), parallel-with-fixed-default = **28 min 20 s** (job 20182), 27.7 GB → 47.9 GB peak RSS, 188% → 1322% CPU. The fixed default uses `max_threads_hint=0` + `pp.num_threads=1` per call (16 samples × 1 inner = 16 cores active on a 16-core node — full saturation in that environment).

**Further refinement (2026-05-01 PM):** the prior fix saturates cores when `num_samples >= num_cores` but underutilizes when `num_samples < num_cores` (e.g. 16 samples on a 32-core node → 16 cores active). Auto-balance default in `BuildProfileParams` now computes inner rayon as `max(1, db_threads / outer_threads)` so the (outer × inner) product targets full core utilization in either regime. Mirrors sylph CLI's `step = max(t/3+1, min(num_samples, t))` formula in spirit. See the head-to-head benchmark below for measurement.

#### Phase 4.2 head-to-head benchmark scaffolding

`tools/sylph_benchmark/` contains a 4-cell {sylph CLI, miint-sylph} × {1, 16 samples} runtime comparison. Submitted via `sbatch tools/sylph_benchmark/run_all.sbatch` (32 CPUs, 96 GB, 4 h cap). Both tools are run with their respective defaults (sylph CLI: `-t 32`, miint: no thread overrides — auto-balance is the new default). The miint side uses a fully streaming pipeline (`COPY (sylph_profile JOIN tax) TO 'output.tsv'`) with no intermediate tables to genuinely surface DuckDB's pipelining advantage. `compare_results.py` validates row-wise numerical equivalence within 0.5pp.

Plan: see `/home/lpatel/.claude/plans/goofy-growing-neumann.md` for the full design + predictions; `tools/sylph_benchmark/results_template.md` for the results table to fill in once it lands.

#### Threading-model comparison: upstream sylph CLI vs miint integration

To run a fair head-to-head, here's what's actually different between the two parallelism strategies. The plan is to run both side-by-side on the same node + sample set and record results in this section.

**Upstream `sylph profile`:**
- Single rayon global thread pool, sized by `-t/--threads`.
- 2-level rayon: outer `chunk.into_par_iter()` over samples + inner `genome_index_vec.par_iter()` over genomes per sample. Work-stealing balances across both levels.
- "Chunk size" auto-computed: profile-mode `step = max(t/3+1, min(num_samples, t))`. For 164 samples + 16 threads, step = 16 → 16 samples in flight at once, each fully fanning out across genomes.
- `--sample-threads N` overrides the auto-step.
- Database loaded once per process (you give it a list of samples on the command line).

**miint `sylph_profile()`:**
- Two separate thread pools: DuckDB's scheduler at the outer (cross-sample) level + rayon's pool inside each per-sample invocation. Coordinated by params, not by work-stealing.
- Outer parallelism: `max_threads_hint` (0 = DuckDB threads, capped at num_samples).
- Inner parallelism: `threads` named param → `pp.num_threads` → sylph's rayon thread count per call.
- Database loaded once per **session** (not per process), so within a DuckDB connection, multiple subsequent `sylph_profile` calls reuse the loaded DB.
- View-driven input: source can be a parquet, a DuckDB table, a SQL view — no file-list construction needed.

**Expected throughput parity:** for the typical case `num_samples >> num_cores`, both strategies should hit ~the same wall-clock because each sample becomes a CPU-bound work unit and cross-sample parallelism dominates. The miint-specific advantages are non-throughput (composability with SQL, single-session DB reuse, no FASTQ-on-disk requirement); the upstream-specific advantages are 2-level rayon work-stealing for `num_samples < num_cores` workloads.

**Benchmark protocol (todo):**
1. Run upstream `sylph profile -t 16 -1 *_R1.fq.gz -2 *_R2.fq.gz syldb.syldb -o ref.tsv` against PRJNA798058 (or a 16-sample subset). Record wall + RSS.
2. Run our parallel-fixed miint version against the same data via the all-samples sbatch. Record wall + RSS.
3. Run miint with `threads := 16` named param to force 2-level (outer DuckDB × inner rayon=16) — measure oversubscription cost.
4. Smaller-N test: 4 samples + 16 cores. Upstream's 2-level should beat our `outer=4, inner=1`. Quantify the gap.
5. Compare result TSVs row-by-row to confirm bit-equivalence (small numerical drift is expected from rayon reduction-order differences; confirm < 0.5pp).

Results table (filled when measurements land):

| Config | Wall | Peak RSS | CPU% | Notes |
|---|---|---|---|---|
| miint serial baseline (max=1, rayon=default) | (job 20176 in flight) | | | current default |
| miint cross-sample (max=0, rayon=1) | TBD | | | proposed default |
| miint hybrid (max=4, rayon=4) | TBD | | | comparison |
| upstream sylph CLI (-t 16) | TBD | | | reference |
| upstream sylph CLI (-t 16, sample-threads 16) | TBD | | | reference |

**RED test 4.3** — auto-detection of `sample_id` column. **DEFERRED:** keep as a follow-up commit. The current behaviour (explicit `sample_id := 'col'` only) covers the immediate use case.

### Phase 5 — real-data regression

**RED test 5.1** (`test/sql/sylph_profile_realworld.test`, gated by `MIINT_SYLPH_REAL_DATA`) — a ~100k-read corpus profiled against the 5-genome syldb produces the same TSV as the embedded fork's CLI.

**GREEN**: have `run_tests.sh` build `expected_real.tsv` on demand using `tools/sylph_oracle.sh` against the embedded fork. Same pattern as sortmerna's `data/sortmerna/real_oracle.blast.gz` (see `embedded-tools.md` line 86-87).

### Phase 6 — docs and discoverability

- `RegisterDocumented*` registration with `function_type = TABLE`, `categories = ["taxonomic_profiling"]`, full `description` and `example`. Auto-generates the user-facing reference page. Already done in the Phase 3.1 stub registration — review and expand as features land.
- **Extensive docstrings + header documentation per existing conventions.** The user pioneered a docstring style on another duckdb-miint branch — peek at it (e.g. `git log --all --oneline | head` to find the branch, then read recent C++ headers/sources for the style) and bring sylph's C++ files (`SylphDatabase.{hpp,cpp}`, `sylph_profile.{hpp,cpp}`) up to that bar. Every public class, every method, every non-obvious inline block. Include: lifetime/ownership rules, thread-safety notes, error semantics, citation of the Shaw & Yu sylph paper at the top of `sylph_profile.cpp`.
- Update `site/src/content/docs/internals/embedded-tools.md` with a sylph section after rype, documenting: location, purpose, build, cross-compile gotchas, force-rebuild stamp file path, the conda libstdcxx-15 break + recovery recipe.
- Smoke-test entry in `localdocs/cli-smoke-tests.md`.

---

## Build / run commands (for me, next session)

```bash
# Move the user's clone into ext/ as a fork:
git -C ext/sylph init   # only if not already a repo with our changes
# (Detail: discuss with user whether ext/sylph is a submodule, vendored, or
#  a sibling worktree. rype is a submodule; vsearch is a fork-as-submodule.
#  For sylph we'll do the same as rype to keep the build pattern uniform.)

# Build everything:
bash build.sh

# Force a sylph rebuild only (cargo's ExternalProject caches aggressively;
# touching source files does NOT invalidate it — see embedded-tools.md):
touch build/release/extension/miint/sylph_build-prefix/src/sylph_build-stamp/sylph_build-configure
bash build.sh

# Run the SQL tests:
bash run_tests.sh

# Run just sylph SQL tests:
./build/release/test/unittest "test/sql/sylph_profile.test"
./build/release/test/unittest "test/sql/sylph_profile_sample_id.test"

# Run the C++ Catch2 tests:
./build/release/extension/miint/tests "[sylph]"

# Real-data regression (gated):
MIINT_SYLPH_REAL_DATA=1 bash run_tests.sh

# Format/check before commit:
conda run -n duckdb-143 make format-fix
make format-check
```

---

## Open questions / decisions deferred

1. ~~Submodule vs vendored fork~~. **DECIDED 2026-04-30:** submodule, mirroring rype. Will be set up at `ext/sylph` pointing at `the-miint/sylph` fork, branch `v0.9.0-miint`.
2. ~~sylph version~~. **DECIDED 2026-04-30:** stay on 0.9.0 (latest upstream). User's local clone is the source of truth.
3. **Long-read mode.** Not in v1. Sylph's `--reads` vs `--long` distinction matters; I'll add a `read_type` parameter when adding long-read support, default `'short'`.
4. **`.sylsp` files** (sample sketches saved to disk). Not exposed in v1; sketching is always in-memory. Could add `sylph_save_sketch(sketch, path)` later.
5. **Multi-syldb queries.** sylph CLI accepts multiple syldb files. v1 supports one. Trivial to extend later by accepting an array of paths.
6. **Numerical tolerance for SQL tests.** 4 decimal places matches uchime's 2-decimal tolerance with margin. If sylph's float math diverges between platforms (it shouldn't — pure deterministic compute on the same inputs), tighten or loosen here.
7. **Allocator fight with DuckDB jemalloc.** rype's `mm_alloc_stubs.c` story applies. sylph uses Rust's default `System` allocator (after we strip `tikv-jemallocator`); cross-language alloc/free is fine because we never free Rust-allocated memory from C++ except through `sylph_*_free` entry points. Verify by running ASan on the test suite once the integration lands.
8. **Sequence VARCHAR encoding.** DuckDB VARCHAR exports as Arrow `LargeUtf8`. Sequences are ASCII-only and validated, so UTF-8 vs. raw bytes is a non-issue, but the FFI must accept either `LargeUtf8` or `LargeBinary` to keep the API forgiving.

---

## Phasing as commits

Following uchime's commit cadence (one TDD cycle per commit, more or less):

1. `sylph: vendor fork; cargo build flags; lib.rs feature gates` (Phase 1.1, no behaviour change)
2. `sylph: add sketch_pair_slices with byte-slice input` (Phase 1.2)
3. `sylph: extract profile_api::run_profile_compute from contain.rs` (Phase 1.3)
4. `sylph: add c_api opaque database + sketch handles, FFI errors` (Phase 2.1)
5. `sylph: add streaming sketch builder FFI` (Phase 2.2)
6. `sylph: add sylph_profile FFI returning Arrow C Data Interface` (Phase 2.3, 2.4)
7. `sylph: register sylph_profile table function (stub)` (Phase 3.1)
8. `sylph: implement sylph_profile bind/init/exec for single sample` (Phase 3.2, 3.3)
9. `sylph: error paths and view support` (Phase 3.4, 3.5)
10. `sylph: named parameters` (Phase 3.6)
11. `sylph: per-sample partitioning` (Phase 4.1, 4.2)
12. `sylph: real-data regression test and oracle script` (Phase 5.1)
13. `sylph: documented function registration and embedded-tools doc` (Phase 6)

Each commit must pass: `bash run_tests.sh`, `./build/release/extension/miint/tests`, `make format-check`. No commit lands red.

---

## Quick checklist before starting

- [x] ext/sylph as submodule on `the-miint/sylph` fork (decided 2026-04-30)
- [x] sylph 0.9.0 (decided 2026-04-30 — latest upstream)
- [x] 4-decimal tolerance accepted (default per plan)
- [x] sample_id auto-detected from `read_fastx` output schema (decided 2026-04-30)
- [x] threading: prefer per-sample batching, fill leftover cores with rayon (decided 2026-04-30)
- [x] Phase 0 fixtures generated and committed (2026-04-30 — bundled sylph test_files reused)
- [ ] First red test (1.1) written and runs red
- [ ] CI / pre-commit hooks pass on the empty branch before adding source

## Phase 0 follow-ups

- ~~Submodule wiring needs a github fork~~. **DONE 2026-04-30.** Forked to `lucaspatel/sylph`; submodule added at `ext/sylph` on branch `v0.9.0-miint` (commit `cf6ee06...`). Oracle script verified to reproduce `expected_profile.tsv` from the submodule build (same md5 as sibling-clone build). Rename to a `the-miint/sylph` later is a one-line `.gitmodules` URL change. Sibling clone at `../sylph/` is now redundant and can be deleted whenever convenient.
- **Conda toolchain conflict.** Building sylph with conda's PATH active fails: `rust-lld: undefined symbol: __libc_csu_{init,fini}` from conda's Scrt1.o. Phase 1's CMake `ExternalProject_Add` needs to run cargo in an environment that doesn't pick up conda's CRT. Either filter `CMAKE_C_COMPILER`/`CMAKE_LINKER` to system gcc, or invoke cargo with a stripped env (`env -i` worked in Phase 0). Document in `embedded-tools.md` next to rype's cross-compile section.
- **Conda libstdcxx-15 break (2026-04-30).** When conda-forge upgrades `libgcc` / `libstdcxx` to 15.x (typically arrives via a transitive dep of `duckdb-cli` etc.), vcpkg-built `.a` files (rocksdb, hdf5, expat, vsearch) start emitting `__isoc23_strtoul`, `__libc_single_threaded`, `getrandom`, `fcntl64` references that the link path's libc.so doesn't export, and the link fails. Fix: pin to libgcc 14.x or earlier:
  ```bash
  conda install -c conda-forge "libgcc<15" "libstdcxx<15" "libgomp<15" "_openmp_mutex"
  rm -rf build/release/vcpkg_installed build/release/CMakeCache.txt build/release/CMakeFiles
  bash build.sh
  ```
  After the downgrade, `bash build.sh` succeeds in ~2-3 min on a warm cache. Worth adding an explicit version constraint to the project's conda env file once one exists.

---

*Plan author: implementation agent (this conversation).*
*Last updated: 2026-04-30. Bump this date on substantive edits.*
