// Umbrella staticlib for rype + sylph + st3. The crate body is intentionally empty:
// cargo treats the upstream crates as rlib dependencies, rustc emits one
// staticlib that links them through a single Rust std, and the `#[no_mangle]
// extern "C"` symbols in rype/sylph survive into the archive's export table.
// `extern crate ... as _` forces the linkage without exposing names.

extern crate rype as _;

#[cfg(feature = "with-sylph")]
extern crate sylph as _;

#[cfg(feature = "with-st3")]
extern crate st3 as _;
