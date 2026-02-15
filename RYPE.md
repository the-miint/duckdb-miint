# RYpe Integration Notes

## rype.h Header vs Rust Implementation — `use_merge_join` removed (RESOLVED)

The `use_merge_join` parameter was removed from the Rust implementation of
`rype_classify_arrow` and `rype_classify_arrow_best_hit` but the rype.h header
was not updated. This caused a silent ABI mismatch: the `0` value for
`use_merge_join` was interpreted as a NULL `out_stream` pointer.

**Fixed** by updating the rype submodule (commit d129d47) which corrected rype.h.

Our `src/rype_classify.cpp` was updated to call with 5 params (no `use_merge_join`).
