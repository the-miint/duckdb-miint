# microtar (vendored)

A lightweight tar library written in ANSI C.

- **Upstream:** https://github.com/rxi/microtar
- **License:** MIT (see `LICENSE`)
- **Version:** 0.1.0

Vendored verbatim (`microtar.c` + `microtar.h`, unmodified) and compiled directly
into the miint extension. Used by `src/taxdump_archive.cpp` to extract members from
NCBI taxonomy `taxdump.tar.gz` archives (the `.gz` layer is handled separately by
zlib). Only the read side of the API is exercised, driven from an in-memory buffer
via microtar's stream callbacks.
