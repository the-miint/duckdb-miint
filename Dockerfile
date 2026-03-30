ARG DUCKDB_VERSION=1.5.1

FROM duckdb/duckdb:${DUCKDB_VERSION} AS duckdb

# Install extension in a shell-capable image
FROM ubuntu:24.04 AS builder
COPY --from=duckdb /duckdb /duckdb
RUN /duckdb -c "INSTALL miint FROM community;"

# Final image: DuckDB with pre-installed extension
FROM duckdb
COPY --from=builder /root/.duckdb /root/.duckdb

WORKDIR /data

ENTRYPOINT ["/duckdb", "-cmd", "LOAD miint;"]
