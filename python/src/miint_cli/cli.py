"""miint CLI — translates command-line arguments into SQL against the miint DuckDB extension."""

import argparse
import os
import sys

import duckdb

# ---------------------------------------------------------------------------
# Relation names used by the CLI.  Every CREATE and every SQL reference to
# an intermediate relation must go through these constants so a rename in
# one place cannot silently break another.
# ---------------------------------------------------------------------------

_V_REF_LENGTHS = "_ref_lengths"
_V_GENOME_LENGTHS = "_genome_lengths"
_V_GENOME_MAP = "_genome_map"
_V_ALIGNMENTS = "_alignments"
_V_PASSING_GENOMES = "_passing_genomes"
_V_FILTERED_ALIGNMENTS = "_filtered_alignments"
_V_QUERIES = "_queries"
_V_SUBJECTS = "_subjects"
_V_SEQUENCES = "_sequences"
_T_READ_TO_SHARD = "_read_to_shard"


# ---------------------------------------------------------------------------
# Connection
# ---------------------------------------------------------------------------


def get_connection(extension_path=None):
    """Create a DuckDB connection with the miint extension loaded."""
    if extension_path:
        con = duckdb.connect(config={"allow_unsigned_extensions": "true"})
        con.execute(f"LOAD '{_sq(extension_path)}'")
    else:
        con = duckdb.connect()
        con.execute("INSTALL miint FROM community")
        con.execute("LOAD miint")
    return con


# ---------------------------------------------------------------------------
# SQL escaping helpers
#
# DuckDB DDL statements (CREATE VIEW/TABLE) do not support prepared-statement
# parameters ($name).  DML statements (COPY, SELECT) do.  Throughout this
# file, prepared parameters are used wherever DuckDB allows them.  For DDL,
# we fall back to _sq() / _dq() escaped interpolation.  Every such site is
# marked with a "-- DDL: no prepared params" comment.
# ---------------------------------------------------------------------------

PARQUET_OPTS = "FORMAT PARQUET, parquet_version 'v2', compression 'zstd'"

# Extensions recognized as sequence input (FASTQ/FASTA variants).
_SEQUENCE_EXTENSIONS = frozenset([
    ".fastq", ".fastq.gz", ".fq", ".fq.gz",
    ".fasta", ".fasta.gz", ".fa", ".fa.gz", ".fna", ".fna.gz",
])


def _sq(s):
    """Escape a string for a SQL single-quoted literal: ' → ''."""
    return s.replace("'", "''")


def _dq(s):
    """Escape a string for a SQL double-quoted identifier: \" → \"\"."""
    return s.replace('"', '""')


def _check_file(path, parser, label="file"):
    """Abort via parser.error() if *path* does not exist."""
    if not os.path.exists(path):
        parser.error(f"{label} not found: {path}")


def detect_format(path):
    """Detect output format and compression from file extension.

    Returns (format_str, compressed) tuple.
    """
    lower = path.lower()

    # Check compound extensions first (.gz variants)
    gz_formats = {
        ".fastq.gz": ("FASTQ", True),
        ".fq.gz": ("FASTQ", True),
        ".fasta.gz": ("FASTA", True),
        ".fa.gz": ("FASTA", True),
        ".fna.gz": ("FASTA", True),
        ".sam.gz": ("SAM", True),
    }
    for ext, result in gz_formats.items():
        if lower.endswith(ext):
            return result

    simple_formats = {
        ".parquet": ("PARQUET", False),
        ".fastq": ("FASTQ", False),
        ".fq": ("FASTQ", False),
        ".fasta": ("FASTA", False),
        ".fa": ("FASTA", False),
        ".fna": ("FASTA", False),
        ".sam": ("SAM", False),
        ".bam": ("BAM", False),
        ".biom": ("BIOM", False),
    }
    for ext, result in simple_formats.items():
        if lower.endswith(ext):
            return result

    raise ValueError(f"Cannot detect format from extension: {path}")


def detect_input_type(path):
    """Classify an input path as 'parquet' or 'sequence'.

    Only .parquet and known sequence extensions are accepted.  Anything else
    raises ValueError so the user gets a clear message instead of a cryptic
    read_fastx error.
    """
    lower = path.lower()
    if lower.endswith(".parquet"):
        return "parquet"
    for ext in _SEQUENCE_EXTENSIONS:
        if lower.endswith(ext):
            return "sequence"
    raise ValueError(
        f"Unrecognized input extension for '{path}'. "
        f"Expected .parquet or a FASTQ/FASTA extension."
    )


def build_copy_options(fmt, compressed):
    """Build COPY option string for the given format."""
    if fmt == "PARQUET":
        return PARQUET_OPTS
    parts = [f"FORMAT {fmt.lower()}"]
    if compressed and fmt in ("FASTQ", "FASTA", "SAM"):
        parts.append("COMPRESSION 'gzip'")
    return ", ".join(parts)


def load_lengths(con, path):
    """Load a 2-column TSV into a non-temp view named _V_REF_LENGTHS.

    Supports both headerless and header-containing TSVs by using DuckDB's
    auto-detect, then normalizing column names by ordinal position.

    Non-temp view because COPY FORMAT SAM reads REFERENCE_LENGTHS via a
    separate connection that cannot see the temp catalog.
    """
    # DESCRIBE supports prepared params
    result = con.execute(
        "DESCRIBE SELECT * FROM read_csv($p, delim='\t')", {"p": path}
    )
    cols = [row[0] for row in result.fetchall()]
    if len(cols) < 2:
        raise ValueError(
            f"Lengths file must have at least 2 columns, got {len(cols)}"
        )
    # DDL: no prepared params.  _dq() escapes any " in column names.
    con.execute(f"""
        CREATE OR REPLACE VIEW {_V_REF_LENGTHS} AS
        SELECT "{_dq(cols[0])}"::VARCHAR AS name,
               "{_dq(cols[1])}"::BIGINT AS length
        FROM read_csv('{_sq(path)}', delim='\t')
    """)


def load_sequences(con, r1, r2=None, as_table=False):
    """Create a relation named _V_SEQUENCES (view or table).

    Parquet input is read with read_parquet; sequence files with read_fastx.
    Use as_table=True when the downstream function scans multiple times
    (e.g. align_minimap2_sharded needs random access per shard).

    Non-temp: the extension reads named tables via a separate connection
    that cannot see the temp catalog.
    """
    kind = "TABLE" if as_table else "VIEW"
    input_type = detect_input_type(r1)
    # DDL: no prepared params
    if input_type == "parquet":
        con.execute(
            f"CREATE OR REPLACE {kind} {_V_SEQUENCES} AS SELECT * FROM read_parquet('{_sq(r1)}')"
        )
    elif r2:
        con.execute(
            f"CREATE OR REPLACE {kind} {_V_SEQUENCES} AS SELECT * FROM read_fastx('{_sq(r1)}', sequence2='{_sq(r2)}')"
        )
    else:
        con.execute(
            f"CREATE OR REPLACE {kind} {_V_SEQUENCES} AS SELECT * FROM read_fastx('{_sq(r1)}')"
        )


# ---------------------------------------------------------------------------
# Commands: convert
# ---------------------------------------------------------------------------


def cmd_convert_sequence(con, args):
    """Convert FASTQ/FASTA to parquet."""
    params = {"r1": args.r1, "out": args.output}
    if args.r2:
        params["r2"] = args.r2
        sql = f"COPY (SELECT * FROM read_fastx($r1, sequence2=$r2)) TO $out ({PARQUET_OPTS})"
    else:
        sql = f"COPY (SELECT * FROM read_fastx($r1)) TO $out ({PARQUET_OPTS})"
    con.execute(sql, params)


def cmd_convert_alignment(con, args):
    """Convert SAM/BAM to parquet."""
    con.execute(
        f"COPY (SELECT * FROM read_alignments($input)) TO $out ({PARQUET_OPTS})",
        {"input": args.input, "out": args.output},
    )


def cmd_convert_biom(con, args):
    """Convert BIOM to parquet."""
    con.execute(
        f"COPY (SELECT * FROM read_biom($input)) TO $out ({PARQUET_OPTS})",
        {"input": args.input, "out": args.output},
    )


# Columns needed by each COPY format when converting FROM parquet.
_FORMAT_COLUMNS = {
    "FASTQ": "read_id, comment, sequence1, qual1",
    "FASTA": "read_id, comment, sequence1",
    "SAM": "read_id, flags, reference, position, mapq, cigar, "
           "mate_reference, mate_position, template_length",
    "BAM": "read_id, flags, reference, position, mapq, cigar, "
           "mate_reference, mate_position, template_length",
    "BIOM": "sample_id, feature_id, value",
}


def _select_for_format(fmt):
    """Return the SELECT column list appropriate for a COPY format.

    FASTQ/FASTA need specific columns to avoid the INTERLEAVE error when
    the parquet schema includes sequence2/qual2 from read_fastx.
    SAM/BAM/BIOM get their required columns to catch schema mismatches early.
    PARQUET gets '*'.
    """
    return _FORMAT_COLUMNS.get(fmt, "*")


def cmd_convert_parquet(con, args):
    """Convert parquet to a bioinformatics format (detected from -o extension)."""
    fmt, compressed = detect_format(args.output)

    if fmt in ("SAM", "BAM") and not args.lengths:
        args._parser.error(f"--lengths is required for {fmt} output")

    if args.lengths:
        load_lengths(con, args.lengths)

    options = build_copy_options(fmt, compressed)
    if fmt in ("SAM", "BAM"):
        options += f", REFERENCE_LENGTHS '{_V_REF_LENGTHS}'"

    cols = _select_for_format(fmt)
    con.execute(
        f"COPY (SELECT {cols} FROM read_parquet($input)) TO $out ({options})",
        {"input": args.input, "out": args.output},
    )


# ---------------------------------------------------------------------------
# Commands: transform
# ---------------------------------------------------------------------------


def _setup_genome_coverage_views(con, lengths_path):
    """Create the views that genome_coverage() macro needs as bare relation args.

    Creates: _V_REF_LENGTHS, _V_GENOME_LENGTHS, _V_GENOME_MAP
    """
    load_lengths(con, lengths_path)
    con.execute(f"""
        CREATE OR REPLACE VIEW {_V_GENOME_LENGTHS} AS
        SELECT name AS genome_id, length AS total_length FROM {_V_REF_LENGTHS}
    """)
    con.execute(f"""
        CREATE OR REPLACE VIEW {_V_GENOME_MAP} AS
        SELECT name AS contig_id, name AS genome_id FROM {_V_REF_LENGTHS}
    """)


def cmd_transform_genome_coverage(con, args):
    """Compute genome coverage from alignment parquet."""
    _setup_genome_coverage_views(con, args.lengths)
    # DDL: no prepared params
    con.execute(
        f"CREATE OR REPLACE VIEW {_V_ALIGNMENTS} AS SELECT * FROM read_parquet('{_sq(args.input)}')"
    )

    con.execute(
        f"""
        COPY (
            SELECT * FROM genome_coverage({_V_ALIGNMENTS}, {_V_GENOME_LENGTHS}, {_V_GENOME_MAP})
        ) TO $out ({PARQUET_OPTS})
    """,
        {"out": args.output},
    )


def cmd_transform_woltka_ogu(con, args):
    """Run WoLTKA OGU feature table generation with optional filters."""
    has_filters = (
        args.alignment_seq_identity is not None
        or args.alignment_query_coverage is not None
    )
    has_genome_cov = args.genome_coverage is not None

    if has_genome_cov and not args.lengths:
        args._parser.error("--lengths is required when --genome-coverage is specified")

    # Build alignments view (with optional filters) — DDL: no prepared params
    where_parts = []
    if args.alignment_seq_identity is not None:
        where_parts.append(
            f"alignment_seq_identity(cigar, tag_nm, tag_md, 'blast') >= {float(args.alignment_seq_identity)}"
        )
    if args.alignment_query_coverage is not None:
        where_parts.append(
            f"alignment_query_coverage(cigar) >= {float(args.alignment_query_coverage)}"
        )
    where_sql = (" WHERE " + " AND ".join(where_parts)) if where_parts else ""

    con.execute(f"""
        CREATE OR REPLACE VIEW {_V_ALIGNMENTS} AS
        SELECT * FROM read_parquet('{_sq(args.input)}'){where_sql}
    """)

    if has_genome_cov:
        _setup_genome_coverage_views(con, args.lengths)
        con.execute(f"""
            CREATE OR REPLACE VIEW {_V_PASSING_GENOMES} AS
            SELECT genome_id FROM genome_coverage(
                {_V_ALIGNMENTS}, {_V_GENOME_LENGTHS}, {_V_GENOME_MAP}
            ) WHERE proportion_covered >= {float(args.genome_coverage)}
        """)
        con.execute(f"""
            CREATE OR REPLACE VIEW {_V_FILTERED_ALIGNMENTS} AS
            SELECT a.* FROM {_V_ALIGNMENTS} a
            SEMI JOIN {_V_PASSING_GENOMES} g ON a.reference = g.genome_id
        """)
        con.execute(
            f"""
            COPY (
                SELECT * FROM woltka_ogu({_V_FILTERED_ALIGNMENTS}, read_id)
            ) TO $out ({PARQUET_OPTS})
        """,
            {"out": args.output},
        )
    else:
        con.execute(
            f"""
            COPY (
                SELECT * FROM woltka_ogu({_V_ALIGNMENTS}, read_id)
            ) TO $out ({PARQUET_OPTS})
        """,
            {"out": args.output},
        )


# ---------------------------------------------------------------------------
# Commands: align
# ---------------------------------------------------------------------------


def cmd_align_minimap2(con, args):
    """Align sequences with minimap2.

    Uses VIEWs (not TABLEs) since non-sharded minimap2 does a single pass.
    Non-temp because the extension reads tables via a separate connection.
    """
    # DDL: no prepared params
    input_type = detect_input_type(args.r1)
    if args.r2:
        con.execute(
            f"CREATE OR REPLACE VIEW {_V_QUERIES} AS SELECT * FROM read_fastx('{_sq(args.r1)}', sequence2='{_sq(args.r2)}')"
        )
    elif input_type == "parquet":
        con.execute(
            f"CREATE OR REPLACE VIEW {_V_QUERIES} AS SELECT * FROM read_parquet('{_sq(args.r1)}')"
        )
    else:
        con.execute(
            f"CREATE OR REPLACE VIEW {_V_QUERIES} AS SELECT * FROM read_fastx('{_sq(args.r1)}')"
        )

    # .mmi → index_path, otherwise → subject_table
    if args.database.lower().endswith(".mmi"):
        con.execute(
            f"""
            COPY (
                SELECT * FROM align_minimap2('{_V_QUERIES}', index_path=$db)
            ) TO $out ({PARQUET_OPTS})
        """,
            {"db": args.database, "out": args.output},
        )
    else:
        # DDL: no prepared params
        con.execute(
            f"CREATE OR REPLACE VIEW {_V_SUBJECTS} AS SELECT read_id, sequence1 FROM read_fastx('{_sq(args.database)}')"
        )
        con.execute(
            f"""
            COPY (
                SELECT * FROM align_minimap2('{_V_QUERIES}', subject_table='{_V_SUBJECTS}')
            ) TO $out ({PARQUET_OPTS})
        """,
            {"out": args.output},
        )


def cmd_align_minimap2_sharded(con, args):
    """Sharded minimap2 alignment with RYpe classification."""
    # Load sequences as TABLE (scanned multiple times)
    load_sequences(con, args.r1, args.r2, as_table=True)

    # Classify with RYpe to get shard assignments — DDL: no prepared params
    con.execute(f"""
        CREATE OR REPLACE TABLE {_T_READ_TO_SHARD} AS
        SELECT read_id, bucket_name AS shard_name
        FROM rype_classify('{_sq(args.rype_index)}', '{_V_SEQUENCES}', id_column='read_id')
    """)

    # Sharded alignment — COPY supports prepared params
    con.execute(
        f"""
        COPY (
            SELECT * FROM align_minimap2_sharded(
                '{_V_SEQUENCES}',
                shard_directory=$shard_dir,
                read_to_shard='{_T_READ_TO_SHARD}'
            )
        ) TO $out ({PARQUET_OPTS})
    """,
        {"shard_dir": args.shard_dir, "out": args.output},
    )


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def build_parser():
    parser = argparse.ArgumentParser(
        prog="miint",
        description="CLI for the miint DuckDB extension — bioinformatics SQL toolkit",
    )
    parser.add_argument(
        "--extension-path",
        help="Path to miint.duckdb_extension (default: install from community)",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- convert ---
    convert = subparsers.add_parser("convert", help="Convert between file formats")
    convert_sub = convert.add_subparsers(dest="convert_command", required=True)

    # convert sequence
    p = convert_sub.add_parser("sequence", help="Convert FASTQ/FASTA to parquet")
    p.add_argument("-1", dest="r1", required=True, help="Forward reads file")
    p.add_argument("-2", dest="r2", help="Reverse reads file (paired-end)")
    p.add_argument("-o", "--output", required=True, help="Output parquet file")
    p.set_defaults(func=cmd_convert_sequence)

    # convert alignment
    p = convert_sub.add_parser("alignment", help="Convert SAM/BAM to parquet")
    p.add_argument("-i", "--input", required=True, help="Input SAM/BAM file")
    p.add_argument("-o", "--output", required=True, help="Output parquet file")
    p.set_defaults(func=cmd_convert_alignment)

    # convert biom
    p = convert_sub.add_parser("biom", help="Convert BIOM to parquet")
    p.add_argument("-i", "--input", required=True, help="Input BIOM file")
    p.add_argument("-o", "--output", required=True, help="Output parquet file")
    p.set_defaults(func=cmd_convert_biom)

    # convert parquet
    p = convert_sub.add_parser(
        "parquet", help="Convert parquet to bioinformatics format"
    )
    p.add_argument("-i", "--input", required=True, help="Input parquet file")
    p.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output file (format detected from extension)",
    )
    p.add_argument(
        "--lengths", help="Reference lengths TSV (required for SAM/BAM output)"
    )
    p.set_defaults(func=cmd_convert_parquet, _parser=p)

    # --- transform ---
    transform = subparsers.add_parser("transform", help="Transform data")
    transform_sub = transform.add_subparsers(
        dest="transform_command", required=True
    )

    # transform genome-coverage
    p = transform_sub.add_parser("genome-coverage", help="Compute genome coverage")
    p.add_argument("-i", "--input", required=True, help="Input alignment parquet")
    p.add_argument("--lengths", required=True, help="Genome lengths TSV")
    p.add_argument("-o", "--output", required=True, help="Output parquet file")
    p.set_defaults(func=cmd_transform_genome_coverage)

    # transform woltka-ogu
    p = transform_sub.add_parser("woltka-ogu", help="WoLTKA OGU feature table")
    p.add_argument("-i", "--input", required=True, help="Input alignment parquet")
    p.add_argument("-o", "--output", required=True, help="Output parquet file")
    p.add_argument(
        "--alignment-seq-identity",
        type=float,
        help="Min sequence identity threshold (0-1)",
    )
    p.add_argument(
        "--alignment-query-coverage",
        type=float,
        help="Min query coverage threshold (0-1)",
    )
    p.add_argument(
        "--lengths", help="Genome lengths TSV (required with --genome-coverage)"
    )
    p.add_argument(
        "--genome-coverage",
        type=float,
        help="Min genome coverage threshold (0-1)",
    )
    p.set_defaults(func=cmd_transform_woltka_ogu, _parser=p)

    # --- align ---
    align = subparsers.add_parser("align", help="Sequence alignment")
    align_sub = align.add_subparsers(dest="align_command", required=True)

    # align minimap2
    p = align_sub.add_parser("minimap2", help="Align with minimap2")
    p.add_argument("-1", dest="r1", required=True, help="Query sequences")
    p.add_argument("-2", dest="r2", help="Reverse reads (paired-end)")
    p.add_argument(
        "-d", "--database", required=True, help="Subject FASTA or .mmi index"
    )
    p.add_argument("-o", "--output", required=True, help="Output parquet file")
    p.set_defaults(func=cmd_align_minimap2)

    # align minimap2-sharded
    p = align_sub.add_parser(
        "minimap2-sharded", help="Sharded minimap2 alignment"
    )
    p.add_argument("-1", dest="r1", required=True, help="Query sequences")
    p.add_argument("-2", dest="r2", help="Reverse reads (paired-end)")
    p.add_argument(
        "--shard-dir", required=True, help="Directory with .mmi shard indexes"
    )
    p.add_argument(
        "--rype-index", required=True, help="RYpe index for shard assignment"
    )
    p.add_argument("-o", "--output", required=True, help="Output parquet file")
    p.set_defaults(func=cmd_align_minimap2_sharded)

    return parser


def _validate_inputs(args, parser):
    """Check that input files exist before running SQL."""
    for attr in ("r1", "r2", "input", "lengths", "database", "rype_index"):
        path = getattr(args, attr, None)
        if path is not None:
            _check_file(path, parser, label=f"--{attr.replace('_', '-')}")
    if hasattr(args, "shard_dir"):
        path = args.shard_dir
        if path is not None and not os.path.isdir(path):
            parser.error(f"--shard-dir not found: {path}")


def main():
    parser = build_parser()
    args = parser.parse_args()

    _validate_inputs(args, parser)

    con = get_connection(args.extension_path)
    try:
        args.func(con, args)
    except (ValueError, duckdb.Error) as exc:
        parser.exit(1, f"error: {exc}\n")
    finally:
        con.close()


if __name__ == "__main__":
    main()
