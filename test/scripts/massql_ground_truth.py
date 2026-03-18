#!/usr/bin/env python3
"""
MassQL Ground Truth Generator

Runs MassQL queries through the Python reference implementation and records
results to a JSON file. The companion SQL test (massql_ground_truth.test)
compares DuckDB's massql() output against these reference values.

Usage:
    conda run -n massql python3 test/scripts/massql_ground_truth.py

Output:
    test/data/massql_ground_truth.json
"""

import warnings

warnings.filterwarnings("ignore")

import json
import os
import sys

from massql import msql_engine

MZML_FILE = "data/mzml/basic_3spectra.mzML"
OUTPUT_FILE = "test/data/massql_ground_truth.json"

# Each entry: (test_id, massql_query, checks)
# checks is a dict of what to verify:
#   "row_count": expected number of rows
#   "scans": expected sorted list of unique scan numbers (1-based, Python convention)
#   "spectrum_indices": expected sorted list of unique spectrum_index (0-based, DuckDB convention)
#   "sum_intensity": expected sum of intensity column (for scansum/raw data)
#   "max_intensity": expected max intensity (for scanmaxint)
#   "precursor_mz": expected precursor_mz value (for scanmz)
#   "mz_values": expected sorted unique mz values
QUERIES = [
    # --- Group 1: Basic peak matching ---
    {
        "id": "basic_ms2prod_220_raw",
        "query": "QUERY MS2DATA WHERE MS2PROD=220",
        "description": "All peaks from scans containing peak near 220",
    },
    {
        "id": "basic_ms2prod_220_scaninfo",
        "query": "QUERY scaninfo(MS2DATA) WHERE MS2PROD=220",
        "description": "Scan info for scans with peak near 220",
    },
    {
        "id": "basic_ms2prod_220_scannum",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=220",
        "description": "Scan numbers for scans with peak near 220",
    },
    {
        "id": "basic_ms2prec_200",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PREC=200",
        "description": "Scans with precursor near 200",
    },
    {
        "id": "basic_ms2nl_80",
        "query": "QUERY scannum(MS2DATA) WHERE MS2NL=80",
        "description": "Scans with neutral loss near 80",
    },
    {
        "id": "basic_ms1mz_200",
        "query": "QUERY scannum(MS1DATA) WHERE MS1MZ=200",
        "description": "MS1 scans with peak near 200",
    },
    # --- Group 2: Metadata conditions ---
    {
        "id": "meta_rtmin",
        "query": "QUERY scannum(MS2DATA) WHERE RTMIN=1.7",
        "description": "MS2 scans with RT >= 1.7",
    },
    {
        "id": "meta_polarity_positive",
        "query": "QUERY scannum(MS2DATA) WHERE POLARITY=Positive",
        "description": "Positive polarity MS2 scans",
    },
    {
        "id": "meta_charge_2",
        "query": "QUERY scannum(MS2DATA) WHERE CHARGE=2",
        "description": "MS2 scans with charge 2",
    },
    # --- Group 3: Intensity qualifiers ---
    {
        "id": "qual_intensitypercent_gt50",
        "query": "QUERY scannum(MS1DATA) WHERE MS1MZ=200:INTENSITYPERCENT>50",
        "description": "MS1 peak at 200 with normalized intensity > 50%",
    },
    {
        "id": "qual_intensityvalue_gt3000",
        "query": "QUERY scannum(MS1DATA) WHERE MS1MZ=200:INTENSITYVALUE>3000",
        "description": "MS1 peak at 200 with raw intensity > 3000",
    },
    {
        "id": "qual_intensityvalue_gt6000",
        "query": "QUERY scannum(MS1DATA) WHERE MS1MZ=200:INTENSITYVALUE>6000",
        "description": "MS1 peak at 200 with raw intensity > 6000 (none match)",
    },
    # --- Group 4: EXCLUDED qualifier ---
    {
        "id": "qual_excluded",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=150:EXCLUDED",
        "description": "MS2 scans WITHOUT peak near 150",
    },
    # --- Group 5: AND multi-condition ---
    {
        "id": "and_prod150_prec200",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=150 AND MS2PREC=200",
        "description": "MS2 with peak near 150 AND precursor near 200",
    },
    {
        "id": "and_prod220_prec300",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=220 AND MS2PREC=300",
        "description": "MS2 with peak near 220 AND precursor near 300",
    },
    {
        "id": "and_prod150_prec300_nomatch",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=150 AND MS2PREC=300",
        "description": "MS2 with peak near 150 AND precursor near 300 (no match)",
    },
    # --- Group 6: OR + CARDINALITY ---
    {
        "id": "or_150_250",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=(150 OR 250)",
        "description": "MS2 scans with peak near 150 OR 250",
    },
    {
        "id": "card_min2",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=(150 OR 250):CARDINALITY=range(min=2,max=5)",
        "description": "MS2 with at least 2 of (150, 250) matching",
    },
    {
        "id": "card_max1",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=(150 OR 250):CARDINALITY=range(min=1,max=1)",
        "description": "MS2 with exactly 1 of (150, 250) matching",
    },
    # --- Group 7: FILTER clause ---
    {
        "id": "filter_same_peak",
        "query": "QUERY MS2DATA WHERE MS2PROD=220 FILTER MS2PROD=220",
        "description": "WHERE finds scan, FILTER keeps only peak near 220",
    },
    {
        "id": "filter_scansum",
        "query": "QUERY scansum(MS2DATA) WHERE MS2PROD=220 FILTER MS2PROD=220",
        "description": "scansum with FILTER narrowing to single peak",
    },
    {
        "id": "filter_cross_peak",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=220 FILTER MS2PROD=320",
        "description": "FILTER does not re-evaluate scan membership",
    },
    # --- Group 8: Aggregation functions ---
    {
        "id": "agg_scansum",
        "query": "QUERY scansum(MS2DATA) WHERE MS2PROD=220",
        "description": "Total intensity of scan with peak near 220",
    },
    {
        "id": "agg_scanmaxint",
        "query": "QUERY scanmaxint(MS2DATA) WHERE MS2PROD=220",
        "description": "Max intensity of scan with peak near 220",
    },
    {
        "id": "agg_scanmz",
        "query": "QUERY scanmz(MS2DATA) WHERE MS2PROD=220",
        "description": "Precursor mz of scan with peak near 220",
    },
    # --- Group 9: Tolerance boundary ---
    {
        "id": "tol_default_within",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=220.09",
        "description": "Default 0.1 Da tolerance: 220.09 matches peak at 220",
    },
    {
        "id": "tol_default_outside",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=220.11",
        "description": "Default 0.1 Da tolerance: 220.11 does NOT match peak at 220",
    },
    # --- Group 10: X-offset patterns ---
    {
        "id": "xoffset_pair",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=X AND MS2PROD=X+100",
        "description": "X-offset pair: scans with two peaks 100 Da apart",
    },
    {
        "id": "xoffset_triplet",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=X AND MS2PROD=X+100 AND MS2PROD=X+200",
        "description": "X-offset triplet: three peaks spaced by 100",
    },
    {
        "id": "xoffset_plus_rtmin",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=X AND MS2PROD=X+100 AND RTMIN=1.7",
        "description": "X-offset pair with RT filter",
    },
    {
        "id": "xoffset_nomatch",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=X AND MS2PROD=X+500",
        "description": "X-offset pair with no match",
    },
    # --- Group 11: Cross-level X patterns ---
    {
        "id": "crosslevel_ms1_ms2prec",
        "query": "QUERY scannum(MS2DATA) WHERE MS1MZ=X AND MS2PREC=X",
        "description": "MS2 scans whose precursor matches an MS1 peak",
    },
    {
        "id": "crosslevel_plus_charge",
        "query": "QUERY scannum(MS2DATA) WHERE MS1MZ=X AND MS2PREC=X AND CHARGE=2",
        "description": "Cross-level with charge filter",
    },
    # --- Group 12: Cycle 3 — INTENSITYTICPERCENT ---
    {
        "id": "qual_intensityticpercent_gt50",
        "query": "QUERY scannum(MS1DATA) WHERE MS1MZ=200:INTENSITYTICPERCENT>50",
        "description": "MS1 peak at 200 with TIC-normalized intensity > 50%",
    },
    # --- Group 13: Cycle 4 — MASSDEFECT ---
    {
        "id": "massdefect_220",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=220:MASSDEFECT=massdefect(min=0.0,max=0.5)",
        "description": "Peak near 220 with mass defect in [0.0, 0.5) — integer m/z qualifies",
    },
    # --- Group 14: Cycle 6 — aminoaciddelta ---
    {
        "id": "aminoaciddelta_nomatch",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=aminoaciddelta(AG)",
        "description": "Peak at amino acid delta mass of AG (~128.058) — no match expected",
    },
    # --- Group 15: Cycle 13 — OTHERSCAN ---
    # Note: Python MassQL's OTHERSCAN does not expand to neighboring scans (appears buggy).
    # Our implementation correctly expands. We record Python's result but expect our result to differ.
    {
        "id": "otherscan_basic",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=220:OTHERSCAN=rtrange(left=0.5,right=0.5)",
        "description": "Scan with peak at 220 + neighboring scans within ±0.5 min RT",
    },
    # --- Group 16: Explicit TOLERANCEMZ ---
    {
        "id": "explicit_tol_wide",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=220:TOLERANCEMZ=1.0",
        "description": "Wide tolerance 1.0 Da still matches peak at 220",
    },
    {
        "id": "explicit_tol_tight_miss",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=220.06:TOLERANCEMZ=0.05",
        "description": "Tight tolerance 0.05 Da, 220.06 is 0.06 away from 220 — no match",
    },
    # --- Group 17: ANY wildcard ---
    {
        "id": "any_ms2prod",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=ANY",
        "description": "ANY matches all MS2 scans with any peak",
    },
    {
        "id": "any_ms2prod_raw",
        "query": "QUERY MS2DATA WHERE MS2PROD=ANY",
        "description": "ALL MS2 peaks from all scans",
    },
    # --- Group 18: MASSDEFECT exclusion test ---
    {
        "id": "massdefect_exclude_integer",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=220:MASSDEFECT=massdefect(min=0.5,max=0.9)",
        "description": "Integer m/z has defect 0.0, not in (0.5, 0.9) — no match expected",
    },
    # --- Group 19: FILTER with intensity qualifier ---
    {
        "id": "filter_with_intensity",
        "query": "QUERY MS2DATA WHERE MS2PROD=220 FILTER MS2PROD=220:INTENSITYPERCENT>50",
        "description": "WHERE matches scan, FILTER keeps peaks near 220 with high intensity",
    },
    # --- Group 20: peptide() ---
    {
        "id": "peptide_nomatch",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=peptide(ACD,charge=1,ion=b)",
        "description": "Peptide b-ion mass for ACD — no match expected in basic_3spectra",
    },
    # --- Group 21: SCANMIN/SCANMAX (Phase 8) ---
    {
        "id": "meta_scanmin_2",
        "query": "QUERY scannum(MS2DATA) WHERE SCANMIN=2",
        "description": "MS2 scans with scan_number >= 2",
    },
    {
        "id": "meta_scanmax_2",
        "query": "QUERY scannum(MS2DATA) WHERE SCANMAX=2",
        "description": "MS2 scans with scan_number <= 2",
    },
    {
        "id": "meta_scan_range",
        "query": "QUERY scannum(MS2DATA) WHERE SCANMIN=2 AND SCANMAX=3",
        "description": "MS2 scans with scan_number in [2, 3]",
    },
    # --- Group 22: TOLERANCEPPM (Phase 8) ---
    {
        "id": "tol_ppm_ms1",
        "query": "QUERY scannum(MS1DATA) WHERE MS1MZ=200:TOLERANCEPPM=5",
        "description": "MS1 peak at 200 with 5 ppm tolerance",
    },
    # --- Group 23: Phase 9 — PPM on X-offset ---
    {
        "id": "ppm_x_offset",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=X:TOLERANCEPPM=10 AND MS2PROD=X+100:TOLERANCEPPM=10",
        "description": "PPM tolerance on X-offset pattern (integer m/z, same result as Da)",
    },
    # --- Group 24: Phase 9 — MOBILITY ---
    {
        "id": "mobility_wide",
        "query": "QUERY scannum(MS2DATA) WHERE MOBILITY=range(min=0,max=100)",
        "description": "MOBILITY wide range (no-op, returns all MS2 scans)",
    },
    {
        "id": "mobility_narrow",
        "query": "QUERY scannum(MS2DATA) WHERE MOBILITY=range(min=99,max=100)",
        "description": "MOBILITY narrow range still no-op (possible Bug 5)",
    },
    # --- Group 25: Phase 9 — OTHERSCAN+OR ---
    # Note: Python OTHERSCAN is broken on MS2PROD (Bug 2). Our result differs.
    {
        "id": "otherscan_or",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=(150 OR 250):OTHERSCAN=rtrange(left=0.5,right=0.5)",
        "description": "OTHERSCAN with OR value list (Python: Bug 2 ignores OTHERSCAN)",
    },
    # --- Group 26: Phase 9 — X=range() and X=massdefect() ---
    {
        "id": "x_range",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=X AND MS2PROD=X+100 AND X=range(min=100,max=300)",
        "description": "X-offset with X=range constraint",
    },
    {
        "id": "x_massdefect_inclusive",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=X AND MS2PROD=X+100 AND X=massdefect(min=0.0,max=0.5)",
        "description": "X-offset with X=massdefect (inclusive boundary, defect=0.0 passes)",
    },
    {
        "id": "x_massdefect_exclusive",
        "query": "QUERY scannum(MS2DATA) WHERE MS2PROD=X AND MS2PROD=X+100 AND X=massdefect(min=0.1,max=0.9)",
        "description": "X-offset with X=massdefect (defect=0.0 excluded from [0.1, 0.9])",
    },
]


def run_ground_truth():
    if not os.path.exists(MZML_FILE):
        print(f"ERROR: {MZML_FILE} not found. Run from project root.", file=sys.stderr)
        sys.exit(1)

    results = {
        "metadata": {
            "mzml_file": MZML_FILE,
            "massql_version": None,
            "note": "scan numbers are 1-based (Python convention). "
            "spectrum_indices are 0-based (DuckDB convention = scan - 1).",
        },
        "queries": {},
    }

    try:
        import massql

        results["metadata"]["massql_version"] = getattr(massql, "__version__", "unknown")
    except Exception:
        pass

    for entry in QUERIES:
        test_id = entry["id"]
        query = entry["query"]
        description = entry["description"]

        try:
            df = msql_engine.process_query(query, MZML_FILE)
        except Exception as e:
            results["queries"][test_id] = {
                "query": query,
                "description": description,
                "error": str(e),
            }
            print(f"  ERROR {test_id}: {e}", file=sys.stderr)
            continue

        result = {
            "query": query,
            "description": description,
            "row_count": len(df),
        }

        # Scan numbers (1-based Python) and spectrum_indices (0-based DuckDB)
        if "scan" in df.columns:
            scans = sorted(df["scan"].unique().tolist())
            result["scans"] = scans
            result["unique_scan_count"] = len(scans)
            result["spectrum_indices"] = [s - 1 for s in scans]

        # Intensity data
        if "i" in df.columns and len(df) > 0:
            result["sum_intensity"] = round(float(df["i"].sum()), 1)
            result["max_intensity"] = round(float(df["i"].max()), 1)

        # m/z data
        if "mz" in df.columns and len(df) > 0:
            result["mz_values"] = sorted([round(float(v), 4) for v in df["mz"].unique().tolist()])

        # Precursor data (for scanmz)
        if "precmz" in df.columns and len(df) > 0:
            result["precursor_mz_values"] = sorted([round(float(v), 4) for v in df["precmz"].unique().tolist()])

        results["queries"][test_id] = result
        print(
            f"  {test_id:40s} unique_scans={result.get('unique_scan_count', '?'):>2}  "
            f"rows={len(df):3d}  spectrum_indices={result.get('spectrum_indices', '?')}"
        )

    # Write output
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    with open(OUTPUT_FILE, "w") as f:
        json.dump(results, f, indent=2)
        f.write("\n")

    print(f"\nWrote {len(results['queries'])} query results to {OUTPUT_FILE}")


if __name__ == "__main__":
    run_ground_truth()
