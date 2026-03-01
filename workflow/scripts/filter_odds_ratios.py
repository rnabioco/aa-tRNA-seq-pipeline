"""
Filter per-sample odds ratio files to keep only well-observed position pairs.

Reads a (potentially large) odds_ratios.tsv.gz in chunks and writes only rows
where total_obs >= min_obs AND p_adjusted < max_p.  This pre-filter matches the
downstream threshold used by clover::filter_linkages() and drastically reduces
file size / R load time.
"""

import argparse
import sys

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Filter odds ratios by minimum observation count."
    )
    parser.add_argument("--input", required=True, help="Input .odds_ratios.tsv.gz")
    parser.add_argument("--output", required=True, help="Output filtered .tsv.gz")
    parser.add_argument(
        "--min-obs",
        type=int,
        default=100,
        help="Minimum total_obs to keep a row (default: 100)",
    )
    parser.add_argument(
        "--max-p",
        type=float,
        default=0.01,
        help="Maximum p_adjusted (BH FDR) to keep a row (default: 0.01)",
    )
    args = parser.parse_args()

    kept = 0
    total = 0
    first_chunk = True

    for chunk in pd.read_csv(
        args.input, sep="\t", compression="gzip", chunksize=100_000
    ):
        total += len(chunk)
        filtered = chunk[
            (chunk["total_obs"] >= args.min_obs)
            & (chunk["p_adjusted"] < args.max_p)
        ]
        kept += len(filtered)
        filtered.to_csv(
            args.output,
            sep="\t",
            index=False,
            compression="gzip",
            mode="w" if first_chunk else "a",
            header=first_chunk,
        )
        first_chunk = False

    print(
        f"Filtered odds ratios: kept {kept:,} / {total:,} rows "
        f"(total_obs >= {args.min_obs}, p_adjusted < {args.max_p})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
