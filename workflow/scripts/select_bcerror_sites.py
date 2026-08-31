import argparse

import pandas as pd

"""
This script selects candidate modification sites from base-calling error rates.

Given the per-position bcerror tables produced by get_bcerror_freqs.py for one
or more samples, it keeps positions whose error rate and coverage clear the
given thresholds in at least --min-samples samples, and writes the union as a
site list.

Selecting sites once across all samples, rather than per sample, means the
downstream read-level analysis tests the same positions everywhere and so can
be compared between samples without a further intersection step.

The output is a TSV with columns:

    ref  pos  n_samples  max_error  mean_error

which is accepted by get_mismatch_calls.py --sites and by clover's
compute_charging_odds_ratios(sites = ).

Positions are passed through unchanged, so they inherit the coordinate system
of the input bcerror tables (tRNA-only and 1-indexed when those were built
with --offset-5p/--offset-3p).

Example:
    python select_bcerror_sites.py --min-error 0.1 --min-cov 20 \
        -o sites.tsv.gz sampleA.bcerror.tsv.gz sampleB.bcerror.tsv.gz
"""


def select_sites(paths, min_error=0.1, min_cov=20, min_samples=1):
    frames = []

    for path in paths:
        df = pd.read_csv(
            path,
            sep="\t",
            usecols=["Reference", "Position", "Spanning_Reads", "BCErrorFreq"],
        )
        df = df[(df["Spanning_Reads"] >= min_cov) & (df["BCErrorFreq"] >= min_error)]
        frames.append(df)

    if not frames:
        return pd.DataFrame(
            columns=["ref", "pos", "n_samples", "max_error", "mean_error"]
        )

    combined = pd.concat(frames, ignore_index=True)

    sites = (
        combined.groupby(["Reference", "Position"], as_index=False)
        .agg(
            n_samples=("BCErrorFreq", "size"),
            max_error=("BCErrorFreq", "max"),
            mean_error=("BCErrorFreq", "mean"),
        )
        .rename(columns={"Reference": "ref", "Position": "pos"})
    )

    sites = sites[sites["n_samples"] >= min_samples]

    return sites.sort_values(["ref", "pos"]).reset_index(drop=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Select candidate modification sites from base-calling error"
    )

    parser.add_argument("bcerror_tsv", nargs="+", help="Paths to bcerror TSV files")
    parser.add_argument("-o", "--output", required=True, help="Path for output TSV")
    parser.add_argument(
        "--min-error",
        type=float,
        default=0.1,
        help="Minimum BCErrorFreq for a position to be called (default: 0.1)",
    )
    parser.add_argument(
        "--min-cov",
        type=int,
        default=20,
        help="Minimum spanning-read coverage for a position (default: 20)",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=1,
        help="Number of samples a position must clear the thresholds in (default: 1)",
    )
    args = parser.parse_args()

    sites = select_sites(
        args.bcerror_tsv,
        min_error=args.min_error,
        min_cov=args.min_cov,
        min_samples=args.min_samples,
    )

    compression = "gzip" if args.output.endswith(".gz") else None
    sites.to_csv(args.output, sep="\t", index=False, compression=compression)

    print(f"selected {len(sites)} sites from {len(args.bcerror_tsv)} sample(s)")
