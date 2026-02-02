"""
Compute pairwise modification odds ratios across tRNA references.

For each pair of positions (including charging status), builds a 2x2
contingency table across all tRNA genes in the sample and computes
odds ratios with Fisher's exact test and BH-corrected p-values.

Inputs:
  - bcerror TSV (from get_bcerror_freqs.py): per-position error metrics
  - charging probability TSV (from get_charging_table.py): per-read CL tag values

Output:
  - Gzipped TSV with pairwise odds ratios and statistics
"""

import argparse
import sys
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, fisher_exact


def binarize_modifications(bcerror_df, mod_threshold, min_coverage):
    """
    Binarize modification calls from bcerror mismatch frequencies.

    Filters positions by minimum coverage, then marks each
    (reference, position) as modified (1) or unmodified (0).

    Returns a DataFrame with columns: Reference, Position, modified
    """
    df = bcerror_df[bcerror_df["Spanning_Reads"] >= min_coverage].copy()
    df["modified"] = (df["MismatchFreq"] >= mod_threshold).astype(int)
    return df[["Reference", "Position", "modified"]]


def build_modification_matrix(mod_df):
    """
    Build binary matrix: rows = position labels, columns = tRNA references.

    Returns (matrix DataFrame, list of position labels).
    Position labels are formatted as "pos_{Position}" from the bcerror data.
    """
    mod_df = mod_df.copy()
    mod_df["pos_label"] = "pos_" + mod_df["Position"].astype(str)

    matrix = mod_df.pivot_table(
        index="pos_label", columns="Reference", values="modified", aggfunc="max"
    )
    return matrix


def compute_charging_row(charging_df, ml_threshold):
    """
    Aggregate per-read charging probabilities to per-reference binary calls.

    A tRNA reference is called 'charged' (1) if >= 50% of its reads
    have CL tag >= ml_threshold.

    Returns a Series indexed by tRNA reference name.
    """
    df = charging_df.copy()
    df["is_charged"] = (df["charging_likelihood"] >= ml_threshold).astype(int)
    frac_charged = df.groupby("tRNA")["is_charged"].mean()
    return (frac_charged >= 0.5).astype(int)


def compute_pairwise_odds_ratios(matrix):
    """
    Compute odds ratios for all pairs of rows in the binary matrix.

    For each pair (pos_i, pos_j), builds a 2x2 contingency table:
        n00 = both 0, n01 = i=0 & j=1, n10 = i=1 & j=0, n11 = both 1

    Applies Haldane correction (add 0.5) when any cell is zero.
    Computes OR, log(OR), SE, 95% CI, and Fisher's exact test.

    Returns a list of result dicts.
    """
    labels = matrix.index.tolist()
    results = []

    for label_i, label_j in combinations(labels, 2):
        row_i = matrix.loc[label_i].values
        row_j = matrix.loc[label_j].values

        # Only use columns (tRNA refs) where both rows have data
        valid = ~(np.isnan(row_i) | np.isnan(row_j))
        ri = row_i[valid].astype(int)
        rj = row_j[valid].astype(int)

        if len(ri) == 0:
            continue

        n11 = int(np.sum((ri == 1) & (rj == 1)))
        n10 = int(np.sum((ri == 1) & (rj == 0)))
        n01 = int(np.sum((ri == 0) & (rj == 1)))
        n00 = int(np.sum((ri == 0) & (rj == 0)))
        total = n00 + n01 + n10 + n11

        # Haldane correction when any cell is zero
        if n00 == 0 or n01 == 0 or n10 == 0 or n11 == 0:
            a, b, c, d = n11 + 0.5, n10 + 0.5, n01 + 0.5, n00 + 0.5
        else:
            a, b, c, d = n11, n10, n01, n00

        odds_ratio = (a * d) / (b * c)
        log_or = np.log(odds_ratio)
        se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        ci_lower = np.exp(log_or - 1.96 * se)
        ci_upper = np.exp(log_or + 1.96 * se)

        # Fisher's exact test on original counts
        table = np.array([[n11, n10], [n01, n00]])
        fisher_or, p_value = fisher_exact(table)

        results.append(
            {
                "pos1": label_i,
                "pos2": label_j,
                "n00": n00,
                "n01": n01,
                "n10": n10,
                "n11": n11,
                "total_obs": total,
                "odds_ratio": odds_ratio,
                "log_odds_ratio": log_or,
                "se_log_or": se,
                "ci_lower": ci_lower,
                "ci_upper": ci_upper,
                "fisher_or": fisher_or,
                "p_value": p_value,
            }
        )

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Compute pairwise modification odds ratios across tRNA references"
    )
    parser.add_argument("--bcerror", required=True, help="Path to bcerror TSV (.gz)")
    parser.add_argument(
        "--charging", required=True, help="Path to charging probability TSV (.gz)"
    )
    parser.add_argument("--output", required=True, help="Output path (.tsv.gz)")
    parser.add_argument(
        "--mod-threshold",
        type=float,
        default=0.3,
        help="Mismatch frequency threshold for calling modified (default: 0.3)",
    )
    parser.add_argument(
        "--ml-threshold",
        type=int,
        default=200,
        help="CL tag threshold for charged classification (default: 200)",
    )
    parser.add_argument(
        "--min-coverage",
        type=int,
        default=10,
        help="Minimum spanning reads per position (default: 10)",
    )
    args = parser.parse_args()

    # Read input data
    bcerror_df = pd.read_csv(args.bcerror, sep="\t")
    charging_df = pd.read_csv(args.charging, sep="\t")

    # Binarize modifications and build matrix
    mod_df = binarize_modifications(bcerror_df, args.mod_threshold, args.min_coverage)
    matrix = build_modification_matrix(mod_df)

    # Add charging row
    charging_row = compute_charging_row(charging_df, args.ml_threshold)

    # Align charging row to matrix columns (tRNA references present in both)
    common_refs = matrix.columns.intersection(charging_row.index)
    if len(common_refs) == 0:
        print(
            "WARNING: No common tRNA references between bcerror and charging data. "
            "Writing empty output.",
            file=sys.stderr,
        )
        empty_df = pd.DataFrame(
            columns=[
                "pos1",
                "pos2",
                "n00",
                "n01",
                "n10",
                "n11",
                "total_obs",
                "odds_ratio",
                "log_odds_ratio",
                "se_log_or",
                "ci_lower",
                "ci_upper",
                "fisher_or",
                "p_value",
                "p_adjusted",
            ]
        )
        empty_df.to_csv(args.output, sep="\t", index=False, compression="gzip")
        return

    matrix = matrix[common_refs]
    charging_series = charging_row[common_refs]
    charging_frame = pd.DataFrame(
        [charging_series.values], index=["charging"], columns=common_refs
    )
    matrix = pd.concat([matrix, charging_frame])

    # Compute pairwise odds ratios
    results = compute_pairwise_odds_ratios(matrix)

    if len(results) == 0:
        print("WARNING: No valid pairs found. Writing empty output.", file=sys.stderr)
        empty_df = pd.DataFrame(
            columns=[
                "pos1",
                "pos2",
                "n00",
                "n01",
                "n10",
                "n11",
                "total_obs",
                "odds_ratio",
                "log_odds_ratio",
                "se_log_or",
                "ci_lower",
                "ci_upper",
                "fisher_or",
                "p_value",
                "p_adjusted",
            ]
        )
        empty_df.to_csv(args.output, sep="\t", index=False, compression="gzip")
        return

    results_df = pd.DataFrame(results)

    # BH correction for multiple testing
    results_df["p_adjusted"] = false_discovery_control(results_df["p_value"].values)

    # Write output
    results_df.to_csv(args.output, sep="\t", index=False, compression="gzip")

    n_sig = (results_df["p_adjusted"] < 0.05).sum()
    print(
        f"Computed {len(results_df)} pairwise odds ratios, "
        f"{n_sig} significant (p_adj < 0.05)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
