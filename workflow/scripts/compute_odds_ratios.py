"""
Compute per-tRNA pairwise modification odds ratios.

For each tRNA, uses individual reads as the unit of observation to ask:
"among reads of this tRNA, is modification at position X correlated with
modification at position Y (and with charging status)?"

Charging status is represented as position 999.

Inputs:
  - modkit extract calls TSV: per-read, per-position modification calls
  - charging probability TSV: per-read CL tag values

Output:
  - Gzipped TSV with pairwise odds ratios and statistics per tRNA
"""

import argparse
import sys
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, fisher_exact


def load_modkit_calls(path):
    """
    Load modkit extract calls TSV.

    Expected columns include: read_id, chrom, ref_position, call_code
    call_code is "-" for canonical, a letter code for modified.
    """
    df = pd.read_csv(path, sep="\t")
    df["modified"] = (df["call_code"] != "-").astype(int)
    return df


def load_charging(path, ml_threshold):
    """
    Load charging probability TSV and binarize.

    Returns DataFrame with read_id and charged (0/1) columns.
    """
    df = pd.read_csv(path, sep="\t")
    df["charged"] = (df["charging_likelihood"] >= ml_threshold).astype(int)
    return df[["read_id", "charged"]]


def compute_or_for_pair(col_i, col_j):
    """
    Compute odds ratio statistics for a pair of binary columns.

    Drops rows where either value is NaN. Returns None if no valid
    observations, otherwise returns a dict of statistics.
    """
    valid = col_i.notna() & col_j.notna()
    ci = col_i[valid].astype(int).values
    cj = col_j[valid].astype(int).values
    total = len(ci)

    if total == 0:
        return None

    n11 = int(np.sum((ci == 1) & (cj == 1)))
    n10 = int(np.sum((ci == 1) & (cj == 0)))
    n01 = int(np.sum((ci == 0) & (cj == 1)))
    n00 = int(np.sum((ci == 0) & (cj == 0)))

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

    table = np.array([[n11, n10], [n01, n00]])
    fisher_or, p_value = fisher_exact(table)

    return {
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


EMPTY_COLUMNS = [
    "tRNA",
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


def main():
    parser = argparse.ArgumentParser(
        description="Compute per-tRNA pairwise modification odds ratios"
    )
    parser.add_argument(
        "--modkit", required=True, help="Path to modkit extract calls TSV (.gz)"
    )
    parser.add_argument(
        "--charging", required=True, help="Path to charging probability TSV (.gz)"
    )
    parser.add_argument("--output", required=True, help="Output path (.tsv.gz)")
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
        help="Minimum reads per tRNA or position pair (default: 10)",
    )
    args = parser.parse_args()

    # Load data
    modkit_df = load_modkit_calls(args.modkit)
    charging_df = load_charging(args.charging, args.ml_threshold)

    if modkit_df.empty:
        print(
            "WARNING: No modkit calls found. Writing empty output.", file=sys.stderr
        )
        pd.DataFrame(columns=EMPTY_COLUMNS).to_csv(
            args.output, sep="\t", index=False, compression="gzip"
        )
        return

    results = []

    # Process each tRNA separately
    for trna, trna_df in modkit_df.groupby("chrom"):
        # Pivot to per-read matrix: rows=read_id, columns=ref_position, values=modified
        read_matrix = trna_df.pivot_table(
            index="read_id",
            columns="ref_position",
            values="modified",
            aggfunc="max",
        )

        # Skip tRNAs with too few reads
        if len(read_matrix) < args.min_coverage:
            continue

        # Add charging as position 999
        read_matrix = read_matrix.merge(
            charging_df.set_index("read_id"),
            left_index=True,
            right_index=True,
            how="left",
        )
        read_matrix = read_matrix.rename(columns={"charged": 999})

        # Get all column labels (positions + 999 for charging)
        cols = read_matrix.columns.tolist()

        if len(cols) < 2:
            continue

        for col_i, col_j in combinations(cols, 2):
            # Drop reads with NaN at either position
            pair_valid = read_matrix[[col_i, col_j]].dropna()

            if len(pair_valid) < args.min_coverage:
                continue

            stats = compute_or_for_pair(pair_valid[col_i], pair_valid[col_j])
            if stats is None:
                continue

            stats["tRNA"] = trna
            stats["pos1"] = int(col_i)
            stats["pos2"] = int(col_j)
            results.append(stats)

    if len(results) == 0:
        print("WARNING: No valid pairs found. Writing empty output.", file=sys.stderr)
        pd.DataFrame(columns=EMPTY_COLUMNS).to_csv(
            args.output, sep="\t", index=False, compression="gzip"
        )
        return

    results_df = pd.DataFrame(results)

    # BH correction across all results from all tRNAs
    results_df["p_adjusted"] = false_discovery_control(results_df["p_value"].values)

    # Order columns
    results_df = results_df[EMPTY_COLUMNS]

    # Write output
    results_df.to_csv(args.output, sep="\t", index=False, compression="gzip")

    n_sig = (results_df["p_adjusted"] < 0.05).sum()
    n_trnas = results_df["tRNA"].nunique()
    print(
        f"Computed {len(results_df)} pairwise odds ratios across {n_trnas} tRNAs, "
        f"{n_sig} significant (p_adj < 0.05)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
