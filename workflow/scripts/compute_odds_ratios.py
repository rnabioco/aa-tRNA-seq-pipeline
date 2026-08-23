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


def read_fasta_lengths(fasta_path):
    """
    Read FASTA file and return dict of {name: sequence_length}.
    """
    lengths = {}
    name = None
    seq_len = 0

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = seq_len
                name = line[1:].split()[0]
                seq_len = 0
            else:
                seq_len += len(line)

        if name is not None:
            lengths[name] = seq_len

    return lengths


def load_modkit_calls_chunked(path, offset_5p=0, offset_3p=0, ref_lengths=None):
    """
    Load modkit extract calls TSV in chunks, selecting only needed columns.

    Reads in 500K-row chunks, filters adapter positions, converts coordinates,
    and accumulates per-tRNA DataFrames in a dict to avoid holding the full
    dataset in memory.

    Returns dict of {tRNA_name: DataFrame} with columns:
        read_id, chrom, ref_position, modified
    """
    CHUNK_SIZE = 500_000
    USECOLS = ["read_id", "chrom", "ref_position", "call_code"]
    DTYPES = {
        "read_id": str,
        "chrom": "category",
        "ref_position": "int32",
        "call_code": "category",
    }

    trna_chunks = {}

    for chunk in pd.read_csv(
        path, sep="\t", usecols=USECOLS, dtype=DTYPES, chunksize=CHUNK_SIZE
    ):
        chunk["modified"] = (chunk["call_code"] != "-").astype("int8")
        chunk.drop(columns=["call_code"], inplace=True)

        # Filter adapter positions
        if offset_5p > 0:
            chunk = chunk[chunk["ref_position"] >= offset_5p]

        if offset_3p > 0 and ref_lengths:
            max_pos = chunk["chrom"].map(ref_lengths)
            chunk = chunk[chunk["ref_position"] < max_pos - offset_3p]

        if chunk.empty:
            continue

        # Convert to 1-indexed tRNA coordinates
        if offset_5p > 0:
            chunk["ref_position"] = chunk["ref_position"] - offset_5p + 1

        # Accumulate per-tRNA
        for trna, trna_df in chunk.groupby("chrom", observed=True):
            if trna in trna_chunks:
                trna_chunks[trna].append(trna_df)
            else:
                trna_chunks[trna] = [trna_df]

    # Concatenate accumulated chunks per tRNA
    result = {}
    for trna, chunks in trna_chunks.items():
        result[trna] = pd.concat(chunks, ignore_index=True)

    return result


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
    parser.add_argument(
        "--offset-5p",
        type=int,
        default=0,
        help="Number of bases to skip at 5' end of each reference (adapter + N). "
        "Positions are renumbered to 1-indexed tRNA coordinates.",
    )
    parser.add_argument(
        "--offset-3p",
        type=int,
        default=0,
        help="Number of bases to skip at 3' end of each reference (3' adapter).",
    )
    parser.add_argument(
        "--reference",
        default=None,
        help="Reference FASTA file (required when --offset-3p > 0 to determine "
        "per-tRNA reference lengths).",
    )
    args = parser.parse_args()

    if args.offset_3p > 0 and args.reference is None:
        parser.error("--reference is required when --offset-3p > 0")

    # Load reference lengths for 3' adapter filtering
    ref_lengths = {}
    if args.reference:
        ref_lengths = read_fasta_lengths(args.reference)

    # Load data: chunked reader returns dict of {tRNA: DataFrame}
    trna_dfs = load_modkit_calls_chunked(
        args.modkit, args.offset_5p, args.offset_3p, ref_lengths
    )
    charging_df = load_charging(args.charging, args.ml_threshold)

    if not trna_dfs:
        print("WARNING: No modkit calls found. Writing empty output.", file=sys.stderr)
        pd.DataFrame(columns=EMPTY_COLUMNS).to_csv(
            args.output, sep="\t", index=False, compression="gzip"
        )
        return

    results = []

    # Process each tRNA separately; discard after processing to free memory
    for trna, trna_df in sorted(trna_dfs.items()):
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
