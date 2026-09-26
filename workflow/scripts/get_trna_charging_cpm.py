#! /usr/bin/env python

"""
Collapses the output of running the CCA charging model and extracting
per-read information on charging likelihood into an ML tag into
per-isodecoder counts (and CPM-normalized counts) of charged and uncharged
tRNAs as determined by the model with a ML >= 200 threshold.

TODO: we should compute this directly from the final BAM file, which
    contains the same charging information in the `CL` tag; no need to write out
    the intermediate per-read charging info. Also, `per_read_charging()`
    does not reflect what this script actuall does. Should be `aggregate_trna_charging()`
    or similar.

tRNA-AA-anticodon-family-species-ref are all preserved from the alignment,
and can be further collapsed as desired in downstream analysis.

A read tied between several references (`tie_refs`, from escpod align's `XA`)
counts 1/(n+1) toward each of its n+1 references rather than 1 toward the
primary alone, so counts are fractional; each read still contributes exactly 1
in total, so a sample's counts sum to its scored read count. A table without a
`tie_refs` column (written before issue #200) counts each read once, as before.

CPM normalization reflects counts per million reads that passed alignment and
the filtering parameters for charging classification; these are full length tRNA
"""

import gzip

import pandas as pd
from get_charging_table import reference_weights


def per_read_charging(input, output, threshold):
    # Read the TSV file into a DataFrame
    df = pd.read_csv(input, sep="\t")

    # Categorize tRNAs as charged or uncharged
    df["status"] = df["charging_likelihood"].apply(
        lambda x: "counts_charged" if x >= threshold else "counts_uncharged"
    )

    # Spread each read over its tie set: one row per (read, reference), each
    # carrying that reference's share of the read.
    if "tie_refs" in df.columns:
        ties = df["tie_refs"].fillna("").astype(str)
    else:
        ties = pd.Series([""] * len(df), index=df.index)
    rows = []
    for ref, tie, status in zip(df["tRNA"], ties, df["status"], strict=True):
        tie_refs = [t for t in tie.split(";") if t]
        for share_ref, weight in reference_weights(ref, tie_refs).items():
            rows.append((share_ref, status, weight))
    shares = pd.DataFrame(rows, columns=["tRNA", "status", "weight"])

    # Group by tRNA and status to get (weighted) counts
    count_data = (
        shares.groupby(["tRNA", "status"])["weight"].sum().unstack(fill_value=0)
    )

    # Ensure both columns exist (handles case where all reads are same status)
    if "counts_charged" not in count_data.columns:
        count_data["counts_charged"] = 0
    if "counts_uncharged" not in count_data.columns:
        count_data["counts_uncharged"] = 0

    # Get total number of reads in the file
    total_reads = len(df)

    # Normalize counts by CPM
    count_data["cpm_charged"] = (count_data["counts_charged"] / total_reads) * 1e6
    count_data["cpm_uncharged"] = (count_data["counts_uncharged"] / total_reads) * 1e6

    # Was opened and never closed: on the gzip branch that risks a truncated
    # file, since the trailer is only written on close.
    open_output = gzip.open if output.endswith(".gz") else open
    with open_output(output, "wt") as output_file:
        count_data.to_csv(output_file, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Process a TSV file to categorize tRNAs as charged or uncharged."
    )
    parser.add_argument(
        "--input", type=str, help="Path to the input TSV file", required=True
    )
    parser.add_argument(
        "--output", type=str, help="Path to the output TSV file", required=True
    )
    parser.add_argument(
        "--ml-threshold",
        type=int,
        default=200,
        help="Threshold for classifying tRNAs as charged (default: 200)",
    )

    # Parse arguments
    args = parser.parse_args()

    # Process the selected file with the given threshold
    per_read_charging(args.input, args.output, args.ml_threshold)
