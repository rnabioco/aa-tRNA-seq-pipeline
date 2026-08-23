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

tRNA-AA-anticodon-family-species-ref are all preserved from BWA alignment,
and can be further collapsed as desired in downstream analysis

CPM normalization reflects counts per million reads that passed alignment and
the filtering parameters for charging classification; these are full length tRNA
"""

import gzip

import pandas as pd


def per_read_charging(input, output, threshold):
    # Read the TSV file into a DataFrame
    df = pd.read_csv(input, sep="\t")

    # Categorize tRNAs as charged or uncharged
    df["status"] = df["charging_likelihood"].apply(
        lambda x: "counts_charged" if x >= threshold else "counts_uncharged"
    )

    # Group by tRNA and status to get counts
    count_data = df.groupby(["tRNA", "status"]).size().unstack(fill_value=0)

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
