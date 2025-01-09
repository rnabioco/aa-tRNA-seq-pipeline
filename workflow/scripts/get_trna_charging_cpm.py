#! /usr/bin/env python

"""
and collapses the output of running a Remora CCA model and extracting
per-read information on charging likelihood into an ML tag into
per-isodecoder counts (and CPM-normalized counts) of charged and uncharged
tRNAs as determined by the model with a ML >= 200 threshold.

tRNA-AA-anticodon-family-species-chargingref are all preserved from BWA alignment,
and can be further collapsed as desired in downstream analysis

CPM normalization reflects counts per million reads that passed alingnment and
the filtering parameters for Remora classification; these are full length tRNA
"""

import pandas as pd
import gzip


def per_read_charging(input, output, threshold):
    # Read the TSV file into a DataFrame
    df = pd.read_csv(input, sep="\t")

    # Categorize tRNAs as charged or uncharged
    df["status"] = df["charging_likelihood"].apply(
        lambda x: "charged" if x >= threshold else "uncharged"
    )

    # Group by tRNA and status to get counts
    count_data = df.groupby(["tRNA", "status"]).size().unstack(fill_value=0)

    # Get total number of reads in the file
    total_reads = len(df)

    # Normalize counts by CPM
    count_data["charged_CPM"] = (count_data["charged"] / total_reads) * 1e6
    count_data["uncharged_CPM"] = (count_data["uncharged"] / total_reads) * 1e6

    if output.endswith(".gz"):
        output_file = gzip.open(output, "wt")
    else:
        output_file = open(output, "w")

    # Write the results to a new file
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
