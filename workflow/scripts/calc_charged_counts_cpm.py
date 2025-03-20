#! /usr/bin/env python

"""
Summary table for charging classification.

Collapses the output of running a Remora CCA model and extracting
per-read information on charging likelihood into an ML tag into
per-isodecoder counts (and CPM-normalized counts) of charged and uncharged
tRNAs as determined by the model with a ML >= 200 threshold.

CPM normalization reflects counts per million reads that passed alignment and
the filtering parameters for Remora classification; these are full length tRNA.
"""

from collections import Counter, defaultdict
from pysam import AlignmentFile

def calc_charging_table(bam_file, tag, threshold):

    charge_counts = defaultdict(Counter)

    with AlignmentFile(bam_file, "rb") as bam:

       for read in bam:
            tag_val = read.get_tag(tag)
            if tag_val >= threshold:
                charge_counts[read.reference]["charged"] += 1
            else:
                charge_counts[read.reference]["uncharged"] += 1

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Calculate charing information from a BAM file."
    )
    parser.add_argument(
        "--bam-file", type=str, help="Path to BAM file", required=True
    )

    parser.add_argument(
        "-t", "--charging-tag",
        type=str,
        default="CL",
        help="Tag containing charing probability"
    )

    parser.add_argument(
        "-c", "--charging-threshold",
        type=int,
        default=200,
        help="Threshold for classifying tRNAs as charged",
    )

    # Parse arguments
    args = parser.parse_args()

    # Process the selected file with the given threshold
    calc_charging_table(args.bam_file, args.charging_tag, args.charging_threshold)
