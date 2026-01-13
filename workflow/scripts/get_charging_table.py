#! /usr/bin/env python

"""
Generate table of read id, ref, value of charging tag
"""

import pysam
import argparse
import csv
import gzip


def extract_tag(bam_file, output_tsv, tag):
    open_func = gzip.open if output_tsv.endswith(".gz") else open
    mode = "wt" if output_tsv.endswith(".gz") else "w"

    with (
        pysam.AlignmentFile(bam_file, "rb") as bam,
        open_func(output_tsv, mode) as tsvfile,
    ):
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["read_id", "tRNA", "charging_likelihood"])

        for read in bam.fetch():
            read_id = read.query_name
            reference = read.reference_name if read.reference_name else "*"
            tag_array = dict(read.tags).get(tag, None)

            # XXX: handle case where there are more than 1 tag value
            # not clear why this is, but we skip for now as it's a small
            # number of reads affected
            if len(tag_array) > 1:
                continue

            tag_value = tag_array[0]

            if tag_value and reference != "*":
                writer.writerow([read_id, reference, tag_value])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract a specified tag from a BAM file and write to TSV."
    )
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument(
        "output_tsv", help="Output TSV file (can be .gz for compression)"
    )
    parser.add_argument("--tag", default="ML", help="BAM tag to extract (default: ML)")

    args = parser.parse_args()
    extract_tag(args.bam_file, args.output_tsv, args.tag)
