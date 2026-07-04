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
            tags_dict = dict(read.tags)
            tag_raw = tags_dict.get(tag, None)

            # Fallback to uppercase tag for backward compat with older BAMs
            # TODO: remove fallback once all BAMs have been reprocessed
            if tag_raw is None and tag.islower():
                tag_raw = tags_dict.get(tag.upper(), None)

            if tag_raw is None:
                continue

            # Handle both scalar (cl:i:200) and array (CL:B:C:200) tag values
            if hasattr(tag_raw, "__len__") and not isinstance(tag_raw, str):
                if len(tag_raw) > 1:
                    continue
                tag_value = tag_raw[0]
            else:
                tag_value = tag_raw

            # Write on tag PRESENCE, not truthiness: a charging tag of 0 is a
            # valid, maximally-confident *uncharged* call (ML score range is
            # 0-255, >=200 = charged). `if tag_value` would silently drop
            # ML==0 reads, biasing charging fraction upward and shrinking the
            # CPM denominator downstream.
            if tag_value is not None and reference != "*":
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
