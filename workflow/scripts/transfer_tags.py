#! /usr/bin/env python

"""
Transfer tags from one BAM file to another based on read IDs.

Output only primary alignments with transferred tags.

Use `--rename` to rename tags during transfer.
"""

from collections import defaultdict
from pysam import AlignmentFile


def transfer_tags(tags, rename, source_bam, target_bam, output_bam):
    renamed_tags = parse_tag_items(rename)

    with (
        AlignmentFile(source_bam, "rb") as source,
        AlignmentFile(target_bam, "rb") as target,
        AlignmentFile(output_bam, "wb", template=target) as output,
    ):
        # Store tags from the source BAM based on read ID
        source_tags = defaultdict(dict)

        for read in source:
            if not read.is_unmapped:
                for tag in tags:
                    if read.has_tag(tag):
                        source_tags[read.query_name][tag] = read.get_tag(tag)

        # Transfer tags and write only primary alignments with transferred tags
        for read in target:
            # Skip secondary alignments (those marked with the 0x100 flag)
            if read.is_secondary or read.is_supplementary:
                continue

            if read.query_name in source_tags:
                for tag, tag_val in source_tags[read.query_name].items():
                    if tag in renamed_tags:
                        read.set_tag(renamed_tags[tag], tag_val)
                    else:
                        read.set_tag(tag, tag_val)

                # Write read only if tags were transferred
                output.write(read)


def parse_tag_items(rename):
    ret = {}
    for item in rename:
        key, val = map(str.strip, item.split("="))
        ret[key] = val
    return ret


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Transfer tags from one BAM file to another based on read IDs, and output only primary alignments with transferred tags."
    )
    parser.add_argument(
        "-t", "--tags", metavar="MM", nargs="+", required=True, help="Tags to transfer"
    )

    parser.add_argument(
        "--rename",
        nargs="+",
        metavar="OLD=NEW",
        help="tags to rename during transfer",
    )

    parser.add_argument("--source", required=True, help="Source BAM file (with tags)")

    parser.add_argument(
        "--target", required=True, help="Target BAM file (without tags)"
    )
    parser.add_argument(
        "--output", required=True, help="Output BAM file with transferred tags"
    )

    args = parser.parse_args()

    transfer_tags(args.tags, args.rename, args.source, args.target, args.output)
