#! /usr/bin/env python

"""
Transfer ML/MM tags from one BAM file to another based on read IDs.
Output only primary alignments with transferred tags."
"""

import pysam


def transfer_tags(source_bam, target_bam, output_bam):
    with (
        pysam.AlignmentFile(source_bam, "rb") as source,
        pysam.AlignmentFile(target_bam, "rb") as target,
        pysam.AlignmentFile(output_bam, "wb", template=target) as output,
    ):
        # Create a dictionary to store ML and MM tags from the source BAM based on read IDs
        source_tags = {}
        for read in source:
            if not read.is_unmapped:
                # Extract ML and MM tags
                ml_tag = read.get_tag("ML") if read.has_tag("ML") else None
                mm_tag = read.get_tag("MM") if read.has_tag("MM") else None
                if ml_tag or mm_tag:
                    source_tags[read.query_name] = (ml_tag, mm_tag)

        # Transfer tags and write only primary alignments with transferred tags
        for read in target:
            # Skip secondary alignments (those marked with the 0x100 flag)
            if read.is_secondary or read.is_supplementary:
                continue

            if read.query_name in source_tags:
                ml_tag, mm_tag = source_tags[read.query_name]
                # Add ML and MM tags to the target read if they exist
                if mm_tag is not None:
                    read.set_tag("MM", mm_tag)
                if ml_tag is not None:
                    read.set_tag("ML", ml_tag)
                # Write the read to the output BAM only if tags were transferred
                output.write(read)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Transfer ML/MM tags from one BAM file to another based on read IDs, and output only primary alignments with transferred tags."
    )
    parser.add_argument(
        "-s", "--source", required=True, help="Source BAM file (with ML/MM tags)"
    )
    parser.add_argument(
        "-t", "--target", required=True, help="Target BAM file (without ML/MM tags)"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output BAM file with transferred tags"
    )

    args = parser.parse_args()

    # Call the function to transfer tags
    transfer_tags(args.source, args.target, args.output)
