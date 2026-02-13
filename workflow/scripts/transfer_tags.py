#! /usr/bin/env python

"""
Transfer tags from one BAM file to another based on read IDs.

Output only primary alignments with transferred tags.

Use `--rename` to rename tags during transfer.
Use `--all-tags` to transfer every tag from source reads.
"""

from pysam import AlignmentFile, IndexedReads


def transfer_tags(tags, rename, source_bam, target_bam, output_bam, all_tags=False):
    renamed_tags = parse_tag_items(rename)

    with (
        AlignmentFile(source_bam, "rb", check_sq=False) as source,
        AlignmentFile(target_bam, "rb") as target,
        AlignmentFile(output_bam, "wb", template=target) as output,
    ):
        # Build read-name index for random access into source BAM
        source_idx = IndexedReads(source)
        source_idx.build()

        for read in target:
            if read.is_secondary or read.is_supplementary:
                continue

            # Look up source read by name
            try:
                source_reads = source_idx.find(read.query_name)
            except KeyError:
                if all_tags:
                    output.write(read)
                continue

            # Get tags from the first matching source read
            source_read = next(source_reads)
            if not all_tags and source_read.is_unmapped:
                continue

            if all_tags:
                read_tags = dict(source_read.get_tags())
            else:
                read_tags = {}
                for tag in tags:
                    if source_read.has_tag(tag):
                        read_tags[tag] = source_read.get_tag(tag)

            if read_tags:
                for tag, tag_val in read_tags.items():
                    if tag in renamed_tags:
                        read.set_tag(renamed_tags[tag], tag_val)
                    else:
                        read.set_tag(tag, tag_val)

            if all_tags or read_tags:
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
        "-t", "--tags", metavar="MM", nargs="+", help="Tags to transfer"
    )

    parser.add_argument(
        "--all-tags",
        action="store_true",
        help="Transfer all tags from source reads (ignores --tags)",
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

    if not args.all_tags and not args.tags:
        parser.error("either --tags or --all-tags is required")

    transfer_tags(
        args.tags or [],
        args.rename or [],
        args.source,
        args.target,
        args.output,
        all_tags=args.all_tags,
    )
