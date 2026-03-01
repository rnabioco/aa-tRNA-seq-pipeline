#! /usr/bin/env python

"""
Transfer tags from one BAM file to another based on read IDs.

Output only primary alignments with transferred tags.

Use `--rename` to rename tags during transfer.
Use `--all-tags` to transfer every tag from source reads.
"""

from pysam import AlignmentFile


def transfer_tags(
    tags, rename, source_bam, target_bam, output_bam, all_tags=False, to_scalar=None, threads=1
):
    renamed_tags = parse_tag_items(rename)
    scalar_tags = set(to_scalar) if to_scalar else set()

    # Collect target read names (primary only) so we only cache matching source reads
    target_names = set()
    with AlignmentFile(target_bam, "rb", threads=threads) as target:
        for read in target:
            if not read.is_secondary and not read.is_supplementary:
                target_names.add(read.query_name)

    # Single sequential pass through source BAM to cache tags for target reads only
    source_tags = {}
    source_unmapped = set()
    with AlignmentFile(source_bam, "rb", check_sq=False, threads=threads) as source:
        for source_read in source:
            name = source_read.query_name
            if name not in target_names or name in source_tags:
                continue  # skip non-target reads; first match wins
            if all_tags:
                source_tags[name] = dict(source_read.get_tags())
            else:
                if source_read.is_unmapped:
                    source_unmapped.add(name)
                tag_dict = {}
                for tag in tags:
                    if source_read.has_tag(tag):
                        tag_dict[tag] = source_read.get_tag(tag)
                source_tags[name] = tag_dict

    with (
        AlignmentFile(target_bam, "rb", threads=threads) as target,
        AlignmentFile(output_bam, "wb", template=target, threads=threads) as output,
    ):
        for read in target:
            if read.is_secondary or read.is_supplementary:
                continue

            read_tags = source_tags.get(read.query_name)

            if read_tags is None:
                if all_tags:
                    output.write(read)
                continue

            if not all_tags and read.query_name in source_unmapped:
                continue

            if read_tags:
                for tag, tag_val in read_tags.items():
                    out_tag = renamed_tags.get(tag, tag)
                    if out_tag in scalar_tags and hasattr(tag_val, "__len__") and len(tag_val) == 1:
                        tag_val = int(tag_val[0])
                    read.set_tag(out_tag, tag_val)

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

    parser.add_argument(
        "--to-scalar",
        nargs="+",
        metavar="TAG",
        help="convert single-element array tags to scalar integers",
    )

    parser.add_argument("--source", required=True, help="Source BAM file (with tags)")

    parser.add_argument(
        "--target", required=True, help="Target BAM file (without tags)"
    )
    parser.add_argument(
        "--output", required=True, help="Output BAM file with transferred tags"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for BAM compression/decompression (default: 1)",
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
        to_scalar=args.to_scalar,
        threads=args.threads,
    )
