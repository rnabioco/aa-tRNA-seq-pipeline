#! /usr/bin/env python

"""
Transfer tags from one BAM file to another based on read IDs.

Output only primary alignments with transferred tags.

Use `--rename` to rename tags during transfer.
Use `--all-tags` to transfer every tag from source reads.
Use `--set-tag` to stamp a constant tag onto every output read.

The output header is built from the TARGET plus the SOURCE's @RG lines, rather
than from the target alone. That is not cosmetic: the target is the aligned BAM,
whose header bwa builds fresh from the reference and which therefore declares no
read groups — while `--all-tags` copies dorado's per-read `RG:Z:` straight back
on. Writing with the target's header alone produced a file where every read
referenced an @RG that did not exist: invalid SAM, accepted by `samtools
quickcheck` but rejected by Picard ValidateSamFile and GATK, and with dorado's
basecall-model provenance destroyed at the header while a dangling pointer
remained on each read.
"""

from pysam import AlignmentFile, AlignmentHeader


def build_output_header(
    target, source_header, rg_sample, rg_library, rg_barcode, comments
):
    """Target header, plus the source's @RG lines with our identity spliced in.

    dorado's @RG `ID` is preserved untouched, because that is what the per-read
    `RG:Z:` values point at — rewriting it would re-break the reference this
    exists to repair. `PU`/`PM`/`DT`/`PL`/`DS` (flowcell, device, basecall
    model) ride along. Only SM/LB/BC are overwritten, replacing the sequencing
    run's names with the pipeline's sample, run and barcode.
    """
    header = target.header.to_dict()

    read_groups = [dict(rg) for rg in (source_header or {}).get("RG", [])]
    for read_group in read_groups:
        if rg_sample:
            read_group["SM"] = rg_sample
        if rg_library:
            read_group["LB"] = rg_library
        if rg_barcode:
            read_group["BC"] = rg_barcode
    if read_groups:
        header["RG"] = read_groups

    if comments:
        header.setdefault("CO", []).extend(comments)

    return AlignmentHeader.from_dict(header)


def transfer_tags(
    tags,
    rename,
    source_bam,
    target_bam,
    output_bam,
    all_tags=False,
    threads=1,
    set_tags=None,
    rg_sample=None,
    rg_library=None,
    rg_barcode=None,
    comments=None,
):
    renamed_tags = parse_tag_items(rename)
    constant_tags = parse_set_tags(set_tags or [])

    # Collect target read names (primary only) so we only cache matching source reads
    target_names = set()
    with AlignmentFile(target_bam, "rb", threads=threads) as target:
        for read in target:
            if not read.is_secondary and not read.is_supplementary:
                target_names.add(read.query_name)

    # Single sequential pass through source BAM to cache tags for target reads only
    source_tags = {}
    source_unmapped = set()
    source_header = None
    with AlignmentFile(source_bam, "rb", check_sq=False, threads=threads) as source:
        source_header = source.header.to_dict()
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

    with AlignmentFile(target_bam, "rb", threads=threads) as target:
        out_header = build_output_header(
            target, source_header, rg_sample, rg_library, rg_barcode, comments
        )

    with (
        AlignmentFile(target_bam, "rb", threads=threads) as target,
        AlignmentFile(output_bam, "wb", header=out_header, threads=threads) as output,
    ):
        for read in target:
            if read.is_secondary or read.is_supplementary:
                continue

            read_tags = source_tags.get(read.query_name)

            if read_tags is None:
                if all_tags:
                    apply_constant_tags(read, constant_tags)
                    output.write(read)
                continue

            if not all_tags and read.query_name in source_unmapped:
                continue

            if read_tags:
                for tag, tag_val in read_tags.items():
                    # Unwrap single-element arrays to scalar values
                    # (e.g., ML:B:C:200 → cl:i:200)
                    if (
                        hasattr(tag_val, "__len__")
                        and not isinstance(tag_val, str)
                        and len(tag_val) == 1
                    ):
                        tag_val = tag_val[0]
                    if tag in renamed_tags:
                        read.set_tag(renamed_tags[tag], tag_val)
                    else:
                        read.set_tag(tag, tag_val)

            if all_tags or read_tags:
                apply_constant_tags(read, constant_tags)
                output.write(read)


def parse_tag_items(rename):
    ret = {}
    for item in rename:
        key, val = map(str.strip, item.split("="))
        ret[key] = val
    return ret


def parse_set_tags(items):
    """Parse `TAG:TYPE:VALUE` into (tag, value, value_type) triples.

    Split on the first two colons only: a value may legitimately contain them.
    """
    parsed = []
    for item in items:
        parts = item.split(":", 2)
        if len(parts) != 3:
            raise ValueError(f"--set-tag expects TAG:TYPE:VALUE, got {item!r}")
        tag, value_type, value = parts
        if value_type == "i":
            value = int(value)
        elif value_type == "f":
            value = float(value)
        parsed.append((tag, value, value_type))
    return parsed


def apply_constant_tags(read, constant_tags):
    for tag, value, value_type in constant_tags:
        read.set_tag(tag, value, value_type=value_type)


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
        "--set-tag",
        nargs="+",
        default=[],
        metavar="TAG:TYPE:VALUE",
        help=(
            "Stamp a constant tag onto every output read, e.g. BC:Z:ldx04. "
            "Repeatable. Unlike --tags, the tag need not exist on any source read."
        ),
    )

    parser.add_argument(
        "--rg-sample",
        help="Overwrite SM on the source's @RG lines (the pipeline's sample name)",
    )
    parser.add_argument(
        "--rg-library",
        help="Overwrite LB on the source's @RG lines (typically the run id)",
    )
    parser.add_argument(
        "--rg-barcode",
        help="Overwrite BC on the source's @RG lines (the sample's barcode)",
    )
    parser.add_argument(
        "--comment",
        nargs="+",
        default=[],
        metavar="TEXT",
        help="Add an @CO header comment. Repeatable.",
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
        threads=args.threads,
        set_tags=args.set_tag,
        rg_sample=args.rg_sample,
        rg_library=args.rg_library,
        rg_barcode=args.rg_barcode,
        comments=args.comment,
    )
