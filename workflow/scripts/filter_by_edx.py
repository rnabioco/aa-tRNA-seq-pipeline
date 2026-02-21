#!/usr/bin/env python3
"""
Filter a BAM file by EDX (3' adapter barcode) identity.

Reads the PT tag from each alignment and keeps only reads whose
3' adapter name matches the expected EDX value.
"""

import argparse
import sys

import pysam


def parse_3p_adapter_from_pt(pt_tag):
    """
    Extract 3' adapter name from a PT tag string.

    PT tag format: "start;end;strand;type|start;end;strand;type"
    3' adapter entries look like:
      - "3p_adapter"       -> "default"
      - "3p_adapter_edx1"  -> "edx1"
      - "3p_adapter_v2"    -> "v2"

    Returns the adapter name or None if no 3' adapter found.
    """
    if not pt_tag:
        return None

    for segment in pt_tag.split("|"):
        fields = segment.split(";")
        if len(fields) < 4:
            continue
        entry_type = fields[3]
        if entry_type.startswith("3p_adapter"):
            suffix = entry_type[len("3p_adapter"):]
            if suffix.startswith("_"):
                return suffix[1:]  # e.g., "edx1", "v2"
            elif suffix == "":
                return "default"
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bam",
        required=True,
        help="Input BAM file (with PT tags from add_adapter_tags)",
    )
    parser.add_argument(
        "--edx",
        required=True,
        help="Expected EDX adapter name (e.g., 'edx1', 'v2')",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output BAM file (filtered)",
    )
    args = parser.parse_args()

    total = 0
    kept = 0
    no_pt = 0

    with pysam.AlignmentFile(args.bam, "rb") as bam_in:
        with pysam.AlignmentFile(args.output, "wb", header=bam_in.header) as bam_out:
            for read in bam_in.fetch(until_eof=True):
                total += 1
                try:
                    pt_tag = read.get_tag("PT")
                except KeyError:
                    no_pt += 1
                    continue

                adapter = parse_3p_adapter_from_pt(pt_tag)
                if adapter == args.edx:
                    bam_out.write(read)
                    kept += 1

    filtered = total - kept
    print(
        f"EDX filter (expected={args.edx}): "
        f"{total} total, {kept} kept, {filtered} filtered "
        f"({no_pt} had no PT tag)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
