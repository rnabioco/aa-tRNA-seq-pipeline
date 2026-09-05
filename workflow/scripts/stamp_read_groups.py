#!/usr/bin/env python
"""
Emit a uBAM's @RG lines, with the pipeline's identity stamped on, as SAM header
text for `bwa mem -H`.

bwa builds the aligned header fresh from the reference and declares no read
groups, while `bwa mem -C` copies dorado's per-read `RG:Z:` straight through
from the FASTQ comment. Without this, every read in the aligned BAM would point
at an @RG that does not exist: invalid SAM, accepted by `samtools quickcheck`
but rejected by Picard ValidateSamFile and GATK, with dorado's basecall-model
provenance destroyed at the header while a dangling pointer remained on each
read. (That was the bug PR #121 fixed in the transfer_tags.py this replaces.)

dorado's @RG `ID` is preserved untouched, because that is what the per-read
`RG:Z:` values point at -- rewriting it would re-break the reference this
exists to repair. `PU`/`PM`/`DT`/`PL`/`DS` (flowcell, device, basecall model)
ride along. Only SM/LB/BC are overwritten, replacing the sequencing run's names
with the pipeline's sample, run and barcode. An unbarcoded sample gets no BC:
absence means "no demultiplexing", not "unknown barcode".
"""

import argparse
import sys

import pysam


def stamp_read_groups(read_groups, sample=None, library=None, barcode=None):
    """Return copies of `read_groups` with SM/LB/BC overwritten where given."""
    stamped = []
    for read_group in read_groups:
        read_group = dict(read_group)
        if sample:
            read_group["SM"] = sample
        if library:
            read_group["LB"] = library
        if barcode:
            read_group["BC"] = barcode
        stamped.append(read_group)
    return stamped


def format_header_lines(read_groups, comments=()):
    """SAM header text: one @RG line per read group, then one @CO per comment."""
    lines = []
    for read_group in read_groups:
        fields = ["@RG"] + [f"{tag}:{value}" for tag, value in read_group.items()]
        lines.append("\t".join(fields))
    for comment in comments:
        lines.append(f"@CO\t{comment}")
    return lines


def read_groups_from_bam(path):
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
        return bam.header.to_dict().get("RG", [])


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("ubam", help="Unaligned BAM whose @RG lines to carry over")
    parser.add_argument("--sample", help="SM: the pipeline's sample name")
    parser.add_argument("--library", help="LB: typically the run id")
    parser.add_argument("--barcode", help="BC: the sample's barcode")
    parser.add_argument(
        "--comment",
        action="append",
        default=[],
        metavar="TEXT",
        help="Add an @CO header line. Repeatable.",
    )
    parser.add_argument("--output", help="Write here instead of stdout", default=None)
    args = parser.parse_args()

    read_groups = read_groups_from_bam(args.ubam)
    if not read_groups:
        print(
            f"WARNING: {args.ubam} declares no @RG; the aligned BAM will carry "
            "none either, and any per-read RG tag on it will dangle.",
            file=sys.stderr,
        )
    lines = format_header_lines(
        stamp_read_groups(read_groups, args.sample, args.library, args.barcode),
        args.comment,
    )

    text = "".join(line + "\n" for line in lines)
    if args.output:
        with open(args.output, "w") as out:
            out.write(text)
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
