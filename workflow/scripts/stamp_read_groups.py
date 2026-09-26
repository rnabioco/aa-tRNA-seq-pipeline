#!/usr/bin/env python
"""
Rewrite a uBAM with the pipeline's identity stamped on, in one pass: the @RG
lines get SM/LB/BC, and on a barcoded sample every record gets a constant
`BC:Z:` tag.

`escpod align` copies the input BAM's @RG/@PG/@CO lines through unchanged and
every record's tags through byte for byte, but it has no way to add a header
line or a tag that was not already on the input. So the sample's identity has
to be on the uBAM before it is aligned, and this is where it is put. It is the
first place both demux backends have converged (see the note at the top of
demux.smk): until here the barcode lives only in the output path.

dorado's @RG `ID` is preserved untouched, because that is what the per-read
`RG:Z:` values point at -- rewriting it would leave every read pointing at an
@RG that does not exist (the dangling-@RG bug PR #121 fixed). `PU`/`PM`/`DT`/
`PL`/`DS` (flowcell, device, basecall model) ride along. Only SM/LB/BC are
overwritten, replacing the sequencing run's names with the pipeline's sample,
run and barcode. An unbarcoded sample gets no BC, on the @RG or on any record:
absence means "no demultiplexing", not "unknown barcode".

The header and the tags are rewritten in the same pass, reader to writer,
because the full-uBAM rewrite is the expensive part. It cannot be streamed
into `escpod align`, which sniffs its input's format and then reopens the path
(a pipe, `-` or `/dev/stdin` all fail on 0.30.0), so escpod_align writes it to a
transient file beside the aligned BAM and deletes it afterwards.
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


def stamp_header(header, sample=None, library=None, barcode=None, comments=()):
    """A copy of a header dict with its @RG stamped and `comments` appended.

    Every other line (@HD, @SQ, @PG, dorado's own @CO) is kept as it was.
    """
    header = dict(header)
    header["RG"] = stamp_read_groups(header.get("RG", []), sample, library, barcode)
    if not header["RG"]:
        del header["RG"]
    if comments:
        header["CO"] = list(header.get("CO", [])) + list(comments)
    return header


def stamp_bam(
    in_path,
    out_path,
    sample=None,
    library=None,
    barcode=None,
    comments=(),
    compress=False,
):
    """Copy `in_path` to `out_path` with a stamped header and, when `barcode`
    is given, a constant `BC:Z:<barcode>` on every record. Returns the number
    of records written.

    `out_path` may be "-" for stdout. Output is uncompressed BAM unless
    `compress` is set, which writes ordinary BGZF. The file is transient
    (escpod_align deletes it once aligned); compressing it is only to cut the
    bytes written to and read back from shared storage, ~4x on dorado uBAM.
    """
    mode = "wb" if compress else "wbu"
    n = 0
    with pysam.AlignmentFile(in_path, "rb", check_sq=False) as src:
        header_dict = src.header.to_dict()
        if not header_dict.get("RG"):
            print(
                f"WARNING: {in_path} declares no @RG; the aligned BAM will carry "
                "none either, and any per-read RG tag on it will dangle.",
                file=sys.stderr,
            )
        header = pysam.AlignmentHeader.from_dict(
            stamp_header(header_dict, sample, library, barcode, comments)
        )
        with pysam.AlignmentFile(out_path, mode, header=header) as dst:
            for read in src.fetch(until_eof=True):
                # Written against dst's header as-is, not re-parsed: htslib
                # serialises the binary record, so the move table and MM/ML
                # arrays are never converted to text and back. A uBAM has no
                # @SQ, so there are no reference ids for the headers to
                # disagree on.
                if barcode:
                    read.set_tag("BC", barcode, value_type="Z")
                dst.write(read)
                n += 1
    return n


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("ubam", help="Unaligned BAM to stamp")
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Output BAM, or - for stdout (uncompressed unless --compress)",
    )
    parser.add_argument("--sample", help="SM: the pipeline's sample name")
    parser.add_argument("--library", help="LB: typically the run id")
    parser.add_argument(
        "--barcode",
        help="BC: the sample's barcode, stamped on the @RG and on every record",
    )
    parser.add_argument(
        "--comment",
        action="append",
        default=[],
        metavar="TEXT",
        help="Add an @CO header line. Repeatable.",
    )
    parser.add_argument(
        "--compress", action="store_true", help="Write BGZF-compressed BAM"
    )
    args = parser.parse_args()

    n = stamp_bam(
        args.ubam,
        args.output,
        sample=args.sample,
        library=args.library,
        barcode=args.barcode,
        comments=args.comment,
        compress=args.compress,
    )
    print(f"stamped {n} records", file=sys.stderr)


if __name__ == "__main__":
    main()
