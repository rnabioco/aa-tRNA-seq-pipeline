"""
Convert modkit/coverage outputs from full-reference coordinates to tRNA-only coordinates.

Filters out positions in 5' and 3' adapter regions and shifts remaining
positions so that the first tRNA base is position 0 (BED formats) or
position 1 (TSV formats, matching bcerror/odds_ratios conventions).

Supports four output formats via --format:
  bedgraph    - 4-column BED: chrom, start, end, value
  bedmethyl   - modkit pileup bedMethyl (BED-like, 0-based half-open)
  modkit_calls - modkit extract calls TSV (header, ref_position col 2)
  modkit_full  - modkit extract full TSV (header, ref_position col 2)
"""

import argparse
import gzip
import sys


def read_fasta_lengths(fasta_path):
    """Read FASTA file and return dict of {name: sequence_length}."""
    lengths = {}
    name = None
    seq_len = 0

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = seq_len
                name = line[1:].split()[0]
                seq_len = 0
            else:
                seq_len += len(line)

        if name is not None:
            lengths[name] = seq_len

    return lengths


def process_bedgraph(infile, outfile, offset_5p, offset_3p, ref_lengths):
    """Convert bedGraph: filter adapter positions, shift to tRNA coords (0-based)."""
    for line in infile:
        fields = line.rstrip("\n").split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])

        ref_len = ref_lengths.get(chrom)
        if ref_len is None:
            continue

        max_pos = ref_len - offset_3p
        if start < offset_5p or end > max_pos:
            continue

        fields[1] = str(start - offset_5p)
        fields[2] = str(end - offset_5p)
        outfile.write(("\t".join(fields) + "\n").encode())


def process_bedmethyl(infile, outfile, offset_5p, offset_3p, ref_lengths):
    """Convert bedMethyl: filter adapter positions, shift start/end/thickStart/thickEnd."""
    for line in infile:
        fields = line.rstrip("\n").split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])

        ref_len = ref_lengths.get(chrom)
        if ref_len is None:
            continue

        max_pos = ref_len - offset_3p
        if start < offset_5p or end > max_pos:
            continue

        fields[1] = str(start - offset_5p)
        fields[2] = str(end - offset_5p)
        # thickStart (col 6) and thickEnd (col 7)
        if len(fields) > 7:
            fields[6] = str(int(fields[6]) - offset_5p)
            fields[7] = str(int(fields[7]) - offset_5p)
        outfile.write(("\t".join(fields) + "\n").encode())


def process_modkit_tsv(infile, outfile, offset_5p, offset_3p, ref_lengths):
    """Convert modkit extract calls/full TSV: filter and shift ref_position to 1-indexed."""
    header = next(infile)
    outfile.write(header.encode() if isinstance(header, str) else header)

    for line in infile:
        fields = line.rstrip("\n").split("\t")
        chrom = fields[3]  # chrom is col 3
        ref_pos = int(fields[2])  # ref_position is col 2

        ref_len = ref_lengths.get(chrom)
        if ref_len is None:
            continue

        max_pos = ref_len - offset_3p
        if ref_pos < offset_5p or ref_pos >= max_pos:
            continue

        # Convert to 1-indexed tRNA coordinate
        fields[2] = str(ref_pos - offset_5p + 1)
        outfile.write(("\t".join(fields) + "\n").encode())


def main():
    parser = argparse.ArgumentParser(
        description="Convert modkit/coverage outputs to tRNA-only coordinates"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input file path or '-' for stdin",
    )
    parser.add_argument("--output", required=True, help="Output gzipped file path")
    parser.add_argument(
        "--format",
        required=True,
        choices=["bedgraph", "bedmethyl", "modkit_calls", "modkit_full"],
        help="Input format",
    )
    parser.add_argument(
        "--offset-5p",
        type=int,
        required=True,
        help="5' adapter + N offset to subtract",
    )
    parser.add_argument(
        "--offset-3p",
        type=int,
        required=True,
        help="3' adapter length to trim",
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA for per-tRNA lengths",
    )
    args = parser.parse_args()

    ref_lengths = read_fasta_lengths(args.reference)

    if args.input == "-":
        infile = sys.stdin
    else:
        infile = open(args.input, "r")

    outfile = gzip.open(args.output, "wb")

    try:
        processors = {
            "bedgraph": process_bedgraph,
            "bedmethyl": process_bedmethyl,
            "modkit_calls": process_modkit_tsv,
            "modkit_full": process_modkit_tsv,
        }
        processors[args.format](
            infile, outfile, args.offset_5p, args.offset_3p, ref_lengths
        )
    finally:
        outfile.close()
        if infile is not sys.stdin:
            infile.close()


if __name__ == "__main__":
    main()
