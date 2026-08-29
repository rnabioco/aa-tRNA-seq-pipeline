import argparse
import csv
import gzip
import sys
from collections import defaultdict

import pysam
from get_charging_table import charging_tag_value

"""
This script emits base-calling error calls from a BAM file, per read and per
position, considering only alignments on the positive strand.

It is the read-level counterpart of get_bcerror_freqs.py: that script reports
the frequency of mismatches, insertions and deletions at each position, while
this one keeps the reads separate. Read-level information is what a
co-occurrence analysis needs -- for example, testing whether a read being
mis-called at a modification site is associated with that same read being
charged.

Two outputs are available from the same pass over the BAM.

--counts writes the per-site summary that such a test actually consumes: for
each position, how many reads fall in each cell of the error-by-charging 2x2
table. Charging comes from the `cl` tag, thresholded at --ml-threshold.

    ref  pos  err_charged  err_uncharged  match_charged  match_uncharged

The positional output is the per-read table, whose schema matches
`modkit extract calls` so the two can be used interchangeably downstream:

    read_id  ref_position  chrom  within_alignment  call_code

`call_code` is one of:

    X   mismatch
    D   deletion
    I   insertion immediately preceding this position
    -   the read matches the reference at this position

Anything other than "-" counts as an error, mirroring the BCErrorFreq column
of get_bcerror_freqs.py (mismatches + insertions + deletions).

Matching positions are omitted by default, because they are ~70% of rows and
carry nothing individually -- only their count per site matters, and --counts
already records that. Pass --include-matches to emit them, which is only
necessary if a downstream consumer needs to know that a specific read matched
at a specific position rather than that it was not covered there.

Emitting every position of every read produces a very large table even so.
Pass --sites with the positions of interest -- typically those called from
base-calling error rates -- to restrict both outputs.

When --offset-5p and --offset-3p are provided, adapter positions are excluded
and remaining positions are reported in tRNA-only coordinates (1-indexed from
the first tRNA nucleotide), matching get_bcerror_freqs.py.

Example:
    python get_mismatch_calls.py --counts counts.tsv.gz --sites sites.tsv \
        --offset-5p 24 --offset-3p 40 sample.bam reference.fasta calls.tsv.gz
"""

MATCH = "-"
MISMATCH = "X"
DELETION = "D"
INSERTION = "I"


def load_sites(path):
    """
    Load a site list restricting which positions are reported.

    Expects a TSV with a header containing reference and position columns,
    named either `ref`/`pos` or `Reference`/`Position`. Positions are in
    tRNA-only coordinates, matching the output of this script.

    Returns dict of {reference: set of positions}, or None if path is None.
    """
    if path is None:
        return None

    opener = gzip.open if path.endswith(".gz") else open
    sites = {}

    with opener(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Site file {path} is empty")

        ref_col = next(
            (c for c in ("ref", "Reference", "chrom") if c in reader.fieldnames), None
        )
        pos_col = next(
            (c for c in ("pos", "Position", "ref_position") if c in reader.fieldnames),
            None,
        )

        if ref_col is None or pos_col is None:
            raise ValueError(
                f"Site file {path} needs reference and position columns; "
                f"found {reader.fieldnames}"
            )

        for row in reader:
            sites.setdefault(row[ref_col], set()).add(int(row[pos_col]))

    return sites


def call_read_errors(read, faidx, ref, trna_start, trna_end):
    """
    Walk one read's CIGAR and return {reference position: call code}.

    Positions are 0-based reference coordinates and are restricted to
    [trna_start, trna_end). Only positions the alignment actually spans are
    included, so a caller can distinguish "no error" from "not covered".
    """
    calls = {}

    ref_pos = read.reference_start
    read_pos = 0
    read_seq = read.query_sequence

    if read_seq is None:
        return calls

    # Track the last insertion position so that an insertion followed
    # immediately by a mismatch is not counted twice, as in get_bcerror_freqs.py.
    ins_pos = -1

    for cigar_op, cigar_len in read.cigartuples:
        if cigar_op in (0, 7, 8):  # Match or mismatch
            for i in range(cigar_len):
                pos = ref_pos + i
                if not trna_start <= pos < trna_end:
                    continue

                ref_base = faidx.fetch(ref, pos, pos + 1).upper()
                read_base = read_seq[read_pos + i].upper()

                if (read_base != ref_base or cigar_op == 8) and ins_pos != pos:
                    calls[pos] = MISMATCH
                else:
                    calls.setdefault(pos, MATCH)

            ref_pos += cigar_len
            read_pos += cigar_len
        elif cigar_op == 1:  # Insertion, does not consume the reference
            if trna_start <= ref_pos < trna_end:
                calls[ref_pos] = INSERTION
            read_pos += cigar_len
            ins_pos = ref_pos
        elif cigar_op == 2:  # Deletion
            for i in range(cigar_len):
                pos = ref_pos + i
                if trna_start <= pos < trna_end:
                    calls[pos] = DELETION
            ref_pos += cigar_len
        elif cigar_op == 3:  # Reference skip
            ref_pos += cigar_len
        elif cigar_op == 4:  # Soft clip
            read_pos += cigar_len
        # Hard clips and padding (5 and 6) consume neither reference nor read

    return calls


def write_mismatch_calls(
    bam_file,
    fasta_file,
    output_tsv,
    sites=None,
    trim_5p=0,
    trim_3p=0,
    counts_tsv=None,
    charging_tag="cl",
    ml_threshold=200,
    include_matches=False,
):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    faidx = pysam.FastaFile(fasta_file)

    n_rows = 0
    n_no_tag = 0

    # (ref, pos) -> [err_charged, err_uncharged, match_charged, match_uncharged]
    counts = defaultdict(lambda: [0, 0, 0, 0]) if counts_tsv else None

    with gzip.open(output_tsv, "wt", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["read_id", "ref_position", "chrom", "within_alignment", "call_code"]
        )

        for ref in faidx.references:
            ref_sites = None
            if sites is not None:
                ref_sites = sites.get(ref)
                if not ref_sites:
                    continue

            ref_len = faidx.get_reference_length(ref)
            trna_start = trim_5p
            trna_end = ref_len - trim_3p

            for read in samfile.fetch(ref):
                if read.is_unmapped or read.is_reverse:
                    continue

                charged = None
                if counts is not None:
                    try:
                        tag_value = charging_tag_value(read, charging_tag)
                    except ValueError:
                        tag_value = None
                    # A tag of 0 is a confident uncharged call, so test for
                    # None rather than truthiness.
                    if tag_value is None:
                        n_no_tag += 1
                    else:
                        charged = tag_value >= ml_threshold

                calls = call_read_errors(read, faidx, ref, trna_start, trna_end)

                for pos in sorted(calls):
                    # Report in tRNA-only coordinates (1-indexed)
                    trna_pos = pos - trim_5p + 1

                    if ref_sites is not None and trna_pos not in ref_sites:
                        continue

                    code = calls[pos]
                    is_error = code != MATCH

                    if charged is not None:
                        cell = (0 if is_error else 2) + (0 if charged else 1)
                        counts[(ref, trna_pos)][cell] += 1

                    if is_error or include_matches:
                        writer.writerow([read.query_name, trna_pos, ref, "TRUE", code])
                        n_rows += 1

    samfile.close()
    faidx.close()

    if counts_tsv:
        write_counts(counts, counts_tsv)

    if n_no_tag:
        print(
            f"WARNING: {n_no_tag} read(s) carried no usable '{charging_tag}' tag "
            f"and are absent from the counts table",
            file=sys.stderr,
        )

    return n_rows


def write_counts(counts, counts_tsv):
    """
    Write the per-site error-by-charging 2x2 counts.

    This is the summary a charging odds ratio consumes, so a downstream
    analysis needs neither the per-read table nor a coverage table alongside
    it: every cell of every site's contingency table is here.
    """
    with gzip.open(counts_tsv, "wt", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "ref",
                "pos",
                "err_charged",
                "err_uncharged",
                "match_charged",
                "match_uncharged",
            ]
        )

        for ref, pos in sorted(counts):
            writer.writerow([ref, pos] + counts[(ref, pos)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Emit per-read base-calling error calls from a BAM file"
    )

    parser.add_argument("bam_file", help="Path to the BAM file")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("output_tsv", help="Path for the gzipped output TSV")
    parser.add_argument(
        "--sites",
        default=None,
        help="Optional TSV of positions to report, with reference and position "
        "columns (ref/pos or Reference/Position). Positions are in tRNA-only "
        "coordinates. Without this, every aligned position of every read is "
        "reported, which is very large.",
    )
    parser.add_argument(
        "--offset-5p",
        type=int,
        default=0,
        help="Number of bases to skip at 5' end of each reference (adapter + N). "
        "Output positions are renumbered starting at 1 after this offset.",
    )
    parser.add_argument(
        "--offset-3p",
        type=int,
        default=0,
        help="Number of bases to skip at 3' end of each reference (3' adapter).",
    )
    parser.add_argument(
        "--counts",
        default=None,
        help="Also write the per-site error-by-charging 2x2 counts to this path. "
        "This is the summary a charging odds ratio consumes.",
    )
    parser.add_argument(
        "--charging-tag",
        default="cl",
        help="BAM tag holding the charging score (default: cl)",
    )
    parser.add_argument(
        "--ml-threshold",
        type=int,
        default=200,
        help="Charging score at or above which a read is called charged (default: 200)",
    )
    parser.add_argument(
        "--include-matches",
        action="store_true",
        help="Emit rows for positions where the read matches the reference. "
        "These are about 70%% of positions and carry nothing that --counts "
        "does not already record.",
    )
    args = parser.parse_args()

    if args.sites is None:
        print(
            "warning: no --sites given, reporting every aligned position of every read",
            file=sys.stderr,
        )

    n_rows = write_mismatch_calls(
        args.bam_file,
        args.fasta_file,
        args.output_tsv,
        sites=load_sites(args.sites),
        trim_5p=args.offset_5p,
        trim_3p=args.offset_3p,
        counts_tsv=args.counts,
        charging_tag=args.charging_tag,
        ml_threshold=args.ml_threshold,
        include_matches=args.include_matches,
    )

    print(f"wrote {n_rows} calls to {args.output_tsv}", file=sys.stderr)
    if args.counts:
        print(f"wrote per-site counts to {args.counts}", file=sys.stderr)
