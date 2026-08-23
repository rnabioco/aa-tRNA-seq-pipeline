#!/usr/bin/env python3
"""
Detect 3' adapter identity per read on an unaligned BAM.

Lightweight script that reuses adapter detection logic from add_adapter_tags.py
to classify which 3' adapter each read matches. Outputs a gzipped TSV of
read_id, adapter_3p, score_best, score_second and margin for downstream EDX
filtering and concordance.

Reads with no adapter match get adapter_3p = "none".

`margin` is the alignment-score gap between the best and second-best adapter.
It exists because the best hit alone cannot distinguish a confident call from a
near-tie between two similar adapters, and the two mean opposite things when the
output is used to measure adapter crosstalk: a near-tie is a read the assay
could not assign, not evidence that one adapter ended up in another's library.
The upstream barcode training corpus gated on margin > 5 for the same reason.
With a single adapter configured there is no second-best, and margin is left
empty rather than being reported as an unbounded gap.
"""

import argparse
import gzip
import sys
from collections import Counter

import pysam
from add_adapter_tags import create_scoring_matrix, find_3p_adapter


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--bam", required=True, help="Input unaligned BAM file")
    parser.add_argument("--output", required=True, help="Output gzipped TSV file")
    parser.add_argument(
        "--adapter-3p",
        action="append",
        default=[],
        dest="adapter_3p_list",
        help="3' adapter as 'name:sequence' (can specify multiple)",
    )
    parser.add_argument(
        "--min-score-3p",
        type=int,
        default=20,
        help="Minimum alignment score for 3' adapter (default: 20)",
    )
    args = parser.parse_args()

    # Parse adapter list
    adapters = []
    for spec in args.adapter_3p_list:
        if ":" in spec:
            name, seq = spec.split(":", 1)
            adapters.append((name, seq))
        else:
            adapters.append(("default", spec))

    if not adapters:
        sys.exit("ERROR: at least one --adapter-3p must be specified")

    matrix = create_scoring_matrix()
    gap_open = 2
    gap_extend = 1
    min_score = args.min_score_3p

    counts = Counter()
    total = 0

    with (
        pysam.AlignmentFile(args.bam, "rb", check_sq=False) as bam_in,
        gzip.open(args.output, "wt") as out,
    ):
        out.write("read_id\tadapter_3p\tscore_best\tscore_second\tmargin\n")
        for read in bam_in.fetch(until_eof=True):
            total += 1
            seq = read.query_sequence
            if seq is None:
                out.write(f"{read.query_name}\tnone\t\t\t\n")
                counts["none"] += 1
                continue

            # Score every adapter rather than short-circuiting on the best,
            # so the runner-up is available for the margin.
            scores = []
            for name, adapter_seq in adapters:
                hit = find_3p_adapter(
                    seq, adapter_seq, matrix, gap_open, gap_extend, min_score
                )
                if hit:
                    scores.append((hit[2], name))  # (start, end, score)
            # Stable sort on score alone: ties keep config order, which is
            # what the previous first-wins comparison did. Sorting on the
            # tuple would break ties by adapter name instead.
            scores.sort(key=lambda s: -s[0])

            if not scores:
                out.write(f"{read.query_name}\tnone\t\t\t\n")
                counts["none"] += 1
                continue

            best_score, adapter_name = scores[0]
            second = f"{scores[1][0]}" if len(scores) > 1 else ""
            margin = f"{best_score - scores[1][0]}" if len(scores) > 1 else ""
            out.write(
                f"{read.query_name}\t{adapter_name}\t{best_score}\t{second}\t{margin}\n"
            )
            counts[adapter_name] += 1

    # Summary to stderr
    print(f"total_reads\t{total}", file=sys.stderr)
    for name in sorted(counts.keys()):
        pct = 100.0 * counts[name] / total if total > 0 else 0.0
        print(f"{name}\t{counts[name]}\t({pct:.1f}%)", file=sys.stderr)


if __name__ == "__main__":
    main()
