#!/usr/bin/env python3
"""
Detect 3' adapter identity per read on an unaligned BAM.

Lightweight script that reuses adapter detection logic from add_adapter_tags.py
to classify which 3' adapter each read matches. Outputs a gzipped two-column TSV
(read_id, adapter_3p) for downstream EDX filtering.

Reads with no adapter match get adapter_3p = "none".
"""

import argparse
import gzip
import sys
from collections import Counter

from add_adapter_tags import create_scoring_matrix, find_best_3p_adapter

import pysam


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

    with pysam.AlignmentFile(args.bam, "rb", check_sq=False) as bam_in:
        with gzip.open(args.output, "wt") as out:
            out.write("read_id\tadapter_3p\n")
            for read in bam_in.fetch(until_eof=True):
                total += 1
                seq = read.query_sequence
                if seq is None:
                    out.write(f"{read.query_name}\tnone\n")
                    counts["none"] += 1
                    continue

                result = find_best_3p_adapter(
                    seq, adapters, matrix, gap_open, gap_extend, min_score
                )
                if result:
                    adapter_name = result[3]  # (start, end, score, name)
                else:
                    adapter_name = "none"

                out.write(f"{read.query_name}\t{adapter_name}\n")
                counts[adapter_name] += 1

    # Summary to stderr
    print(f"total_reads\t{total}", file=sys.stderr)
    for name in sorted(counts.keys()):
        pct = 100.0 * counts[name] / total if total > 0 else 0.0
        print(f"{name}\t{counts[name]}\t({pct:.1f}%)", file=sys.stderr)


if __name__ == "__main__":
    main()
