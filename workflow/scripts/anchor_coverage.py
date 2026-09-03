#!/usr/bin/env python
"""Does each aligned read actually span the CCA anchor the charging model reads?

`classify_charging` anchors on the CCA-adapter junction, which only exists in
reference coordinates, so it can only call a read whose alignment reaches it. A
read that aligns but stops short is uncallable for a structural reason — it is
not that the model looked and declined.

That distinction used to be invisible: under Remora the two were indis-
tinguishable in the outputs, and the only trace of either was that
`align_stats`'s `classified` row came out smaller than its `aligned` row.
`escpod classify --tsv` now names the reads it saw and did not score, so
this script covers the other half — the reads that never reached the model at
all — and the two are read together in `read_attrition.tsv.gz`.

This runs against the ALIGNED bam, which is `temp()` under
`cleanup_intermediates`, so it has to happen inside the pipeline: after a run
finishes the evidence is gone. Reconstructing it for the 2026-08-06 LDX run meant
re-basecalling from POD5.

The anchor is located per reference rather than assumed at a fixed offset:
references differ in tRNA length, and the 3' adapter differs per EDX barcode, so
a hardcoded coordinate would be silently wrong for most of them. It is the three
bases immediately preceding the 3' adapter's constant prefix.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pysam


def locate_anchors(fasta: Path, adapter_prefix: str) -> dict[str, tuple[int, int]]:
    """0-based half-open span of the CCA preceding the 3' adapter, per reference."""
    fa = pysam.FastaFile(str(fasta))
    spans: dict[str, tuple[int, int]] = {}
    for name in fa.references:
        seq = fa.fetch(name).upper()
        i = seq.rfind(adapter_prefix.upper())
        if i >= 3:
            spans[name] = (i - 3, i)
    return spans


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", type=Path, required=True, help="aligned BAM")
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument(
        "--adapter-prefix",
        default="GGCTTCTTCTTGCTCTTAGGAAGGC",
        help="constant 5' part of the 3' adapter; the CCA sits immediately before it",
    )
    args = ap.parse_args()

    spans = locate_anchors(args.reference, args.adapter_prefix)
    counts = {
        "scored": 0,
        "covers_anchor": 0,
        "ends_before_anchor": 0,
        "starts_after_anchor": 0,
        "reverse_strand": 0,
        "no_anchor_in_reference": 0,
    }

    with pysam.AlignmentFile(str(args.bam), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            span = spans.get(rec.reference_name)
            if span is None:
                counts["no_anchor_in_reference"] += 1
                continue
            counts["scored"] += 1
            if rec.is_reverse:
                counts["reverse_strand"] += 1
            if rec.reference_start <= span[0] and rec.reference_end >= span[1]:
                counts["covers_anchor"] += 1
            elif rec.reference_end < span[1]:
                counts["ends_before_anchor"] += 1
            elif rec.reference_start > span[0]:
                counts["starts_after_anchor"] += 1

    scored = counts["scored"] or 1
    uncallable = counts["ends_before_anchor"] + counts["starts_after_anchor"]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output, "wt") as out:
        out.write(
            "sample\tscored\tcovers_anchor\tends_before_anchor\tstarts_after_anchor\t"
            "reverse_strand\tno_anchor_in_reference\tpct_covers\tpct_uncallable\n"
        )
        out.write(
            f"{args.sample}\t{counts['scored']}\t{counts['covers_anchor']}\t"
            f"{counts['ends_before_anchor']}\t{counts['starts_after_anchor']}\t"
            f"{counts['reverse_strand']}\t{counts['no_anchor_in_reference']}\t"
            f"{100 * counts['covers_anchor'] / scored:.4f}\t"
            f"{100 * uncallable / scored:.4f}\n"
        )

    print(
        f"{args.sample}: {counts['scored']:,} scored, "
        f"{100 * counts['covers_anchor'] / scored:.2f}% cover the CCA anchor, "
        f"{100 * uncallable / scored:.2f}% cannot be anchored"
    )


if __name__ == "__main__":
    main()
