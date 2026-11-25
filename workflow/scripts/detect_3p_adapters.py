#!/usr/bin/env python3
"""
Detect 3' adapter identity per read on an unaligned BAM.

Lightweight script that reuses adapter detection logic from add_adapter_tags.py
to classify which 3' adapter each read matches. Outputs a gzipped TSV, one row
per read, for downstream EDX filtering and concordance.

Columns
-------
read_id, adapter_3p
    `adapter_3p` is the call: an adapter name, `none` (no adapter scored above
    --min-score-3p) or `unassigned` (a nearest adapter exists but the gate
    below refused it). `unassigned` is a label of its own on purpose: a read
    the assay could not place is not evidence for the nearest adapter.
score_best, score_second, margin
    Alignment score of the best and second-best adapter by raw score, and
    their gap. With one adapter configured there is no second-best and the
    last two are empty.
adapter_nearest
    The raw best-scoring adapter, before the gate. Equal to `adapter_3p`
    unless the read was gated, which keeps the gate's effect measurable.
hit_start, hit_end, n_aligned, read_len
    Where the best hit sits on the read (0-based, end exclusive) and how many
    adapter bases aligned. The adapter is aligned end to end, so a read cut
    short of the adapter's 3' end shows as n_aligned well under its length and
    hit_end short of read_len, which is how truncation is told from a clean
    miss.
score_frac
    score_best divided by the best adapter's self-score (2 x length). Scores
    scale with adapter length, so the same absolute score means a full-length
    33-nt adapter and a quarter of a 60-nt one; the fraction does not.
reason
    Why a read is `unassigned`: `low_score`, `low_margin`, `stem_only`.

The gate
--------
Off unless asked for (--min-score-frac / --min-margin-frac > 0). A read is
assigned only if its best hit reaches --min-score-frac of that adapter's
self-score, and beats the runner-up by --min-margin-frac of it. Both are
fractions so one setting suits adapters of 33, 45 and 60 nt.

--residual names adapters that are another adapter's stem with the barcode
removed (the plain `edx00`: the shared stem plus a poly-A tail). Every other
adapter, truncated, falls toward it, so it does not compete for the best hit.
It is assigned only when it aligns near full length (--residual-min-frac of its
self-score, which includes the poly-A tail that a barcoded adapter has the
barcode where) AND it beats the best barcoded adapter by the margin. Otherwise
the read is `unassigned`, never `edx00` by default.
"""

import argparse
import gzip
import sys
from collections import Counter

import pysam
from add_adapter_tags import create_scoring_matrix, find_3p_adapter

MATCH_SCORE = 2  # the matrix create_scoring_matrix builds; self-score = 2 x length

COLUMNS = (
    "read_id",
    "adapter_3p",
    "score_best",
    "score_second",
    "margin",
    "adapter_nearest",
    "hit_start",
    "hit_end",
    "n_aligned",
    "read_len",
    "score_frac",
    "reason",
)

UNASSIGNED = "unassigned"


def self_score(adapter_seq):
    """Score of the adapter aligned to itself with no errors."""
    return MATCH_SCORE * len(adapter_seq)


def classify_read(
    seq,
    adapters,
    matrix,
    gap_open=2,
    gap_extend=1,
    min_score=20,
    min_score_frac=0.0,
    min_margin_frac=0.0,
    residual=(),
    residual_min_frac=0.8,
):
    """Call one read's 3' adapter. Returns a dict keyed by COLUMNS[1:].

    `adapters` is a list of (name, sequence). Values are strings/ints ready to
    write; empty string means not applicable.
    """
    out = dict.fromkeys(COLUMNS[1:], "")
    out["adapter_3p"] = "none"
    if not seq:
        return out
    out["read_len"] = len(seq)

    seqs = dict(adapters)
    hits = []  # (score, name, start, end, n_aligned)
    for name, adapter_seq in adapters:
        hit = find_3p_adapter(
            seq, adapter_seq, matrix, gap_open, gap_extend, min_score, detail=True
        )
        if hit:
            hits.append((hit[2], name, hit[0], hit[1], hit[3]))
    if not hits:
        return out

    # Stable sort on score alone: ties keep config order, which is what the
    # original first-wins comparison did.
    hits.sort(key=lambda h: -h[0])
    score, name, start, end, n_aligned = hits[0]
    out.update(
        adapter_nearest=name,
        score_best=score,
        hit_start=start,
        hit_end=end,
        n_aligned=n_aligned,
        score_frac=f"{score / self_score(seqs[name]):.3f}",
    )
    if len(hits) > 1:
        out["score_second"] = hits[1][0]
        out["margin"] = score - hits[1][0]

    gate_on = min_score_frac > 0 or min_margin_frac > 0 or residual
    if not gate_on:
        out["adapter_3p"] = name
        return out

    residual = set(residual)
    barcoded = [h for h in hits if h[1] not in residual]
    stem = [h for h in hits if h[1] in residual]

    # A stem-only adapter is positively called only when it is near full length
    # and clearly beats every barcoded adapter.
    if stem:
        r_score, r_name, r_start, r_end, r_n = stem[0]
        r_self = self_score(seqs[r_name])
        c_score = barcoded[0][0] if barcoded else 0
        if r_score >= residual_min_frac * r_self and r_score - c_score >= max(
            1, round(min_margin_frac * r_self)
        ):
            out.update(
                adapter_3p=r_name,
                score_best=r_score,
                hit_start=r_start,
                hit_end=r_end,
                n_aligned=r_n,
                score_frac=f"{r_score / r_self:.3f}",
            )
            return out

    if not barcoded:
        # Only a stem-only adapter scored, and not convincingly.
        out["adapter_3p"] = UNASSIGNED
        out["reason"] = "stem_only"
        return out

    # Gate the nearest barcoded adapter. Its score and margin are taken within
    # the barcoded set, so a stem-only runner-up does not read as a near-tie.
    c_score, c_name, c_start, c_end, c_n = barcoded[0]
    c_self = self_score(seqs[c_name])
    out.update(
        adapter_nearest=c_name,
        score_best=c_score,
        hit_start=c_start,
        hit_end=c_end,
        n_aligned=c_n,
        score_frac=f"{c_score / c_self:.3f}",
    )
    second = barcoded[1][0] if len(barcoded) > 1 else None
    out["score_second"] = "" if second is None else second
    out["margin"] = "" if second is None else c_score - second

    if c_score < min_score_frac * c_self:
        out["adapter_3p"] = UNASSIGNED
        out["reason"] = "low_score" if not stem else "stem_only"
        return out
    if second is not None and c_score - second < min_margin_frac * c_self:
        out["adapter_3p"] = UNASSIGNED
        out["reason"] = "low_margin"
        return out
    out["adapter_3p"] = c_name
    return out


def parse_adapters(specs):
    adapters = []
    for spec in specs:
        if ":" in spec:
            name, seq = spec.split(":", 1)
            adapters.append((name, seq))
        else:
            adapters.append(("default", spec))
    return adapters


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
    parser.add_argument(
        "--min-score-frac",
        type=float,
        default=0.0,
        help=(
            "Gate: minimum best score as a fraction of the adapter's self-score "
            "(default: 0, off)"
        ),
    )
    parser.add_argument(
        "--min-margin-frac",
        type=float,
        default=0.0,
        help=(
            "Gate: minimum best-vs-second margin as a fraction of the best "
            "adapter's self-score (default: 0, off)"
        ),
    )
    parser.add_argument(
        "--residual",
        action="append",
        default=[],
        help=(
            "Name of a stem-only adapter (repeatable). It does not compete for "
            "the best hit and is assigned only on a near-full-length hit."
        ),
    )
    parser.add_argument(
        "--residual-min-frac",
        type=float,
        default=0.8,
        help="Minimum self-score fraction to assign a --residual adapter (default: 0.8)",
    )
    args = parser.parse_args()

    adapters = parse_adapters(args.adapter_3p_list)
    if not adapters:
        sys.exit("ERROR: at least one --adapter-3p must be specified")
    names = {n for n, _ in adapters}
    unknown = set(args.residual) - names
    if unknown:
        sys.exit(f"ERROR: --residual names no configured adapter: {sorted(unknown)}")

    matrix = create_scoring_matrix()
    counts = Counter()
    total = 0

    with (
        pysam.AlignmentFile(args.bam, "rb", check_sq=False) as bam_in,
        gzip.open(args.output, "wt") as out,
    ):
        out.write("\t".join(COLUMNS) + "\n")
        for read in bam_in.fetch(until_eof=True):
            total += 1
            call = classify_read(
                read.query_sequence,
                adapters,
                matrix,
                min_score=args.min_score_3p,
                min_score_frac=args.min_score_frac,
                min_margin_frac=args.min_margin_frac,
                residual=args.residual,
                residual_min_frac=args.residual_min_frac,
            )
            out.write(
                read.query_name
                + "\t"
                + "\t".join(str(call[c]) for c in COLUMNS[1:])
                + "\n"
            )
            counts[call["adapter_3p"]] += 1

    # Summary to stderr
    print(f"total_reads\t{total}", file=sys.stderr)
    for name in sorted(counts.keys()):
        pct = 100.0 * counts[name] / total if total > 0 else 0.0
        print(f"{name}\t{counts[name]}\t({pct:.1f}%)", file=sys.stderr)


if __name__ == "__main__":
    main()
