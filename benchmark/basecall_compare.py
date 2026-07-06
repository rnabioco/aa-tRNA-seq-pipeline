#! /usr/bin/env python
"""
Compare dorado basecalls across versions.

Given the unmapped BAMs produced by basecalling the SAME pod5 with different
dorado versions, quantify how much the basecalls changed: read-id overlap,
per-read sequence agreement (exact match + mean similarity), and shifts in read
length and mean quality.

Basecalls are EXPECTED to differ across versions (that is what this measures) —
this is a report, not a pass/fail gate. The first --bam given is the reference
that the others are compared against.

Usage:
    python benchmark/basecall_compare.py \
        --bam 1.4.0=s.1.4.0.ubam --bam 2.0.1=s.2.0.1.ubam [--json out.json]
"""

from __future__ import annotations

import argparse
import difflib
import json
from pathlib import Path

import pysam

try:
    import parasail

    _DNA_MATRIX = parasail.dnafull
except Exception:  # parasail optional; fall back to difflib
    parasail = None


def _identity(a: str, b: str) -> float:
    """Global-alignment identity (matches / alignment length) of two reads.

    Uses parasail (Needleman-Wunsch) when available; otherwise difflib with
    autojunk disabled. autojunk MUST be off: on 4-letter sequences >200 nt it
    flags every base as junk and collapses the ratio to a meaningless value.
    """
    if not a or not b:
        return 0.0
    if parasail is not None:
        res = parasail.nw_stats_striped_16(a, b, 10, 1, _DNA_MATRIX)
        return res.matches / res.length if res.length else 0.0
    return difflib.SequenceMatcher(None, a, b, autojunk=False).ratio()


def _read_bam(path: Path) -> dict[str, tuple[str, float, int]]:
    """read_id -> (sequence, mean_quality, length) for primary records."""
    out: dict[str, tuple[str, float, int]] = {}
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary:
                continue
            seq = r.query_sequence or ""
            quals = r.query_qualities
            mq = float(sum(quals) / len(quals)) if quals else 0.0
            out[r.query_name] = (seq, mq, len(seq))
    return out


def _pair_stats(ref: dict, other: dict, sim_cap: int) -> dict:
    shared = sorted(set(ref) & set(other))
    n = len(shared)
    stats = {
        "n_ref": len(ref),
        "n_other": len(other),
        "n_shared": n,
        "read_overlap": n / len(set(ref) | set(other)) if (ref or other) else 0.0,
    }
    if not n:
        return stats

    exact = 0
    idents = []
    dlen = dq = 0.0
    for i, rid in enumerate(shared):
        s_ref, q_ref, l_ref = ref[rid]
        s_oth, q_oth, l_oth = other[rid]
        if s_ref == s_oth:
            exact += 1
            if i < sim_cap:
                idents.append(1.0)
        elif i < sim_cap:
            idents.append(_identity(s_ref, s_oth))
        dlen += l_oth - l_ref
        dq += q_oth - q_ref

    stats["exact_seq_match_frac"] = exact / n
    stats["mean_identity"] = (sum(idents) / len(idents)) if idents else None
    stats["identity_reads_sampled"] = min(n, sim_cap)
    stats["mean_len_delta"] = dlen / n
    stats["mean_qual_delta"] = dq / n
    stats["mean_qual_ref"] = sum(v[1] for v in ref.values()) / len(ref)
    stats["mean_qual_other"] = sum(v[1] for v in other.values()) / len(other)
    return stats


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--bam",
        action="append",
        required=True,
        metavar="LABEL=PATH",
        help="version-labeled uBAM; repeat (>=2). First is the reference.",
    )
    ap.add_argument(
        "--sim-cap",
        type=int,
        default=2000,
        help="max shared reads to compute similarity on (default 2000)",
    )
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()

    if len(args.bam) < 2:
        ap.error("need at least two --bam LABEL=PATH entries")

    labeled = []
    for spec in args.bam:
        if "=" not in spec:
            ap.error(f"--bam expects LABEL=PATH, got {spec!r}")
        label, path = spec.split("=", 1)
        labeled.append((label, Path(path)))

    reads = {label: _read_bam(path) for label, path in labeled}
    ref_label = labeled[0][0]
    ref = reads[ref_label]

    print("=" * 74)
    print(f"dorado basecall comparison (reference = {ref_label})")
    print("=" * 74)
    print("read counts:")
    for label, _ in labeled:
        print(f"  {label:<12} {len(reads[label])} reads")
    print()

    report = {
        "reference": ref_label,
        "counts": {k: len(v) for k, v in reads.items()},
        "pairs": {},
    }
    for label, _ in labeled[1:]:
        s = _pair_stats(ref, reads[label], args.sim_cap)
        report["pairs"][label] = s
        print(f"[{ref_label} vs {label}]")
        print(
            f"  read overlap        : {s['read_overlap']:.4f} ({s['n_shared']} shared)"
        )
        if s["n_shared"]:
            print(f"  exact seq match     : {s['exact_seq_match_frac']:.4f}")
            mi = s["mean_identity"]
            print(
                f"  mean identity       : {mi:.4f} (n={s['identity_reads_sampled']})"
                if mi is not None
                else "  mean identity       : n/a"
            )
            print(f"  mean length delta   : {s['mean_len_delta']:+.2f} bp")
            print(
                f"  mean qual delta     : {s['mean_qual_delta']:+.3f} "
                f"({s['mean_qual_ref']:.2f} -> {s['mean_qual_other']:.2f})"
            )
        print()

    print("-" * 74)
    print(
        "note: cross-version basecall differences are expected; this quantifies them."
    )

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        with open(args.json, "w") as fh:
            json.dump(report, fh, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
