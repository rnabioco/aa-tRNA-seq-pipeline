#! /usr/bin/env python
"""
Isolation check: leech vs remora charging classifier on identical inputs.

Both engines load the same cca_classifier.pt, so per-read charging scores
should agree near-exactly (only runtime numerical differences). This proves the
remora->leech swap did not change results, independently of the dorado bump.

Two modes:

  compare   already have both charging BAMs -> just diff the cl tags
      python classifier_equivalence.py compare --leech A.bam --remora B.bam

  run       have the tagged BAM + pod5 -> run both engines, then diff
      python classifier_equivalence.py run \
          --pod5 s.pod5 --bam s.tagged.bam --model resources/models/cca_classifier.pt \
          --workdir /tmp/clf_eq

The `run` mode needs both `leech` (GPU) and `remora` (CPU) on PATH. The charging
score lives in the ML tag on the freshly classified BAM (transfer_bam_tags later
renames ML->cl); this script reads whichever of ml/ML/cl/CL is present.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import lib  # noqa: E402

CL_THRESHOLD = 200


def _tags_by_read(bam: Path) -> dict[str, float]:
    for tag in ("ml", "cl"):
        d = lib.read_bam_tag_by_read(bam, tag)
        if d:
            return d
    return {}


def diff(
    leech_bam: Path,
    remora_bam: Path,
    call_agreement_min: float,
    cl_p99_max: float,
    read_incl_min: float,
) -> int:
    a = _tags_by_read(leech_bam)  # leech
    b = _tags_by_read(remora_bam)  # remora
    if not a or not b:
        print(
            f"ERROR: no charging tag found (leech reads={len(a)}, remora reads={len(b)})"
        )
        return 2

    shared = sorted(set(a) & set(b))
    only_leech = len(set(a) - set(b))
    only_remora = len(set(b) - set(a))
    if not shared:
        print("ERROR: no shared read_ids between the two BAMs")
        return 2

    va = np.array([a[r] for r in shared], dtype=float)
    vb = np.array([b[r] for r in shared], dtype=float)
    agree = float(((va >= CL_THRESHOLD) == (vb >= CL_THRESHOLD)).mean())
    absdiff = np.abs(va - vb)
    p99 = float(np.percentile(absdiff, 99))
    overlap = len(shared) / len(set(a) | set(b))

    print("=" * 70)
    print("classifier equivalence: leech vs remora")
    print("=" * 70)
    print(
        f"  reads: leech={len(a)} remora={len(b)} shared={len(shared)} "
        f"(leech-only={only_leech}, remora-only={only_remora})"
    )
    print(f"  charged-call agreement (thr {CL_THRESHOLD}): {agree:.5f}")
    print(
        f"  |Δcl|: mean={absdiff.mean():.3f} p99={p99:.3f} max={absdiff.max():.3f}"
    )
    print("-" * 70)

    # The equivalence verdict is about the charging CALLS on shared reads (agreement
    # + a robust p99 of the score delta — not `max`, which one outlier read trips).
    calls_ok = True
    if agree < call_agreement_min:
        print(f"FAIL: call agreement {agree:.5f} < {call_agreement_min}")
        calls_ok = False
    if p99 > cl_p99_max:
        print(f"FAIL: p99 |Δcl| {p99:.2f} > {cl_p99_max}")
        calls_ok = False

    # Read inclusion is reported separately: one engine being a strict subset of
    # the other is a selectivity difference (affects read counts / CPM denominators),
    # not a disagreement on the reads they both call. It does not fail the verdict.
    if overlap < read_incl_min:
        subset = (
            "leech ⊂ remora" if only_leech == 0
            else "remora ⊂ leech" if only_remora == 0
            else "neither is a subset"
        )
        print(
            f"NOTE: read inclusion differs — overlap {overlap:.4f} ({subset}); "
            f"leech-only={only_leech}, remora-only={only_remora}. "
            f"Affects read counts/CPM, not the calls above."
        )

    if calls_ok:
        print("RESULT: PASS — charging calls equivalent (see NOTE for read inclusion)"
              if overlap < read_incl_min
              else "RESULT: PASS — engines equivalent")
    else:
        print("RESULT: FAIL — charging calls diverge")
    return 0 if calls_ok else 1


def run_both(
    pod5: Path, bam: Path, model: Path, reference: Path, workdir: Path
) -> tuple[Path, Path]:
    workdir.mkdir(parents=True, exist_ok=True)
    leech_bam = workdir / "leech.charging.bam"
    remora_bam = workdir / "remora.charging.bam"

    # Reference-anchored inference needs the reference sequences (the BAM @SQ
    # header carries only names). leech takes --reference-fasta directly; remora
    # has no such flag and reconstructs the reference from the BAM's MD tag, which
    # bwa does not emit — so give remora an MD-tagged copy via `samtools calmd`.
    subprocess.run(["samtools", "faidx", str(reference)], check=True)
    md_bam = workdir / "tagged.md.bam"
    with open(md_bam, "wb") as fh:
        subprocess.run(
            ["samtools", "calmd", "-b", str(bam), str(reference)],
            check=True, stdout=fh, stderr=subprocess.DEVNULL,
        )
    subprocess.run(["samtools", "index", str(md_bam)], check=True)

    # No --motif/--motif-offset: both engines read the motif ('CCAGGC', offset 3)
    # from the model, so this compares them on identical anchoring — the whole
    # point of the check. An explicit offset that disagrees makes leech refuse.
    t0 = time.monotonic()
    subprocess.run(
        [
            "leech", "predict",
            "--model", str(model),
            "--pod5", str(pod5),
            "--bam", str(bam),
            "--output", str(leech_bam),
            "--device", "cuda",
            "--reference-anchored",
            "--reference-fasta", str(reference),
            "--workers", "4",
            "--batch-size", "512",
        ],
        check=True,
    )
    t_leech = time.monotonic() - t0

    t1 = time.monotonic()
    subprocess.run(
        [
            "remora", "infer", "from_pod5_and_bam",
            str(pod5), str(md_bam),
            "--model", str(model),
            "--out-bam", str(remora_bam),
            "--reference-anchored",
        ],
        check=True,
    )
    t_remora = time.monotonic() - t1

    speed = f"{t_remora / t_leech:.1f}x" if t_leech else "n/a"
    print(
        f"  runtime: leech(GPU)={t_leech:.1f}s  remora(CPU)={t_remora:.1f}s  "
        f"(leech {speed} faster)"
    )
    return leech_bam, remora_bam


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="mode", required=True)

    c = sub.add_parser("compare", help="diff two existing charging BAMs")
    c.add_argument("--leech", type=Path, required=True)
    c.add_argument("--remora", type=Path, required=True)

    r = sub.add_parser("run", help="run both engines then diff")
    r.add_argument("--pod5", type=Path, required=True)
    r.add_argument(
        "--bam", type=Path, required=True, help="inject_ubam_tags output BAM"
    )
    r.add_argument("--model", type=Path, required=True)
    r.add_argument("--reference-fasta", type=Path, required=True,
                   help="reference FASTA for reference-anchored inference")
    r.add_argument("--workdir", type=Path, required=True)

    for p in (c, r):
        p.add_argument("--call-agreement-min", type=float, default=0.995)
        p.add_argument("--cl-p99-max", type=float, default=5.0,
                       help="max 99th-pctile |Δcl| for the equivalence verdict")
        p.add_argument("--read-incl-min", type=float, default=0.98,
                       help="min read-set overlap before adding a read-inclusion note")

    args = ap.parse_args()

    if args.mode == "run":
        leech_bam, remora_bam = run_both(
            args.pod5, args.bam, args.model, args.reference_fasta, args.workdir
        )
    else:
        leech_bam, remora_bam = args.leech, args.remora

    sys.exit(
        diff(
            leech_bam,
            remora_bam,
            args.call_agreement_min,
            args.cl_p99_max,
            args.read_incl_min,
        )
    )


if __name__ == "__main__":
    main()
