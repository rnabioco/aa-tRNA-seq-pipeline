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
    leech_bam: Path, remora_bam: Path, call_agreement_min: float, cl_abs_diff_max: int
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
    call_a = va >= CL_THRESHOLD
    call_b = vb >= CL_THRESHOLD
    agree = float((call_a == call_b).mean())
    absdiff = np.abs(va - vb)

    print("=" * 70)
    print("classifier equivalence: leech vs remora")
    print("=" * 70)
    print(
        f"  reads: leech={len(a)} remora={len(b)} shared={len(shared)} "
        f"(leech-only={only_leech}, remora-only={only_remora})"
    )
    print(f"  charged-call agreement (thr {CL_THRESHOLD}): {agree:.5f}")
    print(
        f"  |Δcl|: mean={absdiff.mean():.3f} p99={np.percentile(absdiff, 99):.3f} "
        f"max={absdiff.max():.3f}"
    )
    print("-" * 70)

    ok = True
    if agree < call_agreement_min:
        print(f"FAIL: call agreement {agree:.5f} < {call_agreement_min}")
        ok = False
    if absdiff.max() > cl_abs_diff_max:
        print(f"FAIL: max |Δcl| {absdiff.max():.1f} > {cl_abs_diff_max}")
        ok = False
    read_overlap = len(shared) / len(set(a) | set(b))
    if read_overlap < 0.99:
        print(f"FAIL: read-id overlap {read_overlap:.4f} < 0.99")
        ok = False
    print(
        "RESULT:", "PASS — engines are equivalent" if ok else "FAIL — engines diverged"
    )
    return 0 if ok else 1


def run_both(pod5: Path, bam: Path, model: Path, workdir: Path) -> tuple[Path, Path]:
    workdir.mkdir(parents=True, exist_ok=True)
    leech_bam = workdir / "leech.charging.bam"
    remora_bam = workdir / "remora.charging.bam"

    # invocations mirror rules classify_charging_leech / classify_charging.
    # No --motif/--motif-offset: both engines read the motif ('CCAGGC', offset 3)
    # from the model, so this compares them on identical anchoring — the whole
    # point of the check. Passing an explicit offset that disagrees with the
    # model makes leech refuse to run.
    subprocess.run(
        [
            "leech",
            "predict",
            "--model",
            str(model),
            "--pod5",
            str(pod5),
            "--bam",
            str(bam),
            "--output",
            str(leech_bam),
            "--device",
            "cuda",
            "--reference-anchored",
            "--workers",
            "4",
            "--batch-size",
            "512",
        ],
        check=True,
    )

    subprocess.run(
        [
            "remora",
            "infer",
            "from_pod5_and_bam",
            str(pod5),
            str(bam),
            "--model",
            str(model),
            "--out-bam",
            str(remora_bam),
            "--reference-anchored",
        ],
        check=True,
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
    r.add_argument("--workdir", type=Path, required=True)

    for p in (c, r):
        p.add_argument("--call-agreement-min", type=float, default=0.995)
        p.add_argument("--cl-abs-diff-max", type=int, default=5)

    args = ap.parse_args()

    if args.mode == "run":
        leech_bam, remora_bam = run_both(args.pod5, args.bam, args.model, args.workdir)
    else:
        leech_bam, remora_bam = args.leech, args.remora

    sys.exit(diff(leech_bam, remora_bam, args.call_agreement_min, args.cl_abs_diff_max))


if __name__ == "__main__":
    main()
