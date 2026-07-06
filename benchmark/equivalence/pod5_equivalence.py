#! /usr/bin/env python
"""
Isolation check: escpod vs the pod5 CLI for merge/filter losslessness.

The migration replaced `pod5 merge` / `pod5 filter` with `escpod merge` /
`escpod filter`. These must be bit-lossless at the level that matters for
basecalling: identical read-id set and identical raw signal per read. This
script runs both tools on the same input and compares their outputs read-for-read.

Uses the `pod5` python package as the reference oracle to read BOTH outputs, so
it must be installed (it ships in the `benchmark` pixi environment). `escpod`
and (for the run modes) the `pod5` CLI must be on PATH.

Modes:
  merge    run `escpod merge` and `pod5 merge` on the same inputs, then compare
      python pod5_equivalence.py merge --workdir /tmp/eq -- in1.pod5 in2.pod5 ...
  filter   run `escpod filter` and `pod5 filter` with the same --ids, then compare
      python pod5_equivalence.py filter --ids ids.txt --workdir /tmp/eq -- in.pod5
  compare  compare two already-produced pod5 files
      python pod5_equivalence.py compare --a escpod.pod5 --b pod5.pod5
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def _load_pod5():
    try:
        import pod5  # noqa: F401

        return pod5
    except ImportError:
        sys.exit(
            "ERROR: the 'pod5' package is required; it is installed into the "
            "default env by `pixi run setup`. Run via `pixi run python ...`."
        )


def _read_signals(path: Path, limit: int | None):
    """read_id (str) -> (signal ndarray, sample_rate). Optionally cap count."""
    pod5 = _load_pod5()
    out = {}
    with pod5.Reader(str(path)) as reader:
        for i, read in enumerate(reader.reads()):
            if limit is not None and i >= limit:
                break
            out[str(read.read_id)] = (read.signal, float(read.run_info.sample_rate))
    return out


def compare(a: Path, b: Path, limit: int | None) -> int:
    import numpy as np

    sa = _read_signals(a, limit)
    sb = _read_signals(b, limit)
    ids_a, ids_b = set(sa), set(sb)

    print("=" * 70)
    print("pod5 equivalence: escpod (a) vs pod5 (b)")
    print("=" * 70)
    print(f"  a={a.name}: {len(sa)} reads    b={b.name}: {len(sb)} reads")

    ok = True
    only_a, only_b = ids_a - ids_b, ids_b - ids_a
    if only_a or only_b:
        ok = False
        print(f"FAIL: read-id set differs (a-only={len(only_a)}, b-only={len(only_b)})")
        for r in list(only_a)[:3]:
            print(f"        a-only example: {r}")
        for r in list(only_b)[:3]:
            print(f"        b-only example: {r}")

    n_sig_mismatch = 0
    n_rate_mismatch = 0
    for rid in ids_a & ids_b:
        siga, ratea = sa[rid]
        sigb, rateb = sb[rid]
        if ratea != rateb:
            n_rate_mismatch += 1
        if siga.shape != sigb.shape or not np.array_equal(siga, sigb):
            n_sig_mismatch += 1
    if n_sig_mismatch:
        ok = False
        print(f"FAIL: {n_sig_mismatch} shared read(s) have differing raw signal")
    if n_rate_mismatch:
        ok = False
        print(f"FAIL: {n_rate_mismatch} shared read(s) have differing sample_rate")

    print("-" * 70)
    print(
        f"  shared reads checked: {len(ids_a & ids_b)}  "
        f"signal-identical: {len(ids_a & ids_b) - n_sig_mismatch}"
    )
    print(
        "RESULT:",
        "PASS — escpod output is lossless vs pod5"
        if ok
        else "FAIL — escpod output differs from pod5",
    )
    return 0 if ok else 1


def run_merge(inputs: list[Path], workdir: Path) -> tuple[Path, Path]:
    workdir.mkdir(parents=True, exist_ok=True)
    esc = workdir / "escpod.merged.pod5"
    p5 = workdir / "pod5.merged.pod5"
    esc.unlink(missing_ok=True)
    p5.unlink(missing_ok=True)
    subprocess.run(["escpod", "merge", "-o", str(esc), *map(str, inputs)], check=True)
    subprocess.run(["pod5", "merge", *map(str, inputs), "-o", str(p5)], check=True)
    return esc, p5


def run_filter(inp: Path, ids: Path, workdir: Path) -> tuple[Path, Path]:
    workdir.mkdir(parents=True, exist_ok=True)
    esc = workdir / "escpod.filtered.pod5"
    p5 = workdir / "pod5.filtered.pod5"
    esc.unlink(missing_ok=True)
    p5.unlink(missing_ok=True)
    subprocess.run(
        ["escpod", "filter", "--ids", str(ids), "-o", str(esc), str(inp)], check=True
    )
    subprocess.run(
        [
            "pod5",
            "filter",
            str(inp),
            "--ids",
            str(ids),
            "--force-overwrite",
            "-o",
            str(p5),
        ],
        check=True,
    )
    return esc, p5


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="mode", required=True)

    m = sub.add_parser("merge")
    m.add_argument("--workdir", type=Path, required=True)
    m.add_argument("inputs", nargs="+", type=Path)

    f = sub.add_parser("filter")
    f.add_argument("--ids", type=Path, required=True)
    f.add_argument("--workdir", type=Path, required=True)
    f.add_argument("inputs", nargs=1, type=Path)

    c = sub.add_parser("compare")
    c.add_argument("--a", type=Path, required=True, help="escpod output")
    c.add_argument("--b", type=Path, required=True, help="pod5 output")

    for p in (m, f, c):
        p.add_argument(
            "--limit", type=int, default=None, help="cap reads compared (default: all)"
        )

    args = ap.parse_args()
    if args.mode == "merge":
        a, b = run_merge(args.inputs, args.workdir)
    elif args.mode == "filter":
        a, b = run_filter(args.inputs[0], args.ids, args.workdir)
    else:
        a, b = args.a, args.b
    sys.exit(compare(a, b, args.limit))


if __name__ == "__main__":
    main()
