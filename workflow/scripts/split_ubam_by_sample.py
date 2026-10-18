#!/usr/bin/env python
"""
Cut every sample's unaligned BAM out of a run-level basecall in ONE pass.

The run BAM is unindexed and unaligned, so selecting a sample's reads means
decompressing the whole file; doing that once per sample (`samtools view -N`)
costs N full scans of a flowcell's basecalls. This reads it once and routes each
record to the output of every sample whose read-id list names it.

Each `--sample NAME=READ_IDS=OUT_BAM` gives a sample's read-id list (one id per
line) and where to write its BAM. Header and record order are preserved, matching
what `samtools view -b -N` produced. A sample that ends up with no reads fails the
run rather than leaving an empty BAM for the aligner to choke on later.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pysam


def load_routes(specs: list[str]) -> tuple[dict[str, list[str]], dict[str, Path]]:
    """read id -> [sample, ...] and sample -> output path."""
    routes: dict[str, list[str]] = {}
    outputs: dict[str, Path] = {}
    for spec in specs:
        try:
            name, ids, out = spec.split("=", 2)
        except ValueError:
            sys.exit(f"--sample expects NAME=READ_IDS=OUT_BAM, got {spec!r}")
        if name in outputs:
            sys.exit(f"sample {name!r} given twice")
        outputs[name] = Path(out)
        with open(ids) as fh:
            for line in fh:
                rid = line.strip()
                if rid:
                    routes.setdefault(rid, []).append(name)
    return routes, outputs


def split(bam: Path, routes, outputs, threads: int) -> dict[str, int]:
    counts = {name: 0 for name in outputs}
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False, threads=threads) as src:
        handles = {}
        try:
            for name, path in outputs.items():
                path.parent.mkdir(parents=True, exist_ok=True)
                handles[name] = pysam.AlignmentFile(
                    str(path), "wb", header=src.header, threads=threads
                )
            for rec in src.fetch(until_eof=True):
                for name in routes.get(rec.query_name, ()):
                    handles[name].write(rec)
                    counts[name] += 1
        finally:
            for h in handles.values():
                h.close()
    return counts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", type=Path, required=True, help="run-level uBAM")
    ap.add_argument(
        "--sample", action="append", required=True, metavar="NAME=READ_IDS=OUT_BAM"
    )
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    routes, outputs = load_routes(args.sample)
    counts = split(args.bam, routes, outputs, args.threads)
    for name, n in counts.items():
        print(f"{name}: {n:,} reads -> {outputs[name]}")
    empty = [n for n, c in counts.items() if c == 0]
    if empty:
        sys.exit(f"no reads in the run BAM for: {', '.join(empty)}")


if __name__ == "__main__":
    main()
