#! /usr/bin/env python
"""
Snapshot the result-bearing outputs of a pipeline run into a comparable
"fingerprint".

A fingerprint is a directory:

    <fingerprint>/
      fingerprint.json          aggregate metrics per sample + tool manifest
      perread/<sample>.parquet  read_id, tRNA, cl        (for per-read joins)
      perpos/<sample>.bcerror.parquet
      perpos/<sample>.pileup.parquet

Aggregates live in JSON so a fingerprint is human-diffable at a glance; the
big per-read / per-position tables are Parquet sidecars so `compare.py` can do
exact joins (per-read charging-call agreement, per-position error/mod deltas)
without re-reading the whole pipeline output.

Usage:
    python benchmark/fingerprint.py OUTPUT_DIR FINGERPRINT_DIR [--label NAME] [--git-ref REF]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import lib


CL_CHARGED_THRESHOLD = 200  # matches get_trna_charging_cpm.py default


def _charging_summary(df) -> dict:
    cl = df["cl"].to_numpy(dtype=float)
    n = int(len(cl))
    charged = int((cl >= CL_CHARGED_THRESHOLD).sum())
    # coarse histogram so aggregate drift is visible without the sidecar
    hist, _ = np.histogram(cl, bins=16, range=(0, 256))
    return {
        "n_classified": n,
        "charged_fraction": (charged / n) if n else None,
        "cl_mean": float(cl.mean()) if n else None,
        "cl_hist_0_256_16bins": [int(x) for x in hist],
    }


def _cpm_inline(cpm_df) -> dict:
    """Per-tRNA counts/CPM, kept inline (hundreds of rows)."""
    out = {}
    for trna, row in cpm_df.iterrows():
        out[str(trna)] = {k: float(row[k]) for k in cpm_df.columns if k in row}
    return out


def fingerprint(outdir: Path, dest: Path, label: str, git_ref: str | None) -> dict:
    samples = lib.discover_samples(outdir)
    fp: dict = {
        "schema_version": 1,
        "label": label,
        "git_ref": git_ref,
        "output_directory": str(outdir),
        "charged_threshold": CL_CHARGED_THRESHOLD,
        "manifest_tools": (lib.read_manifest(outdir) or {}).get("tools"),
        "samples": {},
    }

    (dest / "perread").mkdir(parents=True, exist_ok=True)
    (dest / "perpos").mkdir(parents=True, exist_ok=True)

    for sample in samples:
        s: dict = {}

        align = lib.read_align_stats(outdir, sample)
        if align:
            s["align_stats"] = align

        prob = lib.read_charging_prob(outdir, sample)
        if prob is not None and not prob.empty:
            s["charging"] = _charging_summary(prob)
            prob.to_parquet(dest / "perread" / f"{sample}.parquet", index=False)

        cpm = lib.read_charging_cpm(outdir, sample)
        if cpm is not None and not cpm.empty:
            s["cpm"] = {
                "n_trna": int(len(cpm)),
                "total_counts": float(
                    cpm.get("counts_charged", 0).sum()
                    + cpm.get("counts_uncharged", 0).sum()
                ),
                "per_trna": _cpm_inline(cpm),
            }

        bcerr = lib.read_bcerror(outdir, sample)
        if bcerr is not None and not bcerr.empty:
            s["bcerror"] = {
                "n_positions": int(len(bcerr)),
                "mean_bcerror_freq": float(bcerr["BCErrorFreq"].mean())
                if "BCErrorFreq" in bcerr
                else None,
            }
            bcerr.to_parquet(dest / "perpos" / f"{sample}.bcerror.parquet", index=False)

        pileup = lib.read_pileup(outdir, sample)
        if pileup is not None and not pileup.empty:
            s["modkit"] = {
                "n_sites": int(len(pileup)),
                "mean_frac_mod": float(pileup["frac_mod"].mean()),
            }
            pileup.to_parquet(dest / "perpos" / f"{sample}.pileup.parquet", index=False)

        fp["samples"][sample] = s

    lib.write_json(dest / "fingerprint.json", fp)
    return fp


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("output_dir", type=Path, help="pipeline output_directory")
    ap.add_argument("fingerprint_dir", type=Path, help="destination fingerprint dir")
    ap.add_argument("--label", default="run", help="name for this run in reports")
    ap.add_argument("--git-ref", default=None, help="git ref the run was built from")
    args = ap.parse_args()

    fp = fingerprint(args.output_dir, args.fingerprint_dir, args.label, args.git_ref)
    n = len(fp["samples"])
    print(
        f"wrote fingerprint '{args.label}' for {n} sample(s) -> {args.fingerprint_dir}"
    )


if __name__ == "__main__":
    main()
