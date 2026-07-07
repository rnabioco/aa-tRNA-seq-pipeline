#! /usr/bin/env python
"""
Compare two fingerprints (see fingerprint.py) under a named tolerance profile
and report whether the migration preserved results.

    python benchmark/compare.py BASELINE_FP CANDIDATE_FP \
        [--profile strict|aggregate] [--tolerances benchmark/tolerances.yml] \
        [--json report.json] [--report-only]

Exit status is 0 when every checked metric is within tolerance, 1 otherwise
(unless --report-only). Metrics absent from either fingerprint are SKIPPED, not
failed, so a partial run still yields a useful partial report.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

import lib

CL_THRESHOLD = 200


@dataclass
class Check:
    sample: str
    metric: str
    status: str  # PASS | FAIL | SKIP
    detail: str
    observed: float | None = None
    limit: float | None = None


@dataclass
class Report:
    checks: list[Check] = field(default_factory=list)

    def add(self, *a, **k):
        self.checks.append(Check(*a, **k))

    @property
    def failed(self) -> list[Check]:
        return [c for c in self.checks if c.status == "FAIL"]

    @property
    def passed(self) -> list[Check]:
        return [c for c in self.checks if c.status == "PASS"]

    @property
    def skipped(self) -> list[Check]:
        return [c for c in self.checks if c.status == "SKIP"]


def _perread(fp_dir: Path, sample: str) -> pd.DataFrame | None:
    p = fp_dir / "perread" / f"{sample}.parquet"
    return pd.read_parquet(p) if p.exists() else None


def _perpos(fp_dir: Path, sample: str, kind: str) -> pd.DataFrame | None:
    p = fp_dir / "perpos" / f"{sample}.{kind}.parquet"
    return pd.read_parquet(p) if p.exists() else None


# --- per-metric comparisons ----------------------------------------------------


def check_align(rep, sample, a, b, tol):
    if not a or not b:
        rep.add(sample, "align_stats", "SKIP", "align_stats missing on one side")
        return
    if (
        "pct_mapped" in a
        and "pct_mapped" in b
        and tol.get("pct_mapped_abs_tol") is not None
    ):
        d = abs(a["pct_mapped"] - b["pct_mapped"])
        lim = tol["pct_mapped_abs_tol"]
        rep.add(
            sample,
            "pct_mapped",
            "PASS" if d <= lim else "FAIL",
            f"|Δ|={d:.3f} pp (base={a['pct_mapped']:.2f}, cand={b['pct_mapped']:.2f})",
            d,
            lim,
        )
    if "n_reads" in a and "n_reads" in b and tol.get("n_reads_rel_tol") is not None:
        base = a["n_reads"] or 1
        d = abs(a["n_reads"] - b["n_reads"]) / base
        lim = tol["n_reads_rel_tol"]
        rep.add(
            sample,
            "n_reads",
            "PASS" if d <= lim else "FAIL",
            f"rel Δ={d:.4f} (base={a['n_reads']:.0f}, cand={b['n_reads']:.0f})",
            d,
            lim,
        )


def check_charging(rep, sample, base_dir, cand_dir, s_a, s_b, tol):
    da, db = _perread(base_dir, sample), _perread(cand_dir, sample)
    if da is None or db is None:
        rep.add(
            sample, "charging", "SKIP", "per-read charging sidecar missing on one side"
        )
    else:
        merged = da.merge(db, on="read_id", suffixes=("_a", "_b"))
        n_union = len(set(da["read_id"]) | set(db["read_id"]))
        frac = len(merged) / n_union if n_union else 0.0
        if tol.get("read_join_min_fraction") is not None:
            lim = tol["read_join_min_fraction"]
            rep.add(
                sample,
                "read_join_fraction",
                "PASS" if frac >= lim else "FAIL",
                f"shared={len(merged)}/{n_union} ({frac:.4f})",
                frac,
                lim,
            )
        if len(merged):
            call_a = merged["cl_a"] >= CL_THRESHOLD
            call_b = merged["cl_b"] >= CL_THRESHOLD
            agree = float((call_a == call_b).mean())
            if tol.get("call_agreement_min") is not None:
                lim = tol["call_agreement_min"]
                rep.add(
                    sample,
                    "charge_call_agreement",
                    "PASS" if agree >= lim else "FAIL",
                    f"{agree:.4f} over {len(merged)} shared reads",
                    agree,
                    lim,
                )
            if tol.get("cl_abs_diff_p99_max") is not None:
                p99 = float(np.percentile((merged["cl_a"] - merged["cl_b"]).abs(), 99))
                lim = tol["cl_abs_diff_p99_max"]
                rep.add(
                    sample,
                    "cl_abs_diff_p99",
                    "PASS" if p99 <= lim else "FAIL",
                    f"p99|Δcl|={p99:.2f}",
                    p99,
                    lim,
                )
    # aggregate charged fraction (works even without sidecars)
    ca = (s_a or {}).get("charging", {}).get("charged_fraction")
    cb = (s_b or {}).get("charging", {}).get("charged_fraction")
    if (
        ca is not None
        and cb is not None
        and tol.get("charged_fraction_abs_tol") is not None
    ):
        d = abs(ca - cb)
        lim = tol["charged_fraction_abs_tol"]
        rep.add(
            sample,
            "charged_fraction",
            "PASS" if d <= lim else "FAIL",
            f"|Δ|={d:.4f} (base={ca:.3f}, cand={cb:.3f})",
            d,
            lim,
        )


def check_cpm(rep, sample, s_a, s_b, tol):
    pa = (s_a or {}).get("cpm", {}).get("per_trna")
    pb = (s_b or {}).get("cpm", {}).get("per_trna")
    if not pa or not pb:
        rep.add(sample, "cpm", "SKIP", "cpm table missing on one side")
        return
    min_reads = tol.get("min_reads_for_cpm_check", 0)
    rel_tol = tol.get("cpm_rel_tol")
    if rel_tol is None:
        return
    worst_trna, worst = None, 0.0
    n_checked = 0
    for trna in set(pa) & set(pb):
        ra, rb = pa[trna], pb[trna]
        reads = ra.get("counts_charged", 0) + ra.get("counts_uncharged", 0)
        if reads < min_reads:
            continue
        for col in ("cpm_charged", "cpm_uncharged"):
            va, vb = ra.get(col), rb.get(col)
            if va is None or vb is None:
                continue
            denom = max(abs(va), 1e-9)
            rel = abs(va - vb) / denom
            n_checked += 1
            if rel > worst:
                worst, worst_trna = rel, f"{trna}/{col}"
    if n_checked == 0:
        rep.add(sample, "cpm", "SKIP", f"no tRNA with >= {min_reads} reads to check")
        return
    rep.add(
        sample,
        "cpm_rel",
        "PASS" if worst <= rel_tol else "FAIL",
        f"worst rel Δ={worst:.4f} at {worst_trna} ({n_checked} values checked)",
        worst,
        rel_tol,
    )


def check_bcerror(rep, sample, base_dir, cand_dir, tol):
    da, db = _perpos(base_dir, sample, "bcerror"), _perpos(cand_dir, sample, "bcerror")
    if da is None or db is None:
        rep.add(sample, "bcerror", "SKIP", "bcerror sidecar missing on one side")
        return
    lim = tol.get("mean_abs_freq_diff_max")
    if lim is None:
        return
    m = da.merge(db, on=["Reference", "Position"], suffixes=("_a", "_b"))
    if m.empty:
        rep.add(sample, "bcerror", "SKIP", "no shared positions")
        return
    diffs = []
    for c in ("MismatchFreq", "InsertionFreq", "DeletionFreq", "BCErrorFreq"):
        if f"{c}_a" in m and f"{c}_b" in m:
            diffs.append((m[f"{c}_a"] - m[f"{c}_b"]).abs())
    mad = float(pd.concat(diffs).mean()) if diffs else 0.0
    rep.add(
        sample,
        "bcerror_mad",
        "PASS" if mad <= lim else "FAIL",
        f"mean|Δfreq|={mad:.5f} over {len(m)} shared positions",
        mad,
        lim,
    )


def check_modkit(rep, sample, base_dir, cand_dir, tol):
    da, db = _perpos(base_dir, sample, "pileup"), _perpos(cand_dir, sample, "pileup")
    if da is None or db is None:
        rep.add(sample, "modkit", "SKIP", "pileup sidecar missing on one side")
        return
    lim = tol.get("frac_mod_mean_abs_diff_max")
    if lim is None:
        return
    keys = ["chrom", "pos", "mod_code"]
    m = da.merge(db, on=keys, suffixes=("_a", "_b"))
    if m.empty:
        rep.add(sample, "modkit", "SKIP", "no shared modified sites")
        return
    mad = float((m["frac_mod_a"] - m["frac_mod_b"]).abs().mean())
    rep.add(
        sample,
        "modkit_frac_mad",
        "PASS" if mad <= lim else "FAIL",
        f"mean|Δfrac_mod|={mad:.5f} over {len(m)} shared sites",
        mad,
        lim,
    )


# --- reporting -----------------------------------------------------------------


def render_text(rep, base_fp, cand_fp, profile) -> str:
    L = []
    L.append("=" * 78)
    L.append(f"benchmark comparison   profile={profile}")
    L.append(f"  baseline : {base_fp.get('label')}  (ref={base_fp.get('git_ref')})")
    L.append(f"  candidate: {cand_fp.get('label')}  (ref={cand_fp.get('git_ref')})")
    L.append("=" * 78)

    ta, tb = base_fp.get("manifest_tools"), cand_fp.get("manifest_tools")
    if ta or tb:
        L.append("tool versions (baseline -> candidate):")
        keys = sorted(set(ta or {}) | set(tb or {}))
        for k in keys:
            va = (ta or {}).get(k, "-")
            vb = (tb or {}).get(k, "-")
            flag = "" if va == vb else "   <-- changed"
            L.append(f"  {k:<20} {va} -> {vb}{flag}")
        L.append("")

    ra, rb = base_fp.get("runtimes") or {}, cand_fp.get("runtimes") or {}
    if ra or rb:
        L.append("runtime — wall seconds, baseline -> candidate (informational):")
        for rule in sorted(set(ra) | set(rb)):
            sa = ra.get(rule, {}).get("wall_s")
            sb = rb.get(rule, {}).get("wall_s")
            if sa is None or sb is None:
                L.append(f"  {rule:<26} {sa if sa is not None else '-'} -> "
                         f"{sb if sb is not None else '-'}  (one side only)")
                continue
            pct = ((sb - sa) / sa * 100) if sa else float("inf")
            arrow = "faster" if sb < sa else "slower" if sb > sa else "same"
            L.append(f"  {rule:<26} {sa:>8.1f} -> {sb:>8.1f}  "
                     f"({pct:+.0f}% {arrow})")
        L.append("")

    samples = sorted({c.sample for c in rep.checks})
    for sample in samples:
        L.append(f"[{sample}]")
        for c in [c for c in rep.checks if c.sample == sample]:
            mark = {"PASS": "ok  ", "FAIL": "FAIL", "SKIP": "skip"}[c.status]
            L.append(f"  {mark}  {c.metric:<22} {c.detail}")
        L.append("")

    n_fail, n_pass, n_skip = len(rep.failed), len(rep.passed), len(rep.skipped)
    L.append("-" * 78)
    L.append(f"SUMMARY: {n_pass} passed, {n_fail} failed, {n_skip} skipped")
    if n_fail:
        L.append(
            "RESULT: FAIL — results diverged beyond tolerance (see FAIL rows above)"
        )
    else:
        L.append("RESULT: PASS — results preserved within the '%s' profile" % profile)
    L.append("-" * 78)
    return "\n".join(L)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("baseline", type=Path, help="baseline fingerprint dir")
    ap.add_argument("candidate", type=Path, help="candidate fingerprint dir")
    ap.add_argument("--profile", default="strict", help="tolerance profile name")
    ap.add_argument(
        "--tolerances", type=Path, default=Path(__file__).parent / "tolerances.yml"
    )
    ap.add_argument("--json", type=Path, default=None, help="write machine report here")
    ap.add_argument(
        "--report-only",
        action="store_true",
        help="always exit 0 (print report, do not gate)",
    )
    args = ap.parse_args()

    profiles = yaml.safe_load(open(args.tolerances))
    if args.profile not in profiles:
        ap.error(
            f"profile '{args.profile}' not in {args.tolerances} "
            f"(have: {', '.join(profiles)})"
        )
    tol = profiles[args.profile]

    base_fp = lib.read_json(args.baseline / "fingerprint.json")
    cand_fp = lib.read_json(args.candidate / "fingerprint.json")

    rep = Report()
    samples = sorted(set(base_fp["samples"]) | set(cand_fp["samples"]))
    for sample in samples:
        s_a = base_fp["samples"].get(sample)
        s_b = cand_fp["samples"].get(sample)
        if s_a is None or s_b is None:
            rep.add(
                sample,
                "sample_presence",
                "FAIL",
                "sample only present in %s"
                % ("baseline" if s_b is None else "candidate"),
            )
            continue
        check_align(
            rep,
            sample,
            s_a.get("align_stats"),
            s_b.get("align_stats"),
            tol.get("align_stats", {}),
        )
        check_charging(
            rep,
            sample,
            args.baseline,
            args.candidate,
            s_a,
            s_b,
            tol.get("charging", {}),
        )
        check_cpm(rep, sample, s_a, s_b, tol.get("cpm", {}))
        check_bcerror(
            rep, sample, args.baseline, args.candidate, tol.get("bcerror", {})
        )
        check_modkit(rep, sample, args.baseline, args.candidate, tol.get("modkit", {}))

    text = render_text(rep, base_fp, cand_fp, args.profile)
    print(text)

    if args.json:
        lib.write_json(
            args.json,
            {
                "profile": args.profile,
                "baseline": base_fp.get("label"),
                "candidate": cand_fp.get("label"),
                "n_pass": len(rep.passed),
                "n_fail": len(rep.failed),
                "n_skip": len(rep.skipped),
                "checks": [c.__dict__ for c in rep.checks],
            },
        )

    if rep.failed and not args.report_only:
        sys.exit(1)


if __name__ == "__main__":
    main()
