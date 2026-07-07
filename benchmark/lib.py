"""
Shared helpers for the aa-tRNA-seq benchmark harness.

Readers here know the exact on-disk schema of the pipeline's summary outputs
(see workflow/rules/*.smk and workflow/scripts/*). They are deliberately
tolerant of missing files: a fingerprint is a best-effort snapshot of whatever
a run produced, and comparison only diffs metrics present in *both* sides.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd
import pysam


# --- layout of a pipeline output_directory -------------------------------------


def tables_dir(outdir: Path) -> Path:
    return outdir / "summary" / "tables"


def modkit_dir(outdir: Path) -> Path:
    return outdir / "summary" / "modkit"


def discover_samples(outdir: Path) -> list[str]:
    """Sample IDs are the subdirectory names under summary/tables/."""
    td = tables_dir(outdir)
    if not td.is_dir():
        return []
    return sorted(p.name for p in td.iterdir() if p.is_dir())


def _table_path(outdir: Path, sample: str, suffix: str) -> Path:
    return tables_dir(outdir) / sample / f"{sample}.{suffix}"


# --- individual output readers -------------------------------------------------


def read_charging_prob(outdir: Path, sample: str) -> pd.DataFrame | None:
    """Per-read charging table: columns read_id, tRNA, charging_likelihood."""
    p = _table_path(outdir, sample, "charging_prob.tsv.gz")
    if not p.exists():
        return None
    df = pd.read_csv(p, sep="\t")
    df = df.rename(columns={"charging_likelihood": "cl"})
    return df[["read_id", "tRNA", "cl"]]


def read_charging_cpm(outdir: Path, sample: str) -> pd.DataFrame | None:
    """Per-isodecoder counts/CPM table, indexed by tRNA."""
    p = _table_path(outdir, sample, "charging.cpm.tsv.gz")
    if not p.exists():
        return None
    return pd.read_csv(p, sep="\t", index_col=0)


def read_align_stats(outdir: Path, sample: str) -> dict | None:
    """Headline alignment metrics.

    align_stats.tsv.gz has one row per pipeline stage (`info` column:
    unmapped/aligned/classified). Total read throughput comes from the
    `unmapped` row (all basecalled reads); mapping rate from the final
    (`classified`, else `aligned`) row — the `unmapped` row is the raw uBAM and
    always reports 0% mapped.
    """
    p = _table_path(outdir, sample, "align_stats.tsv.gz")
    if not p.exists():
        return None
    df = pd.read_csv(p, sep="\t")
    if df.empty:
        return None

    if "info" in df.columns:
        by_stage = {str(r["info"]): r for _, r in df.iterrows()}
        final = by_stage.get("classified")
        if final is None:
            final = by_stage.get("aligned")
        if final is None:
            final = df.iloc[-1]
        total = by_stage.get("unmapped")
        n_reads_row = total if total is not None else final
    else:
        final = df.iloc[0]
        n_reads_row = df.iloc[0]

    out = {}
    if "n_reads" in n_reads_row:
        out["n_reads"] = float(n_reads_row["n_reads"])
    for k in ("mapped_reads", "pct_mapped"):
        if k in final:
            out[k] = float(final[k])
    return out


def read_bcerror(outdir: Path, sample: str) -> pd.DataFrame | None:
    """Per-position base-calling error frequencies."""
    p = _table_path(outdir, sample, "bcerror.tsv.gz")
    if not p.exists():
        return None
    cols = [
        "Reference",
        "Position",
        "MismatchFreq",
        "InsertionFreq",
        "DeletionFreq",
        "BCErrorFreq",
    ]
    df = pd.read_csv(p, sep="\t")
    keep = [c for c in cols if c in df.columns]
    return df[keep]


def read_pileup(outdir: Path, sample: str) -> pd.DataFrame | None:
    """modkit pileup bedMethyl -> (chrom, pos, mod_code, valid_cov, frac_mod).

    bedMethyl columns (0-based): 0 chrom, 1 start, 3 mod_code, 4 score,
    9 valid_cov, 10 percent_modified. `frac_mod` is normalized to 0-1 so
    tolerances read on a fraction scale.
    """
    p = modkit_dir(outdir) / sample / f"{sample}.pileup.bed.gz"
    if not p.exists():
        return None
    df = pd.read_csv(p, sep="\t", header=None)
    if df.shape[1] < 11:
        return None
    out = pd.DataFrame(
        {
            "chrom": df[0],
            "pos": df[1],
            "mod_code": df[3],
            "valid_cov": df[9],
            "frac_mod": df[10].astype(float) / 100.0,
        }
    )
    return out


def read_manifest(outdir: Path) -> dict | None:
    p = outdir / "manifest.json"
    if not p.exists():
        return None
    with open(p) as fh:
        return json.load(fh)


def read_benchmarks(outdir: Path) -> dict:
    """Aggregate Snakemake `benchmark:` TSVs under {outdir}/benchmarks/.

    Layout is benchmarks/<rule>/<name>.tsv, each with a header row and one data
    row (cols: s, h:m:s, max_rss, ...). Returns
    {rule: {"wall_s": total_seconds, "max_rss_mb": peak, "n": count}}.
    """
    bdir = outdir / "benchmarks"
    if not bdir.is_dir():
        return {}
    out: dict[str, dict] = {}
    for rule_dir in sorted(p for p in bdir.iterdir() if p.is_dir()):
        wall = 0.0
        rss = 0.0
        n = 0
        for tsv in rule_dir.glob("*.tsv"):
            try:
                df = pd.read_csv(tsv, sep="\t")
            except Exception:
                continue
            if df.empty:
                continue
            row = df.iloc[0]
            if "s" in row:
                wall += float(row["s"])
                n += 1
            if "max_rss" in row and pd.notna(row["max_rss"]):
                rss = max(rss, float(row["max_rss"]))
        if n:
            out[rule_dir.name] = {
                "wall_s": round(wall, 2),
                "max_rss_mb": round(rss, 1),
                "n": n,
            }
    return out


# --- BAM tag extraction (used by classifier equivalence) -----------------------


def read_bam_tag_by_read(bam_path: Path, tag: str) -> dict[str, float]:
    """Map read_id -> scalar value of `tag` for every primary alignment.

    Mirrors get_charging_table.py: a charging tag of 0 is valid and kept;
    multi-element array tags (e.g. dorado ML mod-base tags) are skipped.
    Falls back to the uppercase tag name for older BAMs.
    """
    out: dict[str, float] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            tags = dict(read.tags)
            val = tags.get(tag)
            if val is None and tag.islower():
                val = tags.get(tag.upper())
            if val is None:
                continue
            if hasattr(val, "__len__") and not isinstance(val, str):
                if len(val) != 1:
                    continue
                val = val[0]
            out[read.query_name] = float(val)
    return out


# --- small utilities -----------------------------------------------------------


def write_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)
        fh.write("\n")


def read_json(path: Path) -> dict:
    with open(path) as fh:
        return json.load(fh)


def openmaybe_gz(path: Path, mode: str = "rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)
