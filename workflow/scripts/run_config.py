#!/usr/bin/env python
"""
Write and check a run's config pair: `config/config-<name>.yml` plus its samples
file. The one place that knows every rule a samples file has to satisfy, so a
config can be validated before a GPU-hour is spent finding out.

    pixi run new-run-config -- --name mouse-liver \\
        --output-dir results/mouse-liver \\
        --run /data/runs/20260910_liver/20260910_1200_P2S-01617-A_PBK12345_abcdef \\
        --sample wt_rep1=ldx01+fdx01 --sample wt_rep2=ldx02+fdx01 \\
        --sample ko_rep1=ldx04+fdx02 \\
        --reference-raw resources/ref/mm39-mature-tRNAs-collapsed.fa \\
        --three-prime plain=GGCTTCTTCTTGCTCTTAGGAAAAAAAAAA \\
        --dry-run

    pixi run check-run-config -- config/config-mouse-liver.yml

Sample specs name the barcode(s) a library carries, on whichever axes apply:

    ldx01              3' LDX code
    ldx01+fdx01        3' LDX code and 5' FDX code (dual index)
    ldx01/edx01        3' LDX code, keep only the edx01 3' adapter
    ldx01+fdx01/edx07  all three
    wdx:barcode03      WarpDemuX barcode (the retired backend)
    -                  unbarcoded: the whole run is this sample

`new` writes two files the pipeline can run, with every override of
config-base.yml explained in a comment beside it, and refuses a plan the
pipeline would refuse. `check` applies the same rules to an existing pair and,
with --dry-run, builds the DAG. Both print a plain report, or JSON with --json
for a caller that is a program.

What is checked, and why each is here:

  * every code against the bundle that has to emit it (barcode_names): a
    misspelt code surfaces otherwise as "no reads were assigned" after demux
  * fdx needs ldx; samples sharing an ldx code all name an fdx or none do;
    tuples are unique per run -- the rules select_demux_reads.py enforces
  * edx names must be declared 3' adapters, and 3' adapters must share one
    length and start with GGC (the CCA|GGC junction the charging model anchors
    on); this is the shape of failure #120 was
  * run directories exist and hold pod5_pass/, pod5_fail/ or pod5/; the raw
    reference exists and ends in CCA
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import subprocess
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve()
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
from barcode_names import bundle_barcode_names  # noqa: E402

POD5_DIRS = ("pod5_pass", "pod5_fail", "pod5")
DEFAULT_FIVE_PRIME = "CCTAAGAGCAAGAAGAAGCCTGG"


class PlanError(ValueError):
    pass


# ------------------------------------------------------------------ specs ---


def parse_spec(spec):
    """`ldx01+fdx01/edx07` -> {"ldx": "ldx01", "fdx": "fdx01", "edx": "edx07"}.

    `-` (or empty) is an unbarcoded sample. `wdx:barcode03` is a WarpDemuX
    barcode. Order of `+` parts does not matter; each axis at most once.
    """
    spec = (spec or "").strip()
    if spec in ("", "-", "~", "null", "none"):
        return {}
    if spec.startswith("wdx:"):
        code = spec[len("wdx:") :]
        if not code:
            raise PlanError(f"{spec!r}: wdx: needs a barcode name")
        return {"wdx": code}
    codes, _, edx = spec.partition("/")
    out = {}
    for part in codes.split("+"):
        if not part:
            raise PlanError(f"{spec!r}: empty code")
        m = re.fullmatch(r"(ldx|fdx)(\d+)", part)
        if not m:
            raise PlanError(
                f"{spec!r}: {part!r} is not an ldxNN or fdxNN code (use wdx:NAME for "
                "WarpDemuX barcodes, or - for an unbarcoded sample)"
            )
        if m.group(1) in out:
            raise PlanError(f"{spec!r}: {m.group(1)} given twice")
        out[m.group(1)] = part
    if edx:
        out["edx"] = edx
    return out


# ------------------------------------------------------------------- plan ---


def load_base_config(repo_root=REPO_ROOT):
    return yaml.safe_load((repo_root / "config" / "config-base.yml").read_text())


def three_prime_list(adapters_cfg):
    """[(name, seq)] for `adapters.three_prime` in either of its shapes."""
    tp = (adapters_cfg or {}).get("three_prime")
    if tp is None:
        return []
    if isinstance(tp, str):
        return [("default", tp)]
    return [(a["name"], a["seq"]) for a in tp]


def validate_plan(plan, base, repo_root=REPO_ROOT):
    """Problems with a plan, as (level, message) with level `error` or `warning`.

    `plan`: {"runs": [{"path": str, "samples": {name: codes}}], "adapters": {...}
    or None, "reference": {"mode","raw_fasta"|"fasta"} or None, "ldx_gpu": bool}
    """
    problems = []
    err = lambda m: problems.append(("error", m))  # noqa: E731
    warn = lambda m: problems.append(("warning", m))  # noqa: E731

    all_codes = [c for run in plan["runs"] for c in run["samples"].values()]
    uses = {axis for c in all_codes for axis in c}
    if "wdx" in uses and ("ldx" in uses or "fdx" in uses):
        err(
            "a run mixes WarpDemuX (wdx:) and escpod (ldx/fdx) samples; the backends are mutually exclusive"
        )
    if "fdx" in uses and "ldx" not in uses:
        err(
            "fdx codes without any ldx code: the FDX axis is joined against the LDX call, so dual-index samples need both"
        )

    # codes against the bundles that emit them
    ldx_model = repo_root / base["ldx"]["model"]
    fdx_model = repo_root / base["fdx"]["model"]
    ldx_names = bundle_barcode_names(ldx_model) if ldx_model.is_dir() else set()
    fdx_names = bundle_barcode_names(fdx_model) if fdx_model.is_dir() else set()
    adapters = plan.get("adapters") or base.get("adapters", {})
    edx_names = {name for name, _ in three_prime_list(adapters)}

    seen_names = set()
    for run in plan["runs"]:
        path = Path(run["path"])
        if not path.is_dir():
            err(f"run directory does not exist: {path}")
        elif not any((path / d).is_dir() for d in POD5_DIRS):
            err(
                f"{path} holds none of {', '.join(POD5_DIRS)}; find_raw_inputs would find no POD5"
            )
        by_ldx = {}
        for name, codes in run["samples"].items():
            # An unbarcoded sample may span several runs (the TSV form pools a
            # sample's paths); a barcoded name must be unique.
            if name in seen_names and codes:
                err(f"sample name {name!r} used twice")
            seen_names.add(name)
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
                err(
                    f"sample name {name!r}: use letters, digits, _ . - only (it becomes a path and a wildcard)"
                )
            if "ldx" in codes and ldx_names and codes["ldx"] not in ldx_names:
                err(
                    f"{name}: {codes['ldx']} is not emitted by {ldx_model.name} ({', '.join(sorted(ldx_names))})"
                )
            if "fdx" in codes and fdx_names and codes["fdx"] not in fdx_names:
                err(
                    f"{name}: {codes['fdx']} is not emitted by {fdx_model.name} ({', '.join(sorted(fdx_names))})"
                )
            if "fdx" in codes and "ldx" not in codes:
                err(f"{name}: fdx without ldx")
            if "edx" in codes and codes["edx"] not in edx_names:
                err(
                    f"{name}: edx {codes['edx']!r} is not a declared 3' adapter ({', '.join(sorted(edx_names)) or 'none'}); EDX names must match adapters.three_prime"
                )
            if "ldx" in codes:
                by_ldx.setdefault(codes["ldx"], []).append((name, codes))
        for code, group in by_ldx.items():
            with_fdx = [n for n, c in group if "fdx" in c]
            if with_fdx and len(with_fdx) != len(group):
                err(
                    f"samples {', '.join(n for n, _ in group)} share {code} but only {', '.join(with_fdx)} name an fdx; all or none"
                )
            tuples = [(c.get("fdx"), c.get("edx")) for _, c in group]
            if len(set(tuples)) != len(tuples):
                err(
                    f"samples {', '.join(n for n, _ in group)} have identical barcode assignments"
                )
        if len(run["samples"]) > 1 and any(not c for c in run["samples"].values()):
            err(
                f"{path}: an unbarcoded sample (-) claims the whole run, so it cannot share a run with barcoded samples"
            )

    # adapters
    tp = three_prime_list(adapters)
    lengths = {len(seq) for _, seq in tp}
    if len(lengths) > 1:
        err(
            f"3' adapters differ in length ({sorted(lengths)}); reference validation requires one length"
        )
    for name, seq in tp:
        if not seq.upper().startswith("GGC"):
            err(
                f"3' adapter {name} does not start with GGC: the charging model anchors on the CCA|GGC junction"
            )
        if not re.fullmatch(r"[ACGTN]+", seq.upper()):
            err(f"3' adapter {name} has non-ACGTN characters")
    fp = adapters.get("five_prime", DEFAULT_FIVE_PRIME)
    if not re.fullmatch(r"[ACGTN]+", str(fp).upper()):
        err("5' adapter has non-ACGTN characters")

    # reference
    ref = plan.get("reference")
    if ref:
        if ref.get("mode") == "build":
            raw = Path(ref["raw_fasta"])
            if not raw.is_file():
                err(f"raw reference does not exist: {raw}")
            else:
                seqs = list(raw.read_text().split(">")[1:])
                bodies = ["".join(s.split("\n")[1:]).upper() for s in seqs]
                if not bodies:
                    err(f"{raw}: no sequences")
                if any(DEFAULT_FIVE_PRIME in b for b in bodies):
                    err(
                        f"{raw}: sequences already carry the 5' adapter; use --fasta (validate mode), not --reference-raw"
                    )
                no_cca = sum(1 for b in bodies if not b.endswith("CCA"))
                if no_cca:
                    warn(
                        f"{raw}: {no_cca} of {len(bodies)} sequences do not end in CCA; build mode appends it with a warning"
                    )
        else:
            fa = Path(ref["fasta"])
            if not fa.is_file():
                err(f"reference does not exist: {fa}")
    return problems


# ----------------------------------------------------------------- render ---


def render_samples(plan):
    """The samples file text and its extension (`yml`, or `tsv` when unbarcoded)."""
    barcoded = any(c for run in plan["runs"] for c in run["samples"].values())
    if not barcoded:
        lines = [
            f"{name}\t{run['path']}" for run in plan["runs"] for name in run["samples"]
        ]
        return "\n".join(lines) + "\n", "tsv"
    out = [
        f"# Samples for `{plan['name']}`, written by run_config.py on {dt.date.today()}.",
        "# One entry per sample: the barcode(s) its library carries. `ldx` is the 3'",
        "# LDX code, `fdx` the 5' FDX code (dual index), `edx` keeps only reads with",
        "# that 3' adapter. Check with: pixi run check-run-config -- config/config-"
        f"{plan['name']}.yml",
        "",
        "runs:",
    ]
    for run in plan["runs"]:
        out.append(f"  - path: {run['path']}")
        out.append("    samples:")
        for name, codes in run["samples"].items():
            if not codes:
                out.append(f"      {name}: ~")
            elif "wdx" in codes:
                fields = [f'wdx: "{codes["wdx"]}"']
                if "edx" in codes:
                    fields.append(f'edx: "{codes["edx"]}"')
                out.append(f"      {name}: {{ {', '.join(fields)} }}")
            else:
                fields = [
                    f'{axis}: "{codes[axis]}"'
                    for axis in ("ldx", "fdx", "edx")
                    if axis in codes
                ]
                out.append(f"      {name}: {{ {', '.join(fields)} }}")
    return "\n".join(out) + "\n", "yml"


def render_config(plan, samples_file):
    uses = {axis for run in plan["runs"] for c in run["samples"].values() for axis in c}
    out = [
        f"# `{plan['name']}`: pipeline config, written by run_config.py on {dt.date.today()}.",
        "#",
        "# Only what differs from config/config-base.yml is set here; everything else",
        "# (basecalling model, charging bundle, escpod version, demux gates) is",
        "# inherited and should be changed THERE, deliberately, not per project.",
        "#",
        f"#   pixi run check-run-config -- config/config-{plan['name']}.yml   # validate + dry-run",
        f"#   sbatch ... snakemake --configfile=config/config-{plan['name']}.yml --profile=cluster/slurm",
        "",
        f"samples: {samples_file}",
        f'output_directory: "{plan["output_dir"]}"',
    ]
    if "wdx" in uses:
        out += [
            "",
            "# WarpDemuX (WDX) barcodes: the retired backend; see config/README.md.",
            "warpdemux:",
            "    enabled: true",
        ]
    if "ldx" in uses:
        out += [
            "",
            "# 3' LDX index, called by escpod on the raw signal. Gates, boundary margin and",
            "# the bundle come from config-base.yml.",
            "ldx:",
            "    enabled: true",
            f"    gpu: {'true' if plan.get('ldx_gpu', True) else 'false'}",
        ]
        if plan.get("threads"):
            out.append(f"    threads: {plan['threads']}")
    if "fdx" in uses:
        out += [
            "",
            "# 5' FDX index, the second axis. A dual-index sample owns the reads on which",
            "# both calls agree. `fused` is off until escpod accepts the LDX boundary",
            "# flags in a fused run; see the fdx block in config-base.yml.",
            "fdx:",
            "    enabled: true",
        ]
    if "edx" in uses:
        out += [
            "",
            "# `edx:` on a sample keeps only reads carrying that 3' adapter, and",
            "# `edx.enabled` emits summary/edx/edx_concordance.tsv.gz.",
            "edx:",
            "    enabled: true",
        ]
    if plan.get("adapters"):
        ad = plan["adapters"]
        out += [
            "",
            "# Adapters actually in these libraries. Overrides the base set: names must",
            "# match `edx:` assignments, all 3' adapters must be one length and start",
            "# with GGC (the CCA|GGC junction the charging model anchors on).",
            "adapters:",
            f'    five_prime: "{ad.get("five_prime", DEFAULT_FIVE_PRIME)}"',
        ]
        tp = three_prime_list(ad)
        if len(tp) == 1 and tp[0][0] == "default":
            out.append(f'    three_prime: "{tp[0][1]}"')
        else:
            out.append("    three_prime:")
            for name, seq in tp:
                out += [f'        - name: "{name}"', f'          seq: "{seq}"']
    ref = plan.get("reference")
    if ref and ref.get("mode") == "build":
        out += [
            "",
            "# Raw tRNA sequences (no adapters, ending in CCA); the pipeline builds the",
            "# adapted reference and its index under output_directory/reference/.",
            "reference:",
            '    mode: "build"',
            f'    raw_fasta: "{ref["raw_fasta"]}"',
        ]
    elif ref:
        out += [
            "",
            "# An already-adapted reference (5' adapter + N + tRNA + CCA + 3' adapter);",
            "# validated, not rebuilt.",
            f'fasta: "{ref["fasta"]}"',
        ]
    if plan.get("cleanup"):
        out += [
            "",
            "# Delete large regenerable intermediates as the run goes; see config/README.md.",
            f"cleanup_intermediates: [{', '.join(plan['cleanup'])}]",
        ]
    return "\n".join(out) + "\n"


# ------------------------------------------------------------------ check ---


def load_pair(config_path, repo_root=REPO_ROOT):
    """A plan reconstructed from an existing config + samples file."""
    cfg = yaml.safe_load(Path(config_path).read_text())
    base = load_base_config(repo_root)
    merged = {**base, **cfg}
    for key in ("ldx", "fdx", "edx", "warpdemux", "adapters", "reference"):
        if isinstance(base.get(key), dict) and isinstance(cfg.get(key), dict):
            merged[key] = {**base[key], **cfg[key]}
    samples_path = Path(merged["samples"])
    if not samples_path.is_absolute():
        samples_path = repo_root / samples_path
    runs = []
    if samples_path.suffix in (".yml", ".yaml"):
        data = yaml.safe_load(samples_path.read_text())
        for run in data.get("runs", []):
            samples = {}
            for name, val in (run.get("samples") or {}).items():
                if val is None:
                    samples[name] = {}
                elif isinstance(val, str):
                    samples[name] = (
                        {"wdx": val}
                        if merged.get("warpdemux", {}).get("enabled")
                        else {"ldx": val}
                    )
                else:
                    samples[name] = {k: v for k, v in val.items() if v is not None}
            runs.append({"path": run["path"], "samples": samples})
    else:
        by_path = {}
        for line in samples_path.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            name, path = line.split()
            by_path.setdefault(path, {})[name] = {}
        runs = [{"path": p, "samples": s} for p, s in by_path.items()]
    ref = None
    if merged.get("reference", {}).get("mode") == "build":
        ref = {"mode": "build", "raw_fasta": merged["reference"]["raw_fasta"]}
    elif merged.get("fasta"):
        ref = {"mode": "validate", "fasta": merged["fasta"]}
    plan = {
        "name": Path(config_path).stem.replace("config-", "", 1),
        "output_dir": merged.get("output_directory"),
        "runs": runs,
        "adapters": cfg.get("adapters"),
        "reference": ref,
        "ldx_gpu": merged.get("ldx", {}).get("gpu", True),
    }
    flags = []
    uses = {axis for run in runs for c in run["samples"].values() for axis in c}
    if "ldx" in uses and not merged.get("ldx", {}).get("enabled"):
        flags.append(("error", "samples name ldx codes but ldx.enabled is off"))
    if "fdx" in uses and not merged.get("fdx", {}).get("enabled"):
        flags.append(("error", "samples name fdx codes but fdx.enabled is off"))
    if "wdx" in uses and not merged.get("warpdemux", {}).get("enabled"):
        flags.append(
            ("error", "samples name wdx barcodes but warpdemux.enabled is off")
        )
    if merged.get("fdx", {}).get("enabled") and not merged.get("ldx", {}).get(
        "enabled"
    ):
        flags.append(("error", "fdx.enabled without ldx.enabled"))
    return plan, merged, flags


def dry_run(config_path, repo_root=REPO_ROOT):
    """Build the DAG with snakemake -n; returns (ok, job table text)."""
    proc = subprocess.run(
        ["snakemake", "-n", "--configfile", str(config_path)],
        capture_output=True,
        text=True,
        cwd=repo_root,
    )
    text = proc.stdout + proc.stderr
    m = re.search(r"^job\s+count.*?^total\s+\d+", text, re.DOTALL | re.MULTILINE)
    table = m.group(0) if m else text[-3000:]
    return proc.returncode == 0, table


# ------------------------------------------------------------------- main ---


class _Run(argparse.Action):
    def __call__(self, parser, ns, value, option_string=None):
        ns.runs.append({"path": value, "samples": {}})


class _Sample(argparse.Action):
    def __call__(self, parser, ns, value, option_string=None):
        if not ns.runs:
            parser.error("--sample must follow the --run it belongs to")
        name, sep, spec = value.partition("=")
        if not sep:
            parser.error(f"--sample expects NAME=SPEC, got {value!r}")
        ns.runs[-1]["samples"][name] = spec


def report(problems, as_json, extra=None):
    if as_json:
        print(
            json.dumps(
                {
                    "problems": [{"level": lv, "message": m} for lv, m in problems],
                    **(extra or {}),
                },
                indent=2,
            )
        )
    else:
        for level, msg in problems:
            print(f"{level.upper()}: {msg}")
        for key, value in (extra or {}).items():
            print(
                f"{key}: {value}"
                if not isinstance(value, str) or "\n" not in value
                else f"{key}:\n{value}"
            )
    return 1 if any(lv == "error" for lv, _ in problems) else 0


def cmd_new(args):
    if not args.runs:
        sys.exit("at least one --run (with its --sample entries) is required")
    try:
        runs = [
            {
                "path": run["path"],
                "samples": {n: parse_spec(s) for n, s in run["samples"].items()},
            }
            for run in args.runs
        ]
    except PlanError as exc:
        sys.exit(f"ERROR: {exc}")
    adapters = None
    if args.three_prime or args.five_prime:
        adapters = {"five_prime": args.five_prime or DEFAULT_FIVE_PRIME}
        if args.three_prime:
            named = [t.partition("=") for t in args.three_prime]
            if len(named) == 1 and not named[0][1]:
                adapters["three_prime"] = named[0][0]
            else:
                adapters["three_prime"] = [
                    {"name": n or "default", "seq": s} for n, _, s in named
                ]
    reference = None
    if args.reference_raw:
        reference = {"mode": "build", "raw_fasta": args.reference_raw}
    elif args.fasta:
        reference = {"mode": "validate", "fasta": args.fasta}
    plan = {
        "name": args.name,
        "output_dir": args.output_dir or f"results/{args.name}",
        "runs": runs,
        "adapters": adapters,
        "reference": reference,
        "ldx_gpu": not args.no_gpu,
        "threads": args.threads,
        "cleanup": args.cleanup,
    }
    problems = validate_plan(plan, load_base_config(args.repo_root), args.repo_root)
    samples_text, ext = render_samples(plan)
    samples_file = f"config/samples-{args.name}.{ext}"
    config_file = f"config/config-{args.name}.yml"
    if any(lv == "error" for lv, _ in problems) and not args.force:
        return report(problems, args.json, {"written": []})
    for rel, text in (
        (samples_file, samples_text),
        (config_file, render_config(plan, samples_file)),
    ):
        path = args.repo_root / rel
        if path.exists() and not args.force:
            problems.append(("error", f"{rel} exists; pass --force to overwrite"))
            return report(problems, args.json, {"written": []})
        path.write_text(text)
    extra = {"written": [samples_file, config_file]}
    if args.dry_run:
        ok, table = dry_run(args.repo_root / config_file, args.repo_root)
        extra["dry_run_ok"] = ok
        extra["dag"] = table
        if not ok:
            problems.append(("error", "snakemake dry-run failed; see dag"))
    return report(problems, args.json, extra)


def cmd_check(args):
    plan, merged, flags = load_pair(args.config, args.repo_root)
    problems = flags + validate_plan(plan, merged, args.repo_root)
    extra = {
        "samples": sum(len(r["samples"]) for r in plan["runs"]),
        "runs": len(plan["runs"]),
        "output_directory": plan["output_dir"],
    }
    if args.dry_run:
        ok, table = dry_run(args.config, args.repo_root)
        extra["dry_run_ok"] = ok
        extra["dag"] = table
        if not ok:
            problems.append(("error", "snakemake dry-run failed; see dag"))
    return report(problems, args.json, extra)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--repo-root", type=Path, default=REPO_ROOT, help=argparse.SUPPRESS
    )
    parser.add_argument("--json", action="store_true", help="machine-readable report")
    sub = parser.add_subparsers(dest="cmd", required=True)

    new = sub.add_parser(
        "new", help="write config/config-<name>.yml and its samples file"
    )
    new.add_argument(
        "--name",
        required=True,
        help="project name; files are config/{config,samples}-<name>",
    )
    new.add_argument("--output-dir", help="output_directory (default results/<name>)")
    new.add_argument(
        "--run",
        action=_Run,
        dest="runs",
        default=[],
        metavar="DIR",
        help="a sequencing run directory (holding pod5_pass/ pod5_fail/ pod5/); repeatable, followed by its --sample entries",
    )
    new.add_argument(
        "--sample",
        action=_Sample,
        dest="runs",
        metavar="NAME=SPEC",
        help="a sample of the preceding --run; SPEC as in the module docstring",
    )
    new.add_argument(
        "--reference-raw", help="raw tRNA FASTA (no adapters) -> reference.mode build"
    )
    new.add_argument("--fasta", help="already-adapted reference FASTA -> validate mode")
    new.add_argument(
        "--three-prime",
        action="append",
        metavar="[NAME=]SEQ",
        help="3' adapter(s) in these libraries; repeatable",
    )
    new.add_argument(
        "--five-prime", metavar="SEQ", help=f"5' adapter (default {DEFAULT_FIVE_PRIME})"
    )
    new.add_argument(
        "--no-gpu", action="store_true", help="run escpod demux on the CPU"
    )
    new.add_argument("--threads", type=int, help="ldx.threads")
    new.add_argument(
        "--cleanup", nargs="+", metavar="TIER", help="cleanup_intermediates tiers"
    )
    new.add_argument(
        "--dry-run", action="store_true", help="build the DAG after writing"
    )
    new.add_argument(
        "--force",
        action="store_true",
        help="overwrite existing files / write despite errors",
    )
    new.set_defaults(func=cmd_new)

    check = sub.add_parser("check", help="validate an existing config + samples pair")
    check.add_argument("config", type=Path)
    check.add_argument("--dry-run", action="store_true", help="also build the DAG")
    check.set_defaults(func=cmd_check)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
