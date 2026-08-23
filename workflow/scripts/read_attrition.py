#!/usr/bin/env python
"""Where a run's reads are lost, as one table, every run.

Each gate's loss is already implicit in some artifact, but only as a difference
between two rows in different files, so nobody looks. On the 2026-08-06 LDX run
the `aligned -> classified` gate was quietly discarding 12.37% and had been for
every prior run; it took reconstructing the cascade by hand to see it. Emitting
this by default is the cheap way to make the next such loss self-announcing.

Sources, all of which survive `cleanup_intermediates`:

  demux_summary.tsv.gz    barcode assignment, `unclassified` among the barcodes
  align_stats.tsv.gz      per sample: `unmapped` (post-basecall), `aligned`,
                          `classified` (carrying a `cl` tag)
  charging_calls.tsv.gz   per sample: one row per read the charging model saw,
                          with a `reason` on every read it did not score
  anchor_coverage.tsv.gz  per sample: how many aligned reads span the CCA
  edx_concordance.tsv.gz  per adapter, `none` among them

Two properties of this cascade are not obvious and are handled explicitly:

**Basecalling can emit MORE reads than it consumes** — dorado splits concatemers
into sub-reads with fresh IDs, so `basecalled` routinely exceeds `assigned` by
~0.5%. A negative loss is reported as `split_gain`, not clamped silently. It is
also the alarm for a POD5 written with non-uniform signal batches, where dorado
skips unreadable reads without erroring and this figure goes sharply negative
(escapepod-rs#195 cost one run 10.6% that way, and it looked like biology).

**The charge-calling gate says why, and it is not all one thing.** Reads the
model declines are split between ones it never saw (no signal, or an alignment
that never reaches the junction) and ones it saw and deliberately abstained on
(`no_aligned_arm` — bwa placed no base of the common arm, where the model scores
balanced accuracy 0.4993 and calls everything charged). Only the second kind is
charging-correlated, and it is the one that biases a charging fraction low, so
the two are reported apart rather than as a single unexplained difference.

**EDX detection is not always a gate.** When no sample carries an `edx:` key it
annotates rather than filters, so it is reported alongside the cascade rather
than inside it — folding it in would double-count reads that also survive
alignment.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--align-stats", type=Path, nargs="+", required=True)
    ap.add_argument("--anchor-coverage", type=Path, nargs="*", default=[])
    ap.add_argument("--charging-calls", type=Path, nargs="*", default=[])
    ap.add_argument("--demux-summary", type=Path, nargs="*", default=[])
    ap.add_argument("--edx-concordance", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    stage = {"unmapped": "basecalled", "aligned": "aligned", "classified": "classified"}
    totals = {"basecalled": 0, "aligned": 0, "classified": 0}
    for p in args.align_stats:
        df = _read_tsv(p)
        for _, row in df.iterrows():
            key = stage.get(row["info"])
            if key:
                totals[key] += int(row["n_reads"])

    sequenced = assigned = None
    for p in args.demux_summary:
        with gzip.open(p, "rt") as fh:
            fh.readline()
            counts = {}
            for line in fh:
                bc, n, _ = line.rstrip("\n").split("\t")
                counts[bc] = int(n)
        sequenced = (sequenced or 0) + sum(counts.values())
        assigned = (assigned or 0) + sum(
            v for k, v in counts.items() if k != "unclassified"
        )

    anchor = {
        "scored": 0,
        "covers_anchor": 0,
        "ends_before_anchor": 0,
        "starts_after_anchor": 0,
    }
    for p in args.anchor_coverage:
        df = _read_tsv(p)
        for k in anchor:
            anchor[k] += int(df[k].sum())

    # The charging model's own account of what it did not score. A blank
    # `reason` is a call; anything else names the population.
    n_called = 0
    no_call: dict[str, int] = {}
    for p in args.charging_calls:
        df = _read_tsv(p)
        if "reason" not in df.columns:
            continue
        reasons = df["reason"].fillna("")
        n_called += int((reasons == "").sum())
        for reason, n in reasons[reasons != ""].value_counts().items():
            no_call[reason] = no_call.get(reason, 0) + int(n)
    n_no_call = sum(no_call.values())
    n_anchored = n_called + n_no_call

    rows = []
    if sequenced is not None:
        rows.append(
            (
                "sequenced",
                "barcode assigned",
                sequenced,
                assigned,
                sequenced - assigned,
                "no barcode decoded",
            )
        )
        rows.append(
            (
                "barcode assigned",
                "basecalled",
                assigned,
                totals["basecalled"],
                max(0, assigned - totals["basecalled"]),
                "basecaller could not read the signal",
            )
        )
    rows.append(
        (
            "basecalled",
            "aligned",
            totals["basecalled"],
            totals["aligned"],
            totals["basecalled"] - totals["aligned"],
            "no alignment to the reference",
        )
    )
    rows.append(
        (
            "aligned",
            "charge-called",
            totals["aligned"],
            totals["classified"],
            totals["aligned"] - totals["classified"],
            "not emitted by the charging model",
        )
    )

    cascade = pd.DataFrame(
        rows,
        columns=["from_stage", "to_stage", "entered", "retained", "lost", "reason"],
    )
    cascade["pct_lost"] = (
        100 * cascade["lost"] / cascade["entered"].where(cascade["entered"] > 0)
    )

    # What the charging model itself reports, where it reported anything. This
    # is measured rather than inferred, so it takes precedence over the
    # structural estimate below.
    if n_anchored:
        if n_no_call:
            parts = ", ".join(
                f"{r} {n:,} ({100 * n / n_anchored:.2f}%)"
                for r, n in sorted(no_call.items(), key=lambda kv: -kv[1])
            )
            measured = (
                f"{n_no_call:,} of {n_anchored:,} anchored reads not scored "
                f"by the model — {parts}"
            )
        else:
            measured = (
                f"all {n_anchored:,} anchored reads were scored; this loss is "
                f"entirely upstream of the model (never anchored)"
            )
        cascade.loc[cascade["to_stage"] == "charge-called", "reason"] = measured

    # How much of the final gate's loss is structurally explained: a read whose
    # alignment never reaches the CCA cannot be anchored, whatever else is true.
    if anchor["scored"]:
        uncallable = anchor["ends_before_anchor"] + anchor["starts_after_anchor"]
        pct_uncallable = 100 * uncallable / anchor["scored"]
        drop = totals["aligned"] - totals["classified"]
        pct_drop = 100 * drop / totals["aligned"] if totals["aligned"] else 0.0
        if not n_anchored:
            cascade.loc[cascade["to_stage"] == "charge-called", "reason"] = (
                f"not emitted by the charging model; {pct_uncallable:.2f}% of "
                f"aligned reads cannot be anchored (do not span the CCA)"
            )
    else:
        pct_uncallable = pct_drop = None

    args.output.parent.mkdir(parents=True, exist_ok=True)
    cascade.to_csv(args.output, sep="\t", index=False)

    print("read attrition")
    for _, r in cascade.iterrows():
        print(
            f"  {r.from_stage:>18} -> {r.to_stage:<16} "
            f"lost {r.lost:>9,} ({r.pct_lost:5.2f}%)"
        )
    if sequenced is not None:
        gain = totals["basecalled"] - assigned
        if gain > 0:
            print(f"  basecalling emitted {gain:+,} sub-reads from split concatemers")
        elif gain < 0:
            print(
                f"  WARNING basecalling LOST {-gain:,} reads outright — check the "
                f"input POD5s for non-uniform signal batches (escapepod-rs#195)"
            )
    if pct_uncallable is not None:
        print(
            f"  of the {pct_drop:.2f}% dropped at charge-calling, {pct_uncallable:.2f}% "
            f"of aligned reads are structurally uncallable (no CCA in the alignment)"
        )
    if n_anchored:
        print(
            f"  charging model: {n_called:,} of {n_anchored:,} anchored reads "
            f"scored, {n_no_call:,} no-called "
            f"({100 * n_no_call / n_anchored:.2f}%)"
        )
        for reason, n in sorted(no_call.items(), key=lambda kv: -kv[1]):
            print(f"    {reason:<18} {n:>9,} ({100 * n / n_anchored:5.2f}%)")
        if n_no_call:
            print(
                "  NOTE no-calls are charging-correlated: report this rate beside "
                "any charging fraction, which is otherwise biased LOW"
            )


if __name__ == "__main__":
    main()
