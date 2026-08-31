#!/usr/bin/env python3
"""Find a `--min-crf-margin` operating point for the WDX4 CTC-CRF bundle.

The concordance study left the CRF close to usable but not there: on the
held-out run it agreed with WarpDemuX on 98.85% of the teacher's confident
calls, with 0.51% cross-barcode -- against a gate of 99.0% and 0.20%. Its
failure is not uniform, though, and the CRF emits a per-read lattice score, so
the question is whether a confidence gate buys the remaining distance at an
acceptable cost in yield.

This is the same dial the LDX path already ships (`ldx.min_crf_margin: 1.0`,
measured there as removing 47.69% of all error at 84.56% discard precision).
Nothing is invented here: the operating point is DERIVED from this run, in the
way `config-base.yml` insists an operating point must be, rather than borrowed
from the LDX panel, which is a different model over a different oligo set.

What the columns mean:

    retained_pct   of escpod's calls that survive the gate. This is the yield
                   cost, and it is charged against ALL calls -- including the
                   reads WarpDemuX could not fingerprint, which is where the
                   yield gain lives.
    agreement      among WDX-confident reads still called by both.
    cross_pct_max  worst barcode's rate of being handed a read WarpDemuX
                   assigned to a DIFFERENT barcode. This is the only flow that
                   biases a charging fraction, so it is the number the gate is
                   really about.
    discard_prec   of the calls this gate discards, the fraction that were
                   disagreements with WarpDemuX. High is good: it means the
                   gate is removing error rather than burning yield.
"""

from __future__ import annotations

import argparse
import glob as globmod
import gzip
import sys

import pandas as pd

UNCLASSIFIED = "unclassified"


def eprint(*a):
    print(*a, file=sys.stderr)


def load_wdx_confident(pred_glob: str, min_conf: float) -> pd.DataFrame:
    frames = []
    for f in sorted(globmod.glob(pred_glob)):
        with gzip.open(f, "rt") as fh:
            df = pd.read_csv(fh)
        df.columns = [c.lstrip("#") for c in df.columns]
        frames.append(df[["read_id", "predicted_barcode", "confidence_score"]])
    pred = pd.concat(frames, ignore_index=True)
    pred = pred[
        (pred["predicted_barcode"] != -1) & (pred["confidence_score"] >= min_conf)
    ]
    pred["wdx"] = pred["predicted_barcode"].map(lambda n: f"barcode{int(n):02d}")
    eprint(f"wdx confident calls (>= {min_conf}): {len(pred):,}")
    return pred[["read_id", "wdx"]]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--crf-calls", required=True)
    ap.add_argument("--wdx-predictions", required=True)
    ap.add_argument("--min-wdx-conf", type=float, default=0.9)
    ap.add_argument("--score", default="crf_margin", help="column to gate on")
    ap.add_argument(
        "--grid",
        default="0,0.5,0.7,1.0,1.5,2.0,2.3,3.0,3.5,4.0,4.6,5.0,6.0",
    )
    ap.add_argument(
        "--scored-barcodes",
        help=(
            "Comma-separated barcodes that actually carry a sample on this "
            "run (e.g. barcode04,barcode05,barcode07). Only these set the "
            "headline cross-barcode rate: a barcode no sample claims is a "
            "discard channel, not a contamination route. Default: all."
        ),
    )
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    crf = pd.read_csv(args.crf_calls)
    labels = {f"bc{n:02d}": f"barcode{n:02d}" for n in (3, 4, 5, 7)}
    crf["escpod"] = crf["barcode"].map(lambda b: labels.get(b, UNCLASSIFIED))
    barcodes = sorted(labels.values())
    eprint(f"crf calls: {len(crf):,}")

    scored = set(args.scored_barcodes.split(",") if args.scored_barcodes else barcodes)
    eprint(f"barcodes carrying samples: {sorted(scored)}")

    wdx = load_wdx_confident(args.wdx_predictions, args.min_wdx_conf)
    j = crf.merge(wdx, on="read_id", how="left")
    j["wdx"] = j["wdx"].fillna("NA")

    n_called_all = int((j["escpod"] != UNCLASSIFIED).sum())
    rows = []
    for t in (float(x) for x in args.grid.split(",")):
        # A read below the gate becomes unclassified; it is not reassigned.
        passes = j[args.score] >= t
        called = (j["escpod"] != UNCLASSIFIED) & passes
        comparable = called & j["wdx"].isin(barcodes)
        same = comparable & (j["escpod"] == j["wdx"])

        cross_max, cross_bc = 0.0, None
        per_bc = {}
        for b in barcodes:
            got = called & (j["escpod"] == b)
            foreign = got & j["wdx"].isin([x for x in barcodes if x != b])
            pct = 100.0 * foreign.sum() / got.sum() if got.sum() else 0.0
            per_bc[f"cross_{b}"] = pct
            per_bc[f"n_{b}"] = int(got.sum())
            # A barcode no sample claims is a discard channel, not a
            # contamination route -- reads landing there reach no analysis. Its
            # rate is also computed on a tiny denominator and swings wildly
            # (run A: 0.198% -> 0.342% -> 0.653% -> 1.23% as the gate TIGHTENS),
            # so letting it set the headline would pick an operating point off
            # noise in a channel that does not matter.
            if b in scored and pct > cross_max:
                cross_max, cross_bc = pct, b

        # What the gate throws away, and whether it deserved to go.
        dropped = (j["escpod"] != UNCLASSIFIED) & ~passes
        drop_cmp = dropped & j["wdx"].isin(barcodes)
        drop_err = drop_cmp & (j["escpod"] != j["wdx"])

        rows.append(
            {
                args.score: t,
                "n_called": int(called.sum()),
                "retained_pct": 100.0 * called.sum() / n_called_all,
                "n_comparable": int(comparable.sum()),
                "agreement": float(same.sum() / comparable.sum())
                if comparable.sum()
                else float("nan"),
                "cross_pct_max": cross_max,
                "cross_worst_barcode": cross_bc,
                "n_dropped": int(dropped.sum()),
                "discard_prec": float(drop_err.sum() / drop_cmp.sum())
                if drop_cmp.sum()
                else float("nan"),
                **per_bc,
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(out.to_string(index=False))

    ok = out[(out["agreement"] >= 0.99) & (out["cross_pct_max"] <= 0.20)]
    print()
    if ok.empty:
        print(
            "NO operating point on this grid reaches agreement>=99.0% and cross<=0.20%."
        )
        best = out.loc[out["cross_pct_max"].idxmin()]
        print(
            f"closest on cross: {args.score}>={best[args.score]} -> "
            f"cross {best['cross_pct_max']:.3f}%, agreement {best['agreement']:.4%}, "
            f"retaining {best['retained_pct']:.1f}% of calls"
        )
    else:
        b = ok.iloc[0]
        print(
            f"USABLE at {args.score} >= {b[args.score]}: agreement "
            f"{b['agreement']:.4%}, cross {b['cross_pct_max']:.3f}%, "
            f"retaining {b['retained_pct']:.1f}% of calls "
            f"({b['n_called']:,} reads)"
        )


if __name__ == "__main__":
    sys.exit(main())
