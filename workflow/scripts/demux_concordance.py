#!/usr/bin/env python3
"""
Compare WarpDemuX's barcode routing against the escpod CRF's, read by read.

This is the measurement that has to pass before WarpDemuX can be retired: the
two backends do not have to agree perfectly, but the ways they disagree have to
be accounted for, because which reads land in which sample is a scientific
property of the run and not an implementation detail.

WHAT AGREEMENT IS AND IS NOT WORTH
----------------------------------
The CRF was trained on WarpDemuX's own calls, so it is a student and WarpDemuX
is its teacher. High agreement is the EXPECTED outcome of a successful
distillation and is evidence about the training, not about which tool is right.
Nothing in this script can establish correctness. What it can do is quantify
the disagreement and split it by cause, so the independent arbiters (3' adapter
purity, tRNA alignment rate, charging fraction) have a defined population to run
on. `--charging` computes the last of those here, since it is the acceptance
criterion and the input is already per-read.

WHY THE DISAGREEMENT IS SPLIT
-----------------------------
Two very different things look like "the CRF disagreed":

  refused   the CRF never decoded the read -- its adapter ended too early for
            the model's window, so it routed to `unclassified` with confidence
            exactly 0. This is a yield question, tunable with --boundary-margin
            and --clamp-max-shift, and says nothing about the classifier.
  differed  the CRF decoded the read and matched it to a different barcode.
            This is a genuine classification disagreement.

Pooling them hides the entire boundary-window question inside a single
concordance percentage, which is how a tunable 14% yield loss stayed invisible
on the LDX panel until someone went looking.

RECOVERING WARPDEMUX'S CALLS
----------------------------
`demux/read_ids/<run>/barcode_mapping.tsv.gz` is the direct record, but it is
`demux_scratch` tier and has been cleaned on most completed runs. The fallback
is the per-sample charging tables: a read earns a row in
`summary/tables/<sample>/<sample>.charging_prob.tsv.gz` only if WarpDemuX routed
it to that sample, so sample membership IS the barcode call. That is how the
training corpus was recovered upstream, and it works on every completed run.

Its one limitation, which matters when reading the output: a read must also have
ALIGNED to a tRNA to earn a row, so this view cannot see reads WarpDemuX called
but that failed alignment, and cannot see WarpDemuX's own `unclassified` at all.
Rows attributed to WarpDemuX are therefore "called AND aligned"; use
--mapping for the unrestricted comparison when the file still exists.
"""

import argparse
import csv
import gzip
import sys
from collections import Counter, defaultdict


UNCLASSIFIED = "unclassified"


def read_crf_calls(path):
    """read_id -> (barcode, confidence) from escpod's classifications CSV."""
    calls = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            read_id = row.get("read_id") or row.get("#read_id")
            barcode = row.get("barcode", UNCLASSIFIED)
            try:
                confidence = float(row.get("confidence", "nan"))
            except ValueError:
                confidence = float("nan")
            calls[read_id] = (barcode, confidence)
    return calls


def read_wdx_from_mapping(path):
    """read_id -> barcode from WarpDemuX's barcode_mapping.tsv.gz."""
    calls = {}
    with gzip.open(path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            calls[row["read_id"]] = row["predicted_barcode"]
    return calls


def read_wdx_from_tables(specs):
    """read_id -> (barcode, charging_likelihood, trna) from charging tables.

    `specs` is a list of "barcode=path" strings. A read appearing under two
    barcodes would mean the tables disagree, which cannot happen for a single
    run -- report it rather than silently keeping one.
    """
    calls, extra, clashes = {}, {}, 0
    for spec in specs:
        barcode, _, path = spec.partition("=")
        if not path:
            sys.exit(f"--tables entries must be BARCODE=PATH, got: {spec!r}")
        with gzip.open(path, "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                read_id = row["read_id"]
                if read_id in calls and calls[read_id] != barcode:
                    clashes += 1
                    continue
                calls[read_id] = barcode
                cl = row.get("charging_likelihood")
                extra[read_id] = (
                    int(cl) if cl not in (None, "", "NA") else None,
                    row.get("tRNA"),
                )
    if clashes:
        print(
            f"warning: {clashes} reads appeared under more than one barcode; "
            "kept the first. Are these tables all from one run?",
            file=sys.stderr,
        )
    return calls, extra


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--classifications",
        required=True,
        help="escpod demux --classifications CSV (the CRF's calls)",
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--mapping",
        help="WarpDemuX barcode_mapping.tsv.gz (direct, includes unclassified)",
    )
    src.add_argument(
        "--tables",
        nargs="+",
        help=(
            "BARCODE=path/to/<sample>.charging_prob.tsv.gz, one per sample. "
            "Fallback for runs whose mapping was cleaned; restricted to reads "
            "that aligned."
        ),
    )
    parser.add_argument(
        "--ml-threshold",
        type=int,
        default=200,
        help="charging_likelihood at or above this is charged (default: 200)",
    )
    parser.add_argument(
        "--output",
        help="Write the contingency table here as TSV (default: stdout summary only)",
    )
    args = parser.parse_args()

    crf = read_crf_calls(args.classifications)
    charging = {}
    if args.mapping:
        wdx = read_wdx_from_mapping(args.mapping)
        wdx_source = "barcode_mapping.tsv.gz (all called reads)"
    else:
        wdx, charging = read_wdx_from_tables(args.tables)
        wdx_source = "charging tables (called AND aligned reads only)"

    shared = wdx.keys() & crf.keys()

    # The contingency table, and the three-way split of every shared read.
    table = Counter()
    buckets = defaultdict(list)
    for read_id in shared:
        w = wdx[read_id]
        c, conf = crf[read_id]
        table[(w, c)] += 1
        if c == w:
            bucket = "agree"
        elif c == UNCLASSIFIED and conf == 0.0:
            # Never decoded: the boundary window, not the classifier.
            bucket = "refused"
        elif c == UNCLASSIFIED:
            bucket = "unclassified_decoded"
        else:
            bucket = "differed"
        buckets[bucket].append(read_id)

    total = len(shared)
    print(f"WarpDemuX source : {wdx_source}")
    print(f"reads in WarpDemuX: {len(wdx)}")
    print(f"reads in CRF      : {len(crf)}")
    print(f"reads in both     : {total}")
    if not total:
        sys.exit("No reads in common — are these from the same run?")
    print()
    print(f"{'outcome':22} {'reads':>10} {'pct':>7}")
    for name in ("agree", "refused", "unclassified_decoded", "differed"):
        n = len(buckets[name])
        print(f"{name:22} {n:>10} {100 * n / total:>6.2f}%")
    print()
    print(
        "NOTE: `agree` measures student-teacher agreement, not accuracy — the\n"
        "      CRF was trained on these very calls. `refused` is a yield knob\n"
        "      (--boundary-margin / --clamp-max-shift), not a classifier error."
    )

    # The acceptance criterion: if the disagreement sets carry the same charging
    # distribution as the agreement set, the swap cannot move a conclusion.
    if charging:
        print()
        print(f"{'outcome':22} {'n_scored':>10} {'charged':>9} {'frac':>7}")
        for name in ("agree", "refused", "unclassified_decoded", "differed"):
            vals = [
                charging[r][0]
                for r in buckets[name]
                if charging.get(r) and charging[r][0] is not None
            ]
            if not vals:
                print(f"{name:22} {0:>10} {'-':>9} {'-':>7}")
                continue
            charged = sum(1 for v in vals if v >= args.ml_threshold)
            print(f"{name:22} {len(vals):>10} {charged:>9} {charged / len(vals):>6.3f}")
        print()
        print(
            "      Charging fraction across the buckets is the acceptance\n"
            "      criterion: if the disagreement sets match the agreement set,\n"
            "      the routing change cannot move a biological conclusion."
        )

    if args.output:
        opener = gzip.open if args.output.endswith(".gz") else open
        with opener(args.output, "wt", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["wdx_barcode", "crf_barcode", "n_reads"])
            for (wb, cb), n in sorted(table.items(), key=lambda kv: -kv[1]):
                w.writerow([wb, cb, n])
        print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
