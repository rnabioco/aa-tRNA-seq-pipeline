#!/usr/bin/env python
"""
Prove a WarpDemuX read -> barcode table now lives in the run's `.p5s` sidecar.

Run after `escpod sidecar annotate --from-table`, and before anything deletes
the table: the table is the only copy of the calls until this passes.

The check is by label counts, because that is what can be compared without
re-deriving escpod's read lookup: every label the table assigns (all but
`unclassified`) must appear in the sidecar exactly as many times as in the table,
and the sidecar must hold at least as many `unclassified` reads as the table
does (it also holds the reads WarpDemuX never saw, which escpod labels
`unclassified` too). A read that failed to land shows up as a short count.

A POD5 directory `<run>/pod5/` has its collection sidecar at `<run>/pod5.p5s`.

usage:
  verify_wdx_sidecar.py --table classifications.csv [--column barcode] DIR [DIR ...]
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

UNASSIGNED = "unclassified"


def sidecar_path(pod5_dir):
    return Path(str(pod5_dir).rstrip("/") + ".p5s")


def table_counts(path, column):
    with open(path, newline="") as fh:
        return Counter(row[column] for row in csv.DictReader(fh))


def sidecar_counts(paths, column):
    import pyarrow as pa

    counts = Counter()
    for p in paths:
        if not p.is_file():
            raise FileNotFoundError(f"no sidecar at {p}")
        with pa.memory_map(str(p)) as src:
            table = pa.ipc.open_file(src).read_all()
        if column not in table.column_names:
            raise KeyError(f"{p} has no {column!r} column (has {table.column_names})")
        counts.update(table.column(column).to_pylist())
    return counts


def discrepancies(expected, found):
    """Human-readable differences; empty when the sidecar carries the table."""
    problems = []
    for label, n in sorted(expected.items()):
        if label == UNASSIGNED:
            if found.get(label, 0) < n:
                problems.append(f"{label}: table {n}, sidecar {found.get(label, 0)}")
        elif found.get(label, 0) != n:
            problems.append(f"{label}: table {n}, sidecar {found.get(label, 0)}")
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--table", required=True, help="classifications table (CSV)")
    ap.add_argument("--column", default="barcode", help="label column, in both")
    ap.add_argument("dirs", nargs="+", help="POD5 directories that were annotated")
    args = ap.parse_args(argv)

    expected = table_counts(args.table, args.column)
    if not expected:
        print(f"{args.table} holds no calls", file=sys.stderr)
        return 1
    try:
        found = sidecar_counts([sidecar_path(d) for d in args.dirs], args.column)
    except (FileNotFoundError, KeyError) as e:
        print(f"sidecar not verified: {e}", file=sys.stderr)
        return 1

    problems = discrepancies(expected, found)
    if problems:
        print(
            "sidecar does not carry the table:", *problems, sep="\n  ", file=sys.stderr
        )
        return 1
    total = sum(expected.values())
    print(
        f"verified: {total} table rows, {len(expected)} labels present in the sidecar"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
