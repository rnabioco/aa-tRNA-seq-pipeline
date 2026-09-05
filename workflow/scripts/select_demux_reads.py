#!/usr/bin/env python
"""
Select the reads demux assigned to one or more samples, joining barcode axes.

A sample is a tuple of barcode calls, one per axis it names: `ldx=ldx01` alone
on a single-index library, `ldx=ldx01,fdx=fdx01` on a dual-index one (a 3' LDX
code and a 5' FDX code). A read belongs to a sample when EVERY axis the sample
names agrees with the call for that read.

The calls come from `escpod demux --classifications` in either of its shapes:

  * one CSV per axis, from one pass per model, columns
    `read_id,barcode,confidence,crf_logp,crf_margin,...`
  * one CSV for all axes, from a fused `--model ldx=... --model fdx=...` pass,
    columns `read_id,ldx,ldx_confidence,ldx_crf_margin,...,fdx,fdx_confidence,...`

An axis's call is read from the column named after the axis when the CSV has
one, else from `barcode`; its lattice margin from `<axis>_crf_margin`, else
`crf_margin`. So the two shapes are joined by the same code, and the same CSV
may be given for several axes.

Gates. `--gate AXIS=NATS` calls a read `unclassified` on that axis when its
margin is below NATS (or missing). escpod applies one `--min-crf-margin` to
every model of a fused pass, while the axes want different operating points --
the LDX bundle's 1.0 against the FDX bundle's declared 3.5 -- so the per-axis
gate lives here, where it can differ. Re-applying a gate escpod already applied
is a no-op, which is what makes the two shapes interchangeable.

Memory is bounded by the reads the FIRST axis assigns to any configured code,
not by the run: its CSV is scanned once keeping only wanted codes, and each
further axis is scanned once keeping only those reads. Samples sharing a
first-axis code are indexed by it, so per-read work does not grow with the
number of samples.

Split children: dorado splits a concatenated read into children with NEW ids,
recorded in `pi:Z`. Demux only ever saw the parent, so a child inherits its
parent's assignment when a parent map (`--parents`, `child<TAB>parent`) is given.

usage:
  select_demux_reads.py --axis ldx=run/classifications.csv \\
                        --axis fdx=run/fdx/classifications.csv \\
                        --gate fdx=3.5 \\
                        --sample fdx01_ldx01:ldx=ldx01,fdx=fdx01 \\
                        --sample fdx02_ldx04:ldx=ldx04,fdx=fdx02 \\
                        [--parents run/split_parents.tsv] \\
                        --output read_ids.txt [--summary summary.tsv]

`unclassified` never matches a configured code, so it needs no special handling.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict


class SelectionError(ValueError):
    """The arguments do not describe a selection that can be made."""


def parse_axes(items):
    """`NAME=PATH` items -> [(name, path), ...] in the order given; first is primary."""
    axes = []
    for item in items:
        name, sep, path = item.partition("=")
        if not sep or not name or not path:
            raise SelectionError(f"--axis expects NAME=PATH, got {item!r}")
        if name in dict(axes):
            raise SelectionError(f"--axis {name!r} given twice")
        axes.append((name, path))
    if not axes:
        raise SelectionError("at least one --axis is required")
    return axes


def parse_gates(items, axes):
    """`AXIS=NATS` items -> {axis: float}."""
    names = [name for name, _ in axes]
    gates = {}
    for item in items:
        axis, sep, value = item.partition("=")
        if not sep or axis not in names:
            raise SelectionError(
                f"--gate expects AXIS=NATS for a declared axis, got {item!r}"
            )
        try:
            gates[axis] = float(value)
        except ValueError:
            raise SelectionError(f"--gate {axis}: {value!r} is not a number") from None
    return gates


def parse_samples(items, axes):
    """`NAME:axis=code,axis=code` items -> {name: {axis: code}}.

    Every sample must name the primary (first) axis: that is the scan the whole
    selection is bounded by, and a sample without it would have to be matched
    against every read of the run.
    """
    axis_names = [name for name, _ in axes]
    primary = axis_names[0]
    samples = {}
    for item in items:
        name, sep, spec = item.partition(":")
        if not sep or not name or not spec:
            raise SelectionError(
                f"--sample expects NAME:axis=code[,axis=code], got {item!r}"
            )
        if name in samples:
            raise SelectionError(f"--sample {name!r} given twice")
        codes = {}
        for pair in spec.split(","):
            axis, sep, code = pair.partition("=")
            if not sep or not axis or not code:
                raise SelectionError(f"--sample {name!r}: bad axis=code {pair!r}")
            if axis not in axis_names:
                raise SelectionError(
                    f"--sample {name!r} names axis {axis!r}, but no --axis {axis}=... was given"
                )
            if axis in codes:
                raise SelectionError(f"--sample {name!r} names axis {axis!r} twice")
            codes[axis] = code
        if primary not in codes:
            raise SelectionError(
                f"--sample {name!r} does not name the primary axis {primary!r}"
            )
        samples[name] = codes
    if not samples:
        raise SelectionError("at least one --sample is required")

    # Two samples with identical tuples would each receive the same reads; a
    # sample whose tuple is a strict prefix of another's (ldx01 against
    # ldx01+fdx01) would swallow the other's reads. Both are config errors and
    # both are cheaper to refuse here than to discover in a BAM.
    by_primary = defaultdict(list)
    for name, codes in samples.items():
        by_primary[codes[primary]].append((name, codes))
    for code, group in by_primary.items():
        axis_sets = {tuple(sorted(codes)) for _, codes in group}
        if len(axis_sets) > 1:
            names = ", ".join(sorted(n for n, _ in group))
            raise SelectionError(
                f"samples {names} share {primary}={code} but name different axis "
                f"sets; every sample sharing a {primary} code must name the same axes"
            )
        tuples = [tuple(sorted(codes.items())) for _, codes in group]
        if len(set(tuples)) != len(tuples):
            names = ", ".join(sorted(n for n, _ in group))
            raise SelectionError(f"samples {names} have identical barcode tuples")
    return samples


def _columns_for(axis, fieldnames, path):
    """(call column, margin column or None) for `axis` in a CSV with `fieldnames`."""
    if axis in fieldnames:
        call = axis
    elif "barcode" in fieldnames:
        call = "barcode"
    else:
        raise SelectionError(
            f"{path}: neither a {axis!r} nor a `barcode` column (columns: "
            f"{fieldnames}); this is not an escpod classifications CSV"
        )
    for margin in (f"{axis}_crf_margin", "crf_margin"):
        if margin in fieldnames:
            return call, margin
    return call, None


def _passes(row, margin_col, gate):
    if gate is None:
        return True
    if margin_col is None:
        raise SelectionError(
            "a --gate was given but the CSV carries no crf_margin column; "
            "run escpod demux with --ref-scores"
        )
    try:
        return float(row[margin_col]) >= gate
    except (TypeError, ValueError):
        return False  # unclassified rows carry no margin


def read_calls(path, axis, wanted, gate=None, restrict=None):
    """read_id -> code, for rows whose (gated) call on `axis` is in `wanted`.

    With `restrict`, only reads in that set are kept, so a secondary axis costs
    memory proportional to the primary axis's candidates rather than to the run.
    """
    calls = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise SelectionError(f"{path}: empty file")
        call_col, margin_col = _columns_for(axis, reader.fieldnames, path)
        for row in reader:
            code = row[call_col]
            if code not in wanted:
                continue
            read_id = row["read_id"]
            if restrict is not None and read_id not in restrict:
                continue
            if not _passes(row, margin_col, gate):
                continue
            calls[read_id] = code
    return calls


def select_reads(axes, samples, gates=None):
    """{sample: set(read_id)} for the reads every named axis agrees on."""
    gates = gates or {}
    primary, primary_path = axes[0]
    wanted_primary = {codes[primary] for codes in samples.values()}
    primary_calls = read_calls(
        primary_path, primary, wanted_primary, gates.get(primary)
    )
    candidates = set(primary_calls)

    secondary_calls = {}
    for axis, path in axes[1:]:
        wanted = {codes[axis] for codes in samples.values() if axis in codes}
        if not wanted:
            continue
        secondary_calls[axis] = read_calls(
            path, axis, wanted, gates.get(axis), restrict=candidates
        )

    by_primary = defaultdict(list)
    for name, codes in samples.items():
        by_primary[codes[primary]].append((name, codes))

    assigned = {name: set() for name in samples}
    for read_id, code in primary_calls.items():
        for name, codes in by_primary[code]:
            if all(
                secondary_calls.get(axis, {}).get(read_id) == want
                for axis, want in codes.items()
                if axis != primary
            ):
                assigned[name].add(read_id)
    return assigned


def add_split_children(assigned, parents_path):
    """Give each split child its parent's assignment. Returns children added per sample."""
    parent_of_sample = {}
    for name, reads in assigned.items():
        for read_id in reads:
            parent_of_sample[read_id] = name
    added = defaultdict(int)
    with open(parents_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            child, parent = line.rstrip("\n").split("\t")
            name = parent_of_sample.get(parent)
            if name is not None:
                assigned[name].add(child)
                added[name] += 1
    return added


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--axis",
        action="append",
        default=[],
        metavar="NAME=CSV",
        help="A barcode axis and its classifications CSV. Repeatable; the first is primary.",
    )
    parser.add_argument(
        "--gate",
        action="append",
        default=[],
        metavar="AXIS=NATS",
        help="Call a read unclassified on AXIS when its crf_margin is below NATS.",
    )
    parser.add_argument(
        "--sample",
        action="append",
        default=[],
        metavar="NAME:axis=code[,axis=code]",
        help="A sample and the code it carries on each axis it names. Repeatable.",
    )
    parser.add_argument(
        "--parents",
        help="split_parents.tsv (child<TAB>parent): children inherit their parent's sample",
    )
    parser.add_argument(
        "--output", required=True, help="Read ids, one per line, sorted"
    )
    parser.add_argument(
        "--summary", help="Per-sample counts as TSV: sample, assigned, split_children"
    )
    args = parser.parse_args(argv)

    try:
        axes = parse_axes(args.axis)
        gates = parse_gates(args.gate, axes)
        samples = parse_samples(args.sample, axes)
        assigned = select_reads(axes, samples, gates)
    except SelectionError as exc:
        parser.error(str(exc))

    children = add_split_children(assigned, args.parents) if args.parents else {}

    union = set().union(*assigned.values())
    with open(args.output, "w") as fh:
        fh.writelines(f"{read_id}\n" for read_id in sorted(union))
    if args.summary:
        with open(args.summary, "w") as fh:
            fh.write("sample\tassigned\tsplit_children\n")
            for name in samples:
                n_children = children.get(name, 0)
                fh.write(f"{name}\t{len(assigned[name]) - n_children}\t{n_children}\n")

    for name, codes in samples.items():
        tuple_str = ",".join(f"{a}={c}" for a, c in codes.items())
        print(
            f"{name} ({tuple_str}): {len(assigned[name]) - children.get(name, 0)} "
            f"assigned reads, {children.get(name, 0)} split children",
            file=sys.stderr,
        )
    empty = [name for name, reads in assigned.items() if not reads]
    if empty:
        axes_str = ", ".join(f"{a} <- {p}" for a, p in axes)
        sys.exit(
            f"ERROR: no reads were assigned to sample(s) {', '.join(empty)}.\n"
            f"Check the per-barcode counts in the demux_summary.tsv.gz beside each "
            f"axis's classifications ({axes_str})."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
