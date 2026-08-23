#!/usr/bin/env python3
"""
Build a contingency table of signal-barcode sample vs EDX (3' adapter) identity.

Reads adapter detection TSVs (from detect_3p_adapters.py), which contain every
read with its detected 3' adapter, and tabulates counts per sample.

For a sample whose library carries a single adapter this is a purity check. For
an adapter-ligation QC run, where one signal barcode deliberately holds a
mixture of adapters, the same table is the measurement itself.

Calls are split into CLEAN and ambiguous by the alignment-score margin between
the best and second-best adapter. Without that split, a read that scored 67 for
one adapter and 66 for another is counted as an adapter-swap event
indistinguishable from a read that matched at 66 with the runner-up at 20 —
which inflates apparent crosstalk with reads the assay simply could not assign.
Both are reported so the gate's effect stays visible: `n_reads` is every call,
`n_clean` only those at or above the margin.
"""

import argparse
import gzip
import sys
from collections import Counter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tsvs",
        nargs="+",
        required=True,
        help="Adapter detection TSV files (one per sample, gzipped)",
    )
    parser.add_argument(
        "--samples",
        nargs="+",
        required=True,
        help="Sample names corresponding to TSV files (same order)",
    )
    parser.add_argument(
        "--min-margin",
        type=int,
        default=5,
        help=(
            "Minimum best-vs-second-best score margin for a call to count as "
            "CLEAN (default: 5; 0 disables the gate)"
        ),
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file (gzipped)",
    )
    args = parser.parse_args()

    if len(args.tsvs) != len(args.samples):
        sys.exit(
            f"Number of TSVs ({len(args.tsvs)}) must match "
            f"number of samples ({len(args.samples)})"
        )

    rows = []

    for tsv_path, sample_name in zip(args.tsvs, args.samples, strict=True):
        counts = Counter()
        clean_counts = Counter()

        print(f"Processing {sample_name}: {tsv_path}", file=sys.stderr)

        with gzip.open(tsv_path, "rt") as f:
            f.readline()  # skip header
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    continue
                adapter = fields[1] or "unknown"
                counts[adapter] += 1

                # An empty margin means there was no second-best to compare
                # against (a single configured adapter, or no match at all), so
                # the gate cannot apply and the call passes unchanged.
                margin_str = fields[4] if len(fields) > 4 else ""
                if adapter == "none":
                    continue
                if not margin_str or int(margin_str) >= args.min_margin:
                    clean_counts[adapter] += 1

        total = sum(counts.values())
        total_clean = sum(clean_counts.values())
        print(
            f"  {sample_name}: {total} reads, {len(counts)} adapter categories, "
            f"{total_clean} clean at margin >= {args.min_margin}",
            file=sys.stderr,
        )

        for adapter, n in sorted(counts.items()):
            pct = 100.0 * n / total if total > 0 else 0.0
            n_clean = clean_counts[adapter]
            # Percentages are within-sample, so each sample's column sums to
            # 100 and samples of different depth stay comparable.
            pct_clean = 100.0 * n_clean / total_clean if total_clean > 0 else 0.0
            rows.append(
                (sample_name, adapter, n, f"{pct:.1f}", n_clean, f"{pct_clean:.1f}")
            )

    with gzip.open(args.output, "wt") as f:
        f.write("sample\tedx_adapter\tn_reads\tpct\tn_clean\tpct_clean\n")
        for sample, adapter, n, pct, n_clean, pct_clean in rows:
            f.write(f"{sample}\t{adapter}\t{n}\t{pct}\t{n_clean}\t{pct_clean}\n")

    print(f"Wrote {len(rows)} rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
