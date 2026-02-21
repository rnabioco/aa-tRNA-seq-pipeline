#!/usr/bin/env python3
"""
Build concordance table of WDX sample assignment vs EDX (3' adapter barcode) identity.

Reads adapter detection TSVs (from detect_3p_adapters.py) that contain all reads
with their detected 3' adapter identity, then tabulates counts per WDX sample.
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
        help="Adapter detection TSV files (one per WDX sample, gzipped)",
    )
    parser.add_argument(
        "--samples",
        nargs="+",
        required=True,
        help="Sample names corresponding to TSV files (same order)",
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

    for tsv_path, sample_name in zip(args.tsvs, args.samples):
        counts = Counter()

        print(f"Processing {sample_name}: {tsv_path}", file=sys.stderr)

        with gzip.open(tsv_path, "rt") as f:
            f.readline()  # skip header
            for line in f:
                parts = line.rstrip("\n").split("\t", 1)
                if len(parts) == 2:
                    adapter = parts[1]
                else:
                    adapter = "unknown"
                counts[adapter] += 1

        total = sum(counts.values())
        print(
            f"  {sample_name}: {total} reads, "
            f"{len(counts)} adapter categories",
            file=sys.stderr,
        )

        for adapter, n in sorted(counts.items()):
            pct = 100.0 * n / total if total > 0 else 0.0
            rows.append((sample_name, adapter, n, f"{pct:.1f}"))

    with gzip.open(args.output, "wt") as f:
        f.write("sample\tedx_adapter\tn_reads\tpct\n")
        for sample, adapter, n, pct in rows:
            f.write(f"{sample}\t{adapter}\t{n}\t{pct}\n")

    print(f"Wrote {len(rows)} rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
