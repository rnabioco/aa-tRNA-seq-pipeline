#!/usr/bin/env python3
"""
Build concordance table of WDX sample assignment vs EDX (3' adapter barcode) identity.

Reads PT tags from final BAMs to determine which 3' adapter each read matched,
then tabulates counts per WDX sample.
"""

import argparse
import gzip
import re
import sys

import pysam


def parse_3p_adapter_from_pt(pt_tag):
    """
    Extract 3' adapter name from a PT tag string.

    PT tag format: "start;end;strand;type|start;end;strand;type"
    3' adapter entries look like:
      - "3p_adapter"       -> "default"
      - "3p_adapter_edx1"  -> "edx1"
      - "3p_adapter_v2"    -> "v2"

    Returns the adapter name or None if no 3' adapter found.
    """
    if not pt_tag:
        return None

    for segment in pt_tag.split("|"):
        fields = segment.split(";")
        if len(fields) < 4:
            continue
        entry_type = fields[3]
        if entry_type.startswith("3p_adapter"):
            suffix = entry_type[len("3p_adapter"):]
            if suffix.startswith("_"):
                return suffix[1:]  # e.g., "edx1", "v2"
            elif suffix == "":
                return "default"
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bams",
        nargs="+",
        required=True,
        help="Final BAM files (one per WDX sample)",
    )
    parser.add_argument(
        "--samples",
        nargs="+",
        required=True,
        help="Sample names corresponding to BAM files (same order)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file (gzipped)",
    )
    args = parser.parse_args()

    if len(args.bams) != len(args.samples):
        sys.exit(
            f"Number of BAMs ({len(args.bams)}) must match "
            f"number of samples ({len(args.samples)})"
        )

    rows = []

    for bam_path, sample_name in zip(args.bams, args.samples):
        counts = {}
        total = 0

        print(f"Processing {sample_name}: {bam_path}", file=sys.stderr)

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                total += 1
                try:
                    pt_tag = read.get_tag("PT")
                except KeyError:
                    pt_tag = None

                adapter = parse_3p_adapter_from_pt(pt_tag)
                if adapter is None:
                    adapter = "no_3p_adapter"

                counts[adapter] = counts.get(adapter, 0) + 1

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
