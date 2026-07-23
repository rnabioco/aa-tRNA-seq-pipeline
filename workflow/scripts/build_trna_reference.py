#!/usr/bin/env python
"""
Validate or build tRNA reference FASTA with adapters.

This script ensures that tRNA reference sequences have the proper structure
required for the aa-tRNA-seq pipeline:
  1. 5' adapter prefix
  2. tRNA sequence ending with CCA
  3. 3' adapter suffix (must start with GGC to form CCAGGC junction)

The CCAGGC 6-mer junction (CCA from tRNA + GGC from adapter) is critical
for the charging classification model, which analyzes the nanopore
signal over this region to distinguish charged vs uncharged tRNAs.

Modes:
  validate: Check existing adapted reference has correct structure
  build: Create adapted reference from raw tRNA sequences (adds CCA if missing)

Both modes deduplicate sequences with identical content, keeping the first
occurrence and reporting collapsed duplicates in the validation/build report.

Requirements:
  - 3' adapter must start with GGC (for CCAGGC junction)
  - No duplicate sequence names
  - In build mode: CCA will be added if missing (with warning)
"""

import argparse
import sys
from collections import defaultdict


def normalize_sequence(seq):
    """Uppercase and convert RNA (U) to DNA (T).

    Reference sources such as GtRNAdb provide RNA-alphabet mature tRNA
    sequences; downstream tooling (BWA index/align) requires DNA, so U->T
    normalization happens on read.
    """
    return seq.upper().replace("U", "T")


def read_fasta(fasta_path):
    """
    Read FASTA file and yield (name, sequence) tuples.
    Handles multi-line sequences. Sequences are normalized to uppercase DNA.
    """
    name = None
    seq_parts = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, normalize_sequence("".join(seq_parts))
                name = line[1:].split()[0]  # Get first word after >
                seq_parts = []
            else:
                seq_parts.append(line)

        if name is not None:
            yield name, normalize_sequence("".join(seq_parts))


def write_fasta(sequences, output_path):
    """
    Write sequences to FASTA format (single line per sequence).
    """
    with open(output_path, "w") as f:
        for name, seq in sequences:
            f.write(f">{name}\n")
            f.write(seq + "\n")


def adapter_matches(expected, actual):
    """
    Check if actual sequence matches expected adapter.
    Handles N as wildcard (matches any base).
    """
    if len(expected) != len(actual):
        return False
    # strict=True is free here: the length guard above already returned.
    return all(e == "N" or e == a for e, a in zip(expected, actual, strict=True))


def deduplicate_sequences(sequences):
    """
    Remove sequences with identical content, keeping the first occurrence.

    Groups sequences by their actual nucleotide content and retains one
    representative per group. Returns the deduplicated list plus a mapping
    of kept name -> list of duplicate names that were collapsed into it.

    Args:
        sequences: List of (name, sequence) tuples.

    Returns:
        Tuple of (deduped_sequences, collapsed_map) where:
          - deduped_sequences: List of (name, seq) with duplicates removed
          - collapsed_map: Dict mapping kept_name -> [dropped_name, ...]
    """
    seen = {}  # seq -> name (first occurrence)
    collapsed = {}  # kept_name -> [dropped_names]
    deduped = []

    for name, seq in sequences:
        if seq in seen:
            kept_name = seen[seq]
            collapsed.setdefault(kept_name, []).append(name)
        else:
            seen[seq] = name
            deduped.append((name, seq))

    return deduped, collapsed


def _within_mismatches(seq_a, seq_b, max_mismatch):
    """Whether two equal-length sequences differ by at most max_mismatch bases.

    Counts incrementally and gives up as soon as the budget is exceeded, which
    is what keeps the all-against-leaders scan cheap: unrelated tRNAs diverge
    within the first few bases.
    """
    mismatches = 0
    for a, b in zip(seq_a, seq_b):
        if a != b:
            mismatches += 1
            if mismatches > max_mismatch:
                return False
    return True


def collapse_similar_sequences(keys, max_mismatch):
    """
    Merge near-identical sequences by Hamming distance.

    Uses greedy leader clustering: each sequence joins the first leader it is
    within max_mismatch substitutions of, otherwise becomes a new leader. Unlike
    single-linkage this bounds the cluster radius, so every member is within
    max_mismatch of its leader rather than merely chained to it.

    Hamming distance is only defined for equal-length sequences, so molecules
    differing by an indel are never merged regardless of max_mismatch.

    Args:
        keys: List of (name, sequence) tuples to cluster on.
        max_mismatch: Maximum substitutions between a member and its leader.

    Returns:
        Tuple of (kept_names, collapsed_map) where:
          - kept_names: List of leader names, in input order
          - collapsed_map: Dict mapping leader_name -> [merged_name, ...]
    """
    if max_mismatch <= 0:
        return [name for name, _ in keys], {}

    # Only equal-length sequences are comparable, so bucket by length first
    by_length = defaultdict(list)
    for name, seq in keys:
        by_length[len(seq)].append((name, seq))

    leaders = {}  # name -> sequence, for leaders only
    collapsed = {}
    order = []

    for bucket in by_length.values():
        bucket_leaders = []  # (name, seq)
        for name, seq in bucket:
            for leader_name, leader_seq in bucket_leaders:
                if _within_mismatches(seq, leader_seq, max_mismatch):
                    collapsed.setdefault(leader_name, []).append(name)
                    break
            else:
                bucket_leaders.append((name, seq))
                leaders[name] = seq

    # Preserve the original input ordering of the leaders
    for name, _ in keys:
        if name in leaders:
            order.append(name)

    return order, collapsed


def validate_reference(input_fasta, output_fasta, report_path, adapter_5p, adapters_3p):
    """
    Validate an existing adapted reference FASTA.

    Checks:
    1. All sequences have correct 5' adapter (allowing N wildcards)
    2. All tRNA portions end with CCA
    3. All sequences have correct 3' adapter matching any provided adapter (must start with GGC)
    4. No duplicate sequence names

    Args:
        adapters_3p: List of 3' adapter sequences. A sequence passes if it
            matches any adapter in the list.
    """
    errors = []
    warnings = []
    stats = defaultdict(int)
    validated_sequences = []
    seen_names = set()

    len_5p = len(adapter_5p)

    # All adapters must be the same length and start with GGC
    adapter_lengths = {len(a) for a in adapters_3p}
    if len(adapter_lengths) != 1:
        errors.append(
            f"All 3' adapters must have the same length. Got lengths: {sorted(adapter_lengths)}"
        )
        # Use the first adapter's length as fallback
        len_3p = len(adapters_3p[0])
    else:
        len_3p = adapter_lengths.pop()

    for adapter_3p in adapters_3p:
        if not adapter_3p.startswith("GGC"):
            errors.append(
                f"3' adapter must start with GGC for CCAGGC junction. Got: {adapter_3p[:3]}"
            )

    for name, seq in read_fasta(input_fasta):
        stats["total_sequences"] += 1

        # Check for duplicates
        if name in seen_names:
            errors.append(f"Duplicate sequence name: {name}")
            stats["duplicate_names"] += 1
            continue
        seen_names.add(name)

        # Check minimum length
        min_len = len_5p + 3 + len_3p  # adapter + CCA + adapter
        if len(seq) < min_len:
            errors.append(
                f"{name}: Sequence too short ({len(seq)} bp). "
                f"Minimum expected: {min_len} bp"
            )
            stats["too_short"] += 1
            continue

        # Extract parts
        actual_5p = seq[:len_5p]
        actual_3p = seq[-len_3p:] if len_3p > 0 else ""
        trna_seq = seq[len_5p:-len_3p] if len_3p > 0 else seq[len_5p:]

        # Check 5' adapter
        if not adapter_matches(adapter_5p, actual_5p):
            errors.append(
                f"{name}: Invalid 5' adapter.\n"
                f"    Expected: {adapter_5p}\n"
                f"    Got:      {actual_5p}"
            )
            stats["invalid_5p_adapter"] += 1
        else:
            stats["valid_5p_adapter"] += 1

        # Check 3' adapter - passes if it matches ANY of the provided adapters
        matched_adapter = None
        for adapter_3p in adapters_3p:
            if adapter_matches(adapter_3p, actual_3p):
                matched_adapter = adapter_3p
                break

        if matched_adapter is None:
            errors.append(
                f"{name}: Invalid 3' adapter.\n"
                f"    Expected one of: {adapters_3p}\n"
                f"    Got:             {actual_3p}"
            )
            stats["invalid_3p_adapter"] += 1
        else:
            stats["valid_3p_adapter"] += 1

        # Check CCA ending (tRNA portion should end with CCA)
        if not trna_seq.endswith("CCA"):
            errors.append(
                f"{name}: tRNA portion does not end with CCA. "
                f"Ends with: '{trna_seq[-3:] if len(trna_seq) >= 3 else trna_seq}'"
            )
            stats["missing_cca"] += 1
        else:
            stats["valid_cca"] += 1

        # Check that CCAGGC junction is present
        junction_start = len_5p + len(trna_seq) - 3  # Position of CCA
        junction = seq[junction_start : junction_start + 6]
        if junction != "CCAGGC":
            warnings.append(
                f"{name}: CCAGGC junction not found. Got: {junction}. "
                "Charging classification may not work correctly."
            )

        validated_sequences.append((name, seq))

    # Deduplicate sequences with identical content
    deduped_sequences, collapsed_map = deduplicate_sequences(validated_sequences)
    n_collapsed = len(validated_sequences) - len(deduped_sequences)
    stats["duplicate_sequences"] = n_collapsed
    stats["unique_sequences"] = len(deduped_sequences)

    # Write validation report
    with open(report_path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("tRNA Reference Validation Report\n")
        f.write("=" * 70 + "\n\n")

        f.write("Configuration:\n")
        f.write(f"  Input file: {input_fasta}\n")
        f.write(f"  5' adapter: {adapter_5p} ({len_5p} bp)\n")
        for i, adapter_3p in enumerate(adapters_3p):
            f.write(f"  3' adapter [{i + 1}]: {adapter_3p} ({len(adapter_3p)} bp)\n")
        f.write("\n")

        f.write("Statistics:\n")
        f.write(f"  Total sequences: {stats['total_sequences']}\n")
        f.write(f"  Valid 5' adapters: {stats['valid_5p_adapter']}\n")
        f.write(f"  Valid 3' adapters: {stats['valid_3p_adapter']}\n")
        f.write(f"  Valid CCA endings: {stats['valid_cca']}\n")
        f.write(f"  Duplicate names: {stats['duplicate_names']}\n")
        f.write(f"  Unique sequences: {stats['unique_sequences']}\n")
        f.write(f"  Duplicate sequences collapsed: {stats['duplicate_sequences']}\n\n")

        if collapsed_map:
            f.write(f"Collapsed duplicates ({n_collapsed} sequences removed):\n")
            f.writelines(
                f"  {kept} <- {', '.join(dropped)}\n"
                for kept, dropped in sorted(collapsed_map.items())
            )
            f.write("\n")

        if warnings:
            f.write(f"WARNINGS ({len(warnings)}):\n")
            f.writelines(f"  - {warn}\n" for warn in warnings)
            f.write("\n")

        if errors:
            f.write(f"ERRORS ({len(errors)}):\n")
            f.writelines(f"  - {err}\n" for err in errors)
            f.write("\n")
            f.write("VALIDATION FAILED\n")
        else:
            f.write("VALIDATION PASSED\n")

    # Fail if errors
    if errors:
        print(
            f"Validation FAILED with {len(errors)} errors. See {report_path}",
            file=sys.stderr,
        )
        for err in errors[:5]:  # Print first 5 errors
            print(f"  ERROR: {err}", file=sys.stderr)
        if len(errors) > 5:
            print(f"  ... and {len(errors) - 5} more errors", file=sys.stderr)
        sys.exit(1)

    # Write deduplicated validated FASTA
    write_fasta(deduped_sequences, output_fasta)

    print(f"Validation PASSED: {stats['total_sequences']} sequences validated")
    if n_collapsed > 0:
        print(
            f"  Deduplicated: {stats['total_sequences']} -> {stats['unique_sequences']} unique sequences"
        )
    print(f"Report written to: {report_path}")
    print(f"Validated reference written to: {output_fasta}")

    return True


def build_reference(
    input_fasta, output_fasta, report_path, adapter_5p, adapter_3p, max_mismatch=0
):
    """
    Build adapted reference from raw tRNA sequences.

    Process:
    1. Read raw tRNA sequences
    2. Check for CCA ending - if missing, add it (with warning)
    3. Prepend 5' adapter (including first tRNA base as variable position)
    4. Append 3' adapter after CCA
    5. Write adapted FASTA

    The resulting structure is:
      5' adapter (23bp) + first_base + tRNA + 3' adapter (40bp)
    """
    errors = []
    warnings = []
    stats = defaultdict(int)
    adapted_sequences = []
    seen_names = set()
    # Original input sequence per name, used as the clustering key when
    # max_mismatch > 0: grouping is on the molecule as provided, before CCA is
    # appended, so added CCA does not shift the equal-length buckets
    raw_by_name = {}

    # Verify 3' adapter starts with GGC
    if not adapter_3p.startswith("GGC"):
        errors.append(
            f"3' adapter must start with GGC for CCAGGC junction. Got: {adapter_3p[:3]}"
        )

    for name, seq in read_fasta(input_fasta):
        stats["total_sequences"] += 1

        # Check for duplicates
        if name in seen_names:
            errors.append(f"Duplicate sequence name: {name}")
            stats["duplicate_names"] += 1
            continue
        seen_names.add(name)
        raw_by_name[name] = seq

        # Check for CCA ending - add if missing
        if seq.endswith("CCA"):
            stats["had_cca"] += 1
            trna_seq = seq
        else:
            # Add CCA to the sequence
            trna_seq = seq + "CCA"
            stats["cca_added"] += 1
            warnings.append(
                f"{name}: CCA added to sequence "
                f"(original ended with '{seq[-3:] if len(seq) >= 3 else seq}')"
            )

        # Build adapted sequence
        # The 5' adapter is 23bp, then we add the first tRNA base
        # (which becomes the variable N position in the adapter scheme)
        first_base = trna_seq[0]
        adapted_seq = adapter_5p + first_base + trna_seq + adapter_3p

        # Verify the CCAGGC junction was created correctly
        # CCA should be at position: len(adapter_5p) + 1 + len(trna_seq) - 3
        junction_start = len(adapter_5p) + 1 + len(trna_seq) - 3
        junction = adapted_seq[junction_start : junction_start + 6]
        if junction != "CCAGGC":
            errors.append(f"{name}: Failed to create CCAGGC junction. Got: {junction}")
            continue

        adapted_sequences.append((name, adapted_seq))
        stats["sequences_built"] += 1

    # Deduplicate adapted sequences with identical content
    deduped_sequences, collapsed_map = deduplicate_sequences(adapted_sequences)
    n_collapsed = len(adapted_sequences) - len(deduped_sequences)
    stats["duplicate_sequences"] = n_collapsed
    stats["unique_sequences"] = len(deduped_sequences)

    # Optionally merge near-identical molecules. Clustering is on the original
    # input sequences, so equal-length buckets reflect the lengths as provided
    # rather than lengths after CCA was appended to some sequences.
    if max_mismatch > 0:
        n_before = len(deduped_sequences)
        keys = [(name, raw_by_name[name]) for name, _ in deduped_sequences]
        kept_names, similar_map = collapse_similar_sequences(keys, max_mismatch)

        adapted_by_name = dict(deduped_sequences)
        deduped_sequences = [(name, adapted_by_name[name]) for name in kept_names]

        # Fold merged names into the collapsed map, carrying along any exact
        # duplicates that had already been folded into a merged name
        for leader, merged in similar_map.items():
            for merged_name in merged:
                collapsed_map.setdefault(leader, []).append(merged_name)
                collapsed_map[leader].extend(collapsed_map.pop(merged_name, []))

        stats["similar_sequences"] = n_before - len(deduped_sequences)
        stats["unique_sequences"] = len(deduped_sequences)
        stats["max_mismatch"] = max_mismatch

    # Write build report
    with open(report_path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("tRNA Reference Build Report\n")
        f.write("=" * 70 + "\n\n")

        f.write("Configuration:\n")
        f.write(f"  Input file: {input_fasta}\n")
        f.write(f"  5' adapter: {adapter_5p} ({len(adapter_5p)} bp)\n")
        f.write(f"  3' adapter: {adapter_3p} ({len(adapter_3p)} bp)\n\n")

        f.write("Input Statistics:\n")
        f.write(f"  Total raw sequences: {stats['total_sequences']}\n")
        f.write(f"  Sequences with CCA: {stats['had_cca']}\n")
        f.write(f"  Sequences with CCA added: {stats['cca_added']}\n")
        f.write(f"  Duplicate names: {stats['duplicate_names']}\n\n")

        f.write("Output Statistics:\n")
        f.write(f"  Sequences built: {stats['sequences_built']}\n")
        f.write(f"  Unique sequences: {stats['unique_sequences']}\n")
        f.write(f"  Duplicate sequences collapsed: {stats['duplicate_sequences']}\n")
        if max_mismatch > 0:
            f.write(
                f"  Near-identical sequences merged: {stats['similar_sequences']}"
                f" (hamming <= {max_mismatch}, equal-length only,"
                f" measured on the input sequences)\n"
            )
        f.write("\n")

        if collapsed_map:
            n_removed = stats["duplicate_sequences"] + stats.get("similar_sequences", 0)
            f.write(f"Collapsed sequences ({n_removed} sequences removed):\n")
            f.writelines(
                f"  {kept} <- {', '.join(dropped)}\n"
                for kept, dropped in sorted(collapsed_map.items())
            )
            f.write("\n")

        if warnings:
            f.write(f"WARNINGS ({len(warnings)}):\n")
            f.writelines(f"  - {warn}\n" for warn in warnings)
            f.write("\n")

        if errors:
            f.write(f"ERRORS ({len(errors)}):\n")
            f.writelines(f"  - {err}\n" for err in errors)
            f.write("\n")
            f.write("BUILD FAILED\n")
        else:
            f.write("BUILD SUCCESSFUL\n")

    # Fail if errors
    if errors:
        print(
            f"Build FAILED with {len(errors)} errors. See {report_path}",
            file=sys.stderr,
        )
        for err in errors[:5]:
            print(f"  ERROR: {err}", file=sys.stderr)
        if len(errors) > 5:
            print(f"  ... and {len(errors) - 5} more errors", file=sys.stderr)
        sys.exit(1)

    # Write deduplicated adapted FASTA
    write_fasta(deduped_sequences, output_fasta)

    print(f"Build SUCCESSFUL: {stats['sequences_built']} sequences created")
    n_unique = stats["sequences_built"] - n_collapsed
    if n_collapsed > 0:
        print(
            f"  Deduplicated: {stats['sequences_built']} -> {n_unique} unique sequences"
        )
    if max_mismatch > 0:
        print(
            f"  Merged near-identical (hamming <= {max_mismatch}, equal-length): "
            f"{n_unique} -> {stats['unique_sequences']}"
        )
    if stats["cca_added"] > 0:
        print(f"  Note: CCA was added to {stats['cca_added']} sequences")
    print(f"Output written to: {output_fasta}")
    print(f"Report written to: {report_path}")

    return True


def trim_reference(input_fasta, output_fasta, adapter_5p, adapter_3p):
    """
    Trim adapter sequences from an adapted reference to produce tRNA-only FASTA.

    The adapted reference structure is:
      5' adapter + N (first tRNA base) + tRNA sequence + 3' adapter

    This function strips the 5' prefix (adapter + N variable position) and
    the 3' adapter to produce sequences containing only the tRNA portion.

    Args:
        input_fasta: Path to adapted reference FASTA
        output_fasta: Path for tRNA-only output FASTA
        adapter_5p: 5' adapter sequence (used to compute prefix length)
        adapter_3p: 3' adapter sequence (used to compute suffix length)
    """
    offset_5p = len(adapter_5p) + 1  # adapter + N variable position
    offset_3p = len(adapter_3p)

    trimmed_sequences = []

    for name, seq in read_fasta(input_fasta):
        if len(seq) <= offset_5p + offset_3p:
            print(
                f"WARNING: {name} too short to trim ({len(seq)} bp, "
                f"need > {offset_5p + offset_3p} bp). Skipping.",
                file=sys.stderr,
            )
            continue

        trna_seq = seq[offset_5p:-offset_3p] if offset_3p > 0 else seq[offset_5p:]
        trimmed_sequences.append((name, trna_seq))

    write_fasta(trimmed_sequences, output_fasta)

    print(
        f"Trimmed {len(trimmed_sequences)} sequences "
        f"(removed {offset_5p} bp 5' prefix, {offset_3p} bp 3' suffix)"
    )
    print(f"Output written to: {output_fasta}")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Validate or build tRNA reference FASTA with adapters",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate existing adapted reference
  python build_trna_reference.py --mode validate \\
      --input reference.fa --output validated.fa --report report.txt

  # Build adapted reference from raw tRNAs
  python build_trna_reference.py --mode build \\
      --input raw_trnas.fa --output adapted.fa --report report.txt

  # Trim adapters to produce tRNA-only FASTA
  python build_trna_reference.py --mode trim \\
      --input adapted.fa --output trna_only.fa

Build mode will:
  - Add CCA to sequences that don't already have it (with warning)
  - Prepend 5' adapter sequence
  - Append 3' adapter sequence after CCA

The CCAGGC junction (tRNA CCA + adapter GGC) is required for the
charging classification model to work correctly.
        """,
    )

    parser.add_argument(
        "--mode",
        choices=["validate", "build", "trim"],
        required=True,
        help="Operation mode: validate existing reference, build new one, or trim adapters",
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Input FASTA file (adapted reference for validate/trim, raw tRNAs for build)",
    )
    parser.add_argument("--output", "-o", required=True, help="Output FASTA file")
    parser.add_argument(
        "--max-mismatch",
        type=int,
        default=0,
        help="Build mode: also merge sequences within this many substitutions of "
        "a cluster representative. 0 (default) collapses exact duplicates only. "
        "Distance is measured on the input sequences, before CCA is appended, "
        "and only between equal-length sequences, so molecules differing by an "
        "indel are never merged.",
    )
    parser.add_argument(
        "--report",
        "-r",
        default=None,
        help="Output validation/build report file (not used for trim mode)",
    )
    parser.add_argument(
        "--adapter-5p",
        default="CCTAAGAGCAAGAAGAAGCCTGG",
        help="5' adapter sequence (default: %(default)s)",
    )
    parser.add_argument(
        "--adapter-3p",
        action="append",
        default=None,
        help="3' adapter sequence (may be specified multiple times for multi-adapter references)",
    )

    args = parser.parse_args()

    # Default 3' adapter if none specified
    adapters_3p = args.adapter_3p or ["GGCTTCTTCTTGCTCTTCCAACCTTGCCTTAAAAAAAAAA"]

    if args.mode == "validate":
        if not args.report:
            parser.error("--report is required for validate mode")
        validate_reference(
            args.input, args.output, args.report, args.adapter_5p, adapters_3p
        )
    elif args.mode == "build":
        if not args.report:
            parser.error("--report is required for build mode")
        # Build mode uses only the first adapter
        build_reference(
            args.input,
            args.output,
            args.report,
            args.adapter_5p,
            adapters_3p[0],
            max_mismatch=args.max_mismatch,
        )
    else:
        # Trim mode uses the first adapter for length calculation
        trim_reference(args.input, args.output, args.adapter_5p, adapters_3p[0])


if __name__ == "__main__":
    main()
