#!/usr/bin/env python
"""
Add adapter position tags to BAM file using parasail semi-global alignment.

Uses SAM-spec pt:Z: tag format for read annotations:
    pt:Z:start;end;strand;type|start;end;strand;type

Example:
    pt:Z:0;24;+;5p_adapter|118;135;+;3p_adapter

Positions are 0-based, relative to the read sequence.
If an adapter is not found, that annotation is omitted.
"""

import argparse
import re
import sys
from collections import defaultdict

import parasail
import pysam

# Default adapter sequences
DEFAULT_ADAPTER_5P = "CCTAAGAGCAAGAAGAAGCCTGG"  # 23bp (excluding trailing N)
DEFAULT_ADAPTER_3P = "GGCTTCTTCTTGCTCTTCCAACCTTGCCTTAAAAAAAAAA"  # 40bp full 3' adapter

# Default scoring parameters for ~15% error tolerance
DEFAULT_MATCH = 2
DEFAULT_MISMATCH = -1
DEFAULT_GAP_OPEN = 2
DEFAULT_GAP_EXTEND = 1

# Default minimum scores
# Semi-global alignment allows partial adapter matches (truncated at either end)
# 5' adapter (23bp): 23 * 2 * 0.85 ≈ 39, use 30 for more tolerance
# 3' adapter (40bp): use 20 to allow partial matches (at least ~12bp with good alignment)
DEFAULT_MIN_SCORE_5P = 30
DEFAULT_MIN_SCORE_3P = 20

# Search regions
DEFAULT_SEARCH_5P = 60  # Search first N bp for 5' adapter
DEFAULT_SEARCH_3P = 80  # Search last N bp for 3' adapter


def create_scoring_matrix(match=DEFAULT_MATCH, mismatch=DEFAULT_MISMATCH):
    """Create parasail scoring matrix."""
    return parasail.matrix_create("ACGT", match, mismatch)


def get_adapter_bounds_from_cigar(cigar_decode, beg_ref):
    """
    Get the actual adapter match bounds from cigar, excluding leading/trailing deletions.

    Semi-global alignment with free query end gaps can produce leading/trailing D ops
    representing parts of the reference (read region) that don't align to the adapter.

    Args:
        cigar_decode: Decoded cigar string (bytes)
        beg_ref: Starting reference position from cigar

    Returns:
        (start, end) tuple with 0-based positions, end is exclusive
    """
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar_decode.decode())

    # Skip leading deletions to find where adapter match actually starts
    ref_pos = beg_ref
    leading_idx = 0
    for i, (length, op) in enumerate(ops):
        length = int(length)
        if op == "D":
            ref_pos += length
            leading_idx = i + 1
        else:
            break
    start = ref_pos

    # Count trailing deletions
    trailing_d_len = 0
    for length, op in reversed(ops):
        length = int(length)
        if op == "D":
            trailing_d_len += length
        else:
            break

    # Walk through remaining ops to find end
    for length, op in ops[leading_idx:]:
        length = int(length)
        if op in "MDN=X":
            ref_pos += length

    end = ref_pos - trailing_d_len
    return start, end


def find_adapter(
    read_seq,
    adapter_seq,
    search_start,
    search_end,
    matrix,
    gap_open,
    gap_extend,
    min_score,
):
    """
    Find adapter in a region of the read using semi-global alignment.

    Args:
        read_seq: Full read sequence
        adapter_seq: Adapter sequence to find
        search_start: Start position in read to search (0-based)
        search_end: End position in read to search (exclusive)
        matrix: Parasail scoring matrix
        gap_open: Gap open penalty
        gap_extend: Gap extension penalty
        min_score: Minimum alignment score to accept

    Returns:
        (start, end, score) tuple with positions relative to full read,
        or None if adapter not found
    """
    region = read_seq[search_start:search_end]
    if len(region) < len(adapter_seq) // 2:
        return None

    # Semi-global alignment: free end gaps on the reference (read region)
    # sg_dx: free gaps at both ends of database (reference/read region)
    # This allows the adapter to be found anywhere within the search region
    result = parasail.sg_dx_trace_scan_sat(
        adapter_seq, region, gap_open, gap_extend, matrix
    )

    if result.score < min_score:
        return None

    # Get alignment coordinates from cigar, excluding leading/trailing deletions
    cigar = result.cigar
    ref_start, ref_end = get_adapter_bounds_from_cigar(cigar.decode, cigar.beg_ref)

    # Convert to full read coordinates
    start = search_start + ref_start
    end = search_start + ref_end

    return (start, end, result.score)


def find_5p_adapter(read_seq, adapter_5p, matrix, gap_open, gap_extend, min_score):
    """Find 5' adapter in the beginning of the read."""
    search_len = min(DEFAULT_SEARCH_5P, len(read_seq))
    return find_adapter(
        read_seq, adapter_5p, 0, search_len, matrix, gap_open, gap_extend, min_score
    )


def find_3p_adapter(read_seq, adapter_3p, matrix, gap_open, gap_extend, min_score):
    """Find 3' adapter in the end of the read."""
    search_len = min(DEFAULT_SEARCH_3P, len(read_seq))
    search_start = max(0, len(read_seq) - search_len)
    return find_adapter(
        read_seq,
        adapter_3p,
        search_start,
        len(read_seq),
        matrix,
        gap_open,
        gap_extend,
        min_score,
    )


def find_best_3p_adapter(read_seq, adapters, matrix, gap_open, gap_extend, min_score):
    """Find best matching 3' adapter from list of (name, seq) tuples.

    Args:
        read_seq: Full read sequence
        adapters: List of (name, sequence) tuples for 3' adapters
        matrix: Parasail scoring matrix
        gap_open: Gap open penalty
        gap_extend: Gap extension penalty
        min_score: Minimum alignment score to accept

    Returns:
        (start, end, score, name) tuple if found, or None if no adapter matches
    """
    best_result = None
    best_name = None
    best_score = min_score - 1

    for name, adapter_seq in adapters:
        result = find_3p_adapter(
            read_seq, adapter_seq, matrix, gap_open, gap_extend, min_score
        )
        if result and result[2] > best_score:
            best_result = result
            best_name = name
            best_score = result[2]

    if best_result:
        return (*best_result, best_name)  # (start, end, score, name)
    return None


def format_pt_tag(adapter_5p_result, adapter_3p_result):
    """
    Format adapter positions as SAM-spec pt tag.

    Format: pt:Z:start;end;strand;type|start;end;strand;type

    For 3' adapters with names, the type will be "3p_adapter_<name>"
    (e.g., "3p_adapter_v1" or "3p_adapter_v2").
    For single/default adapters, the type remains "3p_adapter".
    """
    annotations = []

    if adapter_5p_result:
        start, end, score = adapter_5p_result
        annotations.append(f"{start};{end};+;5p_adapter")

    if adapter_3p_result:
        # Handle both old format (start, end, score) and new format (start, end, score, name)
        if len(adapter_3p_result) == 4:
            start, end, score, name = adapter_3p_result
            if name == "default":
                annotations.append(f"{start};{end};+;3p_adapter")
            else:
                annotations.append(f"{start};{end};+;3p_adapter_{name}")
        else:
            start, end, score = adapter_3p_result
            annotations.append(f"{start};{end};+;3p_adapter")

    if not annotations:
        return None

    return "|".join(annotations)


class Stats:
    """Track adapter detection statistics."""

    def __init__(self):
        self.total = 0
        self.with_5p = 0
        self.with_3p = 0
        self.with_both = 0
        self.with_neither = 0
        self.adapter_3p_counts = defaultdict(int)  # counts by adapter name

    def update(self, has_5p, has_3p, adapter_3p_name=None):
        self.total += 1
        if has_5p:
            self.with_5p += 1
        if has_3p:
            self.with_3p += 1
            if adapter_3p_name:
                self.adapter_3p_counts[adapter_3p_name] += 1
        if has_5p and has_3p:
            self.with_both += 1
        if not has_5p and not has_3p:
            self.with_neither += 1

    def summary(self):
        lines = [
            f"total_reads {self.total}",
            f"with_5p_adapter {self.with_5p}",
            f"with_3p_adapter {self.with_3p}",
            f"with_both_adapters {self.with_both}",
            f"with_neither_adapter {self.with_neither}",
        ]
        # Add per-adapter 3' counts if there are multiple adapters
        if self.adapter_3p_counts:
            for name in sorted(self.adapter_3p_counts.keys()):
                lines.append(f"with_3p_adapter_{name} {self.adapter_3p_counts[name]}")
        return "\n".join(lines)


def process_bam(
    input_bam,
    output_bam,
    adapter_5p,
    adapter_3p_list,
    min_score_5p,
    min_score_3p,
    match,
    mismatch,
    gap_open,
    gap_extend,
    infer_5p_from_alignment=False,
    max_ref_start_for_5p=20,
):
    """Process BAM file and add adapter position tags.

    Args:
        adapter_3p_list: List of (name, sequence) tuples for 3' adapters
    """
    matrix = create_scoring_matrix(match, mismatch)
    stats = Stats()
    adapter_5p_len = len(adapter_5p)

    with pysam.AlignmentFile(input_bam, "rb") as inbam:
        with pysam.AlignmentFile(output_bam, "wb", header=inbam.header) as outbam:
            for read in inbam:
                seq = read.query_sequence

                if seq is None:
                    # Unmapped read without sequence
                    outbam.write(read)
                    continue

                # Find adapters using sequence-based detection
                result_5p = find_5p_adapter(
                    seq, adapter_5p, matrix, gap_open, gap_extend, min_score_5p
                )
                result_3p = find_best_3p_adapter(
                    seq, adapter_3p_list, matrix, gap_open, gap_extend, min_score_3p
                )

                # Fallback: infer 5' adapter from alignment position for truncated reads
                if result_5p is None and infer_5p_from_alignment:
                    if not read.is_unmapped and read.reference_start < max_ref_start_for_5p:
                        # Adapter end position in read = adapter_len - ref_start
                        adapter_end_in_read = adapter_5p_len - read.reference_start
                        if adapter_end_in_read > 0:
                            # Estimate score based on how much adapter is present
                            estimated_score = int(adapter_end_in_read * match * 0.85)
                            result_5p = (0, adapter_end_in_read, estimated_score)

                # Update stats (extract adapter name if present)
                adapter_3p_name = result_3p[3] if result_3p else None
                stats.update(result_5p is not None, result_3p is not None, adapter_3p_name)

                # Add PT tag if any adapter found
                pt_value = format_pt_tag(result_5p, result_3p)
                if pt_value:
                    read.set_tag("pt", pt_value, "Z")

                outbam.write(read)

    return stats


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file")

    parser.add_argument(
        "--adapter-5p",
        default=DEFAULT_ADAPTER_5P,
        help=f"5' adapter sequence (default: {DEFAULT_ADAPTER_5P})",
    )
    parser.add_argument(
        "--adapter-3p",
        action="append",
        default=[],
        dest="adapter_3p_list",
        help=(
            "3' adapter as 'name:sequence' or just 'sequence' (can specify multiple). "
            f"Default if none specified: {DEFAULT_ADAPTER_3P}"
        ),
    )

    parser.add_argument(
        "--min-score-5p",
        type=int,
        default=DEFAULT_MIN_SCORE_5P,
        help=f"Minimum alignment score for 5' adapter (default: {DEFAULT_MIN_SCORE_5P})",
    )
    parser.add_argument(
        "--min-score-3p",
        type=int,
        default=DEFAULT_MIN_SCORE_3P,
        help=f"Minimum alignment score for 3' adapter (default: {DEFAULT_MIN_SCORE_3P})",
    )

    parser.add_argument(
        "--match",
        type=int,
        default=DEFAULT_MATCH,
        help=f"Match score (default: {DEFAULT_MATCH})",
    )
    parser.add_argument(
        "--mismatch",
        type=int,
        default=DEFAULT_MISMATCH,
        help=f"Mismatch penalty (default: {DEFAULT_MISMATCH})",
    )
    parser.add_argument(
        "--gap-open",
        type=int,
        default=DEFAULT_GAP_OPEN,
        help=f"Gap open penalty (default: {DEFAULT_GAP_OPEN})",
    )
    parser.add_argument(
        "--gap-extend",
        type=int,
        default=DEFAULT_GAP_EXTEND,
        help=f"Gap extension penalty (default: {DEFAULT_GAP_EXTEND})",
    )

    parser.add_argument(
        "--infer-5p-from-alignment",
        action="store_true",
        help="Infer 5' adapter presence from alignment position for truncated reads",
    )
    parser.add_argument(
        "--max-ref-start-for-5p",
        type=int,
        default=20,
        help="Max reference start position to infer 5' adapter (default: 20)",
    )

    args = parser.parse_args()

    # Parse 3' adapter list
    adapter_3p_list = []
    if args.adapter_3p_list:
        for adapter_spec in args.adapter_3p_list:
            if ":" in adapter_spec:
                name, seq = adapter_spec.split(":", 1)
                adapter_3p_list.append((name, seq))
            else:
                # No name provided, use "default"
                adapter_3p_list.append(("default", adapter_spec))
    else:
        # No adapters specified, use default
        adapter_3p_list = [("default", DEFAULT_ADAPTER_3P)]

    stats = process_bam(
        args.input,
        args.output,
        args.adapter_5p,
        adapter_3p_list,
        args.min_score_5p,
        args.min_score_3p,
        args.match,
        args.mismatch,
        args.gap_open,
        args.gap_extend,
        args.infer_5p_from_alignment,
        args.max_ref_start_for_5p,
    )

    print(stats.summary(), file=sys.stderr)


if __name__ == "__main__":
    main()
