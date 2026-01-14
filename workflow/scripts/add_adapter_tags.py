#!/usr/bin/env python
"""
Add adapter position tags to BAM file using parasail semi-global alignment.

Uses SAM-spec PT:Z: tag format for read annotations:
    PT:Z:start;end;strand;type|start;end;strand;type

Example:
    PT:Z:0;24;+;5p_adapter|118;135;+;3p_adapter

Positions are 0-based, relative to the read sequence.
If an adapter is not found, that annotation is omitted.
"""

import argparse
import re
import sys

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
    # sg_dx: free gaps at end of database (reference/read)
    # sg_qx: free gaps at end of query (adapter)
    # We want adapter to be able to be truncated, so use sg_qe (free query end)
    # and sg_db (free database begin) for flexibility
    result = parasail.sg_qx_trace_scan_sat(
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


def format_pt_tag(adapter_5p_result, adapter_3p_result):
    """
    Format adapter positions as SAM-spec PT tag.

    Format: PT:Z:start;end;strand;type|start;end;strand;type
    """
    annotations = []

    if adapter_5p_result:
        start, end, score = adapter_5p_result
        annotations.append(f"{start};{end};+;5p_adapter")

    if adapter_3p_result:
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

    def update(self, has_5p, has_3p):
        self.total += 1
        if has_5p:
            self.with_5p += 1
        if has_3p:
            self.with_3p += 1
        if has_5p and has_3p:
            self.with_both += 1
        if not has_5p and not has_3p:
            self.with_neither += 1

    def summary(self):
        return (
            f"total_reads {self.total}\n"
            f"with_5p_adapter {self.with_5p}\n"
            f"with_3p_adapter {self.with_3p}\n"
            f"with_both_adapters {self.with_both}\n"
            f"with_neither_adapter {self.with_neither}"
        )


def process_bam(
    input_bam,
    output_bam,
    adapter_5p,
    adapter_3p,
    min_score_5p,
    min_score_3p,
    match,
    mismatch,
    gap_open,
    gap_extend,
):
    """Process BAM file and add adapter position tags."""
    matrix = create_scoring_matrix(match, mismatch)
    stats = Stats()

    with pysam.AlignmentFile(input_bam, "rb") as inbam:
        with pysam.AlignmentFile(output_bam, "wb", header=inbam.header) as outbam:
            for read in inbam:
                seq = read.query_sequence

                if seq is None:
                    # Unmapped read without sequence
                    outbam.write(read)
                    continue

                # Find adapters
                result_5p = find_5p_adapter(
                    seq, adapter_5p, matrix, gap_open, gap_extend, min_score_5p
                )
                result_3p = find_3p_adapter(
                    seq, adapter_3p, matrix, gap_open, gap_extend, min_score_3p
                )

                # Update stats
                stats.update(result_5p is not None, result_3p is not None)

                # Add PT tag if any adapter found
                pt_value = format_pt_tag(result_5p, result_3p)
                if pt_value:
                    read.set_tag("PT", pt_value, "Z")

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
        default=DEFAULT_ADAPTER_3P,
        help=f"3' adapter sequence (default: {DEFAULT_ADAPTER_3P})",
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

    args = parser.parse_args()

    stats = process_bam(
        args.input,
        args.output,
        args.adapter_5p,
        args.adapter_3p,
        args.min_score_5p,
        args.min_score_3p,
        args.match,
        args.mismatch,
        args.gap_open,
        args.gap_extend,
    )

    print(stats.summary(), file=sys.stderr)


if __name__ == "__main__":
    main()
