#!/usr/bin/env python
"""
Compute pairwise sequence similarity matrix for reference FASTA sequences.

This script calculates all-vs-all pairwise alignments using Needleman-Wunsch
global alignment to identify potential cross-mapping issues due to homologous
tRNA sequences.

Output: TSV file with square similarity matrix (percent identity values).
"""

import argparse
import sys

import numpy as np
import parasail
import pysam


def compute_similarity_matrix(fasta_path):
    """
    Compute all-vs-all pairwise similarity matrix for sequences in FASTA.

    Uses global alignment (Needleman-Wunsch) with scoring parameters matching
    the adapter detection in add_adapter_tags.py.

    Args:
        fasta_path: Path to reference FASTA file

    Returns:
        Tuple of (similarity_matrix, sequence_names)
        - similarity_matrix: numpy array of percent identity values
        - sequence_names: list of sequence IDs
    """
    # Read all sequences from FASTA using pysam
    faidx = pysam.FastaFile(fasta_path)
    names = list(faidx.references)
    n = len(names)

    if n == 0:
        sys.exit(f"No sequences found in {fasta_path}")

    # Fetch all sequences
    sequences = [faidx.fetch(name).upper() for name in names]
    faidx.close()

    # Initialize similarity matrix
    matrix = np.zeros((n, n))

    # Scoring parameters (matching add_adapter_tags.py)
    gap_open = 2
    gap_extend = 1

    # Process each pair
    for i in range(n):
        seq_i = sequences[i]
        len_i = len(seq_i)

        # Diagonal is always 100% identity
        matrix[i, i] = 100.0

        for j in range(i + 1, n):
            seq_j = sequences[j]
            len_j = len(seq_j)

            # Global alignment with statistics
            result = parasail.nw_stats_scan_sat(
                seq_i, seq_j, gap_open, gap_extend, parasail.dnafull
            )

            # Percent identity: matches / max(len_seq1, len_seq2) * 100
            max_len = max(len_i, len_j)
            pct_id = result.matches / max_len * 100

            # Fill symmetric matrix
            matrix[i, j] = pct_id
            matrix[j, i] = pct_id

    return matrix, names


def write_matrix_tsv(matrix, names, output_path):
    """
    Write similarity matrix as TSV file.

    Args:
        matrix: numpy array of similarity values
        names: list of sequence names
        output_path: path for output TSV file
    """
    with open(output_path, "w") as f:
        # Header row
        f.write("\t" + "\t".join(names) + "\n")

        # Data rows
        for i, name in enumerate(names):
            row_values = [f"{matrix[i, j]:.2f}" for j in range(len(names))]
            f.write(name + "\t" + "\t".join(row_values) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("fasta", help="Path to reference FASTA file")
    parser.add_argument("output", help="Path for output TSV file")

    args = parser.parse_args()

    # Compute similarity matrix
    print(f"Computing similarity matrix for {args.fasta}...", file=sys.stderr)
    matrix, names = compute_similarity_matrix(args.fasta)
    print(f"Processed {len(names)} sequences", file=sys.stderr)

    # Write output
    write_matrix_tsv(matrix, names, args.output)
    print(f"Wrote similarity matrix to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
