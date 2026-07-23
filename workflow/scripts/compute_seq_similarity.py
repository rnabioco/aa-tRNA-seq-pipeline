#!/usr/bin/env python
"""
Compute pairwise sequence similarity matrix for reference FASTA sequences.

This script calculates all-vs-all pairwise alignments using Needleman-Wunsch
global alignment to identify potential cross-mapping issues due to homologous
tRNA sequences.

Alignment is quadratic in the number of sequences, so references are collapsed
before aligning:

  - Exact duplicates are always collapsed. This is lossless: identical
    sequences have identical similarity to everything else. Multi-copy tRNA
    gene families make this a large win (danRer11 mature tRNAs: 8879 records,
    3315 distinct sequences).

  - With --max-mismatch > 0, near-identical sequences are additionally
    collapsed by Hamming distance. This is lossy, so the matrix is reported
    over cluster representatives rather than expanded back to every input
    name, and the cluster membership table records what was merged.

Output: TSV file with square similarity matrix (percent identity values),
plus a sidecar TSV of cluster membership.
"""

import argparse
import sys
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import parasail
import pysam

# Scoring parameters (matching add_adapter_tags.py)
GAP_OPEN = 2
GAP_EXTEND = 1


def load_sequences(fasta_path):
    """Read all sequences from a FASTA, returning (names, sequences)."""
    # Regenerate .fai index to ensure it matches the current FASTA content
    pysam.faidx(fasta_path)

    faidx = pysam.FastaFile(fasta_path)
    names = list(faidx.references)
    sequences = [faidx.fetch(name).upper() for name in names]
    faidx.close()

    if not names:
        sys.exit(f"No sequences found in {fasta_path}")

    return names, sequences


def group_exact(sequences):
    """
    Collapse exact duplicate sequences.

    Returns a list of clusters, each a list of input indices sharing an
    identical sequence. The first index of each cluster is its representative.
    Input order is preserved.
    """
    groups = OrderedDict()
    for i, seq in enumerate(sequences):
        groups.setdefault(seq, []).append(i)
    return list(groups.values())


def group_hamming(clusters, sequences, max_mismatch):
    """
    Further merge clusters whose representatives are within max_mismatch
    substitutions of each other.

    Uses greedy leader clustering: each representative joins the first existing
    leader it is within max_mismatch of, otherwise becomes a new leader. Unlike
    single-linkage this bounds the cluster radius, so every member is within
    max_mismatch of its leader rather than merely chained to it.

    Hamming distance is only defined for equal-length sequences, so sequences
    differing by an indel are never merged regardless of max_mismatch.
    """
    if max_mismatch <= 0:
        return clusters

    # Bucket clusters by representative length; only equal lengths are comparable
    by_length = OrderedDict()
    for cluster in clusters:
        seq = sequences[cluster[0]]
        by_length.setdefault(len(seq), []).append(cluster)

    merged = []
    for length, bucket in by_length.items():
        leaders = np.empty((0, length), dtype=np.uint8)
        leader_clusters = []

        for cluster in bucket:
            row = np.frombuffer(sequences[cluster[0]].encode(), dtype=np.uint8)
            if len(leader_clusters):
                dist = (leaders != row).sum(axis=1)
                hit = np.flatnonzero(dist <= max_mismatch)
                if hit.size:
                    leader_clusters[int(hit[0])].extend(cluster)
                    continue
            leaders = np.vstack([leaders, row])
            leader_clusters.append(list(cluster))

        merged.extend(leader_clusters)

    return merged


def _fill_rows(rep_seqs, rows, matrix):
    """Align each row index in `rows` against all higher-indexed sequences.

    Writes into the upper triangle of `matrix`. Rows are disjoint across
    workers, so no locking is needed. parasail releases the GIL, so threads
    give real parallelism here.
    """
    n = len(rep_seqs)
    for i in rows:
        seq_i = rep_seqs[i]
        len_i = len(seq_i)
        # Encode the query once and reuse it across the whole row
        profile = parasail.profile_create_stats_sat(seq_i, parasail.dnafull)
        for j in range(i + 1, n):
            seq_j = rep_seqs[j]
            result = parasail.nw_stats_scan_profile_sat(
                profile, seq_j, GAP_OPEN, GAP_EXTEND
            )
            # Percent identity: matches / max(len_seq1, len_seq2) * 100
            matrix[i, j] = result.matches / max(len_i, len(seq_j)) * 100


def align_all_pairs(rep_seqs, threads=1):
    """Compute the percent-identity matrix for a list of sequences."""
    n = len(rep_seqs)
    matrix = np.zeros((n, n))

    # Row i costs n-i-1 alignments, so deal rows round-robin to balance workers
    threads = max(1, min(threads, n))
    stripes = [range(t, n, threads) for t in range(threads)]

    if threads == 1:
        _fill_rows(rep_seqs, stripes[0], matrix)
    else:
        with ThreadPoolExecutor(max_workers=threads) as pool:
            list(pool.map(lambda r: _fill_rows(rep_seqs, r, matrix), stripes))

    # Mirror the upper triangle and set the diagonal. Not `matrix += matrix.T`:
    # that reads and writes overlapping memory through a view of itself.
    matrix = matrix + matrix.T
    np.fill_diagonal(matrix, 100.0)
    return matrix


def compute_similarity(names, sequences, threads=1, max_mismatch=0, expand=True):
    """
    Compute the pairwise similarity matrix for already-loaded sequences.

    Uses global alignment (Needleman-Wunsch) with scoring parameters matching
    the adapter detection in add_adapter_tags.py.

    Args:
        names: list of sequence IDs
        sequences: list of sequence strings, parallel to names
        threads: Number of worker threads for alignment
        max_mismatch: Additionally collapse sequences within this Hamming
            distance. 0 (default) collapses only exact duplicates.
        expand: Expand the matrix back to every input sequence name. Forced
            off when max_mismatch > 0, where clusters are not interchangeable.

    Returns:
        Tuple of (similarity_matrix, labels, clusters)
        - similarity_matrix: numpy array of percent identity values
        - labels: list of sequence IDs labelling the matrix rows/columns
        - clusters: list of clusters, each a list of input indices
    """
    clusters = group_exact(sequences)
    n_exact = len(clusters)
    clusters = group_hamming(clusters, sequences, max_mismatch)

    print(
        f"{len(names)} sequences -> {n_exact} distinct"
        + (
            f" -> {len(clusters)} clusters (hamming <= {max_mismatch})"
            if max_mismatch > 0
            else ""
        ),
        file=sys.stderr,
    )

    rep_seqs = [sequences[c[0]] for c in clusters]
    matrix = align_all_pairs(rep_seqs, threads=threads)

    # Lossy collapse: representatives are not interchangeable with their members
    if max_mismatch > 0:
        expand = False

    if not expand:
        return matrix, [names[c[0]] for c in clusters], clusters

    # Exact duplicates are interchangeable, so project back to every input name
    cluster_of = np.empty(len(names), dtype=np.intp)
    for k, cluster in enumerate(clusters):
        cluster_of[cluster] = k
    return matrix[np.ix_(cluster_of, cluster_of)], names, clusters


def compute_similarity_matrix(fasta_path, threads=1, max_mismatch=0, expand=True):
    """
    Compute all-vs-all pairwise similarity matrix for sequences in FASTA.

    Returns:
        Tuple of (similarity_matrix, sequence_names)
    """
    names, sequences = load_sequences(fasta_path)
    matrix, labels, _ = compute_similarity(
        names, sequences, threads=threads, max_mismatch=max_mismatch, expand=expand
    )
    return matrix, labels


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

        # Data rows; format in C rather than per-cell Python f-strings
        for name, row in zip(names, np.asarray(matrix)):
            f.write(name + "\t")
            np.savetxt(f, row[None], fmt="%.2f", delimiter="\t")


def write_clusters_tsv(clusters, names, output_path):
    """Write cluster membership so collapsed sequences remain traceable."""
    with open(output_path, "w") as f:
        f.write("cluster_id\trepresentative\tn_members\tmembers\n")
        for k, cluster in enumerate(clusters):
            members = ",".join(names[i] for i in cluster)
            f.write(f"{k}\t{names[cluster[0]]}\t{len(cluster)}\t{members}\n")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("fasta", help="Path to reference FASTA file")
    parser.add_argument("output", help="Path for output TSV file")
    parser.add_argument(
        "--threads", type=int, default=1, help="Worker threads for alignment"
    )
    parser.add_argument(
        "--max-mismatch",
        type=int,
        default=0,
        help="Additionally collapse sequences within this Hamming distance. "
        "0 (default) collapses exact duplicates only, which is lossless.",
    )
    parser.add_argument(
        "--clusters",
        help="Path for cluster membership TSV "
        "(default: <output> with a .clusters.tsv suffix)",
    )

    args = parser.parse_args()

    # Compute similarity matrix
    print(f"Computing similarity matrix for {args.fasta}...", file=sys.stderr)
    names, sequences = load_sequences(args.fasta)
    matrix, labels, clusters = compute_similarity(
        names, sequences, threads=args.threads, max_mismatch=args.max_mismatch
    )
    print(f"Processed {matrix.shape[0]} sequences", file=sys.stderr)

    # Write output
    write_matrix_tsv(matrix, labels, args.output)
    print(f"Wrote similarity matrix to {args.output}", file=sys.stderr)

    clusters_path = args.clusters or (
        args.output.removesuffix(".tsv") + ".clusters.tsv"
    )
    write_clusters_tsv(clusters, names, clusters_path)
    print(f"Wrote cluster membership to {clusters_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
