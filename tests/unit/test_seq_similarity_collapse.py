"""Tests for sequence collapsing in compute_seq_similarity.py.

The similarity matrix is quadratic in the number of sequences, so identical
sequences are collapsed before aligning. These tests pin down that the
exact-duplicate collapse is lossless and that the lossy Hamming collapse
reports what it merged.
"""

import numpy as np
import pysam
import pytest

from compute_seq_similarity import (
    compute_similarity,
    compute_similarity_matrix,
    group_exact,
    group_hamming,
    load_sequences,
    write_clusters_tsv,
)

# Two pairs of duplicates (A at 0/3, B at 1/4) plus a unique C
DUP_FASTA = (
    ">a1\nACGTACGTAA\n"
    ">b1\nTTGGTTGGTT\n"
    ">c1\nACCCACCCAC\n"
    ">a2\nACGTACGTAA\n"
    ">b2\nTTGGTTGGTT\n"
)


@pytest.fixture
def dup_fasta(temp_dir):
    fa = temp_dir / "dups.fa"
    fa.write_text(DUP_FASTA)
    pysam.faidx(str(fa))
    return fa


class TestGroupExact:
    def test_groups_duplicates(self, dup_fasta):
        _, seqs = load_sequences(str(dup_fasta))
        clusters = group_exact(seqs)
        assert clusters == [[0, 3], [1, 4], [2]]

    def test_all_unique_is_identity(self):
        clusters = group_exact(["AAAA", "CCCC", "GGGG"])
        assert clusters == [[0], [1], [2]]

    def test_covers_every_input_exactly_once(self, dup_fasta):
        _, seqs = load_sequences(str(dup_fasta))
        clusters = group_exact(seqs)
        assert sorted(i for c in clusters for i in c) == list(range(len(seqs)))


class TestGroupHamming:
    def test_zero_is_a_noop(self):
        clusters = [[0], [1]]
        assert group_hamming(clusters, ["AAAA", "AAAT"], 0) is clusters

    def test_merges_within_distance(self):
        seqs = ["AAAA", "AAAT"]
        merged = group_hamming([[0], [1]], seqs, 1)
        assert merged == [[0, 1]]

    def test_does_not_merge_beyond_distance(self):
        seqs = ["AAAA", "AATT"]
        assert group_hamming([[0], [1]], seqs, 1) == [[0], [1]]

    def test_never_merges_different_lengths(self):
        # One substitution apart in content, but different lengths, so Hamming
        # distance is undefined and they must stay separate at any threshold.
        seqs = ["AAAA", "AAAAA"]
        assert group_hamming([[0], [1]], seqs, 99) == [[0], [1]]

    def test_bounded_radius_not_single_linkage(self):
        # A-B and B-C are each 1 apart, A-C is 2 apart. Single linkage would
        # chain all three; leader clustering must not put C with A.
        seqs = ["AAAA", "AAAT", "AATT"]
        merged = group_hamming([[0], [1], [2]], seqs, 1)
        assert merged == [[0, 1], [2]]

    def test_preserves_all_members(self):
        seqs = ["AAAA", "AAAT", "CCCC"]
        merged = group_hamming([[0, 7], [1], [2]], seqs, 1)
        assert sorted(i for c in merged for i in c) == [0, 1, 2, 7]


class TestDedupIsLossless:
    def test_matches_bruteforce_all_vs_all(self, dup_fasta):
        """Collapsing duplicates must reproduce the naive all-vs-all matrix."""
        import parasail

        from compute_seq_similarity import GAP_EXTEND, GAP_OPEN

        names, seqs = load_sequences(str(dup_fasta))
        n = len(seqs)

        expected = np.zeros((n, n))
        for i in range(n):
            expected[i, i] = 100.0
            for j in range(i + 1, n):
                r = parasail.nw_stats_scan_sat(
                    seqs[i], seqs[j], GAP_OPEN, GAP_EXTEND, parasail.dnafull
                )
                pct = r.matches / max(len(seqs[i]), len(seqs[j])) * 100
                expected[i, j] = expected[j, i] = pct

        actual, labels = compute_similarity_matrix(str(dup_fasta))
        assert labels == names
        np.testing.assert_allclose(actual, expected)

    def test_expanded_shape_is_full_input(self, dup_fasta):
        matrix, labels = compute_similarity_matrix(str(dup_fasta))
        assert matrix.shape == (5, 5)
        assert len(labels) == 5

    def test_duplicates_are_exactly_100(self, dup_fasta):
        matrix, _ = compute_similarity_matrix(str(dup_fasta))
        assert matrix[0, 3] == 100.0
        assert matrix[1, 4] == 100.0

    @pytest.mark.parametrize("threads", [1, 2, 4, 8])
    def test_threading_is_deterministic(self, dup_fasta, threads):
        serial, _ = compute_similarity_matrix(str(dup_fasta), threads=1)
        parallel, _ = compute_similarity_matrix(str(dup_fasta), threads=threads)
        np.testing.assert_array_equal(serial, parallel)


class TestHammingCollapseReporting:
    def test_collapse_reports_representatives_not_all_names(self, dup_fasta):
        names, seqs = load_sequences(str(dup_fasta))
        matrix, labels, clusters = compute_similarity(
            names, seqs, max_mismatch=1
        )
        # Lossy collapse must not silently expand back to all input names
        assert len(labels) == matrix.shape[0]
        assert len(labels) < len(names)
        assert labels == [names[c[0]] for c in clusters]

    def test_every_input_appears_in_some_cluster(self, dup_fasta):
        names, seqs = load_sequences(str(dup_fasta))
        _, _, clusters = compute_similarity(names, seqs, max_mismatch=2)
        assert sorted(i for c in clusters for i in c) == list(range(len(names)))

    def test_clusters_tsv_records_membership(self, dup_fasta, temp_dir):
        names, seqs = load_sequences(str(dup_fasta))
        _, _, clusters = compute_similarity(names, seqs, max_mismatch=0)
        out = temp_dir / "clusters.tsv"
        write_clusters_tsv(clusters, names, str(out))

        lines = out.read_text().strip().split("\n")
        assert lines[0] == "cluster_id\trepresentative\tn_members\tmembers"
        assert len(lines) == 1 + len(clusters)

        rows = [line.split("\t") for line in lines[1:]]
        assert rows[0][1] == "a1"
        assert rows[0][2] == "2"
        assert rows[0][3] == "a1,a2"
        # Every input name is accounted for
        listed = [n for r in rows for n in r[3].split(",")]
        assert sorted(listed) == sorted(names)
