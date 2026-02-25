"""Unit tests for compute_seq_similarity.py."""

import numpy as np
import pysam
import pytest

from compute_seq_similarity import compute_similarity_matrix, write_matrix_tsv


@pytest.fixture
def identical_fasta(temp_dir):
    """FASTA with two identical sequences."""
    fa = temp_dir / "identical.fa"
    fa.write_text(">seq1\nACGTACGT\n>seq2\nACGTACGT\n")
    pysam.faidx(str(fa))
    return fa


@pytest.fixture
def different_fasta(temp_dir):
    """FASTA with two very different sequences."""
    fa = temp_dir / "different.fa"
    fa.write_text(">seq1\nAAAAAAAA\n>seq2\nTTTTTTTT\n")
    pysam.faidx(str(fa))
    return fa


@pytest.fixture
def single_fasta(temp_dir):
    """FASTA with one sequence."""
    fa = temp_dir / "single.fa"
    fa.write_text(">seq1\nACGTACGT\n")
    pysam.faidx(str(fa))
    return fa


class TestComputeSimilarityMatrix:
    def test_identical_sequences(self, identical_fasta):
        matrix, names = compute_similarity_matrix(str(identical_fasta))
        assert names == ["seq1", "seq2"]
        assert matrix[0, 1] == 100.0
        assert matrix[1, 0] == 100.0

    def test_different_sequences(self, different_fasta):
        matrix, names = compute_similarity_matrix(str(different_fasta))
        assert matrix[0, 1] < 50.0

    def test_single_sequence(self, single_fasta):
        matrix, names = compute_similarity_matrix(str(single_fasta))
        assert matrix.shape == (1, 1)
        assert matrix[0, 0] == 100.0

    def test_symmetric(self, identical_fasta):
        matrix, _ = compute_similarity_matrix(str(identical_fasta))
        np.testing.assert_array_equal(matrix, matrix.T)

    def test_diagonal_100(self, different_fasta):
        matrix, _ = compute_similarity_matrix(str(different_fasta))
        for i in range(matrix.shape[0]):
            assert matrix[i, i] == 100.0


class TestWriteMatrixTsv:
    def test_output_format(self, temp_dir):
        matrix = np.array([[100.0, 50.0], [50.0, 100.0]])
        names = ["seq1", "seq2"]
        out_path = temp_dir / "matrix.tsv"
        write_matrix_tsv(matrix, names, str(out_path))
        lines = out_path.read_text().strip().split("\n")
        assert len(lines) == 3  # header + 2 data rows
        assert "seq1" in lines[0]
        assert "seq2" in lines[0]

    def test_round_trip_values(self, temp_dir):
        matrix = np.array([[100.0, 75.55], [75.55, 100.0]])
        names = ["a", "b"]
        out_path = temp_dir / "matrix.tsv"
        write_matrix_tsv(matrix, names, str(out_path))
        lines = out_path.read_text().strip().split("\n")
        row1_vals = lines[1].split("\t")[1:]
        assert float(row1_vals[0]) == 100.0
        assert float(row1_vals[1]) == 75.55
