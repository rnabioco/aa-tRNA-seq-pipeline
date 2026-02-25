"""Unit tests for compute_odds_ratios.py."""

import numpy as np
import pandas as pd
import pytest

from compute_odds_ratios import compute_or_for_pair, load_charging, read_fasta_lengths


class TestComputeOrForPair:
    def test_perfect_positive_correlation(self):
        """When both columns are identical, OR should be very high (with Haldane correction)."""
        col = pd.Series([1, 1, 1, 0, 0, 0])
        result = compute_or_for_pair(col, col)
        assert result is not None
        assert result["odds_ratio"] > 10  # Strong positive
        assert result["n11"] == 3
        assert result["n00"] == 3
        assert result["n10"] == 0
        assert result["n01"] == 0

    def test_no_correlation(self):
        """Independent columns should give OR ≈ 1."""
        np.random.seed(42)
        n = 1000
        col_i = pd.Series(np.random.binomial(1, 0.5, n))
        col_j = pd.Series(np.random.binomial(1, 0.5, n))
        result = compute_or_for_pair(col_i, col_j)
        assert result is not None
        assert 0.5 < result["odds_ratio"] < 2.0

    def test_negative_correlation(self):
        """Anti-correlated columns should give OR < 1."""
        col_i = pd.Series([1, 1, 1, 0, 0, 0, 1, 0])
        col_j = pd.Series([0, 0, 0, 1, 1, 1, 0, 1])
        result = compute_or_for_pair(col_i, col_j)
        assert result is not None
        assert result["odds_ratio"] < 1.0

    def test_haldane_correction_zero_cell(self):
        """When a cell is zero, Haldane +0.5 correction applied."""
        col_i = pd.Series([1, 1, 1, 0, 0])
        col_j = pd.Series([1, 1, 1, 0, 0])
        result = compute_or_for_pair(col_i, col_j)
        # n10 and n01 are 0, so Haldane should kick in
        assert result["n10"] == 0
        assert result["n01"] == 0
        assert np.isfinite(result["odds_ratio"])

    def test_all_same_values(self):
        """All-same column should still return result (all ones or all zeros)."""
        col_i = pd.Series([1, 1, 1, 1])
        col_j = pd.Series([0, 1, 0, 1])
        result = compute_or_for_pair(col_i, col_j)
        assert result is not None
        assert result["total_obs"] == 4

    def test_empty_input(self):
        """Empty series should return None."""
        col_i = pd.Series([], dtype=float)
        col_j = pd.Series([], dtype=float)
        result = compute_or_for_pair(col_i, col_j)
        assert result is None

    def test_nan_handling(self):
        """NaN values should be dropped before computation."""
        col_i = pd.Series([1, 0, np.nan, 1, 0])
        col_j = pd.Series([1, 0, 1, np.nan, 0])
        result = compute_or_for_pair(col_i, col_j)
        assert result is not None
        assert result["total_obs"] == 3  # Only rows where both are non-NaN

    def test_known_2x2_table(self):
        """Verify against a known 2×2 table: [[10,5],[3,12]]."""
        # Build columns: 10 (1,1), 5 (1,0), 3 (0,1), 12 (0,0)
        col_i = pd.Series([1]*10 + [1]*5 + [0]*3 + [0]*12)
        col_j = pd.Series([1]*10 + [0]*5 + [1]*3 + [0]*12)
        result = compute_or_for_pair(col_i, col_j)
        assert result["n11"] == 10
        assert result["n10"] == 5
        assert result["n01"] == 3
        assert result["n00"] == 12
        expected_or = (10 * 12) / (5 * 3)
        assert abs(result["odds_ratio"] - expected_or) < 0.01

    def test_output_dict_keys(self):
        """Result should have all expected keys."""
        col_i = pd.Series([1, 0, 1, 0])
        col_j = pd.Series([0, 1, 1, 0])
        result = compute_or_for_pair(col_i, col_j)
        expected_keys = {
            "n00", "n01", "n10", "n11", "total_obs",
            "odds_ratio", "log_odds_ratio", "se_log_or",
            "ci_lower", "ci_upper", "fisher_or", "p_value",
        }
        assert set(result.keys()) == expected_keys

    def test_ci_contains_or(self):
        """Confidence interval should contain the odds ratio."""
        col_i = pd.Series([1]*20 + [0]*20 + [1]*5 + [0]*5)
        col_j = pd.Series([1]*20 + [0]*5 + [0]*20 + [1]*5)
        result = compute_or_for_pair(col_i, col_j)
        assert result["ci_lower"] <= result["odds_ratio"] <= result["ci_upper"]


class TestLoadCharging:
    def test_binarize_above_threshold(self, temp_dir):
        df = pd.DataFrame({
            "read_id": ["r1", "r2", "r3"],
            "charging_likelihood": [250, 100, 200],
        })
        path = temp_dir / "charging.tsv"
        df.to_csv(path, sep="\t", index=False)
        result = load_charging(str(path), ml_threshold=200)
        assert list(result["charged"]) == [1, 0, 1]

    def test_output_columns(self, temp_dir):
        df = pd.DataFrame({
            "read_id": ["r1"],
            "charging_likelihood": [150],
        })
        path = temp_dir / "charging.tsv"
        df.to_csv(path, sep="\t", index=False)
        result = load_charging(str(path), ml_threshold=200)
        assert list(result.columns) == ["read_id", "charged"]


class TestReadFastaLengths:
    def test_single_sequence(self, temp_dir):
        fa = temp_dir / "test.fa"
        fa.write_text(">seq1\nACGTACGT\n")
        result = read_fasta_lengths(str(fa))
        assert result == {"seq1": 8}

    def test_multiple_sequences(self, temp_dir):
        fa = temp_dir / "test.fa"
        fa.write_text(">seq1\nACGT\n>seq2\nACGTACGTAC\n")
        result = read_fasta_lengths(str(fa))
        assert result == {"seq1": 4, "seq2": 10}

    def test_multiline_sequence(self, temp_dir):
        fa = temp_dir / "test.fa"
        fa.write_text(">seq1\nACGT\nACGT\n")
        result = read_fasta_lengths(str(fa))
        assert result == {"seq1": 8}

    def test_empty_file(self, temp_dir):
        fa = temp_dir / "test.fa"
        fa.write_text("")
        result = read_fasta_lengths(str(fa))
        assert result == {}
