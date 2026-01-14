"""
Unit tests for get_trna_charging_cpm.py

Tests CPM calculation from per-read charging likelihood data.
"""

import gzip
import pandas as pd
import pytest

from get_trna_charging_cpm import per_read_charging


class TestPerReadCharging:
    """Tests for per_read_charging function."""

    def test_basic_cpm_calculation(self, temp_dir):
        """Should correctly calculate CPM for charged/uncharged."""
        input_tsv = temp_dir / "input.tsv"
        output_tsv = temp_dir / "output.tsv"

        # Create input with 10 reads: 6 charged (>=200), 4 uncharged (<200)
        with open(input_tsv, "w") as f:
            f.write("read_id\ttRNA\tcharging_likelihood\n")
            # 6 charged reads for tRNA-Ala
            for i in range(6):
                f.write(f"read{i}\ttRNA-Ala-AGC-1-1\t{220 + i}\n")
            # 4 uncharged reads for tRNA-Ala
            for i in range(4):
                f.write(f"read{i+6}\ttRNA-Ala-AGC-1-1\t{100 + i}\n")

        per_read_charging(str(input_tsv), str(output_tsv), threshold=200)

        # Read output
        df = pd.read_csv(output_tsv, sep="\t", index_col=0)

        assert "tRNA-Ala-AGC-1-1" in df.index
        row = df.loc["tRNA-Ala-AGC-1-1"]

        # Check counts
        assert row["counts_charged"] == 6
        assert row["counts_uncharged"] == 4

        # Check CPM (should be counts / 10 * 1e6)
        assert row["cpm_charged"] == 600000  # 6/10 * 1e6
        assert row["cpm_uncharged"] == 400000  # 4/10 * 1e6

    def test_threshold_boundary(self, temp_dir):
        """Reads at exactly threshold should be classified as charged."""
        input_tsv = temp_dir / "input.tsv"
        output_tsv = temp_dir / "output.tsv"

        with open(input_tsv, "w") as f:
            f.write("read_id\ttRNA\tcharging_likelihood\n")
            f.write("read1\ttRNA-Ala-AGC-1-1\t200\n")  # Exactly at threshold
            f.write("read2\ttRNA-Ala-AGC-1-1\t199\n")  # Below threshold

        per_read_charging(str(input_tsv), str(output_tsv), threshold=200)

        df = pd.read_csv(output_tsv, sep="\t", index_col=0)
        row = df.loc["tRNA-Ala-AGC-1-1"]

        assert row["counts_charged"] == 1  # 200 is >= 200
        assert row["counts_uncharged"] == 1  # 199 is < 200

    def test_multiple_trnas(self, temp_dir):
        """Should handle multiple different tRNAs."""
        input_tsv = temp_dir / "input.tsv"
        output_tsv = temp_dir / "output.tsv"

        with open(input_tsv, "w") as f:
            f.write("read_id\ttRNA\tcharging_likelihood\n")
            # tRNA-Ala: 2 charged, 1 uncharged
            f.write("read1\ttRNA-Ala-AGC-1-1\t220\n")
            f.write("read2\ttRNA-Ala-AGC-1-1\t230\n")
            f.write("read3\ttRNA-Ala-AGC-1-1\t100\n")
            # tRNA-Gly: 1 charged, 2 uncharged
            f.write("read4\ttRNA-Gly-GCC-1-1\t250\n")
            f.write("read5\ttRNA-Gly-GCC-1-1\t50\n")
            f.write("read6\ttRNA-Gly-GCC-1-1\t150\n")

        per_read_charging(str(input_tsv), str(output_tsv), threshold=200)

        df = pd.read_csv(output_tsv, sep="\t", index_col=0)

        # Check tRNA-Ala
        ala = df.loc["tRNA-Ala-AGC-1-1"]
        assert ala["counts_charged"] == 2
        assert ala["counts_uncharged"] == 1

        # Check tRNA-Gly
        gly = df.loc["tRNA-Gly-GCC-1-1"]
        assert gly["counts_charged"] == 1
        assert gly["counts_uncharged"] == 2

    def test_gzip_output(self, temp_dir):
        """Should write gzipped output when filename ends in .gz."""
        input_tsv = temp_dir / "input.tsv"
        output_tsv = temp_dir / "output.tsv.gz"

        with open(input_tsv, "w") as f:
            f.write("read_id\ttRNA\tcharging_likelihood\n")
            f.write("read1\ttRNA-Ala-AGC-1-1\t220\n")
            f.write("read2\ttRNA-Ala-AGC-1-1\t100\n")

        per_read_charging(str(input_tsv), str(output_tsv), threshold=200)

        # Should be able to read gzipped file
        df = pd.read_csv(output_tsv, sep="\t", index_col=0)
        assert "tRNA-Ala-AGC-1-1" in df.index

    def test_custom_threshold(self, temp_dir):
        """Should respect custom threshold parameter."""
        input_tsv = temp_dir / "input.tsv"
        output_tsv = temp_dir / "output.tsv"

        with open(input_tsv, "w") as f:
            f.write("read_id\ttRNA\tcharging_likelihood\n")
            f.write("read1\ttRNA-Ala-AGC-1-1\t150\n")
            f.write("read2\ttRNA-Ala-AGC-1-1\t100\n")

        # With threshold=150, first read should be charged
        per_read_charging(str(input_tsv), str(output_tsv), threshold=150)

        df = pd.read_csv(output_tsv, sep="\t", index_col=0)
        row = df.loc["tRNA-Ala-AGC-1-1"]

        assert row["counts_charged"] == 1
        assert row["counts_uncharged"] == 1

    def test_all_charged(self, temp_dir):
        """Should handle case where all reads are charged."""
        input_tsv = temp_dir / "input.tsv"
        output_tsv = temp_dir / "output.tsv"

        with open(input_tsv, "w") as f:
            f.write("read_id\ttRNA\tcharging_likelihood\n")
            f.write("read1\ttRNA-Ala-AGC-1-1\t255\n")
            f.write("read2\ttRNA-Ala-AGC-1-1\t220\n")

        per_read_charging(str(input_tsv), str(output_tsv), threshold=200)

        df = pd.read_csv(output_tsv, sep="\t", index_col=0)
        row = df.loc["tRNA-Ala-AGC-1-1"]

        assert row["counts_charged"] == 2
        assert row["counts_uncharged"] == 0

    def test_all_uncharged(self, temp_dir):
        """Should handle case where all reads are uncharged."""
        input_tsv = temp_dir / "input.tsv"
        output_tsv = temp_dir / "output.tsv"

        with open(input_tsv, "w") as f:
            f.write("read_id\ttRNA\tcharging_likelihood\n")
            f.write("read1\ttRNA-Ala-AGC-1-1\t50\n")
            f.write("read2\ttRNA-Ala-AGC-1-1\t100\n")

        per_read_charging(str(input_tsv), str(output_tsv), threshold=200)

        df = pd.read_csv(output_tsv, sep="\t", index_col=0)
        row = df.loc["tRNA-Ala-AGC-1-1"]

        assert row["counts_uncharged"] == 2
        assert row["counts_charged"] == 0
