"""Unit tests for get_align_stats.py."""

import pysam
import pytest

from conftest import create_bam_with_reads
from get_align_stats import CountArray, ReadStats, get_read_stats


class TestCountArray:
    def test_push(self):
        ca = CountArray(max=10)
        ca.push(3)
        ca.push(3)
        ca.push(5)
        assert ca.counts[3] == 2
        assert ca.counts[5] == 1
        assert ca.total == 3

    def test_push_overflow(self):
        """Pushing beyond max should be silently ignored."""
        ca = CountArray(max=5)
        ca.push(10)
        assert ca.total == 0

    def test_nth_single(self):
        ca = CountArray(max=10)
        ca.push(5)
        assert ca.nth(0) == 5

    def test_nth_multiple(self):
        ca = CountArray(max=10)
        ca.push(2)
        ca.push(2)
        ca.push(7)
        assert ca.nth(0) == 2
        assert ca.nth(1) == 2
        assert ca.nth(2) == 7

    def test_nth_out_of_range(self):
        ca = CountArray(max=10)
        ca.push(3)
        assert ca.nth(5) is None

    def test_mean_empty(self):
        ca = CountArray(max=10)
        assert ca.mean() == 0

    def test_mean_single(self):
        ca = CountArray(max=10)
        ca.push(4)
        assert ca.mean() == 4.0

    def test_mean_multiple(self):
        ca = CountArray(max=10)
        ca.push(2)
        ca.push(4)
        ca.push(6)
        assert ca.mean() == pytest.approx(4.0)

    def test_quantile_odd(self):
        """Median of [2, 4, 6] = 4."""
        ca = CountArray(max=10)
        ca.push(2)
        ca.push(4)
        ca.push(6)
        assert ca.quantile(0.5) == 4.0

    def test_quantile_even(self):
        """Median of [2, 4] = 3.0."""
        ca = CountArray(max=10)
        ca.push(2)
        ca.push(4)
        assert ca.quantile(0.5) == 3.0

    def test_quantile_empty(self):
        ca = CountArray(max=10)
        assert ca.quantile(0.5) == 0

    def test_quantile_invalid_p(self):
        ca = CountArray(max=10)
        ca.push(1)
        with pytest.raises(ValueError):
            ca.quantile(1.5)


class TestReadStats:
    def test_initial_state(self):
        rs = ReadStats()
        assert rs.mapped_reads == 0
        assert rs.n_pos_reads == 0

    def test_summary_keys(self):
        rs = ReadStats()
        summary = rs.summary()
        expected = {
            "mapped_reads", "pos_reads", "mapq0_reads", "mapq_pass_reads",
            "mean_length", "median_length", "mean_base_quality",
            "median_base_quality", "mean_MAPQ", "median_MAPQ",
        }
        assert set(summary.keys()) == expected


class TestGetReadStats:
    def test_basic_bam(self, temp_dir):
        """Test with mapped, unmapped, and reverse reads."""
        bam_path = temp_dir / "test.bam"
        reads = [
            {"name": "r1", "seq": "ACGTACGT", "flag": 0, "mapq": 60},  # mapped fwd
            {"name": "r2", "seq": "ACGTACGT", "flag": 16, "mapq": 30},  # mapped rev
            {"name": "r3", "seq": "ACGTACGT", "flag": 4},  # unmapped
        ]
        create_bam_with_reads(bam_path, reads)
        result = get_read_stats(str(bam_path))
        assert result["n_reads"] == 3
        assert result["mapped_reads"] == 2

    def test_flag_filtering(self, temp_dir):
        """When flag is set, only matching reads should be counted."""
        bam_path = temp_dir / "test.bam"
        reads = [
            {"name": "r1", "seq": "ACGTACGT", "flag": 0, "mapq": 60},
            {"name": "r2", "seq": "ACGTACGT", "flag": 16, "mapq": 30},
        ]
        create_bam_with_reads(bam_path, reads)
        # Flag 16 = reverse strand
        result = get_read_stats(str(bam_path), flag=16)
        assert result["n_reads"] == 1

    def test_qname_deduplication(self, temp_dir):
        """Duplicate qnames should only be counted once."""
        bam_path = temp_dir / "test.bam"
        reads = [
            {"name": "r1", "seq": "ACGTACGT", "flag": 0, "mapq": 60},
            {"name": "r1", "seq": "ACGTACGT", "flag": 256, "mapq": 60},  # secondary
        ]
        create_bam_with_reads(bam_path, reads)
        result = get_read_stats(str(bam_path))
        assert result["n_reads"] == 1

    def test_supplementary_records_are_not_counted(self, temp_dir):
        """A chimeric read gets a supplementary record; it is still one read."""
        bam_path = temp_dir / "test.bam"
        reads = [
            {"name": "r1", "seq": "ACGTACGT", "flag": 0, "mapq": 60},
            {"name": "r1", "seq": "ACGT", "flag": 2048, "mapq": 60},
        ]
        create_bam_with_reads(bam_path, reads)
        result = get_read_stats(str(bam_path))
        assert result["n_reads"] == 1

    def test_the_primary_record_supplies_the_stats(self, temp_dir):
        """
        File order must not decide which record a read is measured by.

        Counting the first record per query name -- which is what holding every
        seen name in a set amounted to -- measures whichever alignment the file
        happens to list first. Here that is a four-base supplementary at MAPQ
        10, while the read itself is eight bases at MAPQ 60.
        """
        bam_path = temp_dir / "test.bam"
        reads = [
            {"name": "r1", "seq": "ACGT", "flag": 2048, "mapq": 10},
            {"name": "r1", "seq": "ACGTACGT", "flag": 0, "mapq": 60},
        ]
        create_bam_with_reads(bam_path, reads)
        result = get_read_stats(str(bam_path))
        assert result["n_reads"] == 1
        assert result["mean_length"] == 8
        assert result["mean_MAPQ"] == 60
