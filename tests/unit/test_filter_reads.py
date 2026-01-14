"""
Unit tests for filter_reads.py

Tests filtering logic for truncation detection, strand filtering,
mapping quality, and adapter edit distance.
"""

import argparse
import pysam
import pytest

from filter_reads import (
    FILTER_CODES,
    FILTER_TAG,
    FilterStats,
    count_adapter_edits,
    filter_bam,
)


class TestFilterStats:
    """Tests for FilterStats class."""

    def test_initial_state(self):
        """Fresh FilterStats should have zero counts."""
        stats = FilterStats()
        assert stats.n_align == 0
        assert stats.n_filtered == 0
        for reason in FILTER_CODES:
            assert stats.filter_counts[reason] == 0

    def test_log_passed_read(self):
        """Read with tag=0 should not increment filtered count."""
        stats = FilterStats()
        stats.log(0)

        assert stats.n_align == 1
        assert stats.n_filtered == 0

    def test_log_5p_truncation(self):
        """Read with 5p_trunc flag should be counted."""
        stats = FilterStats()
        stats.log(FILTER_CODES["5p_trunc"])

        assert stats.n_align == 1
        assert stats.n_filtered == 1
        assert stats.filter_counts["5p_trunc"] == 1

    def test_log_multiple_flags(self):
        """Read with multiple flags should increment all counters."""
        stats = FilterStats()
        tag = FILTER_CODES["5p_trunc"] | FILTER_CODES["low_mapq"]
        stats.log(tag)

        assert stats.n_align == 1
        assert stats.n_filtered == 1
        assert stats.filter_counts["5p_trunc"] == 1
        assert stats.filter_counts["low_mapq"] == 1

    def test_summary_format(self):
        """Summary should contain expected fields."""
        stats = FilterStats()
        stats.log(0)  # Passed
        stats.log(FILTER_CODES["unmapped"])  # Failed

        summary = stats.summary()

        assert "total_alignments 2" in summary
        assert "failed_alignments 1" in summary
        assert "passed_alignments 1" in summary


class TestCountAdapterEdits:
    """Tests for count_adapter_edits function."""

    def test_perfect_alignment(self, temp_dir):
        """Perfect alignment in adapter region should return 0 edits."""
        bam_path = temp_dir / "test.bam"

        # Reference: 100bp, adapter region 78-100 (last 22bp)
        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref", "LN": 100}]}

        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 100)]  # 100M
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            # Add MD tag for perfect match (100 matches)
            read.set_tag("MD", "100")
            outf.write(read)

        with pysam.AlignmentFile(str(bam_path), "rb") as inbam:
            read = next(inbam)
            # Adapter region from position 78 to 100
            edits = count_adapter_edits(read, 78, 100)

        assert edits is not None
        assert edits == 0

    def test_alignment_not_reaching_adapter(self, temp_dir):
        """Alignment not reaching adapter region should return None."""
        bam_path = temp_dir / "test.bam"

        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref", "LN": 100}]}

        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A" * 50
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 50)]  # Only aligns to first 50bp
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            read.set_tag("MD", "50")
            outf.write(read)

        with pysam.AlignmentFile(str(bam_path), "rb") as inbam:
            read = next(inbam)
            # Adapter region from position 78 to 100
            edits = count_adapter_edits(read, 78, 100)

        assert edits is None


class TestFilterBam:
    """Integration tests for filter_bam function."""

    def test_filter_unmapped_reads(self, temp_dir):
        """Unmapped reads should be filtered."""
        input_bam = temp_dir / "input.bam"
        output_bam = temp_dir / "output.bam"

        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref", "LN": 100}]}

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            # Mapped read
            read1 = pysam.AlignedSegment()
            read1.query_name = "mapped"
            read1.query_sequence = "A" * 100
            read1.flag = 0
            read1.reference_id = 0
            read1.reference_start = 0
            read1.cigartuples = [(0, 100)]
            read1.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read1)

            # Unmapped read
            read2 = pysam.AlignedSegment()
            read2.query_name = "unmapped"
            read2.query_sequence = "A" * 50
            read2.flag = 4  # Unmapped flag
            read2.reference_id = -1
            read2.reference_start = -1
            read2.query_qualities = pysam.qualitystring_to_array("I" * 50)
            outf.write(read2)

        args = argparse.Namespace(
            input_bam=str(input_bam),
            output_bam=str(output_bam),
            five_p_truncation=-1,  # Disable
            three_p_truncation=-1,  # Disable
            min_mapq=0,
            only_positive=False,
            rescue_multi_mappers=False,
            trna_table=None,
            max_edit_dist=None,
            failed_bam=None,
        )

        filter_bam(args)

        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            reads = list(bam)

        # Only mapped read should be in output
        assert len(reads) == 1
        assert reads[0].query_name == "mapped"

    def test_filter_negative_strand(self, temp_dir):
        """Negative strand reads should be filtered when only_positive=True."""
        input_bam = temp_dir / "input.bam"
        output_bam = temp_dir / "output.bam"

        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref", "LN": 100}]}

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            # Positive strand read
            read1 = pysam.AlignedSegment()
            read1.query_name = "positive"
            read1.query_sequence = "A" * 100
            read1.flag = 0
            read1.reference_id = 0
            read1.reference_start = 0
            read1.cigartuples = [(0, 100)]
            read1.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read1)

            # Negative strand read
            read2 = pysam.AlignedSegment()
            read2.query_name = "negative"
            read2.query_sequence = "A" * 100
            read2.flag = 16  # Reverse strand flag
            read2.reference_id = 0
            read2.reference_start = 0
            read2.cigartuples = [(0, 100)]
            read2.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read2)

        args = argparse.Namespace(
            input_bam=str(input_bam),
            output_bam=str(output_bam),
            five_p_truncation=-1,
            three_p_truncation=-1,
            min_mapq=0,
            only_positive=True,  # Filter negative strand
            rescue_multi_mappers=False,
            trna_table=None,
            max_edit_dist=None,
            failed_bam=None,
        )

        filter_bam(args)

        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            reads = list(bam)

        # Only positive strand read should be in output
        assert len(reads) == 1
        assert reads[0].query_name == "positive"

    def test_filter_5p_truncation(self, temp_dir):
        """Reads with excessive 5' truncation should be filtered."""
        input_bam = temp_dir / "input.bam"
        output_bam = temp_dir / "output.bam"

        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref", "LN": 100}]}

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            # Good read - starts at position 0
            read1 = pysam.AlignedSegment()
            read1.query_name = "full_length"
            read1.query_sequence = "A" * 100
            read1.flag = 0
            read1.reference_id = 0
            read1.reference_start = 0
            read1.cigartuples = [(0, 100)]
            read1.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read1)

            # Truncated read - starts at position 30
            read2 = pysam.AlignedSegment()
            read2.query_name = "truncated"
            read2.query_sequence = "A" * 70
            read2.flag = 0
            read2.reference_id = 0
            read2.reference_start = 30  # Beyond 5' truncation threshold
            read2.cigartuples = [(0, 70)]
            read2.query_qualities = pysam.qualitystring_to_array("I" * 70)
            outf.write(read2)

        args = argparse.Namespace(
            input_bam=str(input_bam),
            output_bam=str(output_bam),
            five_p_truncation=24,  # Allow up to 24bp truncation
            three_p_truncation=-1,  # Disable
            min_mapq=0,
            only_positive=False,
            rescue_multi_mappers=False,
            trna_table=None,
            max_edit_dist=None,
            failed_bam=None,
        )

        filter_bam(args)

        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            reads = list(bam)

        # Only full length read should pass
        assert len(reads) == 1
        assert reads[0].query_name == "full_length"

    def test_adds_zf_tag(self, temp_dir):
        """All reads should get the zf filter tag."""
        input_bam = temp_dir / "input.bam"
        output_bam = temp_dir / "output.bam"
        failed_bam = temp_dir / "failed.bam"

        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref", "LN": 100}]}

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "test"
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read)

        args = argparse.Namespace(
            input_bam=str(input_bam),
            output_bam=str(output_bam),
            five_p_truncation=-1,
            three_p_truncation=-1,
            min_mapq=0,
            only_positive=False,
            rescue_multi_mappers=False,
            trna_table=None,
            max_edit_dist=None,
            failed_bam=str(failed_bam),
        )

        filter_bam(args)

        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            read = next(bam)
            assert read.has_tag(FILTER_TAG)
            assert read.get_tag(FILTER_TAG) == 0  # Passed

    def test_low_mapq_filter(self, temp_dir):
        """Reads below min MAPQ should be filtered."""
        input_bam = temp_dir / "input.bam"
        output_bam = temp_dir / "output.bam"

        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref", "LN": 100}]}

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            # High MAPQ read
            read1 = pysam.AlignedSegment()
            read1.query_name = "high_mapq"
            read1.query_sequence = "A" * 100
            read1.flag = 0
            read1.reference_id = 0
            read1.reference_start = 0
            read1.mapping_quality = 60
            read1.cigartuples = [(0, 100)]
            read1.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read1)

            # Low MAPQ read
            read2 = pysam.AlignedSegment()
            read2.query_name = "low_mapq"
            read2.query_sequence = "A" * 100
            read2.flag = 0
            read2.reference_id = 0
            read2.reference_start = 0
            read2.mapping_quality = 5
            read2.cigartuples = [(0, 100)]
            read2.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read2)

        args = argparse.Namespace(
            input_bam=str(input_bam),
            output_bam=str(output_bam),
            five_p_truncation=-1,
            three_p_truncation=-1,
            min_mapq=20,  # Require MAPQ >= 20
            only_positive=False,
            rescue_multi_mappers=False,
            trna_table=None,
            max_edit_dist=None,
            failed_bam=None,
        )

        filter_bam(args)

        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            reads = list(bam)

        assert len(reads) == 1
        assert reads[0].query_name == "high_mapq"
