"""
Unit tests for add_adapter_tags.py

Tests CIGAR parsing, adapter finding with parasail alignment,
PT tag formatting, and BAM processing.
"""

import tempfile
from array import array
from pathlib import Path

import pysam
import pytest

from add_adapter_tags import (
    DEFAULT_ADAPTER_3P,
    DEFAULT_ADAPTER_5P,
    DEFAULT_GAP_EXTEND,
    DEFAULT_GAP_OPEN,
    DEFAULT_MIN_SCORE_3P,
    DEFAULT_MIN_SCORE_5P,
    Stats,
    create_scoring_matrix,
    find_3p_adapter,
    find_5p_adapter,
    find_adapter,
    format_pt_tag,
    get_adapter_bounds_from_cigar,
    process_bam,
)


class TestGetAdapterBoundsFromCigar:
    """Tests for get_adapter_bounds_from_cigar function."""

    def test_simple_match(self):
        """Simple match without leading/trailing deletions."""
        # 23M - 23 matches
        cigar = b"23M"
        beg_ref = 0
        start, end = get_adapter_bounds_from_cigar(cigar, beg_ref)
        assert start == 0
        assert end == 23

    def test_leading_deletion(self):
        """Leading deletion should be skipped."""
        # 5D10M - 5 deletions then 10 matches
        cigar = b"5D10M"
        beg_ref = 0
        start, end = get_adapter_bounds_from_cigar(cigar, beg_ref)
        # Start should skip the 5D
        assert start == 5
        assert end == 15

    def test_trailing_deletion(self):
        """Trailing deletion should be excluded from end."""
        # 10M5D - 10 matches then 5 deletions
        cigar = b"10M5D"
        beg_ref = 0
        start, end = get_adapter_bounds_from_cigar(cigar, beg_ref)
        assert start == 0
        assert end == 10  # Excludes trailing D

    def test_both_leading_and_trailing_deletions(self):
        """Both leading and trailing deletions should be handled."""
        # 3D15M4D - 3 deletions, 15 matches, 4 deletions
        cigar = b"3D15M4D"
        beg_ref = 0
        start, end = get_adapter_bounds_from_cigar(cigar, beg_ref)
        assert start == 3  # Skip leading 3D
        assert end == 18  # 3 + 15 = 18, excludes trailing 4D

    def test_insertion_in_middle(self):
        """Insertions don't consume reference positions."""
        # 10M2I10M - 10 matches, 2 insertions, 10 matches
        cigar = b"10M2I10M"
        beg_ref = 0
        start, end = get_adapter_bounds_from_cigar(cigar, beg_ref)
        assert start == 0
        assert end == 20  # 10 + 10, insertions don't add

    def test_with_nonzero_beg_ref(self):
        """Non-zero beg_ref should be respected."""
        cigar = b"20M"
        beg_ref = 5
        start, end = get_adapter_bounds_from_cigar(cigar, beg_ref)
        assert start == 5
        assert end == 25


class TestFindAdapter:
    """Tests for find_adapter and related functions."""

    @pytest.fixture
    def scoring_matrix(self):
        return create_scoring_matrix()

    def test_perfect_match_5p(self, scoring_matrix):
        """Perfect match of 5' adapter should return correct positions."""
        # Read starts with exact 5' adapter
        read_seq = DEFAULT_ADAPTER_5P + "GCGGCTATAGCTCAGTTGGTA" + "CCA"

        result = find_5p_adapter(
            read_seq,
            DEFAULT_ADAPTER_5P,
            scoring_matrix,
            DEFAULT_GAP_OPEN,
            DEFAULT_GAP_EXTEND,
            min_score=20,  # Use lower threshold for testing
        )

        assert result is not None
        start, end, score = result
        # Start should be near 0 (may have small offset due to alignment)
        assert start <= 5
        assert score > 20

    def test_perfect_match_3p(self, scoring_matrix):
        """Perfect match of 3' adapter should return correct positions."""
        read_seq = "GCGGCTATAGCTCAGTTGGTACCA" + DEFAULT_ADAPTER_3P

        result = find_3p_adapter(
            read_seq,
            DEFAULT_ADAPTER_3P,
            scoring_matrix,
            DEFAULT_GAP_OPEN,
            DEFAULT_GAP_EXTEND,
            DEFAULT_MIN_SCORE_3P,
        )

        assert result is not None
        start, end, score = result
        # 3' adapter should be at the end
        assert start == 24  # After tRNA sequence
        assert end == len(read_seq)
        assert score > DEFAULT_MIN_SCORE_3P

    def test_adapter_with_mismatches(self, scoring_matrix):
        """Adapter with ~10% mismatches should still be found."""
        # Introduce 2 mismatches in 23bp adapter (~9% error rate)
        adapter_mutated = "CCTAAGAGCAAGTAGAAGCCTGG"  # Changed AA to TA at pos 13-14
        read_seq = adapter_mutated + "GCGGCTATAGCTCAGTTGGTA"

        result = find_5p_adapter(
            read_seq,
            DEFAULT_ADAPTER_5P,
            scoring_matrix,
            DEFAULT_GAP_OPEN,
            DEFAULT_GAP_EXTEND,
            min_score=20,  # Use lower threshold for testing
        )

        assert result is not None
        start, end, score = result
        # Start should be near 0
        assert start <= 5
        # Score should still be reasonable
        assert score > 20

    def test_adapter_too_many_mismatches(self, scoring_matrix):
        """Adapter with >20% mismatches should not be found."""
        # Introduce 5 mismatches in 23bp adapter (~22% error rate)
        adapter_bad = "CCTTTGAGCTTGTAGTGGCCTGG"
        read_seq = adapter_bad + "GCGGCTATAGCTCAGTTGGTA"

        result = find_5p_adapter(
            read_seq,
            DEFAULT_ADAPTER_5P,
            scoring_matrix,
            DEFAULT_GAP_OPEN,
            DEFAULT_GAP_EXTEND,
            DEFAULT_MIN_SCORE_5P,
        )

        # Should return None due to low score
        assert result is None

    def test_adapter_not_present(self, scoring_matrix):
        """Random sequence should not match adapter."""
        read_seq = "ACGTACGTACGTACGTACGTACGTACGT"

        result = find_5p_adapter(
            read_seq,
            DEFAULT_ADAPTER_5P,
            scoring_matrix,
            DEFAULT_GAP_OPEN,
            DEFAULT_GAP_EXTEND,
            DEFAULT_MIN_SCORE_5P,
        )

        assert result is None

    def test_short_sequence(self, scoring_matrix):
        """Very short sequence should return None."""
        read_seq = "ACGT"

        result = find_5p_adapter(
            read_seq,
            DEFAULT_ADAPTER_5P,
            scoring_matrix,
            DEFAULT_GAP_OPEN,
            DEFAULT_GAP_EXTEND,
            DEFAULT_MIN_SCORE_5P,
        )

        assert result is None

    def test_partial_3p_adapter_match(self, scoring_matrix):
        """Partial 3' adapter (truncated polyA tail) should still match."""
        # 3' adapter without full polyA tail
        partial_3p = DEFAULT_ADAPTER_3P[:30]  # First 30 bp
        read_seq = "GCGGCTATAGCTCAGTTGGTACCA" + partial_3p

        result = find_3p_adapter(
            read_seq,
            DEFAULT_ADAPTER_3P,
            scoring_matrix,
            DEFAULT_GAP_OPEN,
            DEFAULT_GAP_EXTEND,
            DEFAULT_MIN_SCORE_3P,
        )

        # Should still find adapter with lower threshold
        assert result is not None

    def test_adapter_with_polyt_prefix(self, scoring_matrix):
        """5' adapter after poly-T prefix should be found (regression test for sg_dx fix).

        Real nanopore reads often have poly-T at the start followed by the adapter.
        This test ensures the adapter can be found anywhere within the search region.
        """
        # Real read pattern: poly-T + adapter + tRNA
        read_seq = "TTTTTTTTTTTT" + DEFAULT_ADAPTER_5P + "GCGGCTATAGCTCAGTTGGTA"

        result = find_5p_adapter(
            read_seq,
            DEFAULT_ADAPTER_5P,
            scoring_matrix,
            DEFAULT_GAP_OPEN,
            DEFAULT_GAP_EXTEND,
            DEFAULT_MIN_SCORE_5P,
        )

        assert result is not None
        start, end, score = result
        # Adapter should be found after the poly-T prefix
        assert start >= 10  # After poly-T
        assert start <= 15
        assert score >= DEFAULT_MIN_SCORE_5P

    def test_adapter_with_errors_embedded(self, scoring_matrix):
        """Adapter with errors after poly-T should still be found.

        Tests the realistic case from actual data where adapter has ~10% errors.
        """
        # Pattern from real data: poly-T + adapter with errors
        # Real: CCTAAGAGCAAGGGGAAGCCTGG vs config: CCTAAGAGCAAGAAGAAGCCTGG
        adapter_with_errors = "CCTAAGAGCAAGGGGAAGCCTGG"
        read_seq = "TTTTTTTTTTTT" + adapter_with_errors + "TGGAGGATGCGGGCAGCGAGTCCCG"

        result = find_5p_adapter(
            read_seq,
            DEFAULT_ADAPTER_5P,
            scoring_matrix,
            DEFAULT_GAP_OPEN,
            DEFAULT_GAP_EXTEND,
            min_score=20,  # Lower threshold due to errors
        )

        assert result is not None
        start, end, score = result
        assert start >= 10  # After poly-T
        assert score >= 20


class TestFormatPtTag:
    """Tests for format_pt_tag function."""

    def test_both_adapters(self):
        """Both adapters present should create proper PT tag."""
        result_5p = (0, 23, 46)  # start, end, score
        result_3p = (100, 140, 80)

        pt_tag = format_pt_tag(result_5p, result_3p)

        assert pt_tag == "0;23;+;5p_adapter|100;140;+;3p_adapter"

    def test_only_5p_adapter(self):
        """Only 5' adapter should create single annotation."""
        result_5p = (0, 23, 46)
        result_3p = None

        pt_tag = format_pt_tag(result_5p, result_3p)

        assert pt_tag == "0;23;+;5p_adapter"

    def test_only_3p_adapter(self):
        """Only 3' adapter should create single annotation."""
        result_5p = None
        result_3p = (100, 140, 80)

        pt_tag = format_pt_tag(result_5p, result_3p)

        assert pt_tag == "100;140;+;3p_adapter"

    def test_neither_adapter(self):
        """No adapters should return None."""
        result_5p = None
        result_3p = None

        pt_tag = format_pt_tag(result_5p, result_3p)

        assert pt_tag is None


class TestStats:
    """Tests for Stats class."""

    def test_initial_counts(self):
        """Fresh Stats should have zero counts."""
        stats = Stats()
        assert stats.total == 0
        assert stats.with_5p == 0
        assert stats.with_3p == 0
        assert stats.with_both == 0
        assert stats.with_neither == 0

    def test_update_with_5p(self):
        """Update with only 5' adapter."""
        stats = Stats()
        stats.update(has_5p=True, has_3p=False)

        assert stats.total == 1
        assert stats.with_5p == 1
        assert stats.with_3p == 0
        assert stats.with_both == 0
        assert stats.with_neither == 0

    def test_update_with_3p(self):
        """Update with only 3' adapter."""
        stats = Stats()
        stats.update(has_5p=False, has_3p=True)

        assert stats.total == 1
        assert stats.with_5p == 0
        assert stats.with_3p == 1
        assert stats.with_both == 0
        assert stats.with_neither == 0

    def test_update_with_both(self):
        """Update with both adapters."""
        stats = Stats()
        stats.update(has_5p=True, has_3p=True)

        assert stats.total == 1
        assert stats.with_5p == 1
        assert stats.with_3p == 1
        assert stats.with_both == 1
        assert stats.with_neither == 0

    def test_update_with_neither(self):
        """Update with no adapters."""
        stats = Stats()
        stats.update(has_5p=False, has_3p=False)

        assert stats.total == 1
        assert stats.with_5p == 0
        assert stats.with_3p == 0
        assert stats.with_both == 0
        assert stats.with_neither == 1

    def test_summary_format(self):
        """Summary should have expected format."""
        stats = Stats()
        stats.update(has_5p=True, has_3p=True)
        stats.update(has_5p=True, has_3p=False)

        summary = stats.summary()
        assert "total_reads 2" in summary
        assert "with_5p_adapter 2" in summary
        assert "with_3p_adapter 1" in summary
        assert "with_both_adapters 1" in summary
        assert "with_neither_adapter 0" in summary


class TestProcessBam:
    """Integration tests for process_bam function."""

    def test_adds_pt_tags(self, sample_bam, temp_dir):
        """process_bam should add PT tags to reads with adapters."""
        output_bam = temp_dir / "output.bam"

        stats = process_bam(
            str(sample_bam),
            str(output_bam),
            DEFAULT_ADAPTER_5P,
            [("default", DEFAULT_ADAPTER_3P)],
            min_score_5p=20,  # Lower threshold for test data
            min_score_3p=20,
            match=2,
            mismatch=-1,
            gap_open=2,
            gap_extend=1,
        )

        # Check output BAM
        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            reads = list(bam)

        # First read has adapters - should have PT tag
        read1 = [r for r in reads if r.query_name == "read1"][0]
        assert read1.has_tag("PT")
        pt_tag = read1.get_tag("PT")
        # At least one adapter should be found
        assert "adapter" in pt_tag

        # Verify stats
        assert stats.total >= 2

    def test_handles_unmapped_reads(self, temp_dir):
        """Unmapped reads without sequence should pass through."""
        input_bam = temp_dir / "input.bam"
        output_bam = temp_dir / "output.bam"

        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref", "LN": 100}]}

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "unmapped"
            read.flag = 4  # Unmapped
            read.reference_id = -1
            read.reference_start = -1
            # No sequence set
            outf.write(read)

        stats = process_bam(
            str(input_bam),
            str(output_bam),
            DEFAULT_ADAPTER_5P,
            [("default", DEFAULT_ADAPTER_3P)],
            DEFAULT_MIN_SCORE_5P,
            DEFAULT_MIN_SCORE_3P,
            2,
            -1,
            2,
            1,
        )

        # Should complete without error
        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            reads = list(bam)
        assert len(reads) == 1

    def test_preserves_existing_tags(self, temp_dir):
        """Existing tags should be preserved."""
        input_bam = temp_dir / "input.bam"
        output_bam = temp_dir / "output.bam"

        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref", "LN": 200}]}

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "with_tags"
            read.query_sequence = (
                DEFAULT_ADAPTER_5P + "ACGTACGTACGTACGTACGTCCA" + DEFAULT_ADAPTER_3P
            )
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, len(read.query_sequence))]
            read.query_qualities = pysam.qualitystring_to_array(
                "I" * len(read.query_sequence)
            )
            read.set_tag("XY", "existing")
            read.set_tag("ML", array("B", [200]))
            outf.write(read)

        process_bam(
            str(input_bam),
            str(output_bam),
            DEFAULT_ADAPTER_5P,
            [("default", DEFAULT_ADAPTER_3P)],
            min_score_5p=20,  # Lower threshold for test data
            min_score_3p=20,
            match=2,
            mismatch=-1,
            gap_open=2,
            gap_extend=1,
        )

        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            read = next(bam)
            assert read.has_tag("XY")
            assert read.get_tag("XY") == "existing"
            assert read.has_tag("ML")
            assert read.has_tag("PT")

    def test_alignment_based_5p_detection(self, temp_dir):
        """Alignment-based detection should infer 5' adapter from ref position."""
        input_bam = temp_dir / "input.bam"
        output_bam = temp_dir / "output.bam"

        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref", "LN": 200}]}

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            # Read starting at ref position 10 (truncated, missing first 10bp of adapter)
            # This read would NOT be detected by sequence-based detection
            # but SHOULD be detected by alignment-based detection
            read = pysam.AlignedSegment()
            read.query_name = "truncated_read"
            # Sequence without the full 5' adapter (starts mid-adapter)
            read.query_sequence = "AAGAAGCCTGG" + "ACGTACGTACGTACGT" + DEFAULT_ADAPTER_3P
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 10  # Starts at position 10 (mid-adapter)
            read.cigartuples = [(0, len(read.query_sequence))]
            read.query_qualities = pysam.qualitystring_to_array(
                "I" * len(read.query_sequence)
            )
            outf.write(read)

        # Without alignment-based detection
        stats_no_infer = process_bam(
            str(input_bam),
            str(output_bam),
            DEFAULT_ADAPTER_5P,
            [("default", DEFAULT_ADAPTER_3P)],
            min_score_5p=30,  # High threshold - won't detect partial
            min_score_3p=20,
            match=2,
            mismatch=-1,
            gap_open=2,
            gap_extend=1,
            infer_5p_from_alignment=False,
        )

        # Should NOT detect 5' adapter without inference
        assert stats_no_infer.with_5p == 0

        # With alignment-based detection
        stats_with_infer = process_bam(
            str(input_bam),
            str(output_bam),
            DEFAULT_ADAPTER_5P,
            [("default", DEFAULT_ADAPTER_3P)],
            min_score_5p=30,
            min_score_3p=20,
            match=2,
            mismatch=-1,
            gap_open=2,
            gap_extend=1,
            infer_5p_from_alignment=True,
            max_ref_start_for_5p=20,
        )

        # SHOULD detect 5' adapter with inference (ref_start=10 < 20)
        assert stats_with_infer.with_5p == 1

        # Check the PT tag
        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            read = next(bam)
            assert read.has_tag("PT")
            pt_tag = read.get_tag("PT")
            assert "5p_adapter" in pt_tag
