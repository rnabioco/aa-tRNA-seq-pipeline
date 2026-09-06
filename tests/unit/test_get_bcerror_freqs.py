"""Unit tests for get_bcerror_freqs.py."""

import pysam
import pytest

from conftest import create_bam_with_reads
from get_bcerror_freqs import calculate_error_frequencies


@pytest.fixture
def simple_fasta(temp_dir):
    """Create a simple indexed FASTA."""
    fa_path = temp_dir / "ref.fa"
    fa_path.write_text(">ref1\nACGTACGT\n")
    pysam.faidx(str(fa_path))
    return fa_path


@pytest.fixture
def two_reference_fasta(temp_dir):
    """Two references that differ at every coordinate."""
    fa_path = temp_dir / "two.fa"
    fa_path.write_text(">ref1\nACGTACGT\n>ref2\nTTTTTTTT\n")
    pysam.faidx(str(fa_path))
    return fa_path


@pytest.fixture
def adapted_fasta(temp_dir):
    """Create FASTA with adapter regions for trimming tests."""
    fa_path = temp_dir / "ref.fa"
    # 4bp 5' adapter + 8bp tRNA + 4bp 3' adapter = 16bp total
    fa_path.write_text(">ref1\nAAAAACGTACGTTTTT\n")
    pysam.faidx(str(fa_path))
    return fa_path


class TestCalculateErrorFrequencies:
    def test_perfect_alignment(self, temp_dir, simple_fasta):
        """Perfect match should have all freq = 0."""
        bam_path = temp_dir / "test.bam"
        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref1", "LN": 8}]}
        reads = [
            {"name": "r1", "seq": "ACGTACGT", "flag": 0, "ref_id": 0, "ref_start": 0, "mapq": 60,
             "cigar": [(0, 8)]},
        ]
        create_bam_with_reads(bam_path, reads, header=header)
        df = calculate_error_frequencies(str(bam_path), str(simple_fasta))
        assert len(df) == 8
        assert all(df["MismatchFreq"] == 0)
        assert all(df["InsertionFreq"] == 0)
        assert all(df["DeletionFreq"] == 0)

    def test_single_mismatch(self, temp_dir, simple_fasta):
        """One mismatch at position 0: ref=A, read=T."""
        bam_path = temp_dir / "test.bam"
        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref1", "LN": 8}]}
        reads = [
            {"name": "r1", "seq": "TCGTACGT", "flag": 0, "ref_id": 0, "ref_start": 0, "mapq": 60,
             "cigar": [(0, 8)]},
        ]
        create_bam_with_reads(bam_path, reads, header=header)
        df = calculate_error_frequencies(str(bam_path), str(simple_fasta))
        row0 = df[df["Position"] == 1].iloc[0]
        assert row0["MismatchFreq"] == 1.0
        # Other positions should be 0
        assert all(df[df["Position"] != 1]["MismatchFreq"] == 0)

    def test_insertion(self, temp_dir, simple_fasta):
        """Read with insertion at position 2."""
        bam_path = temp_dir / "test.bam"
        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref1", "LN": 8}]}
        # CIGAR: 2M1I6M = match 2, insert 1, match 6 → aligns to full 8bp ref
        reads = [
            {"name": "r1", "seq": "ACNGTACGT", "flag": 0, "ref_id": 0, "ref_start": 0, "mapq": 60,
             "cigar": [(0, 2), (1, 1), (0, 6)]},
        ]
        create_bam_with_reads(bam_path, reads, header=header)
        df = calculate_error_frequencies(str(bam_path), str(simple_fasta))
        # Insertion recorded at ref_pos=2 (Position 3 in 1-indexed)
        row = df[df["Position"] == 3].iloc[0]
        assert row["InsertionFreq"] == 1.0

    def test_deletion(self, temp_dir, simple_fasta):
        """Read with deletion at position 3."""
        bam_path = temp_dir / "test.bam"
        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref1", "LN": 8}]}
        # CIGAR: 3M1D4M → matches 3bp, deletes 1bp (pos 3), matches 4bp
        reads = [
            {"name": "r1", "seq": "ACGACGT", "flag": 0, "ref_id": 0, "ref_start": 0, "mapq": 60,
             "cigar": [(0, 3), (2, 1), (0, 4)]},
        ]
        create_bam_with_reads(bam_path, reads, header=header)
        df = calculate_error_frequencies(str(bam_path), str(simple_fasta))
        row = df[df["Position"] == 4].iloc[0]  # 1-indexed position 4 = 0-indexed 3
        assert row["DeletionFreq"] > 0

    def test_coverage_counting(self, temp_dir, simple_fasta):
        """Multiple reads should accumulate coverage."""
        bam_path = temp_dir / "test.bam"
        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref1", "LN": 8}]}
        reads = [
            {"name": f"r{i}", "seq": "ACGTACGT", "flag": 0, "ref_id": 0, "ref_start": 0, "mapq": 60,
             "cigar": [(0, 8)]}
            for i in range(5)
        ]
        create_bam_with_reads(bam_path, reads, header=header)
        df = calculate_error_frequencies(str(bam_path), str(simple_fasta))
        assert all(df["Spanning_Reads"] == 5)

    def test_reverse_strand_excluded(self, temp_dir, simple_fasta):
        """Reverse strand reads should be excluded."""
        bam_path = temp_dir / "test.bam"
        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref1", "LN": 8}]}
        reads = [
            {"name": "r1", "seq": "ACGTACGT", "flag": 0, "ref_id": 0, "ref_start": 0, "mapq": 60, "cigar": [(0, 8)]},
            {"name": "r2", "seq": "ACGTACGT", "flag": 16, "ref_id": 0, "ref_start": 0, "mapq": 60, "cigar": [(0, 8)]},
        ]
        create_bam_with_reads(bam_path, reads, header=header)
        df = calculate_error_frequencies(str(bam_path), str(simple_fasta))
        assert all(df["Spanning_Reads"] == 1)

    def test_adapter_trimming_with_offsets(self, temp_dir, adapted_fasta):
        """With offsets, adapter positions should be excluded and coords shifted."""
        bam_path = temp_dir / "test.bam"
        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref1", "LN": 16}]}
        reads = [
            {"name": "r1", "seq": "AAAACGTACGTTTTT", "flag": 0, "ref_id": 0, "ref_start": 0, "mapq": 60,
             "cigar": [(0, 15)]},
        ]
        create_bam_with_reads(bam_path, reads, header=header)
        df = calculate_error_frequencies(str(bam_path), str(adapted_fasta), trim_5p=4, trim_3p=4)
        # Should only have 8 positions (16 - 4 - 4) renamed to 1..8
        assert len(df) == 8
        assert df["Position"].min() == 1
        assert df["Position"].max() == 8

    def test_each_reference_is_scored_against_its_own_sequence(
        self, temp_dir, two_reference_fasta
    ):
        """
        Both reads match their own reference perfectly, so neither has a
        mismatch.

        The reference is read once per contig and then indexed by position,
        rather than fetched a base at a time inside the CIGAR loop. Hoisting
        that out of the loop is only safe if the sequence is re-read for every
        contig; getting it wrong scores all of them against whichever sequence
        was loaded first, and that is invisible until two references differ at
        the same coordinate.
        """
        bam_path = temp_dir / "two.bam"
        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "ref1", "LN": 8}, {"SN": "ref2", "LN": 8}],
        }
        reads = [
            {"name": "r1", "seq": "ACGTACGT", "flag": 0, "ref_id": 0},
            {"name": "r2", "seq": "TTTTTTTT", "flag": 0, "ref_id": 1},
        ]
        create_bam_with_reads(bam_path, reads, header=header)

        df = calculate_error_frequencies(str(bam_path), str(two_reference_fasta))

        for ref in ("ref1", "ref2"):
            rows = df[df["Reference"] == ref]
            assert rows["Bases_Mapped"].sum() == 8, ref
            assert rows["MismatchFreq"].sum() == 0, ref

    def test_nucleotide_freqs_sum(self, temp_dir, simple_fasta):
        """Nucleotide frequencies should sum to approximately 1.0 where bases are mapped."""
        bam_path = temp_dir / "test.bam"
        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "ref1", "LN": 8}]}
        reads = [
            {"name": "r1", "seq": "ACGTACGT", "flag": 0, "ref_id": 0, "ref_start": 0, "mapq": 60, "cigar": [(0, 8)]},
        ]
        create_bam_with_reads(bam_path, reads, header=header)
        df = calculate_error_frequencies(str(bam_path), str(simple_fasta))
        freq_sum = df["A_Freq"] + df["T_Freq"] + df["G_Freq"] + df["C_Freq"] + df["N_freq"]
        assert all(abs(freq_sum - 1.0) < 0.01)
