"""
Unit tests for build_trna_reference.py

Tests FASTA reading/writing, adapter matching, reference validation,
and reference building with CCA handling.
"""

import pytest
from pathlib import Path

from build_trna_reference import (
    read_fasta,
    write_fasta,
    adapter_matches,
    validate_reference,
    build_reference,
    collapse_similar_sequences,
)

# Default adapters
ADAPTER_5P = "CCTAAGAGCAAGAAGAAGCCTGG"
ADAPTER_3P = "GGCTTCTTCTTGCTCTTCCAACCTTGCCTTAAAAAAAAAA"


class TestReadFasta:
    """Tests for read_fasta function."""

    def test_single_sequence(self, temp_dir):
        """Read single-line FASTA."""
        fasta = temp_dir / "test.fa"
        fasta.write_text(">seq1\nACGTACGT\n")

        sequences = list(read_fasta(str(fasta)))

        assert len(sequences) == 1
        assert sequences[0] == ("seq1", "ACGTACGT")

    def test_multiple_sequences(self, temp_dir):
        """Read multiple sequences."""
        fasta = temp_dir / "test.fa"
        fasta.write_text(">seq1\nACGT\n>seq2\nTGCA\n")

        sequences = list(read_fasta(str(fasta)))

        assert len(sequences) == 2
        assert sequences[0] == ("seq1", "ACGT")
        assert sequences[1] == ("seq2", "TGCA")

    def test_multiline_sequence(self, temp_dir):
        """Read sequence split across multiple lines."""
        fasta = temp_dir / "test.fa"
        fasta.write_text(">seq1\nACGT\nTGCA\nAAAA\n")

        sequences = list(read_fasta(str(fasta)))

        assert len(sequences) == 1
        assert sequences[0] == ("seq1", "ACGTTGCAAAAA")

    def test_uppercase_conversion(self, temp_dir):
        """Sequences should be converted to uppercase."""
        fasta = temp_dir / "test.fa"
        fasta.write_text(">seq1\nacgtACGT\n")

        sequences = list(read_fasta(str(fasta)))

        assert sequences[0] == ("seq1", "ACGTACGT")

    def test_rna_to_dna_conversion(self, temp_dir):
        """RNA (U) sequences should be normalized to DNA (T).

        GtRNAdb mature-tRNA FASTAs use the RNA alphabet; downstream BWA
        requires DNA.
        """
        fasta = temp_dir / "test.fa"
        fasta.write_text(">seq1\nACGUacgu\n")

        sequences = list(read_fasta(str(fasta)))

        assert sequences[0] == ("seq1", "ACGTACGT")

    def test_header_parsing(self, temp_dir):
        """Only first word of header should be used as name."""
        fasta = temp_dir / "test.fa"
        fasta.write_text(">seq1 description here\nACGT\n")

        sequences = list(read_fasta(str(fasta)))

        assert sequences[0][0] == "seq1"

    def test_empty_lines_ignored(self, temp_dir):
        """Empty lines should be skipped."""
        fasta = temp_dir / "test.fa"
        fasta.write_text(">seq1\n\nACGT\n\n>seq2\nTGCA\n")

        sequences = list(read_fasta(str(fasta)))

        assert len(sequences) == 2
        assert sequences[0] == ("seq1", "ACGT")


class TestWriteFasta:
    """Tests for write_fasta function."""

    def test_write_single_sequence(self, temp_dir):
        """Write single sequence."""
        output = temp_dir / "out.fa"
        sequences = [("seq1", "ACGTACGT")]

        write_fasta(sequences, str(output))

        content = output.read_text()
        assert ">seq1\n" in content
        assert "ACGTACGT\n" in content

    def test_write_multiple_sequences(self, temp_dir):
        """Write multiple sequences."""
        output = temp_dir / "out.fa"
        sequences = [("seq1", "ACGT"), ("seq2", "TGCA")]

        write_fasta(sequences, str(output))

        content = output.read_text()
        assert ">seq1\n" in content
        assert ">seq2\n" in content

    def test_roundtrip(self, temp_dir):
        """Read and write should preserve sequences."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        input_fa.write_text(">seq1\nACGTACGT\n>seq2\nTGCATGCA\n")

        sequences = list(read_fasta(str(input_fa)))
        write_fasta(sequences, str(output_fa))
        result = list(read_fasta(str(output_fa)))

        assert sequences == result


class TestAdapterMatches:
    """Tests for adapter_matches function."""

    def test_exact_match(self):
        """Exact match should return True."""
        assert adapter_matches("ACGT", "ACGT") is True

    def test_no_match(self):
        """Mismatch should return False."""
        assert adapter_matches("ACGT", "ACGA") is False

    def test_wildcard_n(self):
        """N should match any base."""
        assert adapter_matches("ACNT", "ACGT") is True
        assert adapter_matches("ACNT", "ACAT") is True
        assert adapter_matches("ACNT", "ACCT") is True

    def test_multiple_wildcards(self):
        """Multiple N wildcards should all match."""
        assert adapter_matches("NNNN", "ACGT") is True

    def test_length_mismatch(self):
        """Different lengths should return False."""
        assert adapter_matches("ACGT", "ACG") is False
        assert adapter_matches("ACG", "ACGT") is False


class TestValidateReference:
    """Tests for validate_reference function."""

    def test_valid_reference(self, temp_dir):
        """Valid adapted reference should pass validation."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        report = temp_dir / "report.txt"

        # Create valid adapted sequence
        trna = "GCGGCTATAGCTCAGTTGGTACCA"  # Ends with CCA
        adapted = ADAPTER_5P + "G" + trna + ADAPTER_3P
        input_fa.write_text(f">tRNA-Test\n{adapted}\n")

        result = validate_reference(
            str(input_fa),
            str(output_fa),
            str(report),
            ADAPTER_5P,
            [ADAPTER_3P],
        )

        assert result is True
        assert output_fa.exists()
        assert "VALIDATION PASSED" in report.read_text()

    def test_invalid_5p_adapter(self, temp_dir):
        """Wrong 5' adapter should fail validation."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        report = temp_dir / "report.txt"

        # Wrong 5' adapter
        bad_5p = "AAAAAAAAAAAAAAAAAAAAAA" + "G"
        adapted = bad_5p + "GCGGCTATAGCTCAGTTGGTACCA" + ADAPTER_3P
        input_fa.write_text(f">tRNA-Test\n{adapted}\n")

        with pytest.raises(SystemExit):
            validate_reference(
                str(input_fa),
                str(output_fa),
                str(report),
                ADAPTER_5P,
                [ADAPTER_3P],
            )

        assert "VALIDATION FAILED" in report.read_text()

    def test_invalid_3p_adapter(self, temp_dir):
        """Wrong 3' adapter should fail validation."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        report = temp_dir / "report.txt"

        # Wrong 3' adapter
        bad_3p = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        adapted = ADAPTER_5P + "G" + "GCGGCTATAGCTCAGTTGGTACCA" + bad_3p
        input_fa.write_text(f">tRNA-Test\n{adapted}\n")

        with pytest.raises(SystemExit):
            validate_reference(
                str(input_fa),
                str(output_fa),
                str(report),
                ADAPTER_5P,
                [ADAPTER_3P],
            )

        assert "VALIDATION FAILED" in report.read_text()

    def test_missing_cca(self, temp_dir):
        """Sequence without CCA should fail validation."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        report = temp_dir / "report.txt"

        # tRNA without CCA ending
        trna = "GCGGCTATAGCTCAGTTGGTATTT"  # Ends with TTT, not CCA
        adapted = ADAPTER_5P + "G" + trna + ADAPTER_3P
        input_fa.write_text(f">tRNA-Test\n{adapted}\n")

        with pytest.raises(SystemExit):
            validate_reference(
                str(input_fa),
                str(output_fa),
                str(report),
                ADAPTER_5P,
                [ADAPTER_3P],
            )

        assert "does not end with CCA" in report.read_text()

    def test_duplicate_names(self, temp_dir):
        """Duplicate sequence names should fail validation."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        report = temp_dir / "report.txt"

        trna = "GCGGCTATAGCTCAGTTGGTACCA"
        adapted = ADAPTER_5P + "G" + trna + ADAPTER_3P
        # Same name twice
        input_fa.write_text(f">tRNA-Test\n{adapted}\n>tRNA-Test\n{adapted}\n")

        with pytest.raises(SystemExit):
            validate_reference(
                str(input_fa),
                str(output_fa),
                str(report),
                ADAPTER_5P,
                [ADAPTER_3P],
            )

        assert "Duplicate sequence name" in report.read_text()

    def test_sequence_too_short(self, temp_dir):
        """Sequence shorter than minimum should fail."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        report = temp_dir / "report.txt"

        input_fa.write_text(">tRNA-Test\nACGT\n")

        with pytest.raises(SystemExit):
            validate_reference(
                str(input_fa),
                str(output_fa),
                str(report),
                ADAPTER_5P,
                [ADAPTER_3P],
            )

        assert "too short" in report.read_text()

    def test_adapter_3p_must_start_ggc(self, temp_dir):
        """3' adapter not starting with GGC should fail."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        report = temp_dir / "report.txt"

        trna = "GCGGCTATAGCTCAGTTGGTACCA"
        bad_3p = "AACTTCTTCTTGCTCTTCCAACCTTGCCTTAAAAAAAAAA"  # Starts with AAC
        adapted = ADAPTER_5P + "G" + trna + bad_3p
        input_fa.write_text(f">tRNA-Test\n{adapted}\n")

        with pytest.raises(SystemExit):
            validate_reference(
                str(input_fa),
                str(output_fa),
                str(report),
                ADAPTER_5P,
                [bad_3p],
            )

        assert "must start with GGC" in report.read_text()


class TestBuildReference:
    """Tests for build_reference function."""

    def test_build_with_cca(self, temp_dir):
        """Build reference from tRNA already having CCA."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        report = temp_dir / "report.txt"

        # Raw tRNA with CCA
        input_fa.write_text(">tRNA-Test\nGCGGCTATAGCTCAGTTGGTACCA\n")

        result = build_reference(
            str(input_fa),
            str(output_fa),
            str(report),
            ADAPTER_5P,
            ADAPTER_3P,
        )

        assert result is True
        assert output_fa.exists()

        # Check output structure
        sequences = list(read_fasta(str(output_fa)))
        assert len(sequences) == 1
        name, seq = sequences[0]

        # Verify adapters present
        assert seq.startswith(ADAPTER_5P)
        assert seq.endswith(ADAPTER_3P)

        # Verify CCAGGC junction
        assert "CCAGGC" in seq

        report_text = report.read_text()
        assert "BUILD SUCCESSFUL" in report_text
        assert "Sequences with CCA: 1" in report_text

    def test_build_adds_cca(self, temp_dir):
        """Build should add CCA when missing."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        report = temp_dir / "report.txt"

        # Raw tRNA without CCA
        input_fa.write_text(">tRNA-Test\nGCGGCTATAGCTCAGTTGGTA\n")

        result = build_reference(
            str(input_fa),
            str(output_fa),
            str(report),
            ADAPTER_5P,
            ADAPTER_3P,
        )

        assert result is True

        # Check output has CCA added
        sequences = list(read_fasta(str(output_fa)))
        name, seq = sequences[0]

        # CCAGGC junction should exist
        assert "CCAGGC" in seq

        report_text = report.read_text()
        assert "Sequences with CCA added: 1" in report_text

    def test_build_multiple_sequences(self, temp_dir):
        """Build should handle multiple sequences."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        report = temp_dir / "report.txt"

        input_fa.write_text(
            ">tRNA-Ala\nGCGGCTATAGCTCAGTTGGTACCA\n"
            ">tRNA-Gly\nGCATGCATGCATCCA\n"
        )

        result = build_reference(
            str(input_fa),
            str(output_fa),
            str(report),
            ADAPTER_5P,
            ADAPTER_3P,
        )

        assert result is True

        sequences = list(read_fasta(str(output_fa)))
        assert len(sequences) == 2

        report_text = report.read_text()
        assert "Sequences built: 2" in report_text

    def test_build_rejects_duplicates(self, temp_dir):
        """Build should reject duplicate sequence names."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        report = temp_dir / "report.txt"

        input_fa.write_text(
            ">tRNA-Test\nGCGGCTATAGCTCAGTTGGTACCA\n"
            ">tRNA-Test\nGCATGCATGCATCCA\n"
        )

        with pytest.raises(SystemExit):
            build_reference(
                str(input_fa),
                str(output_fa),
                str(report),
                ADAPTER_5P,
                ADAPTER_3P,
            )

        assert "Duplicate sequence name" in report.read_text()

    def test_build_rejects_bad_3p_adapter(self, temp_dir):
        """Build should reject 3' adapter not starting with GGC."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        report = temp_dir / "report.txt"

        input_fa.write_text(">tRNA-Test\nGCGGCTATAGCTCAGTTGGTACCA\n")

        bad_3p = "AACTTCTTCTTGCTCTTCCAACCTTGCCTTAAAAAAAAAA"  # Starts with AAC

        with pytest.raises(SystemExit):
            build_reference(
                str(input_fa),
                str(output_fa),
                str(report),
                ADAPTER_5P,
                bad_3p,
            )

        assert "must start with GGC" in report.read_text()


class TestCollapseSimilarSequences:
    """Merging near-identical reference sequences by Hamming distance."""

    def test_zero_keeps_everything(self):
        keys = [("a", "ACGTACGT"), ("b", "ACGTACGA")]
        kept, collapsed = collapse_similar_sequences(keys, 0)
        assert kept == ["a", "b"]
        assert collapsed == {}

    def test_merges_within_distance(self):
        keys = [("a", "ACGTACGT"), ("b", "ACGTACGA")]
        kept, collapsed = collapse_similar_sequences(keys, 1)
        assert kept == ["a"]
        assert collapsed == {"a": ["b"]}

    def test_does_not_merge_beyond_distance(self):
        keys = [("a", "ACGTACGT"), ("b", "ACGTAAAA")]
        kept, collapsed = collapse_similar_sequences(keys, 2)
        assert kept == ["a", "b"]
        assert collapsed == {}

    def test_never_merges_different_lengths(self):
        """Hamming distance is undefined across lengths, so an indel never merges."""
        keys = [("a", "ACGTACGT"), ("b", "ACGTACGTA")]
        kept, collapsed = collapse_similar_sequences(keys, 99)
        assert kept == ["a", "b"]
        assert collapsed == {}

    def test_bounded_radius_not_single_linkage(self):
        """a-b and b-c are each 1 apart but a-c is 2, so c must not join a."""
        keys = [("a", "AAAA"), ("b", "AAAT"), ("c", "AATT")]
        kept, collapsed = collapse_similar_sequences(keys, 1)
        assert kept == ["a", "c"]
        assert collapsed == {"a": ["b"]}

    def test_every_name_accounted_for(self):
        keys = [("a", "AAAA"), ("b", "AAAT"), ("c", "CCCC"), ("d", "CCCG")]
        kept, collapsed = collapse_similar_sequences(keys, 1)
        seen = set(kept) | {m for v in collapsed.values() for m in v}
        assert seen == {"a", "b", "c", "d"}

    def test_leader_is_never_also_merged(self):
        keys = [("a", "AAAA"), ("b", "AAAT"), ("c", "AAAA")]
        kept, collapsed = collapse_similar_sequences(keys, 1)
        merged = {m for v in collapsed.values() for m in v}
        assert not (set(kept) & merged)

    def test_input_order_preserved(self):
        keys = [("z", "AAAA"), ("y", "CCCC"), ("x", "GGGG")]
        kept, _ = collapse_similar_sequences(keys, 1)
        assert kept == ["z", "y", "x"]


class TestBuildReferenceCollapse:
    """max_mismatch wiring through build_reference."""

    def test_collapses_near_identical_in_build(self, tmp_path):
        # Two tRNAs differing by one base, plus one clearly distinct
        fasta = tmp_path / "raw.fa"
        fasta.write_text(
            ">t1\nGGGGGGGGGGGGGGGGGGGGCCA\n"
            ">t2\nGGGGGGGGGGGGGGGGGGGACCA\n"
            ">t3\nAAAAAAAAAAAAAAAAAAAACCA\n"
        )
        out = tmp_path / "adapted.fa"
        report = tmp_path / "report.txt"

        build_reference(
            str(fasta), str(out), str(report), ADAPTER_5P, ADAPTER_3P, max_mismatch=1
        )

        names = [n for n, _ in read_fasta(str(out))]
        assert names == ["t1", "t3"]
        assert "t1 <- t2" in report.read_text()

    def test_default_does_not_collapse_near_identical(self, tmp_path):
        fasta = tmp_path / "raw.fa"
        fasta.write_text(
            ">t1\nGGGGGGGGGGGGGGGGGGGGCCA\n>t2\nGGGGGGGGGGGGGGGGGGGACCA\n"
        )
        out = tmp_path / "adapted.fa"
        report = tmp_path / "report.txt"

        build_reference(str(fasta), str(out), str(report), ADAPTER_5P, ADAPTER_3P)

        names = [n for n, _ in read_fasta(str(out))]
        assert names == ["t1", "t2"]

    def test_grouping_uses_input_lengths_not_post_cca(self, tmp_path):
        """Clustering keys are the input sequences, before CCA is appended.

        t1 already ends in CCA (23 nt); t2 lacks it (20 nt). As provided they
        differ in length, so Hamming distance is undefined and they must not
        merge. After CCA is appended both are 23 nt and differ by a single
        base, so clustering on the post-CCA sequence would wrongly merge them.
        """
        fasta = tmp_path / "raw.fa"
        fasta.write_text(
            ">t1\nGGGGGGGGGGGGGGGGGGGGCCA\n>t2\nGGGGGGGGGGGGGGGGGGGA\n"
        )
        out = tmp_path / "adapted.fa"
        report = tmp_path / "report.txt"

        build_reference(
            str(fasta), str(out), str(report), ADAPTER_5P, ADAPTER_3P, max_mismatch=1
        )

        names = [n for n, _ in read_fasta(str(out))]
        assert names == ["t1", "t2"]
