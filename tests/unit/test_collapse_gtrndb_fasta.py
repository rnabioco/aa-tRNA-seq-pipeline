"""
Unit tests for collapse_gtrndb_fasta.py

Tests header parsing, CCA stripping, sequence collapsing, and end-to-end
FASTA processing for the GtRNAdb redundancy collapsing tool.
"""

import pytest

from collapse_gtrndb_fasta import (
    parse_trna_name,
    strip_trailing_cca,
    collapse_sequences,
    read_fasta,
    write_fasta,
    write_mapping,
)


class TestParseTrnaName:
    """Tests for parse_trna_name function."""

    def test_with_species_prefix(self):
        """Standard GtRNAdb header with species prefix."""
        result = parse_trna_name(
            "Escherichia_coli_str_K-12_substr_MG1655_tRNA-Ala-GGC-1-1"
        )
        assert result is not None
        assert result["trna_type"] == "tRNA"
        assert result["amino_acid"] == "Ala"
        assert result["anticodon"] == "GGC"
        assert result["family_num"] == "1"
        assert result["copy_num"] == "1"
        assert result["short_name"] == "tRNA-Ala-GGC-1-1"
        assert "K-12" in result["prefix"]

    def test_no_prefix(self):
        """Header with no species prefix."""
        result = parse_trna_name("tRNA-Gly-GCC-2-3")
        assert result is not None
        assert result["trna_type"] == "tRNA"
        assert result["amino_acid"] == "Gly"
        assert result["anticodon"] == "GCC"
        assert result["family_num"] == "2"
        assert result["copy_num"] == "3"
        assert result["short_name"] == "tRNA-Gly-GCC-2-3"
        assert result["prefix"] == ""

    def test_pretrna(self):
        """pretRNA type should be parsed."""
        result = parse_trna_name("Ecoli_pretRNA-Arg-CCG-1-1")
        assert result is not None
        assert result["trna_type"] == "pretRNA"
        assert result["amino_acid"] == "Arg"
        assert result["short_name"] == "pretRNA-Arg-CCG-1-1"

    def test_ile2(self):
        """Ile2 amino acid variant."""
        result = parse_trna_name("Ecoli_tRNA-Ile2-CAU-1-1")
        assert result is not None
        assert result["amino_acid"] == "Ile2"

    def test_fmet(self):
        """fMet initiator tRNA."""
        result = parse_trna_name("Ecoli_tRNA-fMet-CAU-1-1")
        assert result is not None
        assert result["amino_acid"] == "fMet"

    def test_imet(self):
        """iMet initiator tRNA."""
        result = parse_trna_name("Human_tRNA-iMet-CAU-1-1")
        assert result is not None
        assert result["amino_acid"] == "iMet"

    def test_sec(self):
        """SeC (selenocysteine) tRNA."""
        result = parse_trna_name("Human_tRNA-SeC-TCA-1-1")
        assert result is not None
        assert result["amino_acid"] == "SeC"

    def test_und(self):
        """Und (undetermined) tRNA."""
        result = parse_trna_name("Ecoli_tRNA-Und-NNN-1-1")
        assert result is not None
        assert result["amino_acid"] == "Und"
        assert result["anticodon"] == "NNN"

    def test_sup(self):
        """Sup (suppressor) tRNA."""
        result = parse_trna_name("Ecoli_tRNA-Sup-CTA-1-1")
        assert result is not None
        assert result["amino_acid"] == "Sup"

    def test_hyphenated_species(self):
        """Species name with hyphens (e.g. K-12)."""
        result = parse_trna_name(
            "Escherichia_coli_str_K-12_substr_MG1655_tRNA-Lys-TTT-3-2"
        )
        assert result is not None
        assert result["amino_acid"] == "Lys"
        assert result["anticodon"] == "TTT"
        assert result["family_num"] == "3"
        assert result["copy_num"] == "2"

    def test_unparseable(self):
        """Non-GtRNAdb headers should return None."""
        assert parse_trna_name("random_sequence_1") is None
        assert parse_trna_name("chr1:1000-2000") is None
        assert parse_trna_name("") is None

    def test_multi_digit_family_copy(self):
        """Multi-digit family and copy numbers."""
        result = parse_trna_name("Human_tRNA-Ala-AGC-12-34")
        assert result is not None
        assert result["family_num"] == "12"
        assert result["copy_num"] == "34"
        assert result["short_name"] == "tRNA-Ala-AGC-12-34"


class TestStripTrailingCca:
    """Tests for strip_trailing_cca function."""

    def test_with_cca(self):
        """Sequence ending with CCA should be stripped."""
        assert strip_trailing_cca("GCGGCTATAGCTCAGTTGGTACCA") == "GCGGCTATAGCTCAGTTGGTA"

    def test_without_cca(self):
        """Sequence not ending with CCA should be unchanged."""
        assert strip_trailing_cca("GCGGCTATAGCTCAGTTGGTA") == "GCGGCTATAGCTCAGTTGGTA"

    def test_partial_cc(self):
        """Partial CCA (CC only) should not be stripped."""
        assert strip_trailing_cca("GCGGCTATAGCTCAGTTGGTACC") == "GCGGCTATAGCTCAGTTGGTACC"

    def test_short_sequence(self):
        """Short sequences should be handled."""
        assert strip_trailing_cca("CCA") == ""
        assert strip_trailing_cca("CA") == "CA"
        assert strip_trailing_cca("A") == "A"

    def test_lowercase_cca(self):
        """Lowercase CCA should also be stripped (via upper())."""
        assert strip_trailing_cca("ACGTcca") == "ACGT"


class TestCollapseSequences:
    """Tests for collapse_sequences function."""

    def test_identical_copies_collapse(self):
        """Identical gene copies within an isodecoder should collapse."""
        records = [
            ("Ecoli_tRNA-Ala-GGC-1-1", "GCGGCTATAGCTCAGTTGGTACCA"),
            ("Ecoli_tRNA-Ala-GGC-1-2", "GCGGCTATAGCTCAGTTGGTACCA"),
            ("Ecoli_tRNA-Ala-GGC-1-3", "GCGGCTATAGCTCAGTTGGTACCA"),
        ]
        output_seqs, mapping = collapse_sequences(records)

        assert len(output_seqs) == 1
        assert output_seqs[0][0] == "tRNA-Ala-GGC-1-1"
        assert output_seqs[0][1] == "GCGGCTATAGCTCAGTTGGTACCA"

        assert len(mapping) == 3
        assert mapping[0]["is_representative"] is True
        assert mapping[1]["is_representative"] is False
        assert mapping[2]["is_representative"] is False

        # All map to same collapsed name
        for row in mapping:
            assert row["collapsed_name"] == "tRNA-Ala-GGC-1-1"

    def test_different_families_kept_separate(self):
        """Different family numbers with different sequences should be kept."""
        records = [
            ("Ecoli_tRNA-Ala-GGC-1-1", "GCGGCTATAGCTCAGTTGGTACCA"),
            ("Ecoli_tRNA-Ala-GGC-2-1", "AAAGCTATAGCTCAGTTGGTACCA"),
        ]
        output_seqs, mapping = collapse_sequences(records)

        assert len(output_seqs) == 2
        names = {s[0] for s in output_seqs}
        assert "tRNA-Ala-GGC-1-1" in names
        assert "tRNA-Ala-GGC-2-1" in names

    def test_cca_stripping_for_comparison(self):
        """Sequences differing only in CCA should collapse."""
        records = [
            ("Ecoli_tRNA-Ala-GGC-1-1", "GCGGCTATAGCTCAGTTGGTACCA"),
            ("Ecoli_tRNA-Ala-GGC-1-2", "GCGGCTATAGCTCAGTTGGTA"),
        ]
        output_seqs, mapping = collapse_sequences(records)

        assert len(output_seqs) == 1
        # Representative keeps its original sequence (with CCA)
        assert output_seqs[0][0] == "tRNA-Ala-GGC-1-1"
        assert output_seqs[0][1] == "GCGGCTATAGCTCAGTTGGTACCA"

    def test_mixed_identical_and_different(self):
        """Mix of identical and different sequences within isodecoder."""
        records = [
            ("Ecoli_tRNA-Lys-TTT-1-1", "AAAAACCA"),
            ("Ecoli_tRNA-Lys-TTT-1-2", "AAAAACCA"),
            ("Ecoli_tRNA-Lys-TTT-2-1", "BBBBBCCA"),
        ]
        output_seqs, mapping = collapse_sequences(records)

        assert len(output_seqs) == 2
        names = {s[0] for s in output_seqs}
        assert "tRNA-Lys-TTT-1-1" in names
        assert "tRNA-Lys-TTT-2-1" in names

    def test_single_copy_passthrough(self):
        """Single-copy families should pass through unchanged."""
        records = [
            ("Ecoli_tRNA-Trp-CCA-1-1", "GCGGCTATAGCTCAGTTGGTACCA"),
        ]
        output_seqs, mapping = collapse_sequences(records)

        assert len(output_seqs) == 1
        assert output_seqs[0][0] == "tRNA-Trp-CCA-1-1"
        assert len(mapping) == 1
        assert mapping[0]["is_representative"] is True

    def test_different_isodecoders_independent(self):
        """Different isodecoders should be handled independently."""
        records = [
            ("Ecoli_tRNA-Ala-GGC-1-1", "AAAAACCA"),
            ("Ecoli_tRNA-Gly-GCC-1-1", "AAAAACCA"),
        ]
        output_seqs, mapping = collapse_sequences(records)

        # Same sequence but different isodecoders => both kept
        assert len(output_seqs) == 2

    def test_unparseable_error_by_default(self):
        """Unparseable headers should raise ValueError by default."""
        records = [("random_seq", "ACGT")]
        with pytest.raises(ValueError, match="Could not parse"):
            collapse_sequences(records)

    def test_unparseable_kept_with_flag(self):
        """Unparseable headers should be kept with --keep-unparsed."""
        records = [
            ("random_seq", "ACGT"),
            ("Ecoli_tRNA-Ala-GGC-1-1", "GCGGCTATAGCTCAGTTGGTACCA"),
        ]
        output_seqs, mapping = collapse_sequences(records, keep_unparsed=True)

        assert len(output_seqs) == 2
        names = {s[0] for s in output_seqs}
        assert "random_seq" in names
        assert "tRNA-Ala-GGC-1-1" in names

        unparsed_rows = [r for r in mapping if r["isodecoder"] == "unparsed"]
        assert len(unparsed_rows) == 1

    def test_duplicate_input_names_error(self):
        """Duplicate input names should raise ValueError."""
        records = [
            ("Ecoli_tRNA-Ala-GGC-1-1", "AAAAACCA"),
            ("Ecoli_tRNA-Ala-GGC-1-1", "BBBBBCCA"),
        ]
        with pytest.raises(ValueError, match="Duplicate input name"):
            collapse_sequences(records)

    def test_pretrna_and_trna_not_collapsed(self):
        """pretRNA and tRNA of same anticodon should not collapse together."""
        records = [
            ("Ecoli_tRNA-Ala-GGC-1-1", "AAAAACCA"),
            ("Ecoli_pretRNA-Ala-GGC-1-1", "AAAAACCA"),
        ]
        output_seqs, mapping = collapse_sequences(records)

        assert len(output_seqs) == 2
        names = {s[0] for s in output_seqs}
        assert "tRNA-Ala-GGC-1-1" in names
        assert "pretRNA-Ala-GGC-1-1" in names

    def test_mapping_fields(self):
        """Mapping rows should have all expected fields."""
        records = [("Ecoli_tRNA-Ala-GGC-1-1", "GCGGCTATAGCTCAGTTGGTACCA")]
        _, mapping = collapse_sequences(records)

        row = mapping[0]
        assert row["collapsed_name"] == "tRNA-Ala-GGC-1-1"
        assert row["original_name"] == "Ecoli_tRNA-Ala-GGC-1-1"
        assert row["isodecoder"] == "tRNA-Ala-GGC"
        assert row["amino_acid"] == "Ala"
        assert row["anticodon"] == "GGC"
        assert row["family_num"] == "1"
        assert row["copy_num"] == "1"
        assert row["is_representative"] is True
        assert row["sequence_length"] == 24


class TestEndToEnd:
    """End-to-end tests with FASTA files."""

    def test_full_pipeline(self, temp_dir):
        """Full FASTA in -> collapsed FASTA + mapping out."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        mapping_tsv = temp_dir / "mapping.tsv"

        # Write input FASTA with redundant copies
        input_fa.write_text(
            ">Ecoli_tRNA-Ala-GGC-1-1\n"
            "GCGGCTATAGCTCAGTTGGTACCA\n"
            ">Ecoli_tRNA-Ala-GGC-1-2\n"
            "GCGGCTATAGCTCAGTTGGTACCA\n"
            ">Ecoli_tRNA-Ala-GGC-2-1\n"
            "AAAGCTATAGCTCAGTTGGTACCA\n"
            ">Ecoli_tRNA-Gly-GCC-1-1\n"
            "TTTTTCCA\n"
        )

        # Read, collapse, write
        records = list(read_fasta(str(input_fa)))
        output_seqs, mapping_rows = collapse_sequences(records)
        write_fasta(output_seqs, str(output_fa))
        write_mapping(mapping_rows, str(mapping_tsv))

        # Verify output FASTA
        result_seqs = list(read_fasta(str(output_fa)))
        assert len(result_seqs) == 3
        names = [s[0] for s in result_seqs]
        assert names == ["tRNA-Ala-GGC-1-1", "tRNA-Ala-GGC-2-1", "tRNA-Gly-GCC-1-1"]

        # Verify mapping TSV
        mapping_text = mapping_tsv.read_text()
        lines = mapping_text.strip().split("\n")
        assert len(lines) == 5  # header + 4 data rows
        assert "collapsed_name" in lines[0]

        # Check representative flags
        assert "True" in lines[1]  # Ala-GGC-1-1 is representative
        assert "False" in lines[2]  # Ala-GGC-1-2 is not

    def test_roundtrip_single_copy(self, temp_dir):
        """Single-copy families should round-trip with stripped prefix."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        mapping_tsv = temp_dir / "mapping.tsv"

        input_fa.write_text(
            ">Ecoli_tRNA-Trp-CCA-1-1\n"
            "GCGGCTATAGCTCAGTTGGTACCA\n"
        )

        records = list(read_fasta(str(input_fa)))
        output_seqs, mapping_rows = collapse_sequences(records)
        write_fasta(output_seqs, str(output_fa))
        write_mapping(mapping_rows, str(mapping_tsv))

        result_seqs = list(read_fasta(str(output_fa)))
        assert len(result_seqs) == 1
        assert result_seqs[0][0] == "tRNA-Trp-CCA-1-1"
        assert result_seqs[0][1] == "GCGGCTATAGCTCAGTTGGTACCA"

    def test_empty_fasta_detected(self, temp_dir):
        """Empty input should produce empty records list."""
        input_fa = temp_dir / "empty.fa"
        input_fa.write_text("")

        records = list(read_fasta(str(input_fa)))
        assert len(records) == 0

    def test_keep_unparsed_end_to_end(self, temp_dir):
        """Unparsed sequences should pass through in output."""
        input_fa = temp_dir / "input.fa"
        output_fa = temp_dir / "output.fa"
        mapping_tsv = temp_dir / "mapping.tsv"

        input_fa.write_text(
            ">adapter_sequence\n"
            "CCTAAGAGCAAGAAGAAGCCTGG\n"
            ">Ecoli_tRNA-Ala-GGC-1-1\n"
            "GCGGCTATAGCTCAGTTGGTACCA\n"
        )

        records = list(read_fasta(str(input_fa)))
        output_seqs, mapping_rows = collapse_sequences(records, keep_unparsed=True)
        write_fasta(output_seqs, str(output_fa))
        write_mapping(mapping_rows, str(mapping_tsv))

        result_seqs = list(read_fasta(str(output_fa)))
        assert len(result_seqs) == 2
        names = {s[0] for s in result_seqs}
        assert "adapter_sequence" in names
        assert "tRNA-Ala-GGC-1-1" in names
