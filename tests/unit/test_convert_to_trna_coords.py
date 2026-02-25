"""Unit tests for convert_to_trna_coords.py."""

import io
import pytest

from convert_to_trna_coords import (
    process_bedgraph,
    process_bedmethyl,
    process_modkit_tsv,
)


REF_LENGTHS = {"tRNA-Ala": 100, "tRNA-Gly": 80}
OFFSET_5P = 24
OFFSET_3P = 40


class TestProcessBedgraph:
    def test_coordinate_shift(self):
        """Positions should be shifted by offset_5p."""
        infile = io.StringIO("tRNA-Ala\t30\t31\t50\n")
        outfile = io.BytesIO()
        process_bedgraph(infile, outfile, OFFSET_5P, OFFSET_3P, REF_LENGTHS)
        result = outfile.getvalue().decode()
        fields = result.strip().split("\t")
        assert fields[0] == "tRNA-Ala"
        assert fields[1] == "6"  # 30 - 24
        assert fields[2] == "7"  # 31 - 24
        assert fields[3] == "50"

    def test_adapter_filtering_5p(self):
        """Positions within 5' adapter should be filtered."""
        infile = io.StringIO("tRNA-Ala\t10\t11\t50\n")
        outfile = io.BytesIO()
        process_bedgraph(infile, outfile, OFFSET_5P, OFFSET_3P, REF_LENGTHS)
        assert outfile.getvalue() == b""

    def test_adapter_filtering_3p(self):
        """Positions within 3' adapter should be filtered."""
        # ref_len=100, offset_3p=40, max_pos=60, so end>60 is filtered
        infile = io.StringIO("tRNA-Ala\t55\t61\t50\n")
        outfile = io.BytesIO()
        process_bedgraph(infile, outfile, OFFSET_5P, OFFSET_3P, REF_LENGTHS)
        assert outfile.getvalue() == b""

    def test_missing_chrom_skipped(self):
        """Unknown chromosome should be skipped."""
        infile = io.StringIO("unknown_chrom\t30\t31\t50\n")
        outfile = io.BytesIO()
        process_bedgraph(infile, outfile, OFFSET_5P, OFFSET_3P, REF_LENGTHS)
        assert outfile.getvalue() == b""


class TestProcessBedmethyl:
    def test_thick_start_end_shifted(self):
        """thickStart and thickEnd (cols 6,7) should also be shifted."""
        # bedMethyl: chrom start end name score strand thickStart thickEnd ...
        line = "tRNA-Ala\t30\t31\tmod\t100\t+\t30\t31\t255,0,0\t1\n"
        infile = io.StringIO(line)
        outfile = io.BytesIO()
        process_bedmethyl(infile, outfile, OFFSET_5P, OFFSET_3P, REF_LENGTHS)
        result = outfile.getvalue().decode()
        fields = result.strip().split("\t")
        assert fields[1] == "6"   # start shifted
        assert fields[2] == "7"   # end shifted
        assert fields[6] == "6"   # thickStart shifted
        assert fields[7] == "7"   # thickEnd shifted

    def test_adapter_filtering(self):
        """Adapter positions should be filtered for bedMethyl too."""
        line = "tRNA-Ala\t10\t11\tmod\t100\t+\t10\t11\t255,0,0\t1\n"
        infile = io.StringIO(line)
        outfile = io.BytesIO()
        process_bedmethyl(infile, outfile, OFFSET_5P, OFFSET_3P, REF_LENGTHS)
        assert outfile.getvalue() == b""


class TestProcessModkitTsv:
    def test_header_preserved(self):
        """Header line should be passed through unchanged."""
        header = "read_id\tfwd_pos\tref_position\tchrom\tmod_strand\n"
        data = "r1\t5\t30\ttRNA-Ala\t+\n"
        infile = io.StringIO(header + data)
        outfile = io.BytesIO()
        process_modkit_tsv(infile, outfile, OFFSET_5P, OFFSET_3P, REF_LENGTHS)
        lines = outfile.getvalue().decode().strip().split("\n")
        assert lines[0] == header.rstrip("\n")

    def test_ref_position_to_1_indexed(self):
        """ref_position (col 2) should be converted to 1-indexed tRNA coords."""
        header = "read_id\tfwd_pos\tref_position\tchrom\tmod_strand\n"
        data = "r1\t5\t30\ttRNA-Ala\t+\n"
        infile = io.StringIO(header + data)
        outfile = io.BytesIO()
        process_modkit_tsv(infile, outfile, OFFSET_5P, OFFSET_3P, REF_LENGTHS)
        lines = outfile.getvalue().decode().strip().split("\n")
        fields = lines[1].split("\t")
        assert fields[2] == "7"  # 30 - 24 + 1

    def test_adapter_filtering(self):
        """Adapter positions should be filtered."""
        header = "read_id\tfwd_pos\tref_position\tchrom\tmod_strand\n"
        data = "r1\t5\t10\ttRNA-Ala\t+\n"  # pos 10 < offset_5p=24
        infile = io.StringIO(header + data)
        outfile = io.BytesIO()
        process_modkit_tsv(infile, outfile, OFFSET_5P, OFFSET_3P, REF_LENGTHS)
        lines = outfile.getvalue().decode().strip().split("\n")
        assert len(lines) == 1  # Only header
