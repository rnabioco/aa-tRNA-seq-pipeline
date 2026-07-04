"""
Unit tests for get_charging_table.py

Tests ML tag extraction from BAM files.
"""

import gzip
from array import array

import pysam
import pytest

from get_charging_table import extract_tag


class TestExtractTag:
    """Tests for extract_tag function."""

    def test_extracts_ml_tag(self, temp_dir):
        """Should extract ML tag values to TSV."""
        input_bam = temp_dir / "input.bam"
        output_tsv = temp_dir / "output.tsv"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "tRNA-Ala-AGC-1-1", "LN": 100}],
        }

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read.set_tag("ML", array("B", [220]))
            outf.write(read)

        pysam.index(str(input_bam))

        extract_tag(str(input_bam), str(output_tsv), "ML")

        # Check output
        with open(output_tsv) as f:
            lines = f.readlines()

        assert len(lines) == 2  # Header + 1 read
        assert "read_id\ttRNA\tcharging_likelihood\n" == lines[0]
        assert "read1\ttRNA-Ala-AGC-1-1\t220\n" == lines[1]

    def test_handles_gzip_output(self, temp_dir):
        """Should write gzipped output when filename ends in .gz."""
        input_bam = temp_dir / "input.bam"
        output_tsv = temp_dir / "output.tsv.gz"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "tRNA-Ala-AGC-1-1", "LN": 100}],
        }

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read.set_tag("ML", array("B", [200]))
            outf.write(read)

        pysam.index(str(input_bam))

        extract_tag(str(input_bam), str(output_tsv), "ML")

        # Check gzipped output
        with gzip.open(str(output_tsv), "rt") as f:
            lines = f.readlines()

        assert len(lines) == 2
        assert "read1" in lines[1]

    def test_multiple_reads(self, temp_dir):
        """Should handle multiple reads."""
        input_bam = temp_dir / "input.bam"
        output_tsv = temp_dir / "output.tsv"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [
                {"SN": "tRNA-Ala-AGC-1-1", "LN": 100},
                {"SN": "tRNA-Gly-GCC-1-1", "LN": 100},
            ],
        }

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            for i, (ref_id, ml_val) in enumerate([(0, 220), (0, 150), (1, 255)]):
                read = pysam.AlignedSegment()
                read.query_name = f"read{i}"
                read.query_sequence = "A" * 100
                read.flag = 0
                read.reference_id = ref_id
                read.reference_start = 0
                read.cigartuples = [(0, 100)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 100)
                read.set_tag("ML", array("B", [ml_val]))
                outf.write(read)

        pysam.index(str(input_bam))

        extract_tag(str(input_bam), str(output_tsv), "ML")

        with open(output_tsv) as f:
            lines = f.readlines()

        assert len(lines) == 4  # Header + 3 reads

    def test_skips_unmapped_reads(self, temp_dir):
        """Unmapped reads should be excluded."""
        input_bam = temp_dir / "input.bam"
        output_tsv = temp_dir / "output.tsv"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "tRNA-Ala-AGC-1-1", "LN": 100}],
        }

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
            read1.set_tag("ML", array("B", [200]))
            outf.write(read1)

            # Unmapped read with ML tag
            read2 = pysam.AlignedSegment()
            read2.query_name = "unmapped"
            read2.query_sequence = "A" * 50
            read2.flag = 4
            read2.reference_id = -1
            read2.reference_start = -1
            read2.query_qualities = pysam.qualitystring_to_array("I" * 50)
            read2.set_tag("ML", array("B", [150]))
            outf.write(read2)

        pysam.index(str(input_bam))

        extract_tag(str(input_bam), str(output_tsv), "ML")

        with open(output_tsv) as f:
            lines = f.readlines()

        # Only mapped read should appear (plus header)
        assert len(lines) == 2
        assert "mapped" in lines[1]

    def test_retains_ml_zero_reads(self, temp_dir):
        """ML==0 is a valid, maximally-confident *uncharged* call and must be
        written to the table. Regression test for a bug where the write gate
        used `if tag_value` (falsy at 0), silently dropping ML==0 reads and
        biasing charging fraction upward."""
        input_bam = temp_dir / "input.bam"
        output_tsv = temp_dir / "output.tsv"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "tRNA-Ala-AGC-1-1", "LN": 100}],
        }

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            # One confidently-uncharged read (ML=0) and one charged (ML=255)
            for name, ml_val in [("uncharged_zero", 0), ("charged_max", 255)]:
                read = pysam.AlignedSegment()
                read.query_name = name
                read.query_sequence = "A" * 100
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 0
                read.cigartuples = [(0, 100)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 100)
                read.set_tag("ML", array("B", [ml_val]))
                outf.write(read)

        pysam.index(str(input_bam))

        extract_tag(str(input_bam), str(output_tsv), "ML")

        with open(output_tsv) as f:
            lines = f.readlines()

        # Header + both reads (the ML=0 read must NOT be dropped)
        assert len(lines) == 3
        assert any(
            line.startswith("uncharged_zero\t") and line.rstrip().endswith("\t0")
            for line in lines[1:]
        )

    def test_skips_and_warns_on_multielement_tag(self, temp_dir, capsys):
        """A multi-element tag array (e.g. the dorado mod-base ML tag) is not a
        single charging score and is skipped, but the skip must be reported on
        stderr rather than dropped silently (silent drops bias the charging
        fraction like the ML==0 bug did)."""
        input_bam = temp_dir / "input.bam"
        output_tsv = temp_dir / "output.tsv"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "tRNA-Ala-AGC-1-1", "LN": 100}],
        }

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            # One scalar-usable read and one multi-element read (2 mod probs)
            scalar_read = pysam.AlignedSegment()
            scalar_read.query_name = "scalar"
            scalar_read.query_sequence = "A" * 100
            scalar_read.flag = 0
            scalar_read.reference_id = 0
            scalar_read.reference_start = 0
            scalar_read.cigartuples = [(0, 100)]
            scalar_read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            scalar_read.set_tag("ML", array("B", [200]))
            outf.write(scalar_read)

            multi_read = pysam.AlignedSegment()
            multi_read.query_name = "multi"
            multi_read.query_sequence = "A" * 100
            multi_read.flag = 0
            multi_read.reference_id = 0
            multi_read.reference_start = 0
            multi_read.cigartuples = [(0, 100)]
            multi_read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            multi_read.set_tag("ML", array("B", [10, 250]))
            outf.write(multi_read)

        pysam.index(str(input_bam))

        extract_tag(str(input_bam), str(output_tsv), "ML")

        with open(output_tsv) as f:
            lines = f.readlines()

        # Header + only the scalar read; the multi-element read is skipped
        assert len(lines) == 2
        assert lines[1].startswith("scalar\t")
        assert not any(line.startswith("multi\t") for line in lines[1:])

        # ...but the skip is reported, not silent
        err = capsys.readouterr().err
        assert "skipped 1" in err

    def test_different_tag(self, temp_dir):
        """Should work with different tag names."""
        input_bam = temp_dir / "input.bam"
        output_tsv = temp_dir / "output.tsv"

        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"SN": "tRNA-Ala-AGC-1-1", "LN": 100}],
        }

        with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A" * 100
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.cigartuples = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read.set_tag("CL", array("B", [180]))
            outf.write(read)

        pysam.index(str(input_bam))

        extract_tag(str(input_bam), str(output_tsv), "CL")

        with open(output_tsv) as f:
            lines = f.readlines()

        assert len(lines) == 2
        assert "180" in lines[1]
